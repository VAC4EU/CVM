# Program Information  ----------------------------------------------------
#
# Program:      step_12_2_run_scri.R 
# Author:       Svetlana Belitser
# Description:  calls functions which runs an SCRI on specified datasets 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival" 
#               input:  g_output/scri/*  
#               output:  g_output/scri/*  
#
#               parameters: in 07_scri_inputs.R 
#  
#               function: p_macro/scri_tools.R                    
#               function: p_macro/table1.R  
#
#
#               function: scri_create_input
#               function: add_to_report_list
#               function: add_to_models_list
#
#  

####       print_during_running <- F    #  T or F


# Housekeeping  -----------------------------------------------------------
# install and load packages

if(!any(ls()=="thisdir"))   thisdir   <- getwd()
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"/g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"/g_output/")

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp, "scri")),           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,"/log_files/scri")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(dirtemp, "models")),         showWarnings = FALSE, recursive = TRUE)

for (subpop in subpopulations_non_empty) {
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE, direxp, direxpsubpop[[subpop]])
  
  # SCCS output_directory to export  
  sdr0 <- paste0(thisdirexp, "scri/")
  dir.create(sdr0, showWarnings = FALSE, recursive = TRUE)
  
  # SCCS output_directory for models 
  sdr_models0 <- paste0(dirtemp, "models/", "scri/") 																							 
  dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
  
  
  # Import Data -------------------------------------------------------------
  load(paste0(dirtemp, "scri/", intermediate_data, suffix[[subpop]], ".RData"))
  temp_name<-get(paste0(intermediate_data, suffix[[subpop]]))
  rm(list=paste0(intermediate_data, suffix[[subpop]]))
  assign(intermediate_data, as.data.frame(temp_name))
  rm(temp_name)
  
  
  if(F){  # this is only for testing
    
    if(!any(ls()=="thisdirexp")) thisdirexp <- paste0(thisdir,"/g_export/")
    dir.create(file.path(paste0(thisdirexp, "scri")), showWarnings = FALSE, recursive = TRUE)

    # SCCS output_directory  
    sdr0 <- paste0(thisdirexp, "scri/")
    dir.create(sdr0, showWarnings = FALSE, recursive = TRUE)
    
    # SCCS output_directory for models 
    sdr_models0 <- paste0(diroutput, "scri/")
    dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
    
    # Import Data -------------------------------------------------------------
    load(paste0(dirtemp, "scri/", intermediate_data, ".RData"))
    scri_data_extract <- eval(parse(text = intermediate_data))
    
    
    library("survival")
    library("MASS")
    library("qpdf")
    
  } 
  
  scri_input <- scri_data_extract
  
  dap <- ifelse( any(names(scri_input)=="DAP"), scri_input$DAP[1], "")
  
  
  
  #############   SCRI models ############################
  #
  #
  #
  
  ##########
  #              (vax2 takes precedence over vax1)
  #
  #   the risk window of dose 2 takes precedence over the risk window of dose 1
  
  old_width = options(width=300)
  
  paral   <- T    # if T ==> library(parallel) is started in function 'scri_strata'
  n_cores <- NA   # Don't define it! Define it only if it doesn't work on its own! number of cores/threads to use.
  
  print_during_running <- F
  plot_during_running  <- F  
  CI_draw              <- T
  
  # which part to run:
  lmain <- T
  ldist <- T
  
  onlyPfizer <- F
  
  
  paral_vars <- c("lab_pre", "vax1_end_before_vax2")
  
  
  ##########################################
  #
  #      definitions of risk windows:
  #
  
  lab_pre  <- c( "pre-exposure",  "buffer" )
  
  
  vax1_end_before_vax2 <- substitute(pmin(days_vax2-1, study_exit_days, na.rm=T))
  vax2_end_before_vax3 <- substitute(pmin(days_vax3-1, study_exit_days, na.rm=T))
  
  rws_def_2vax_7  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),  lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                       v1      = list( t0 = days_vax1, cuts = c( 0, 1, 8, 15, 29, 61 ),    tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                       v2      = list( t0 = days_vax2, cuts = c( 0, 1, 8, 15, 29, 61 ),    tend = study_exit_days,       lab = "dose 2 "   ))
                                #all      = list( t0 = study_entry_days, tend=study_exit_days,     lab = "all", lab_add_interval=F   )
  )
  
  rws_def_2vax_14  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),       lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                        v1      = list( t0 = days_vax1, cuts = c( 0, 1, 15, 29, 61 ),  tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                        v2      = list( t0 = days_vax2, cuts = c( 0, 1, 15, 29, 61 ),  tend = study_exit_days,       lab = "dose 2 "   ))
                                 #all      = list( t0 = study_entry_days, tend=study_exit_days,     lab = "all", lab_add_interval=F   )
  )
  rws_def_2vax_28 <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),  lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                       v1      = list( t0 = days_vax1, cuts = c( 0, 1, 29 ),    tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                       v2      = list( t0 = days_vax2, cuts = c( 0, 1, 29 ),    tend = study_exit_days,       lab = "dose 2 "   ))
                                #all      = list( t0 = study_entry_days, tend=study_exit_days,     lab = "all", lab_add_interval=F   )
  )
  
  rws_def_3vax_7  <- substitute(list(prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),               lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                     v1      = list( t0 = days_vax1, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                     v2      = list( t0 = days_vax2, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = vax2_end_before_vax3,  lab = "dose 2 "   ), 
                                     v3      = list( t0 = days_vax3, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,       lab = "dose 3 "   ))
  )
  rws_def_3vax_14  <- substitute(list(prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),     lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                      v1      = list( t0 = days_vax1, cuts = c( 0, 1, 15, 29, 61 ), tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                      v2      = list( t0 = days_vax2, cuts = c( 0, 1, 15, 29, 61 ), tend = vax2_end_before_vax3,  lab = "dose 2 "   ), 
                                      v3      = list( t0 = days_vax3, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,       lab = "dose 3 "   ))
  )
  rws_def_3vax_28 <- substitute(list(prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),               lab  = lab_pre , lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
                                     v1      = list( t0 = days_vax1, cuts = c( 0, 1, 29 ), tend = vax1_end_before_vax2,  lab = "dose 1 "   ), 
                                     v2      = list( t0 = days_vax2, cuts = c( 0, 1, 29 ), tend = vax2_end_before_vax3,  lab = "dose 2 "   ), 
                                     v3      = list( t0 = days_vax3, cuts = c( 0, 1, 29 ), tend = study_exit_days,       lab = "dose 3 "   ))
  )
  
  
  # with overlap & the last intervals last till the end of observation !!!
  rws_def_3vax_overlap_7  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),          lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                               v1      = list( t0 = days_vax1, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                               v2      = list( t0 = days_vax2, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    ), 
                                               v3      = list( t0 = days_vax3, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 3 ", no_last_interval=F    ))
  )
  rws_def_3vax_overlap_14  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),      lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                                v1      = list( t0 = days_vax1, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                                v2      = list( t0 = days_vax2, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    ), 
                                                v3      = list( t0 = days_vax3, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 3 ", no_last_interval=F    ))
  )
  rws_def_3vax_overlap_28 <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),          lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                               v1      = list( t0 = days_vax1, cuts = c( 0, 1, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                               v2      = list( t0 = days_vax2, cuts = c( 0, 1, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    ), 
                                               v3      = list( t0 = days_vax3, cuts = c( 0, 1, 29, 61 ), tend = study_exit_days,  lab = "dose 3 ", no_last_interval=F    ))
  )
  rws_def_2vax_overlap_7  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),          lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                               v1      = list( t0 = days_vax1, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                               v2      = list( t0 = days_vax2, cuts = c( 0, 1, 8, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    )) 
  )
  rws_def_2vax_overlap_14  <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),      lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                                v1      = list( t0 = days_vax1, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                                v2      = list( t0 = days_vax2, cuts = c( 0, 1, 15, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    )) 
  )
  rws_def_2vax_overlap_28 <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),          lab  = lab_pre , no_last_interval=T   ), # also can be added:    add_interval = T 
                                               v1      = list( t0 = days_vax1, cuts = c( 0, 1, 29, 61 ), tend = study_exit_days,  lab = "dose 1 ", no_last_interval=F    ), 
                                               v2      = list( t0 = days_vax2, cuts = c( 0, 1, 29, 61 ), tend = study_exit_days,  lab = "dose 2 ", no_last_interval=F    )) 
  )
  
  # with overlap & with cutting off the last intervals, i.e., not until the end of observation !!!
  rws_def_3vax_overlap_7_bounded <- substitute(list(  prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),          lab = lab_pre ,  no_last_interval=T   ), # also can be added:    add_interval = T 
                                                      v1      = list( t0 = days_vax1, cuts = c( 0, 1, 8, 15, 29, 61 ), lab = "dose 1 ", no_last_interval=T    ), 
                                                      v2      = list( t0 = days_vax2, cuts = c( 0, 1, 8, 15, 29, 61 ), lab = "dose 2 ", no_last_interval=T    ), 
                                                      v3      = list( t0 = days_vax3, cuts = c( 0, 1, 8, 15, 29, 61 ), lab = "dose 3 ", no_last_interval=T    ))
  )
  
  
  
  ##########################################
  #
  # define distance categories for number of days between vax1 and vax2 ('dist12'), between vax2 and vax3 ('dist32')
  #
  distances <- list( list( dist12 = c(-Inf,-1,    3*7,  Inf),
                           dist23 = c(-Inf,-1,   10*7,  Inf),
                           name   = "_dist3w"                     ),
                     list( dist12 = c(-Inf,-1,    4*7,  Inf),
                           dist23 = c(-Inf,-1,   12*7,  Inf),
                           name   = "_dist4w"                    ),
                     list( dist12 = c(-Inf,-1,    5*7,  Inf),
                           dist23 = c(-Inf,-1,   16*7,  Inf),
                           name   = "_dist5w"                    ),
                     list( dist12 = c(-Inf,-1,    6*7,  Inf),
                           dist23 = c(-Inf,-1,   20*7,  Inf),
                           name   = "_dist6w"                    ),
                     list( dist12 = c(-Inf,-1,   10*7,  Inf),
                           dist23 = c(-Inf,-1,   24*7,  Inf),
                           name   = "_dist10w"                    ),
                     list( dist12 = c(-Inf,-1,   11*7,  Inf),
                           dist23 = c(-Inf,-1,   25*7,  Inf),
                           name   = "_dist11w"                    ),
                     list( dist12 = c(-Inf,-1,   12*7,  Inf),
                           dist23 = c(-Inf,-1,   26*7,  Inf),
                           name   = "_dist12w"                    )
  )

  
 
  
  ##########################################
  #
  #      definitions of time_dependent variables :
  #
  
  ########## variable 'brand'   ##############
  
  
  brands <- unique(unlist(scri_input[,paste0("type_vax",1:nvax)]))
  brands <- brands[!is.na(brands)]
  tmp <- match(brands,c("Pfizer","Moderna","AstraZeneca","J&J","JJ"))
  brands <- c( brands[tmp[!is.na(tmp)]],  brands[is.na(tmp)] )  
  brands
  
  # create new type_vax1_short, type_vax2_short, ... : as substring(..., 1, 5 )
  for(iv in 1:nvax)
    scri_input[paste0("type_vax",iv,"_short")] <- format(substring(scri_input[,paste0("type_vax",iv)], 1, 5), width=5)
  
  
  brand_2v_def <- substitute(list( splits_names  = "vax", 
                                splits        = cbind.data.frame( days_vax1, days_vax2),  
                                lab           = c("no_vax","dose 1","dose 2"),
                                ref           = "pre-", 
                                change        = list( new_name         = "br",
                                                      #add_begin_sep = "_",
                                                      from             = "vax",
                                                      condition_var    = "vax",
                                                      replace       = list( c(value="dose 1", var_name="type_vax1_short"),
                                                                            c(value="dose 2", var_name="type_vax2_short"),
                                                                            c(value="dose 3", var_name="type_vax3_short")  ),
                                                      separate_variables = c("d1","d2","d3")
                                )
  ))
  
  
  brand_distance_2v_def <- substitute(list( splits_names  = "vax", 
                                         splits        = cbind.data.frame( days_vax1, days_vax2),  
                                         lab           = c("no_vax","dose 1","dose 2"),
                                         ref           = "pre-", 
                                         change        = list( 
                                           list( new_name         = "br_dist",
                                                 from             = "vax",
                                                 condition_var    = "vax",
                                                 #add_begin_sep    = "_",
                                                 replace       = list( c( value="dose 1", var_name="type_vax1_short" ),
                                                                       c( value="dose 2", var_name="type_vax2_short" ),
                                                                       c( value="dose 3", var_name="type_vax3_short" )  ),
                                                 separate_variables = c("d1","d2","d3"),
                                                 add_dist =list( c( value="dose 1", var_name="vax01_distance" ),
                                                                 c( value="dose 2", var_name="vax12_distance" ),
                                                                 c( value="dose 3", var_name="vax23_distance" )  )
                                           ))
  ))
  
  
  distance_2v_def <- substitute(list( splits_names  = "vax", 
                                   splits        = cbind.data.frame( days_vax1, days_vax2),  
                                   lab           = c("no_vax","dose 1","dose 2"),
                                   ref           = "pre-", 
                                   change        = list( new_name         = "dist",
                                                         from             = "vax",
                                                         condition_var    = "vax",
                                                         #add_begin_sep    = "_",
                                                         replace       = list( c( value="dose 1", var_name="vax01_distance" ),
                                                                               c( value="dose 2", var_name="vax12_distance" ),
                                                                               c( value="dose 3", var_name="vax23_distance" )  
                                                         ) )
  ))
  
  
  
  brand_3v_def <- substitute(list( splits_names  = "vax", 
                                splits        = cbind.data.frame( days_vax1, days_vax2, days_vax3),  
                                lab           = c("no_vax","dose 1","dose 2","dose 3"),
                                ref           = "pre-", 
                                change        = list( new_name         = "br",
                                                      #add_begin_sep = "_",
                                                      from             = "vax",
                                                      condition_var    = "vax",
                                                      replace       = list( c(value="dose 1", var_name="type_vax1_short"),
                                                                            c(value="dose 2", var_name="type_vax2_short"),
                                                                            c(value="dose 3", var_name="type_vax3_short")  ),
                                                      separate_variables = c("d1","d2","d3")
                                )
  ))
  
  
  brand_distance_3v_def <- substitute(list( splits_names  = "vax", 
                                         splits        = cbind.data.frame( days_vax1, days_vax2, days_vax3),  
                                         lab           = c("no_vax","dose 1","dose 2","dose 3"),
                                         ref           = "pre-", 
                                         change        = list( 
                                           list( new_name         = "br_dist",
                                                 from             = "vax",
                                                 condition_var    = "vax",
                                                 #add_begin_sep    = "_",
                                                 replace       = list( c( value="dose 1", var_name="type_vax1_short" ),
                                                                       c( value="dose 2", var_name="type_vax2_short" ),
                                                                       c( value="dose 3", var_name="type_vax3_short" )  ),
                                                 separate_variables = c("d1","d2","d3"),
                                                 add_dist =list( c( value="dose 1", var_name="vax01_distance" ),
                                                                 c( value="dose 2", var_name="vax12_distance" ),
                                                                 c( value="dose 3", var_name="vax23_distance" )  )
                                           ))
  ))
  
  
  distance_3v_def <- substitute(list( splits_names  = "vax", 
                                   splits        = cbind.data.frame( days_vax1, days_vax2, days_vax3),  
                                   lab           = c("no_vax","dose 1","dose 2","dose 3"),
                                   ref           = "pre-", 
                                   change        = list( new_name         = "dist",
                                                         from             = "vax",
                                                         condition_var    = "vax",
                                                         #add_begin_sep    = "_",
                                                         replace       = list( c( value="dose 1", var_name="vax01_distance" ),
                                                                               c( value="dose 2", var_name="vax12_distance" ),
                                                                               c( value="dose 3", var_name="vax23_distance" )  
                                                         ) )
  ))
  
  
  
  
  ########## variable 'cumulative brand'   ##############
  
  if(F){
    scri_input$type_vax1_d1  <- paste0( "d1:",substring(scri_input$type_vax1_short,1,4) )
    scri_input$type_vax12[ !is.na(scri_input$type_vax2_short)] <- paste0( scri_input$type_vax1_d1, "_d2:", substring(scri_input$type_vax2_short,1,4))[!is.na(scri_input$type_vax2_short)]
    scri_input$type_vax123[!is.na(scri_input$type_vax3_short)] <- paste0( scri_input$type_vax12,   "_d3:", substring(scri_input$type_vax3_short,1,4))[!is.na(scri_input$type_vax3_short)]
    
    brand_123cum_def <- substitute(list( 
      splits_names  = "vax123", 
      splits        = cbind.data.frame( days_vax1, days_vax2, days_vax3),  
      lab           = c("no_vax","dose 1","dose 2","dose 3"),
      change        = list( name          = "brand_123cum",
                            #add_begin_sep = "_",
                            replace       = list( c(value="dose 1", var_name="type_vax1_d1"),
                                                  c(value="dose 2", var_name="type_vax12"),
                                                  c(value="dose 3", var_name="type_vax123") ),
                            separate_variables = c("d1","d2","d3")
      ),     
      ref           = "pre-"))
    
    brand_12cum_def <- substitute(list( 
      splits_names  = "vax12", 
      splits        = cbind.data.frame( days_vax1, days_vax2 ),  
      lab           = c("no_vax","dose 1","dose 2"),
      change        = list( name          = "brand_12cum",
                            #add_begin_sep = "_",
                            replace       = list( c(value="dose 1", var_name="type_vax1_d1"),
                                                  c(value="dose 2", var_name="type_vax12") ),
                            separate_variables = c("d1","d2","d3")
      ),     
      ref           = "pre-"))
  }
  
  #
  ########################################################
  
  ########################################
  #
  # specify calendar time interval for adjusting:
  #
  lengths <- c( 10, 14, 21,30)
  #lengths <- c( 7,14,10, 21,30)
  #lengths <- c( 3,7,14,10, 21,30)
  starts  <- c(-1,-5,-8,-10,-1, -3)
  time_seq <- vector("list",length=length(lengths))
  names(time_seq) <- paste0("period",lengths,"d")
  for(i in 1:length(lengths))
    time_seq[[i]] <- seq(min(scri_input$study_entry_days,na.rm=T)-starts[i],max(scri_input$study_exit_days,na.rm=T)+lengths[i]-1,by=lengths[i])
  
  ########################################
  #
  #    define order of coefficient labels:
  #
  
  
  lab_orders <- list(  
    c("F", "M" ),
    c("sex0","sex:0","sexF","sex:F", "sex1", "sex:1" , "sexM", "sex:M" ),
    paste0("age",c("(-1,30]","(30,40]","(30,50]","(30,60]","(30,Inf]",">=30",">30","(40,50]","(50,60]","(50,65]","(50,Inf]",">=50",">50",">=60","(60,Inf]", ">65" )),
    #c("(-1,30]","(30,40]","(30,50]","(30,60]","(30,Inf]",">=30",">30","(40,50]","(50,60]","(50,65]","(50,Inf]",">=50",">50",">=60","(60,Inf]", ">65" ),
    c("d1:","d2:","d3:" ),
    #c("Pfi", "Mod", "Ast",  "JJ","J&J" ),
    #c("Pfizer", "Moderna", "AstraZeneca", "JJ","J&J" ),
    c("Pfi","no_Pf","no Pf", "Mod","no_Mo","no Mo", "Ast","no_AZ","no AZ",  "JJ","J&J","no_JJ","no_J&","no JJ","no_J&" ),
    c("Pfizer","no_Pfizer","no Pfizer", "Moderna","no_Moderna","no Moderna", "AstraZeneca","no_AZ", "no AZ", "JJ","J&J", "no_JJ","no_J&J", "no JJ","no J&J" ),
    c("pre-","buf", "dose 1","dose 2", "dose 3" ),
    c("buf", "dose 1 pre-", paste0( rep(paste0("dose ",1:3),each=10), rep(c("pre-"," pre-","<"," <","("," (","["," [",">"," >"),3) ) ),
    c("[-90;-30],[-29;-1],[0;0]","[1;7]","[1;14]","[1;28]","[8;14]","[15;28]",">28","[29;60]",">60","[61;180]", ">180"),
    paste0(":",c("(-1,30]","(30,60]","(60,Inf]","(60, Inf]")),
    paste0(":",c("(-Inf,-1]","(-1,70]","(70,Inf]")),
    paste0(":",c("(-\U221E,-1]","(-Inf,-1]","(-1,21]","(-1,30]","(21,35]","(30,60]","(35,56]","(56,84]","(60,Inf]","(60, Inf]","(84, Inf]","(84, \U221E]")),
    paste0(":",c("(-\U221E,-1]","(-Inf,-1]","(-1,175]","(175,Inf]","(175, \U221E]"))
  )
  
  
  # during testing: may use only one vector to adjust for calendar time, for example, time_seq[5]:
  # time_seq <- time_seq[5]

  # colors:
  col_list <- c("red", "green3", "orange",  "deepskyblue", "magenta2", "gray", "cyan2","chocolate1" ) 
  
  
  # outcome variables:
  ae_events <- c("myocarditis","pericarditis")
  #ae_events <- c("myocarditis","pericarditis","myopericarditis")
  
  
  # create covid selection variables: no_covid_before_myocard_30d, no_covid_before_pericar_30d, no_covid_before_myoperi_30d, ...
  for(iae in ae_events){ 
    scri_input[,paste0("no_covid_before_",substring(iae,1,7),"_30d")] <- scri_input$covid19_date > ( scri_input[,paste0(iae,"_date")] - 30 )  &  !is.na(scri_input$covid19_date) & !is.na(scri_input[,paste0(iae,"_date")])
    scri_input[is.na(scri_input$covid19_date),paste0("no_covid_before_",substring(iae,1,7),"_30d")] <- T
  }                                      
  {
    cat("\n\ntable for variables:\n")
    print(c("covid",paste0("no_covid_before_",substring(ae_events,1,7),"_30d")))
    table1(scri_input[,c("covid",paste0("no_covid_before_",substring(ae_events,1,7),"_30d"))] )
  }
  
  
  # select only those who didnot have COVID before the event and 30 days extra to make sure that 
  # if someone has COVID then there is no connection to the event.
  # create covid selection variables: no_covid_before_myocard_30d, no_covid_before_pericar_30d, no_covid_before_myoperi_30d, ...
  for(iae in ae_events){ 
    scri_input[,paste0("no_covid_before_",substring(iae,1,7),"_plus30d")] <- scri_input$covid19_date > ( scri_input[,paste0(iae,"_date")] + 30 )  &  !is.na(scri_input$covid19_date) & !is.na(scri_input[,paste0(iae,"_date")])
    scri_input[is.na(scri_input$covid19_date),paste0("no_covid_before_",substring(iae,1,7),"_plus30d")] <- T
  }                                      
  {
    cat("\n\ntable for variables:\n")
    print(c("covid",paste0("no_covid_before_",substring(ae_events,1,7),"_plus30d")))
    table1(scri_input[,c("covid",paste0("no_covid_before_",substring(ae_events,1,7),"_plus30d"))] )
  }
  
  
  
  
  
  #################################################
  #################################################
  #
  #  run scri analysis
  # 
  
  scri_input0 <- scri_input
  
  
  
  if(onlyPfizer){   # run analysis only if at least one of the doses is Pfizer
    scri_input0 <- scri_input0[ scri_input0$type_vax1=="Pfizer" | is.na(scri_input0$type_vax2) | (!is.na(scri_input0$type_vax2) & scri_input0$type_vax2=="Pfizer"),]
    scri_input0$type_vax1[!is.na(scri_input0$type_vax1) & scri_input0$type_vax1!="Pfizer"] <- "no_Pfizer"
    scri_input0$type_vax2[!is.na(scri_input0$type_vax2) & scri_input0$type_vax2!="Pfizer"] <- "no_Pfizer"
    
    # create new type_vax1_short, type_vax2_short, ... : as substring(..., 1, 5 )
    for(iv in 1:nvax)
      scri_input0[paste0("type_vax",iv,"_short")] <- format(substring(scri_input0[,paste0("type_vax",iv)], 1, 5), width=5)
  }
  
  
  
  for(icovid in c( "all_data", "no_covid_before_event_plus30d", "no_covid_start_control_rw" ) ){
  #for(icovid in c( "all_data", "no_covid_before_event_plus30d", "no_covid_before_event_30d", "no_covid_start_control_rw" ) ){
      
    if( length(grep("covid",tolower(icovid)))>0 & all(tolower(names(scri_input0)) != c("covid19_date")) ) next 
    
    # SCCS output_directory  
    sdr <- paste0(sdr0, icovid,"/")
    dir.create(sdr, showWarnings = FALSE, recursive = TRUE)
    
    # SCCS output_directory  
    sdr_models <- paste0(sdr_models0, icovid,"/")
    dir.create(sdr_models, showWarnings = FALSE, recursive = TRUE)
    
    scri_input <- scri_input0
    
    for(iae in ae_events){
      if( icovid =="all_data"                      ) scri_input[,paste0("cond_covid_",substring(iae,1,7))] <- rep(T,nrow(scri_input))
      if( icovid =="no_covid_before_event_plus30d" ) scri_input[,paste0("cond_covid_",substring(iae,1,7))] <- scri_input[,paste0("no_covid_before_",substring(iae,1,7),"_plus30d") ]
      if( icovid =="no_covid_before_event_30d"     ) scri_input[,paste0("cond_covid_",substring(iae,1,7))] <- scri_input[,paste0("no_covid_before_",substring(iae,1,7),"_30d") ]
      if( icovid =="no_covid_start_control_rw"     ) scri_input[,paste0("cond_covid_",substring(iae,1,7))] <- scri_input$covid==0
    }  
    
    covid_select_text <- switch (icovid,
                                 "all_data"                      = "alldata",
                                 "no_covid_before_event_plus30d" = "covidplus30d",
                                 "no_covid_before_event_30d"     = "covid30d",
                                 "no_covid_start_control_rw"     = "covidstart"
    )
    
    vax_priority <- paste0("_vax2prior_",covid_select_text)
    
    
    
    if(lmain){
      ########### no strata: #############
      #
      
      # for(iae in ae_events)
      #   brand_images( plot_data, ae_event=iae, tit="")
      
      cat("\n\nAnalysis for ",icovid," together.\n\n")
      
      models_list <- list()
      report_list <- list()
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        
        first_plot <- T
        ae_event_first <- T  
        
        for(iii in ifelse(nvax>=3,6,3):1){
          
          if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
          if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
          if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
          
          if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
          if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
          if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
          
          #vax_priority <- paste0("_vax2prior_",covid_select_text)
          specif_name  <- "_no_split" 
          
          global_name  <- paste0( vax_priority, specif_name )
          output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
          cat(paste(output_name, format(Sys.time()),"\n"))
          
          formula_text <-  "~ lab"
          
          res <- scri_strata( output_name  = output_name, 
                              formula_text = formula_text,       time_seq = time_seq, 
                              event_time = paste0(iae,"_days"), event = iae, id="person_id",
                              rws          = rws_def,
                              start_obs    = "study_entry_days", end_obs = "study_exit_days",
                              data         = scri_input[cond_iae,],
                              image_plots = ae_event_first,
                              lab_orders = lab_orders,
                              paral = T, paral_vars = paral_vars, n_cores = n_cores,
                              lprint = print_during_running,
                              global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                              CI = CI_draw
          )
          first_plot <- F
          ae_event_first <- F
          
          report_list <- add_to_report_list(res$tabs,     output_name)
          models_list <- add_to_models_list(res$scri_all, output_name)
        }
      }
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      
      ###########
      #  save report_list and model_list
      save_results(global_name, report_list, models_list)
      #
      #######################################################################################
      
      
      
      ########### per brand:  all brands in one model #############
      
      cat("\n\nAnalysis for ",icovid," per brand.\n\n")
      
      models_list <- list()
      report_list <- list()
      
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        first_plot <- T
        ae_event_first <- T
        
        for(iii in ifelse(nvax>=3,6,3):1){
          
          
          if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
          if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
          if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
          
          if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
          if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
          if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
          
          specif_name  <-"_brands" 
          
          global_name  <- paste0( vax_priority, specif_name )
          output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
          cat(paste(output_name, format(Sys.time()),"\n"))
          
          formula_text <-  "~ br:lab"
          
          if(nvax>=3) time_dep <- list( brand_3v_def )
          else        time_dep <- list( brand_2v_def )
          
          res <- scri_strata( output_name  = output_name, 
                              formula_text = formula_text,       time_seq = time_seq, 
                              event_time = paste0(iae,"_days"), event = iae, id="person_id",
                              rws          = rws_def,
                              time_dep     = time_dep,
                              start_obs    = "study_entry_days", end_obs = "study_exit_days",
                              data         = scri_input[cond_iae,],
                              image_plots  = ae_event_first, image_brand=T,
                              lab_orders   = lab_orders,
                              paral = T, paral_vars = paral_vars, n_cores = n_cores,
                              lprint       = print_during_running,
                              global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                              CI = CI_draw
          )
          first_plot <- F
          ae_event_first <- F
          
          report_list <- add_to_report_list(res$tabs,     output_name)
          models_list <- add_to_models_list(res$scri_all, output_name)
        }
      }  
      
      
      
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste( report_list[[istr]][[1]][[1]][1,"event"], formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      ####
      #  save report_list and model_list
      save_results(global_name, report_list, models_list)
      #
      #######################################################################################
      
      
      
      
      #######################################################################################
      
      ########### per sex & age30: one model #############
      #
      
      cat("\n\nAnalysis for ",icovid," per sex, age30, brand.\n\n")
      
      if(all(names(scri_input) !="age30"))
        scri_input$age30 <- paste0("age",as.character(cut(scri_input$age_at_study_entry, c(-1,30,Inf))))
      #levels(scri_input$age30) <- gsub( "Inf","\U221E", levels(scri_input$age30),fixed=T)
      
      scri_input$sex_age30 <- paste0("sex:",scri_input$sex, " ", scri_input$age30)
      
      models_list <- list()
      report_list <- list()
      
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        for(istr in unique(scri_input$sex_age30[cond_iae]) ){
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <- istr   #  "_sex_age30" 
            
            global_name0 <- paste0( vax_priority, "_sex_age30" )
            global_name  <- paste0( vax_priority, specif_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ sex_age30:br:lab"
            
            if(nvax>=3) time_dep <- list( brand_3v_def )
            else        time_dep <- list( brand_2v_def )
            
            res <- scri_strata(  output_name  = output_name, 
                                 formula_text = formula_text,       time_seq = time_seq, 
                                 event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                 rws          = rws_def,
                                 time_dep     = time_dep, 
                                 #combine_vars =  c("sex","age4"), 
                                 start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                 data         = scri_input[scri_input$sex_age30==istr & cond_iae, ],
                                 nvax         = min(3,nvax),
                                 image_plots  = ae_event_first, image_brand=T, image_tit=istr,
                                 lab_orders   = lab_orders,
                                 paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                 lprint       = print_during_running,
                                 global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                 CI = CI_draw
            )
            ae_event_first <- F
            first_plot <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
      }  
      
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      ####
      #  save report_list and model_list
      save_results(global_name0, report_list, models_list)
      #
      #######################################################################################
      
      
      
      
      
      
      
      ########### per age30_50: one model #############
      #
      
      cat("\n\nAnalysis for ",icovid," per age30_50, brand.\n\n")
      
      if(all(names(scri_input) !="age30_50"))
        scri_input$age30_50 <- paste0("age",as.character(cut(scri_input$age_at_study_entry, c(-1,30,50,Inf))))
      #levels(scri_input$age30_50) <- gsub( "Inf","\U221E", levels(scri_input$age30_50),fixed=T)
      
      
      models_list <- list()
      report_list <- list()
      
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        for(istr in unique(scri_input$age30_50[cond_iae]) ){
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <-  istr  # "_age30_50" # istr
            
            global_name0 <- paste0( vax_priority, "_age30_50" )
            global_name  <- paste0( vax_priority, specif_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ age30_50:br:lab"
            
            if(nvax>=3) time_dep <- list( brand_3v_def )
            else        time_dep <- list( brand_2v_def )
            
            res <- scri_strata(  output_name  = output_name, 
                                 formula_text = formula_text,       time_seq = time_seq, 
                                 event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                 rws          = rws_def,
                                 time_dep     = time_dep,  
                                 #combine_vars =  c("sex","age4"), 
                                 start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                 data         = scri_input[scri_input$age30_50==istr & cond_iae,],
                                 image_plots = ae_event_first, image_brand=T, image_tit=istr,
                                 lab_orders   = lab_orders,
                                 paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                 lprint       = print_during_running,
                                 global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                 CI = CI_draw
            )
            ae_event_first <- F
            first_plot <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
      }  
      
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      ####
      #  save report_list and model_list
      save_results(global_name0, report_list, models_list)
      #
      #######################################################################################
      
      
      
      
      
      
      
      ########### per age30  #############
      #
      
      cat("\n\nAnalysis for ",icovid," per age30, brand.\n\n")
      
      if(all(names(scri_input) !="age30"))
        scri_input$age30 <- paste0("age",as.character(cut(scri_input$age_at_study_entry, c(-1,30,Inf))))
      #levels(scri_input$age30) <- gsub( "Inf","\U221E", levels(scri_input$age30),fixed=T)
      
      
      models_list <- list()
      report_list <- list()
      
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        for(istr in unique(scri_input$age30[cond_iae]) ){
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <- istr  #  "_age30" 
            
            global_name0 <- paste0( vax_priority, "_age30" )
            global_name  <- paste0( vax_priority, specif_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ age30:br:lab"
            
            if(nvax>=3) time_dep <- list( brand_3v_def )
            else        time_dep <- list( brand_2v_def )
            
            res <- scri_strata(  output_name  = output_name,  
                                 formula_text = formula_text,       time_seq = time_seq, 
                                 event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                 rws          = rws_def,
                                 time_dep     = time_dep,
                                 #combine_vars =  c("sex","age4"), 
                                 start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                 data         = scri_input[scri_input$age30==istr & cond_iae,],
                                 image_plots = ae_event_first, image_brand=T, image_tit=istr,
                                 lab_orders   = lab_orders,
                                 paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                 lprint       = print_during_running,
                                 global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                 CI = CI_draw
            )
            ae_event_first <- F
            first_plot <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
      }  
      
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      ####
      #  save report_list and model_list
      save_results(global_name0, report_list, models_list)
      #
      #######################################################################################
      
      
      
      ########### per sex  #############
      #
      
      cat("\n\nAnalysis for ",icovid," per sex, brand.\n\n")
      
      scri_input$sexc <- paste0("sex:",scri_input$sex )
      
      models_list <- list()
      report_list <- list()
      
      for(iae in ae_events){
        
        cat("\n",paste(iae, format(Sys.time()),"\n\n"))
        cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
        
        for(istr in unique(scri_input$sexc[cond_iae]) ){
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <- istr  # "_sex" # istr
            
            global_name0 <- paste0( vax_priority, "_sex" )
            global_name  <- paste0( vax_priority, specif_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ sexc:br:lab"
            
            if(nvax>=3) time_dep <- list( brand_3v_def )
            else        time_dep <- list( brand_2v_def )
            
            res <- scri_strata(  output_name  = output_name,
                                 formula_text = formula_text,       time_seq = time_seq, 
                                 event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                 rws          = rws_def,
                                 time_dep     = time_dep, 
                                 #combine_vars =  c("sex","age4"), 
                                 start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                 data         = scri_input[scri_input$sexc==istr & cond_iae,],
                                 image_plots = ae_event_first, image_brand=T, image_tit=istr,
                                 lab_orders   = lab_orders,
                                 paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                 lprint       = print_during_running,
                                 global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                 CI = CI_draw
            )
            ae_event_first <- F
            first_plot <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
      }  
      
      # plots:
      if(plot_during_running) 
        for(istr in names(report_list) )
          if(!is.null(report_list[[istr]][[1]]))
            plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
      gc()
      
      ####
      #  save report_list and model_list
      save_results(global_name0, report_list, models_list)
      #
      #######################################################################################
      
      
    }  # end of if(lmain)   
    
    
    #  days_to_week_infinity <- function( var, 
    #                                     infinity = T, 
    #                                     change_to_weeks = F, round_full_week = T, week_start=1, 
    #                                     add_before = "", add_after = "w" ,
    #                                     changes = list(list("-1]","<0"))
    #                                     ){
    #    browser()
    #    # change days to week
    #    if(change_to_weeks){
    #      if(round_full_week) var <- 7 * ( (var + week_start) %/% 7 )
    #      for(iw in c(-100:-1,100:1))
    #        var <- gsub( as.character(iw*7), paste0(iw,"w"), var, fixed=T)
    #    }
    #    
    #    if(infinity)
    #      var <- gsub( "Inf", "\U221E", var, fixed=T)
    #    var
    #    
    #  }
    #  days_to_week_infinity( levels(scri_input$dose12_diff_cat) )
    
    
    
    
    
    
    
    #########################################################################  
    
    #                      distance between doses:
    
    if(ldist){
    
      #########################################################################  
      ######################################################################### 
      #
      #            Preparation for distance between doses analyses:
      #
      
      
      scri_input$dose12_diff  <-  as.numeric(difftime(scri_input$date_vax2, scri_input$date_vax1 ,units="days"))
      
      if(F){
        # hist( scri_input$dose12_diff,seq(0,max( scri_input$dose12_diff,na.rm=T)+9,by=10))
        # plot( (scri_input$myopericarditis_days - scri_input$days_vax1), scri_input$dose12_diff)
        # plot( (scri_input$myopericarditis_days - scri_input$days_vax2), scri_input$dose12_diff)
        plot( (scri_input$myopericarditis_days - scri_input$days_vax2), scri_input$dose12_diff, ylim=c(0,100))
        abline(v=c(0,7,14,21,28,60))
        
        par(mfcol=c(2,3))
        for(ibr in brands){ 
          ylim <- c(0,max(scri_input$dose12_diff, na.rm=T))
          xlim <- range((scri_input$myopericarditis_days - scri_input$days_vax2), na.rm=T)
          
          xlim <- c(-120,150); ylim <- c(0,120)
          
          cond <- scri_input$myocarditis==1 & scri_input$type_vax2 == ibr
          plot( (scri_input$myopericarditis_days - scri_input$days_vax2)[cond], scri_input$dose12_diff[cond], ylim=ylim, xlim=xlim)
          abline(v=c(0,7,14,21,28,60)); title(paste("myocarditis",ibr))
          
          cond <- scri_input$pericarditis==1 & scri_input$type_vax2 == ibr
          plot( (scri_input$myopericarditis_days - scri_input$days_vax2)[cond], scri_input$dose12_diff[cond], ylim=ylim, xlim=xlim)
          abline(v=c(0,7,14,21,28,60)); title(paste("pericarditis",ibr))
        }
      }
      
      #   # define distance categories for number of days between vax1 and vax2 ('dist12'), between vax2 and vax3 ('dist32')
      #   distances <- list( list( dist12 = c(-Inf,-1,      10*7,  Inf),
      #                            dist23 = c(-Inf,-1,      25*7,  Inf),
      #                            name   = "_dist10w"                        ),
      #                      list( dist12 = c(-Inf,-1, 4*7,  8*7,  Inf),
      #                            dist23 = c(-Inf,-1,      25*7,  Inf),
      #                            name   = "_dist4_8w"                       )
      #                 )
      
      
      # loop for different categorizing of the distance between vaccines
      for(idist in 1:length(distances)){
        
        distance_12 <- distances[[idist]][["dist12"]]
        dist_name   <- distances[[idist]][["name"]]
        
        cat("\n\nAnalysis with distance between vaccines: ",dist_name ," \n\n")
        
        scri_input$dose12_diff_cat <- scri_input$dose12_diff
        scri_input$dose12_diff_cat[is.na(scri_input$dose12_diff_cat)] <- -999999999
        scri_input$dose12_diff_cat <- paste0("dist12:",cut( scri_input$dose12_diff_cat, distance_12 ))
        #scri_input$dose12_diff_cat <- paste0("dist12:",cut( scri_input$dose12_diff_cat, c(-Inf,-1, 10*7    ,Inf) ))   #10 weeks
        
        scri_input$dose12_diff_cat <- gsub( " Inf","Inf", scri_input$dose12_diff_cat,fixed=T)
        #scri_input$dose12_diff_cat <- gsub( "Inf","\U221E", scri_input$dose12_diff_cat,fixed=T)
        
        table1(scri_input$dose12_diff_cat)
        
        
        if(nvax>2){
          
          distance_23 <- distances[[idist]][["dist23"]]
          
          scri_input$dose23_diff  <-  as.numeric(difftime(scri_input$date_vax3, scri_input$date_vax2 ,units="days"))
          scri_input$dose23_diff_cat <- scri_input$dose23_diff
          scri_input$dose23_diff_cat[is.na(scri_input$dose23_diff_cat)] <- -999999999
          scri_input$dose23_diff_cat <- paste0("dist23:",cut( scri_input$dose23_diff_cat, distance_23 ))
          #scri_input$dose23_diff_cat <- paste0("dist23:",cut( scri_input$dose23_diff_cat, c(-Inf,-1,  25*7, Inf) ))
          
          scri_input$dose23_diff_cat <- gsub( " Inf","Inf", scri_input$dose23_diff_cat,fixed=T)
          #scri_input$dose23_diff_cat <- gsub( "Inf","\U221E", scri_input$dose23_diff_cat,fixed=T)
          table1(scri_input$dose23_diff_cat)
        }
        
        
        scri_input$vax01_distance <- paste0(  "d1:", "_dist:first_dose" )
        scri_input$vax12_distance <- paste0(  "d2:", "_", scri_input$dose12_diff_cat )
        scri_input$vax23_distance <- paste0(  "d3:", "_", scri_input$dose23_diff_cat )
        table1(scri_input$vax01_distance)
        table1(scri_input$vax12_distance)
        table1(scri_input$vax23_distance)
        
        
        
        ########### no strata  interaction with distance between vaccines:     #############
        #
        
        # for(iae in ae_events)
        #   brand_images( plot_data, ae_event=iae, tit="")
        
        cat("\n\nAnalysis for ",icovid," per distance between vaccines.\n\n")
        
        models_list <- list()
        report_list <- list()
        
        for(iae in ae_events){
          
          cat("\n",paste(iae, format(Sys.time()),"\n\n"))
          cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <- "_no_split" 
            
            global_name  <- paste0( vax_priority, specif_name, dist_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ lab:dist"
            
            if(nvax>=3) time_dep <- list( distance_3v_def )
            else        time_dep <- list( distance_2v_def )
            
            res <- scri_strata( output_name  = output_name, 
                                formula_text = formula_text,  time_seq = time_seq, 
                                event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                rws          = rws_def,
                                time_dep     = time_dep,              
                                start_obs    = "study_entry_days", end_obs = "study_exit_days",
                                data         = scri_input[cond_iae,],
                                image_plots = ae_event_first,
                                lab_orders = lab_orders,
                                paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                lprint = print_during_running,
                                global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                CI = CI_draw
            )
            first_plot <- F
            ae_event_first <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
        
        # plots:
        if(plot_during_running) 
          for(istr in names(report_list) )
            if(!is.null(report_list[[istr]][[1]]))
              plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
        gc()
        
        
        ###########
        #  save report_list and model_list
        save_results(global_name, report_list, models_list)
        #
        #######################################################################################
        
        
        
        
        ########### brand  interaction with distance between vaccines:     #############
        #
        
        # for(iae in ae_events)
        #   brand_images( plot_data, ae_event=iae, tit="")
        
        cat("\n\nAnalysis for ",icovid," per brand, distance between vaccines.\n\n")
        
        models_list <- list()
        report_list <- list()
        
        for(iae in ae_events){
          
          cat("\n",paste(iae, format(Sys.time()),"\n\n"))
          cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
          
          first_plot <- T
          ae_event_first <- T
          
          for(iii in ifelse(nvax>=3,6,3):1){
            
            
            if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
            if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
            if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
            
            if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
            if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
            if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
            
            specif_name  <- "_brand" 
            
            global_name  <- paste0( vax_priority, specif_name, dist_name )
            output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
            cat(paste(output_name, format(Sys.time()),"\n"))
            
            formula_text <-  "~ br_dist : lab "
            
            if(nvax>=3) time_dep <- list( brand_distance_3v_def )
            else        time_dep <- list( brand_distance_2v_def )
            
            res <- scri_strata(   
              output_name  = output_name, 
              formula_text = formula_text,       time_seq = time_seq, 
              event_time = paste0(iae,"_days"), event = iae, id="person_id",
              rws          = rws_def,
              time_dep     = time_dep,     
              start_obs    = "study_entry_days", end_obs = "study_exit_days",
              data         = scri_input[cond_iae,],
              image_plots = ae_event_first, image_brand=T, 
              lab_orders = lab_orders,
              paral = T, paral_vars = paral_vars, n_cores = n_cores,
              lprint = print_during_running,
              global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
              CI = CI_draw
            )
            first_plot <- F
            ae_event_first <- F
            
            report_list <- add_to_report_list(res$tabs,     output_name)
            models_list <- add_to_models_list(res$scri_all, output_name)
          }
        }
        
        # plots:
        if(plot_during_running) 
          for(istr in names(report_list) )
            if(!is.null(report_list[[istr]][[1]]))
              plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
        gc()
        
        ###########
        #  save report_list and model_list
        save_results(global_name, report_list, models_list)
        #
        #######################################################################################
        
        
        
        
        ########### per age30 per brand  interaction with distance between vaccines #############
        #
        
        cat("\n\nAnalysis for ",icovid," per age30, brand, distance between vaccines.\n\n")
        
        if(all(names(scri_input) !="age30"))
          scri_input$age30 <- paste0("age",as.character(cut(scri_input$age_at_study_entry, c(-1,30,Inf))))
        #levels(scri_input$age30) <- gsub( "Inf","\U221E", levels(scri_input$age30),fixed=T)
        
        
        models_list <- list()
        report_list <- list()
        
        for(iae in ae_events){
          
          cat("\n",paste(iae, format(Sys.time()),"\n\n"))
          cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
          
          for(istr in unique(scri_input$age30[cond_iae]) ){
            
            first_plot <- T
            ae_event_first <- T
            
            for(iii in ifelse(nvax>=3,6,3):1){
              
              
              if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
              if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
              if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
              
              if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
              if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
              if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
              
              specif_name  <- paste0(istr,"_brand")  #  "_age30" 
              
              global_name0 <- paste0( vax_priority, "_age30"   , dist_name )
              global_name  <- paste0( vax_priority, specif_name, dist_name )
              output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
              cat(paste(output_name, format(Sys.time()),"\n"))
              
              formula_text <-  "~ age30 : br_dist : lab"
              
              if(nvax>=3) time_dep <- list( brand_distance_3v_def )
              else        time_dep <- list( brand_distance_2v_def )
              
              res <- scri_strata(  output_name  = output_name,  
                                   formula_text = formula_text,       time_seq = time_seq, 
                                   event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                   rws          = rws_def,
                                   time_dep     = time_dep,
                                   #combine_vars =  c("sex","age4"), 
                                   start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                   data         = scri_input[scri_input$age30==istr & cond_iae,],
                                   image_plots = ae_event_first, image_brand=T, image_tit=istr,
                                   lab_orders   = lab_orders,
                                   paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                   lprint       = print_during_running,
                                   global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                   CI = CI_draw
              )
              ae_event_first <- F
              first_plot <- F
              
              report_list <- add_to_report_list(res$tabs,     output_name)
              models_list <- add_to_models_list(res$scri_all, output_name)
            }
          }
        }  
        
        # plots:
        if(plot_during_running) 
          for(istr in names(report_list) )
            if(!is.null(report_list[[istr]][[1]]))
              plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
        gc()
        
        ####
        #  save report_list and model_list
        save_results(global_name0, report_list, models_list)
        #
        #######################################################################################
        
        
        
        
        
        
        ########### per sex per brand  interaction with distance between vaccines #############
        #
        
        cat("\n\nAnalysis for ",icovid," per sex, brand, distance between vaccines.\n\n")
        
        if(all(names(scri_input) !="sexc"))
          scri_input$sexc <- paste0("sex:",scri_input$sex )
        
        
        models_list <- list()
        report_list <- list()
        
        for(iae in ae_events){
          
          cat("\n",paste(iae, format(Sys.time()),"\n\n"))
          cond_iae <- scri_input[,paste0("cond_covid_",substring(iae,1,7))]
          
          for(istr in unique(scri_input$sexc[cond_iae]) ){
            
            first_plot <- T
            ae_event_first <- T
            
            for(iii in ifelse(nvax>=3,6,3):1){
              
              
              if(iii==1) { rws_def <- rws_def_2vax_28; irw <- 28; idose <- 2 }
              if(iii==2) { rws_def <- rws_def_2vax_14; irw <- 14; idose <- 2 }
              if(iii==3) { rws_def <- rws_def_2vax_7;  irw <- 7 ; idose <- 2 }
              
              if(iii==4) { rws_def <- rws_def_3vax_28; irw <- 28; idose <- 3 }
              if(iii==5) { rws_def <- rws_def_3vax_14; irw <- 14; idose <- 3 }
              if(iii==6) { rws_def <- rws_def_3vax_7;  irw <- 7 ; idose <- 3 }
              
              specif_name  <- paste0(istr,"_brand")  #  "_sex" 
              
              global_name0 <- paste0( vax_priority, "_sex"     , dist_name )
              global_name  <- paste0( vax_priority, specif_name, dist_name )
              output_name  <- paste0( "_",substring(iae,1,7), global_name,"_",idose,"v","_",irw)
              cat(paste(output_name, format(Sys.time()),"\n"))
              
              formula_text <-  "~ sexc : br_dist : lab"
              
              if(nvax>=3) time_dep <- list( brand_distance_3v_def )
              else        time_dep <- list( brand_distance_2v_def )
              
              res <- scri_strata(  output_name  = output_name,  
                                   formula_text = formula_text,       time_seq = time_seq, 
                                   event_time = paste0(iae,"_days"), event = iae, id="person_id",
                                   rws          = rws_def,
                                   time_dep     = time_dep,
                                   #combine_vars =  c("sex","age4"), 
                                   start_obs    = "study_entry_days", end_obs = "study_exit_days", 
                                   data         = scri_input[scri_input$sexc==istr & cond_iae,],
                                   image_plots = ae_event_first, image_brand=T, image_tit=istr,
                                   lab_orders   = lab_orders,
                                   paral = T, paral_vars = paral_vars, n_cores = n_cores,
                                   lprint       = print_during_running,
                                   global_plot_name = paste0( substring(iae,1,7),global_name), add_global_plot = !first_plot,
                                   CI = CI_draw
              )
              ae_event_first <- F
              first_plot <- F
              
              report_list <- add_to_report_list(res$tabs,     output_name)
              models_list <- add_to_models_list(res$scri_all, output_name)
            }
          }
        }  
        gc()
        
        # plots:
        if(plot_during_running) 
          for(istr in names(report_list) )
            if(!is.null(report_list[[istr]][[1]]))
              plot_res(report_list[[istr]][[1]], main=paste(formula_text," + cal_time_cat"), col=col_list)
        
        ####
        #  save report_list and model_list
        save_results(global_name0, report_list, models_list)
        #
        #######################################################################################
        
      } # end of 'distance' loop
    }   # end of   if(ldist)    
      
  }  # end of 'covid' loop
  
  
  ##########
  # restore options:
  options(old_width)
  
  
}  # end of 'subpop' loop



