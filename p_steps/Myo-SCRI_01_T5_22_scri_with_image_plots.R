# Program Information  ----------------------------------------------------
#
# Program:      Myo-SCRI-01-T5_20_.R # in the 'steps' directory 
# Author:       Svetlana Belitser
# Description:  calls functions to run SCRI analyses on a specified dataset 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival" 
#               input:   g_intermediate/data_vax_SCRI.RData
#               output:  g_output/scri/*  
#
#               parameters: in 08_SCRI_parameters.R 
#  
#               function: p_macro/scri_tools.R                    
#
#

# Housekeeping  -----------------------------------------------------------
# install and load packages

if(!any(ls()=="thisdir"))   thisdir   <- getwd()
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"/g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"/g_output/")

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp,   "scri")),            showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,   "/log_files/scri")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "scri")),            showWarnings = FALSE, recursive = TRUE)


for (subpop in subpopulations_non_empty) {
  
  # scri output directory for export:  
  sdr0 <- paste0(direxpsubpop[[subpop]], "scri_with_plots/")
  dir.create(file.path(sdr0), showWarnings = FALSE, recursive = TRUE)
  
  sdr0_flowchart <- paste0(direxpsubpop[[subpop]], "scri_with_plots/flowchart_for_vax/")
  dir.create(file.path(sdr0_flowchart), showWarnings = FALSE, recursive = TRUE)
  
  
  # scri output_directory for models NOT for export:
  sdr_models0 <- paste0(diroutput, "scri_with_plots/")
  if(length(subpopulations_non_empty)>1) sdr_models0 <- paste0(diroutput, "scri_with_plots/",subpop,"/")
  dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
  
  
  cat(paste0('\n\t"',subpop,'":\n\n'))
  
  
  # copy 'g_export/.../scri_small/' into "g_export/.../scri_with_plots":
  if(dir.exists(paste0(direxpsubpop[[subpop]], "scri_small/")))
    copyDirectory(from=paste0(direxpsubpop[[subpop]], "scri_small/"), to=sdr0, recursive=TRUE)
  
  # copy 'g_output/.../scri_small/' into "g_output/.../scri_with_plots":
  if(length(subpopulations_non_empty)>1){  
    if(dir.exists(paste0(diroutput, "scri_small/",subpop,"/")))
      copyDirectory(from=paste0(diroutput, "scri_small/",subpop,"/"), to=sdr_models0, recursive=TRUE)
  }
  else {
    if(dir.exists(paste0(diroutput, "scri_small/")))
      copyDirectory(from=paste0(diroutput, "scri_small/"), to=sdr_models0, recursive=TRUE)
  }
  
  #detach("package:data.table", unload = TRUE)
  
  
  # Import Data -------------------------------------------------------------
  
  # Load dataset 'data_vax' by loading file 'data_vax_SCRI.RData'
  load(paste0(dirtemp, "scri/data_vax_SCRI", suffix[[subpop]], ".RData"))
  data_vax <- as.data.frame(data_vax)
  
  if("pat_n" %in% names(data_vax)) id <- "pat_n"
  if(any(tolower(names(data_vax))=="datasource")) 
    if(length(unique(data_vax[,tolower(names(data_vax))=="datasource"])))
      thisdatasource <- data_vax[1,tolower(names(data_vax))=="datasource"]
  if(!any(ls()=="thisdatasource")) thisdatasource <- ""
  
  gc()
  
  #############   SCRI models ############################
  #
  #
  old_width = options(width=300)
  
  ########################################
  #
  # create a list with calendar time interval:
  #
  time_seq <- vector("list",length=length(time_interval_width))
  for(i in 1:length(time_interval_width))
    time_seq[[i]] <- seq(min(data_vax[,"study_entry_days"],na.rm=T)-time_interval_starts[i],max(data_vax[,"study_exit_days"],na.rm=T)+time_interval_width[i]-1,by=time_interval_width[i])
  names(time_seq) <- paste0("period",time_interval_width,"d_start_minimum_minus_",-time_interval_starts,"d_or_",sapply(time_seq,function(x)ifelse(x[1]<0,"minus_","")),sapply(time_seq,function(x)abs(x[1])),"d")
  
  days_1jan2020 <- as.numeric(difftime(as.Date("2020-08-31"),as.Date("2020-01-01"),units="days"))
  
  ########################################
  
  ##########################################################################################
  #
  #   these parameters are used for all analyses:
  #
  extra_options <- list( lparal               = lparal,    # if T ==> library(parallel) is started in function 'scri'
                         n_cores              = n_cores,   # Don't define it! Define it only if it doesn't work on its own! number of cores/threads to use.
                         time_seq             = time_seq, 
                         time_seq_ref         = time_seq_ref, 
                         data_source          = thisdatasource,
                         print_during_running = F,
                         lprint               = F,
                         plot_during_running  = F, 
                         leventplot           = leventplot,
                         max_n_points         = max_n_points,  # ?1000
                         lplot                = lplot,
                         CI_draw              = CI_draw,
                         lplot_hist           = F,
                         lforest              = lforest,
                         #lplots               = T,
                         col                  = col_list,
                         performance          = T
  )
  lmain      <- T
  scri_small <- F
  lcovid     <- F
  ldist      <- T
  
  # default in function 'define_rws': (vax2 takes precedence over vax1) the risk window of dose 2 takes precedence over the risk window of dose 1
  
  
  
  ###################################################
  #  baseline tables without events info
  if(!file.exists(paste0(sdr0_flowchart, thisdatasource,"_characteristics_vax_number.RData")))
      characteristics(data=data_vax, vax_name="vax_number", vax_part=T, event_vax_part=F, path=sdr0_flowchart, id="pat_n", condition_value="", age="age_at_study_entry", lparal=lparal_flowchart, n_cores=n_cores_flowchart )
  
  
  for(iae in ae_events){
    
    cat(paste("\n",iae,":\n\n"))
    
    
    if(!(paste0(iae,"_days") %in% names(data_vax))) { cat(paste0('\nevent "',iae,'" not found.\n\n')); next }
    
    data_vax[,paste0(iae,"_death_days")] <- pmin(data_vax[,paste0(iae,"_days")],data_vax$death_days,na.rm=T)
    data_vax[,paste0(iae,"_death_date")] <- pmin(data_vax[,paste0(iae,"_date")],data_vax$death_date,na.rm=T)
    data_vax[,paste0(iae,"_death")]      <- as.numeric(!is.na(data_vax[,paste0(iae,"_death_days")]))
    
    if(lmain){
      
      print("main part")
      
      # SCCS output_directory for the event:  EXPORT 
      sdr <- paste0(sdr0, iae,"/")
      dir.create(sdr, showWarnings = FALSE, recursive = TRUE)
      
      # SCCS output_directory for the event:  LOCAL 
      sdr_models <- paste0(sdr_models0, iae,"/")
      dir.create(sdr_models, showWarnings = FALSE, recursive = TRUE)
      
      # SCCS output_directory for the event:  EXPORT 
      sdr_cohort <- paste0(sdr,"cohort/")
      dir.create(sdr_cohort, showWarnings = FALSE, recursive = TRUE)
      
      ##########################################################################################
      #
      #   add 'sdr_tabs' and 'sdr_models' parameters
      #
      extra_options$sdr_tabs    <- sdr              
      extra_options$sdr_models  <- sdr_models 
      extra_options$path_cohort <- sdr_cohort
      
      
      ###########################################################################################
      ###################################  event  ############################################### 
      #   
      event_info        <- list( event=iae,                  event_time =paste0(iae,      "_days"), event_date =paste0(iae,      "_date") )
      
      ###########################################################################################
      #
      #             vax_name="vax_number": dose1, dose2, dose3, ...
      #
      vax_def0        <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax1="dose 1",  vax_time = "vax_days",        vax_date     = "vax_date", 
                                               id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      
      ###################################################
      #  baseline tables for event and 'vax_number':
      if(!file.exists(paste0(sdr, thisdatasource,"_characteristics_",iae,"_vax_number.RData")))
        characteristics(data=data_vax, vax_name="vax_number", vax_part=F, event_vax_part=T, event=iae, path=sdr, id="pat_n", condition_value="", age="age_at_study_entry" ,
                        lparal=lparal_flowchart, n_cores=n_cores_flowchart )
      
      data_vax$inclusion <- !is.na(data_vax[,event_info$event_time]) & !is.na(data_vax$vax_days_v1) &   
        ( data_vax[,vax_def0$data_parameters$start_obs] <= pmax(data_vax$vax_days_v1-90-365,days_1jan2020) | (data_vax$age_at_study_entry<=2 & data_vax[,vax_def0$data_parameters$start_obs]<=data_vax$vax_days_v1-90) ) & 
        ( data_vax$vax_days_v1-90 <= data_vax[,event_info$event_time]  | data_vax$vax_days_v1-90 <= data_vax$death_days  & !is.na(data_vax$death_days)  )
     
    
      data_vax_strata             <- data_vax[data_vax$inclusion, ]
      
      ##############  vax_number & no split  ######################
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info,  extra_parameters = extra_options, add_to_itself=F, lplots=F  )        
      
      ## cut_points_name="7d" { [-91;-30], [-29;-1], [0;0], [1;7], [8;14], [15;21], [22;28] } ; with 'çohort' plots
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info,  extra_parameters = extra_options, add_to_itself=T, lplots=T   )       
      
 
      #############  vax_number & brand ( no distance):  ####################### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=F)
      
      ## cut_points_name="7d" ; with 'çohort' plots
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options, add_to_itself=T, lplots=T )
      
       #######################################################################
      #
      #       strata analyse with brand for: age30, age30_50, sex, sex_age30 
      #
      
      # strata variables:
      for( strata_var in c( "age30","age30_50", "sexc", "sex_age30") ){
        
        print(strata_var)
        
        # values of the current strata variable
        strata_values <- unique(data_vax[,strata_var]); strata_values <- strata_values[!is.na(strata_values)]
        
        # delete strata "age(-1,30]" if variable 'age30_50'  because it is also in 'age30'
        if(strata_var=="age30_50") strata_values <- strata_values[strata_values!="age(-1,30]"]
        
        for(strata_value in strata_values){ 
          
          data_vax$strata_cond        <- data_vax[,strata_var]==strata_value
          data_vax_strata             <- data_vax[data_vax$strata_cond & data_vax$inclusion            ,]  
          data_vax_strata_also_deaths <- data_vax[data_vax$strata_cond & data_vax$inclusion_also_deaths,]  
          #if(nrow(data_vax_strata)==0) next
          
          ###########  vax_number & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options, add_to_itself=F, lplots=F  )        
          
          ## cut_points_name="7d" { [-91;-30], [-29;-1], [0;0], [1;7], [8;14], [15;21], [22;28] } ;  with 'cohort' plots
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options, add_to_itself=T, lplots=T   )        
          
           ###########  vax_number & brand ( no distance) per stratum  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"; with 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options, add_to_itself=T, lplots=T )
          
         } # end for strata_value
        gc()
      } # end for strata_var
    } # end lmain
    
    
    ################################################################################
    #
    #            for covid subsets
    #
    if(lcovid){
      if(any(tolower(names(data_vax)) != "covid19_date") ){
        
        print("covid")
        
        if(!lmain)
          data_vax$inclusion <- !is.na(data_vax[,event_info$event_time]) & !is.na(data_vax$vax_days_v1) &   
            ( data_vax[,vax_def0$data_parameters$start_obs] <= pmax(data_vax$vax_days_v1-90-365,days_1jan2020) | (data_vax$age_at_study_entry<=2 & data_vax[,vax_def0$data_parameters$start_obs]<=data_vax$vax_days_v1-90) ) & 
            ( data_vax$vax_days_v1-90 <= data_vax[,event_info$event_time]  | data_vax$vax_days_v1-90 <= data_vax$death_days  & !is.na(data_vax$death_days)  )
        
        # SCCS output_subdirectory 'covid' in 'event' directory 
        sdr_covid <- paste0(sdr0, iae, "/covid/")
        dir.create(sdr_covid, showWarnings = FALSE, recursive = TRUE)
        
        # SCCS output_subdirectory 'covid' in the event directory: 
        sdr_covid_models <- paste0(sdr_models0, iae, "/covid/")
        dir.create(sdr_covid_models, showWarnings = FALSE, recursive = TRUE)
        
        # SCCS output_directory for the event:  LOCAL 
        sdr_covid_cohort <- paste0(sdr,"cohort/")
        dir.create(sdr_covid_cohort, showWarnings = FALSE, recursive = TRUE)
        
        #########################
        #   copy extra_options  and 
        #   change  'sdr_tabs' and 'sdr_models' parameters    
        #
        extra_options_covid             <- extra_options
        extra_options_covid$sdr_tabs    <- sdr_covid
        extra_options_covid$sdr_models  <- sdr_covid_models
        extra_options_covid$path_cohort <- sdr_covid_cohort
        ##############
        
        # all data:
        # cat("\n\nALL DATA: table for variables:\n")
        # print(substring(iae,1,7))
        # table1( unique(data_vax[ ,c(id, iae)]) [ , iae]  )
        
        
        for(icovid in c( "covidplus30d","nocovid")){
          
          ############# create covid variables: ##############.
          data_vax$covid <- as.numeric(!is.na(data_vax[,"covid19_date"]))
          if(icovid == "covidplus30d"){
            #####
            # select only those who did not have COVID before (the event date + 30 days). These 30 days added extra to make sure that 
            # if someone had COVID then there was no connection to the event.
            # create covid selection variables: no_covid_before_myocard_30d, no_covid_before_pericar_30d, no_covid_before_myoperi_30d, ...
            data_vax[,"no_covid_before_event_plus30d"] <- ( (as.Date(data_vax[,paste0(iae,"_date")]) + 30) < as.Date(data_vax$covid19_date) )  &  !is.na(data_vax$covid19_date) & !is.na(data_vax[,paste0(iae,"_date")])
            data_vax[is.na(data_vax$covid19_date) | is.na(data_vax[,paste0(iae,"_date")]),"no_covid_before_event_plus30d"] <- T
            
            data_vax[                                          , "covid_selection_name"] <- ""
            data_vax[data_vax[,"no_covid_before_event_plus30d"], "covid_selection_name"] <- paste0("no_covid_before_",iae,"_plus30d")
            covid_value <- paste0("no_covid_before_",iae,"_plus30d")
            
            #cat(paste0("\n\nno_covid_before ",iae," plus30d: table for variables:\n"))
            #print(c("covid","no_covid_before_event_plus30d"))
            #print(table1( unique( data_vax[ ,c(id, "covid","covid_selection_name",iae)] )[ ,c(iae,"covid","covid_selection_name")] ))
          }
          
          if(icovid =="nocovid"){
            #####
            # no covid in [start_of_observation; stop_of_observation]:
            data_vax[                                , "covid_selection_name"] <- ""
            data_vax[is.na(data_vax[,"covid19_date"]), "covid_selection_name"] <- "no_covid_observed"
            covid_value <- "no_covid_observed"
            
            #cat(paste0("\n\nno covid in [start_of_observation; stop_of_observation]:\n"))
            #print(c("covid","no_covid_observed"))
            #print(table1( unique( data_vax[ ,c(id, "covid","covid_selection_name",iae)] )[ ,c(iae,"covid","covid_selection_name")] ))
          }
          
          
          
          ###########################################################################################
          ###################################  event  ############################################### 
          #   
          event_info        <- list( event=iae,                  event_time =paste0(iae,      "_days"), event_date =paste0(iae,      "_date") )
          
          ###########################################################################################
          #
          #             vax_name="vax_number": dose1, dose2, dose3, ...
          #
          vax_def0        <- scri_data_parameters( data = data_vax, vax_name  = "vax_number",       vax1="dose 1",  vax_time = "vax_days",        vax_date     = "vax_date",
                                                   id   = "pat_n",  start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          
          ###################################################
          #  baseline tables without event info (only one time to run, i.e., for the first event).
          if(iae==ae_events[1]){
            if(file.exists(paste0(sdr0_flowchart, "characteristics_vax_number_",icovid,".RData"))) file.rename(from=paste0(sdr0_flowchart, "characteristics_vax_number_",icovid,".RData"), to=paste0(sdr0_flowchart, thisdatasource,"_characteristics_vax_number_",icovid,".RData")) 
            if(file.exists(paste0(sdr0_flowchart, "characteristics_vax_number_",icovid,".txt"  ))) file.rename(from=paste0(sdr0_flowchart, "characteristics_vax_number_",icovid,".txt"  ), to=paste0(sdr0_flowchart, thisdatasource,"_characteristics_vax_number_",icovid,".txt"  )) 

            if(!file.exists(paste0(sdr0_flowchart, thisdatasource,"_characteristics_vax_number_",icovid,".RData")))
              characteristics(data=data_vax[data_vax$covid_selection_name!="",], vax_name="vax_number", vax_part=T, event_vax_part=F, path=sdr0_flowchart, id="pat_n", condition_value=icovid, age="age_at_study_entry", lparal=lparal_flowchart, n_cores=n_cores_flowchart )
          }
          ###################################################
          #  baseline tables for event and 'vax_number':
   
          if(!file.exists(paste0(sdr_covid, thisdatasource,"_characteristics_",iae,"_vax_number_",icovid,".RData")))
            characteristics(data=data_vax[data_vax$covid_selection_name!="",], vax_name="vax_number", vax_part=F, event_vax_part=T, event=iae, path=sdr_covid, subpopulations=c(1,2), id="pat_n", condition_value=icovid, age="age_at_study_entry" ,
                            lparal=lparal_flowchart, n_cores=n_cores_flowchart )
          
          ###########  vax_number & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options_covid, add_to_itself=F, lplots=F  ) 
          ## cut_points_name="7d"; no 'cohort' plots
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options_covid, add_to_itself=T, lplots=T, leventplot=F )        
          
          
          ###########  vax_number & brand ( no distance):  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options_covid, add_to_itself=F, lplots=F )
          ## cut_points_name="7d"; no 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options_covid, add_to_itself=T, lplots=T, leventplot=F )
          
          
        }   # end  covid 'for'
      }# end codid variable
      gc()
    }  # end lcovid
    
    
    
    
    
    
    
    #########################################################################
    #
    #  this part:  without time adjustment
    #
    if(ldist){
      
      print("additional")
      
      if(!lmain)      
        data_vax$inclusion <-  !is.na(data_vax[,event_info$event_time])  & !is.na(data_vax$vax_days_v1) &   
          ( data_vax[,vax_def0$data_parameters$start_obs] <= pmax(data_vax$vax_days_v1-90-365,days_1jan2020) | (data_vax$age_at_study_entry<=2 & data_vax[,vax_def0$data_parameters$start_obs]<=data_vax$vax_days_v1-90) ) & 
          ( data_vax$vax_days_v1-90 <= data_vax[,event_info$event_time]  | data_vax$vax_days_v1-90 <= data_vax$death_days  & !is.na(data_vax$death_days)  )
      
      ##########################################
      #
      #   create additional variables:
      #
      if( !all( c("type_with_prev","type_history","type_history_sorted") %in% names(data_vax)) ){
        
        # create variable with tow brands: from the previous and current doses :
        data_vax$type_with_prev <- format(data_vax[,"vax_brand_short"])
        cond_prev_exists <- c( F, data_vax[-1,id] == data_vax[-nrow(data_vax),id] ) 
        cond_next_exists <- c( data_vax[-nrow(data_vax),id] == data_vax[-1,id], F ) 
        if(any(cond_prev_exists)) 
          data_vax$type_with_prev[ cond_prev_exists ] <- paste0( data_vax$type_with_prev[cond_next_exists],"-", format(data_vax[,"vax_brand_short"])[cond_prev_exists] )
        
        # create variable with history of brands:
        data_vax$type_history <- format(data_vax[,"vax_brand_short"])
        prev_steps <- 1
        while(T){
          cond_prev_exists <- c( F, data_vax[-1,id] == data_vax[-nrow(data_vax),id] ) & data_vax$vax_n==prev_steps+1
          cond_next_exists <- c( data_vax[-nrow(data_vax),id] == data_vax[-1,id], F )  & data_vax$vax_n==prev_steps  
          if(!any(cond_prev_exists)) break
          data_vax$type_history[ cond_prev_exists ] <- paste0( data_vax$type_history[cond_next_exists],"-", format(data_vax[,"vax_brand_short"])[cond_prev_exists] )
          prev_steps <- prev_steps + 1
        }
        # table1(data_vax[,,"type_history"])
        
        # create variable with sorted history of brands:
        data_vax$type_history_sorted <- unlist(lapply( strsplit(data_vax$type_history, "[ |-]+"), function(x)paste0(sort(x), collapse ="-") ))
        # table1(data_vax[,c("type_history_sorted","type_history")])
      }
      #
      ##################################
      
      
      # SCCS output_subdirectory 'distance_combi' in 'event' directory 
      sdr_dist <- paste0(sdr0, iae, "/distance_combi/")
      dir.create(sdr_dist, showWarnings = FALSE, recursive = TRUE)
      
      # SCCS output_subdirectory 'distance_combi' in the event directory: 
      sdr_dist_models <- paste0(sdr_models0, iae,"/distance_combi/")
      dir.create(sdr_dist_models, showWarnings = FALSE, recursive = TRUE)
      
      # SCCS output_subdirectory 'distance_combi' in the event directory: 
      sdr_dist_cohort <- paste0(sdr_dist,"cohort/")
      dir.create(sdr_dist_cohort, showWarnings = FALSE, recursive = TRUE)
      
      #########################
      #   copy extra_options  and 
      #   change  'sdr_tabs' and 'sdr_models' parameters    
      #
      extra_options_dist             <- extra_options
      extra_options_dist$sdr_tabs    <- sdr_dist
      extra_options_dist$sdr_models  <- sdr_dist_models
      extra_options_dist$path_cohort <- sdr_dist_cohort
      extra_options_dist$time_seq    <- NULL
      
      data_vax_strata <- data_vax[data_vax$inclusion,]
      
      ###########################################################################################
      ###################################  event  ############################################### 
      #   
      event_info        <- list( event=iae,                  event_time =paste0(iae,      "_days"), event_date =paste0(iae,      "_date") )
      event_info_deaths <- list( event=paste0(iae,"_death"), event_time =paste0(iae,"_death_days"), event_date =paste0(iae,"_death_date") )
      
      ###########################################################################################
      #
      #             vax_name="vax_number": dose1, dose2, dose3, ...
      #
      vax_def0        <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax1="dose 1",  vax_time = "vax_days",        vax_date     = "vax_date", 
                                               id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      vax_def0_deaths <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax1="dose 1",  vax_time = "vax_days",        vax_date     = "vax_date", 
                                               id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days"  )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      
      ###########  vax_number & dist  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61] } ; with 'cohort' plots
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d",   no_last_interval_after=T, 
                             data=data_vax, vax_dep = c( after="dist_gt_60" ) )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=T  )        
      
      ###########  vax_number & brand with distance:  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d"; with 'cohort' plots 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
      
      extra_options_dist$extra_name <- ""
      
      ###########################################################################################
      #
      #             vax_name="vax_name": dose1.1, dose1.2, boost1, boost2, ...
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_name",         vax1="dose 1.1",  vax_time = "vax_days",        vax_date      = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      
      ###################################################
      #  baseline tables without event info
      if(iae==ae_events[1])
        characteristics(data=data_vax, vax_name="vax_name",   vax_part=T, event_vax_part=F, path=sdr0_flowchart, id="pat_n", condition_value=vax_def0$data_parameters$vax_name, age="age_at_study_entry", lparal=lparal_flowchart, n_cores=n_cores_flowchart )

      ###################################################
      #  baseline tables for event and 'vax_number':
      characteristics(data=data_vax, vax_name="vax_name", vax_part=F, event_vax_part=T, event=iae, path=sdr_dist, subpopulations=c(1,2), id="pat_n", condition_value=vax_def0$data_parameters$vax_name, age="age_at_study_entry", 
                      lparal=lparal_flowchart, n_cores=n_cores_flowchart )
    
      ###########  vax_name & no split  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 }
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T,  data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
      
      ## cut_points_name="7d" ; with 'cohort' plots
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d",  no_last_interval_after=T,  data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T  )        
      
      ###########  vax_name & brand ( without distance ):  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d" ; with 'cohort' plots
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
      
      
      ###########  vax_name & brand wih distance:  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d"; no 'cohort'  plots
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F  )
      
      extra_options_dist$extra_name <- ""
      
      gc()
      
      ###########################################################################################
      #
      #             for combination of brands (historical)
      #
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",  vax1="dose 1",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      
      ###########  vax_number & combination of previous and current brand :  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_with_prev"  ))
      # without formula; no 'cohort' plots:
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
      if(class(res)[[1]]== "try-error") forest_plots_tab(res[[1]][[1]][[1]])      
      
      ###########  vax_number & brand history sorted, or ignore order, i.e. you don't know which one was the first, which one was the second, ...  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } ; no 'cohort' plots:
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_history_sorted"  ))
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
      
      
      ###########  vax_number & brand history  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } ; no 'cohort' plots:
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_history"  ))
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
      
      extra_options_dist$extra_name <- ""
      
      
      
      
      ########################################################
      #
      #    analyse for "age30","age30_50","sexc","sex_age30"
      #
      
      # strata variables:
      for( strata_var in c( "age30","age30_50", "sexc", "sex_age30") ){
        
        # SCCS output_subdirectory 'distance_combi' in 'event' directory 
        sdr_dist_stratum <- paste0(sdr0, iae, "/distance_combi/",ifelse(strata_var%in%c("age30","age30_50"), "age", strata_var), "/" )
        dir.create(sdr_dist_stratum, showWarnings = FALSE, recursive = TRUE)
        
        # SCCS output_subdirectory 'distance_combi' in the event directory: 
        sdr_dist_stratum_models <- paste0(sdr_models0, iae,"/distance_combi/",ifelse(strata_var%in%c("age30","age30_50"), "age", strata_var), "/" )
        dir.create(sdr_dist_stratum_models, showWarnings = FALSE, recursive = TRUE)
        
        #########################
        #   copy extra_options  and 
        #   change  'sdr_tabs' and 'sdr_models' parameters    
        #
        extra_options_dist$sdr_tabs   <- sdr_dist_stratum
        extra_options_dist$sdr_models <- sdr_dist_stratum_models
        
        
        # values of the current strata variable
        strata_values <- unique(data_vax[,strata_var]); strata_values <- strata_values[!is.na(strata_values)]
        
        # delete strata "age(-1,30]" if variable 'age30_50'  because it is also in 'age30'
        if(strata_var=="age30_50") strata_values <- strata_values[strata_values!="age(-1,30]"]
        
        for(strata_value in strata_values){ 
          
          data_vax$strata_cond <- data_vax[,strata_var]==strata_value
          data_vax_strata <- data_vax[data_vax$strata_cond & data_vax$inclusion,]
          
          ###########################################################################################
          ###################################  event  ############################################### 
          #   
          event_info <- list( event=iae, event_time =paste0(iae,"_days"), event_date =paste0(iae,"_date") )
          
          
          ###########################################################################################
          #
          #             vax_name="vax_number": dose1, dose2, dose3, ...
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",     vax1="dose 1",    vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_number & dist  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61] } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d",   no_last_interval_after=T, 
                                 data=data_vax, vax_dep = c( after="dist_gt_60" ) )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
          ## cut_points_name="7d"; no 'cohort' plots
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d",   no_last_interval_after=T, 
                                 data=data_vax, vax_dep = c( after="dist_gt_60" ) )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value, 
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F  )        
          
          ###########  vax_number & brand wih distance:  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"; no 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value, 
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F )
          
          extra_options_dist$extra_name <- ""
          
          
          ###########################################################################################
          #
          #             vax_name="vax_name": dose1.1, dose1.2, boost1, boost2, ...
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_name",     vax1="dose 1.1",      vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_name & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 }
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T,  data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
          
          ## cut_points_name="7d"; no 'cohort' plots:
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d",  no_last_interval_after=T,  data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F  )        
          
          ###########  vax_name & brand ( wihout distance ):  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"; no 'cohort' plots:
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F )
          
          
          ###########  vax_name & brand wih distance:  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"; no 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T, leventplot=F )
          
          extra_options_dist$extra_name <- ""
          
          gc()
          
          ###########################################################################################
          #
          #             for combination of brands (historical)
          #
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",   vax1="dose 1",      vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_number & combination of previous and current brand :  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_with_prev"  ))
          # without formula; no 'cohort' plots:
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
          if(class(res)[[1]]== "try-error") forest_plots_tab(res[[1]][[1]][[1]])      
          
          ###########  vax_number & brand history sorted, or ignore order, i.e. you don't know which one was the first, which one was the second, ...  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } ; no 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_history_sorted"  ))
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
          
          
          ###########  vax_number & brand history  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } ; no 'cohort' plots
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_history"  ))
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, inclusion="inclusion", strata_var="strata_cond", strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T, leventplot=F ) )
          
          extra_options_dist$extra_name <- ""
          
          gc()
          
        } # end for strata_value
      } # end for strata_var
      
    } # end of ldist
    
  } # end iae      
  ##########
  # restore options:
  options(old_width)
  
  
}  # end of 'subpop' loop

