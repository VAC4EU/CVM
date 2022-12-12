# Program Information  ----------------------------------------------------
#
# Program:      step_12_2_run_scri.R 
# Author:       Svetlana Belitser
# Description:  calls functions which runs an SCRI on specified datasets 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival" 
#               input:   data_vax_SCRI 
#               output:  g_output/scri/*  
#
#               parameters: in 07_scri_inputs.R 
#  
#               function: p_macro/scri_tools.R                    
#               function: p_macro/table1.R  
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
  
  # scri output_directory for export:  
  sdr0 <- paste0(direxpsubpop[[subpop]], "scri/")
  dir.create(file.path(sdr0), showWarnings = FALSE, recursive = TRUE)
  
  # scri output_directory for models not for export:
  sdr_models0 <- paste0(diroutput, "scri/")
  if(length(subpopulations_non_empty)>1) sdr_models0 <- paste0(diroutput, "scri/",subpop,"/")
  dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
  
  
  cat(paste0('\n\t"',subpop,'":\n\n'))
  
  
  memory.limit(10^13)
  #detach("package:data.table", unload = TRUE)
  
  
  # Import Data -------------------------------------------------------------
  
  # Load dataset 'data_vax' by loading file 'data_vax_SCRI.RData'
  load(paste0(dirtemp, "data_vax_SCRI", suffix[[subpop]], ".RData"))
  
  if("pat_n" %in% names(data_vax)) id <- "pat_n"
  
  gc()
  
  ########################################
  #
  # create a list with calendar time interval:
  #
  time_seq <- vector("list",length=length(time_interval_width))
  for(i in 1:length(time_interval_width))
    time_seq[[i]] <- seq(min(data_vax[,"study_entry_days"],na.rm=T)-time_interval_starts[i],max(data_vax[,"study_exit_days"],na.rm=T)+time_interval_width[i]-1,by=time_interval_width[i])
  names(time_seq) <- paste0("period",time_interval_width,"d_start_",-time_interval_starts,"d")
  
  ########################################
  
  ##########################################################################################
  #
  #   these parameters are used for all analyses:
  #
  extra_options <- list( paral                = lparal,    # if T ==> library(parallel) is started in function 'scri'
                         n_cores              = n_cores,   # Don't define it! Define it only if it doesn't work on its own! number of cores/threads to use.
                         time_seq             = time_seq, 
                         print_during_running = F,
                         lprint               = F,
                         plot_during_running  = F, 
                         leventplot           = leventplot, 
                         max_n_points         = max_n_points,  # ?1000
                         lplot                = lplot,
                         CI_draw              = CI_draw,
                         lforest              = lforest,
                         lplots               = T,
                         col                  = col_list,
                         performance          = T
  )
  
  # default in function 'define_rws': (vax2 takes precedence over vax1) the risk window of dose 2 takes precedence over the risk window of dose 1
  
  for(iae in ae_events){
    
    print(iae)
    
    if(!(paste0(iae,"_days") %in% names(data_vax))) { cat(paste0('\nevent "',iae,'" not found.\n\n')); next }

    if(lmain){
      
      print("main part")
      
      # SCCS output_directory for the event:  EXPORT 
      sdr <- paste0(sdr0, iae,"/")
      dir.create(sdr, showWarnings = FALSE, recursive = TRUE)
      
      # SCCS output_directory for the event:  LOCAL 
      sdr_models <- paste0(sdr_models0, iae,"/")
      dir.create(sdr_models, showWarnings = FALSE, recursive = TRUE)
      
      ##########################################################################################
      #
      #   add 'sdr_tabs' and 'sdr_models' parameters
      #
      extra_options$sdr_tabs   <- sdr              
      extra_options$sdr_models <- sdr_models 
      
      
      
      
      ###########################################################################################
      ###################################  event  ############################################### 
      #   
      event_info <- list( event=iae, event_time =paste0(iae,"_days"), event_date =paste0(iae,"_date") )
      
      
      ###########################################################################################
      #
      #             vax_name="vax_number": dose1, dose2, dose3, ...
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      
      ###################################################
      #  baseline tables
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
      #characteristics(data=data_vax, event=iae, path_file_name=paste0(sdr,"baseline.txt"), vax_name="vax_number", condition_value="", age="age_at_study_entry", lab_orders=vax_def$lab_orders )
      
      ###########  vax_number & no split  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, event_info=event_info,  extra_parameters = extra_options, add_to_itself=F, lplots=F  )        
      
      ## cut_points_name="7d" { [-91;-30], [-29;-1], [0;0], [1;7], [8;14], [15;21], [22;28] } 
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, event_info=event_info,  extra_parameters = extra_options, add_to_itself=T, lplots=T   )        
      
      
      ###########  vax_number & brand ( no distance):  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=F)
      
      ## cut_points_name="7d"
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options, add_to_itself=T, lplots=T )
      
      
      
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
          
          data_vax_strata <- data_vax[data_vax[,strata_var]==strata_value,]
          
          ###########  vax_number & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options, add_to_itself=F, lplots=F  )        
          
          ## cut_points_name="7d" { [-91;-30], [-29;-1], [0;0], [1;7], [8;14], [15;21], [22;28] } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options, add_to_itself=T, lplots=T   )        
          
          
          ###########  vax_number & brand ( no distance) per stratum  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          # in models with calendar time intervals use only events from the current stratum:
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, use_all_events=F,
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
        
        # SCCS output_subdirectory 'covid' in 'event' directory 
        sdr_covid <- paste0(sdr, "covid","/")
        dir.create(sdr_covid, showWarnings = FALSE, recursive = TRUE)
        
        # SCCS output_subdirectory 'covid' in the event directory: 
        sdr_covid_models <- paste0(sdr_models, "covid","/")
        dir.create(sdr_covid_models, showWarnings = FALSE, recursive = TRUE)
        
        #########################
        #   copy extra_options  and 
        #   change  'sdr_tabs' and 'sdr_models' parameters    
        #
        extra_options_covid            <- extra_options
        extra_options_covid$sdr_tabs   <- sdr_covid
        extra_options_covid$sdr_models <- sdr_covid_models
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
          #
          #             vax_name="vax_number": dose1, dose2, dose3, ...
          #
          vax_def0 <- scri_data_parameters(data=data_vax, vax_name="vax_number", vax_time="vax_days", vax_date="vax_date" )
          
          ###################################################
          #  baseline tables
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
          #characteristics(data=data_vax, event=iae, path_file_name=paste0(sdr_covid,"baseline_",icovid,".txt"), vax_name="vax_number", condition_value=icovid, age="age_at_study_entry", lab_orders=vax_def$lab_orders )
          
          ###########  vax_number & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options_covid, add_to_itself=F, lplots=F  ) 
          ## cut_points_name="7d"
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info,  extra_parameters = extra_options_covid, add_to_itself=T, lplots=T )        
          
          
          ###########  vax_number & brand ( no distance):  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options_covid, add_to_itself=F, lplots=F )
          ## cut_points_name="7d"
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, strata_var="covid_selection_name", strata_value=covid_value, use_all_events=F,
                       event_info=event_info, extra_parameters = extra_options_covid, add_to_itself=T, lplots=T )
          
          
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
      
      #########################
      #   copy extra_options  and 
      #   change  'sdr_tabs' and 'sdr_models' parameters    
      #
      extra_options_dist            <- extra_options
      extra_options_dist$sdr_tabs   <- sdr_dist
      extra_options_dist$sdr_models <- sdr_dist_models
      extra_options_dist$time_seq   <- NULL
      
      
      ###########################################################################################
      ###################################  event  ############################################### 
      #   
      event_info <- list( event=iae, event_time =paste0(iae,"_days"), event_date =paste0(iae,"_date") )
      
      
      ###########################################################################################
      #
      #             vax_name="vax_number": dose1, dose2, dose3, ...
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      
      ###################################################
      #  baseline tables
      #vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, data=data_vax, vax_dep = c( after="dist_gt_60" ) )
      ##characteristics(data=data_vax, event=iae, path_file_name=paste0(sdr_dist,"baseline_vax_number.txt"), vax_name="vax_number", condition_value=vax_def0$data_parameters$vax_name, age="age_at_study_entry", lab_orders=vax_def$lab_orders )
      
      ###########  vax_number & dist  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61] } 
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d",   no_last_interval_after=T, 
                             data=data_vax, vax_dep = c( after="dist_gt_60" ) )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F  )        
      
      ###########  vax_number & brand wih distance:  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d"
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
      
      extra_options_dist$extra_name <- ""
      
      ###########################################################################################
      #
      #             vax_name="vax_name": dose1.1, dose1.2, boost1, boost2, ...
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_name",         vax_time = "vax_days",        vax_date     = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      ###################################################
      #  baseline tables
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T,  data=data_vax )
      #characteristics(data=data_vax, event=iae, path_file_name=paste0(sdr_dist,"baseline_vax_name.txt"), vax_name="vax_name", condition_value=vax_def0$data_parameters$vax_name, age="age_at_study_entry", lab_orders=vax_def$lab_orders )
      
      ###########  vax_name & no split  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 }
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T,  data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
      
      ## cut_points_name="7d" : 
      vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d",  no_last_interval_after=T,  data=data_vax )
      res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax, event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T  )        
      
      ###########  vax_name & brand ( wihout distance ):  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d"
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short" ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
      
      
      ###########  vax_name & brand wih distance:  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
      
      ## cut_points_name="7d"
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
      res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
      
      extra_options_dist$extra_name <- ""
      
      gc()
      
      ###########################################################################################
      #
      #             for combination of brands (historical)
      #
      #
      vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                        id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
      extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
      
      ###########  vax_number & combination of previous and current brand :  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_with_prev"  ))
      # without formula:
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T ) )
      if(class(res)[[1]]== "try-error") forest_plots_tab(res[[1]][[1]][[1]])      
      
      ###########  vax_number & brand history sorted, or ignore order, i.e. you don't know which one was the first, which one was the second, ...  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_history_sorted"  ))
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T ) )
      
      
      ###########  vax_number & brand history  ##### 
      # 
      ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
      vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                            data=data_vax, vax_dep = c( before="type_history"  ))
      res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax, event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=T ) )
      
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
          
          data_vax_strata <- data_vax[data_vax[,strata_var]==strata_value,]
          
          ###########################################################################################
          ###################################  event  ############################################### 
          #   
          event_info <- list( event=iae, event_time =paste0(iae,"_days"), event_date =paste0(iae,"_date") )
          
          
          ###########################################################################################
          #
          #             vax_name="vax_number": dose1, dose2, dose3, ...
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_number & dist  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61] } 
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d",   no_last_interval_after=T, 
                                 data=data_vax, vax_dep = c( after="dist_gt_60" ) )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
          ## cut_points_name="7d"
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d",   no_last_interval_after=T, 
                                 data=data_vax, vax_dep = c( after="dist_gt_60" ) )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T  )        
          
          ###########  vax_number & brand wih distance:  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value, 
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
          
          extra_options_dist$extra_name <- ""
          
          
          ###########################################################################################
          #
          #             vax_name="vax_name": dose1.1, dose1.2, boost1, boost2, ...
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_name",         vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_name & no split  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181], >181 }
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d",  no_last_interval_after=T,  data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=F, lplots=F  )        
          
          ## cut_points_name="7d"
          vax_def  <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="^7d",  no_last_interval_after=T,  data=data_vax )
          res <- scri( formula = "~ lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info,  extra_parameters = extra_options_dist, add_to_itself=T, lplots=T  )        
          
          ###########  vax_name & brand ( wihout distance ):  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61], [62;181]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62,182), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short" ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
          
          
          ###########  vax_name & brand wih distance:  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;0], [1;28], [28;61]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,29,62), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=F )
          
          ## cut_points_name="7d"
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,1,8,15,22,29), cut_points_name="7d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="vax_brand_short", after="dist_gt_60"  ))
          res <- scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                       event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=T, lplots=T )
          
          extra_options_dist$extra_name <- ""
          
          gc()
          
          ###########################################################################################
          #
          #             for combination of brands (historical)
          #
          #
          vax_def0 <- scri_data_parameters( data =  data_vax, vax_name  = "vax_number",       vax_time = "vax_days",        vax_date     = "vax_date", 
                                            id   = "pat_n",   start_obs = "study_entry_days", end_obs  = "study_exit_days", censored_vars = "death_days" )
          extra_options_dist$extra_name <- vax_def0$data_parameters$vax_name
          
          ###########  vax_number & combination of previous and current brand :  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_with_prev"  ))
          # without formula:
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T ) )
          if(class(res)[[1]]== "try-error") forest_plots_tab(res[[1]][[1]][[1]])      
          
          ###########  vax_number & brand history sorted, or ignore order, i.e. you don't know which one was the first, which one was the second, ...  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_history_sorted"  ))
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options_dist, add_to_itself=F, lplots=T ) )
          
          
          ###########  vax_number & brand history  ##### 
          # 
          ## cut_points_name="28d" :  { [-91;-30], [-29;-1], [0;28]  } 
          vax_def <- define_rws(vax_def0,  cut_points_before = c(-90,-29,0), cut_points_after = c(0,29), cut_points_name="28d", no_last_interval_after=T, 
                                data=data_vax, vax_dep = c( before="type_history"  ))
          res <- try( scri( formula = "~ brand:lab", vax_def = vax_def, data = data_vax_strata, strata_value=strata_value,
                            event_info=event_info, extra_parameters = extra_options, add_to_itself=F, lplots=T ) )
          
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

