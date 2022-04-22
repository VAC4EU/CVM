# Program Information  ----------------------------------------------------

# Program:      step_12_1_create_scri_dataset.R 
# Author:       Svetlana Belitser; Anna Schultze, Ema Alsina, Sophie Bots, Ivonne Martens 
# Description:  calls a function which creates SCRI dataset in wide format 
#               takes user provided brands and outcomes
#               outputs datasets provided N is >5. 
#               if N is < 5 no output data is provided and a warning printed to console 
# Requirements: 
#               dependencies: 07_scri_inputs.R
#               input:  D3_study_variables_for_SCRI (g_intermediate) in RData format 
#               parameters: in 07_scri_inputs.R
#               output: g_output/scri/scri_input_[brand]_[outcome]_flowchart.txt
#                       g_intermediate/scri/scri_input_[brand]_[outcome]_input.csv               

# Housekeeping  -----------------------------------------------------------
# load functions:
#for(ifunc in func_names)
#  source(paste0(func_dir,ifunc))

if(!any(ls()=="thisdir"))   thisdir   <- getwd()
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"/g_intermediate/")
#if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"/g_output/")


# ensure required folders are created 
for (subpop in subpopulations_non_empty) {
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE, direxp, direxpsubpop[[subpop]])
  
  dir.create(file.path(paste0(dirtemp,    "scri")),            showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(paste0(thisdirexp, "scri")),            showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(paste0(thisdirexp, "scri/flowcharts")), showWarnings = FALSE, recursive = TRUE)
  
  # Import Data -------------------------------------------------------------
  load(paste0(dirtemp, raw_data_name, suffix[[subpop]], ".RData"))
  temp_name<-get(paste0(raw_data_name, suffix[[subpop]]))
  rm(list=paste0(raw_data_name, suffix[[subpop]]))
  assign(intermediate_data, as.data.frame(temp_name))
  rm(temp_name)
  
  #load(paste0(dirtemp, raw_data))
  #scri_data_extract <- eval(parse(text = raw_data_name))
  
  
  
  if(F){   # for Meritxell
    date_vars <- c("date_vax1", "date_vax2","date_vax_1", "date_vax_2","covid19_date", "covid_date", "covid_19_date", 
                   "study_entry_date", "study_exit_date", "myocarditis_date", "pericarditis_date", "myopericarditis_date", "death_date" ) 
    for(date_var in date_vars)
      if( any(names(scri_data_extract)==date_var) )
        if(mode(scri_data_extract[,date_var])=="numeric")  scri_data_extract[,date_var] = as.Date(scri_data_extract[,date_var], origin="1970-01-01")
  }
  
  
  
  
  # get max number of vaccine doses:
  nvax <- names(scri_data_extract)[ substring(names(scri_data_extract),1,8)=="date_vax" ]
  nvax <- max(as.numeric(substring(nvax,9)))
  
  # change vax names to standard names: ( "type_vax_2"  ==>   "type_vax2" )
  for(iv in 1:nvax){
    names(scri_data_extract)[names(scri_data_extract)==paste0("type_vax_",iv)] <- paste0("type_vax",iv)
    scri_data_extract[,paste0("vax",iv)] <- as.numeric(!is.na(scri_data_extract[,paste0("date_vax",iv)]))
  }
  #table1(scri_data_extract[,paste0("type_vax",1:3)])
  
  # covid data:
  names(scri_data_extract)[names(scri_data_extract) %in% c("covid_19_date","covid_date")] <- "covid19_date"
  scri_data_extract$covid <- as.numeric(!is.na(scri_data_extract$covid19_date))
  
  if(F){  # for ancient or silly computers and VDI    ==> first delete all no_vax1 or no_myoperi 
    dim(scri_data_extract)
    scri_data_extract <- scri_data_extract[ !is.na(scri_data_extract$myocarditis_date) |  !is.na(scri_data_extract$pericarditis_date), ]
    dim(scri_data_extract)
    scri_data_extract <- scri_data_extract[ !is.na(scri_data_extract$date_vax1), ]
    dim(scri_data_extract)
  }
  
  # create variable 'myopericarditis_date'  
  print("Freq table for non-missing(T) values:")
  scri_data_extract$myopericarditis_date <- pmin(scri_data_extract$myocarditis_date,scri_data_extract$pericarditis_date,na.rm=T)
  print(table1( paste0(  "myoperi:", !is.na(scri_data_extract$myopericarditis_date),
                         " myo:",!is.na(scri_data_extract$myocarditis_date),
                         " peri:", !is.na(scri_data_extract$pericarditis_date) )
  ))
  #print(table( myo_perio=paste( !is.na(scri_data_extract$myocarditis_date),!is.na(scri_data_extract$pericarditis_date)), myoperi= !is.na(scri_data_extract$myopericarditis_date)))
  
  # create 'myocarditis', 'pericarditis' and 'myopericarditis'  variables
  scri_data_extract$myocarditis     <- as.numeric(!is.na(scri_data_extract$myocarditis_date ))
  scri_data_extract$pericarditis    <- as.numeric(!is.na(scri_data_extract$pericarditis_date))
  scri_data_extract$myopericarditis <- as.numeric(!is.na(scri_data_extract$myopericarditis_date))
  
  # distribution of vaccine dates between vax1 and vax2. The minimum should not be small !!!!
  scri_data_extract$dose_diff = as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1,units="days"))
  summary(scri_data_extract$dose_diff)
  if(any(min(scri_data_extract$dose_diff[!is.na(scri_data_extract$dose_dif)])<14)) 
    warning(paste("There are ", sum(scri_data_extract$dose_diff<14, na.rm=T), " persons with (date_vax2 - date_vax1) < 14 days."))
  
  
  # change time scale to days since study_entry_date (or 01-09-2020)  
  start_scri_date <- as.Date("2020-09-01")
  scri_data_extract$start_study_date <- start_scri_date
  scri_data_extract$start_study_date <- pmax( scri_data_extract$study_entry_date, scri_data_extract$start_study_date)
  # age at start_study study_entry_date (or 01-09-2020)
  #scri_data_extract$age_at_study_entry <- scri_data_extract$age_at_study_entry
  #scri_data_extract$age_at_study_entry <- as.numeric(difftime(scri_data_extract$start_study_date,scri_data_extract$birth_date,units="days")) / 365.25 
  print("age at study entry date:")
  summary(scri_data_extract$age_at_study_entry)
  
  #if( any(min(scri_data_extract$age_at_study_entry,na.rm=T) < as.numeric(difftime(scri_data_extract$start_study_date,Sys.Date(),units="days")) / 365.25) ) 
  #  warning("check the range of variable 'age_at_study_entry'")
  
  
  #####  create variables with days from 'start_study_date'
  # calculate #days from 'start_study_date'  for the following variables:
  date_vars <- c("study_entry_date", "study_exit_date", paste0("date_vax",1:nvax),  # "date_vax1","date_vax2", # "date_vax3"
                 "myocarditis_date", "pericarditis_date", "myopericarditis_date","death_date", "covid19_date" )  #, "birth_date")
  if(any( !(date_vars %in% names(scri_data_extract)) ))  # check if there are such colname in scri_data_extract
    warning( paste0("Variable[s] '", paste0(date_vars[!(date_vars %in% names(scri_data_extract))],collapse=","),"' is not in the dataset."))
  date_vars <- date_vars[date_vars %in% names(scri_data_extract)]   # founded column names
  date_vars_days <- gsub("date","days",date_vars )                  # new var-names with 'days' instead of 'date'
  for(icol in 1:length(date_vars)){
    scri_data_extract[,date_vars_days[icol]] <- round(as.numeric( difftime( scri_data_extract[,date_vars[icol]], start_scri_date, units="days") ))  
    #scri_data_extract[,date_vars_days[icol]] <- round(as.numeric( difftime( scri_data_extract[,date_vars[icol]], scri_data_extract$start_study_date, units="days") ))  
    gc()
  }
  # check that myoperi after the start of observations or study_entry_date (???)
  #if( any(scri_data_extract$myopericarditis_days<0))
  #  warning(paste("There are ", sum(scri_data_extract$myopericarditis_days<0, na.rm=T), " persons with  'myopericarditis_date' before 'study_entry_date' ."))
  
  # check that myoperi after the start of observations or study_entry_date (???)
  if( any( tmp <-(as.numeric(difftime(scri_data_extract$myopericarditis_date,as.Date("2020-09-01"),units="days"))>=0 & 
                  scri_data_extract$myopericarditis_days<0),na.rm=T ) )
    warning(paste("There are ", sum(tmp, na.rm=T), " persons with  'myopericarditis_date' after '31-08-2020' and before 'study_entry_date' ."))
   
  
  
  
  # check!!!
  if(T){   
    
    ######################################################
    #  check the 'study_entry_days' and 'study_exit_days'
    #
    if(any( scri_data_extract$study_entry_days > scri_data_extract$days_vax1 & !is.na(scri_data_extract$days_vax1) & !is.na(scri_data_extract$study_entry_days), na.rm=T )){
      warning(paste("'study_entry_days' after the vax1 for ", sum( scri_data_extract$start > scri_data_extract$days_vax1, na.rm=T), "rows. They are deleted!"))
      nrow0<-nrow(scri_data_extract)  
      scri_data_extract <- scri_data_extract[ scri_data_extract$study_entry_days <= scri_data_extract$days_vax1 & scri_data_extract$vax1==1 & !is.na(scri_data_extract$study_entry_days),]  # if study_entry_days > vax1 ==>  delete rows
      print(c(new_nrow=nrow(scri_data_extract), old_nrow=nrow0)) 
    }
    
    if(any( scri_data_extract$study_exit_days < scri_data_extract$days_vax1, na.rm=T )){
      warning(paste0("There are ",sum( scri_data_extract$study_exit_days < scri_data_extract$days_vax1,na.rm=T)," persons with 'study_exit_time' before vax1. They are deleted!"))
      nrow0<-nrow(scri_data_extract)  
      scri_data_extract <- scri_data_extract[ scri_data_extract$days_vax1<=scri_data_extract$study_exit_days & scri_data_extract$vax1==1 & !is.na(scri_data_extract$study_entry_days),]  # if study_exit_days < vax1 ==>  delete rows
      print(c(new_nrow=nrow(scri_data_extract), old_nrow=nrow0)) 
    }  
    
    if(any( (cond <- !is.na(scri_data_extract$days_vax2) & scri_data_extract$study_exit_days < scri_data_extract$days_vax2 & !is.na(scri_data_extract$dstudy_exit_days) ) )){
      warning(paste0("There are ",sum( cond, na.rm=T)," persons with 'study_exit_time' before vax2 ==> study_exit_days <- days_vax2 !"))
      scri_data_extract$study_exit_days[cond]  <- scri_data_extract$days_vax2[cond] 
    } 
    
    
    #delete vax2 if vax2_days - vax1_days < 5
    if(nvax>2){
      while( any( cond <- !is.na(scri_data_extract$days_vax2) & (scri_data_extract$days_vax2 - scri_data_extract$days_vax1) < 5) ){
        warning(paste(sum(cond),"'dose2' were replace with next dose because dose2 is less than 5 days after dose1."))
        for(iv in 3:nvax){
          scri_data_extract[cond, paste0("date_vax",iv-1) ] <- scri_data_extract[cond, paste0("date_vax",iv) ]
          scri_data_extract[cond, paste0("days_vax",iv-1) ] <- scri_data_extract[cond, paste0("days_vax",iv) ]
          scri_data_extract[cond, paste0("type_vax",iv-1) ] <- scri_data_extract[cond, paste0("type_vax",iv) ]
          scri_data_extract[cond, paste0("vax"     ,iv-1) ] <- scri_data_extract[cond, paste0("vax"     ,iv) ]
        } 
      }  
    }
    scri_data_extract$dose_diff  <-  as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1 ,units="days"))
    
    # start of obs >29 days before the vax 1:
    nrow(scri_data_extract)
    scri_data_extract <- scri_data_extract[scri_data_extract$days_vax1- scri_data_extract$study_entry_days >29,]
    nrow(scri_data_extract)
    
    # table of numbers of vaccine doses:
    table1(scri_data_extract[,paste0("vax",1:nvax)])
    scri_data_extract[is.na(scri_data_extract$date_vax2) & !is.na(scri_data_extract$date_vax3),]
    
    #
    # nrow(scri_data_extract)
    # scri_data_extract <- scri_data_extract[  ! ( is.na(scri_data_extract$date_vax2) & !is.na(scri_data_extract$date_vax3)  ),]
    # nrow(scri_data_extract)
  
  }
  
  
  
  
  
  
  
  ############
  #
  # Conditions:
  #
  # condition
  scri_data_extract$cond_gender        <-  !is.na(scri_data_extract$sex)                  # non-missing gender
  scri_data_extract$cond_age           <-  !is.na(scri_data_extract$age_at_study_entry)      # non-missing age
  scri_data_extract$cond_vax1          <-  !is.na(scri_data_extract$date_vax1)            # non-missing vax1 
  scri_data_extract$cond_type_vax1     <-  !is.na(scri_data_extract$type_vax1)            # non-missing type_vax1
  scri_data_extract$cond_myoperi       <-  !is.na(scri_data_extract$myopericarditis_date) # non-missing myoperi
  
  scri_data_extract$cond_vax1_myoperi  <-  scri_data_extract$cond_vax1 & scri_data_extract$cond_myoperi  # non=missing vax1 and myoperi
  
  # condtitons for risk windows:
  prevax_per_start <- -90
  buffer_per_start <- -29
  
  # condition for myoperi before 01-09-2020:
  scri_data_extract$cond_myoperi_from_sept2020    <-  !is.na(scri_data_extract$myopericarditis_date)  & 
    scri_data_extract$myopericarditis_date >= as.Date("2020-09-01") 
  
  # condition for myoperi before vax1-90days:
  scri_data_extract$cond_from_control_window <- scri_data_extract$cond_vax1_myoperi & 
    scri_data_extract$myopericarditis_days >= scri_data_extract$days_vax1 + prevax_per_start 
  
  # condition myoperi during buffer interval
  scri_data_extract$cond_not_buffer <-   scri_data_extract$cond_vax1_myoperi & 
    !( scri_data_extract$days_vax1 + buffer_per_start <=  scri_data_extract$myopericarditis_days  & 
         scri_data_extract$myopericarditis_days       <   scri_data_extract$days_vax1 )
  
  ######
  # condition for between risk windows of vax1 and vax2
  vax1_end         <- 28
  vax2_end         <- 28
  
  # condition for myoperi after last risk window (of vax1 or vax2):
  end_last_rw <- pmax(scri_data_extract$days_vax1 + vax1_end, scri_data_extract$days_vax2 + vax2_end, na.rm=T )
  end_last_rw <- pmin(end_last_rw, scri_data_extract$study_exit_days )
  scri_data_extract$cond_till_last_rw_28d <-  scri_data_extract$myopericarditis_days <= end_last_rw &  
    !is.na(scri_data_extract$myopericarditis_days) & !is.na(scri_data_extract$date_vax1)
  
  scri_data_extract$cond_not_betw_rw1_rw2_28d <-  scri_data_extract$cond_vax1_myoperi & !is.na(scri_data_extract$date_vax2)     & 
    ( scri_data_extract$myopericarditis_days <=  scri_data_extract$days_vax1 + vax1_end   |  
        scri_data_extract$days_vax2            <=  scri_data_extract$myopericarditis_days     ) 
  scri_data_extract$cond_not_betw_rw1_rw2_28d <- scri_data_extract$cond_not_betw_rw1_rw2_28d | is.na(scri_data_extract$date_vax2)
  
  #####
  # condition for between risk windows of vax1 and vax2
  vax1_end <- 60
  vax2_end <- 60
  
  # condition for myoperi not after last risk window (of vax1 or vax2):
  end_last_rw <- pmax(scri_data_extract$days_vax1 + vax1_end, scri_data_extract$days_vax2 + vax2_end, na.rm=T )
  end_last_rw <- pmin(end_last_rw, scri_data_extract$study_exit_days )
  scri_data_extract$cond_till_last_rw_60d <-  scri_data_extract$myopericarditis_days <= end_last_rw &  
    !is.na(scri_data_extract$myopericarditis_days) & !is.na(scri_data_extract$date_vax1)
  
  scri_data_extract$cond_not_betw_rw1_rw2_60d <-  scri_data_extract$cond_vax1_myoperi & !is.na(scri_data_extract$date_vax2)    & 
         ( scri_data_extract$myopericarditis_days <= scri_data_extract$days_vax1 + vax1_end   | 
             scri_data_extract$days_vax2          <= scri_data_extract$myopericarditis_days     )
  scri_data_extract$cond_not_betw_rw1_rw2_60d <- scri_data_extract$cond_not_betw_rw1_rw2_60d | is.na(scri_data_extract$date_vax2)
  
  
  
  #######
  ###  create  inclusion variables 
  
  scri_data_extract$include <- 
    scri_data_extract$cond_gender                   &     
    scri_data_extract$cond_age                      &       
    scri_data_extract$cond_vax1                     &       
    scri_data_extract$cond_type_vax1                &  
    scri_data_extract$cond_myoperi                  &   
    scri_data_extract$cond_myoperi_from_sept2020    &     
    scri_data_extract$cond_from_control_window  
  
  scri_data_extract$end_study_28d <- 
    scri_data_extract$include                       & 
    scri_data_extract$cond_till_last_rw_28d            
  
  scri_data_extract$study_period_28d <- 
    scri_data_extract$include                       & 
    scri_data_extract$cond_not_betw_rw1_rw2_28d     & 
    scri_data_extract$cond_till_last_rw_28d       
  
  
  scri_data_extract$study_period_60d <- 
    scri_data_extract$include                       & 
    scri_data_extract$cond_not_betw_rw1_rw2_60d     &  
    scri_data_extract$cond_till_last_rw_60d       
  
  
  ############
  #
  # flowchart small:
  #
  
  old_width = options(width=200)
  
  # create constant 'dap'
  dap <- ifelse( any(names(scri_data_extract)=="DAP"), scri_data_extract$DAP[1], "")
  # create vector brands with available brands, and "all"
  brands <- c("all", sort(unique(scri_data_extract$type_vax1)))
  flowchart_list     <- list()
  for(ibr in brands){
    # select brand
    cond_brand <- rep(T,nrow(scri_data_extract))
    if(ibr != "all") cond_brand <- scri_data_extract$type_vax1==ibr & !is.na(scri_data_extract$type_vax1)
    
    flowchart_simple <- with(scri_data_extract[cond_brand,],
      c(      all                     = length( cond       <- rep(T,length(sex))                        )      
            , gender                  = sum(    cond       <- cond      &  cond_gender                  )      
            , age                     = sum(    cond       <- cond      &  cond_age                     )      
            , vax1                    = sum(    cond       <- cond      &  cond_vax1                    )      
            , type_vax1               = sum(    cond       <- cond      &  cond_type_vax1               ) 
            , myoperi                 = sum(    cond       <- cond      &  cond_myoperi                 )  
            , myoperi_from_sep2020    = sum(    cond       <- cond      &  cond_myoperi_from_sept2020   ) 
            , myoperi_from_rw1        = sum(  ( cond_incl  <- cond      &  cond_from_control_window )   ) 
            , after_last28            = sum(    cond       <- cond_incl &  cond_till_last_rw_28d        ) 
            , between_rw12_28         = sum(    cond       <- cond      &  cond_not_betw_rw1_rw2_28d    )
            , without_buffer_28       = sum(    cond       <- cond      &  cond_not_buffer              )
            , after_last60            = sum(    cond       <- cond_incl &  cond_till_last_rw_60d        ) 
            , between_rw12_60         = sum(    cond       <- cond      &  cond_not_betw_rw1_rw2_60d    ) 
            , without_buffer_60       = sum(    cond       <- cond      &  cond_not_buffer              )
      )) 
    
    flowchart_simple <- cbind.data.frame( name=names(flowchart_simple), n=flowchart_simple )
    rownames(flowchart_simple) <- NULL
  
    flowchart_simple$label <- ""
    flowchart_simple$label[ flowchart_simple$name == "all"                    ] <- "total"                              
    flowchart_simple$label[ flowchart_simple$name == "gender"                 ] <- "sex not missing"                             
    flowchart_simple$label[ flowchart_simple$name == "age"                    ] <- "age not missing"                             
    flowchart_simple$label[ flowchart_simple$name == "vax1"                   ] <- "received at least one vaccine dose of any brand"                             
    flowchart_simple$label[ flowchart_simple$name == "type_vax1"              ] <- "vaccine brand not missing"                             
    flowchart_simple$label[ flowchart_simple$name == "myoperi"                ] <- "at least one myopericarditis diagnosis ever"                             
    flowchart_simple$label[ flowchart_simple$name == "myoperi_from_sep2020"   ] <- "at least one myopericarditis diagnosis in time period of interest (starts 1sep2020)"                             
    flowchart_simple$label[ flowchart_simple$name == "myoperi_from_rw1"       ] <- "at least one myopericarditis diagnosis within study period (starts at control window)"                             
    flowchart_simple$label[ flowchart_simple$name == "after_last28"           ] <- "at least one myopericarditis diagnosis between start control and end last 28-day risk window"                                   
    flowchart_simple$label[ flowchart_simple$name == "between_rw12_28"        ] <- "excluding myopericarditis diagnoses occuring between 28-day risk window 1 and 2"                             
    flowchart_simple$label[ flowchart_simple$name == "without_buffer_28"      ] <- "excluding myopericarditis diagnoses occuring during the buffer period (28 day risk windows)"                             
    flowchart_simple$label[ flowchart_simple$name == "after_last60"           ] <- "at least one myopericarditis diagnosis between start control and end last 60-day risk window"                                   
    flowchart_simple$label[ flowchart_simple$name == "between_rw12_60"        ] <- "excluding myopericarditis diagnoses occuring between 60-day risk window 1 and 2"                               
    flowchart_simple$label[ flowchart_simple$name == "without_buffer_60"      ] <- "excluding myopericarditis diagnoses occuring during the buffer period (60 day risk windows)"    
    
    flowchart_simple         <- flowchart_simple[,c(3,1,2)]
    flowchart_simple$percent <- 100 * flowchart_simple$n /  flowchart_simple$n[flowchart_simple$name=="all"]
    #print(format(flowchart_simple[,c("label","n","percent")],  justify="left", digits=3))
    
    
    flowchart_simple_interactions <-
        with(scri_data_extract[cond_brand,],
         table1(
            paste(
               c("no_sex",               "   sex"               ) [ 1 + cond_gender                ]
              ,c("no_age",               "   age"               ) [ 1 + cond_age                   ]     
              ,c("no_vax1",              "   vax1"              ) [ 1 + cond_vax1                  ]                         
              ,c("no_vax1_type",         "   vax1_type"         ) [ 1 + cond_type_vax1             ]
              ,c("no_myoperi",           "   myoperi"           ) [ 1 + cond_myoperi               ] 
              ,c("before_sep2020",       "  from_sep2020"       ) [ 1 + cond_myoperi_from_sept2020 ]
              ,c("before_vax1-90",       "  from_vax1-90"       ) [ 1 + cond_from_control_window   ]
              ,c("after_last_rw28d",     " till_last_rw28d"     ) [ 1 + cond_till_last_rw_28d      ]  
              ,c("    betw_rw1_rw2_28d", "not_betw_rw1_rw2_28d" ) [ 1 + cond_not_betw_rw1_rw2_28d  ]
              ,c("after_last_rw60d",     " till_last_rw60d"     ) [ 1 + cond_till_last_rw_60d      ]
              ,c("    betw_rw1_rw2_60d", "not_betw_rw1_rw2_60d" ) [ 1 + cond_not_betw_rw1_rw2_60d  ]
              ,c("   buffer",            "no_buffer"            ) [ 1 + cond_not_buffer            ]
      )))
    
    # 28 days
    flowchart_simple_interactions_28d <-
        with(scri_data_extract[cond_brand,],
         table1(
            paste(
                 c("no_sex",               "   sex"               ) [ 1 + cond_gender                ]
                ,c("no_age",               "   age"               ) [ 1 + cond_age                   ]     
                ,c("no_vax1",              "   vax1"              ) [ 1 + cond_vax1                  ]                         
                ,c("no_vax1_type",         "   vax1_type"         ) [ 1 + cond_type_vax1             ]
                ,c("no_myoperi",           "   myoperi"           ) [ 1 + cond_myoperi               ] 
                ,c("before_sep2020",       "  from_sep2020"       ) [ 1 + cond_myoperi_from_sept2020 ]
                ,c("before_vax1-90",       "  from_vax1-90"       ) [ 1 + cond_from_control_window   ]
                ,c("after_last_rw28d",     " till_last_rw28d"     ) [ 1 + cond_till_last_rw_28d      ]  
                ,c("    betw_rw1_rw2_28d", "not_betw_rw1_rw2_28d" ) [ 1 + cond_not_betw_rw1_rw2_28d  ]
                ,c("   buffer",            "no_buffer"            ) [ 1 + cond_not_buffer            ]
    )))
    # 60 days
    flowchart_simple_interactions_60d <-
        with(scri_data_extract[cond_brand,],
         table1(
            paste(
               c("no_sex",               "   sex"               ) [ 1 + cond_gender                ]
              ,c("no_age",               "   age"               ) [ 1 + cond_age                   ]     
              ,c("no_vax1",              "   vax1"              ) [ 1 + cond_vax1                  ]                         
              ,c("no_vax1_type",         "   vax1_type"         ) [ 1 + cond_type_vax1             ]
              ,c("no_myoperi",           "   myoperi"           ) [ 1 + cond_myoperi               ] 
              ,c("before_sep2020",       "  from_sep2020"       ) [ 1 + cond_myoperi_from_sept2020 ]
              ,c("before_vax1-90",       "  from_vax1-90"       ) [ 1 + cond_from_control_window   ]
              ,c("after_last_rw60d",     " till_last_rw60d"     ) [ 1 + cond_till_last_rw_60d      ]
              ,c("    betw_rw1_rw2_60d", "not_betw_rw1_rw2_60d" ) [ 1 + cond_not_betw_rw1_rw2_60d  ]
              ,c("   buffer",            "no_buffer"            ) [ 1 + cond_not_buffer            ]
      )))
    
    if(ibr=="all"){ 
      sink(paste0(diroutput, 'scri/flowcharts/flowchart_',dap,'.txt'),append = F)
      cat(paste0("\n\nDAP: ", dap, "  Brands: \n\n"))
      print(table1( scri_data_extract[,c("type_vax1","type_vax2")], title="" ))
    }
    else     sink(paste0(diroutput, 'scri/flowcharts/flowchart_',dap,'.txt'),append = T)
    
      #old_width = options(width=200)
    
      cat(paste0("\n\nDAP: ", dap, "  Brand: ", ibr, "\n\n"))
      print(format(flowchart_simple[,c("label","n","percent")],  justify="left", digits=3))
      
      if(ibr=="all"){ 
        cat(paste0("\n\nDAP: ", dap, "  Brands for included persons: \n\n"))
        print(table1( scri_data_extract[scri_data_extract$include,c("type_vax1","type_vax2")], title="" ))
      }
      
      cat(paste0('\n\nDAP: ', dap, '  Brand: ', ibr, '   "Interaction"-Table:\n\n'))
      print(  flowchart_simple_interactions )
      
      # rw is 28 days
      cat(paste0('\n\nDAP: ', dap, '  Brand: ', ibr, '   "Interaction"-Table for riks windows of 28 days:\n\n'))
      print(  flowchart_simple_interactions_28d )
      # rw is 60d
      cat(paste0('\n\nDAP: ', dap, '  Brand: ', ibr, '   "Interaction"-Table for riks windows of 60 days:\n\n'))
      print(  flowchart_simple_interactions_60d )
      
      #options(old_width)
    
    sink()
    
    flowchart_list <- c( 
      flowchart_list, 
      list(list( 
            flowchart                  = format(flowchart_simple,  justify="left"),
            flowchart_interactions     = flowchart_simple_interactions,
            flowchart_interactions_28d = flowchart_simple_interactions_28d,
            flowchart_interactions_60d = flowchart_simple_interactions_60d,
            dap                        = dap,
            brand                      = ibr
    )))
    names(flowchart_list)[length(flowchart_list)] <- paste0(dap,"_",ibr)
  
  }  
  
  save(flowchart_list, file = paste0(diroutput,'scri/flowcharts/flowchart_',dap,'.RData') )
  
  
  #############
  #
  #        select included rows:
  #
  nrow0 <- nrow(scri_data_extract)
  scri_data_extract <- scri_data_extract[ scri_data_extract$include,  ]
  print( c(rows_from_control_window = nrow(scri_data_extract), old_nrow=nrow0, percent=100*nrow(scri_data_extract)/nrow0), digits=2 ) # without period between risk windows of vax1 and vax2 with buffer period
  
  
  scri_data_extract_subpop <- paste0("scri_data_extract", suffix[[subpop]])
  assign(scri_data_extract_subpop, scri_data_extract)
  save(scri_data_extract_subpop, file = paste0(dirtemp, "scri/scri_data_extract", suffix[[subpop]],".RData"), list = scri_data_extract_subpop)
  
  # save(scri_data_extract, file = paste0(dirtemp, "scri/scri_data_extract.RData"))
  
  
  # scri_data_extract <- scri_data_extract[ scri_data_extract$include & scri_data_extract$study_period_28d,  ] #  only in rw28
  # scri_data_extract <- scri_data_extract[ scri_data_extract$include & scri_data_extract$study_period_60d,  ] #  only in rw60
  #
  # scri_data_extract <- scri_data_extract[ scri_data_extract$include & scri_data_extract$study_period_28d & scri_data_extract$cond_not_buffer,  ] #  only in rw28 & without buffer
  # scri_data_extract <- scri_data_extract[ scri_data_extract$include & scri_data_extract$study_period_60d & scri_data_extract$cond_not_buffer,  ] #  only in rw60 & without buffer
  # 
  
  
  ###############################
  #
  #   Baseline characteristics:
  #
  
  for(iae in c("myopericarditis","myocarditis","pericarditis")){
    
    event_data <- scri_data_extract[scri_data_extract[,iae]==1,]
  
    for(i in 1:3){
      if(i==1) {
        bas_name <- "  included at day -90"
        bas_name_short <- paste0("_",iae,"_after_ms90d")
        bas_data <- event_data # inclusion at day -90
      }  
      if(i>1) {
        if(i==2) {                           # inclusion at day -90 and end 28 days after last vaccine
          vax1_end         <- 28
          vax2_end         <- 28
          bas_name <- "  included at day -90 and before last vaccine + 28 days"
          bas_name_short <- paste0("_",iae,"_after_ms90d_28d")
        }  
        if(i==3) {                          # inclusion at day -90 and end 60 days after last vaccine
          vax1_end         <- 60
          vax2_end         <- 60
          bas_name <- "  included at day -90 and before last vaccine + 60 days"
          bas_name_short <- paste0("_",iae,"_after_ms90d_60d")
        }  
        end_last_rw <- pmax(event_data$days_vax1 + vax1_end, event_data$days_vax2 + vax2_end, na.rm=T )
        end_last_rw <- pmin(end_last_rw, event_data$study_exit_days )
        bas_data    <- event_data [ event_data[,paste0(iae,"_days")] <= end_last_rw,]
      }
      
      # sex per brand:
      var <- "sex"
      bas <- c( list(all  = table1(bas_data[,var])),
                tapply(bas_data[,var], bas_data$type_vax1, table1) )
      assign(paste0(dap,bas_name_short,"_basl_sex"),bas); rm(bas)
      
      # age per brand:
      var <- "age_at_study_entry"
      age_cat  <- c(-1,30,120)
      age_cat2 <- c(-1,18,24,29,39,49,55,65,80,120)
      bas <-  list(summary = 
                     t(sapply(c(
                       all=list( c(summary(bas_data[,var],na.rm=T), n_missing=sum(is.na(bas_data[,var])))),
                       tapply(bas_data[,var], bas_data$type_vax1, function(x)c(summary(x,na.rm=T),n_missing=sum(is.na(x))))),c))  ,
                   quartile = t(sapply(c(
                     all=list( quantile(bas_data[,var],c(0,0.25,0.5,0.75,1),na.rm=T) ),
                     tapply(bas_data[,var], bas_data$type_vax1, quantile, c(0,0.25,0.5,0.75,1), na.rm=T )),c))  ,
                   cat1 = c(  all=list( table1(cut(bas_data[,var],age_cat))),
                              tapply(cut(bas_data[,var],age_cat), bas_data$type_vax1, table1)  ) ,
                   cat2 = c(  all=list( table1(cut(bas_data[,var],age_cat2))),
                              tapply(cut(bas_data[,var],age_cat2), bas_data$type_vax1, table1)  ) 
      )
      assign(paste0(dap,bas_name_short,"_basl_age"),bas); rm(bas)
      
      # age per brand and sex:
      bas <-  list(summary_per_sex = 
                     t(sapply(c(
                       tapply(bas_data[,var], paste0("all_",bas_data$sex),  function(x)c(summary(x,na.rm=T),n_missing=sum(is.na(x)))) ,
                       tapply(bas_data[,var], paste0(bas_data$type_vax1,"_",bas_data$sex), function(x)c(summary(x,na.rm=T),n_missing=sum(is.na(x)))) ),c)),
                   quartile_per_sex = 
                     t(sapply(c(
                       tapply(bas_data[,var], paste0("all_",bas_data$sex),   quantile, c(0,0.25,0.5,0.75,1), na.rm=T ) ,
                       tapply(bas_data[,var], paste0(bas_data$type_vax1,"_",bas_data$sex),  quantile, c(0,0.25,0.5,0.75,1), na.rm=T )),c)) ,
                   cat1_per_sex = 
                     c(  all=list( table1(cbind.data.frame(bas_data$sex, cut(bas_data[,var],age_cat)))),
                         tapply(cut(bas_data[,var],age_cat), paste0(bas_data$type_vax1,"_",bas_data$sex), table1)  ) 
                   )
      assign(paste0(dap,bas_name_short,"_basl_age_per_sex"),bas); rm(bas)
      
      var <- "dose_diff"
      days_betwee_cat <- c(-1,7,14,21,28,60,90,10000)
      bas <- c( brands                = list(table1(  bas_data[,c("type_vax1","type_vax2")], title="" )),
                   days_between_cat_all  = list(table1(  cut(bas_data[,var],days_betwee_cat))),
                   tapply(  cut(bas_data[,var],days_betwee_cat), paste0("days_between_cat_",bas_data$type_vax1), table1),
                   days_between_summary  = list(t(sapply(c(
                     all=list(summary( bas_data[,var], na.rm=T ) ),
                     tapply(  bas_data[,var], bas_data$type_vax1, summary, na.rm=T)),c)))  )
      assign(paste0(dap,bas_name_short,"_basl_vax1_vax2"),bas); rm(bas)
      
      
      
      save(list=c( paste0(dap,bas_name_short,"_basl_sex"),
                   paste0(dap,bas_name_short,"_basl_age"),
                   paste0(dap,bas_name_short,"_basl_age_per_sex"),
                   paste0(dap,bas_name_short,"_basl_vax1_vax2")), 
           file = paste0(thisdirexp, "scri/flowcharts/baseline_",dap,bas_name_short,".RData"))
      
      
      sink(paste0(thisdirexp, 'scri/flowcharts/baseline_',dap,bas_name_short,'.txt'),append = F)
      
      cat(paste0('\n\nDAP: ', dap, bas_name, '  sex:\n\n'))
      print(get(paste0(dap,bas_name_short,"_basl_sex")))
      
      cat(paste0('\n\nDAP: ', dap, bas_name, '  age:\n\n'))
      print(get(paste0(dap,bas_name_short,"_basl_age")))
      
      cat(paste0('\n\nDAP: ', dap, bas_name, '  age per sex:\n\n'))
      print(get(paste0(dap,bas_name_short,"_basl_age_per_sex")))
      
      cat(paste0('\n\nDAP: ', dap, bas_name, '  vax1 and vax2:\n\n'))
      print(get(paste0(dap,bas_name_short,"_basl_vax1_vax2")))
      
      sink()
      
    }
  }
  
  ##########
  # restore options:
  options(old_width)
  
  
  
  
  
  ############
  #
  # flowchart LARGE:
  #
  if(F){
    dap <- ifelse( any(names(scri_data_extract)=="DAP"), scri_data_extract$DAP[1], "")
    sink(paste0('./g_output/scri/flowchart_',dap,'.txt'),append = F)
      cat("\n\n\n")
      dim(scri_data_extract)
      # non-missing vax1 dates
      na_date_vax1 <- c("no_vax1","vax1")[ 1 + as.numeric(!is.na(scri_data_extract$date_vax1)) ]
      print(table1(na_date_vax1, title="vaccine: yes or no ('date_vax1' is missing):"))
      cat("___________________________\n\n\n")
      # .. and brand
      flowchart <- cbind.data.frame(na_date_vax1=na_date_vax1, brand=scri_data_extract$type_vax1) 
      print(table1( flowchart$brand, title="Brands:"))
      print(table1( flowchart ))
      cat('___________________________\n\n\n')
      # .. and sex
      flowchart <- cbind.data.frame(flowchart, sex=scri_data_extract$sex) 
      print(table1(flowchart$sex, title="Gender (M,F  or 0=F,1=M)"))
      print(table1(flowchart))
      print(tapply(flowchart$sex, paste0(flowchart$na_date_vax1," & ",flowchart$brand),  table1) )
      cat("___________________________\n\n\n")
      # .. and age at start_study (01-09-2020)
      age_cat_fc <- age_cat
      if( min(scri_data_extract$age_at_study_entry,na.rm=T) < age_cat[1] ) age_cat_fc <- c(min(scri_data_extract$age_at_study_entry,na.rm=T), age_cat_fc)
      if( max(scri_data_extract$age_at_study_entry,na.rm=T) > age_cat[length(age_cat)] ) age_cat_fc <- c(age_cat_fc, max(scri_data_extract$age_at_study_entry,na.rm=T) )
      flowchart <- cbind.data.frame(flowchart, age_start_fu=cut(scri_data_extract$age_at_study_entry, age_cat_fc) ) 
      print(table1(flowchart$age_start_fu, title="age_at_start_study categories:"))
      print(table1(flowchart))
      print(tapply(paste(flowchart$sex," & ",flowchart$age_start_fu), paste0(flowchart$na_date_vax1," & ",flowchart$brand),  table1) )
      cat("___________________________\n\n\n")
      # obs.start after 1-09-2020:
      obs_before_fu_or_baby <- rep("", nrow(scri_data_extract))
      obs_before_fu_or_baby[ is.na(scri_data_extract$study_entry_date)  ] <- "start_obs is missing"
      obs_before_fu_or_baby[ scri_data_extract$study_entry_date <  scri_data_extract$start_study_date ] <- "start_obs before fu_start"
      obs_before_fu_or_baby[ scri_data_extract$study_entry_date >= scri_data_extract$start_study_date ] <- "start_obs after fu_start (no_vax1)"
      obs_before_fu_or_baby[ scri_data_extract$study_entry_date >= scri_data_extract$start_study_date & 
                               scri_data_extract$study_entry_date < scri_data_extract$date_vax1-90 ] <- "start_obs in [start_fu;vax1-90d)"
      obs_before_fu_or_baby[ !is.na(scri_data_extract$age_at_study_entry) & scri_data_extract$age_at_study_entry<=0 ] <- "born after fu_start"
                               
      flowchart <- cbind.data.frame(flowchart, obs_before_fu_or_baby = obs_before_fu_or_baby)
                                    
      print(table1(flowchart$obs_before_fu_or_baby, title="Table for study_entry_date or start_study_date:"))
      print(table1(flowchart))
      print(tapply(paste0(flowchart$obs_before_fu_or_baby," & age:",flowchart$age_start_fu), 
                   paste0(flowchart$na_date_vax1," & ",flowchart$brand, " & sex:",flowchart$sex),  table1) )
      #print(table1(!is.na(scri_data_extract$study_entry_date)))
      print("summary(study_entry_date):")
      print(summary(scri_data_extract$study_entry_date))
      #print(summary(as.numeric(difftime(scri_data_extract$start_study_date, scri_data_extract$study_entry_date, units="days"))))
      cat("___________________________\n\n\n")
      # myoperi dates in fu period:
      myoperi_before_fu <- rep("", nrow(scri_data_extract))
      myoperi_before_fu[ is.na(scri_data_extract$myopericarditis_date)  ] <- "myoperi is missing"
      myoperi_before_fu[ scri_data_extract$myopericarditis_date < scri_data_extract$start_study_date ] <- "myoperi   before fu"
      myoperi_before_fu[ scri_data_extract$myopericarditis_date >= scri_data_extract$start_study_date ] <- "myoperi after start_fu (no_vax1)"
      myoperi_before_fu[ scri_data_extract$myopericarditis_date >= scri_data_extract$date_vax1-90 ] <- "myoperi after (vax1-90d)"
      myoperi_before_fu[ scri_data_extract$myopericarditis_date >= scri_data_extract$start_study_date & 
                           scri_data_extract$myopericarditis_date < scri_data_extract$date_vax1-90 ] <- "myoperi  in [start_fu;vax1-90d)"
      #myoperi_before_fu[ scri_data_extract$myopericarditis_date > scri_data_extract$date_vax1-90 ] <- "myoperi in [start_fu;vax1-90d]"
      flowchart <- cbind.data.frame(flowchart, myoperi_before_fu = myoperi_before_fu)
      
      print(table1(flowchart$myoperi_before_fu,title="Table for myo-peri-carditis:"))
      print(table1(flowchart[,c("obs_before_fu_or_baby","myoperi_before_fu")]))
      print(table1(flowchart))
      print(tapply(paste(flowchart$myoperi_before_fu," & ",flowchart$age_start_fu), paste0(flowchart$na_date_vax1," & ",flowchart$brand, " & ",flowchart$sex),  table1) )
      cat("\nAbout 'study_entry_date':\n\t#non missing:",sum(!is.na(scri_data_extract$study_entry_date)), "\n\tsummary:\n" )
      print(summary(scri_data_extract$study_entry_date))
      #print(summary(as.numeric(difftime(scri_data_extract$start_study_date, scri_data_extract$study_entry_date, units="days"))))
      cat("___________________________\n\n\n")
    
      # condition
      cond_scri <- with( scri_data_extract,
                         !is.na(date_vax1) & 
                           !is.na(type_vax1) &   
                           !is.na(sex) &   
                           !is.na(age_at_study_entry) &   
                           !is.na(myopericarditis_date) & 
                           myopericarditis_date >=  study_entry_date  #  or start_study_date ? # can be improtant for babies!  
      )
      cond_scri_90d <- cond_scri & 
        scri_data_extract$myopericarditis_date >= scri_data_extract$date_vax1-90 
      
      prevax_per_start <- -90
      vax1_end <- 28
      vax2_end <- 28
      
      # condition for myoperi before start:
      data$cond_from_control_window <- data$myopericarditis_days >= data$days_vax1 + prevax_per_start 
      
      # condition myoperi during buffer interval
      data$cond_not_buffer <-  
        !( data$days_vax1 + buffer_per_start <=  data$myopericarditis_days  & 
           data$myopericarditis_days         <   data$days_vax1  )
      
      # condition for between risk windows of vax1 and vax2
      data$cond_not_betw_rw1_rw2_28d <-  
        data$vax2==1     & 
        !( data$days_vax1 + vax1_end < data$myopericarditis_days  & 
           data$myopericarditis_days <   data$days_vax2 )
      
      # condition for myoperi after last risk window (of vax1 or vax2):
      tmp <- pmax(data$days_vax1 + vax1_end, data$days_vax2 + vax2_end )
      data$cond_till_last_rw_28d <- tmp >= data$myopericarditis_days  
      
      
      vax1_end <- 60
      vax2_end <- 60
      # condition for between risk windows of vax1 and vax2
      data$cond_not_betw_rw1_rw2_60d <-  
        data$vax2==1     & 
        !( data$days_vax1 + vax1_end < data$myopericarditis_days  & 
           data$myopericarditis_days <   data$days_vax2 ) 
      
      # condition for myoperi after last risk window (of vax1 or vax2):
      tmp <- pmax(data$days_vax1 + vax1_end, data$days_vax2 + vax2_end )
      data$cond_till_last_rw_60d <- tmp >= data$myopericarditis_days  
    
      
      # condition
      cond_scri <- with( scri_data_extract,
                         !is.na(date_vax1) & 
                           !is.na(type_vax1) &   
                           !is.na(sex) &   
                           !is.na(age_at_study_entry) &   
                           !is.na(myopericarditis_date) & 
                           myopericarditis_date >=  study_entry_date  #  or start_study_date ? # can be improtant for babies!  
      )
      cond_scri_90d <- cond_scri & 
                       scri_data_extract$myopericarditis_date >= scri_data_extract$date_vax1-90 
                         
      flowchart <- cbind.data.frame( cond_scri_90d=c("excluded"," included")[1+as.numeric(cond_scri_90d)], flowchart )   
      
      cat("\n\n\nFull Flowchart:\n\n")
      print(table1( flowchart$cond_scri_90d, title="#persons included till now:" ))
      print(tapply( flowchart$cond_scri_90d, flowchart$brand, table1 ))
      print(table1( flowchart ))
      cat("___________________________\n\n\n")
      
      cat("\n\n\nFull Flowchart only for vaccinated:\n\n")
      print(table1( flowchart$cond_scri_90d[flowchart$na_date_vax1!="no_vax1"], title="#vaccinated persons included till now:" ))
      print(with(flowchart[flowchart$na_date_vax1!="no_vax1",], tapply( cond_scri_90d, brand, table1) ))
      print(table1( flowchart[flowchart$na_date_vax1!="no_vax1",] ))
      cat("___________________________\n\n\n")
      
      
      print(table1(flowchart[,c("obs_before_fu_or_baby","myoperi_before_fu")]))
      cat("___________________________\n\n\n")
      
      # include only 90 days before vax1:
      print("Selection persons with 90 days before vax1:")
      nrow0<-nrow(scri_data_extract0)
      scri_data_extract <- scri_data_extract0[cond_scri_90d, ]
      print(cbind.data.frame( included_rows=nrow(scri_data_extract), all_rows=nrow0, percent_included= 100*nrow(scri_data_extract)/nrow0 ), digits=2)
      
      
      # include persons from start_fu:
      print("Selection persons from the fu start:")
      nrow0<-nrow(scri_data_extract0)
      scri_data_extract_from_fu <- scri_data_extract0[cond_scri, ]
      print(cbind.data.frame( included_rows=nrow(scri_data_extract_from_fu), all_rows=nrow0, percent_included= 100*nrow(scri_data_extract_from_fu)/nrow0 ), digits=2)
      
    sink()
  
    
    save(scri_data_extract_subpop, file = paste0(dirtemp, 'scri/flowchart_',dap, suffix[[subpop]],'.RData'), list = scri_data_extract_subpop)
    
    assign(scri_data_extract_subpop, scri_data_extract)
    save(scri_data_extract_subpop, file = paste0(dirtemp, "scri/scri_data_extract", suffix[[subpop]],".RData"), list = scri_data_extract_subpop)
    rm(list = c(scri_data_extract_subpop))
    
    # save(scri_data_extract, file = paste0(dirtemp, 'scri/flowchart_',dap,'.RData'))
    
    # save(scri_data_extract, file = paste0(dirtemp, "scri/scri_data_extract.RData"))
    
  }
  
  
  
  save(nvax, file = paste0(dirtemp, "nvax", suffix[[subpop]], ".RData"))

}
