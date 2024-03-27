# Program Information  ----------------------------------------------------
#
# Program:      01_sensitivity_clean_data.R 
# Description:  this script does some data cleaning relevant for all designs, and outputs a dataset for 
#               SCRI sensitivity analyses. 
#               It uses some code from Svetlana's original programming to ensure consistency in the data cleaning 
# Requirements: dataset in RData format called D3_study_population_SCRI, note, should be in g_intermediate
# Output:       this outputs a single dataset called sccs_data_extract to g_intermediate 

# 0. HOUSEKEEPING ---------------------------------------------------------
# set in 00_functions file 

# Start the subpopulation
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  # 1. IMPORT DATA ----------------------------------------------------------
  
  # Import data 
  message("Note: Initiating data cleaning")
  scri_data_extract <- as.data.frame(get(load(paste0(dirtemp, raw_data_name, suffix[[subpop]], ".RData"))[[1]]))
  
  # 2. DATA CLEANING ------------------------------------------------------
  ## data cleaning which should apply to all case series, irrespective of DAP and outcome 
  ## QC and sense checks 
  
  # Check key variables exist and are expected format, otherwise throw warning 
  key_vars <- c("date_vax1", "date_vax2", "date_vax3", "study_entry_date", "study_exit_date")
  
  # check the vaccine, entry, and exit dates exist 
  invisible(lapply(key_vars, function(x) if (!(x %in% colnames(scri_data_extract))) {
    stop(paste(x," does not exist"))
  }))
  
  # Harmonise variable names; remove underscores
  oldnames <- c(names(scri_data_extract[grepl("vax", names(scri_data_extract))]),
                names(scri_data_extract[grepl("covid", names(scri_data_extract))]))
  newnames <- gsub("vax_", "vax", oldnames)
  newnames <- gsub("covid_19|covid_date", "covid19", newnames)
  names(scri_data_extract)[match(oldnames, names(scri_data_extract))] <- newnames
  
  # create list of date variable names for later use (except the study start date)
  date_vars <- colnames(scri_data_extract)[grepl("date", colnames(scri_data_extract))==T]
  date_vars <- date_vars[date_vars != "start_study_date"]
  date_vars_days <- gsub("date","days",date_vars)
  
  # convert numeric to date (dates may be numeric in BIFAP)
  scri_data_extract[date_vars] <- lapply(scri_data_extract[date_vars],
                                         function(x) if (is.numeric(x)) as.Date(x, origin="1970-01-01") else x)
  
  # check for characters 
  invisible(lapply(colnames(scri_data_extract[date_vars]), function(x) {
    if (is.character(scri_data_extract[[x]])) {
      message(paste(x, "is character, check the format for converting dates as program assumes YYYY-MM-DD"))
    }
  }))
  
  # convert character to date (assumed format YYYY-MM-DD)
  scri_data_extract[date_vars] <- lapply(scri_data_extract[date_vars],
                                         function(x) if (is.character(x)) as.Date(x) else x)
  
  # create a variable for study start date 
  # only relevant where someone enters late, not applicable to this study as there is no changing eligibility in age which would allow someone late entry 
  start_scri_date <- as.Date("2020-09-01")
  
  # create vaccine indicator variables 
  scri_data_extract$vax1 <- !is.na(scri_data_extract$date_vax1)
  scri_data_extract$vax2 <- !is.na(scri_data_extract$date_vax2)
  scri_data_extract$vax3 <- !is.na(scri_data_extract$date_vax3)
  
  # reset study start date to be the entry date 
  scri_data_extract$start_study_date <- pmax(scri_data_extract$study_entry_date, start_scri_date, na.rm = T)
  # below code is retained for run on simulated data as this only has people entering after the study start date 
  #scri_data_extract$start_study_date <- start_scri_date 
  
  scri_data_extract$start_study_days <- round(difftime(scri_data_extract$start_study_date, as.Date("2020-09-01"), units = "days"),0)

  # get max number of vaccine doses
  nvax <- names(scri_data_extract)[ substring(names(scri_data_extract),1,8)=="date_vax" ]
  nvax <- max(as.numeric(substring(nvax,9)))
  
  # Basic dataset information 
  message(paste("Note: there are", nrow(scri_data_extract), "individuals in this dataset"))
  message(paste("Note: the maximum number of vaccine doses was", nvax))
  
  # flag how common late entry into the dataset is, anticipated to be very low or zero 
  if(any(scri_data_extract$study_entry_date > start_scri_date)) 
    message(paste("Note: there are", sum(scri_data_extract$study_entry_date > start_scri_date), "persons with study entry date after the SCRI start date"))
  
  # change the date variables to days from the calendar start date of our observation period 
  scri_data_extract[date_vars_days] <- lapply(scri_data_extract[date_vars],
                                              function(x){round(difftime(x, start_scri_date, units = "days"), 0)})
  
  # Confirm that the first vaccine dose occurs after the study entry date, remove rows if this does not happen
  if(any(scri_data_extract$start_study_date > scri_data_extract$date_vax1 & !is.na(scri_data_extract$date_vax1) & !is.na(scri_data_extract$start_study_date), na.rm=T )){
    message(paste("Note: 'start_study_date' after the vax1 for ", sum(scri_data_extract$start_study_date > scri_data_extract$date_vax1, na.rm=T), "rows. They are deleted!"))
    nrow0<-nrow(scri_data_extract)  
    scri_data_extract <- scri_data_extract[ scri_data_extract$start_study_date <= scri_data_extract$date_vax1 & scri_data_extract$vax1==1 & !is.na(scri_data_extract$start_study_date),]  # if study_entry_days > vax1 ==>  delete rows
    message(paste("new row", nrow(scri_data_extract), "old row", nrow0)) 
  }
  
  # Check the time difference (in days) between vaccine doses
  # It should be at least 14 days; print warning if this is not the case
  scri_data_extract$dose_diff = as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1,units="days"))
  message("Note: Distribution of time between dose 1 and 2 is printed in log")
  summary(scri_data_extract$dose_diff)
  if(any(min(scri_data_extract$dose_diff[!is.na(scri_data_extract$dose_diff)])<14)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff<14, na.rm=T), " persons with (date_vax2 - date_vax1) < 14 days."))
  
  # number of people with third dose before second dose, reset date to second dose plus 1 if not
  if(any((cond <- !is.na(scri_data_extract$days_vax3) & !is.na(scri_data_extract$days_vax2) & scri_data_extract$days_vax3 < scri_data_extract$days_vax2))){
    message(paste0("Note: There are ",sum( cond, na.rm=T)," persons with vax3 before vax2! date of dose 3 set to dose 2 + 1"))
    scri_data_extract$days_vax3[cond] <- scri_data_extract$days_vax2[cond] + 1 
  }

  # do the same for dose 3 to dose 2
  # add an extra check to check whether the planned control time of 88 days is available 
  scri_data_extract$dose_diff2 = as.numeric(difftime(scri_data_extract$date_vax3, scri_data_extract$date_vax2,units="days"))
  message("Note: Distribution of time between dose 2 and 3 is printed in log")
  summary(scri_data_extract$dose_diff2)
  if(any(min(scri_data_extract$dose_diff2[!is.na(scri_data_extract$dose_diff2)])<14)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff2<14, na.rm=T), " persons with (date_vax3 - date_vax2) < 14 days."))
  if(any(min(scri_data_extract$dose_diff2[!is.na(scri_data_extract$dose_diff2)])<88)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff2<14, na.rm=T), " persons with (date_vax3 - date_vax2) < 88 days (needed control time for post vaccine SCRI)"))
  
  # Check whether someone's study exit date occurs before the vaccine, remove rows if this does happen
  if(any(scri_data_extract$study_exit_days < scri_data_extract$days_vax1 & !is.na(scri_data_extract$days_vax1), na.rm = T)) {
    message(paste0("Note: There are ",sum(scri_data_extract$study_exit_days < scri_data_extract$days_vax1,na.rm=T)," persons with 'study_exit_time' before vax1. They are deleted!"))
    nrow0<-nrow(scri_data_extract)  
    scri_data_extract <- scri_data_extract[scri_data_extract$days_vax1 <= scri_data_extract$study_exit_days & !is.na(scri_data_extract$days_vax1),]  # if study_exit_days < vax1 ==>  delete rows
  }  
  
  # number of people with second vaccine dose before end of study (cond is just used in the check for clarity)
  if(any((cond <- !is.na(scri_data_extract$days_vax2) & scri_data_extract$study_exit_days < scri_data_extract$days_vax2) )){
    message(paste0("Note: There are ",sum( cond, na.rm=T)," persons with 'study_exit_time' before vax2! Exit date is set to second dose date"))
    scri_data_extract$study_exit_days[cond] <- scri_data_extract$days_vax2[cond] 
  }

  # number of people with third dose before end of study (change exit date to third dose)
  if(any((cond <- !is.na(scri_data_extract$days_vax3) & scri_data_extract$study_exit_days < scri_data_extract$days_vax3) )){
    message(paste0("Note: There are ",sum( cond, na.rm=T)," persons with 'study_exit_time' before vax3! Exit date is set to third dose date"))
    scri_data_extract$study_exit_days[cond] <- scri_data_extract$days_vax3[cond] 
  }

  # number of people with myocarditis after end of study (cond is just used in the check for clarity)
  if(any(scri_data_extract$study_exit_date < scri_data_extract$myocarditis_date & !is.na(scri_data_extract$myocarditis_date), na.rm = T)) {
    message(paste0("Note: There are ",sum(scri_data_extract$study_exit_date < scri_data_extract$myocarditis_date,na.rm=T)," persons with 'study_exit_time' before myocarditis_date."))
  }  
  
  # cleaning of vaccine doses (only applies if more than 2)
  # follows conventions set for myo/perio analyses, so have retained 
  if(nvax > 2){
    while(any(cond <- !is.na(scri_data_extract$days_vax2) & (scri_data_extract$days_vax2 - scri_data_extract$days_vax1) < 5)){
      message(paste(sum(cond),"'dose2' were replace with next dose because dose2 is less than 5 days after dose1."))
      for(iv in 3:nvax){
        scri_data_extract[cond, paste0("date_vax",iv-1) ] <- scri_data_extract[cond, paste0("date_vax",iv) ]
        scri_data_extract[cond, paste0("days_vax",iv-1) ] <- scri_data_extract[cond, paste0("days_vax",iv) ]
        scri_data_extract[cond, paste0("type_vax",iv-1) ] <- scri_data_extract[cond, paste0("type_vax",iv) ]
        scri_data_extract[cond, paste0("vax"     ,iv-1) ] <- scri_data_extract[cond, paste0("vax"     ,iv) ]
      } 
    }  
  }
  
  
  # update the basic dataset summaries after the above cleaning 
  scri_data_extract$dose_diff  <-  as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1 ,units="days"))
  
  message("Distribution of time between dose 1 and 2 after data cleaning")
  summary(scri_data_extract$dose_diff)
  if(any(min(scri_data_extract$dose_diff[!is.na(scri_data_extract$dose_diff)])<14)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff<14, na.rm=T), " persons with (date_vax2 - date_vax1) < 14 days after data cleaning"))
  
  scri_data_extract$dose_diff2 = as.numeric(difftime(scri_data_extract$date_vax3, scri_data_extract$date_vax2,units="days"))
  message("Distribution of time between dose 2 and 3 after data cleaning")
  summary(scri_data_extract$dose_diff2)
  if(any(min(scri_data_extract$dose_diff2[!is.na(scri_data_extract$dose_diff2)])<14)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff2<14, na.rm=T), " persons with (date_vax3 - date_vax2) < 14 days after data cleaning"))
  if(any(min(scri_data_extract$dose_diff2[!is.na(scri_data_extract$dose_diff2)])<88)) 
    message(paste("Note: There are ", sum(scri_data_extract$dose_diff2<14, na.rm=T), " persons with (date_vax3 - date_vax2) < 88 days after data cleaning"))
  
  # create a variable for the last possible vaccine date for the post vaccine SCRI analyses (max of dose 1 and 2, because control time will start after dose 2)
  scri_data_extract$days_last_vax <- pmax(scri_data_extract$days_vax1, scri_data_extract$days_vax2, na.rm = T)
  
  # create myopericarditis date variable 
  scri_data_extract$myopericarditis_date <- pmin(scri_data_extract$myocarditis_date, scri_data_extract$pericarditis_date,na.rm=T)
  
  # what is the time distribution between doses? 
  scri_data_extract$diff_d1d2 <- as.numeric(scri_data_extract$date_vax2 - scri_data_extract$date_vax1)
  scri_data_extract$diff_d2d3 <- as.numeric(scri_data_extract$date_vax3 - scri_data_extract$date_vax2)
  
  # number of people with vaccine after event 
  scri_data_extract$dose1_afterevent <- (scri_data_extract$myocarditis_date <= scri_data_extract$date_vax1)
  scri_data_extract$dose2_afterevent <- ((scri_data_extract$myocarditis_date <= scri_data_extract$date_vax2) & !is.na(scri_data_extract$date_vax2)) 
  scri_data_extract$dose3_afterevent <- ((scri_data_extract$myocarditis_date <= scri_data_extract$date_vax3) & !is.na(scri_data_extract$date_vax3)) 
  
  # Time between event and doses?
  scri_data_extract$diff_event_d1 <- as.numeric(scri_data_extract$myocarditis_date - scri_data_extract$date_vax1) 
  scri_data_extract$diff_event_d2 <- as.numeric(scri_data_extract$myocarditis_date - scri_data_extract$date_vax2)
  scri_data_extract$diff_event_d3 <- as.numeric(scri_data_extract$myocarditis_date - scri_data_extract$date_vax3) 
  
  # create a numeric ID variable
  scri_data_extract$numeric_id <- match(scri_data_extract$person_id, unique(scri_data_extract$person_id))
  
  ### throw a warning if duplicate PIDs 
  if(any(duplicated(scri_data_extract$numeric_id))) { 
    warning("there are unexpected duplicate PIDs in the data, check with DAP")
    message("Note: duplicate PIDs dropped to enable function to run")
    scri_data_extract <- scri_data_extract %>% 
      distinct(numeric_id, .keep_all = TRUE) # should really be dropped by DAP 
  }
  
  # 3. SAVE DATA ------------------------------------------------------------
  # saving data with new name to prevent overwriting any prior files
  sccs_data_extract <- scri_data_extract 
  save(sccs_data_extract,file = paste0(dirtemp,"sccs_sensitivity/","sccs_data_extract", suffix[[subpop]], ".RData"))
  
}
