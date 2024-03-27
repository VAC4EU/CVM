# Program Information  ----------------------------------------------------
#
# Program:      04d_standard_sccs.R
# Description:  run standard SCCS
# Dependencies: 00_sensitivity_functions 
#               01_clean_data.R
#               02_select_population.R 
# 
# 0. HOUSEKEEPING ---------------------------------------------------------

# Start the subpopulation
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  # specify the study design
  design <- "extended_sccs"
  
  # 1. READ IN DATA ------------------------------------------------------------
  load(paste0(dirtemp, "sccs_sensitivity/", "sccs_population", suffix[[subpop]], ".RData"))
  dir.create(file.path(paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", design)),         showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(paste0(diroutput, "sccs_sensitivity/", design, suffix[[subpop]])),         showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(sccs_population) < 5) {
    
    message(paste("There are only", nrow(sccs_population), "people, execution of script halted"))
    
  } else {
    
    # STUDY INPUTS ------------------------------------------------------------
    # specify study inputs and create generic variables 
    
    ## study design options
    preexp_start  <- -29
    vax1_end          <- 28
    vax2_end          <- 28
    vax3_end          <- 28
    
    ## outcome 
    sccs_population$outcome_date <- sccs_population$myocarditis_date
    sccs_population$outcome_days <- round(difftime(sccs_population$outcome_date, as.Date("2020-09-01"), units = "days"),0)
    sccs_population$outcome_days <- as.numeric(sccs_population$outcome_days)
    sccs_population$outcome_binary <- as.numeric(!is.na(sccs_population$outcome_date))
    
    # CREATE VARIABLES --------------------------------------------------------
    # create variables (in days) needed to fit the chosen study design  
    
    # pre-exposure period start, end and duration 
    sccs_population$preexp1 <- as.numeric(sccs_population$days_vax1 + preexp_start)
    sccs_population$preexp1_end <- as.numeric(sccs_population$days_vax1)
    # start minus end + 1 as the model fit is inclusive (so if start and end is the same, modelled length is 1 day)
    sccs_population$preexp1_length <- as.numeric(sccs_population$preexp1_end - sccs_population$preexp1) + 1
    
    sccs_population$preexp2 <- as.numeric(sccs_population$days_vax2 + preexp_start)
    sccs_population$preexp2_end <- as.numeric(sccs_population$days_vax2)
    sccs_population$preexp2_length <- as.numeric(sccs_population$preexp2_end - sccs_population$preexp2) + 1
    
    sccs_population$preexp3 <- as.numeric(sccs_population$days_vax3 + preexp_start)
    sccs_population$preexp3_end <- as.numeric(sccs_population$days_vax3)
    sccs_population$preexp3_length <- as.numeric(sccs_population$preexp3_end - sccs_population$preexp3) + 1
    
    # first risk window start, end and duration. Censor risk period at day before 2nd dose if occurs, as day of 2nd dose should be considered separately 
    sccs_population$risk_d1 <- as.numeric(sccs_population$days_vax1 + 1)
    sccs_population$risk_d1_end <- as.numeric(pmin((sccs_population$days_vax1 + vax1_end), sccs_population$study_exit_days, sccs_population$days_vax2-1, na.rm = T))
    # plus 1 to account for inclusive model fit
    sccs_population$risk_d1_length <- as.numeric((sccs_population$risk_d1_end - sccs_population$risk_d1)) + 1
    
    # second risk window start, end and duration. Censor risk period at 3rd dose if occurs. 
    sccs_population$risk_d2 <- as.numeric(sccs_population$days_vax2 + 1)
    # 2nd dose end min of dose 2 risk period, study exit and dose 3. Involved call because want to ignore missing dose 3 but return missing if missing dose 2 
    sccs_population$risk_d2_end <- ifelse(!is.na(sccs_population$days_vax2), as.numeric(pmin((sccs_population$days_vax2 + vax2_end), sccs_population$study_exit_days, sccs_population$days_vax3, na.rm = T)), NA_real_)
    # plus 1 to account for inclusive model fit
    sccs_population$risk_d2_length <- as.numeric((sccs_population$risk_d2_end - sccs_population$risk_d2)) + 1
    
    # third risk window start, end and duration. Censor risk period at 3rd dose if occurs. 
    sccs_population$risk_d3 <- as.numeric(sccs_population$days_vax3 + 1)
    sccs_population$risk_d3_end <- ifelse(!is.na(sccs_population$days_vax3), as.numeric(pmin((sccs_population$days_vax3 + vax3_end), sccs_population$study_exit_days, na.rm = T)), NA_real_)
    sccs_population$risk_d3_length <- as.numeric((sccs_population$risk_d3_end - sccs_population$risk_d3)) + 1
    
    sccs_population$study_length <- as.numeric(sccs_population$study_exit_days)
    
    # start time set in cleaning file as all vars calculate in days relative to this time 
    sccs_population$ref_start <- sccs_population$start_study_days
    sccs_population$ref_end <- sccs_population$study_exit_days
    
    # end at the last of the end of the control window, first and second risk windows (depending on design)
    sccs_population$ref_end <- as.numeric(sccs_population$ref_end)
    
    ## calendar time adjustment 
    # create a calendar time variable for calendar time adjustment (30-day interval between start and end in the data where model is fit)
    max_day = max(sccs_population$ref_end, na.rm = T)
    min_day = min(sccs_population$ref_start, na.rm = T)
    
    # generate a cut-off for each month to allow us to summarise events per month
    caltime_groups <- seq(from = min_day-10, to = max_day-10, by = 60)
    
    # first cut-off is start of second age group, so create the same minus first element
    caltime_cutoffs_for_model <- caltime_groups[-1]
    
    # some groups may not have events. Find out which ones 
    outcome_sums <- numeric(length(caltime_groups) - 1)
    
    for (i in 1:(length(outcome_sums))) {
      outcome_sums[i] <- sum(between(sccs_population$outcome_days, caltime_groups[i], caltime_groups[i+1]), na.rm = T)
    } 
    
    ref_max <- as.numeric(which.max(outcome_sums))
    
    # NUMBER OF EVENTS --------------------------------------------------------
    # function to summarise the number and duration of each period 
    
    table3a.interim <- lapply(split(sccs_population, sccs_population$type_vax1), describe_sccs_periods)
    table3a <- do.call(rbind, table3a.interim )
    
    table3a <- cbind(Names = rownames(table3a), table3a)
    table3a$Names <- gsub("\\..*","",table3a$Names)
    
    write.csv(table3a, paste0(direxpsubpop[[subpop]],"sccs_sensitivity/", design, "/table_3a_describe_period_by_vaccine.csv"), row.names = FALSE)
    
    # FIT MODEL  --------------------------------------------------------------
    
    ## fit model 
    unadjusted <- lapply(split(sccs_population, sccs_population$type_vax1), quietly(fit_unadjusted_extended_sccs))
    adjusted <- lapply(split(sccs_population, sccs_population$type_vax1), quietly(fit_adjusted_extended_sccs))
    
    ## process output
    
    # print all messages and errors to a separate file ('conditions') 
    ## unadjusted 
    conditions <- paste0(diroutput,"sccs_sensitivity/", design, suffix[[subpop]], "/conditions_unadjusted.txt") 
    cat(paste(design, ": unadjusted model conditions\n\n"), file = conditions) 
    invisible(lapply(seq_along(unadjusted), process_conditions, inputfile = unadjusted, filename = conditions))
    
    ## adjusted 
    conditions <- paste0(diroutput,"sccs_sensitivity/", design, suffix[[subpop]], "/conditions_adjusted.txt") 
    cat(paste(design, ": adjusted model conditions\n\n"), file = conditions) 
    invisible(lapply(seq_along(adjusted), process_conditions, inputfile = adjusted, filename = conditions))
    
    # process results
    modeloutput1.list <- lapply(seq_along(unadjusted), process_output, inputfile = unadjusted)
    modeloutput2.list <- lapply(seq_along(adjusted), process_output, inputfile = adjusted)
    
    modeloutput1.list <- lapply(modeloutput1.list, function(x) if(!is.null(x)) ungroup(x))
    modeloutput2.list <- lapply(modeloutput2.list, function(x) if(!is.null(x)) ungroup(x))
    
    modeloutput1 <- do.call(rbind, modeloutput1.list)
    modeloutput2 <- do.call(rbind, modeloutput2.list)
    
    table3b <- rbind(modeloutput1, modeloutput2)
    
    # export results 
    write.csv(table3b, paste0(direxpsubpop[[subpop]],"sccs_sensitivity/", design, "/table_3b_model_results.csv"), row.names = FALSE)
    
  }
}