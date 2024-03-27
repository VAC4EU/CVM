#  MODEL EVALUATION  ----------------------------------------------------
# a series of different models to check that model specification is correct 
# these are intended to be run interactively, and for results to be viewed in the console 

# 1. PREPARATION

  # READ IN DATA ------------------------------------------------------------
  # QUESTION to Davide: the suffix added here doesn't change the name of the loaded dataset, is that correct? 
  # as below I use just the 'sccs_population' dataset to create all variables. I think if there was an issue it would have not run before...
  load(paste0(dirtemp, "sccs_sensitivity/", "sccs_population", suffix[[subpop]], ".RData"))
 
  # STUDY INPUTS ------------------------------------------------------------
  # specify study inputs and create generic variables 
  
  ## study design options
  control_duration  <- -59
  preexp_start  <- -29
  vax1_end          <- 28
  vax2_end          <- 28
  
  ## outcome 
  sccs_population$outcome_date <- sccs_population$myocarditis_date
  sccs_population$outcome_days <- round(difftime(sccs_population$outcome_date, as.Date("2020-09-01"), units = "days"),0)
  sccs_population$outcome_days <- as.numeric(sccs_population$outcome_days)
  sccs_population$outcome_binary <- as.numeric(!is.na(sccs_population$outcome_date))
  
  # CREATE VARIABLES --------------------------------------------------------
  # create variables (in days) needed to fit the chosen study design  
  
  # pre-exposure period start, end and duration 
  sccs_population$preexp <- as.numeric(sccs_population$days_vax1 + preexp_start)
  sccs_population$preexp_end <- as.numeric(sccs_population$days_vax1)
  # start minus end + 1 as the model fit is inclusive (so if start and end is the same, modelled length is 1 day)
  sccs_population$preexp_length <- as.numeric(sccs_population$preexp_end - sccs_population$preexp) + 1
  
  # first risk window start, end and duration. Censor risk period at day before 2nd dose if occurs, as day of 2nd dose should be considered separately 
  sccs_population$risk_d1 <- as.numeric(sccs_population$days_vax1 + 1)
  sccs_population$risk_d1_end <- as.numeric(pmin((sccs_population$days_vax1 + vax1_end), sccs_population$study_exit_days, sccs_population$days_vax2-1, na.rm = T))
  # plus 1 to account for inclusive model fit
  sccs_population$risk_d1_length <- as.numeric((sccs_population$risk_d1_end - sccs_population$risk_d1)) + 1
  
  # second risk window start, end and duration. Censor risk period at 3rd dose if occurs. 
  sccs_population$risk_d2 <- as.numeric(sccs_population$days_vax2 + 1)
  # if else needed here to not accidentally populate if dose 2 is missing 
  sccs_population$risk_d2_end <- ifelse(!is.na(sccs_population$days_vax2), as.numeric(pmin((sccs_population$days_vax2 + vax2_end), sccs_population$study_exit_days, sccs_population$days_vax3, na.rm = T)), NA_real_)
  # plus 1 to account for inclusive model fit
  sccs_population$risk_d2_length <- as.numeric((sccs_population$risk_d2_end - sccs_population$risk_d2)) + 1
  
  # time between doses, starts day 1 after dose one risk end or dose 2 whichever is first, ends at dose2 
  # na.rm is F because we want this to be missing if days_vax2 is ever missing 
  sccs_population$between_start <- as.numeric(pmin(sccs_population$risk_d1_end + 1), sccs_population$days_vax2)
  sccs_population$between_start <- ifelse(!is.na(sccs_population$days_vax2), as.numeric(sccs_population$between_start), NA_real_)
  sccs_population$between_end <- as.numeric(pmin(sccs_population$days_vax2, sccs_population$study_exit_days))
  # plus 1 to account for inclusive model fit
  sccs_population$between_length <- as.numeric((sccs_population$between_end - (sccs_population$between_start))) + 1
  
  # control time  
  sccs_population$c_start <- as.numeric(sccs_population$days_vax1 + preexp_start + control_duration)
  sccs_population$c_end <- as.numeric(sccs_population$days_vax1 + preexp_start)
  
  # length of study and control period 
  sccs_population$c_length <- as.numeric(sccs_population$c_end - sccs_population$c_start) + 1 
  sccs_population$study_length <- as.numeric(pmin(sccs_population$study_exit_days, sccs_population$days_last_vax + vax2_end, sccs_population$days_vax3, na.rm = T)) + 1
  
  # start at the first of control window or risk window (depending on design)
  sccs_population$ref_start <- as.numeric(sccs_population$c_start)
  
  # end at the last of the end of the control window, first and second risk windows (depending on design)
  sccs_population$ref_end <- as.numeric(pmax(sccs_population$risk_d1_end, sccs_population$risk_d2_end,  na.rm = T))
  sccs_population$ref_end <- as.numeric(sccs_population$ref_end)
  
  ## calendar time adjustment 
  # create a calendar time variable for calendar time adjustment (30-day interval between start and end in the data where model is fit)
  max_day = max(as.numeric(sccs_population$study_exit_days), na.rm = T)
  min_day = min(as.numeric(sccs_population$start_study_days), na.rm = T)
  
  # generate a cut-off for each month to allow us to summarise events per month
  caltime_groups <- seq(from = min_day-1, to = max_day, by = 30)
  
  # first cut-off is start of second age group, so create the same minus first element
  caltime_cutoffs_for_model <- caltime_groups[-1]
  
  # if outside of study period remove 
  max_ref_day = max(as.numeric(sccs_population$ref_start), na.rm = T)
  min_ref_day = min(as.numeric(sccs_population$ref_end), na.rm = T)
  
  caltime_groups <- subset(caltime_groups, caltime_groups>min_ref_day)
  caltime_groups <- subset(caltime_groups, caltime_groups<max_ref_day)
  caltime_cutoffs_for_model <- subset(caltime_cutoffs_for_model, caltime_cutoffs_for_model>min_ref_day)
  caltime_cutoffs_for_model <- subset(caltime_cutoffs_for_model, caltime_cutoffs_for_model<max_ref_day)
  
  # some groups may not have events. Find out which ones 
  outcome_sums <- numeric(length(caltime_groups) - 1)
  
  for (i in 1:(length(outcome_sums))) {
    outcome_sums[i] <- sum(between(sccs_population$outcome_days, caltime_groups[i], caltime_groups[i+1]), na.rm = T)
  } 
  
  ref_max <- as.numeric(which.max(outcome_sums))
  
# 2: EVALUATE DIFFERENT MODELS 

  # UNADJUSTED  
  
    # standard option from package 
    standardsccs(
      formula = event ~ risk_d1, 
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
      aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = sccs_population 
    )
    
    # edited option 
    standardsccs2(
      formula = event ~ risk_d1, 
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
      aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = sccs_population 
    )

  # ADJUSTED  
    
    
    

