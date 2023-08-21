# 1. FUNCTIONS ---------------------------------------------------------------
# read functions required by the pipeline 

# UTILITY FUNCTIONS

# 1. Function to summarise a vector of binary variables
# Returns a tibble of counts when the variable == 1, and percentages 

indicator_summary <- function(data, ...) {
  
  {{data}} %>% 
    select(all_of(...)) %>% 
    summarise(
      across(
        .cols  = all_of(...),
        list(count = sum, percent = mean), 
        na.rm = TRUE, 
        .names = "{.fn}-{.col}")) %>% 
    mutate(across(contains("percent"), ~round(.x*100, digits = 2))) %>% 
    pivot_longer(everything(), 
                 names_to = c(".value", "variable"), 
                 names_sep = "-") 
  
}

# 2. Function to summarise a vector of continuous variables 
# Returns a tibble of median and IQR (called counts and percentages to enable binding)

cont_summary <- function(data, ...) {
  
  {{data}} %>% 
    select(all_of(...)) %>% 
    summarise(
      across(
        .cols  = all_of(...),
        .fns = list(median = ~median(., na.rm = T),
                    q1 = ~quantile(.,probs =0.25, na.rm = T), 
                    q3 = ~quantile(.,probs =0.75, na.rm = T), 
                    min = ~ifelse(length(.[!is.na(.)]) == 0, NA_real_, min(., na.rm = TRUE)), 
                    max = ~ifelse(length(.[!is.na(.)]) == 0, NA_real_, max(., na.rm = TRUE))), 
        .names = "{.fn}-{.col}")) %>% 
    mutate(across(everything(), ~as.numeric(.x))) %>% 
    mutate(across(everything(), ~round(.x, digits = 0))) %>% 
    pivot_longer(everything(), 
                 names_to = c(".value", "variable"), 
                 names_sep = "-") 
  
}

# 3. Wrap the SCCS calls used 
# wrap the sccs functions in purrrs safely (ensures each call returns results and errors)

safe_sccs2 <- safely(standardsccs2)
safe_eventdepenexp <- safely(eventdepenexp)

# PROJECT SPECIFIC FUNCTIONS 

## TABULATE EVENT-DEPENDENCY INFORMATION 
#' @description Takes clean data and tabulates basic information of interest for assessing event-dependent exposures 
#' Should be run on the SCCS dataset (i.e., don't restrict control time)
#' @param data = Input dataframe, should be loaded into R before invoking, as function does not do this. 
#' @returns A dataframe with basic summaries of the vaccine/event relationship  

event_dependency_table <- function(data) {
  
  vaccine_type <- unique(data$type_vax1)
  message(paste("summarising data for", vaccine_type))
  
  # list the variables of interest
  var_list <- c("vax1", "vax2", "vax3", "dose1_afterevent", "dose2_afterevent", "dose3_afterevent")
  # use a utility function from 00_functions to save the N and % of each dose 
  vaccine_summary <- indicator_summary(data, var_list)
  message("Note: summarised categorical variables succesfully")
  
  # list continous variables of interest
  cont_list <- c("diff_d1d2", "diff_d2d3", "diff_event_d1", "diff_event_d2", "diff_event_d3")
  # use a utility function from 00_functions to save the median and IQR for each continous variable 
  vaccine_duration_summary <- cont_summary(data, cont_list) 
  message("Note: summarised continous variables succesfully")
  
  # Make this into a table 
  vaccine_table <- as.data.frame(bind_rows(vaccine_summary, vaccine_duration_summary))
  names(vaccine_table) <- c("Variable", "N", "%", "Median", "Q1", "Q3", "Min", "Max")
  message("Note: created vaccine table succesfully")
  
 vaccine_table <- vaccine_table %>% 
  # mutate(across(everything(), replace_na, replace = "")) %>% 
    mutate(across(everything(), ~ifelse(.x == "Inf", "", .x))) %>%   
    mutate(across(everything(), ~ifelse(.x == "-Inf", "", .x)))
  
  message("Note: cleaned vaccine table succesfully")
  
  return(vaccine_table)  
  
} 

## PLOT EVENT-DEPENDENCY INFORMATION 
#' @description generate exposure-centered interval plot
#' @param data = Input dataframe, should be loaded into R before invoking, as function does not do this. 
#' @param dose = Dose Number, 1, 2 etc. 
#' @returns ggplot object

event_dependency_plot <- function(data, dose) {
  
  xvar <- paste0("diff_event_d", dose)
  vaxvar <- paste0("date_vax", dose)
  vaccine_type <- unique(data$type_vax1)
  message(paste("plotting histogram of", vaccine_type, xvar))
  
  if (sum(!is.na(data[[vaxvar]])) < 1) { 
    
    message(paste(sum(!is.na(data[[vaxvar]])), "people with ", vaccine_type, "dose", dose, ", plot not created"))
    
  } else {
    
  plot <- ggplot(data, aes(x=eval(parse(text = xvar)))) + 
    geom_histogram(colour = "lightslategrey", fill = "lightsteelblue", na.rm = T) + 
    theme_minimal() + 
    labs(y = "Number of people", 
         x = "Difference (days)", 
         title = paste("Exposure-Centered Interval Plot", "Dose", dose, "-", vaccine_type),  
         subtitle = "Distribution of time between the event and vaccine dose") + 
    geom_vline(xintercept=0, colour="darkslategrey") 
  
  return(plot)
  
  } 
  
} 

#' SCRI period description
#' @description tabulate number of events and duration of each period
#' @param data = input dataframe
#' @returns dataframe 

describe_scri_periods <- function(data) {
  
  # number of events in risk windows 
  n_risk1 <- sum(as.numeric((data$days_vax1 < data$outcome_days) & (data$outcome_days <= data$risk_d1_end)), na.rm = T) 
  n_risk2 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  
  # number of events in pre-exposure period 
  n_preexp <- sum(as.numeric((data$preexp_days < data$outcome_days) & (data$outcome_days <= data$preexp_end)), na.rm = T) 
  
  # number of events in "control" period 
  n_control <- sum(as.numeric((data$c_start < data$outcome_days) & (data$outcome_days <= data$c_end)), na.rm = T) 
  
  # number of events in "between" period 
  n_between <- sum(as.numeric((data$between_start < data$outcome_days) & (data$outcome_days <= data$between_end)), na.rm = T) 
  
  # describe total number of events in each period 
  n_total <- n_risk1 + n_risk2 + n_preexp + n_control + n_between
  
  # Duration
  duration_vars <- c("study_length", "c_length", "preexp_length", "risk_d1_length", "between_length", "risk_d2_length")
  # note, this will throw errors if applied to an empty vector even though it is expected behavior. Warning can be ignored. 
  duration_summary <- cont_summary(data, duration_vars)
  duration_summary <- duration_summary[,-1]
  
  count_summary <- rbind(n_total, n_control, n_preexp, n_risk1, n_between, n_risk2) 
  row_headings <- rbind("Overall", "Control", "Pre-Exposure", "Dose 1 Risk Period", "In between doses", "Dose 2 Risk Period") 
  
  event_summary <- cbind(row_headings,count_summary,duration_summary)
  
  # return the edited data frame 
  return(event_summary)
  
}   

#' SCCS period description
#' @description tabulate number of events and duration of each period
#' @param data = input dataframe
#' @returns dataframe 

describe_sccs_periods <- function(data) {
  
  # describe total number of events in each period 
  n_total <- nrow(data)
  
  # number of events in risk windows 
  n_risk1 <- sum(as.numeric((data$days_vax1 < data$outcome_days) & (data$outcome_days <= data$risk_d1_end)), na.rm = T) 
  n_risk2 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  n_risk3 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  
  # number of events in pre-exposure period 
  n_preexp1 <- sum(as.numeric((data$preexp1_days < data$outcome_days) & (data$outcome_days <= data$preexp1_end)), na.rm = T) 
  n_preexp2 <- sum(as.numeric((data$preexp2_days < data$outcome_days) & (data$outcome_days <= data$preexp2_end)), na.rm = T) 
  n_preexp3 <- sum(as.numeric((data$preexp3_days < data$outcome_days) & (data$outcome_days <= data$preexp3_end)), na.rm = T) 
  
  # Duration
  duration_vars <- c("study_length", "preexp1_length", "risk_d1_length", "preexp2_length", "risk_d2_length", "preexp3_length", "risk_d3_length")
  duration_summary <- cont_summary(data, duration_vars)
  duration_summary <- duration_summary[,-1]
  
  count_summary <- rbind(n_total, n_preexp1, n_risk1, n_preexp2, n_risk2, n_preexp3, n_risk3) 
  row_headings <- rbind("Overall", "Pre-Exposure 1", "Dose 1 Risk Period", "Pre-Exposure 2", "Dose 2 Risk Period", "Pre-Exposure 3", "Dose 3 Risk Period") 
  
  event_summary <- cbind(row_headings,count_summary,duration_summary)
  
  # return the edited data frame 
  return(event_summary)
  
}   

#' Run an unadjusted SCRI
#' @description fit an SCRI using the SCCS package and format the output 
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_unadjusted_scri <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting unadjusted SCRI for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {

    output.1 <- safe_sccs2(
      formula = event ~ risk_d1, 
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
      aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = data 
    )
    
    output.1.result <- output.1$result
    
    if(!is.null(output.1.result)) {
      
      message("unadjusted model generated output")
      
      summary.1 <- as.data.frame(cbind(output.1.result$n, 
                                       output.1.result$nevent, 
                                       output.1.result$conf.int, 
                                       output.1.result$coefficients[,c(1,3)]))
      summary.1$var <- rownames(summary.1)
      summary.1$analysis <- "unadjusted"
      
      rownames(summary.1) <- NULL
      
      summary.1 <- summary.1 %>%
        rename(irr = `exp(coef)`,
               lci = `lower .95`,
               uci = `upper .95`, 
               n = `V2`) %>%
        mutate(study_design = design) %>% 
        select(study_design, analysis, n, var, irr, lci, uci) %>% 
        mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                               var == "risk_d12" ~ "dose 1 risk window", 
                               var == "risk_d13" ~ "between doses", 
                               var == "risk_d14" ~ "dose 2 risk window", 
                               TRUE ~ var)) 
      
    } else {
      
      summary.1 <- data.frame(study_design = design, 
                              analysis = "unadjusted", 
                              n = "n/a", 
                              var = "empty output", 
                              irr = "-", 
                              lci = "-", 
                              uci = "-")
      
      message("empty summary output table for unadjusted analyses created")  
      
    }
  
  output.list <- list(table = summary.1, 
                      error = output.1$error) 
  
  # return object 
  return(output.list)
  } 
}

#' Run an adjusted SCRI
#' @description fit an adjusted SCRI using the SCCS package and format the output 
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_adjusted_scri <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("Fitting adjusted SCRI for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {
    
    message(paste("using ref_max = ", ref_max, "as the reference group for the age adjustment"))
    
    output.2 <- safe_sccs2(
      formula = event ~ risk_d1 + relevel(age, ref = ref_max), 
      indiv = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
      aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
      aevent = outcome_days, 
      agegrp = caltime_cutoffs_for_model, 
      dataformat = "multi",
      sameexpopar = F,  
      data = data
    )
    
    output.2.result <- output.2$result
  
    # adjusted output 
    
    if(!is.null(output.2.result)) {
      
      message("adjusted model generated output")
      
      summary.2 <- as.data.frame(cbind(output.2.result$n, 
                                       output.2.result$nevent, 
                                       output.2.result$conf.int, 
                                       output.2.result$coefficients[,c(1,3)]))
      summary.2$var <- rownames(summary.2)
      summary.2$analysis <- "adjusted"
      rownames(summary.2) <- NULL
      
      summary.2 <- summary.2 %>%
        rename(irr = `exp(coef)`,
               lci = `lower .95`,
               uci = `upper .95`, 
               n = `V2`) %>%
        mutate(study_design = design) %>% 
        select(study_design, analysis, n, var, irr, lci, uci) 
      
    } else {
      
      summary.2 <- data.frame(study_design = design, 
                              analysis = "adjusted", 
                              n = "n/a", 
                              var = "empty output", 
                              irr = "-", 
                              lci = "-", 
                              uci = "-")
      
      message("empty summary output table for adjusted analyses created")  
      
    }
    
    # add labels 
      summary.2 <- summary.2 %>%  
      mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                             var == "risk_d12" ~ "dose 1 risk window", 
                             var == "risk_d13" ~ "between doses", 
                             var == "risk_d14" ~ "dose 2 risk window", 
                             TRUE ~ var)) 
    
    output.list <- list(table = summary.2, 
                        error = output.2$error) 
    
    # return object 
    return(output.list)
    
  } 
}

#' Run unadjusted SCCS
#' @description fit an unadjusted standard SCCS
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_unadjusted_standard_sccs <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting unadjusted SCCS for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {

    output.1 <- safe_sccs2(
      formula = event ~ risk_d1, 
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = data 
    )
    
    output.1.result <- output.1$result
    
    if(!is.null(output.1.result)) {
    
      message("unadjusted model generated output")
      
      summary.1 <- as.data.frame(cbind(output.1.result$n, 
                                       output.1.result$nevent, 
                                       output.1.result$conf.int, 
                                       output.1.result$coefficients[,c(1,3)]))
      summary.1$var <- rownames(summary.1)
      summary.1$analysis <- "unadjusted"
      
      rownames(summary.1) <- NULL
      
      summary.1 <- summary.1 %>%
        rename(irr = `exp(coef)`,
               lci = `lower .95`,
               uci = `upper .95`, 
               n = `V2`) %>%
        mutate(study_design = design) %>% 
        select(study_design, analysis, n, var, irr, lci, uci) %>%  
        mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                               var == "risk_d12" ~ "dose 1 risk window", 
                               var == "risk_d13" ~ "dose 2 pre-exposure", 
                               var == "risk_d14" ~ "dose 2 risk window", 
                               var == "risk_d15" ~ "dose 3 pre-exposure", 
                               var == "risk_d16" ~ "dose 3 risk window", 
                               TRUE ~ var)) 
    } else {
      
      summary.1 <- data.frame(study_design = design, 
                              analysis = "unadjusted", 
                              n = "n/a", 
                              var = "did not fit", 
                              irr = "-", 
                              lci = "-", 
                              uci = "-")
      
      message("empty summary output table for unadjusted analyses created")  
    
    } 
  
  output.list <- list(table = summary.1, 
                      error = output.1$error) 
  
  # return object 
  return(output.list)
  
  } 
}

#' Run adjusted SCCS
#' @description fit an adjusted standard SCCS
#' requires input parameters defining start and end of study period in the dataset
#' uses a wrapped and amended version of standardsccs from the SCCS package 
#' amendment ensures that a userspecified reference level can be used 
#' @param data = input dataframe
#' @returns dataframe 

fit_adjusted_standard_sccs <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting adjusted SCCS for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {
    # specify the reference level 
    # note 
    message(paste("using ref_max = ", ref_max, "as the reference group for the age adjustment"))
      
      output.2 <- safe_sccs2(
        formula = event ~ risk_d1 + relevel(age, ref = ref_max), 
        indiv = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
        aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
        aevent = outcome_days, 
        agegrp = caltime_cutoffs_for_model, 
        dataformat = "multi",
        sameexpopar = F,  
        data = data
      )
      
    output.2.result <- output.2$result
    
    if(!is.null(output.2.result)) {
        
        message("adjusted model generated output")
    
        summary.2 <- as.data.frame(cbind(output.2.result$n, 
                                         output.2.result$nevent, 
                                         output.2.result$conf.int, 
                                         output.2.result$coefficients[,c(1,3)]))
        summary.2$var <- rownames(summary.2)
        summary.2$analysis <- "adjusted"
        rownames(summary.2) <- NULL
          
        summary.2 <- summary.2 %>%
          rename(irr = `exp(coef)`,
                 lci = `lower .95`,
                 uci = `upper .95`, 
                 n = `V2`) %>%
          mutate(study_design = design) %>% 
          select(study_design, analysis, n, var, irr, lci, uci) %>%  
          mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                                 var == "risk_d12" ~ "dose 1 risk window", 
                                 var == "risk_d13" ~ "dose 2 pre-exposure", 
                                 var == "risk_d14" ~ "dose 2 risk window", 
                                 var == "risk_d15" ~ "dose 3 pre-exposure", 
                                 var == "risk_d16" ~ "dose 3 risk window", 
                                 TRUE ~ var)) 
        
      } else {
        
        summary.2 <- data.frame(study_design = design, 
                                analysis = "adjusted", 
                                n = "n/a", 
                                var = "did not fit", 
                                irr = "-", 
                                lci = "-", 
                                uci = "-")
          
        message("empty summary output table for adjusted analyses created")  
        
      } 
      
    output.list <- list(table = summary.2, 
                        error = output.2$error) 
      
    # return object 
    return(output.list)
    
  } 
}

#' Run unadjusted extended SCCS
#' @description fit an unadjusted extended SCCS
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_unadjusted_extended_sccs <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting unadjusted extended SCCS for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {
  
    output.1 <- safe_eventdepenexp(
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = data 
    )
    
  output.1.result <- output.1$result
  
  if(!is.null(output.1.result)) {
    
    summary.1 <- as.data.frame(cbind(output.1.result$n, 
                                     output.1.result$nevent, 
                                     output.1.result$conf.int, 
                                     output.1.result$coefficients[,c(1,3)]))
    summary.1$var <- rownames(summary.1)
    summary.1$analysis <- "unadjusted"
    
    rownames(summary.1) <- NULL
    
    summary.1 <- summary.1 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci) %>% 
      mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                             var == "risk_d12" ~ "dose 1 risk window", 
                             var == "risk_d13" ~ "dose 2 pre-exposure", 
                             var == "risk_d14" ~ "dose 2 risk window", 
                             var == "risk_d15" ~ "dose 3 pre-exposure", 
                             var == "risk_d16" ~ "dose 3 risk window", 
                             TRUE ~ var)) 
  } else {
    
    message("unadjusted model generated output")
    
    summary.1 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not generate output", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for unadjusted analyses created")  
  }
  
    output.list <- list(table = summary.1, 
                        error = output.1$error) 
    
    # return object 
    return(output.list)
    
  
  } 
}

#' Run adjusted extended SCCS
#' @description fit an adjusted extended SCCS
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_adjusted_extended_sccs <- function(data) {
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting adjusted extended SCCS for", vaccine_type))
  
  # stop execution if <5 individuals 
  if (nrow(data) < 5) {
    
    message(paste("<5 people with", vaccine_type,", execution of script halted"))
    
  } else {
     
      output.2 <- safe_eventdepenexp(
        indiv = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
        aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
        aevent = outcome_days, 
        agegrp = caltime_cutoffs_for_model, 
        dataformat = "multi",
        sameexpopar = F,  
        data = data
      )
      
    output.2.result <- output.2$result
      
    if(!is.null(output.2.result)) {
      
      summary.2 <- as.data.frame(cbind(output.2.result$n, 
                                       output.2.result$nevent, 
                                       output.2.result$conf.int, 
                                       output.2.result$coefficients[,c(1,3)]))
      summary.2$var <- rownames(summary.2)
      summary.2$analysis <- "adjusted"
      rownames(summary.2) <- NULL
      
      summary.2 <- summary.2 %>%
        rename(irr = `exp(coef)`,
               lci = `lower .95`,
               uci = `upper .95`, 
               n = `V2`) %>%
        mutate(study_design = design) %>% 
        select(study_design, analysis, n, var, irr, lci, uci) %>%  
        mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                               var == "risk_d12" ~ "dose 1 risk window", 
                               var == "risk_d13" ~ "dose 2 pre-exposure", 
                               var == "risk_d14" ~ "dose 2 risk window", 
                               var == "risk_d15" ~ "dose 3 pre-exposure", 
                               var == "risk_d16" ~ "dose 3 risk window", 
                               TRUE ~ var)) 
      
    } else {
      
      summary.2 <- data.frame(study_design = design, 
                              analysis = "adjusted", 
                              n = "n/a", 
                              var = "did not fit", 
                              irr = "-", 
                              lci = "-", 
                              uci = "-")
      
      message("empty summary output table for adjusted analyses created")  
      
    }
    
    output.list <- list(table = summary.2, 
                        error = output.2$error) 
    
    # return object 
    return(output.list)
    
  }
}

#' Process messages from models
#' @description this function prints a "log" file to a user specified destination
#' the file contains only model messages, warnings and outputs 
#' function designed to be applied using seq_along for easy list subsetting
#' @param index = looping index 
#' @param inputfile = output from model, as a list created by safely_n_quietly
#' @param filename = filename and destination 

process_conditions <- function(index, inputfile, filename) {
  
  i <- index 
  data <- inputfile 
  condition_file_name <- filename
  
  cat(paste(toupper(names(data[i])), "\n\n"), file = condition_file_name, append = TRUE)
  
  cat("messages\n", file = condition_file_name, append = TRUE)
  cat(as.character(data[[i]]$messages), file = condition_file_name, append = TRUE)
  cat("\n", file = condition_file_name, append = TRUE)
  
  cat("warnings\n", file = condition_file_name, append = TRUE)
  cat(as.character(data[[i]]$warnings), file = condition_file_name, append = TRUE)
  cat("\n", file = condition_file_name, append = TRUE)
  
  cat("errors\n", file = condition_file_name, append = TRUE)
  cat(as.character(data[[i]]$result$error), file = condition_file_name, append = TRUE)
  cat("\n", file = condition_file_name, append = TRUE)
  
}

#' Process output from models
#' @description this adds a column flagging if there are errors or warnings for each model
#' function designed to be applied using seq_along for easy list subsetting
#' @param index = looping index 
#' @param inputfile = output from model, as a list created by safely_n_quietly
#' @returns df = data frame with results with status column added 

process_output <- function(index,inputfile) {
  
  i <- index
  data <- inputfile
  
  # create logical tests to check whether results exist with warnings or errors 
  warning_exists <- as.logical(length(data[[i]]$warning)>0)
  error_exists <- as.logical(length(data[[i]]$result$error)>0)
  result_exists <- as.logical(length(data[[i]]$result$table)>0)
  
  # basic formatting where results exist 
  if (result_exists) {
    
    df <- as.data.frame((data[[i]]$result$table)) %>% 
      mutate(model_status = "", 
             Names = names(data[i])) %>% 
      relocate(Names)
    
    # add a note to the csv of result to flag the model status 
    if ((result_exists) && (!warning_exists) && (!error_exists)) {
      
      df <- df %>% 
        mutate(model_status = replace(model_status,1,"no warnings or errors")) 
      
    }
    
    if ((result_exists) && (warning_exists)) {
      
      df <- df %>% 
        mutate(model_status = replace(model_status,1,"warnings generated, check condition file")) 
      
    }
    
    if ((result_exists) && (error_exists)) {
      
      df <- df %>% 
        mutate(model_status = replace(model_status,1,"errors generated, check condition file")) 
      
    }
    
    return(df)
    
  } 
  
} 