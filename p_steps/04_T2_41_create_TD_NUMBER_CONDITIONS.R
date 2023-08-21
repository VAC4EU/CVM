#-------------------------------
# create the dataset with number of conditiosn

# input
# D3_TD_variable_condition (condition is a parameter, values: DIAB CANCER PULMON OBES CKD HIV IMMUNOSUP SICKLE CVD)

# output
# D3_TD_NUMBER_CONDITIONS 

# use variable_definition[[element]] to retrieve all the conceptsets associated to one of the elements; see 3_30 for a similar thing

for (subpop in subpopulations_non_empty) {
  
  print(subpop)
  # 
  # cov_ALL <- data.table()
  
  # Load all TD datasets for covariates
  D3_TD_list <- lapply(list_of_covariates_for_cohort,
                       function(x) {
                         print(paste("Loading", x))
                         temp <- get(load(paste0(dirTD, "D3_TD_variable_", x,
                                                   suffix[[subpop]], ".RData"))[[1]])
                         
                         # When person enter the study, can start from 1 by not changing value_of_variable to -1
                         temp[, min_date := min(date), by = "person_id"]
                         temp[, min_date := fifelse(min_date == date, T, F)]
                         temp[!min_date & value_of_variable == 0, value_of_variable := -1]
                         temp[, min_date := NULL]
                         setnames(temp, "value_of_variable", "contribution")
                         })
  
  D3_TD_list <- rbindlist(D3_TD_list)
  
  # Sum contribution/value_of_variable in each date for each person
  D3_TD_list <- D3_TD_list[, .(total_contribution = sum(contribution)), by = c("person_id", "date")]
  
  # Order by person and date the take the cumulative sum
  setorder(D3_TD_list, "person_id", "date")
  D3_TD_list <- D3_TD_list[, .(date, cumulative = cumsum(total_contribution)), by = "person_id"]
  setnames(D3_TD_list, "cumulative", "value_of_variable")
  
  # Save the dataset
  name_export_df <- paste0("D3_TD_variable_NUMBER_CONDITIONS", suffix[[subpop]])
  assign(name_export_df, D3_TD_list)
  save(name_export_df, file = paste0(dirTD, "/", name_export_df, ".RData"), list = name_export_df)
  rm(list = name_export_df)
  
  rm(D3_TD_list)
}

