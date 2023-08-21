#------------------------------------------------------------------
# create events and components of OUTCOMES 

# input: concept set datasets of outcomes (narrow and possible), D4_study_population
# output: for each outcome OUTCOME, D3_components_OUTCOME.RData and D3_events_OUTCOME_type.RData, for type = narrow, possible

### STEP SETUP
# Vocabulary for value of new variable type_of_date and old variable which contains dates
vocabulary_covariate_times <- vector()
vocabulary_covariate_times[["baseline"]] <- "study_entry_date"
vocabulary_covariate_times[["vax1"]] <- "date_vax1"
# vocabulary_covariate_times[["first_covid"]] <- "b"

# final names for covariates and drug proxies
final_COV_variables <- paste0(COV_variables, "_in_365gg")
final_DP_variables <- paste0(DP_variables, "_in_90gg")
covid_variables <- "COVID"
pregnancy_variables <- "PREGNANCY"

# Create an empty covariate dataset
empty_covariate_df <- data.table(person_id = character(),
                                 date = as.Date(as.POSIXct(character())),
                                 type_covariate = character(),
                                 meaning_renamed = character(),
                                 codvar = character(),
                                 type_of_date = character())

# Same for the pregnancy one
empty_pregnancy_df <- data.table(person_id = character(),
                                 pregnancy_start_date = as.Date(as.POSIXct(character())),
                                 pregnancy_end_date = as.Date(as.POSIXct(character())),
                                 type_covariate = character(),
                                 meaning_renamed = character(),
                                 codvar = character(),
                                 type_of_date = character())

# Create an empty D3_covariates_ALL with the correct column types
base_matrix <- matrix(nrow = 0, ncol = length(c(final_COV_variables, final_DP_variables,
                                                covid_variables, pregnancy_variables)) + 2)
empty_covariates_ALL <- setnames(data.table(base_matrix), c("person_id", "type_of_date", final_COV_variables,
                                                            final_DP_variables, covid_variables, pregnancy_variables))
empty_covariates_ALL[, person_id := as.character(person_id)][, type_of_date := as.character(type_of_date)]
empty_covariates_ALL <- empty_covariates_ALL[, (c(final_COV_variables, final_DP_variables,
                                                  covid_variables, pregnancy_variables)) := lapply(.SD,function(x) as.integer(x)), .SDcols = sapply(empty_covariates_ALL, is.logical)]

print('create events and create components of COVARIATES')

### STEP START
for (subpop in subpopulations_non_empty) {
  
  print(subpop)
  
  cov_ALL <- data.table()
  
  # Import the study population
  name_D4_study_population <- paste0("D4_study_population", suffix[[subpop]])
  load(paste0(diroutput, name_D4_study_population, ".RData")) 
  study_population <- get(name_D4_study_population)[, .(person_id, study_entry_date)]
  
  # Create another empty D3_covariates_ALL with all the combinations of covariate_time and person_id
  empty_row_covariates_ALL <- as.data.table(expand.grid(person_id = unlist(study_population[, .(person_id)]),
                                                        type_of_date = names(vocabulary_covariate_times)))
  empty_row_covariates_ALL <- rbind(empty_covariates_ALL, empty_row_covariates_ALL, fill = T)
  
  # Import doses dataset and retain only date of first vaccination
  load(paste0(dirtemp,"D3_vaccines_curated.RData"))
  D3_vaccines_curated <- D3_vaccines_curated[dose_curated == 1, ][, .(person_id, date_curated)]
  setnames(D3_vaccines_curated, "date_curated", "date_vax1")
  
  # Merge population and first vaccinations. Left join since we want all the persons in the study population
  study_population <- merge(study_population, D3_vaccines_curated, all.x = T, by = "person_id")
  
  for (COVARIATE in c(COV_variables, DP_variables, covid_variables, pregnancy_variables)) {
    print(paste("NOW CALCULATING THE COVARIATE:", COVARIATE))
    
    # Create the list of dataset generated after createconceptsetdataset
    nameconceptset <- variable_definition[[COVARIATE]]
    nameconceptset <- nameconceptset[nameconceptset %in% sub('\\.RData$', '', list.files(dirconceptsets))]
    
    if (COVARIATE %in% covid_variables) {
      nameconceptset <- paste0("D3_covid_episodes", suffix[[subpop]])
    } else if (COVARIATE %in% pregnancy_variables) {
      if (skip_pregnancy) {
        nameconceptset <- c()
      } else {
        nameconceptset <- "D3_pregnancy_final"
      }
    }
    
    # # Load, then select/create only variable of interest. If original concept empty use a base df.
    # if (COVARIATE %in% covid_variables) {
    #   conceptsets_list <- lapply(c("D3_covid_episodes"), function(x) {
    #     df <- get(load(paste0(dirtemp, x, suffix[[subpop]], ".RData"))[[1]])
    #     if(nrow(df) == 0) {df <- empty_covariate_df} else {df <- df[, .(person_id, date, meaning_renamed = "covid",
    #                                                                     codvar = "covid")]}
    #     return(df)
    #   })
    # } else if (COVARIATE %in% pregnancy_variables) {
    #   conceptsets_list <- lapply(c("D3_pregnancy_final"), function(x) {
    #     if (skip_pregnancy) {
    #       df <- empty_pregnancy_df
    #     } else {
    #       df <- as.data.table(get(load(paste0(dirpregnancy, x, ".RData"))[[1]]))
    #       if(nrow(df) == 0) {df <- empty_pregnancy_df} else {
    #         df <- df[, .(person_id, pregnancy_start_date, pregnancy_end_date,
    #                      meaning_renamed = meaning_of_principal_record, codvar = "pregnancy")]}
    #     }
    #     return(df)
    #   })
    # } else {
    #   conceptsets_list <- lapply(nameconceptset, function(x) {
    #     df <- get(load(paste0(dirconceptsets, x, ".RData"))[[1]])
    #     if(nrow(df) == 0) {df <- empty_covariate_df} else {df <- df[, .(person_id, date, meaning_renamed, codvar)]}
    #     return(df)
    #   })
    # }
    
    # Create additionalvar parameter for MergeFilterAndCollapse
    additionalvar <- list()
    summarystat <- list()
    
    for (covariate_time in names(vocabulary_covariate_times)) {
      
      # new_var <- c(paste0("binary_at_", covariate_time), paste0("fifelse(", selectionOUTCOME, ", 1, 0)"))
      # additionalvar[[covariate_time]] <- new_var
      # summarystat[[covariate_time]] <- c("max", paste0("binary_at_", covariate_time), paste0("binary_at_", covariate_time))
      # 
      # Merge the concepts, join with study population and the take the maximum date for each person
      
      full_components <- if (!is.null(nameconceptset)) data.table() else data.table(person_id = character(),
                                                                                    date_exists = integer())
      
      for (conceptset in nameconceptset) {
        
        # Add to selectionOUTCOME the correct variable to select
        date_var_name <- vocabulary_covariate_times[[covariate_time]]
        lookback <- if (conceptset %in% DP_variables) 90 else 365
        
        if (COVARIATE %in% covid_variables) {
          selectionOUTCOME <- sprintf("!is.na(date) & !is.na(%s) & date < %s", date_var_name, date_var_name)
        } else if (COVARIATE %in% pregnancy_variables) {
          selectionOUTCOME <- sprintf("!is.na(pregnancy_start_date) & !is.na(pregnancy_end_date) & !is.na(%s) &
                                    pregnancy_start_date <= %s & pregnancy_end_date >= %s",
                                      date_var_name, date_var_name, date_var_name)
        } else {
          selectionOUTCOME <- sprintf("!is.na(date) & !is.na(%s) & date >= %s - %s & date < %s",
                                      date_var_name, date_var_name, lookback, date_var_name)
          # delete records that are not observed in this whole subpopulation
          if (this_datasource_has_subpopulations){
            selectionOUTCOME <- paste0(selectionOUTCOME, ' & ', select_in_subpopulationsEVENTS[[subpop]])
          }
        }
        
        if (COVARIATE %in% covid_variables) {
          conceptset_rdata <- paste0(dirtemp, conceptset, ".RData")
        } else if (COVARIATE %in% pregnancy_variables) {
          conceptset_rdata <- paste0(dirpregnancy, conceptset, ".RData")
        } else {
          conceptset_rdata <- paste0(dirconceptsets, conceptset, ".RData")
        }
        
        # Che existence and recode person_id as character
        components <- MergeFilterAndCollapse(listdatasetL = list(as.data.table(get(load(conceptset_rdata)[[1]])
        )[, person_id := as.character(person_id)]),
        datasetS = study_population,
        key = "person_id",
        condition = selectionOUTCOME,
        strata = c("person_id"),
        summarystat = list(c("exist", "person_id", "date_exists")))
        
        if (COVARIATE %not in% c(covid_variables, pregnancy_variables)) {
          rm(list = conceptset)
        }
        full_components <- list(full_components, components)
        rm(components)
        full_components <- unique(rbindlist(full_components))
        
      }
      
      # add the variable name
      components <- full_components[, type_covariate := COVARIATE][, type_of_date := covariate_time]
      
      # Append to cov_ALL to create a single file with all the covariates
      cov_ALL <- rbind(cov_ALL, components)
    }
  }
  
  cov_ALL <- unique(cov_ALL)
  cov_ALL <- cov_ALL[, type_covariate := paste0(type_covariate, 
                                                fcase(type_covariate %in% DP_variables, "_in_90gg",
                                                      type_covariate %in% COV_variables, "_in_365gg",
                                                      default = ""))]
  
  # Covariates as columns, 1 as values
  cov_ALL <- dcast(cov_ALL, person_id + type_of_date ~ type_covariate, value.var = "date_exists")
  
  # Join with the empty_covariates_ALL the previous dataset to add missing covariates(columns) 
  half_full_covariates_ALL <- rbind(empty_covariates_ALL, cov_ALL, fill = T)
  
  # Update join with dataset with all rows/columns but empty and the one calculated above with all the observations
  full_covariates_ALL <- empty_row_covariates_ALL[half_full_covariates_ALL, on = .(person_id, type_of_date),
                                                  names(half_full_covariates_ALL) := mget(
                                                    paste0("i.", names(half_full_covariates_ALL)))]
  
  for (i in names(full_covariates_ALL)) full_covariates_ALL[is.na(get(i)), (i):=0]
  
  # ADD any_risk_factors as the maximum value between each covariate
  full_covariates_ALL <- full_covariates_ALL[, any_risk_factors := do.call(pmax, c(.SD, list(na.rm = T))),
                                             .SDcols = c(final_COV_variables, final_DP_variables)]
  
  # Save D3_covariates_ALL
  nameobjectCOVARIATES <- paste0("D3_covariates_ALL", suffix[[subpop]])
  assign(nameobjectCOVARIATES, full_covariates_ALL)
  save(nameobjectCOVARIATES, file = paste0(dirtemp, nameobjectCOVARIATES, ".RData"), list = nameobjectCOVARIATES)
}



