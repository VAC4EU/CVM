# Create D3 for SCRI. 
#-----------------------------------------------
# input: D3_Total_study_population, D3_events_ALL_OUTCOMES, D3_covid_episodes
# output: D3_study_population_SCRI

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  # Load study population, vaccines and covariates
  load(paste0(dirtemp, "D3_Total_study_population", suffix[[subpop]], ".RData"))
  study_population <- get(paste0("D3_Total_study_population", suffix[[subpop]]))
  rm(list = paste0("D3_Total_study_population", suffix[[subpop]]))
  
  # Remove not used columns
  study_population[, spell_start_date := NULL]
  
  # Change names for vaccination and dates
  setnames(study_population, do.call(paste0, expand.grid(c("type_vax_", "date_vax_"), seq_len(max_number_doses))),
           do.call(paste0, expand.grid(c("type_vax", "date_vax"), seq_len(max_number_doses))))
  
  # Calculate age at start follow-up and remove date_of_birth
  study_population[, age_at_study_entry := age_fast(date_of_birth, start_followup_study)]
  study_population[, date_of_birth := NULL]
  
  # Rename start_followup_study to study_entry_date
  setnames(study_population, "start_followup_study", "study_entry_date")
  
  # Create column with name of datasource
  study_population[, datasource := thisdatasource]
  
  # Load the covid episodes
  load(paste0(dirtemp, "D3_covid_episodes", suffix[[subpop]], ".RData"))
  covid_episodes <- get(paste0("D3_covid_episodes", suffix[[subpop]]))
  rm(list = paste0("D3_covid_episodes", suffix[[subpop]]))
  
  covid_episodes <- covid_episodes[date >= start_COVID_diagnosis_date, ]
  
  # Select only the first diagnosis for each person
  covid_episodes <- covid_episodes[covid_episodes[,.I[which.min(date)], by = c("person_id")]$V1]
  
  # Change column name to covid19_date
  setnames(covid_episodes, "date", "covid19_date")
  
  # Merge population and covid episodes
  study_population_covid <- merge(study_population, covid_episodes, all.x = T, by = "person_id")
  
  # Save the study_entry for each persons to be used in selecting the events
  study_population_entry <- study_population[, .(person_id, study_entry_date)]
  
  # Load events and keep only relevant columns
  load(paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"))
  events_ALL_OUTCOMES <- get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list = paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, .(person_id, date, type_outcome)]
  
  # Filter for events of interest
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[type_outcome %in% SCRI_variables_vocabulary$vac4eu]
  
  # Update join name of events for SCRI
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[type_outcome %in% SCRI_variables_vocabulary$vac4eu]
  
  temp_SCRI_variables_vocabulary <- copy(SCRI_variables_vocabulary)
  temp_SCRI_variables_vocabulary[ , scri := paste0(scri, "_date")]
  
  setnames(temp_SCRI_variables_vocabulary, "vac4eu", "type_outcome")
  events_ALL_OUTCOMES[temp_SCRI_variables_vocabulary, on = .(type_outcome), type_outcome := i.scri]
  
  # Merge to get the study_entry_date and filter only events occuring after it. Then remove study_entry_date
  events_ALL_OUTCOMES <- merge(events_ALL_OUTCOMES, study_population_entry, by = "person_id")
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[date >= study_entry_date, ][, study_entry_date := NULL]
  
  # Select only the first event for each person
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[events_ALL_OUTCOMES[,.I[which.min(date)], by = c("person_id", "type_outcome")]$V1]
  
  # From long to wide
  events_ALL_OUTCOMES <- dcast(events_ALL_OUTCOMES, person_id ~ type_outcome, value.var = "date")
  scri_variables <- paste0(SCRI_variables_vocabulary$scri, "_date")
  missing_columns <- setdiff(scri_variables,
                             check_columns_exist(events_ALL_OUTCOMES, scri_variables))
  
  if (length(missing_columns) != 0) {
    events_ALL_OUTCOMES[, (missing_columns) := as.Date(as.POSIXct(character()))]
  }
  
  # Merge population and events
  study_population_covid_events <- merge(study_population_covid, events_ALL_OUTCOMES, all.x = T, by = "person_id")
  
  # Save dataset for weekly countpersontime
  nameobject <- paste0("D3_study_population_SCRI", suffix[[subpop]])
  assign(nameobject, study_population_covid_events)
  save(nameobject, file = paste0(dirtemp, nameobject, ".RData"), list = nameobject)
  rm(list = nameobject)
  
}
