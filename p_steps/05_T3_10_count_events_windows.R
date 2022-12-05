# COUNT PERSON TIME PER MYOCARD COHORT
#-----------------------------------------------
# To estimate the incidence rates of risks in periods

# input: D3_myocard_cohort, D3_events_ALL_OUTCOMES.RData ,list_outcomes_observed.RData
# output: D4_persontime_windows (exported to csv)


print("COUNT EVENTS by windows")

source(paste0(dirmacro,"CountPersonTimeV13.8.R"))

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"))
  load(paste0(dirtemp, "D3_study_population_by_window_and_dose", suffix[[subpop]], ".RData"))
  
  events_ALL_OUTCOMES <- get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  study_population <- get(paste0("D3_study_population_by_window_and_dose", suffix[[subpop]]))
  
  rm(list = paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list = paste0("D3_study_population_by_window_and_dose", suffix[[subpop]]))
  
  # missing_OUTCOME_variables <- setdiff(OUTCOME_variables, events_ALL_OUTCOMES[, unique(type_outcome)])
  #if (length(missing_OUTCOME_variables) > 0) {
  vars_to_add <- data.table(person_id = study_population[1, person_id], date = ymd(99991231),
                            type_outcome = c(OUTCOME_variables, CONTROL_variables), meaning_renamed = "DO NOT USE",
                            codvar = "DO NOT USE", event_record_vocabulary = "DO NOT USE")
  events_ALL_OUTCOMES <- rbind(events_ALL_OUTCOMES, vars_to_add)
  #}
  
  max_exit <- study_population[, ceiling_date(max(end_date_of_period), 'year') %m-% days(1)]
  
  not_recurrent_OUTCOME_variables <- setdiff(c(OUTCOME_variables, CONTROL_variables), recurrent_OUTCOME_variables)
  
  print("not recurrent")
  
  persontime_4weeks_not_recurrent <- CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = gsub('-', '', as.character(study_start)),
    End_study_time = gsub('-', '', as.character(max_exit)),
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex", "dose", "type_vax", "period"),
    Name_event = "type_outcome",
    Date_event = "date",
    Age_bands = Agebands_countpersontime,
    Increment = "year",
    Outcomes_nrec = not_recurrent_OUTCOME_variables,
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T
  )
  
  print("recurrent")
  
  persontime_4weeks_recurrent <- CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = gsub('-', '', as.character(study_start)),
    End_study_time = gsub('-', '', as.character(max_exit)),
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex", "dose", "type_vax", "period"),
    Name_event = "type_outcome",
    Date_event = "date",
    Age_bands = Agebands_countpersontime,
    Increment = "year",
    Outcomes_rec = recurrent_OUTCOME_variables, 
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T,
    Rec_period = c(rep(30, length(recurrent_OUTCOME_variables)))
  )
  
  print("Merging")
  
  persontime_4weeks_not_recurrent <- persontime_4weeks_not_recurrent[, .SD,
                                                                     .SDcols = unique(names(
                                                                       persontime_4weeks_not_recurrent))]
  
  persontime_4weeks_recurrent <- persontime_4weeks_recurrent[, .SD,
                                                             .SDcols = unique(names(persontime_4weeks_recurrent))]
  
  persontime_4weeks <- merge(persontime_4weeks_not_recurrent, persontime_4weeks_recurrent,
                              by = c("sex", "dose", "type_vax", "period", "year", "Ageband", "Persontime"))
  
  persontime_4weeks[, (colnames(persontime_4weeks)[grepl("^Persontime", colnames(persontime_4weeks))]) := NULL]

  print("Saving")
  
  nameoutput<-paste0("D4_count_events_windows", suffix[[subpop]])
  assign(nameoutput, persontime_4weeks)
  save(nameoutput, file  =paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
}
