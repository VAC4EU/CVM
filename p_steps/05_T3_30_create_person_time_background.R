# COUNT PERSON TIME PER COVID PER VACCINATED
#-----------------------------------------------
# To estimate the weekly incidence rates of risks in 2020 by data source for vaccinated cohort

# input: D3_events_ALL_OUTCOMES, D3_vaxweeks_vaccin_cohort
# output: D4_persontime_risk_month

print("COUNT PERSON TIME for background")

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"))
  load(paste0(dirtemp, "D3_study_population_by_dose", suffix[[subpop]], ".RData"))
  
  events_ALL_OUTCOMES <- get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  study_population <- get(paste0("D3_study_population_by_dose", suffix[[subpop]]))
  
  rm(list = paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list = paste0("D3_study_population_by_dose", suffix[[subpop]]))
  
  study_population <- study_population[dose == 0, ][, c("dose", "type_vax") := NULL]
  
  # missing_OUTCOME_variables <- setdiff(OUTCOME_variables, events_ALL_OUTCOMES[, unique(type_outcome)])
  #if (length(missing_OUTCOME_variables) > 0) {
    vars_to_add <- data.table(person_id = study_population[1, person_id], date = ymd(99991231),
                              type_outcome = c(OUTCOME_variables, CONTROL_variables, "DEATH"), meaning_renamed = "DO NOT USE",
                              codvar = "DO NOT USE", event_record_vocabulary = "DO NOT USE")
    events_ALL_OUTCOMES <- rbind(events_ALL_OUTCOMES, vars_to_add)
  #}
  
  max_exit <- study_population[, ceiling_date(max(end_date_of_period), 'year') %m-% days(1)]
  
  not_recurrent_OUTCOME_variables <- setdiff(c(OUTCOME_variables, CONTROL_variables, "DEATH"), recurrent_OUTCOME_variables)
  
  print("not recurrent")
  
  if (thisdatasource %in% c("UOSL")) {
    split_by = NULL
  } else {
    split_by = c("sex", "COVID19")
  }
  
  persontime_monthly_not_recurrent <- CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = gsub('-', '', as.character(study_start)),
    End_study_time = gsub('-', '', as.character(max_exit)),
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex", "COVID19"),
    Name_event = "type_outcome",
    Date_event = "date",
    Age_bands = Agebands_countpersontime,
    Increment = "year",
    Outcomes_nrec = not_recurrent_OUTCOME_variables,
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T,
    intermediate_folder = dirtemp,
    split_by = split_by
    )
  
  print("recurrent")
  
  persontime_monthly_recurrent <- CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = gsub('-', '', as.character(study_start)),
    End_study_time = gsub('-', '', as.character(max_exit)),
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex", "COVID19"),
    Name_event = "type_outcome",
    Date_event = "date",
    Age_bands = Agebands_countpersontime,
    Increment = "year",
    Outcomes_rec = recurrent_OUTCOME_variables, 
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T,
    Rec_period = c(rep(30, length(recurrent_OUTCOME_variables))),
    intermediate_folder = dirtemp,
    split_by = split_by
  )
  
  persontime_monthly_not_recurrent <- persontime_monthly_not_recurrent[, .SD,
                                                                       .SDcols = unique(names(
                                                                         persontime_monthly_not_recurrent))]
  
  persontime_monthly_recurrent <- persontime_monthly_recurrent[, .SD,
                                                               .SDcols = unique(names(persontime_monthly_recurrent))]
  
  print("Merging")
  
  persontime_monthly <- merge(persontime_monthly_not_recurrent, persontime_monthly_recurrent,
                              by = c("sex", "year", "Ageband", "COVID19", "Persontime"))
  
  print("Saving")
  
  nameoutput <- paste0("D4_persontime_background", suffix[[subpop]])
  assign(nameoutput, persontime_monthly)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
}
