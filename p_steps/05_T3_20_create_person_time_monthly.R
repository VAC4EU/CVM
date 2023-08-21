# COUNT PERSON TIME PER COVID PER VACCINATED
#-----------------------------------------------
# To estimate the weekly incidence rates of risks in 2020 by data source for vaccinated cohort

# input: D3_vaxweeks_vaccin_cohort
# output: D4_persontime_risk_month

print("COUNT PERSON TIME by month")

split_by <- if(thisdatasource %in% datasource_needing_split_conceptsets) c("sex", "dose", "type_vax", "COVID19", 1) else NULL

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(dirtemp, "D3_study_population_by_dose", suffix[[subpop]], ".RData"))
  study_population <- get(paste0("D3_study_population_by_dose", suffix[[subpop]]))
  rm(list = paste0("D3_study_population_by_dose", suffix[[subpop]]))
  
  # # missing_OUTCOME_variables <- setdiff(OUTCOME_variables, events_ALL_OUTCOMES[, unique(type_outcome)])
  # #if (length(missing_OUTCOME_variables) > 0) {
  #   vars_to_add <- data.table(person_id = study_population[1, person_id], date = ymd(99991231),
  #                             type_outcome = c(OUTCOME_variables, CONTROL_variables), meaning_renamed = "DO NOT USE",
  #                             codvar = "DO NOT USE", event_record_vocabulary = "DO NOT USE")
  #   events_ALL_OUTCOMES <- rbind(events_ALL_OUTCOMES, vars_to_add)
  # #}
  
  max_exit <- study_population[, ceiling_date(max(end_date_of_period), 'month') %m-% days(1)]
  
  not_recurrent_OUTCOME_variables <- setdiff(c(OUTCOME_variables, CONTROL_variables), recurrent_OUTCOME_variables)
  
  persontime_monthly <- CountPersonTime(
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = gsub('-', '', as.character(study_start)),
    End_study_time = gsub('-', '', as.character(max_exit)),
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex", "dose", "type_vax", "COVID19"),
    Age_bands = Agebands_countpersontime,
    Increment = "month",
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T,
    intermediate_folder = dirtemp,
    split_by = split_by
    )
  
  persontime_monthly <- persontime_monthly[, .SD, .SDcols = unique(names(persontime_monthly))]
  
  print("Saving")
  
  nameoutput <- paste0("D4_persontime_monthly", suffix[[subpop]])
  assign(nameoutput, persontime_monthly)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
}
