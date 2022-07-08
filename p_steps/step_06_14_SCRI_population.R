COVnames <- c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE")

for (subpop in subpopulations_non_empty) {
  
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_outcomes_covid",suffix[[subpop]],".RData")) #L1plus
  load(paste0(dirtemp,"D3_study_population",suffix[[subpop]],".RData"))
  
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  outcomes_covid<-get(paste0("D3_outcomes_covid", suffix[[subpop]]))
  rm(list=paste0("D3_outcomes_covid", suffix[[subpop]]))
  study_population<-get(paste0("D3_study_population", suffix[[subpop]]))
  rm(list=paste0("D3_study_population", suffix[[subpop]]))
  
  SCRI <- study_population[, c("person_id", "sex", "age_at_study_entry", "type_vax_1", "type_vax_2", "study_entry_date",
                               "study_exit_date", "date_vax1", "date_vax2", "immunosuppressants_at_study_entry")]
  rm(study_population)
  
  events_card <- events_ALL_OUTCOMES[name_event %in% c("Myocardalone_narrow", "PERICARD_narrow"), ]
  events_card <- events_card[, c("person_id", "name_event", "date_event")]
  events_card <- events_card[events_card[,.I[which.min(date_event)], by = c("person_id", "name_event")][['V1']]]
  events_card <- dcast(events_card, person_id ~ name_event, value.var = "date_event")
  events_card <- rbind(events_card, data.table(person_id = character(), Myocardalone_narrow = Date(),
                                               PERICARD_narrow = Date()), fill = TRUE)
  setnames(events_card, c("Myocardalone_narrow", "PERICARD_narrow"), c("myocarditis_date", "pericarditis_date"), 
           skip_absent = TRUE)
  SCRI <- merge(SCRI, events_card, all.x = T, by = "person_id")
  rm(events_card, events_ALL_OUTCOMES)
  
  load(paste0(diroutput,"D3_study_population_cov_ALL", suffix[[subpop]],".RData"))
  study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL", suffix[[subpop]]))
  rm(list=paste0("D3_study_population_cov_ALL", suffix[[subpop]]))
  study_population_cov_ALL <- study_population_cov_ALL[, "comorbidity" := fifelse(rowSums(.SD) == 0, 0, 1),
                                                       .SDcols = paste0(COVnames, "_either_DX_or_DP")]
  study_population_cov_ALL <- study_population_cov_ALL[, c("person_id", "comorbidity")]
  SCRI <- merge(SCRI, study_population_cov_ALL, all.x = T, by = "person_id")
  rm(study_population_cov_ALL)
  
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData"))
  study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
  rm(list=paste0("D4_study_population", suffix[[subpop]]))
  study_population <- study_population[, c("person_id", "date_of_death")]
  setnames(study_population, "date_of_death", "death_date")
  SCRI <- merge(SCRI, study_population, all.x = T, by = "person_id")
  rm(study_population)
  
  any_covid <- unique(outcomes_covid[, c("person_id", "date_event")])
  setnames(any_covid, "date_event", "covid_19_date")
  SCRI <- merge(SCRI, any_covid, all.x = T, by = "person_id")
  rm(any_covid, outcomes_covid)
  
  for (col_date in c("myocarditis_date", "pericarditis_date", "death_date", "covid_19_date")) {
    SCRI <- SCRI[get(col_date) > study_exit_date, (col_date) := NA]
  }
  
  SCRI_names <- copy(colnames(SCRI))
  SCRI <- SCRI[, DAP := thisdatasource]
  setcolorder(SCRI, c("DAP", SCRI_names))
  
  tempname<-paste0("D3_study_variables_for_SCRI", suffix[[subpop]])
  assign(tempname, SCRI)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  
  rm(list = paste0("D3_study_variables_for_SCRI", suffix[[subpop]]))
  rm(SCRI)
  
}