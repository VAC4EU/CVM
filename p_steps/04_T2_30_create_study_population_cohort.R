#-------------------------------
# create the input population for the cohort study

# input
# D3_study_population_SCRI, D3_PERSONS

# output
# D3_study_population_cohort


load(paste0(dirtemp, "D3_PERSONS.RData"))

for (subpop in subpopulations_non_empty){
  print(subpop)
  
  load(paste0(dirtemp,"D3_study_population_SCRI", suffix[[subpop]], ".RData"))
 
  study_population <- get(paste0("D3_study_population_SCRI", suffix[[subpop]]))
  
  study_population <- study_population[,.(person_id,sex,age_at_study_entry,study_entry_date,study_exit_date,datasource,date_vax1,date_of_death)]
  
  
  study_population <- merge(study_population,D3_PERSONS,by = "person_id",all.x = TRUE)
  
  study_population <- study_population[,.(person_id,sex,birth_date,age_at_study_entry,study_entry_date,study_exit_date,datasource,date_vax1,date_of_death)]
  
  nameoutput <- paste0("D3_study_population_cohort", suffix[[subpop]])
  assign(nameoutput, study_population)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
  rm(list = nameoutput)
}
