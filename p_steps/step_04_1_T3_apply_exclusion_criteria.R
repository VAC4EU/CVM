#---------------------------------------------------------------
# Apply exclusion criteria to create study population 

# input: D3_selection_criteria
# output: Flowchart_basic_exclusion_criteria, Flowchart_exclusion_criteria, D4_study_population, D4_study_source_population

print('FLOWCHART')

#USE THE FUNCTION CREATEFLOWCHART TO SELECT THE SUBJECTS IN POPULATION

for (subpop in subpopulations_non_empty){
  print(subpop)
  load(paste0(dirtemp,"D3_selection_criteria",suffix[[subpop]],".RData"))
  selection_criteria<-get(paste0("D3_selection_criteria", suffix[[subpop]]))


selected_population <- CreateFlowChart(
  dataset = selection_criteria,
  listcriteria = c("sex_or_birth_date_missing", "birth_date_absurd", "no_observation_period", "death_before_study_entry", "no_observation_period_including_study_start"),
  flowchartname = paste0("Flowchart_basic_exclusion_criteria",suffix[[subpop]]))

suppressWarnings(
  if(!(file.exists(direxp))) {dir.create(file.path(direxp))}
)
suppressWarnings(
  if(!(file.exists(diroutput))) {dir.create(file.path(diroutput))}
)


to_study_pop <- unique(selected_population)[, .(person_id, sex, date_of_birth, date_of_death,study_entry_date,start_follow_up, study_exit_date, insufficient_run_in)]

D4_study_source_population <- to_study_pop[, .(person_id, sex, date_of_birth, date_of_death, study_entry_date,start_follow_up, study_exit_date)]

thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])

fwrite(get(paste0("Flowchart_basic_exclusion_criteria",suffix[[subpop]])), paste0(thisdirexp,"Flowchart_basic_exclusion_criteria.csv"))


selected_population_all_filter <- CreateFlowChart(
  dataset = selected_population[,.(person_id, sex, date_of_birth, date_of_death, study_entry_date, start_follow_up, study_exit_date, insufficient_run_in)],
  listcriteria = c("insufficient_run_in"),
  flowchartname = paste0("Flowchart_exclusion_criteria",suffix[[subpop]]))

D4_study_population <- unique(selected_population_all_filter)[, .(person_id, sex, date_of_birth, date_of_death, study_entry_date,
                                                                  start_follow_up, study_exit_date)]


fwrite(get(paste0("Flowchart_exclusion_criteria",suffix[[subpop]])), paste0(thisdirexp,"Flowchart_exclusion_criteria.csv"))

assign(paste0("D4_study_population",suffix[[subpop]]),D4_study_population)
save(list=paste0("D4_study_population",suffix[[subpop]]),file = paste0(diroutput, "D4_study_population",suffix[[subpop]],".RData"))
rm(list=(paste0("D4_study_population",suffix[[subpop]])))

assign(paste0("D4_study_source_population",suffix[[subpop]]),D4_study_source_population)
save(list=paste0("D4_study_source_population",suffix[[subpop]]),file = paste0(diroutput, "D4_study_source_population",suffix[[subpop]],".RData"))
rm(list=(paste0("D4_study_source_population",suffix[[subpop]])))

rm(list=paste0("Flowchart_basic_exclusion_criteria",suffix[[subpop]]))
rm(list=paste0("Flowchart_exclusion_criteria",suffix[[subpop]] ))
rm(list=paste0("D3_selection_criteria", suffix[[subpop]]))
}

if(this_datasource_has_subpopulations==T) rm(D4_study_population,D4_study_source_population)
rm( selection_criteria,selected_population, to_study_pop, selected_population_all_filter)



