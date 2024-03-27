##%######################################################%##
#                                                          #
####         SELECT PREGNANT POPULATION INSIDE          ####
####            THE STUDY BEFORE COVID START            ####
#                                                          #
##%######################################################%##

for (subpop in subpopulations_non_empty){
  print(subpop)
  
  # Create flowchart for adults and save D4_study_population
  load(paste0(diroutput, "D4_study_population", suffix[[subpop]], ".RData"))
  selected_population <- get(paste0("D4_study_population", suffix[[subpop]]))
  rm(list = paste0("D4_study_population", suffix[[subpop]]))
  
  # Load pregnancies
  pregnancies_df <- as.data.table(get(load(paste0(dirpregnancy, "D3_pregnancy_final.RData"))[[1]]))
  pregnancies_df <- pregnancies_df[, .(pregnancy_id , person_id, pregnancy_start_date, pregnancy_end_date)]
  pregnancies_df <- pregnancies_df[pregnancy_start_date >= start_lookback & pregnancy_end_date < start_COVID_diagnosis_date, ]
  
  # Keep only pregnancy overlapping the study period and before start of covid
  setkey(selected_population, person_id, spell_start_date, study_exit_date)
  setkey(pregnancies_df, person_id, pregnancy_start_date, pregnancy_end_date)
  selected_pregnancies <- data.table::foverlaps(pregnancies_df, selected_population, type = "within", nomatch = NULL)
  selected_pregnancies <- selected_pregnancies[!is.na(spell_start_date), ]
  
  # Define trimesters
  selected_pregnancies[, start_preg_period_during_1 := pregnancy_start_date]
  selected_pregnancies[, end_preg_period_during_1 := pmin(pregnancy_start_date %m+% days(97), pregnancy_end_date)]
  selected_pregnancies[, start_preg_period_during_2 := pregnancy_start_date %m+% days(98)]
  selected_pregnancies[start_preg_period_during_2 > pregnancy_end_date, start_preg_period_during_2 := NA]
  selected_pregnancies[, end_preg_period_during_2 := pmin(pregnancy_start_date %m+% days(195), pregnancy_end_date)]
  selected_pregnancies[, start_preg_period_during_3 := pregnancy_start_date %m+% days(196)]
  selected_pregnancies[start_preg_period_during_3 > pregnancy_end_date, start_preg_period_during_3 := NA]
  selected_pregnancies[, end_preg_period_during_3 := pregnancy_end_date]
  
  # Transform the dataset to long format
  colA = paste("start_preg_period_during", 1:3, sep = "_")
  colB = paste("end_preg_period_during", 1:3, sep = "_")
  selected_pregnancies = data.table::melt(selected_pregnancies, measure = list(colA, colB),
                                          value.name = c("period_entry_date", "period_exit_date"),
                                          variable.name = "pregnancy_period")
  selected_pregnancies[, pregnancy_period := as.character(pregnancy_period)]
  selected_pregnancies[.(pregnancy_period = as.character(seq_along(colA)),
                         to = c("first", "second", "third")), on = "pregnancy_period", pregnancy_period := i.to]
  
  # Clean dataset
  selected_pregnancies[, c("pregnancy_start_date", "pregnancy_end_date") := NULL]
  selected_pregnancies <- selected_pregnancies[!is.na(period_entry_date), ]
  
  nameoutput <- paste0("D4_pregnancy_study_population", suffix[[subpop]])
  assign(nameoutput, selected_pregnancies)
  save(nameoutput, file = paste0(diroutput, nameoutput, ".RData"), list = nameoutput)
  rm(list = nameoutput)
}
