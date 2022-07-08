for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  dataset_name <- paste0("D3_study_population", suffix[[subpop]])
  load(paste0(dirtemp, dataset_name, ".RData"))
  study_population <- get(dataset_name)
  rm(list=dataset_name)
  
  dataset_name <- paste0("D3_outcomes_covid", suffix[[subpop]])
  load(paste0(dirtemp, dataset_name, ".RData"))
  outcomes_covid <- get(dataset_name)
  rm(list = dataset_name)
  
  outcomes_covid <- outcomes_covid[, .(person_id, date_event)]
  outcomes_covid <- outcomes_covid[outcomes_covid[,.I[which.min(date_event)], by = c("person_id")]$V1]
  
  setnames(study_population,
           c("CV_at_study_entry", "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
             "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
             "immunosuppressants_at_study_entry", "at_risk_at_study_entry"),
           c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE",
             "IMMUNOSUPPR", "any_risk_factors"))
  n_doses <- 2
  study_population <- study_population %>%
    select(person_id, sex, study_entry_date, study_exit_date, ageband_at_study_entry, type_vax_1, type_vax_2,
           study_entry_date_vax1, study_exit_date_vax1, study_entry_date_vax2, study_exit_date_vax2, CV, COVCANCER,
           COVCOPD, COVHIV, COVCKD, COVDIAB, COVOBES, COVSICKLE, IMMUNOSUPPR, any_risk_factors)
  
  study_population <- study_population[, study_exit_date := pmin(study_entry_date_vax1 - 1, study_exit_date, na.rm = T)]
  colA = c(paste("study_entry_date_vax", 1:n_doses, sep = ""), "study_entry_date")
  colB = c(paste("study_exit_date_vax", 1:n_doses, sep = ""), "study_exit_date")
  
  colD = paste("type_vax", 1:n_doses, sep = "_")
  study_population <- data.table::melt(study_population, measure = list(colA, colB, colD),
                           variable.name = "Dose", variable.factor = F,
                           value.name = c("study_entry_date", "study_exit_date", "type_vax"),
                           na.rm = F)
  
  study_population <- study_population[!is.na(study_entry_date), ][Dose == as.character(n_doses + 1), Dose := "0"]
  study_population <- merge(study_population, outcomes_covid, all.x = T, by = "person_id")
  
  rm(outcomes_covid)
  
  end_2020 <- ymd(20201231)
  
  study_population <- study_population[Dose == "0", study_exit_date := pmin(study_exit_date, end_2020, na.rm = T)]
  
  study_population_0 <- study_population[Dose == "0", ]
  
  divide_period_per_event <- function(.data, .x, start_period, end_period) {

    temp_0 <- .data[get(end_period) < get(.x) | is.na(get(.x)), ][, COVID19 := 0]
    temp_3 <- .data[get(start_period) > get(.x), ][, COVID19 := 1]
    temp_1 <- .data[data.table::between(get(.x), get(start_period), get(end_period)), ]
    
    temp_2 <- copy(temp_1)[, COVID19 := 0][, (end_period) := get(.x) - 1][get(end_period) >= get(start_period), ]
    temp_1 <- temp_1[, COVID19 := 1][, (start_period) := get(.x)]
    
    temp_0 <- rbind(temp_0, temp_1, temp_2, temp_3)
    return(temp_0)
  }
  
  study_population_0 <- divide_period_per_event(study_population_0, "date_event", "study_entry_date", "study_exit_date")
  study_population_0 <- study_population_0[, date_event := NULL][, Period := "0"]
  setnames(study_population_0, c("study_entry_date", "study_exit_date"), c("start_date_of_period", "end_date_of_period"))
  
  study_population_1 <- study_population[Dose != "0", ]
  
  vector_periods <- c(2, 28, 60, 180, 365)
  effective_starts <- c(0, head(vector_periods, -1))
  
  study_entry_names <- paste0("study_entry_date_", effective_starts)
  study_exit_names <- paste0("study_exit_date_", vector_periods)
  study_population_1 <- study_population_1[, (study_entry_names) := lapply(effective_starts, `+`, study_entry_date)]
  study_population_1 <- study_population_1[, (study_exit_names) := lapply(vector_periods - 1, `+`, study_entry_date)]
  setnames(study_population_1, c("study_entry_date", "study_exit_date"), c("orig_study_entry_date", "orig_study_exit_date"))
  
  study_population_1 <- data.table::melt(study_population_1, measure = list(study_entry_names, study_exit_names),
                                         variable.name = "Period", variable.factor = F,
                                         value.name = c("study_entry_date", "study_exit_date"))
  study_population_1[.(Period = as.character(seq_along(vector_periods)), to = as.character(vector_periods)), on = "Period", Period := i.to]
  study_population_1 <- study_population_1[study_entry_date <= orig_study_exit_date, ]
  study_population_1 <- study_population_1[study_exit_date > orig_study_exit_date, study_exit_date := orig_study_exit_date]
  
  study_population_1 <- divide_period_per_event(study_population_1, "date_event", "study_entry_date", "study_exit_date")
  
  study_population_1 <- study_population_1[, c("orig_study_entry_date", "orig_study_exit_date", "date_event") := NULL]
  setnames(study_population_1, c("study_entry_date", "study_exit_date"), c("start_date_of_period", "end_date_of_period"))
  
  study_population <- rbind(study_population_0, study_population_1)
  
  rm(study_population_0, study_population_1)
  
  tempname <- paste0("D3_myocard_cohort", suffix[[subpop]])
  assign(tempname, study_population)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  
  rm(list=paste0("D3_myocard_cohort", suffix[[subpop]]))
  rm(study_population)
  
}