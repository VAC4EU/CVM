
for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  namedataset3<-paste0("D4_persontime_windows",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_windows",suffix[[subpop]],".RData"))
  persontime_windows<-get(paste0("D4_persontime_windows", suffix[[subpop]]))
  rm(list = namedataset3)
  
  persontime_windows <- persontime_windows[, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                                               "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                                               "any_risk_factors") := NULL]
  
  persontime_windows <- persontime_windows[Dose != 0 | year == 2020,][, year := NULL]
  cols_to_sums <- names(persontime_windows)[8:length(persontime_windows)]
  
  persontime_windows <- persontime_windows[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("sex", "Ageband", "Dose", "type_vax", "Period"),
                                         .SDcols = cols_to_sums]
  
  all_sex <- copy(persontime_windows)[, lapply(.SD, sum, na.rm=TRUE),
                                     by = c("Ageband", "Dose", "type_vax", "Period"),
                                     .SDcols = cols_to_sums]
  all_sex <- all_sex[, sex := "both_sexes"]
  
  persontime_windows <- rbind(persontime_windows, all_sex)
  
  all_ages <- copy(persontime_windows)[, lapply(.SD, sum, na.rm=TRUE),
                                      by = c("sex", "Dose", "type_vax", "Period"),
                                      .SDcols = cols_to_sums]
  all_ages <- unique(all_ages[, Ageband := "all_birth_cohorts"])
  
  persontime_windows <- rbind(persontime_windows, all_ages)
  
  setorder(persontime_windows, "Period")
  
  all_dose <- copy(persontime_windows)[Dose %in% c(1, 2, 3), lapply(.SD, sum, na.rm=TRUE),
                                     by = c("Ageband", "sex", "type_vax", "Period"),
                                     .SDcols = cols_to_sums]
  all_dose <- all_dose[, Dose := "both_doses"]
  
  persontime_windows <- rbind(persontime_windows, all_dose)
  
  all_fup <- copy(persontime_windows)[Period %in% c(2, 28), lapply(.SD, sum, na.rm=TRUE),
                                    by = c("Ageband", "sex", "type_vax", "Dose"),
                                    .SDcols = cols_to_sums]
  all_fup <- all_fup[, Period := "fup_until_4"]
  
  all_man <- copy(persontime_windows)[, lapply(.SD, sum, na.rm=TRUE),
                                    by = c("Ageband", "sex", "Period", "Dose"),
                                    .SDcols = cols_to_sums]
  all_man <- all_man[, type_vax := "all_manufacturer"]
  
  persontime_windows <- rbind(persontime_windows, all_man)
  
  for (i in names(persontime_windows)){
    persontime_windows[is.na(get(i)), (i) := 0]
  }
  
  # persontime_windows <- persontime_windows[, (cols_to_sums) := lapply(.SD, cumsum),
  #                                      by = c("Dose", "type_vax", "Ageband", "sex"),
  #                                      .SDcols = cols_to_sums]
  
  persontime_windows <- rbind(persontime_windows, all_fup)
  
  nameoutput3<-paste0("D4_persontime_windows_aggregated",suffix[[subpop]])
  assign(nameoutput3,rbind(persontime_windows, all_fup))
  
  save(nameoutput3,file=paste0(diroutput,nameoutput3,".RData"),list=nameoutput3)
  rm(list=nameoutput3)
  rm(nameoutput3)
  rm(all_ages, all_sex, all_man, all_fup, all_dose, cols_to_sums)
  
}
