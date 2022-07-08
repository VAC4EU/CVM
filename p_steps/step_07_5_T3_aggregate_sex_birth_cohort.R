#-----------------------------------------------
# Aggregate persontime datasets by birth cohorts

# input: D4_persontime_risk_week, D4_persontime_benefit_week, D4_persontime_risk_year, D4_persontime_benefit_year
# output: D4_persontime_risk_week_BC, D4_persontime_benefit_week_BC, D4_persontime_risk_year_BC, D4_persontime_benefit_year_BC

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  namedataset1<-paste0("D4_persontime_risk_week",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_risk_week",suffix[[subpop]],".RData"))
  
  get(namedataset1)[, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                        "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                        "any_risk_factors") := NULL]
  
  cols_to_sums <- names(get(namedataset1))[6:length(get(namedataset1))]
  
  assign(namedataset1, get(namedataset1)[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("sex", "ageband_at_study_entry", "Dose", "type_vax", "week"),
                                         .SDcols = cols_to_sums])
  
  all_sex <- copy(get(namedataset1))[, lapply(.SD, sum, na.rm=TRUE),
                                     by = c("ageband_at_study_entry", "Dose", "type_vax", "week"),
                                     .SDcols = cols_to_sums]
  all_sex <- all_sex[, sex := "both_sexes"]
  
  assign(namedataset1,rbind(get(namedataset1), all_sex))
  
  all_ages <- copy(get(namedataset1))[, lapply(.SD, sum, na.rm=TRUE),
                                      by = c("sex", "Dose", "type_vax", "week"),
                                      .SDcols = cols_to_sums]
  all_ages <- unique(all_ages[, ageband_at_study_entry := "all_birth_cohorts"])
  
  assign(namedataset1,rbind(get(namedataset1), all_ages))
  
  nameoutput1<-paste0("D4_persontime_risk_week_BC",suffix[[subpop]])
  assign(nameoutput1 , bc_divide_60(get(namedataset1), c("sex", "Dose", "type_vax", "week"), cols_to_sums))
  
  save(nameoutput1,file=paste0(diroutput,nameoutput1,".RData"),list=nameoutput1)
  rm(list=nameoutput1)
  rm(list=namedataset1)
  rm(namedataset1,nameoutput1)
  
  
  
  namedataset2<-paste0("D4_persontime_benefit_week",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_benefit_week",suffix[[subpop]],".RData"))
  
 assign(namedataset2,get(namedataset2)[, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                                                               "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                                                               "any_risk_factors") := NULL]
 )
  
  cols_to_sums <- names(get(namedataset2))[6:length(get(namedataset2))]
  
  assign(namedataset2, get(namedataset2)[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("sex", "ageband_at_study_entry", "Dose", "type_vax", "week"),
                                         .SDcols = cols_to_sums])

  all_sex <- copy(get(namedataset2))[, lapply(.SD, sum, na.rm=TRUE),
                                              by = c("ageband_at_study_entry", "Dose", "type_vax", "week"),
                                              .SDcols = cols_to_sums]
  all_sex <- all_sex[, sex := "both_sexes"]
  
  assign(namedataset2, rbind(get(namedataset2), all_sex))
  
  all_ages <- copy(get(namedataset2))[, lapply(.SD, sum, na.rm=TRUE),
                                               by = c("sex", "Dose", "type_vax", "week"),
                                               .SDcols = cols_to_sums]
  all_ages <- unique(all_ages[, ageband_at_study_entry := "all_birth_cohorts"])
  
  assign(namedataset2,rbind(get(namedataset2), all_ages))
  
  nameoutput2<-paste0("D4_persontime_benefit_week_BC",suffix[[subpop]])
  assign(nameoutput2,bc_divide_60(get(namedataset2), c("sex", "Dose", "type_vax", "week"), cols_to_sums))
  
  save(nameoutput2,file=paste0(diroutput,nameoutput2,".RData"),list=nameoutput2)
  rm(list=nameoutput2)
  rm(list=namedataset2)
  rm(namedataset2,nameoutput2)
  
  
  namedataset3<-paste0("D4_persontime_risk_year",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_risk_year",suffix[[subpop]],".RData"))
  
  assign(namedataset3,get(namedataset3)[, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                                                         "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                                                         "any_risk_factors") := NULL])
  
  assign(namedataset3, get(namedataset3)[Dose != 0 | year == 2020,][, year := NULL])
  cols_to_sums <- names(get(namedataset3))[6:length(get(namedataset3))]
  
  assign(namedataset3, get(namedataset3)[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("sex", "ageband_at_study_entry", "Dose", "type_vax", "week_fup"),
                                         .SDcols = cols_to_sums])
  
  all_sex <- copy(get(namedataset3))[, lapply(.SD, sum, na.rm=TRUE),
                                           by = c("ageband_at_study_entry", "Dose", "type_vax", "week_fup"),
                                           .SDcols = cols_to_sums]
  all_sex <- all_sex[, sex := "both_sexes"]
  
  assign(namedataset3,rbind(get(namedataset3), all_sex))
  
  all_ages <- copy(get(namedataset3))[, lapply(.SD, sum, na.rm=TRUE),
                                            by = c("sex", "Dose", "type_vax", "week_fup"),
                                            .SDcols = cols_to_sums]
  all_ages <- unique(all_ages[, ageband_at_study_entry := "all_birth_cohorts"])
  
  assign(namedataset3, rbind(get(namedataset3), all_ages))

  assign(namedataset3, bc_divide_60(get(namedataset3), c("sex", "Dose", "type_vax", "week_fup"), cols_to_sums))
  
  setorder(get(namedataset3), "ageband_at_study_entry")
  
  vax_dose <- unique(copy(get(namedataset3))[, c("Dose", "type_vax", "week_fup")])
  vax_dose <- vax_dose[, .SD[which.max(as.integer(week_fup))], by = c("Dose", "type_vax")]
  vax_dose_0 <- vax_dose[week_fup == 0,]
  week_vax_dose <- vax_dose[, lapply(.SD, rep, vax_dose[, week_fup])][, week_fup := unlist(
    lapply(vax_dose[, week_fup], seq_len))]
  week_vax_dose <- rbind(week_vax_dose, vax_dose_0)
  
  sex_vect <- rep(c("0", "1", "both_sexes"), each = nrow(week_vax_dose))
  week_vax_dose <- week_vax_dose[, lapply(.SD, rep, 3)][, sex := sex_vect]
  observed_agebands <- unique(copy(get(namedataset3))[, ageband_at_study_entry])
  Agebands_vect <- rep(observed_agebands, each = nrow(week_vax_dose))
  empty_risk_year <- week_vax_dose[, lapply(.SD, rep, length(observed_agebands))][, ageband_at_study_entry := Agebands_vect]
  
  nameoutput3<-paste0("D4_persontime_risk_year_BC",suffix[[subpop]])
  assign(nameoutput3, merge(empty_risk_year, get(namedataset3), all.x = T,
                                      by = c("Dose", "type_vax", "week_fup", "sex", "ageband_at_study_entry")))
  
  all_dose <- copy(get(nameoutput3))[Dose %in% c(1, 2), lapply(.SD, sum, na.rm=TRUE),
                                      by = c("ageband_at_study_entry", "sex", "type_vax", "week_fup"),
                                      .SDcols = cols_to_sums]
  all_dose <- all_dose[, Dose := "both_doses"]
  
  assign(nameoutput3,rbind(get(nameoutput3), all_dose))
  
  all_fup <- copy(get(nameoutput3))[week_fup %in% c(1, 2, 3, 4), lapply(.SD, sum, na.rm=TRUE),
                                     by = c("ageband_at_study_entry", "sex", "type_vax", "Dose"),
                                     .SDcols = cols_to_sums]
  all_fup <- all_fup[, week_fup := "fup_until_4"]
  
  all_man <- copy(get(nameoutput3))[, lapply(.SD, sum, na.rm=TRUE),
                                     by = c("ageband_at_study_entry", "sex", "week_fup", "Dose"),
                                     .SDcols = cols_to_sums]
  all_man <- all_man[, type_vax := "all_manufacturer"]
  
  assign(nameoutput3,rbind(get(nameoutput3), all_man))
  
  for (i in names(get(nameoutput3))){
    get(nameoutput3)[is.na(get(i)), (i) := 0]
  }
  
  assign(nameoutput3, get(nameoutput3)[, (cols_to_sums) := lapply(.SD, cumsum),
                                       by = c("Dose", "type_vax", "ageband_at_study_entry", "sex"),
                                       .SDcols = cols_to_sums])
  
  assign(nameoutput3,rbind(get(nameoutput3), all_fup))
  
  save(nameoutput3,file=paste0(diroutput,nameoutput3,".RData"),list=nameoutput3)
  rm(list=nameoutput3)
  rm(list=namedataset3)
  rm(namedataset3,nameoutput3)
  
  
  namedataset4<-paste0("D4_persontime_benefit_year",suffix[[subpop]])
  load(paste0(diroutput,namedataset4,".RData"))
  
  assign(namedataset4, get(namedataset4)[, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                                                               "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                                                               "any_risk_factors") := NULL])
  
  assign(namedataset4,get(namedataset4)[Dose != 0 | year == 2020,][, year := NULL])
  cols_to_sums <- names(get(namedataset4))[6:length(get(namedataset4))]
  
  assign(namedataset4, get(namedataset4)[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("sex", "ageband_at_study_entry", "Dose", "type_vax", "week_fup"),
                                         .SDcols = cols_to_sums])
  all_sex <- copy(get(namedataset4))[, lapply(.SD, sum, na.rm=TRUE),
                                              by = c("ageband_at_study_entry", "Dose", "type_vax", "week_fup"),
                                              .SDcols = cols_to_sums]
  all_sex <- all_sex[, sex := "both_sexes"]
  
  assign(namedataset4,rbind(get(namedataset4), all_sex))
  
  all_ages <- copy(get(namedataset4))[, lapply(.SD, sum, na.rm=TRUE),
                                               by = c("sex", "Dose", "type_vax", "week_fup"),
                                               .SDcols = cols_to_sums]
  all_ages <- unique(all_ages[, ageband_at_study_entry := "all_birth_cohorts"])
  assign(namedataset4,rbind(get(namedataset4), all_ages))
  
  nameoutput4<-paste0("D4_persontime_benefit_year_BC",suffix[[subpop]])
  assign(nameoutput4,bc_divide_60(get(namedataset4), c("sex", "Dose", "type_vax", "week_fup")
                                                , cols_to_sums))
  
  save(nameoutput4,file=paste0(diroutput,nameoutput4,".RData"),list=nameoutput4)
  rm(list=nameoutput4)
  rm(list=namedataset4)
  rm(namedataset4,nameoutput4)
  
  rm(all_ages, all_sex, cols_to_sums,empty_risk_year,vax_dose,vax_dose_0,week_vax_dose)
  
}
