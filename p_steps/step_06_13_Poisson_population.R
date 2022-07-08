# ----------------------------------
# Create population for poisson

# input: D3_vaxweeks_including_not_vaccinated, D3_Vaccin_cohort_cov_ALL, D3_outcomes_covid
# output: D3_vaxweeks_poisson

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  load(paste0(dirtemp,"D3_vaxweeks_including_not_vaccinated", suffix[[subpop]],".RData"))
  D3_vaxweeks <- get(paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
  rm(list=paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))

  load(paste0(diroutput,"D3_Vaccin_cohort_cov_ALL", suffix[[subpop]],".RData"))
  D3_Vaccin_cohort <- get(paste0("D3_Vaccin_cohort_cov_ALL", suffix[[subpop]]))
  rm(list=paste0("D3_Vaccin_cohort_cov_ALL", suffix[[subpop]]))
  
  load(paste0(dirtemp,"D3_outcomes_covid",suffix[[subpop]],".RData"))
  outcomes_covid <- get(paste0("D3_outcomes_covid", suffix[[subpop]]))
  rm(list=paste0("D3_outcomes_covid", suffix[[subpop]]))
  
  D3_vaxweeks <- D3_vaxweeks[week_fup <= 4, ][, week_fup := NULL]
  
  D3_vaxweeks <- D3_vaxweeks[, all_risk_factors := rowSums(.SD),
                             .SDcols = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                                         "COVSICKLE", "IMMUNOSUPPR")]
  setnames(D3_vaxweeks, c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                          "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors", "all_risk_factors"),
           c("CV_at_study_entry", "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
             "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
             "IMMUNOSUPPR_at_study_entry", "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry"))
  
  setnames(D3_vaxweeks, c("type_vax", "Dose", "sex"), c("Vaccine1", "Dose1", "Gender"))
  D3_vaxweeks <- D3_vaxweeks[is.na(Vaccine1), Vaccine1 := "none"]
  
  ori_cols <- c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
                "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP",
                "COVSICKLE_either_DX_or_DP", "IMMUNOSUPPR_at_date_vax_1")
  Vaccin_comorbidity <- D3_Vaccin_cohort[, all_risk_factors := rowSums(.SD),
                                         .SDcols = ori_cols]
  setnames(Vaccin_comorbidity, c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP",
                                 "COVHIV_either_DX_or_DP", "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP",
                                 "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP", "IMMUNOSUPPR_at_date_vax_1",
                                 "all_covariates_non_CONTR", "all_risk_factors"),
           c("CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
             "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax", "COVSICKLE_at_date_vax",
             "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax", "all_risk_factors_at_date_vax"))
  Vaccin_comorbidity <- Vaccin_comorbidity[, c("person_id", "type_vax_2", "CV_at_date_vax", "COVCANCER_at_date_vax",
                                               "COVCOPD_at_date_vax", "COVHIV_at_date_vax", "COVCKD_at_date_vax",
                                               "COVDIAB_at_date_vax", "COVOBES_at_date_vax", "COVSICKLE_at_date_vax",
                                               "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                                               "all_risk_factors_at_date_vax")]
  setnames(Vaccin_comorbidity, "type_vax_2", "Vaccine2")
  
  D3_vaxweeks <- merge(D3_vaxweeks, Vaccin_comorbidity, all.x = T, by = "person_id")
  D3_vaxweeks <- D3_vaxweeks[Vaccine2 == "0" | is.na(Vaccine2), Vaccine2 := "none"]
  D3_vaxweeks <- D3_vaxweeks[, Dose2 := 0][Dose1 == 2, Dose2 := 1][Dose1 == 2, Dose1 := 1][Dose2 == 0, Vaccine2 := "none"]
  
  rm(Vaccin_comorbidity)
  
  outcomes_covid <- outcomes_covid[, .(person_id, date_event)]
  outcomes_covid <- outcomes_covid[outcomes_covid[,.I[which.min(date_event)], by = c("person_id")]$V1]
  
  D3_vaxweeks <- merge(D3_vaxweeks, outcomes_covid, all.x = T, by = "person_id")
  rm(outcomes_covid)
  
  divided_rows <- D3_vaxweeks[data.table::between(date_event, start_date_of_period, end_date_of_period), ]
  divided_rows0 <- copy(divided_rows)[, COVID19 := 0][, end_date_of_period := date_event - 1]
  divided_rows0 <- divided_rows0[end_date_of_period >= start_date_of_period, ]
  divided_rows1 <- copy(divided_rows)[, COVID19 := 1][, start_date_of_period := date_event]
  divided_rows <- rbind(divided_rows0, divided_rows1)
  rm(divided_rows0, divided_rows1)
  
  D3_vaxweeks <- D3_vaxweeks[is.na(date_event) | !data.table::between(date_event, start_date_of_period, end_date_of_period), ]
  D3_vaxweeks <- D3_vaxweeks[, COVID19 := fifelse(is.na(date_event) | date_event > end_date_of_period, 0, 1)]
  D3_vaxweeks <- rbind(D3_vaxweeks, divided_rows)
  rm(divided_rows)
  
  D3_vaxweeks <- D3_vaxweeks[, DAP := paste(thisdatasource, subpop, sep = "_")]
  
  D3_vaxweeks <- D3_vaxweeks[, c("person_id", "start_date_of_period", "end_date_of_period", "DAP", "Gender", "ageband_at_study_entry",
                                 "COVID19", "Vaccine1", "Vaccine2", "Dose1", "Dose2", "CV_at_study_entry",
                                 "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                                 "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                                 "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                                 "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                                 "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                                 "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                                 "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                                 "all_risk_factors_at_date_vax")]
  
  tempname<-paste0("D3_vaxweeks_poisson",suffix[[subpop]])
  assign(tempname, D3_vaxweeks)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  
  rm(list=paste0("D3_vaxweeks_poisson", suffix[[subpop]]))
  rm(D3_vaxweeks)
  rm(D3_Vaccin_cohort)
}
