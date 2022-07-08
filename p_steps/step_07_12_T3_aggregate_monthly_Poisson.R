#-----------------------------------------------
# Aggregate persontime datasets for poisson

# input: D4_persontime_risk_month_poisson
# output: D4_persontime_monthly_poisson_RF

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  load(paste0(diroutput,"D4_persontime_risk_month_poisson", suffix[[subpop]],".RData"))
  D4_persontime_poisson <- get(paste0("D4_persontime_risk_month_poisson", suffix[[subpop]]))
  rm(list=paste0("D4_persontime_risk_month_poisson", suffix[[subpop]]))
  
  D4_persontime_poisson$year = as.integer(lapply(strsplit(D4_persontime_poisson$month, split = "-"), "[", 1))
  D4_persontime_poisson$month = as.integer(lapply(strsplit(D4_persontime_poisson$month, split = "-"), "[", 2))
  
  D4_persontime_poisson <- D4_persontime_poisson[, Dose1 := as.integer(Dose1)]
  
  col_names <- copy(colnames(D4_persontime_poisson))[colnames(D4_persontime_poisson) %not in% c("DAP", "Gender", "ageband_at_study_entry",
                                                                                                "month", "year", "COVID19", "Vaccine1", "Vaccine2", "Dose1", "Dose2", "CV_at_study_entry",
                                                                                                "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                                                                                                "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                                                                                                "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                                                                                                "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                                                                                                "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                                                                                                "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                                                                                                "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                                                                                                "all_risk_factors_at_date_vax")]
  D4_persontime_poisson <- D4_persontime_poisson[, c("DAP", "Gender", "ageband_at_study_entry", "month", "year", "COVID19", "Vaccine1", "Vaccine2",
                                                     "Dose1", "Dose2", "CV_at_study_entry",
                                                     "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                                                     "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                                                     "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                                                     "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                                                     "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                                                     "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                                                     "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                                                     "all_risk_factors_at_date_vax", ..col_names)]
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  nameoutput1<-paste0("D4_persontime_monthly_poisson_RF",suffix[[subpop]])
  assign(nameoutput1, D4_persontime_poisson)
  
  save(nameoutput1,file=paste0(thisdirexp,nameoutput1,".RData"), list=nameoutput1)
  rm(list=nameoutput1)
  rm(col_names, D4_persontime_poisson)
  
}