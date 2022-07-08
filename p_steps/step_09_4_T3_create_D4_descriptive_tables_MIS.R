#-----------------------------------------------
# Create D4 descriptive tables for MIS/Myocard

# input: D3_Vaccin_cohort, D3_study_population, D3_study_population_cov_ALL, D4_population_b, D4_population_c, D4_population_d, D3_Vaccin_cohort_cov_ALL
# output: D4_descriptive_dataset_age_studystart_MIS, D4_descriptive_dataset_ageband_studystart_MIS, D4_descriptive_dataset_sex_studystart_MIS, D4_descriptive_dataset_covariate_studystart_MIS, D4_followup_fromstudystart_MIS_c, D4_followup_fromstudystart_MIS, D4_descriptive_dataset_age_vax1_MIS, D4_descriptive_dataset_ageband_vax_MIS, D4_followup_from_vax_MIS_d, D4_descriptive_dataset_covid_studystart_c_MIS, D4_descriptive_dataset_age_studystart_c_MIS, D4_descriptive_dataset_ageband_studystart_c_MIS, D4_descriptive_dataset_covariate_covid_c_MIS, D4_followup_fromstudystart_MIS_c_total, D4_descriptive_dataset_sex_vaccination_MIS

na_to_0 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
  return(DT)
}

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  if(this_datasource_has_subpopulations == T) dirD4tables <-paste0(thisdirexp,"D4 tables/")
  suppressWarnings(if (!file.exists(dirD4tables)) dir.create(file.path(dirD4tables)))
  
  load(paste0(dirtemp,"D3_Vaccin_cohort",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_study_population",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D3_study_population_cov_ALL",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_b",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_c",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_d",suffix[[subpop]],".RData"))
  
  D3_Vaccin_cohort<-get(paste0("D3_Vaccin_cohort",suffix[[subpop]]))
  D3_study_population<-get(paste0("D3_study_population",suffix[[subpop]]))
  D3_study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
  D4_population_b<-get(paste0("D4_population_b",suffix[[subpop]]))
  D4_population_c<-get(paste0("D4_population_c",suffix[[subpop]]))
  D4_population_d<-get(paste0("D4_population_d",suffix[[subpop]]))
  
  
  D4_population_c <- D4_population_c[, .(person_id, ageband_at_1_jan_2021, fup_days)]
  D4_population_d <- D4_population_d[, .(person_id, ageband_at_1_jan_2021, fup_days, type_vax_1)]
  
  D4_population_b_temp <- D4_population_b[, .(person_id, age_at_1_jan_2021, ageband_at_1_jan_2021)]
  setnames(D4_population_b_temp, c("ageband_at_1_jan_2021", "age_at_1_jan_2021"), c("ageband_at_study_entry", "age_at_study_entry"))
  
  D3_study_population <- D3_study_population[, age_at_study_entry := NULL]
  D3_study_population <- merge(D3_study_population, D4_population_b_temp, by = "person_id")
  
  D4_descriptive_dataset_age_studystart <- D3_study_population[, .(person_id, age_at_study_entry, fup_days)]
  D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                                   as.list(round(quantile(age_at_study_entry, probs = c(0.25, 0.50, 0.75)), 0))]
  D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[, c("Age_mean", "Age_min", "Age_max") :=
                                                                                   list(round(mean(age_at_study_entry), 0), min(age_at_study_entry), max(age_at_study_entry))]
  
  D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[, Followup := round(sum(fup_days) / 365.25)][, Datasource := thisdatasource]
  D4_descriptive_dataset_age_studystart <- unique(D4_descriptive_dataset_age_studystart[, .(Datasource, Followup, Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
  
  nameoutput <- paste0("D4_descriptive_dataset_age_studystart_MIS")
  assign(nameoutput, D4_descriptive_dataset_age_studystart)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  load(paste0(dirtemp,"D3_study_population",suffix[[subpop]],".RData"))
  D3_study_population<-get(paste0("D3_study_population",suffix[[subpop]]))
  
  D4_descriptive_dataset_ageband_studystart <- D3_study_population[, .(person_id, ageband_at_study_entry)]
  vect_recode <- paste0("AgeCat_", Agebands_labels)
  names(vect_recode) <- Agebands_labels
  D4_descriptive_dataset_ageband_studystart[, ageband_at_study_entry := vect_recode[ageband_at_study_entry]]
  
  D4_descriptive_dataset_ageband_studystart <- unique(D4_descriptive_dataset_ageband_studystart[, N := .N, by = "ageband_at_study_entry"][, person_id := NULL])
  D4_descriptive_dataset_ageband_studystart <- D4_descriptive_dataset_ageband_studystart[, Datasource := thisdatasource]
  D4_descriptive_dataset_ageband_studystart <- data.table::dcast(D4_descriptive_dataset_ageband_studystart, Datasource ~ ageband_at_study_entry, value.var = "N")
  
  nameoutput <- paste0("D4_descriptive_dataset_ageband_studystart_MIS")
  assign(nameoutput, D4_descriptive_dataset_ageband_studystart)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_descriptive_dataset_sex_studystart <- D3_study_population[, .(person_id, sex)]
  D4_descriptive_dataset_sex_studystart <- unique(D4_descriptive_dataset_sex_studystart[, N := .N, by = "sex"][, person_id := NULL])
  D4_descriptive_dataset_sex_studystart <- D4_descriptive_dataset_sex_studystart[, Datasource := thisdatasource]
  D4_descriptive_dataset_sex_studystart <- D4_descriptive_dataset_sex_studystart[, sex := fifelse(sex == 1, "Sex_male", "Sex_female")]
  D4_descriptive_dataset_sex_studystart <- data.table::dcast(D4_descriptive_dataset_sex_studystart, Datasource ~ sex, value.var = "N")
  
  nameoutput <- paste0("D4_descriptive_dataset_sex_studystart_MIS")
  assign(nameoutput, D4_descriptive_dataset_sex_studystart)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  setnames(D3_study_population_cov_ALL,
           c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
             "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP",
             "IMMUNOSUPPR_at_study_entry", "all_covariates_non_CONTR"),
           c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"))
  
  cols_chosen <- c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")
  D4_descriptive_dataset_covariate_studystart <- D3_study_population_cov_ALL[, lapply(.SD, sum, na.rm=TRUE), .SDcols = cols_chosen]
  D4_descriptive_dataset_covariate_studystart <- D4_descriptive_dataset_covariate_studystart[, Datasource := thisdatasource]
  D4_descriptive_dataset_covariate_studystart <- D4_descriptive_dataset_covariate_studystart[, .(Datasource, CV, Cancer, CLD, HIV, CKD, Diabetes, Obesity, Sicklecell, immunosuppressants, any_risk_factors)]
  
  nameoutput <- paste0("D4_descriptive_dataset_covariate_studystart_MIS")
  assign(nameoutput, D4_descriptive_dataset_covariate_studystart)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_followup_fromstudystart_MIS_c <- D4_population_c[, round(as.integer(sum(fup_days)) / 365), by = "ageband_at_1_jan_2021"]
  D4_followup_fromstudystart_MIS_c <- data.table::dcast(D4_followup_fromstudystart_MIS_c, . ~ ageband_at_1_jan_2021, value.var = c("V1"))
  D4_followup_fromstudystart_MIS_c <- D4_followup_fromstudystart_MIS_c[, . := NULL][, Datasource := thisdatasource]
  setnames(D4_followup_fromstudystart_MIS_c, Agebands_labels, paste0("Followup_", Agebands_labels), skip_absent=TRUE)
  
  nameoutput <- paste0("D4_followup_fromstudystart_MIS_c")
  assign(nameoutput, D4_followup_fromstudystart_MIS_c)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  
  
  
  D4_followup_fromstudystart <- D4_population_b[, .(person_id, sex, ageband_at_1_jan_2021, fup_days)]
  # D4_followup_fromstudystart <- D3_study_population[, .(person_id, sex, age_at_study_entry, study_entry_date, study_exit_date, fup_days)]
  
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, sex := fifelse(sex == 1, "Followup_males", "Followup_females")]
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, sex_value := sum(fup_days), by = "sex"]
  
  setnames(D4_followup_fromstudystart, "ageband_at_1_jan_2021", "ageband_at_study_entry")
  
  vect_recode <- paste0("Followup_", Agebands_labels)
  names(vect_recode) <- Agebands_labels
  D4_followup_fromstudystart[, ageband_at_study_entry := vect_recode[ageband_at_study_entry]]
  
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, cohort_value := sum(fup_days), by = "ageband_at_study_entry"]
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, Followup_total := sum(fup_days)]
  D4_followup_fromstudystart <- unique(D4_followup_fromstudystart[, c("person_id", "fup_days") := NULL][, Datasource := thisdatasource])
  
  D4_followup_sex <- unique(D4_followup_fromstudystart[, .(Datasource, sex, sex_value)])
  D4_followup_sex <- data.table::dcast(D4_followup_sex, Datasource ~ sex, value.var = c("sex_value"))
  
  D4_followup_cohort <- unique(D4_followup_fromstudystart[, .(Datasource, ageband_at_study_entry, cohort_value)])
  D4_followup_cohort <- data.table::dcast(D4_followup_cohort, Datasource ~ ageband_at_study_entry, value.var = c("cohort_value"))
  
  D4_followup_fromstudystart <- Reduce(merge, list(D4_followup_sex, D4_followup_cohort))
  
  vect_col_to_year <- c("Followup_males", "Followup_females", vect_recode)
  
  tot_col <- c("Datasource", vect_col_to_year)
  
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, lapply(.SD, as.numeric), .SDcols = vect_col_to_year]
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, (vect_col_to_year) := round(.SD / 365.25),
                                                           .SDcols = vect_col_to_year]
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, Datasource := thisdatasource]
  
  D4_followup_fromstudystart <- D4_followup_fromstudystart[, ..tot_col]
  
  nameoutput <- paste0("D4_followup_fromstudystart_MIS")
  assign(nameoutput, D4_followup_fromstudystart)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  rm(D4_population_b,D4_population_b_temp)
  
  
  
  
  
  
  
  
  
  
  D4_descriptive_dataset_age_vax1 <- D3_Vaccin_cohort[, .(person_id, type_vax_1, date_vax1, fup_vax1, age_at_date_vax_1)]
  D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, Month_vax1 := month(date_vax1)]
  setnames(D4_descriptive_dataset_age_vax1, "type_vax_1", "Vax_dose1")
  D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                       as.list(round(quantile(age_at_date_vax_1, probs = c(0.25, 0.50, 0.75)), 0)), by = c("Vax_dose1", "Month_vax1")]
  D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, c("Age_mean", "Age_min", "Age_max") :=
                                                                       list(round(mean(age_at_date_vax_1), 0), min(age_at_date_vax_1), max(age_at_date_vax_1)), by = c("Vax_dose1", "Month_vax1")]
  
  D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, Followup_vax1 := sum(fup_vax1), by = c("Vax_dose1", "Month_vax1")][, Datasource := thisdatasource]
  D4_descriptive_dataset_age_vax1 <- unique(D4_descriptive_dataset_age_vax1[, .(Datasource, Vax_dose1, Month_vax1, Followup_vax1,
                                                                                Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
  
  nameoutput <- paste0("D4_descriptive_dataset_age_vax1_MIS")
  assign(nameoutput, D4_descriptive_dataset_age_vax1)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_descriptive_dataset_ageband_vax <- D3_Vaccin_cohort[, .(person_id, type_vax_1, ageband_at_date_vax_1)]
  
  vect_recode <- paste0("AgeCat_", Agebands_labels)
  names(vect_recode) <- Agebands_labels
  D4_descriptive_dataset_ageband_vax[, ageband_at_date_vax_1 := vect_recode[ageband_at_date_vax_1]]
  
  D4_descriptive_dataset_ageband_vax <- unique(D4_descriptive_dataset_ageband_vax[, N := .N, by = c("ageband_at_date_vax_1", "type_vax_1")][, person_id := NULL])
  
  D4_descriptive_dataset_ageband_vax <- D4_descriptive_dataset_ageband_vax[, Datasource := thisdatasource]
  D4_descriptive_dataset_ageband_vax <- data.table::dcast(D4_descriptive_dataset_ageband_vax, Datasource + type_vax_1 ~ ageband_at_date_vax_1, value.var = "N")
  # D4_descriptive_dataset_ageband_vax <- na_to_0(D4_descriptive_dataset_ageband_vax)
  
  nameoutput <- paste0("D4_descriptive_dataset_ageband_vax_MIS")
  assign(nameoutput, D4_descriptive_dataset_ageband_vax)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_followup_from_vax_MIS_d <- D3_Vaccin_cohort[, .(person_id, type_vax_1, age_at_date_vax_1, fup_vax1)]
  rm(D3_Vaccin_cohort)
  D4_followup_from_vax_MIS_d [,age_at_date_vax_1:=cut(age_at_date_vax_1, breaks = Agebands,  labels = Agebands_labels)]
  D4_followup_from_vax_MIS_d <- D4_followup_from_vax_MIS_d[, sum(fup_vax1), by = "type_vax_1"]
  
  nameoutput <- paste0("D4_followup_from_vax_MIS_d")
  assign(nameoutput, D4_followup_from_vax_MIS_d)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  rm(D3_study_population)
  
  
  
  
  
  
  load(paste0(diroutput,"D3_study_population_cov_ALL",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_c",suffix[[subpop]],".RData"))
  
  
  D3_study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
  D4_population_c<-get(paste0("D4_population_c",suffix[[subpop]]))
  
  
  D4_descriptive_dataset_age_studystart_c <- D4_population_c[, .(person_id, age_at_1_jan_2021, ageband_at_1_jan_2021,cohort_entry_date_MIS_c,fup_days)]
  setnames(D4_descriptive_dataset_age_studystart_c, c("ageband_at_1_jan_2021", "age_at_1_jan_2021"), c("ageband_at_study_entry", "age_at_study_entry"))
  
  D4_descriptive_dataset_covid_studystart_c <- D4_descriptive_dataset_age_studystart_c[, .(person_id, cohort_entry_date_MIS_c)]
  setorder(D4_descriptive_dataset_covid_studystart_c, cohort_entry_date_MIS_c)
  D4_descriptive_dataset_covid_studystart_c[, covid_month:= as.character(substr(cohort_entry_date_MIS_c, 1, 7))][,cohort_entry_date_MIS_c:=NULL]
  D4_descriptive_dataset_covid_studystart_c <- unique(D4_descriptive_dataset_covid_studystart_c[, N := .N, by = "covid_month"][, person_id := NULL])
  D4_descriptive_dataset_covid_studystart_c <- D4_descriptive_dataset_covid_studystart_c[, Datasource := thisdatasource]
  D4_descriptive_dataset_covid_studystart_c <- data.table::dcast(D4_descriptive_dataset_covid_studystart_c, Datasource ~ covid_month, value.var = "N")
  
  
  nameoutput <- paste0("D4_descriptive_dataset_covid_studystart_c_MIS")
  assign(nameoutput, D4_descriptive_dataset_covid_studystart_c)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_descriptive_dataset_age_studystart_c <- D4_descriptive_dataset_age_studystart_c[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                                       as.list(round(quantile(age_at_study_entry, probs = c(0.25, 0.50, 0.75)), 0))]
  D4_descriptive_dataset_age_studystart_c <- D4_descriptive_dataset_age_studystart_c[, c("Age_mean", "Age_min", "Age_max") :=
                                                                                       list(round(mean(age_at_study_entry), 0), min(age_at_study_entry), max(age_at_study_entry))]
  
  D4_descriptive_dataset_age_studystart_c <- D4_descriptive_dataset_age_studystart_c[, Followup := round(sum(fup_days) / 365.25)][, Datasource := thisdatasource]
  D4_descriptive_dataset_age_studystart_c <- unique(D4_descriptive_dataset_age_studystart_c[, .(Datasource, Followup, Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
  
  nameoutput <- paste0("D4_descriptive_dataset_age_studystart_c_MIS")
  assign(nameoutput, D4_descriptive_dataset_age_studystart_c)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  D4_descriptive_dataset_ageband_studystart_c <- D4_population_c[, .(person_id, ageband_at_1_jan_2021)]
  setnames(D4_descriptive_dataset_ageband_studystart_c, "ageband_at_1_jan_2021", "ageband_at_study_entry")
  
  vect_recode <- paste0("AgeCat_", Agebands_labels)
  names(vect_recode) <- Agebands_labels
  D4_descriptive_dataset_ageband_studystart_c[, ageband_at_study_entry := vect_recode[ageband_at_study_entry]]
  
  D4_descriptive_dataset_ageband_studystart_c <- unique(D4_descriptive_dataset_ageband_studystart_c[, N := .N, by = "ageband_at_study_entry"][, person_id := NULL])
  D4_descriptive_dataset_ageband_studystart_c <- D4_descriptive_dataset_ageband_studystart_c[, Datasource := thisdatasource]
  D4_descriptive_dataset_ageband_studystart_c <- data.table::dcast(D4_descriptive_dataset_ageband_studystart_c, Datasource ~ ageband_at_study_entry, value.var = "N")
  cols_to_insert <- setdiff(paste0("AgeCat_", Agebands_labels), colnames(D4_descriptive_dataset_ageband_studystart_c))
  if (length(cols_to_insert) != 0) {
    D4_descriptive_dataset_ageband_studystart_c[, (cols_to_insert) := 0]
  }
  
  nameoutput <- paste0("D4_descriptive_dataset_ageband_studystart_c_MIS")
  assign(nameoutput, D4_descriptive_dataset_ageband_studystart_c)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  
  setnames(D4_population_c,
           c("CV_at_covid", "COVCANCER_at_covid", "COVCOPD_at_covid", "COVHIV_at_covid", "COVCKD_at_covid",
             "COVDIAB_at_covid", "COVOBES_at_covid", "COVSICKLE_at_covid", "immunosuppressants_at_covid", "at_risk_at_covid"),
           c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"))
  
  cols_chosen <- c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")
  D4_descriptive_dataset_covariate_covid_c_MIS<-merge(D4_population_c,D3_study_population_cov_ALL,all.x=T,by="person_id")
  
  D4_descriptive_dataset_covariate_covid_c_MIS <- D4_population_c[, lapply(.SD, sum, na.rm=TRUE), .SDcols = cols_chosen]
  D4_descriptive_dataset_covariate_covid_c_MIS <- D4_descriptive_dataset_covariate_covid_c_MIS[, Datasource := thisdatasource]
  D4_descriptive_dataset_covariate_covid_c_MIS <- D4_descriptive_dataset_covariate_covid_c_MIS[, .(Datasource, CV, Cancer, CLD, HIV, CKD, Diabetes, Obesity, Sicklecell, immunosuppressants, any_risk_factors)]
  
  nameoutput <- paste0("D4_descriptive_dataset_covariate_covid_c_MIS")
  assign(nameoutput, D4_descriptive_dataset_covariate_covid_c_MIS)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  D4_followup_fromstudystart_MIS_c <- D4_population_c[, sum(fup_days)]
  rm(D4_population_c)
  D4_followup_fromstudystart_MIS_c <- data.table(total = D4_followup_fromstudystart_MIS_c)
  
  nameoutput <- paste0("D4_followup_fromstudystart_MIS_c_total")
  assign(nameoutput, D4_followup_fromstudystart_MIS_c)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  
  
  
  load(paste0(diroutput,"D4_population_d",suffix[[subpop]],".RData"))
  D4_population_d<-get(paste0("D4_population_d",suffix[[subpop]]))
  
  
  D4_descriptive_dataset_sex_vaccination <- D4_population_d[, .(person_id, sex, type_vax_1)]
  D4_descriptive_dataset_sex_vaccination <- unique(D4_descriptive_dataset_sex_vaccination[, N := .N, by = c("sex", "type_vax_1")][, person_id := NULL])
  D4_descriptive_dataset_sex_vaccination <- D4_descriptive_dataset_sex_vaccination[, Datasource := thisdatasource]
  D4_descriptive_dataset_sex_vaccination <- D4_descriptive_dataset_sex_vaccination[, sex := fifelse(sex == 1, "Sex_male", "Sex_female")]
  D4_descriptive_dataset_sex_vaccination <- data.table::dcast(D4_descriptive_dataset_sex_vaccination, Datasource + type_vax_1 ~ sex, value.var = "N")
  
  nameoutput <- paste0("D4_descriptive_dataset_sex_vaccination_MIS")
  assign(nameoutput, D4_descriptive_dataset_sex_vaccination)
  fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
  rm(list=nameoutput)
  rm(D4_population_d)
  
}



