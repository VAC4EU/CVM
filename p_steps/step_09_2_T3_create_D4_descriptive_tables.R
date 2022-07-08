#-----------------------------------------------
# Create D4 descriptive tables

# input: D3_Vaccin_cohort, D3_study_population, D3_study_population_cov_ALL, D3_Vaccin_cohort_cov_ALL
# output: D4_descriptive_dataset_age_studystart, D4_descriptive_dataset_ageband_studystart, D4_descriptive_dataset_sex_studystart, D4_descriptive_dataset_covariate_studystart, D4_descriptive_dataset_covariate_vax, D4_followup_fromstudystart, D4_descriptive_dataset_age_vax1, D4_descriptive_dataset_age_vax2, D4_descriptive_dataset_ageband_vax, D4_descriptive_dataset_sex_vaccination, D4_followup_from_vax, Density_plot_distance_doses, Histogram_distance_doses, D4_distance_doses


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
   load(paste0(diroutput,"D3_Vaccin_cohort_cov_ALL",suffix[[subpop]],".RData"))
   
   
   Vaccin_cohort<-get(paste0("D3_Vaccin_cohort",suffix[[subpop]]))
   rm(list=paste0("D3_Vaccin_cohort",suffix[[subpop]]))
   study_population<-get(paste0("D3_study_population",suffix[[subpop]]))
   rm(list=paste0("D3_study_population",suffix[[subpop]]))
   study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
   rm(list=paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
   Vaccin_cohort_cov_ALL<-get(paste0("D3_Vaccin_cohort_cov_ALL",suffix[[subpop]]))
   rm(list=paste0("D3_Vaccin_cohort_cov_ALL",suffix[[subpop]]))
   
   
   D4_descriptive_dataset_age_studystart <- study_population[, .(person_id, age_at_study_entry, fup_days)]
   D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                                     as.list(round(quantile(age_at_study_entry, probs = c(0.25, 0.50, 0.75)), 0))]
   D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[, c("Age_mean", "Age_min", "Age_max") :=
                                                                                     list(round(mean(age_at_study_entry), 0), min(age_at_study_entry), max(age_at_study_entry))]
   
   D4_descriptive_dataset_age_studystart <- D4_descriptive_dataset_age_studystart[, Followup := round(sum(fup_days) / 365.25)][, Datasource := thisdatasource]
   D4_descriptive_dataset_age_studystart <- unique(D4_descriptive_dataset_age_studystart[, .(Datasource, Followup, Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
   
   nameoutput <- paste0("D4_descriptive_dataset_age_studystart")
   assign(nameoutput, D4_descriptive_dataset_age_studystart)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_descriptive_dataset_ageband_studystart <- study_population[, .(person_id, ageband_at_study_entry)]
   vect_recode <- paste0("AgeCat_", Agebands_labels)
   names(vect_recode) <- Agebands_labels
   D4_descriptive_dataset_ageband_studystart[, ageband_at_study_entry := vect_recode[ageband_at_study_entry]]
   
   D4_descriptive_dataset_ageband_studystart <- unique(D4_descriptive_dataset_ageband_studystart[, N := .N, by = "ageband_at_study_entry"][, person_id := NULL])
   older60 <- copy(D4_descriptive_dataset_ageband_studystart)[ageband_at_study_entry %in% paste0("AgeCat_", Agebands60), sum(N)]
   older60 <- data.table::data.table(ageband_at_study_entry = "AgeCat_60+", N = older60)
   D4_descriptive_dataset_ageband_studystart <- rbind(D4_descriptive_dataset_ageband_studystart, older60)
   
   D4_descriptive_dataset_ageband_studystart <- D4_descriptive_dataset_ageband_studystart[, Datasource := thisdatasource]
   D4_descriptive_dataset_ageband_studystart <- data.table::dcast(D4_descriptive_dataset_ageband_studystart, Datasource ~ ageband_at_study_entry, value.var = "N")
   
   nameoutput <- paste0("D4_descriptive_dataset_ageband_studystart")
   assign(nameoutput, D4_descriptive_dataset_ageband_studystart)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   D4_descriptive_dataset_sex_studystart <- study_population[, .(person_id, sex)]
   D4_descriptive_dataset_sex_studystart <- unique(D4_descriptive_dataset_sex_studystart[, N := .N, by = "sex"][, person_id := NULL])
   D4_descriptive_dataset_sex_studystart <- D4_descriptive_dataset_sex_studystart[, Datasource := thisdatasource]
   D4_descriptive_dataset_sex_studystart <- D4_descriptive_dataset_sex_studystart[, sex := fifelse(sex == 1, "Sex_male", "Sex_female")]
   D4_descriptive_dataset_sex_studystart <- data.table::dcast(D4_descriptive_dataset_sex_studystart, Datasource ~ sex, value.var = "N")
   
   nameoutput <- paste0("D4_descriptive_dataset_sex_studystart")
   assign(nameoutput, D4_descriptive_dataset_sex_studystart)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   setnames(study_population_cov_ALL,
            c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
              "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP",
              "IMMUNOSUPPR_at_study_entry", "all_covariates_non_CONTR"),
            c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"))
   
   cols_chosen <- c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")
   D4_descriptive_dataset_covariate_studystart <- study_population_cov_ALL[, lapply(.SD, sum, na.rm=TRUE), .SDcols = cols_chosen]
   D4_descriptive_dataset_covariate_studystart <- D4_descriptive_dataset_covariate_studystart[, Datasource := thisdatasource]
   D4_descriptive_dataset_covariate_studystart <- D4_descriptive_dataset_covariate_studystart[, .(Datasource, CV, Cancer, CLD, HIV, CKD, Diabetes, Obesity, Sicklecell, immunosuppressants, any_risk_factors)]
   
   nameoutput <- paste0("D4_descriptive_dataset_covariate_studystart")
   assign(nameoutput, D4_descriptive_dataset_covariate_studystart)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   setnames(Vaccin_cohort_cov_ALL,
            c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
              "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP",
              "IMMUNOSUPPR_at_date_vax_1", "all_covariates_non_CONTR"),
            c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"))
   
   cols_chosen <- c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")
   D4_descriptive_dataset_covariate_vax <- Vaccin_cohort_cov_ALL[, lapply(.SD, sum, na.rm=TRUE), .SDcols = cols_chosen]
   rm(Vaccin_cohort_cov_ALL)
   D4_descriptive_dataset_covariate_vax <- D4_descriptive_dataset_covariate_vax[, Datasource := thisdatasource]
   D4_descriptive_dataset_covariate_vax <- D4_descriptive_dataset_covariate_vax[, .(Datasource, CV, Cancer, CLD, HIV, CKD, Diabetes, Obesity, Sicklecell, immunosuppressants, any_risk_factors)]
   
   nameoutput <- paste0("D4_descriptive_dataset_covariate_vax")
   assign(nameoutput, D4_descriptive_dataset_covariate_vax)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_followup_fromstudystart <- study_population[, .(person_id, sex, ageband_at_study_entry, study_entry_date, study_exit_date, fup_days)]
   dec31 = ymd(20201231)
   jan1 = ymd(20210101)
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, Followup_2020 := fifelse(study_exit_date > dec31, correct_difftime(dec31, study_entry_date), fup_days)]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[study_exit_date > dec31, Followup_2021 := correct_difftime(study_exit_date, jan1)]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, sex := fifelse(sex == 1, "Followup_males", "Followup_females")]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, sex_value := sum(fup_days), by = "sex"]
   
   vect_recode <- paste0("Followup_", Agebands_labels)
   names(vect_recode) <- Agebands_labels
   D4_followup_fromstudystart[, ageband_at_study_entry := vect_recode[ageband_at_study_entry]]
   
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, cohort_value := sum(fup_days), by = "ageband_at_study_entry"]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, Followup_total := sum(fup_days)]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, c("Followup_2020", "Followup_2021") := list(sum(Followup_2020), sum(Followup_2021, na.rm = T))]
   D4_followup_fromstudystart <- unique(D4_followup_fromstudystart[, c("person_id", "fup_days") := NULL][, Datasource := thisdatasource])
   
   D4_followup_sex <- unique(D4_followup_fromstudystart[, .(Datasource, sex, sex_value)])
   D4_followup_sex <- data.table::dcast(D4_followup_sex, Datasource ~ sex, value.var = c("sex_value"))
   
   D4_followup_cohort <- unique(D4_followup_fromstudystart[, .(Datasource, ageband_at_study_entry, cohort_value)])
   older60 <- copy(D4_followup_cohort)[ageband_at_study_entry %in% paste0("Followup_", Agebands60), sum(cohort_value)]
   older60 <- data.table::data.table(Datasource = thisdatasource, ageband_at_study_entry = "Followup_60+", cohort_value = older60)
   D4_followup_cohort <- rbind(D4_followup_cohort, older60)
   
   D4_followup_cohort <- data.table::dcast(D4_followup_cohort, Datasource ~ ageband_at_study_entry, value.var = c("cohort_value"))
   
   D4_followup_complete <- unique(D4_followup_fromstudystart[, .(Datasource, Followup_total, Followup_2020, Followup_2021)])
   
   D4_followup_fromstudystart <- Reduce(merge, list(D4_followup_sex, D4_followup_cohort, D4_followup_complete))
   
   vect_col_to_year <- c("Followup_males", "Followup_females", "Followup_total",
                         vect_recode, "Followup_60+", "Followup_2020", "Followup_2021")
   
   tot_col <- c("Datasource", vect_col_to_year)
   
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, lapply(.SD, as.numeric), .SDcols = vect_col_to_year]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, (vect_col_to_year) := round(.SD / 365.25),
                                                            .SDcols = vect_col_to_year]
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, Datasource := thisdatasource]
   
   D4_followup_fromstudystart <- D4_followup_fromstudystart[, ..tot_col]
   
   nameoutput <- paste0("D4_followup_fromstudystart")
   assign(nameoutput, D4_followup_fromstudystart)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   
   D4_descriptive_dataset_age_vax1 <- Vaccin_cohort[, .(person_id, type_vax_1, date_vax1, fup_vax1, age_at_date_vax_1)]
   D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, Month_vax1 := month(date_vax1)]
   setnames(D4_descriptive_dataset_age_vax1, "type_vax_1", "Vax_dose1")
   D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                         as.list(round(quantile(age_at_date_vax_1, probs = c(0.25, 0.50, 0.75)), 0)), by = c("Vax_dose1", "Month_vax1")]
   D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, c("Age_mean", "Age_min", "Age_max") :=
                                                                         list(round(mean(age_at_date_vax_1), 0), min(age_at_date_vax_1), max(age_at_date_vax_1)), by = c("Vax_dose1", "Month_vax1")]
   
   D4_descriptive_dataset_age_vax1 <- D4_descriptive_dataset_age_vax1[, Followup_vax1 := sum(fup_vax1), by = c("Vax_dose1", "Month_vax1")][, Datasource := thisdatasource]
   D4_descriptive_dataset_age_vax1 <- unique(D4_descriptive_dataset_age_vax1[, .(Datasource, Vax_dose1, Month_vax1, Followup_vax1,
                                                                                 Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
   
   nameoutput <- paste0("D4_descriptive_dataset_age_vax1")
   assign(nameoutput, D4_descriptive_dataset_age_vax1)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_descriptive_dataset_age_vax2 <- Vaccin_cohort[, .(person_id, type_vax_2, date_vax2, fup_vax2, age_at_date_vax_2)]
   D4_descriptive_dataset_age_vax2 <- D4_descriptive_dataset_age_vax2[!is.na(type_vax_2)]
   D4_descriptive_dataset_age_vax2 <- D4_descriptive_dataset_age_vax2[, Month_vax2 := month(date_vax2)]
   setnames(D4_descriptive_dataset_age_vax2, "type_vax_2", "Vax_dose2")
   D4_descriptive_dataset_age_vax2 <- D4_descriptive_dataset_age_vax2[,c("Age_P25", "Age_P50", "Age_p75") :=
                                                                         as.list(round(quantile(age_at_date_vax_2, probs = c(0.25, 0.50, 0.75)), 0)), by = c("Vax_dose2", "Month_vax2")]
   D4_descriptive_dataset_age_vax2 <- D4_descriptive_dataset_age_vax2[, c("Age_mean", "Age_min", "Age_max") :=
                                                                         list(round(mean(age_at_date_vax_2), 0), min(age_at_date_vax_2), max(age_at_date_vax_2)), by = c("Vax_dose2", "Month_vax2")]
   
   D4_descriptive_dataset_age_vax2 <- D4_descriptive_dataset_age_vax2[, Followup_vax2 := sum(fup_vax2), by = c("Vax_dose2", "Month_vax2")][, Datasource := thisdatasource]
   D4_descriptive_dataset_age_vax2 <- unique(D4_descriptive_dataset_age_vax2[, .(Datasource, Vax_dose2, Month_vax2, Followup_vax2,
                                                                                 Age_P25, Age_P50, Age_p75, Age_mean, Age_min, Age_max)])
   
   nameoutput <- paste0("D4_descriptive_dataset_age_vax2")
   assign(nameoutput, D4_descriptive_dataset_age_vax2)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_descriptive_dataset_ageband_vax <- Vaccin_cohort[, .(person_id, type_vax_1, ageband_at_date_vax_1)]
   
   vect_recode <- paste0("Agecat_", Agebands_labels)
   names(vect_recode) <- Agebands_labels
   D4_descriptive_dataset_ageband_vax[, ageband_at_date_vax_1 := vect_recode[ageband_at_date_vax_1]]
   
   D4_descriptive_dataset_ageband_vax <- unique(D4_descriptive_dataset_ageband_vax[, N := .N, by = c("ageband_at_date_vax_1",
                                                                                                     "type_vax_1")][, person_id := NULL])
   
   older60 <- copy(D4_descriptive_dataset_ageband_vax)[ageband_at_date_vax_1 %in% paste0("Agecat_", Agebands60),
                                                       lapply(.SD, sum, na.rm=TRUE), by = c("type_vax_1"),
                                                       .SDcols = "N"]
   older60 <- unique(older60[, ageband_at_date_vax_1 := "Agecat_60+"])
   D4_descriptive_dataset_ageband_vax <- rbind(D4_descriptive_dataset_ageband_vax, older60)
   
   D4_descriptive_dataset_ageband_vax <- D4_descriptive_dataset_ageband_vax[, Datasource := thisdatasource]
   D4_descriptive_dataset_ageband_vax <- data.table::dcast(D4_descriptive_dataset_ageband_vax, Datasource + type_vax_1 ~ ageband_at_date_vax_1, value.var = "N")
   # D4_descriptive_dataset_ageband_vax <- na_to_0(D4_descriptive_dataset_ageband_vax)
   D4_descriptive_dataset_ageband_vax <- na_to_0(D4_descriptive_dataset_ageband_vax)
   
   nameoutput <- paste0("D4_descriptive_dataset_ageband_vax")
   assign(nameoutput, D4_descriptive_dataset_ageband_vax)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_descriptive_dataset_sex_vaccination <- Vaccin_cohort[, .(person_id, sex, type_vax_1)]
   D4_descriptive_dataset_sex_vaccination <- unique(D4_descriptive_dataset_sex_vaccination[, N := .N, by = c("sex", "type_vax_1")][, person_id := NULL])
   D4_descriptive_dataset_sex_vaccination <- D4_descriptive_dataset_sex_vaccination[, Datasource := thisdatasource]
   D4_descriptive_dataset_sex_vaccination <- D4_descriptive_dataset_sex_vaccination[, sex := fifelse(sex == 1, "Sex_male", "Sex_female")]
   D4_descriptive_dataset_sex_vaccination <- data.table::dcast(D4_descriptive_dataset_sex_vaccination, Datasource + type_vax_1 ~ sex, value.var = "N")
   
   nameoutput <- paste0("D4_descriptive_dataset_sex_vaccination")
   assign(nameoutput, D4_descriptive_dataset_sex_vaccination)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_followup_from_vax <- study_population[, .(person_id, sex, ageband_at_date_vax_1, fup_vax1, fup_vax2)]
   
   followup_vax1 <- D4_followup_from_vax[!is.na(fup_vax1), ][, sex := fifelse(sex == 1, "Followup_males_vax1", "Followup_females_vax1")]
   followup_vax2 <- D4_followup_from_vax[!is.na(fup_vax2), ][, sex := fifelse(sex == 1, "Followup_males_vax2", "Followup_females_vax2")]
   followup_vax1 <- followup_vax1[, sex_value := sum(fup_vax1), by = "sex"]
   followup_vax2 <- followup_vax2[, sex_value := sum(fup_vax2), by = "sex"]
   
   vect_recode <- paste0("Followup_", Agebands_labels)
   names(vect_recode) <- Agebands_labels
   followup_vax1[, ageband_at_date_vax_1 := paste0(vect_recode, "_vax1")[ageband_at_date_vax_1]]
   followup_vax2[, ageband_at_date_vax_1 := paste0(vect_recode, "_vax2")[ageband_at_date_vax_1]]
   
   followup_vax1 <- followup_vax1[, cohort_value := sum(fup_vax1), by = "ageband_at_date_vax_1"]
   followup_vax2 <- followup_vax2[, cohort_value := sum(fup_vax2), by = "ageband_at_date_vax_1"]
   
   followup_vax1 <- followup_vax1[, Followup_total_vax1 := sum(fup_vax1)]
   followup_vax2 <- followup_vax2[, Followup_total_vax2 := sum(fup_vax2)]
   
   followup_vax1 <- unique(followup_vax1[, c("person_id", "fup_vax1") := NULL][, Datasource := thisdatasource])
   followup_vax2 <- unique(followup_vax2[, c("person_id", "fup_vax2") := NULL][, Datasource := thisdatasource])
   
   followup_vax1_sex <- unique(followup_vax1[, .(Datasource, sex, sex_value)])
   followup_vax1_sex <- data.table::dcast(followup_vax1_sex, Datasource ~ sex, value.var = c("sex_value"))
   followup_vax2_sex <- unique(followup_vax2[, .(Datasource, sex, sex_value)])
   followup_vax2_sex <- data.table::dcast(followup_vax2_sex, Datasource ~ sex, value.var = c("sex_value"))
   
   followup_vax1_cohort <- unique(followup_vax1[, .(Datasource, ageband_at_date_vax_1, cohort_value)])
   empty_vax1_cohort <- data.table(Datasource = thisdatasource, ageband_at_date_vax_1 = paste0(vect_recode, "_vax1"))
   followup_vax1_cohort <- merge(empty_vax1_cohort, followup_vax1_cohort, all.x = T)
   older60_vax1 <- copy(followup_vax1_cohort)[ageband_at_date_vax_1 %in% paste0("Followup_", Agebands60, "_vax1"), sum(cohort_value, na.rm = T)]
   older60_vax1 <- data.table::data.table(Datasource = thisdatasource, ageband_at_date_vax_1 = "Followup_60+_vax1",
                                          cohort_value = older60_vax1)
   followup_vax1_cohort <- rbind(followup_vax1_cohort, older60_vax1)
   followup_vax1_cohort <- data.table::dcast(followup_vax1_cohort, Datasource ~ ageband_at_date_vax_1, value.var = c("cohort_value"))
   followup_vax2_cohort <- unique(followup_vax2[, .(Datasource, ageband_at_date_vax_1, cohort_value)])
   empty_vax2_cohort <- data.table(Datasource = thisdatasource,
                                   ageband_at_date_vax_1 = paste0(vect_recode, "_vax2"))
   followup_vax2_cohort <- merge(empty_vax2_cohort, followup_vax2_cohort, all.x = T)
   older60_vax2 <- copy(followup_vax2_cohort)[ageband_at_date_vax_1 %in% paste0("Followup_", Agebands60, "_vax2"), sum(cohort_value, na.rm = T)]
   older60_vax2 <- data.table::data.table(Datasource = thisdatasource, ageband_at_date_vax_1 = "Followup_60+_vax2",
                                          cohort_value = older60_vax2)
   followup_vax2_cohort <- rbind(followup_vax2_cohort, older60_vax2)
   followup_vax2_cohort <- data.table::dcast(followup_vax2_cohort, Datasource ~ ageband_at_date_vax_1, value.var = c("cohort_value"))
   
   followup_vax1_complete <- unique(followup_vax1[, .(Datasource, Followup_total_vax1)])
   followup_vax2_complete <- unique(followup_vax2[, .(Datasource, Followup_total_vax2)])
   
   D4_followup_from_vax <- Reduce(merge, list(followup_vax1_sex, followup_vax2_sex, followup_vax1_cohort,
                                              followup_vax2_cohort, followup_vax1_complete, followup_vax2_complete))
   
   vect_col_to_year <- c("Followup_males_vax1", "Followup_males_vax2", "Followup_females_vax1", "Followup_females_vax2",
                         "Followup_total_vax1", "Followup_total_vax2", paste0(vect_recode, "_vax1"), "Followup_60+_vax1",
                         paste0(vect_recode, "_vax2"), "Followup_60+_vax2")
   
   tot_col <- c("Datasource", vect_col_to_year)
   
   D4_followup_from_vax <- D4_followup_from_vax[, lapply(.SD, as.numeric), .SDcols = vect_col_to_year]
   D4_followup_from_vax <- D4_followup_from_vax[, (vect_col_to_year) := round(.SD / 365.25), .SDcols = vect_col_to_year]
   D4_followup_from_vax <- D4_followup_from_vax[, Datasource := thisdatasource]
   
   D4_followup_from_vax <- D4_followup_from_vax[, ..tot_col]
   D4_followup_from_vax <- na_to_0(D4_followup_from_vax)
   
   nameoutput <- paste0("D4_followup_from_vax")
   assign(nameoutput, D4_followup_from_vax)
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   D4_distance_doses <- Vaccin_cohort[, .(person_id, date_vax1, date_vax2, type_vax_1, type_vax_2)]
   D4_distance_doses <- D4_distance_doses[!is.na(date_vax2) & type_vax_1 == type_vax_2, ]
   D4_distance_doses <- D4_distance_doses[, distance := correct_difftime(date_vax2 - 1, date_vax1)]
   
   
   theme_set(theme_bw())
   
   # Plot
   g <- ggplot(D4_distance_doses, aes(as.numeric(distance), fill = factor(type_vax_1))) +
      geom_histogram(aes(y=..density..), alpha=0.5, 
                     position="identity", binwidth = 1) +
      geom_density(alpha=0.2) + 
      scale_x_continuous(limits = c(as.numeric(D4_distance_doses[, min(distance)]) - 5,
                                    as.numeric(D4_distance_doses[, max(distance)]) + 5),
                         oob = scales::oob_keep) +
      labs(title = "Density plot of distances", 
           subtitle = "Distances between First and Second Doses Grouped by Manufacturers",
           caption = "Source: D4_distance_doses",
           x = "Distance",
           fill = "Manufacturers")
   suppressMessages(ggsave(paste0(direxp, "Density_plot_distance_doses.png"), plot = g,
                           units = c("cm"), dpi = 600))
   
   g <- ggplot(D4_distance_doses, aes(as.numeric(distance), fill = factor(type_vax_1))) +
      geom_histogram(alpha=0.5, 
                     position="identity", binwidth = 1) +
      geom_density(alpha=0.2) + 
      scale_x_continuous(limits = c(as.numeric(D4_distance_doses[, min(distance)]) - 5,
                                    as.numeric(D4_distance_doses[, max(distance)]) + 5),
                         oob = scales::oob_keep) +
      labs(title = "Density plot of distances", 
           subtitle = "Distances between First and Second Doses Grouped by Manufacturers",
           caption = "Source: D4_distance_doses",
           x = "Distance",
           fill = "Manufacturers")
   suppressMessages(ggsave(paste0(direxp, "Histogram_distance_doses.png"), plot = g,
                           units = c("cm"), dpi = 600))
   
   
   D4_distance_doses <- D4_distance_doses[,c("P25", "P50", "p75") :=
                                             as.list(round(quantile(distance, probs = c(0.25, 0.50, 0.75)), 0)), by = "type_vax_1"]
   
   empty_distances <- data.table(Datasource = thisdatasource, type_vax_1 = c("Moderna", "Pfizer", "AstraZeneca"))
   
   D4_distance_doses <- unique(D4_distance_doses[, .(type_vax_1, P25, P50, p75)][, Datasource := thisdatasource])
   D4_distance_doses <- merge(empty_distances, D4_distance_doses, all.x = T)
   D4_distance_doses <- data.table::dcast(D4_distance_doses, Datasource ~ type_vax_1, value.var = c("P25", "P50", "p75"))
   setnames(D4_distance_doses,
            c("P25_AstraZeneca", "P25_Moderna", "P25_Pfizer", "P50_AstraZeneca", "P50_Moderna", "P50_Pfizer",
              "p75_AstraZeneca", "p75_Moderna", "p75_Pfizer"),
            c("Distance_P25_AZ_1_2", "Distance_P25_Moderna1_2", "Distance_P25_Pfizer1_2", "Distance_P50_AZ_1_2",
              "Distance_P50_Moderna1_2", "Distance_P50_Pfizer1_2", "Distance_P75_AZ_1_2", "Distance_P75_Moderna1_2",
              "Distance_P75_Pfizer1_2"))
   D4_distance_doses <- D4_distance_doses[, c("Datasource", "Distance_P25_Pfizer1_2", "Distance_P50_Pfizer1_2", "Distance_P75_Pfizer1_2",
                                              "Distance_P25_Moderna1_2", "Distance_P50_Moderna1_2", "Distance_P75_Moderna1_2",
                                              "Distance_P25_AZ_1_2", "Distance_P50_AZ_1_2", "Distance_P75_AZ_1_2")]
   
   nameoutput <- paste0("D4_distance_doses")
   assign(nameoutput, D4_distance_doses)
   
   
   fwrite(get(nameoutput), file = paste0(dirD4tables, nameoutput,".csv"))
   rm(list=nameoutput)
   
   
   rm(study_population, Vaccin_cohort, dec31, jan1, D4_followup_sex, D4_followup_cohort,
      D4_followup_complete, followup_vax1, followup_vax2, followup_vax1_sex,
      followup_vax2_sex, followup_vax1_cohort, empty_vax1_cohort, followup_vax2_cohort, empty_vax2_cohort,
      followup_vax1_complete, followup_vax2_complete, empty_distances, study_population_cov_ALL,g)
}
