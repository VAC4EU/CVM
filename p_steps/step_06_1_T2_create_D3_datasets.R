# Create D3 cohort datasets. They are temporary since they don't have risks at vaccination
#-----------------------------------------------
# input: D4_study_population, selected_doses, D3_study_population_cov_ALL
# output: D3_study_population_no_risk, D3_Vaccin_cohort_no_risk

for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]])) 
  
load(paste0(dirtemp, "selected_doses.RData"))
load(paste0(diroutput,"D3_study_population_cov_ALL",suffix[[subpop]],".RData"))
study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL", suffix[[subpop]]))


selected_doses <- selected_doses[, .(person_id, date, vx_dose, vx_manufacturer)]
D3_doses <- merge(study_population, selected_doses, all.x = T, by="person_id")
setnames(D3_doses, "start_follow_up", "start_lookback")

D3_doses[, vx_dose := as.character(vx_dose)]

D3_doses <- data.table::dcast(D3_doses, person_id + sex + date_of_birth + date_of_death + study_entry_date + start_lookback +
                                study_exit_date ~ vx_dose, value.var = c("date", "vx_manufacturer"))

setnames(D3_doses, c("date_1", "date_2", "vx_manufacturer_1", "vx_manufacturer_2"),
         c("date_vax1", "date_vax2", "type_vax_1", "type_vax_2"))

D3_doses <- D3_doses[study_exit_date < date_vax1, c("date_vax1", "type_vax_1") := NA]
D3_doses <- D3_doses[study_exit_date < date_vax2, c("date_vax2", "type_vax_2") := NA]
D3_doses <- D3_doses[!is.na(date_vax1), c("study_entry_date_vax1", "study_exit_date_vax1") := list(date_vax1, fifelse(is.na(date_vax2), study_exit_date, date_vax2 - 1))]
D3_doses <- D3_doses[!is.na(date_vax2), c("study_entry_date_vax2", "study_exit_date_vax2") := list(date_vax2, study_exit_date)]

D3_doses <- D3_doses[, age_at_study_entry := floor(lubridate::time_length(correct_difftime(study_entry_date, date_of_birth), "years"))]
D3_doses <- D3_doses[, ageband_at_study_entry := cut(age_at_study_entry, breaks = Agebands, labels = Agebands_labels)]

D3_doses <- D3_doses[, age_at_1_jan_2021 := floor(lubridate::time_length(correct_difftime(firstjan2021, date_of_birth), "years"))]
D3_doses <- D3_doses[, ageband_at_1_jan_2021 := cut(age_at_1_jan_2021, breaks = Agebands, labels = Agebands_labels)]

D3_doses <- D3_doses[, age_at_date_vax_1 := floor(lubridate::time_length(correct_difftime(date_vax1, date_of_birth), "years"))]
D3_doses <- D3_doses[, ageband_at_date_vax_1 := cut(age_at_date_vax_1, breaks = Agebands, labels = Agebands_labels)]

D3_doses <- D3_doses[, fup_days := correct_difftime(study_exit_date, study_entry_date)]
D3_doses <- D3_doses[, fup_no_vax := fifelse(is.na(study_entry_date_vax1), fup_days, correct_difftime(study_entry_date_vax1 - 1, study_entry_date))]
D3_doses <- D3_doses[!is.na(study_entry_date_vax1), fup_vax1 := correct_difftime(study_exit_date_vax1, study_entry_date_vax1)]
D3_doses <- D3_doses[!is.na(study_entry_date_vax2), fup_vax2 := correct_difftime(study_exit_date_vax2, study_entry_date_vax2)]

study_population_cov_ALL <- study_population_cov_ALL[, .(person_id, CV_either_DX_or_DP, COVCANCER_either_DX_or_DP,
                                                               COVCOPD_either_DX_or_DP, COVHIV_either_DX_or_DP,
                                                               COVCKD_either_DX_or_DP, COVDIAB_either_DX_or_DP,
                                                               COVOBES_either_DX_or_DP, COVSICKLE_either_DX_or_DP,
                                                               IMMUNOSUPPR_at_study_entry, all_covariates_non_CONTR)]

setnames(study_population_cov_ALL,
         c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
           "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP",
           "IMMUNOSUPPR_at_study_entry", "all_covariates_non_CONTR"),
         c("CV_at_study_entry", "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
           "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
           "immunosuppressants_at_study_entry", "at_risk_at_study_entry"))

D3_doses <- merge(D3_doses, study_population_cov_ALL, all.x = T, by="person_id")

study_population_no_risk <- D3_doses[, .(person_id, sex, date_of_birth, start_lookback, study_entry_date, study_exit_date,
                                    date_vax1, date_vax2, age_at_study_entry, ageband_at_study_entry, age_at_1_jan_2021,
                                    ageband_at_1_jan_2021, CV_at_study_entry, COVCANCER_at_study_entry,
                                    COVCOPD_at_study_entry, COVHIV_at_study_entry, COVCKD_at_study_entry,
                                    COVDIAB_at_study_entry, COVOBES_at_study_entry, COVSICKLE_at_study_entry,
                                    immunosuppressants_at_study_entry, at_risk_at_study_entry, age_at_date_vax_1,
                                    ageband_at_date_vax_1, type_vax_1, type_vax_2, study_entry_date_vax1, study_exit_date_vax1,
                                    study_entry_date_vax2, study_exit_date_vax2, fup_days, fup_no_vax, fup_vax1, fup_vax2)]

tempname<-paste0("D3_study_population_no_risk",suffix[[subpop]])
assign(tempname,study_population_no_risk)
save(list=tempname,file=paste0(dirtemp,tempname,".RData"))

Vaccin_cohort_no_risk <- get(tempname)[!is.na(date_vax1) & (is.na(date_vax2) | date_vax1 < date_vax2), ]
Vaccin_cohort_no_risk <- Vaccin_cohort_no_risk[, age_at_date_vax_2 := floor(time_length(correct_difftime(date_vax2, date_of_birth), "years"))]
Vaccin_cohort_no_risk <- Vaccin_cohort_no_risk[, .(person_id, sex, date_of_birth, study_entry_date, study_exit_date, date_vax1,
                                         date_vax2, age_at_date_vax_1, ageband_at_date_vax_1, age_at_date_vax_2, type_vax_1, type_vax_2,
                                         study_entry_date_vax1, study_exit_date_vax1, study_entry_date_vax2,
                                         study_exit_date_vax2, fup_vax1, fup_vax2)]

tempname<-paste0("D3_Vaccin_cohort_no_risk",suffix[[subpop]])
assign(tempname,Vaccin_cohort_no_risk)
save(list=tempname,file=paste0(dirtemp,tempname,".RData"))

rm(list=paste0("D3_Vaccin_cohort_no_risk", suffix[[subpop]]))
rm(list=paste0("D3_study_population_no_risk", suffix[[subpop]]))
rm(list=paste0("D4_study_population", suffix[[subpop]]))
rm(list=paste0("D3_study_population_cov_ALL", suffix[[subpop]]))
}

rm(selected_doses, D3_doses,study_population,study_population_cov_ALL,Vaccin_cohort_no_risk,study_population_no_risk)

