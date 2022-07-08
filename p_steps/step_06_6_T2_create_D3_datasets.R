# ----------------------------------
# Create basic D3 datasets based on the population

# input: D3_study_population_no_risk, D3_study_population_cov_ALL , D3_Vaccin_cohort_cov_ALL
# output: D3_study_population, D3_Vaccin_cohort, D3_vaxweeks_vaccin_cohort, D3_studyweeks, D3_vaxweeks, D3_vaxweeks_including_not_vaccinated

for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D3_study_population_no_risk",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D3_study_population_cov_ALL",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D3_Vaccin_cohort_cov_ALL",suffix[[subpop]],".RData"))
  
  study_population_no_risk<-get(paste0("D3_study_population_no_risk", suffix[[subpop]]))
  rm(list=paste0("D3_study_population_no_risk", suffix[[subpop]]))
  
  study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL", suffix[[subpop]]))
  Vaccin_cohort_cov_ALL<-get(paste0("D3_Vaccin_cohort_cov_ALL", suffix[[subpop]]))
  rm(list=paste0("D3_study_population_cov_ALL", suffix[[subpop]]))
  rm(list=paste0("D3_Vaccin_cohort_cov_ALL", suffix[[subpop]]))
  
  cols_risk_factors = c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP",
                   "COVHIV_either_DX_or_DP", "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP",
                   "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP", "all_covariates_non_CONTR")
  
study_population_cov_ALL <- study_population_cov_ALL[, c("person_id", ..cols_risk_factors, "IMMUNOSUPPR_at_study_entry")]

setnames(study_population_cov_ALL, c(cols_risk_factors, "IMMUNOSUPPR_at_study_entry"),
         c("CV_at_study_entry", "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
           "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
           "at_risk_at_study_entry", "immunosuppressants_at_study_entry"))

Vaccin_cohort_cov_ALL <- Vaccin_cohort_cov_ALL[, c("person_id", ..cols_risk_factors, "IMMUNOSUPPR_at_date_vax_1")]

setnames(Vaccin_cohort_cov_ALL, c(cols_risk_factors, "IMMUNOSUPPR_at_date_vax_1"),
         c("CV_at_date_vax_1", "COVCANCER_at_date_vax_1", "COVCOPD_at_date_vax_1", "COVHIV_at_date_vax_1",
           "COVCKD_at_date_vax_1", "COVDIAB_at_date_vax_1", "COVOBES_at_date_vax_1", "COVSICKLE_at_date_vax_1",
           "at_risk_at_date_vax_1", "immunosuppressants_at_date_vax_1"))

D3_study_population <- merge(study_population_no_risk, Vaccin_cohort_cov_ALL, all.x = T, by="person_id")

D3_study_population <- D3_study_population[, .(person_id, sex, date_of_birth, start_lookback, study_entry_date,
                                               study_exit_date, date_vax1, date_vax2, age_at_study_entry,
                                               ageband_at_study_entry, age_at_1_jan_2021, ageband_at_1_jan_2021,
                                               CV_at_study_entry, COVCANCER_at_study_entry, COVCOPD_at_study_entry,
                                               COVHIV_at_study_entry, COVCKD_at_study_entry, COVDIAB_at_study_entry,
                                               COVOBES_at_study_entry, COVSICKLE_at_study_entry,
                                               immunosuppressants_at_study_entry, at_risk_at_study_entry,
                                               age_at_date_vax_1, ageband_at_date_vax_1, CV_at_date_vax_1,
                                               COVCANCER_at_date_vax_1, COVCOPD_at_date_vax_1, COVHIV_at_date_vax_1,
                                               COVCKD_at_date_vax_1, COVDIAB_at_date_vax_1, COVOBES_at_date_vax_1,
                                               COVSICKLE_at_date_vax_1, immunosuppressants_at_date_vax_1,
                                               at_risk_at_date_vax_1, type_vax_1, type_vax_2, study_entry_date_vax1,
                                               study_exit_date_vax1, study_entry_date_vax2, study_exit_date_vax2,
                                               fup_days, fup_no_vax, fup_vax1, fup_vax2)]

changeCols <- names(Filter(is.factor, D3_study_population))
D3_study_population <- D3_study_population[, (changeCols) := lapply(.SD, as.character), .SDcols = changeCols]

nameoutput<-paste0("D3_study_population",suffix[[subpop]])
assign(nameoutput,D3_study_population)
save(nameoutput, file = paste0(dirtemp, "D3_study_population",suffix[[subpop]],".RData"),list=nameoutput)

D3_Vaccin_cohort <- D3_study_population[!is.na(date_vax1) & (is.na(date_vax2) | date_vax1 < date_vax2), ]
D3_Vaccin_cohort <- D3_Vaccin_cohort[, age_at_date_vax_2 := floor(time_length(correct_difftime(date_vax2, date_of_birth), "years"))]
D3_Vaccin_cohort <- D3_Vaccin_cohort[, .(person_id, sex, date_of_birth, study_entry_date, study_exit_date, date_vax1,
                                         date_vax2, age_at_date_vax_1, ageband_at_date_vax_1, CV_at_date_vax_1,
                                         COVCANCER_at_date_vax_1, COVCOPD_at_date_vax_1, COVHIV_at_date_vax_1,
                                         COVCKD_at_date_vax_1, COVDIAB_at_date_vax_1, COVOBES_at_date_vax_1,
                                         COVSICKLE_at_date_vax_1, immunosuppressants_at_date_vax_1,
                                         at_risk_at_date_vax_1, age_at_date_vax_2, type_vax_1, type_vax_2,
                                         study_entry_date_vax1, study_exit_date_vax1, study_entry_date_vax2,
                                         study_exit_date_vax2, fup_vax1, fup_vax2)]

nameoutput<-paste0("D3_Vaccin_cohort",suffix[[subpop]])
assign(nameoutput,D3_Vaccin_cohort)
save(nameoutput, file = paste0(dirtemp, "D3_Vaccin_cohort",suffix[[subpop]],".RData"),list=nameoutput)

rm(list=paste0("D3_Vaccin_cohort", suffix[[subpop]]))

persons_at_risk <- copy(study_population_cov_ALL)[, .(person_id, at_risk_at_study_entry)]
D3_vaxweeks_vaccin_cohort <- D3_study_population[!is.na(date_vax1) & (is.na(date_vax2) | date_vax1 < date_vax2), ]
D3_vaxweeks_vaccin_cohort <- D3_vaxweeks_vaccin_cohort[, .(person_id, sex, ageband_at_study_entry, study_entry_date,
                                                           study_entry_date_vax1)]
D3_vaxweeks_vaccin_cohort <- merge(D3_vaxweeks_vaccin_cohort, persons_at_risk, by = "person_id", all.x = T)
D3_vaxweeks_vaccin_cohort <- D3_vaxweeks_vaccin_cohort[, study_exit_date := study_entry_date_vax1 - 1]
D3_vaxweeks_vaccin_cohort[, "study_entry_date_vax1" := NULL]

nameoutput<-paste0("D3_vaxweeks_vaccin_cohort",suffix[[subpop]])

assign(nameoutput,D3_vaxweeks_vaccin_cohort)
save(nameoutput, file = paste0(dirtemp, "D3_vaxweeks_vaccin_cohort",suffix[[subpop]],".RData"),list=nameoutput)
rm(list=nameoutput)

cohort_to_vaxweeks <- D3_study_population[, .(person_id, sex, study_entry_date, study_exit_date, ageband_at_study_entry,
                                              study_entry_date_vax1, study_exit_date_vax1, study_entry_date_vax2,
                                              study_exit_date_vax2, type_vax_1, type_vax_2, fup_no_vax, fup_vax1,
                                              fup_vax2, at_risk_at_study_entry)]

rm(list=paste0("D3_study_population", suffix[[subpop]]))

cohort_to_vaxweeks <- cohort_to_vaxweeks[!is.na(study_entry_date_vax1), study_exit_date := study_entry_date_vax1 - 1]
cols <- c("fup_no_vax", "fup_vax1", "fup_vax2")
cohort_to_vaxweeks <- cohort_to_vaxweeks[, (cols) := lapply(.SD, calc_precise_week), .SDcols = cols]

colA = c(paste("study_entry_date_vax", 1:2, sep = ""), "study_entry_date")
colB = c(paste("study_exit_date_vax", 1:2, sep = ""), "study_exit_date")
colC = c(paste("fup_vax", 1:2, sep = ""), "fup_no_vax")
colD = paste("type_vax", 1:2, sep = "_")
cohort_to_vaxweeks <- data.table::melt(cohort_to_vaxweeks, measure = list(colA, colB, colC, colD), variable.name = "Dose",
                                       value.name = c("study_entry_date", "study_exit_date", "fup", "type_vax"), na.rm = F)
cohort_to_vaxweeks <- cohort_to_vaxweeks[!is.na(study_entry_date) & !is.na(study_exit_date) & !is.na(fup), ]
cohort_to_vaxweeks <- cohort_to_vaxweeks[, Dose := as.character(Dose)][Dose == "3", Dose := "0"]
vaxweeks <- copy(cohort_to_vaxweeks)
studyweeks <- copy(cohort_to_vaxweeks)
vaxweeks <- vaxweeks[Dose %in% c(1, 2), ]
setnames(vaxweeks, c("study_entry_date", "study_exit_date", "fup"),
         c("study_entry_date_vax", "study_exit_date_vax", "fup_vax"))
vaxweeks <- as.data.table(lapply(vaxweeks, rep, vaxweeks$fup_vax))
vaxweeks[, week := rowid(person_id, Dose, study_entry_date_vax, study_exit_date_vax)]
vaxweeks <- vaxweeks[, week := week - 1][, fup_vax := fup_vax - 1]
vaxweeks <- vaxweeks[, start_date_of_period := study_entry_date_vax + 7 * week]
vaxweeks <- vaxweeks[, end_date_of_period := fifelse(week == fup_vax, study_exit_date_vax, start_date_of_period + 6)]
vaxweeks <- vaxweeks[, month := month(start_date_of_period)]
vaxweeks <- vaxweeks[, .(person_id, start_date_of_period, end_date_of_period, Dose, week, month, sex, type_vax, ageband_at_study_entry)]

vaxweeks_including_not_vaccinated <- copy(vaxweeks)[, month := NULL]

vaxweeks <- vaxweeks[, c("sex", "type_vax", "ageband_at_study_entry") := NULL]
vaxweeks <- vaxweeks[, .(person_id, start_date_of_period, end_date_of_period, Dose, week, month)]

tempname<-paste0("D3_vaxweeks",suffix[[subpop]])
assign(tempname,vaxweeks)
save(tempname, file = paste0(dirtemp, "D3_vaxweeks",suffix[[subpop]],".RData"),list=tempname)
rm(list=paste0("D3_vaxweeks",suffix[[subpop]]))

D3_studyweeks_not_vaccinated <- copy(studyweeks)[Dose == 0, ][, fup := 0]
D3_studyweeks_not_vaccinated <- D3_studyweeks_not_vaccinated[, .(person_id, study_entry_date, study_exit_date, Dose, fup, sex, type_vax, ageband_at_study_entry)]

tempname<-paste0("D3_studyweeks",suffix[[subpop]])
assign(tempname,studyweeks)
save(tempname, file = paste0(dirtemp, "D3_studyweeks",suffix[[subpop]],".RData"),list=tempname)
rm(list=paste0("D3_studyweeks",suffix[[subpop]]))

vaxweeks_including_not_vaccinated <- vaxweeks_including_not_vaccinated[, week := week + 1]

setnames(D3_studyweeks_not_vaccinated, c("study_entry_date", "study_exit_date", "fup"),
         c("start_date_of_period", "end_date_of_period", "week"))

vaxweeks_including_not_vaccinated <- rbind(vaxweeks_including_not_vaccinated, D3_studyweeks_not_vaccinated)
setnames(vaxweeks_including_not_vaccinated, "week", "week_fup")

setnames(study_population_cov_ALL,
         c("CV_at_study_entry", "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
           "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
           "immunosuppressants_at_study_entry", "at_risk_at_study_entry"),
         c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE",
           "IMMUNOSUPPR", "any_risk_factors"))

vaxweeks_including_not_vaccinated <- merge(vaxweeks_including_not_vaccinated, study_population_cov_ALL,
                                              all.x = T, by = "person_id", allow.cartesian = T)


tempname<-paste0("D3_vaxweeks_including_not_vaccinated",suffix[[subpop]])
assign(tempname,vaxweeks_including_not_vaccinated)
save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)

rm(list=paste0("D3_vaxweeks_including_not_vaccinated",suffix[[subpop]]))
}

rm(study_population_no_risk, Vaccin_cohort_cov_ALL, cohort_to_vaxweeks,
   colA, colB, colC, colD, studyweeks,vaxweeks,vaxweeks_including_not_vaccinated,
   D3_studyweeks_not_vaccinated, study_population_cov_ALL, persons_at_risk)

if(this_datasource_has_subpopulations==T) rm(D3_study_population,D3_vaxweeks_vaccin_cohort,D3_Vaccin_cohort)