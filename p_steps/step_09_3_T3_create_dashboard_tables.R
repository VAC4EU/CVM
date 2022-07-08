#-----------------------------------------------
# Aggregate persontime datasets for poisson

# input: D3_vaxweeks, D3_Vaccin_cohort, D3_study_population, D4_doses_weeks, D3_study_population_cov_ALL, RES_IR_benefit_week_BC, RES_IR_benefit_fup_BC, RES_IR_benefit_week_RF, RES_IR_benefit_fup_RF, RES_IR_risk_week_BC, RES_IR_risk_fup_BC, RES_IR_risk_week_RF, RES_IR_risk_fup_RF
# output: DOSES_BIRTHCOHORTS, Intermediate coverage, COVERAGE_BIRTHCOHORTS, DOSES_RISKFACTORS, COVERAGE_RISKFACTORS, BENEFIT_BIRTHCOHORTS_CALENDARTIME, BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION, BENEFIT_RISKFACTORS_CALENDARTIME, BENEFIT_RISKFACTORS_TIMESINCEVACCINATION, RISK_BIRTHCOHORTS_CALENDARTIME, RISK_BIRTHCOHORTS_TIMESINCEVACCINATION, RISK_RISKFACTORS_CALENDARTIME, RISK_RISKFACTORS_TIMESINCEVACCINATION

for (subpop in subpopulations_non_empty) {  
   print(subpop)
   
   thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
   
   
   if(this_datasource_has_subpopulations == T)  dirdashboard <- paste0(thisdirexp,"dashboard tables/")
   
   suppressWarnings(if (!file.exists(dirdashboard)) dir.create(file.path(dirdashboard)))
   
   
   load(paste0(dirtemp,"D3_vaxweeks",suffix[[subpop]],".RData"))
   load(paste0(dirtemp,"D3_Vaccin_cohort",suffix[[subpop]],".RData"))
   load(paste0(dirtemp,"D3_study_population",suffix[[subpop]],".RData"))
   load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))
   load(paste0(thisdirexp,"D4_doses_weeks.RData"))
   
   vaxweeks<-get(paste0("D3_vaxweeks",suffix[[subpop]]))
   rm(list=paste0("D3_vaxweeks",suffix[[subpop]]))
   Vaccin_cohort<-get(paste0("D3_Vaccin_cohort",suffix[[subpop]]))
   rm(list=paste0("D3_Vaccin_cohort",suffix[[subpop]]))
   study_population<-get(paste0("D3_study_population",suffix[[subpop]]))
   rm(list=paste0("D3_study_population",suffix[[subpop]]))
   list_outcomes_obs<-get(paste0("list_outcomes_observed",suffix[[subpop]]))
   rm(list=paste0("list_outcomes_observed",suffix[[subpop]]))
   doses_weeks<-get(paste0("D4_doses_weeks"))
   rm(list=paste0("D4_doses_weeks"))
   
   # Birth Cohort ----------------------------------------------------------------------------------------------------
   
   cohort_to_doses_weeks <- study_population[, .(person_id, sex, type_vax_1, type_vax_2, ageband_at_study_entry)]
   cohort_to_doses_weeks <- cohort_to_doses_weeks[!is.na(type_vax_1), ]
   
   all_mondays <- seq.Date(as.Date("19000101","%Y%m%d"), study_end, by = "week")
   
   monday_week <- seq.Date(from = find_last_monday(study_start, all_mondays), to = find_last_monday(study_end, all_mondays),
                           by = "week")
   double_weeks <- data.table(weeks_to_join = monday_week, monday_week = monday_week)
   all_days_df <- data.table(all_days = seq.Date(from = find_last_monday(study_start, monday_week), to = study_end, by = "days"))
   all_days_df <- merge(all_days_df, double_weeks, by.x = "all_days", by.y = "weeks_to_join", all.x = T)
   all_days_df <- all_days_df[, monday_week := nafill(monday_week, type="locf")]
   all_days_df <- all_days_df[all_days >= study_start,]
   
   
   all_tuesdays <- seq.Date(as.Date("19000102","%Y%m%d"), study_end, by = "week")
   
   tuesday_week <- seq.Date(from = find_last_monday(study_start, all_tuesdays), to = find_last_monday(study_end, all_tuesdays),
                            by = "week")
   double_weeks <- data.table(weeks_to_join = tuesday_week, tuesday_week = tuesday_week)
   all_tuesdays_df <- data.table(all_days = seq.Date(from = find_last_monday(study_start, all_tuesdays), to = study_end, by = "days"))
   all_tuesdays_df <- merge(all_tuesdays_df, double_weeks, by.x = "all_days", by.y = "weeks_to_join", all.x = T)
   all_tuesdays_df <- all_tuesdays_df[, tuesday_week := nafill(tuesday_week, type="locf")]
   all_tuesdays_df <- all_tuesdays_df[all_days >= study_start,]
   
   
   
   
   exited_pop <- copy(vaxweeks)[, .(person_id, end_date_of_period, Dose, week)]
   exited_pop <- exited_pop[exited_pop[, .I[Dose == max(Dose)], by = "person_id"]$V1]
   exited_pop <- exited_pop[exited_pop[, .I[week == max(week)], by = "person_id"]$V1][, week := NULL]
   
   exited_pop_1 <- copy(exited_pop)[Dose == 2, ][, Dose := 1]
   exited_pop <- rbind(exited_pop, exited_pop_1)
   
   exited_pop <- merge(exited_pop, all_tuesdays_df, by.x = "end_date_of_period",
                       by.y = "all_days", all.x = T)[, end_date_of_period := NULL]
   setnames(exited_pop,  "tuesday_week", "week")
   exited_pop <- exited_pop[, week := week + 6]
   
   pop_traits <- copy(study_population)[, .(person_id, ageband_at_study_entry, CV_at_study_entry, COVCANCER_at_study_entry,
                                            COVCOPD_at_study_entry, COVHIV_at_study_entry, COVCKD_at_study_entry,
                                            COVDIAB_at_study_entry, COVOBES_at_study_entry, COVSICKLE_at_study_entry,
                                            immunosuppressants_at_study_entry, at_risk_at_study_entry, type_vax_1,
                                            type_vax_2)]
   pop_traits <- pop_traits[!(is.na(type_vax_1) & is.na(type_vax_2)), ]
   pop_traits <- melt(pop_traits, id.vars = c("person_id", "ageband_at_study_entry", "type_vax_1", "type_vax_2"),
                      measure.vars = c("CV_at_study_entry", "COVCANCER_at_study_entry",
                                       "COVCOPD_at_study_entry", "COVHIV_at_study_entry", "COVCKD_at_study_entry",
                                       "COVDIAB_at_study_entry", "COVOBES_at_study_entry", "COVSICKLE_at_study_entry",
                                       "immunosuppressants_at_study_entry", "at_risk_at_study_entry"),
                      variable.name = "riskfactor")
   pop_traits <- melt(pop_traits, id.vars = c("person_id", "ageband_at_study_entry", "riskfactor", "value"),
                      measure.vars = c("type_vax_1", "type_vax_2"),
                      variable.name = "Dose", value.name = "vx_manufacturer")
   pop_traits <- pop_traits[, Dose := fifelse(Dose == "type_vax_1", "1", "2")]
   
   exited_pop <- merge(exited_pop, pop_traits, by = c("person_id", "Dose"), all.x = T)
   setnames(exited_pop, "Dose", "dose")
   
   exited_pop_birth_cohorts <- unique(copy(exited_pop)[, c("riskfactor", "value") := NULL])
   exited_pop_birth_cohorts <- exited_pop_birth_cohorts[, .N, by = c("dose", "week", "ageband_at_study_entry", "vx_manufacturer")]
   
   exited_pop_risk_factors <- unique(copy(exited_pop)[value == 1, ][, c("ageband_at_study_entry", "value") := NULL])
   exited_pop_risk_factors <- exited_pop_risk_factors[, .N, by = c("dose", "week", "riskfactor", "vx_manufacturer")]
   
   all_pop <- copy(exited_pop_birth_cohorts)[ageband_at_study_entry %in% Agebands_labels, .(N = sum(N)),
                                             by = c("dose", "week", "vx_manufacturer")][, ageband_at_study_entry := "all_agebands"]
   exited_pop_birth_cohorts <- rbind(exited_pop_birth_cohorts, all_pop)
   
   exited_pop_birth_cohorts <- bc_divide_60(exited_pop_birth_cohorts, c("week", "vx_manufacturer", "dose"), "N", only_old = T)
   
   setorder(exited_pop_birth_cohorts, week)
   exited_pop_birth_cohorts <- exited_pop_birth_cohorts[, exited := cumsum(N),
                                                        by = c("dose", "ageband_at_study_entry", "vx_manufacturer")][, N := NULL]
   exited_pop_birth_cohorts <- exited_pop_birth_cohorts[, week := format(week, "%Y%m%d")]
   setnames(exited_pop_birth_cohorts, "ageband_at_study_entry", "ageband")
   
   setorder(exited_pop_risk_factors, week)
   exited_pop_risk_factors <- exited_pop_risk_factors[, exited := cumsum(N),
                                                      by = c("dose", "riskfactor", "vx_manufacturer")][, N := NULL]
   exited_pop_risk_factors <- exited_pop_risk_factors[, week := format(week, "%Y%m%d")]
   
   
   
   vaxweeks <- vaxweeks[week == 0]
   vaxweeks <- merge(vaxweeks, all_days_df, by.x = "start_date_of_period", by.y = "all_days", all.x = T)
   
   vaxweeks_to_dos_bir_cor_base <- merge(vaxweeks, cohort_to_doses_weeks, by = "person_id")
   
   vaxweeks_to_dos_bir_cor <- vaxweeks_to_dos_bir_cor_base[, vx_manufacturer := fifelse(Dose == 1, type_vax_1, type_vax_2)]
   vaxweeks_to_dos_bir_cor <- vaxweeks_to_dos_bir_cor[, Datasource := thisdatasource]
   
   vaxweeks_to_dos_bir_cor <- vaxweeks_to_dos_bir_cor[, .(Datasource, monday_week, vx_manufacturer, Dose, ageband_at_study_entry)]
   
   vaxweeks_to_dos_bir_cor <- vaxweeks_to_dos_bir_cor[, .(N = .N), by = c("Datasource", "monday_week", "vx_manufacturer", "Dose", "ageband_at_study_entry")]
   
   vaxweeks_to_dos_bir_cor <- bc_divide_60(vaxweeks_to_dos_bir_cor, c("Datasource", "monday_week", "vx_manufacturer", "Dose"),
                                           "N", only_old = T)
   
   all_ages <- copy(vaxweeks_to_dos_bir_cor)[ageband_at_study_entry %in% Agebands_labels, .(N = sum(N)), by = c("Datasource", "monday_week", "vx_manufacturer", "Dose")]
   all_ages <- all_ages[, ageband_at_study_entry := "all_agebands"]
   vaxweeks_to_dos_bir_cor <- rbind(vaxweeks_to_dos_bir_cor, all_ages)
   
   setnames(vaxweeks_to_dos_bir_cor, c("Datasource", "monday_week", "Dose", "ageband_at_study_entry"),
            c("datasource", "week", "dose", "ageband"))
   
   complete_df <- expand.grid(datasource = thisdatasource, week = monday_week, vx_manufacturer = c("Moderna", "Pfizer", "AstraZeneca", "J&J", "UKN"),
                              dose = c("1", "2"), ageband = c(Agebands_labels, "all_agebands", "60+"))
   
   vaxweeks_to_dos_bir_cor <- merge(vaxweeks_to_dos_bir_cor, complete_df, all.y = T, by = c("datasource", "week", "vx_manufacturer", "dose", "ageband"))
   DOSES_BIRTHCOHORTS <- vaxweeks_to_dos_bir_cor[is.na(N), N := 0][, week := format(week, "%Y%m%d")]
   
   nameoutput <- paste0("DOSES_BIRTHCOHORTS")
   assign(nameoutput, DOSES_BIRTHCOHORTS)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   
   
   # fwrite(DOSES_BIRTHCOHORTS, file = paste0(dirdashboard, "DOSES_BIRTHCOHORTS.csv"))
   # 
   # tot_pop_cohorts <- study_population[, birth_cohort := findInterval(year(date_of_birth), c(1940, 1950, 1960, 1970, 1980, 1990))]
   # tot_pop_cohorts$birth_cohort <- as.character(tot_pop_cohorts$birth_cohort)
   # tot_pop_cohorts <- tot_pop_cohorts[.(birth_cohort = c("0", "1", "2", "3", "4", "5", "6"),
   #                                      to = c("<1940", "1940-1949", "1950-1959", "1960-1969", "1970-1979",
   #                                             "1980-1989", "1990+")),
   #                                    on = "birth_cohort", birth_cohort := i.to]
   # tot_pop_cohorts <- tot_pop_cohorts[, .(pop_cohorts = .N), by = c("birth_cohort")]
   # all_pop <- unique(copy(tot_pop_cohorts)[, pop_cohorts := sum(pop_cohorts)][, birth_cohort := "all_birth_cohorts"])
   # tot_pop_cohorts <- rbind(tot_pop_cohorts, all_pop)
   # older60 <- copy(tot_pop_cohorts)[birth_cohort %in% c("<1940", "1940-1949", "1950-1959"), sum(pop_cohorts)]
   # older60 <- data.table::data.table(birth_cohort = "<1960", pop_cohorts = older60)
   # tot_pop_cohorts <- rbind(tot_pop_cohorts, older60)
   
   
   tot_pop_cohorts <- unique(doses_weeks[, .(week = format(Week_number, "%Y%m%d"), ageband = ageband_at_study_entry,
                                             Persons_in_week)])
   all_pop <- copy(tot_pop_cohorts)[ageband %in% Agebands_labels, .(Persons_in_week = sum(Persons_in_week)),
                                    by = "week"][, ageband := "all_agebands"]
   tot_pop_cohorts <- rbind(tot_pop_cohorts, all_pop)
   
   
   COVERAGE_BIRTHCOHORTS <- merge(DOSES_BIRTHCOHORTS, tot_pop_cohorts, by = c("week", "ageband"), all.x = T)
   setorder(COVERAGE_BIRTHCOHORTS, week)
   
   COVERAGE_BIRTHCOHORTS <- COVERAGE_BIRTHCOHORTS[, cum_N := cumsum(N), by = c("datasource", "vx_manufacturer", "dose", "ageband")]
   
   COVERAGE_BIRTHCOHORTS <- merge(COVERAGE_BIRTHCOHORTS, exited_pop_birth_cohorts,
                                  by = c("week", "vx_manufacturer", "dose", "ageband"), all.x = T)
   COVERAGE_BIRTHCOHORTS <- COVERAGE_BIRTHCOHORTS[is.na(exited), exited := 0][, cum_N := cum_N - exited][, exited := NULL]
   
   COVERAGE_BIRTHCOHORTS <- COVERAGE_BIRTHCOHORTS[, percentage := round(cum_N / Persons_in_week  * 100, 1)]
   
   rm(list=nameoutput)
   nameoutput <- paste0("Intermediate coverage",suffix[[subpop]])
   assign(nameoutput, COVERAGE_BIRTHCOHORTS)
   save(nameoutput,file=paste0(dirtemp,nameoutput,".RData"),list=nameoutput)
   COVERAGE_BIRTHCOHORTS <- COVERAGE_BIRTHCOHORTS[, .(datasource, week, vx_manufacturer, dose, ageband, percentage)]
   
   rm(list=nameoutput)
   nameoutput <- paste0("COVERAGE_BIRTHCOHORTS")
   assign(nameoutput, COVERAGE_BIRTHCOHORTS)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   
   if (any(COVERAGE_BIRTHCOHORTS[, percentage] > 100)) {
      stop("Percentage in COVERAGE_BIRTHCOHORTS > 100")
   }
   
   rm(list=nameoutput)
   
   # Risk Factors ----------------------------------------------------------------------------------------------------
   
   load(paste0(diroutput,"D3_study_population_cov_ALL",suffix[[subpop]],".RData"))
   
   study_population_cov_ALL<-get(paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
   rm(list=paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
   
   
   setnames(study_population_cov_ALL,
            c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP", "COVHIV_either_DX_or_DP",
              "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP", "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP",
              "IMMUNOSUPPR_at_study_entry", "all_covariates_non_CONTR"),
            c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
              "any_risk_factors"))
   
   study_population_cov_ALL <- study_population_cov_ALL[, .(person_id, CV, COVCANCER, COVCOPD, COVHIV, COVCKD,
                                                            COVDIAB, COVOBES, COVSICKLE, IMMUNOSUPPR, any_risk_factors)]
   
   study_population_cov_ALL <- melt(study_population_cov_ALL,
                                    measure.vars = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                     "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                    variable.name = "riskfactor", value.name = "to_drop")
   
   study_population_cov_ALL <- study_population_cov_ALL[to_drop == 1, ]
   vaxweeks_to_dos_risk <- merge(vaxweeks_to_dos_bir_cor_base, study_population_cov_ALL, by = "person_id")
   
   vaxweeks_to_dos_risk <- vaxweeks_to_dos_risk[, vx_manufacturer := fifelse(Dose == 1, type_vax_1, type_vax_2)]
   vaxweeks_to_dos_risk <- vaxweeks_to_dos_risk[, Datasource := thisdatasource]
   
   vaxweeks_to_dos_risk <- vaxweeks_to_dos_risk[, .(Datasource, monday_week, vx_manufacturer, Dose, riskfactor)]
   
   setnames(vaxweeks_to_dos_risk, c("Datasource", "monday_week", "Dose"), c("datasource", "week", "dose"))
   
   vaxweeks_to_dos_risk <- vaxweeks_to_dos_risk[, .(N = .N), by = c("datasource", "week", "vx_manufacturer",
                                                                    "dose", "riskfactor")]
   
   complete_df <- expand.grid(datasource = thisdatasource, week = monday_week, vx_manufacturer = c("Moderna", "Pfizer", "AstraZeneca", "J&J", "UKN"),
                              dose = c("1", "2"), riskfactor = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                                 "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"))
   
   vaxweeks_to_dos_risk <- merge(vaxweeks_to_dos_risk, complete_df, all.y = T, by = c("datasource", "week", "vx_manufacturer", "dose", "riskfactor"))
   DOSES_RISKFACTORS <- vaxweeks_to_dos_risk[is.na(N), N := 0][, week := format(week, "%Y%m%d")]
   
   nameoutput <- paste0("DOSES_RISKFACTORS")
   assign(nameoutput, DOSES_RISKFACTORS)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   
   
   tot_pop_cohorts <- study_population_cov_ALL[, .(pop_cohorts = .N), by = c("riskfactor")]
   COVERAGE_RISKFACTORS <- merge(DOSES_RISKFACTORS, tot_pop_cohorts, by = "riskfactor", all.x = T)
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[is.na(pop_cohorts), pop_cohorts := 0]
   setorder(COVERAGE_RISKFACTORS, week)
   
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[, cum_N := cumsum(N), by = c("datasource", "vx_manufacturer", "dose", "riskfactor")]
   COVERAGE_RISKFACTORS <- merge(COVERAGE_RISKFACTORS, exited_pop_risk_factors,
                                 by = c("week", "vx_manufacturer", "dose", "riskfactor"), all.x = T)
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[is.na(exited), exited := 0][, cum_N := cum_N - exited][, exited := NULL]
   
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[, percentage := round(cum_N / pop_cohorts * 100, 1)]
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[is.nan(percentage), percentage := 0]
   COVERAGE_RISKFACTORS <- COVERAGE_RISKFACTORS[, .(datasource, week, vx_manufacturer, dose, riskfactor, percentage)]
   
   rm(list=nameoutput)
   nameoutput <- paste0("COVERAGE_RISKFACTORS")
   assign(nameoutput, COVERAGE_RISKFACTORS)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   
   if (any(COVERAGE_RISKFACTORS[, percentage] > 100)) {
      stop("Percentage in COVERAGE_RISKFACTORS > 100")
   }
   
   rm(list=nameoutput)
   
   rm(vaxweeks, cohort_to_doses_weeks, all_mondays, monday_week, double_weeks, all_days_df, vaxweeks_to_dos_bir_cor,
      all_ages, complete_df, study_population, tot_pop_cohorts, all_pop, Vaccin_cohort, study_population_cov_ALL, vaxweeks_to_dos_bir_cor_base, vaxweeks_to_dos_risk)
   
   # Benefit ------------------------------------------------------------------------------------------------------------
   thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
   
   IR_benefit_week<-fread(paste0(thisdirexp,"RES_IR_benefit_week_BC.csv"))
   
   BBC <- IR_benefit_week[, Dose := as.character(Dose)][ageband_at_study_entry != "0_59", ]
   BBC <- BBC[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   colA = paste("COVID_L", 1:5, "plus_b", sep = "")
   colB = paste("IR_COVID_L", 1:5, "plus", sep = "")
   colC = paste("lb_COVID_L", 1:5, "plus", sep = "")
   colD = paste("ub_COVID_L", 1:5, "plus", sep = "")
   
   BBC <- correct_col_type(BBC)
   
   BBC <- data.table::melt(BBC, measure = list(colA, colB, colC, colD), variable.name = "COVID",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   BBC <- BBC[is.na(ub), ub := 0]
   setnames(BBC, c("ageband_at_study_entry", "Dose", "type_vax"), c("ageband", "dose", "vx_manufacturer"))
   BBC <- BBC[, datasource := thisdatasource][sex == "both_sexes", ][, week := format(week, "%Y%m%d")]
   BBC <- BBC[, .(datasource, week, vx_manufacturer, dose, ageband, COVID, Numerator, IR, lb, ub)]
   vect_recode_COVID <- c("1" = "L1", "2" = "L2", "3" = "L3", "4" = "L4", "5" = "L5")
   BBC <- BBC[ , COVID := vect_recode_COVID[COVID]]
   
   nameoutput <- paste0("BENEFIT_BIRTHCOHORTS_CALENDARTIME")
   assign(nameoutput, BBC)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(BBC, IR_benefit_week)
   
   
   IR_benefit_fup<- fread(paste0(thisdirexp,"RES_IR_benefit_fup_BC.csv")) 
   
   
   BBT <- IR_benefit_fup[, Dose := as.character(Dose)][ageband_at_study_entry != "0_59", ]
   BBT <- BBT[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   colA = paste("COVID_L", 1:5, "plus_b", sep = "")
   colB = paste("IR_COVID_L", 1:5, "plus", sep = "")
   colC = paste("lb_COVID_L", 1:5, "plus", sep = "")
   colD = paste("ub_COVID_L", 1:5, "plus", sep = "")
   
   BBT <- correct_col_type(BBT)
   
   BBT <- data.table::melt(BBT, measure = list(colA, colB, colC, colD), variable.name = "COVID",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   BBT <- BBT[is.na(ub), ub := 0]
   setnames(BBT, c("ageband_at_study_entry", "Dose", "type_vax"), c("ageband", "dose", "vx_manufacturer"))
   BBT <- BBT[, datasource := thisdatasource][sex == "both_sexes", ]
   BBT <- BBT[, .(datasource, week_fup, vx_manufacturer, dose, ageband, COVID, Numerator, IR, lb, ub)]
   setnames(BBT, c("week_fup"), c("week_since_vaccination"))
   vect_recode_COVID <- c("1" = "L1", "2" = "L2", "3" = "L3", "4" = "L4", "5" = "L5")
   BBT <- BBT[ , COVID := vect_recode_COVID[COVID]]
   
   nameoutput <- paste0("BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION")
   assign(nameoutput, BBT)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(BBT, IR_benefit_fup)
   
   
   IR_benefit_week<-fread(paste0(thisdirexp,"RES_IR_benefit_week_RF.csv"))
   
   BRC <- IR_benefit_week[, Dose := as.character(Dose)]
   BRC <- BRC[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   colA = paste("COVID_L", 1:5, "plus_b", sep = "")
   colB = paste("IR_COVID_L", 1:5, "plus", sep = "")
   colC = paste("lb_COVID_L", 1:5, "plus", sep = "")
   colD = paste("ub_COVID_L", 1:5, "plus", sep = "")
   
   BRC <- correct_col_type(BRC)
   
   BRC <- data.table::melt(BRC, measure = list(colA, colB, colC, colD), variable.name = "COVID",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   BRC <- BRC[is.na(ub), ub := 0]
   setnames(BRC, c("Dose", "type_vax"), c("dose", "vx_manufacturer"))
   BRC <- BRC[, datasource := thisdatasource][sex == "both_sexes", ][, week := format(week, "%Y%m%d")]
   BRC <- BRC[, .(datasource, week, vx_manufacturer, dose, riskfactor, COVID, Numerator, IR, lb, ub)]
   vect_recode_COVID <- c("1" = "L1", "2" = "L2", "3" = "L3", "4" = "L4", "5" = "L5")
   BRC <- BRC[ , COVID := vect_recode_COVID[COVID]]
   
   nameoutput <- paste0("BENEFIT_RISKFACTORS_CALENDARTIME")
   assign(nameoutput, BRC)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(BRC, IR_benefit_week)
   
   IR_benefit_fup<-fread(paste0(thisdirexp,"RES_IR_benefit_fup_RF.csv"))
   
   
   BRT <- IR_benefit_fup[, Dose := as.character(Dose)]
   BRT <- BRT[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   colA = paste("COVID_L", 1:5, "plus_b", sep = "")
   colB = paste("IR_COVID_L", 1:5, "plus", sep = "")
   colC = paste("lb_COVID_L", 1:5, "plus", sep = "")
   colD = paste("ub_COVID_L", 1:5, "plus", sep = "")
   
   BRT <- correct_col_type(BRT)
   
   BRT <- data.table::melt(BRT, measure = list(colA, colB, colC, colD), variable.name = "COVID",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   BRT <- BRT[is.na(ub), ub := 0]
   setnames(BRT, c("Dose", "type_vax"), c("dose", "vx_manufacturer"))
   BRT <- BRT[, datasource := thisdatasource][sex == "both_sexes", ]
   BRT <- BRT[, .(datasource, week_fup, vx_manufacturer, dose, riskfactor, COVID, Numerator, IR, lb, ub)]
   setnames(BRT, c("week_fup"), c("week_since_vaccination"))
   vect_recode_COVID <- c("1" = "L1", "2" = "L2", "3" = "L3", "4" = "L4", "5" = "L5")
   BRT <- BRT[ , COVID := vect_recode_COVID[COVID]]
   
   nameoutput <- paste0("BENEFIT_RISKFACTORS_TIMESINCEVACCINATION")
   assign(nameoutput, BRT)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(BRT, vect_recode_COVID, IR_benefit_fup)
   
   
   
   
   # Risk ------------------------------------------------------------------------------------------------------------
   
   IR_risk_week<-fread(paste0(thisdirexp,"RES_IR_risk_week_BC.csv")) 
   
   RBC <- IR_risk_week[, Dose := as.character(Dose)][ageband_at_study_entry != "0_59", ]
   RBC <- RBC[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   list_risk <- list_outcomes_obs
   colA = paste0(list_risk, "_b")
   colB = paste0("IR_", list_risk)
   colC = paste0("lb_", list_risk)
   colD = paste0("ub_", list_risk)
   
   RBC <- correct_col_type(RBC)
   
   RBC <- data.table::melt(RBC, measure = list(colA, colB, colC, colD), variable.name = "AESI",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   setnames(RBC, c("ageband_at_study_entry", "Dose", "type_vax"), c("ageband", "dose", "vx_manufacturer"))
   RBC <- RBC[, datasource := thisdatasource][sex == "both_sexes", ][, week := format(week, "%Y%m%d")]
   RBC <- RBC[, .(datasource, week, vx_manufacturer, dose, ageband, AESI, Numerator, IR, lb, ub)]
   vect_recode_AESI <- list_risk
   names(vect_recode_AESI) <- c(as.character(seq_len(length(list_risk))))
   RBC <- RBC[ , AESI := vect_recode_AESI[AESI]]
   
   nameoutput <- paste0("RISK_BIRTHCOHORTS_CALENDARTIME")
   assign(nameoutput, RBC)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(RBC, IR_risk_week)
   
   
   IR_risk_fup<-fread(paste0(thisdirexp,"RES_IR_risk_fup_BC.csv")) 
   
   
   RBT <- IR_risk_fup[, Dose := as.character(Dose)][ageband_at_study_entry != "0_59", ]
   RBT <- RBT[Dose != "both_doses" & week_fup != "fup_until_4" & type_vax != "all_manufacturer", ]
   RBT <- RBT[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   list_risk <- list_outcomes_obs
   
   colA = paste0(list_risk, "_b")
   colB = paste0("IR_", list_risk)
   colC = paste0("lb_", list_risk)
   colD = paste0("ub_", list_risk)
   
   RBT <- correct_col_type(RBT)
   
   RBT <- data.table::melt(RBT, measure = list(colA, colB, colC, colD), variable.name = "AESI",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   RBT <- RBT[is.na(ub), ub := 0]
   setnames(RBT, c("ageband_at_study_entry", "Dose", "type_vax"), c("ageband", "dose", "vx_manufacturer"))
   RBT <- RBT[, datasource := thisdatasource][sex == "both_sexes", ]
   RBT <- RBT[, .(datasource, week_fup, vx_manufacturer, dose, ageband, AESI, Numerator, IR, lb, ub)]
   setnames(RBT, c("week_fup"), c("week_since_vaccination"))
   vect_recode_AESI <- list_risk
   names(vect_recode_AESI) <- c(as.character(seq_len(length(list_risk))))
   RBT <- RBT[ , AESI := vect_recode_AESI[AESI]]
   
   nameoutput <- paste0("RISK_BIRTHCOHORTS_TIMESINCEVACCINATION")
   assign(nameoutput, RBT)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(RBT, IR_risk_fup)
   
   
   IR_risk_week <- fread(paste0(thisdirexp,"RES_IR_risk_week_RF.csv"))
   RRC <- IR_risk_week[, Dose := as.character(Dose)]
   RRC <- RRC[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   list_risk <- list_outcomes_obs
   
   colA = paste0(list_risk, "_b")
   colB = paste0("IR_", list_risk)
   colC = paste0("lb_", list_risk)
   colD = paste0("ub_", list_risk)
   
   RRC <- correct_col_type(RRC)
   
   RRC <- data.table::melt(RRC, measure = list(colA, colB, colC, colD), variable.name = "AESI",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   setnames(RRC, c("Dose", "type_vax"), c("dose", "vx_manufacturer"))
   RRC <- RRC[, datasource := thisdatasource][sex == "both_sexes", ][, week := format(week, "%Y%m%d")]
   RRC <- RRC[, .(datasource, week, vx_manufacturer, dose, riskfactor, AESI, Numerator, IR, lb, ub)]
   vect_recode_AESI <- list_risk
   names(vect_recode_AESI) <- c(as.character(seq_len(length(list_risk))))
   RRC <- RRC[ , AESI := vect_recode_AESI[AESI]]
   
   nameoutput <- paste0("RISK_RISKFACTORS_CALENDARTIME")
   assign(nameoutput, RRC)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   rm(RRC, IR_risk_week)
   
   
   IR_risk_fup <- fread(paste0(thisdirexp,"RES_IR_risk_fup_RF.csv"))
   RRT <- IR_risk_fup[, Dose := as.character(Dose)]
   RRT <- RRT[Dose == 0, c("Dose", "type_vax") := list("no_dose", "none")]
   list_risk <- list_outcomes_obs
   
   colA = paste0(list_risk, "_b")
   colB = paste0("IR_", list_risk)
   colC = paste0("lb_", list_risk)
   colD = paste0("ub_", list_risk)
   
   RRT <- correct_col_type(RRT)
   
   RRT <- data.table::melt(RRT, measure = list(colA, colB, colC, colD), variable.name = "AESI",
                           value.name = c("Numerator", "IR", "lb", "ub"), na.rm = F)
   
   RRT <- RRT[is.na(ub), ub := 0]
   setnames(RRT, c("Dose", "type_vax"), c("dose", "vx_manufacturer"))
   RRT <- RRT[, datasource := thisdatasource][sex == "both_sexes", ]
   RRT <- RRT[, .(datasource, week_fup, vx_manufacturer, dose, riskfactor, AESI, Numerator, IR, lb, ub)]
   setnames(RRT, c("week_fup"), c("week_since_vaccination"))
   vect_recode_AESI <- list_risk
   names(vect_recode_AESI) <- c(as.character(seq_len(length(list_risk))))
   RRT <- RRT[ , AESI := vect_recode_AESI[AESI]]
   
   nameoutput <- paste0("RISK_RISKFACTORS_TIMESINCEVACCINATION")
   assign(nameoutput, RRT)
   fwrite(get(nameoutput), file = paste0(dirdashboard, nameoutput,".csv"))
   rm(list=nameoutput)
   
   rm(RRT, vect_recode_AESI, IR_risk_fup)
   
}
