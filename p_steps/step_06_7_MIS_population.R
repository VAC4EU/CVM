# ----------------------------------
# Calculate the cohort D3 datasets for MIS/myocard, For cohort c the dataset is temporary since it doesn't contains the risks at covid diagnosis.

# input: D3_events_ALL_OUTCOMES, D3_outcomes_covid , D3_study_population
# output: D3_study_variables_for_MIS, D4_population_b, D3_selection_criteria_c, D4_population_c_no_risk, D3_selection_criteria_d, D4_population_d

persontime_benefit_week <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  end_dic_2021<-as_date("20211231")
  first_jan_2020<-as_date("20200101")
  first_jan_2021<-as_date("20210101")
  
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_outcomes_covid",suffix[[subpop]],".RData")) #L1plus
  load(paste0(dirtemp,"D3_study_population",suffix[[subpop]],".RData"))
  
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  outcomes_covid<-get(paste0("D3_outcomes_covid", suffix[[subpop]]))
  rm(list=paste0("D3_outcomes_covid", suffix[[subpop]]))
  study_population<-get(paste0("D3_study_population", suffix[[subpop]]))
  rm(list=paste0("D3_study_population", suffix[[subpop]]))
  
  
  #add date of first covid to the population
  temp_covid<-unique(outcomes_covid[name_event=="COVID_L1plus",.(person_id,date_event)])
  temp_covid<-temp_covid[,min(date_event,na.rm = T),by="person_id"]
  
  setnames(temp_covid,"V1","covid_date")
  D3_study_variables_for_MIS <- merge(study_population, temp_covid, all.x = T, by="person_id")
  rm(temp_covid)
  
  # #add date of MIS broad to the population
  # temp_MIS_broad<-events_ALL_OUTCOMES[name_event=="MIS_broad",][,.(person_id,date_event)]
  # setnames(temp_MIS_broad,"date_event","MIS_date_broad")
  # D3_study_variables_for_MIS <- merge(D3_study_variables_for_MIS, temp_MIS_broad, all.x = T, by="person_id")
  # rm(temp_MIS_broad)
  # 
  # #add date of MIS narrow to the population
  # temp_MIS_narrow<-events_ALL_OUTCOMES[name_event=="MIS_narrow",][,.(person_id,date_event)]
  # setnames(temp_MIS_narrow,"date_event","MIS_date_narrow")
  # D3_study_variables_for_MIS <- merge(D3_study_variables_for_MIS, temp_MIS_narrow, all.x = T, by="person_id")
  # rm(temp_MIS_narrow)
  
  tempname<-paste0("D3_study_variables_for_MIS",suffix[[subpop]])
  assign(tempname,D3_study_variables_for_MIS)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  
  #COHORT B
  #add the study entry date for MIS
  D3_study_variables_for_MIS[,study_entry_date_MIS_b:=max(first_jan_2020,study_entry_date,na.rm = T),by="person_id"]
  #add the cohort entry date for MIS
  D3_study_variables_for_MIS[,cohort_entry_date_MIS_b:=study_entry_date_MIS_b,by="person_id"]
  #add the study exit date for MIS
  D3_study_variables_for_MIS[,study_exit_date_MIS_b:=min(end_dic_2021,study_exit_date,covid_date-1, date_vax1-1,na.rm = T),by="person_id"]
  # calculate correct fup_days
  D3_study_variables_for_MIS[, fup_days := correct_difftime(study_exit_date_MIS_b, cohort_entry_date_MIS_b)]
  #select the variables and save
  D4_population_b<-D3_study_variables_for_MIS[,.(person_id,sex,age_at_1_jan_2021,ageband_at_1_jan_2021,study_entry_date_MIS_b, cohort_entry_date_MIS_b, study_exit_date_MIS_b, fup_days)]
  
  tempname<-paste0("D4_population_b",suffix[[subpop]])
  assign(tempname,D4_population_b)
  save(tempname, file = paste0(diroutput, tempname,".RData"),list=tempname)
  rm(list=tempname)
  if (this_datasource_has_subpopulations==T) rm(D4_population_b)
  
  #---------------------------------
  #COHORT C
  #add the study entry date for MIS
  D3_study_variables_for_MIS[,study_entry_date_MIS_c:=first_jan_2020,by="person_id"]
  #add the cohort entry date for MIS
  D3_study_variables_for_MIS[,cohort_entry_date_MIS_c:=covid_date,by="person_id"]
  #add the study exit date for MIS
  D3_study_variables_for_MIS[,study_exit_date_MIS_c:=min(end_dic_2021,study_exit_date,date_vax1-1,na.rm = T),by="person_id"]
  # calculate correct fup_days
  D3_study_variables_for_MIS[, fup_days := correct_difftime(study_exit_date_MIS_c, cohort_entry_date_MIS_c)]
  #select the variables and save                           
  
  D3_selection_criteria_c <- D3_study_variables_for_MIS[is.na(covid_date) | study_exit_date_MIS_c <= cohort_entry_date_MIS_c, not_in_cohort_c:=1]
  D3_selection_criteria_c <- D3_selection_criteria_c[, not_in_cohort_c := replace(.SD, is.na(.SD), 0), .SDcols = "not_in_cohort_c"]
  
  tempname<-paste0("D3_selection_criteria_c",suffix[[subpop]])
  assign(tempname,D3_selection_criteria_c)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  
  D4_population_c_no_risk <- CreateFlowChart(
    dataset = D3_selection_criteria_c[,.(person_id,sex,age_at_1_jan_2021,ageband_at_1_jan_2021,study_entry_date_MIS_c, cohort_entry_date_MIS_c, study_exit_date_MIS_c,not_in_cohort_c, fup_days, CV_at_date_vax_1, COVCANCER_at_date_vax_1, COVCOPD_at_date_vax_1,
                                         COVHIV_at_date_vax_1, COVCKD_at_date_vax_1, COVDIAB_at_date_vax_1,
                                         COVOBES_at_date_vax_1, COVSICKLE_at_date_vax_1, immunosuppressants_at_date_vax_1,
                                         at_risk_at_date_vax_1)],
    listcriteria = c("not_in_cohort_c"),
    flowchartname = "Flowchart_cohort_c")
  
  tempname<-paste0("D4_population_c_no_risk",suffix[[subpop]])
  assign(tempname,D4_population_c_no_risk)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  rm(list=tempname)
  rm(D3_selection_criteria_c)
  if (this_datasource_has_subpopulations==T) rm(D4_population_c_no_risk)
  #---------------------------------
  #COHORT D
  
  #d crea: covid prima vaccinazione (usa in strato countpersontime) e e covdi dopo(usa per exit)
  
  #covid_d<-unique(outcomes_covid[name_event=="COVID_L1plus",.(person_id,date_event)])
  #add the study entry date for MIS
  D3_study_variables_for_MIS[,study_entry_date_MIS_d:=first_jan_2021]
  #add the cohort entry date for MIS
  D3_study_variables_for_MIS[,cohort_entry_date_MIS_d:=max(date_vax1, first_jan_2021), by="person_id"]
  
  D3_selection_criteria_d <- D3_study_variables_for_MIS[is.na(date_vax1) | study_exit_date <= cohort_entry_date_MIS_d, not_in_cohort_d:=1]
  rm(D3_study_variables_for_MIS)
  D3_selection_criteria_d <- D3_selection_criteria_d[, not_in_cohort_d := replace(.SD, is.na(.SD), 0), .SDcols = "not_in_cohort_d"]
  
  tempname<-paste0("D3_selection_criteria_d",suffix[[subpop]])
  assign(tempname,D3_selection_criteria_d)
  save(tempname, file = paste0(dirtemp, tempname,".RData"),list=tempname)
  
  D4_population_d <- CreateFlowChart(
    dataset = D3_selection_criteria_d,
    listcriteria = c("not_in_cohort_d"),
    flowchartname = "Flowchart_cohort_d")
  
  #add the study exit date for MIS
  #D3_study_variables_for_MIS[covid_date>date_vax1,study_exit_date_MIS_d:=min(covid_date,study_exit_date,na.rm = T)]
  D4_population_d[covid_date<=date_vax1,history_covid:=1][covid_date>date_vax1 | is.na(covid_date),history_covid:=0]
  D4_population_d[covid_date>date_vax1 & !is.na(covid_date),study_exit_date_MIS_d:=min(covid_date-1,study_end,study_exit_date),by="person_id"]
  D4_population_d[covid_date<date_vax1 | is.na(covid_date),study_exit_date_MIS_d:=min(study_end,study_exit_date),by="person_id"]
  
  # calculate correct fup_days
  D4_population_d[, fup_days := correct_difftime(study_exit_date_MIS_d, cohort_entry_date_MIS_d)]
  
  D4_population_d<-D4_population_d[,.(person_id,sex,age_at_1_jan_2021,ageband_at_1_jan_2021,study_entry_date_MIS_d,cohort_entry_date_MIS_d,study_exit_date_MIS_d,date_vax1,covid_date,history_covid,age_at_date_vax_1,type_vax_1,not_in_cohort_d, fup_days, CV_at_date_vax_1, COVCANCER_at_date_vax_1, COVCOPD_at_date_vax_1,
                                      COVHIV_at_date_vax_1, COVCKD_at_date_vax_1, COVDIAB_at_date_vax_1,
                                      COVOBES_at_date_vax_1, COVSICKLE_at_date_vax_1, immunosuppressants_at_date_vax_1,
                                      at_risk_at_date_vax_1)]
  
  D4_population_d<-D4_population_d[study_exit_date_MIS_d>study_entry_date_MIS_d,]
  
  tempname<-paste0("D4_population_d",suffix[[subpop]])
  assign(tempname,D4_population_d)
  save(tempname, file = paste0(diroutput, tempname,".RData"),list=tempname)
  
  rm(list=tempname)
  rm(tempname)
  rm(D3_selection_criteria_d)
  if (this_datasource_has_subpopulations==T) rm(D4_population_d)
  rm(events_ALL_OUTCOMES,study_population,outcomes_covid)
  
}


