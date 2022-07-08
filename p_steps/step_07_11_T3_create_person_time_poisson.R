# COUNT PERSON TIME PER COVID PER POISSON
#-----------------------------------------------
# To estimate the weekly incidence rates of risks in 2020 by data source for poisson

# input: D3_events_ALL_OUTCOMES, D3_vaxweeks_poisson
# output: D4_persontime_risk_month_poisson

print("COUNT PERSON TIME by month for risks (Poisson)")

persontime_risk_year <- vector(mode = 'list')

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = "20200101"
  
  load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_vaxweeks_poisson",suffix[[subpop]],".RData"))
  
  list_outcomes<-get(paste0("list_outcomes_observed", suffix[[subpop]]))
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  study_population<-get(paste0("D3_vaxweeks_poisson", suffix[[subpop]]))
  
  endyear<- substr(study_population[,max(end_date_of_period)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  list_recurrent_outcomes <- list_outcomes[str_detect(list_outcomes, "^GENCONV_") | str_detect(list_outcomes, "^ANAPHYL_")]
  list_outcomes <- setdiff(list_outcomes, list_recurrent_outcomes)
  
  sex_vect <- c(0,1)
  
  for (ageband in Agebands_labels) {
    for (sex in  sex_vect) {
      nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), "_", sex, suffix[[subpop]])
      assign(nameoutput, study_population[ageband_at_study_entry == ageband & Gender == sex, ])
      save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
      rm(list=nameoutput)
    }
  }
  
  df_events_ages <- c(paste0("pop_age_", gsub("-", "_", Agebands_labels), "_", "0"),
                      paste0("pop_age_", gsub("-", "_", Agebands_labels), "_", "1"))
  
  for (events_df_sex in df_events_ages) {
    print(paste("age/sex", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    print("recurrent")
    
    nameoutput <- paste0("Recurrent_output_file",suffix[[subpop]])
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex, suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "start_date_of_period",
      End_date = "end_date_of_period",
      Strata = c("DAP", "Gender", "ageband_at_study_entry",
                 "COVID19", "Vaccine1", "Vaccine2", "Dose1", "Dose2", "CV_at_study_entry",
                 "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                 "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                 "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                 "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                 "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                 "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                 "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                 "all_risk_factors_at_date_vax"),
      Name_event = "name_event",
      Date_event = "date_event",
      Increment = "month",
      Outcomes_rec =   list_recurrent_outcomes,
      Aggregate = T,
      Rec_period = c(rep(30, length(list_recurrent_outcomes)))
    ))
    
    save(nameoutput, file=paste0(dirtemp,"D3_recurrent_year",suffix[[subpop]],".RData"),list=nameoutput)
    rm(list=nameoutput)
    
    print("normal")
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "start_date_of_period",
      End_date = "end_date_of_period",
      Strata = c("DAP", "Gender", "ageband_at_study_entry",
                 "COVID19", "Vaccine1", "Vaccine2", "Dose1", "Dose2", "CV_at_study_entry",
                 "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                 "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                 "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                 "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                 "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                 "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                 "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                 "all_risk_factors_at_date_vax"),
      Name_event = "name_event",
      Date_event = "date_event",
      Increment = "month",
      Outcomes_nrec = list_outcomes,
      Aggregate = T
    ))
    
    load(paste0(dirtemp,"D3_recurrent_year",suffix[[subpop]],".RData"))
    print("Merging")
    nameoutput<-paste0("Output_file",suffix[[subpop]])
    assign(nameoutput,merge(get(paste0("Output_file",suffix[[subpop]])), get(paste0("Recurrent_output_file",suffix[[subpop]])) ,
                            by = c("DAP", "Gender", "ageband_at_study_entry",
                                   "COVID19", "Vaccine1", "Vaccine2", "Dose1", "Dose2", "CV_at_study_entry",
                                   "COVCANCER_at_study_entry", "COVCOPD_at_study_entry", "COVHIV_at_study_entry",
                                   "COVCKD_at_study_entry", "COVDIAB_at_study_entry", "COVOBES_at_study_entry",
                                   "COVSICKLE_at_study_entry", "IMMUNOSUPPR_at_study_entry",
                                   "any_risk_factors_at_study_entry", "all_risk_factors_at_study_entry",
                                   "CV_at_date_vax", "COVCANCER_at_date_vax", "COVCOPD_at_date_vax", "COVHIV_at_date_vax",
                                   "COVCKD_at_date_vax", "COVDIAB_at_date_vax", "COVOBES_at_date_vax",
                                   "COVSICKLE_at_date_vax", "IMMUNOSUPPR_at_date_vax", "any_risk_factors_at_date_vax",
                                   "all_risk_factors_at_date_vax", "month", "Persontime"),
                            all = T)
    )
    print("Saving")
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign(paste0("Output_file",suffix[[subpop]]),rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  for (i in names(get(paste0("Output_file",suffix[[subpop]])))){
    get(paste0("Output_file",suffix[[subpop]]))[is.na(get(i)), (i):=0]
  }
  
  persontime_risk_year <- get(paste0("Output_file",suffix[[subpop]]))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(persontime_risk_year,file=paste0(thisdirexp,"D4_persontime_risk_month_poisson.csv"))
  
  nameoutput<-paste0("D4_persontime_risk_month_poisson", suffix[[subpop]])
  assign(nameoutput,persontime_risk_year)
  save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)
  
  rm(list=paste0("Output_file",suffix[[subpop]]) )
  rm(list=paste0("Recurrent_output_file",suffix[[subpop]]))
  rm(list=nameoutput)
  rm(list=paste0("D3_vaxweeks_poisson", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("list_outcomes_observed", suffix[[subpop]]))
}

for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_risk_month_poisson")
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  assign(tempname,fread(paste0(thisdirexp,tempname,".csv")))
  thisdirsmallcountsremoved <- ifelse(this_datasource_has_subpopulations == FALSE,dirsmallcountsremoved,dirsmallcountsremovedsubpop[[subpop]])
  col<-colnames(get(tempname))[-(1:6)]
  temp<-paste0(col,"=5")
  temp2<-paste("c(",paste(temp, collapse = ','),")")
  suppressWarnings(
    DRE_Treshold(
      Inputfolder = thisdirexp,
      Outputfolder = thisdirsmallcountsremoved,
      Delimiter = ",",
      Varlist = c(eval(parse(text=(temp2)))),
      FileContains = "D4_persontime_risk_month_poisson"
    )
  )
  rm(list=tempname)
}
# rm(list = nameobject)
rm(persontime_risk_year,events_ALL_OUTCOMES,study_population)
