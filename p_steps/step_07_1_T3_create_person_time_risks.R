# COUNT PERSON TIME PER COVID
#-----------------------------------------------
#To estimate the weekly incidence rates of risks in 2020 by data source

# input: D3_vaxweeks_including_not_vaccinated, D3_events_ALL_OUTCOMES.RData ,list_outcomes_observed.RData
# output: D4_persontime_risk_week (exported to csv)

print("COUNT PERSON TIME by week for risks")


D4_persontime_risk_week <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
    
  start_week=seq.Date(as.Date("20200106","%Y%m%d"),Sys.Date(),by = "week")
    
    load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))
    load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
    load(paste0(dirtemp,"D3_vaxweeks_including_not_vaccinated",suffix[[subpop]],".RData"))
    
    list_outcomes<-get(paste0("list_outcomes_observed", suffix[[subpop]]))
    events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
    study_population<-get(paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
  
  max_exit<-study_population[,max(end_date_of_period)]
  last_event<-events_ALL_OUTCOMES[,max(date_event)]
  
  if (last_event<ymd("20200101")) {
    next
  }
  end_persontime_studytime<-min(max_exit,last_event)
  start_week=start_week[start_week<=end_persontime_studytime]
  end_week=start_week+6
  start_week=gsub("-","", start_week)
  end_week=gsub("-","", end_week)
  
  list_recurrent_outcomes <- list_outcomes[str_detect(list_outcomes, "^GENCONV_") | str_detect(list_outcomes, "^ANAPHYL_")]
  list_outcomes <- setdiff(list_outcomes, list_recurrent_outcomes)
  
  for (ageband in Agebands_labels) {
    for (current_sex in c(0, 1)) {
      nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), "_", current_sex, suffix[[subpop]])
      assign(nameoutput, study_population[ageband_at_study_entry == ageband & sex == current_sex, ])
      save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
      rm(list=nameoutput)
    }
  }
  
  df_events_ages <- c(paste0("pop_age_", gsub("-", "_", Agebands_labels), "_0"),
                      paste0("pop_age_", gsub("-", "_", Agebands_labels), "_1"))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    print("recurrent")
  
    for (i in 1:length(start_week)){
      start_persontime_studytime = start_week[i]
      end_persontime_studytime = end_week[i]
      nameoutput <- paste0("Output_file",start_week[i])
      print(nameoutput)
      assign(nameoutput,CountPersonTime(
        Dataset_events = events_ALL_OUTCOMES,
        Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
        Person_id = "person_id",
        Start_study_time = start_persontime_studytime,
        End_study_time = end_persontime_studytime,
        Start_date = "start_date_of_period",
        End_date = "end_date_of_period",
        #Birth_date = "date_of_birth",
        Strata = c("sex","ageband_at_study_entry","Dose","type_vax", "CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                   "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
        Name_event = "name_event",
        Date_event = "date_event",
        #Age_bands = c(0,19,29,39,49,59,69,79),F
        Increment="week",
        Outcomes_rec = list_recurrent_outcomes, 
        Unit_of_age = "year",
        include_remaning_ages = T,
        Aggregate = T,
        Rec_period = c(rep(30, length(list_recurrent_outcomes)))
      )
      )
      if (i==1) {
        persontime_risk_week<-get(paste0("Output_file",start_week[i]))
        rm(list=paste0("Output_file",start_week[i]))
      }else{
        persontime_risk_week<-rbind(persontime_risk_week,get(paste0("Output_file",start_week[i])),fill=TRUE)
        rm(list=paste0("Output_file",start_week[i]))
      }
    }
    recurrent_persontime_risk_week <- persontime_risk_week
    rm(persontime_risk_week)
    for (i in 1:length(start_week)){
      start_persontime_studytime = start_week[i]
      end_persontime_studytime = end_week[i]
      nameoutput <- paste0("Output_file",start_week[i])
      print(nameoutput)
      assign(nameoutput,CountPersonTime(
        Dataset_events = events_ALL_OUTCOMES,
        Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
        Person_id = "person_id",
        Start_study_time = start_persontime_studytime,
        End_study_time = end_persontime_studytime,
        Start_date = "start_date_of_period",
        End_date = "end_date_of_period",
        #Birth_date = "date_of_birth",
        Strata = c("sex","ageband_at_study_entry","Dose","type_vax", "CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
                   "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
        Name_event = "name_event",
        Date_event = "date_event",
        #Age_bands = c(0,19,29,39,49,59,69,79),F
        Increment="week",
        Outcomes_nrec = list_outcomes, 
        Unit_of_age = "year",
        include_remaning_ages = T,
        Aggregate = T
      )
      )
      if (i==1) {
        persontime_risk_week<-get(paste0("Output_file",start_week[i]))
        rm(list=paste0("Output_file",start_week[i]))
      }else{
        persontime_risk_week<-rbind(persontime_risk_week,get(paste0("Output_file",start_week[i])),fill=TRUE)
        rm(list=paste0("Output_file",start_week[i]))
      }
    }
    
    print("Merging")
    nameoutput<-paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, merge(persontime_risk_week, recurrent_persontime_risk_week,
                                  by = c("sex","ageband_at_study_entry","Dose","type_vax", "CV", "COVCANCER", "COVCOPD",
                                         "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                                         "any_risk_factors", "week", "Persontime"), all = T))
    print("Saving")
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"), list=nameoutput)
    rm(list=paste0(events_df_sex, suffix[[subpop]]))
    
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
  
  persontime_risk_week <- get(paste0("Output_file",suffix[[subpop]]))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(persontime_risk_week,file=paste0(thisdirexp,"D4_persontime_risk_week.csv"))
  
  nameoutput<-paste0("D4_persontime_risk_week",suffix[[subpop]])
  assign(nameoutput,persontime_risk_week)
  save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)
  
  rm(list=paste0("D4_persontime_risk_week",suffix[[subpop]]))
  rm(list=paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
  rm(list=paste0("list_outcomes_observed", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
}


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_risk_week")
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  assign(tempname,fread(paste0(thisdirexp,tempname,".csv")))
    thisdirsmallcountsremoved <- ifelse(this_datasource_has_subpopulations == FALSE,dirsmallcountsremoved,dirsmallcountsremovedsubpop[[subpop]])
  col<-colnames(get(tempname))[-(1:3)]
  temp<-paste0(col,"=5")
  temp2<-paste("c(",paste(temp, collapse = ','),")")
  suppressWarnings(
    DRE_Treshold(
      Inputfolder = thisdirexp,
      Outputfolder = thisdirsmallcountsremoved,
      Delimiter = ",",
      Varlist = c(eval(parse(text=(temp2)))),
      FileContains = "D4_persontime_risk_week"
    )
  )
  rm(list=tempname)
}

rm(persontime_risk_week,events_ALL_OUTCOMES,study_population,list_outcomes,recurrent_persontime_risk_week)
