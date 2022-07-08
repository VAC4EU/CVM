# COUNT PERSON TIME PER COVID PER VACCINATED
#-----------------------------------------------
# To estimate the weekly incidence rates of risks in 2020 by data source for vaccinated cohort

# input: D3_events_ALL_OUTCOMES, D3_vaxweeks_vaccin_cohort
# output: D4_persontime_risk_month

print("COUNT PERSON TIME by month for risks")

D4_persontime_risk_month <- vector(mode = 'list')

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  start_month=seq.Date(as.Date("20200101","%Y%m%d"),Sys.Date(),by = "month")

  load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_vaxweeks_vaccin_cohort",suffix[[subpop]],".RData"))
  
  list_outcomes<-get(paste0("list_outcomes_observed", suffix[[subpop]]))
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  study_population<-get(paste0("D3_vaxweeks_vaccin_cohort", suffix[[subpop]]))


  max_exit<-study_population[,max(study_exit_date)]
  last_event<-events_ALL_OUTCOMES[,max(date_event)]
  if (last_event<ymd("20200101")) {
    next
  }
  end_persontime_studytime<-min(max_exit,last_event)
  start_month=start_month[start_month<=end_persontime_studytime]
  end_month=start_month %m+% months(1) - 1
  start_month=gsub("-","", start_month)
  end_month=gsub("-","", end_month)
  
  list_recurrent_outcomes <- list_outcomes[str_detect(list_outcomes, "^GENCONV_") | str_detect(list_outcomes, "^ANAPHYL_")]
  list_outcomes <- setdiff(list_outcomes, list_recurrent_outcomes)
  
  for (i in 1:length(start_month)){
    start_persontime_studytime = start_month[i]
    end_persontime_studytime = end_month[i]
    nameoutput <- paste0("Output_file",start_month[i])
    print(nameoutput)
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = study_population,
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "study_entry_date",
      End_date = "study_exit_date",
      #Birth_date = "date_of_birth",
      Strata = c("sex", "ageband_at_study_entry", "at_risk_at_study_entry"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),F
      Increment="month",
      Outcomes_rec = list_recurrent_outcomes, 
      Unit_of_age = "year",
      include_remaning_ages = T,
      Aggregate = T,
      Rec_period = c(rep(30, length(list_recurrent_outcomes)))
    )
    )
    if (i==1) {
      persontime_risk_month<-get(paste0("Output_file",start_month[i]))
      rm(list=paste0("Output_file",start_month[i]))
    }else{
      persontime_risk_month<-rbind(persontime_risk_month,get(paste0("Output_file",start_month[i])),fill=TRUE)
      rm(list=paste0("Output_file",start_month[i]))
    }
  }
  recurrent_persontime_risk_month <- persontime_risk_month
  rm(persontime_risk_month)
  for (i in 1:length(start_month)){
    start_persontime_studytime = start_month[i]
    end_persontime_studytime = end_month[i]
    nameoutput <- paste0("Output_file",start_month[i])
    print(nameoutput)
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = study_population,
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "study_entry_date",
      End_date = "study_exit_date",
      #Birth_date = "date_of_birth",
      Strata = c("sex", "ageband_at_study_entry", "at_risk_at_study_entry"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),F
      Increment="month",
      Outcomes_nrec = list_outcomes, 
      Unit_of_age = "year",
      include_remaning_ages = T,
      Aggregate = T
    )
    )
    if (i==1) {
      persontime_risk_month<-get(paste0("Output_file",start_month[i]))
      rm(list=paste0("Output_file",start_month[i]))
    }else{
      persontime_risk_month<-rbind(persontime_risk_month,get(paste0("Output_file",start_month[i])),fill=TRUE)
      rm(list=paste0("Output_file",start_month[i]))
    }
  }
  
  persontime_risk_month <- merge(persontime_risk_month, recurrent_persontime_risk_month,
                                by = c("sex","ageband_at_study_entry", "month", "at_risk_at_study_entry", "Persontime"), all = T)
  
  for (i in names(persontime_risk_month)){
    persontime_risk_month[is.na(get(i)), (i):=0]
  }
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(persontime_risk_month,file=paste0(thisdirexp,"D4_persontime_risk_month.csv"))
  
  nameoutput<-paste0("D4_persontime_risk_month",suffix[[subpop]])
  assign(nameoutput,persontime_risk_month)
  save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)
  

  rm(list=nameoutput)
  rm(list=paste0("D3_vaxweeks_vaccin_cohort", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("list_outcomes_observed", suffix[[subpop]]))
  
  # thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  # fwrite(persontime_risk_month,file=paste0(thisdirexp,"D4_persontime_risk_month.csv"))
  # if (this_datasource_has_subpopulations == TRUE){ 
  #   D4_persontime_risk_month[[subpop]] <- persontime_risk_month
  # }else{
  #   D4_persontime_risk_month <- persontime_risk_month
  # }
}



for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_risk_month")
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
      FileContains = "D4_persontime_risk_month"
    )
  )
  rm(list=tempname)
}

rm(persontime_risk_month,recurrent_persontime_risk_month,events_ALL_OUTCOMES,study_population)
#rm(D3_vaxweeks_vaccin_cohort,D4_persontime_risk_month,persontime_risk_month,study_population,events_ALL_OUTCOMES)
