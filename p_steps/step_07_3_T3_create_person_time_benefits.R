# COUNT PERSON TIME PER COVID
#-----------------------------------------------
#To estimate the weekly incidence rates of COVID-19 (overall and by severity level) in 2020 by data source

# input: D3_vaxweeks_including_not_vaccinated, D3_outcomes_covid ,list_outcomes_observed_COVID
# output: D4_persontime_benefit_week (exported to csv)

print("COUNT PERSON TIME PER COVID by week benefits")

persontime_benefit_week <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  start_week=seq.Date(as.Date("20200106","%Y%m%d"),Sys.Date(),by = "week")
  
  load(paste0(dirpargen,"list_outcomes_observed_COVID",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_outcomes_covid",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_vaxweeks_including_not_vaccinated",suffix[[subpop]],".RData"))
  
  list_outcomes<-get(paste0("list_outcomes_observed_COVID", suffix[[subpop]]))
  events_ALL_OUTCOMES<-get(paste0("D3_outcomes_covid", suffix[[subpop]]))
  study_population<-get(paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
  
  max_exit<-study_population[,max(end_date_of_period)]
  if (nrow(events_ALL_OUTCOMES) != 0) {
    last_event<-events_ALL_OUTCOMES[,max(date_event)]
    if (last_event<ymd("20200101")) {
      next
    }
    end_persontime_studytime<-min(max_exit,last_event)
  } else {
    end_persontime_studytime<-max_exit
  }
  start_week=start_week[start_week<=end_persontime_studytime]
  end_week=start_week+6
  start_week=gsub("-","", start_week)
  end_week=gsub("-","", end_week)
  
  for (i in 1:length(start_week)){
    start_persontime_studytime = start_week[i]
    end_persontime_studytime = end_week[i]
    nameoutput <- paste0("Output_file",start_week[i])
    print(nameoutput)
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = study_population,
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
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="week",
      Outcomes_nrec = list_outcomes, 
      # Unit_of_age = "year",
      # include_remaning_ages = T,
      Aggregate = T
    )
    )
    if (i==1) {
      persontime_benefit_week<-get(paste0("Output_file",start_week[i]))
      rm(list=paste0("Output_file",start_week[i]))
    }else{
      persontime_benefit_week<-rbind(persontime_benefit_week,get(paste0("Output_file",start_week[i])),fill=TRUE)
      rm(list=paste0("Output_file",start_week[i]))
    }
  }
  
  
  # fwrite(persontime_benefit_week,file=paste0(thisdirexp,"D4_persontime_benefit_week.csv"))
  # if (this_datasource_has_subpopulations == TRUE){ 
  #   D4_persontime_benefit_week[[subpop]] <- persontime_benefit_week
  # }else{
  #   D4_persontime_benefit_week <- persontime_benefit_week
  # }
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(persontime_benefit_week,file=paste0(thisdirexp,"D4_persontime_benefit_week.csv"))
  
  nameoutput<-paste0("D4_persontime_benefit_week",suffix[[subpop]])
  assign(nameoutput,persontime_benefit_week)
  save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)
  
  rm(list=paste0("D4_persontime_benefit_week",suffix[[subpop]]))
  rm(list=paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
  rm(list=paste0("list_outcomes_observed_COVID", suffix[[subpop]]))
  rm(list=paste0("D3_outcomes_covid", suffix[[subpop]]))
}


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_benefit_week")
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
      FileContains = "D4_persontime_benefit_week"
    )
  )
  rm(list=tempname)
}

rm(persontime_benefit_week,events_ALL_OUTCOMES,study_population,list_outcomes)
