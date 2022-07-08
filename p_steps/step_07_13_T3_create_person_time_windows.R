# COUNT PERSON TIME PER MYOCARD COHORT
#-----------------------------------------------
# To estimate the incidence rates of risks in periods

# input: D3_myocard_cohort, D3_events_ALL_OUTCOMES.RData ,list_outcomes_observed.RData
# output: D4_persontime_windows (exported to csv)


print("COUNT PERSON TIME by year for risks for MYOCARD")

persontime_risk_year <- vector(mode = 'list')

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = format(study_start, "%Y%m%d")
  
  load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_myocard_cohort",suffix[[subpop]],".RData"))
  
  list_outcomes<-get(paste0("list_outcomes_observed", suffix[[subpop]]))
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  study_population<-get(paste0("D3_myocard_cohort", suffix[[subpop]]))
  
  rm(list=paste0("D3_myocard_cohort", suffix[[subpop]]))
  
  endyear<- substr(study_population[,max(end_date_of_period)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  list_recurrent_outcomes <- list_outcomes[str_detect(list_outcomes, "^GENCONV_") | str_detect(list_outcomes, "^ANAPHYL_")]
  list_outcomes <- setdiff(list_outcomes, list_recurrent_outcomes)
  
  print("recurrent")
  
  nameoutput <- paste0("Recurrent_output_file", suffix[[subpop]])
  assign(nameoutput,CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = start_persontime_studytime,
    End_study_time = end_persontime_studytime,
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex","Dose","type_vax","Period", "COVID19", "CV", "COVCANCER", "COVCOPD",
               "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
    Name_event = "name_event",
    Date_event = "date_event",
    Age_bands = c(0, 4, 11, 17, 24, 29, 39, 49, 59, 69, 79),
    Increment="year",
    Outcomes_rec =   list_recurrent_outcomes, 
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T,
    Rec_period = c(rep(30, length(list_recurrent_outcomes)))
  ))
  
  save(nameoutput, file=paste0(dirtemp,"D3_recurrent_year",suffix[[subpop]],".RData"),list=nameoutput)
  rm(list=nameoutput)
    
  print("normal")
  nameoutput <- paste0("Output_file",suffix[[subpop]])
  assign(nameoutput,CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = start_persontime_studytime,
    End_study_time = end_persontime_studytime,
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    Birth_date = "date_of_birth",
    Strata = c("sex","Dose","type_vax","Period", "COVID19", "CV", "COVCANCER", "COVCOPD",
               "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
    Name_event = "name_event",
    Date_event = "date_event",
    Age_bands = c(0, 4, 11, 17, 24, 29, 39, 49, 59, 69, 79),
    Increment="year",
    Outcomes_nrec = list_outcomes, 
    Unit_of_age = "year",
    include_remaning_ages = T,
    Aggregate = T
  ))
    
  load(paste0(dirtemp,"D3_recurrent_year",suffix[[subpop]],".RData"))
  print("Merging")
  
  assign(paste0("Output_file",suffix[[subpop]]), merge(get(paste0("Output_file",suffix[[subpop]])),
                                                       get(paste0("Recurrent_output_file",suffix[[subpop]])) ,
                          by = c("sex","Dose","Ageband","type_vax","Period", "COVID19", "CV", "COVCANCER", "COVCOPD",
                                 "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors", "year", "Persontime"),
                          all = T)
  )
  print("Saving")
  
  for (i in names(get(paste0("Output_file",suffix[[subpop]])))){
    get(paste0("Output_file",suffix[[subpop]]))[is.na(get(i)), (i):=0]
  }
  
  persontime_windows <- get(paste0("Output_file",suffix[[subpop]]))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(persontime_windows,file=paste0(thisdirexp,"D4_persontime_windows.csv"))
  
  nameoutput<-paste0("D4_persontime_windows",suffix[[subpop]])
  assign(nameoutput,persontime_windows)
  save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)
  
  rm(list=paste0("Output_file",suffix[[subpop]]) )
  rm(list=paste0("Recurrent_output_file",suffix[[subpop]]))
  rm(list=nameoutput)
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("list_outcomes_observed", suffix[[subpop]]))
}

for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_windows")
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
      FileContains = "D4_persontime_windows"
    )
  )
  rm(list=tempname)
}
# rm(list = nameobject)
rm(persontime_windows,events_ALL_OUTCOMES,study_population)
