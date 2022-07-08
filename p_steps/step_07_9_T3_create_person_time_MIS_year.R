# COUNT PERSON TIME PER COVID PER MIS/MYOCARD
#-----------------------------------------------
# To estimate the weekly incidence rates of risks in 2020 by data source for MIS/Myocard

# input: D3_events_ALL_OUTCOMES, D4_population_b, D4_population_c, D4_population_d
# output: D4_persontime_b, D4_persontime_monthly_b, D4_persontime_c, D4_persontime_monthly_c, D4_persontime_d, D4_persontime_monthly_d

#COHORT B

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = "20200101"
  
  load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_b",suffix[[subpop]],".RData"))
  
  events_ALL_OUTCOMES<-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  population_b<-get(paste0("D4_population_b", suffix[[subpop]]))
  rm(list=paste0("D4_population_b", suffix[[subpop]]))
  
  endyear<- substr(population_b[,max(study_exit_date_MIS_b)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_b[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_b",
      End_date = "study_exit_date_MIS_b",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="year",
      Outcomes_nrec = list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  nameoutput<-paste0("D4_persontime_b",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_b")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  
  rm(list=nameoutput)
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_b[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput,CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_b",
      End_date = "study_exit_date_MIS_b",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="month",
      Outcomes_nrec =  list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  nameoutput<-paste0("D4_persontime_monthly_b",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_monthly_b")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  
  rm(list=nameoutput)
  
  
}

rm(population_b)


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_b")
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
      FileContains = "D4_persontime_b"
    )
  )
  rm(list=tempname)
}

for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_monthly_b")
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
      FileContains = "D4_persontime_monthly_b"
    )
  )
  rm(list=tempname)
}

#COHORT C
for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = "20200101"
  
  load(paste0(diroutput,"D4_population_c",suffix[[subpop]],".RData"))
  
  population_c<-get(paste0("D4_population_c", suffix[[subpop]]))
  rm(list=paste0("D4_population_c", suffix[[subpop]]))
  
  endyear<- substr(population_c[,max(study_exit_date_MIS_c)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_c[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_c",
      End_date = "study_exit_date_MIS_c",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="year",
      Outcomes_nrec =  list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  nameoutput<-paste0("D4_persontime_c",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_c")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  rm(list=nameoutput)
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_c[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_c",
      End_date = "study_exit_date_MIS_c",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021"),
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="month",
      Outcomes_nrec =   list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  
  nameoutput<-paste0("D4_persontime_monthly_c",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_monthly_c")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  
  rm(list=nameoutput)
  rm(nameoutput)
  rm(population_c)
}


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_c")
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
      FileContains = "D4_persontime_c"
    )
  )
  rm(list=tempname)
}

for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_monthly_c")
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
      FileContains = "D4_persontime_monthly_c"
    )
  )
  rm(list=tempname)
}

#COHORT D
for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = "20200101"
  
  load(paste0(diroutput,"D4_population_d",suffix[[subpop]],".RData"))
  
  population_d<-get(paste0("D4_population_d", suffix[[subpop]]))
  rm(list=paste0("D4_population_d", suffix[[subpop]]))
  
  endyear<- substr(population_d[,max(study_exit_date_MIS_d)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_d[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_d",
      End_date = "study_exit_date_MIS_d",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021","type_vax_1","history_covid"), #add covid before vaccine
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="year",
      Outcomes_nrec = list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  nameoutput<-paste0("D4_persontime_d",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_d")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  rm(list=nameoutput)
  
  for (ageband in Agebands_labels) {
    nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
    assign(nameoutput, population_d[ageband_at_1_jan_2021 == ageband, ])
    save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
    rm(list=nameoutput)
  }
  
  df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))
  
  for (events_df_sex in df_events_ages) {
    print(paste("Age", substring(events_df_sex, 9)))
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    
    nameoutput <- paste0("Output_file",suffix[[subpop]])
    assign(nameoutput, CountPersonTime(
      Dataset_events = events_ALL_OUTCOMES,
      Dataset = get(paste0(events_df_sex,suffix[[subpop]])),
      Person_id = "person_id",
      Start_study_time = start_persontime_studytime,
      End_study_time = end_persontime_studytime,
      Start_date = "cohort_entry_date_MIS_d",
      End_date = "study_exit_date_MIS_d",
      #Birth_date = "date_of_birth",
      Strata = c("sex","ageband_at_1_jan_2021","type_vax_1","history_covid"), #add covid before vaccine
      Name_event = "name_event",
      Date_event = "date_event",
      #Age_bands = c(0,19,29,39,49,59,69,79),
      Increment="month",
      Outcomes_nrec = list_outcomes_MIS, 
      #Unit_of_age = "year",
      #include_remaning_ages = T,
      Aggregate = T
    ))
    
    save(nameoutput, file = paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"),list=nameoutput)
    rm(list=paste0(events_df_sex,suffix[[subpop]]))
    
  }
  
  vect_df_persontime <- list()
  for (events_df_sex in df_events_ages) {
    load(paste0(dirtemp, events_df_sex,suffix[[subpop]], ".RData"))
    vect_df_persontime <- append(vect_df_persontime, list(get(paste0("Output_file",suffix[[subpop]]) )))
  }
  
  assign("Output_file", rbindlist(vect_df_persontime))
  rm(vect_df_persontime)
  
  nameoutput<-paste0("D4_persontime_monthly_d",suffix[[subpop]])
  assign(nameoutput ,Output_file)
  rm(Output_file)
  save(nameoutput, file = paste0(diroutput, nameoutput,".RData"),list=nameoutput)
  
  nameoutput2<-paste0("D4_persontime_monthly_d")
  assign(nameoutput2 ,get(nameoutput))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(get(nameoutput2),file=paste0(thisdirexp,nameoutput2,".csv"))
  rm(list=nameoutput)
  rm(population_d)
  
}


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_d")
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
      FileContains = "D4_persontime_d"
    )
  )
  rm(list=tempname)
}

for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_monthly_d")
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
      FileContains = "D4_persontime_monthly_d"
    )
  )
  rm(list=tempname)
}

rm(events_ALL_OUTCOMES)