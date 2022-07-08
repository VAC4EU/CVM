#-----------------------------------------------
# set D3_events_ALL_OUTCOMES which contains the first outcome per person

# input: D4_study_population, D3_events_OUTCOME_narrow, D3_events_OUTCOME_possible, for all outcomes OUTCOME; conceptset datasets for CONTROL_events, D3_eventsSecondary_SECCOMP
# output: D3_events_ALL_OUTCOMES, list_outcomes_observed.RData

print("CREATE EVENTS PER ALL OUTCOMES")

OUTCOMEnoCOVID <- OUTCOME_events[OUTCOME_events!="COVID"]

list_outcomes_observed_temp <- vector(mode="list")
list_outcomes_observed_only_diagnosis_temp <- vector(mode="list")
list_outcomes_observed_for_QC <- vector(mode="list")
D3_events_ALL_OUTCOMES <- vector(mode="list")

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
  
  running_list_outcomes_observed <- c()
  running_list_outcomes_observed_only_diagnosis <- c()

    study_population <- as.data.table(study_population)
  
  empty_events_ALL_OUTCOMES <- study_population[1,.(person_id)]
  empty_events_ALL_OUTCOMES <- empty_events_ALL_OUTCOMES[,name_event := "test"] 
  empty_events_ALL_OUTCOMES <- empty_events_ALL_OUTCOMES[,date_event := as.Date('20010101',date_format)] 
  empty_events_ALL_OUTCOMES <- empty_events_ALL_OUTCOMES[name_event!="test",] 
  events_ALL_OUTCOMES <- empty_events_ALL_OUTCOMES
  for (OUTCOME in OUTCOMEnoCOVID) {
    print(OUTCOME)
    namedatasetnarrow <- paste0('D3_events',"_",OUTCOME,"_narrow",suffix[[subpop]])
    namedatasetpossible <- paste0('D3_events',"_",OUTCOME,"_possible",suffix[[subpop]])
    load(paste0(dirtemp,namedatasetnarrow,'.RData'))
    load(paste0(dirtemp,namedatasetpossible,'.RData'))
    
    dataset <- vector(mode="list")
    dataset[['narrow']] <- as.data.table(get(namedatasetnarrow))
    dataset[['possible']] <- as.data.table(get(namedatasetpossible))

    dataset[['narrow']] <- dataset[['narrow']][,.(person_id, date, codvar, meaning_of_event, event_record_vocabulary)]
    dataset[['possible']] <- dataset[['possible']][,.(person_id, date, codvar, meaning_of_event, event_record_vocabulary)]

    for (type in c("narrow","possible")) {
      if ( nrow(dataset[[type]]) > 0 ) {
        dataset[[type]] <-  dataset[[type]][!is.na(person_id),]
        dataset[[type]] <-  dataset[[type]][!is.na(date),]
      }
    }  
    if ( nrow(dataset[['narrow']]) == 0){ 
      dataset[['broad']] <- dataset[['possible']]
    }
    if ( nrow(dataset[['possible']]) == 0){ 
      dataset[['broad']] <- dataset[['narrow']]
    }
    if ( nrow(dataset[['narrow']]) > 0 & nrow(dataset[['possible']]) > 0 ){
    dataset[['broad']] <- as.data.table(rbind(dataset[['narrow']], dataset[['possible']]), fill = T)
    }
    
    for (type in c("narrow","broad")) {
      dataset[[type]][,name_event:=paste0(OUTCOME,'_',type)]
      if (nrow(dataset[[type]]) > 0 ){
        setnames(dataset[[type]], c("date", "codvar", "meaning_of_event", "event_record_vocabulary"),
                 c("date_event", "code_first_event", "meaning_of_first_event", "coding_system_of_code_first_event"), 
                 skip_absent = T)
        dataset[[type]] <- dataset[[type]][,.(person_id, date_event, name_event, code_first_event, meaning_of_first_event, coding_system_of_code_first_event)]
        dataset[[type]] <- unique(dataset[[type]])
        events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES, dataset[[type]], fill = T))
        running_list_outcomes_observed <- c(running_list_outcomes_observed, paste0(OUTCOME,'_',type))
      }
    }
    rm(namedatasetnarrow,list = namedatasetnarrow)
    rm(namedatasetpossible,list = namedatasetpossible)
  }
  running_list_outcomes_observed_only_diagnosis <- running_list_outcomes_observed
  running_list_outcomes_for_QC <- running_list_outcomes_observed_only_diagnosis
  # add CONTROL_events
  for (CONTROL in CONTROL_events){
    print(CONTROL)
    load(paste0(dirtemp,CONTROL,".RData"))
    filecovariate <- get(CONTROL)
    if (this_datasource_has_subpopulations == TRUE){
      selection = "!is.na(person_id)"
      for (meaningevent in exclude_meaning_of_event[[thisdatasource]][[subpop]]){
        selection <- paste0(selection," & meaning_of_event!= '",meaningevent,"'")
      }
      filecovariate = filecovariate[eval(parse(text = selection)),]
    }
    temp <- MergeFilterAndCollapse(
      listdatasetL = list(filecovariate),
      condition= "date >= study_entry_date - 365",
      key = "person_id",
      datasetS = study_population,
      saveintermediatedataset=F,
      strata=c("person_id"),
      summarystat = list(
        list(c("min"),"date","date_event")
        )
    )
    if (CONTROL=='CONTRHYPERT'){
      # for CONTRHYPERT the algorithm is 'either DX or at least 2 DP'
      DPconcept <- paste0('DP_',CONTROL)
      load(paste0(dirtemp,DPconcept,".RData"))
      tempDP <- MergeFilterAndCollapse(
        listdatasetL = list(get(DPconcept)),
        condition= "date >= study_entry_date - 365",
        key = "person_id",
        datasetS = study_population,
        saveintermediatedataset=F,
        strata=c("person_id"),
        summarystat = list(
          list(c("second"),"date","date_event_DP")
        )
      )
      tempDP <- tempDP[!is.na(date_event_DP),]
      temp <- merge(temp,tempDP, all = T, by="person_id")
      temp <- temp[is.na(date_event),date_event := date_event_DP]
      temp <- temp[!is.na(date_event_DP),date_event := min(date_event,date_event_DP), by = "person_id"]
      temp <- temp[,-"date_event_DP"]
      rm(list = DPconcept)
      rm(tempDP)
    }
    temp <- temp[,name_event:=CONTROL]
    events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES,temp,fill = T))
    running_list_outcomes_observed <- c(running_list_outcomes_observed,CONTROL)
    rm(list = CONTROL)
    rm(filecovariate)
  }
  #add composite events
  events_MIS_KD_narrow<-events_ALL_OUTCOMES[name_event=="MIS_narrow" | name_event=="KD_narrow"]
  events_MIS_KD_narrow<-events_MIS_KD_narrow[,name_event:="MIS_KD_narrow"]
  events_MIS_KD_broad<-events_ALL_OUTCOMES[name_event=="MIS_broad" | name_event=="KD_broad",]
  events_MIS_KD_broad<-events_MIS_KD_broad[,name_event:="MIS_KD_broad"]
  events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES,events_MIS_KD_narrow,fill = T))
  events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES,events_MIS_KD_broad,fill = T))
  running_list_outcomes_observed <- c(running_list_outcomes_observed,"MIS_KD_narrow","MIS_KD_broad")
  rm(events_MIS_KD_narrow)
  rm(events_MIS_KD_broad)
  
  # add DEATH
  load(paste0(dirtemp,"D3_events_DEATH.RData"))
  temp <- as.data.table(D3_events_DEATH)
  temp <- temp[,name_event:="DEATH"]
  temp <- temp[,date_event:=date]
  temp <- temp[,-"date"]
  events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES,temp,fill = T))
  running_list_outcomes_observed <- c(running_list_outcomes_observed,"DEATH")
  rm(temp)
  
  # add secondary components
  
  for (SECCOMP in SECCOMPONENTS) { 
    print(SECCOMP)
    nameobjectSECCOMP <- paste0('D3_eventsSecondary',"_",SECCOMP,suffix[[subpop]])
    load(paste0(dirtemp,paste0(nameobjectSECCOMP,".RData")))
    temp <- get(nameobjectSECCOMP)
    temp <- temp[,name_event:= SECCOMP]
    if ( nrow(temp) > 0 ){
      setkey(temp,person_id, date)
      temp <- temp[,`:=`(date_event = min(date),
                                                   code_first_event = codvarA[1], meaning_of_first_event = meaning_of_eventA[1], coding_system_of_code_first_event = event_record_vocabularyA[1]),by = c("person_id", "name_event")]
      temp <- temp[,.(person_id, date_event, name_event, code_first_event, meaning_of_first_event, coding_system_of_code_first_event)]
      temp <- unique(temp)
      events_ALL_OUTCOMES <- as.data.table(rbind(events_ALL_OUTCOMES, temp, fill = T))
      running_list_outcomes_observed <- c(running_list_outcomes_observed, SECCOMP)
      running_list_outcomes_for_QC <- c(running_list_outcomes_for_QC,SECCOMP)
    }
    rm(nameobjectSECCOMP, list = nameobjectSECCOMP)
    }
    
    list_outcomes_observed_temp <- running_list_outcomes_observed
    list_outcomes_observed_only_diagnosis_temp <- running_list_outcomes_observed_only_diagnosis
    list_outcomes_observed_for_QC <- running_list_outcomes_for_QC
    events_ALL_OUTCOMES <- events_ALL_OUTCOMES
    
    nametemp<-paste0("list_outcomes_observed",suffix[[subpop]])
    assign(nametemp,list_outcomes_observed_temp)
    save(nametemp,file=paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"),list=nametemp)
    
    nametemp<-paste0("list_outcomes_observed_only_diagnosis",suffix[[subpop]])
    assign(nametemp,list_outcomes_observed_only_diagnosis_temp)
    save(nametemp,file=paste0(dirtemp,"list_outcomes_observed_only_diagnosis",suffix[[subpop]],".RData"),list=nametemp)
    
    nametemp<-paste0("list_outcomes_observed_for_QC",suffix[[subpop]])
    assign(nametemp,list_outcomes_observed_for_QC)
    save(nametemp,file=paste0(dirtemp,"list_outcomes_observed_for_QC",suffix[[subpop]],".RData"),list=nametemp)

    nametemp<-paste0("D3_events_ALL_OUTCOMES",suffix[[subpop]])
    assign(nametemp,events_ALL_OUTCOMES)
    save(nametemp,file=paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"),list=nametemp)
    rm(list=paste0("D3_events_ALL_OUTCOMES",suffix[[subpop]]))
    
    rm(list=paste0("D4_study_population", suffix[[subpop]]))
}


rm(events_ALL_OUTCOMES, dataset,D3_events_DEATH, empty_events_ALL_OUTCOMES,study_population,temp,list_outcomes_observed_temp,list_outcomes_observed_only_diagnosis_temp)
