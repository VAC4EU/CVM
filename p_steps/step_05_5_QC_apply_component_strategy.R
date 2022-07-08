#------------------------------------------------------------------
# split outcomes in components

# input: D4_study_population, CONCEPT, named OUTCOME
# output: QC_all_components_OUTCOME.RData


print('BREAK DOWN OUTCOMES PER MEANINGS')


COHORT_to_be_used <- vector(mode = 'list')
COHORT_dates <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
  
  COHORT_TMP <- study_population[,.(person_id,study_entry_date,study_exit_date)]
  
  FirstJan<-vector(mode="list")
  for (year in c("2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018","2019")) {
    FirstJan[[year]]<-as.Date(as.character(paste0(year,"0101")), date_format)
  }
  
  firstyear = firstYearComponentAnalysis
  secondyear = secondYearComponentAnalysis
  
  COHORT_TMP <- COHORT_TMP[study_entry_date<=as.Date(FirstJan[[secondyear]]) + 365 & study_exit_date >= as.Date(FirstJan[[firstyear]]),]

  COHORT_dates <- COHORT_TMP
  COHORT_to_be_used <- COHORT_TMP[,.(person_id)]
  rm(list=paste0("D4_study_population", suffix[[subpop]]))
}

 
for (OUTCOME in OUTCOME_events) {
  print(OUTCOME)
  OUTCOME_components <- vector(mode = 'list')
  for (subpop in subpopulations_non_empty) {
    namedatasetnarrow <- paste0(OUTCOME,"_narrow")
    namedatasetpossible <- paste0(OUTCOME,"_possible")
    load(paste0(dirtemp,namedatasetnarrow,'.RData'))
    load(paste0(dirtemp,namedatasetpossible,'.RData'))

    COHORT_TMP <- COHORT_to_be_used
    COHORT_DATES <- COHORT_dates
    
    dataset <- vector(mode="list")
    dataset[['narrow']] <- as.data.table(get(namedatasetnarrow))
    dataset[['possible']] <- as.data.table(get(namedatasetpossible))
    dataset[['narrow']] <- dataset[['narrow']][,.(person_id,date,meaning_of_event)]
    dataset[['possible']] <- dataset[['possible']][,.(person_id,date,meaning_of_event)]
    for (type in c("narrow","possible")) {
      if ( nrow(dataset[[type]]) > 0 ) {
        dataset[[type]] <-  dataset[[type]][!is.na(person_id),]
        dataset[[type]][,type_concept_set := paste0(type)]
      }
      else{
        dataset[[type]][,type_concept_set := NA]
      }
    }  
    if ( nrow(dataset[['narrow']]) == 0 ) { 
      dataset[['broad']] <- dataset[['possible']]
    }
    if ( nrow(dataset[['possible']]) == 0 ) { 
      dataset[['broad']] <- dataset[['narrow']]
    }
    if ( nrow(dataset[['narrow']]) > 0 & nrow(dataset[['possible']]) > 0 ){
      dataset[['broad']] <- as.data.table(rbind(dataset[['narrow']],dataset[['possible']]))
    }
    if ('mo_meaning' %in% colnames(dataset[['broad']])  ){
      dataset[['broad']] <- dataset[['broad']][is.na(meaning_of_event) & !(is.na(mo_meaning)),meaning_of_event:= mo_meaning]
    }
    dataset[['broad']] <- dataset[['broad']][is.na(meaning_of_event),meaning_of_event:="survey"]
    dataset[['broad']] <- dataset[['broad']][,year:=year(date)] 
    dataset[['broad']] <- dataset[['broad']][,n:=1] 
  
    OUTCOME_detailed_components <- MergeFilterAndCollapse(listdatasetL =  list(dataset[['broad']]),
                                                          
                                                          condition= paste0("date>=study_entry_date - 365 & year >=",firstyear, "& year <=",secondyear),
                                                          key = "person_id",
                                                          datasetS = COHORT_DATES,
                                                          saveintermediatedataset=F,
                                                          strata=c("person_id","type_concept_set","meaning_of_event","year"),
                                                          summarystat = list(
                                                            list(c("max"),"n","has_component")
                                                          )
    )
  
    OUTCOME_detailed_components <- OUTCOME_detailed_components[year == as.numeric(firstyear),todrop:=1] 
    OUTCOME_todrop <- OUTCOME_detailed_components[todrop==1,]
    OUTCOME_todrop <- OUTCOME_todrop[,.(person_id,todrop)]
    OUTCOME_detailed_components <- OUTCOME_detailed_components[is.na(todrop),]
    OUTCOME_detailed_components <- OUTCOME_detailed_components[,component := paste0(type_concept_set,'_',meaning_of_event)]
    OUTCOME_detailed_components <- OUTCOME_detailed_components[,.(person_id,component,has_component)]
    
    if ( nrow(OUTCOME_detailed_components) > 0 ) {
      OUTCOME_reshaped <- dcast(OUTCOME_detailed_components, person_id ~ component, value.var = "has_component")
    }
    else {
      OUTCOME_reshaped <- COHORT_TMP
    }
    OUTCOME_reshaped[is.na(OUTCOME_reshaped)] <- 0
    OUTCOME_merged <- merge(COHORT_TMP,OUTCOME_reshaped,by = "person_id", all.x =T)
    OUTCOME_merged <- merge(OUTCOME_merged,OUTCOME_todrop,by = "person_id", all.x =T)
    OUTCOME_merged[is.na(OUTCOME_merged)] <- 0
    OUTCOME_merged <- OUTCOME_merged[, person_id := NULL]
    OUTCOME_aggregated <- OUTCOME_merged[, .(.N), by = c(colnames(OUTCOME_merged))]
    nameobject <- paste0("QC_all_components_",OUTCOME)

    OUTCOME_components <- OUTCOME_aggregated
    assign(nameobject, OUTCOME_aggregated)
    thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
    fwrite(get(nameobject),file=paste0(thisdirexp,paste0(nameobject,".csv")))


  assign(nameobject, OUTCOME_components)
  save(nameobject,file=paste0(dirtemp,paste0(nameobject,".RData")),list = nameobject)
  rm(nameobject,list=nameobject)
  }
  
  
  
  
  rm(namedatasetnarrow , list = namedatasetnarrow)
  rm(namedatasetpossible , list = namedatasetpossible)
  rm(OUTCOME_detailed_components, OUTCOME_todrop, OUTCOME_reshaped, OUTCOME_merged, OUTCOME_aggregated,OUTCOME_components)
  
}

rm(study_population, COHORT_TMP, COHORT_to_be_used,COHORT_DATES, COHORT_dates)

# mask small counts
for (subpop in subpopulations_non_empty){
  col <- c('N')
  temp<-paste0(col,"=5")
  temp2<-paste("c(",paste(temp, collapse = ','),")")
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  thisdirsmallcountsremoved <- ifelse(this_datasource_has_subpopulations == FALSE,dirsmallcountsremoved,dirsmallcountsremovedsubpop[[subpop]])
  suppressWarnings(
    DRE_Treshold(
      Inputfolder = thisdirexp,
      Outputfolder = thisdirsmallcountsremoved,
      Delimiter = ",",
      Varlist = c(eval(parse(text=(temp2)))),
      FileContains = "QC_all_components_"
    )
  )
}
