#------------------------------------------------------------------
# create secondary components 

# input: concept set datasets involved in secondary components, D4_study_population
# output: for each secondary component SECCOMPONENT, D3_eventsSecondary_SECCOMP.RData

print('create secondary components SECCOMPONENTS.')

# load(paste0(dirpargen,"subpopulations_non_empty.RData"))

##for each SECCOMP create components

# create emptyconceptsetsataset

load(paste0(dirtemp,"CAD_possible.RData"))
emptyconceptsetsataset <- CAD_possible[is.na(person_id) & !is.na(person_id)][,conceptsetname := '']
rm(CAD_possible)

for (SECCOMP in SECCOMPONENTS) {
  componentsSECCOMP <- vector(mode="list")
  print(SECCOMP)
  for (subpop in subpopulations_non_empty) {
    print(subpop)
    load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 
    study_population<-get(paste0("D4_study_population", suffix[[subpop]]))

    COHORT_TMP <- as.data.table(study_population)  

    COHORT_TMP <- COHORT_TMP[,.(person_id,study_entry_date)]
    #create datasets A and B to be merged for the secondary component
    datasets_to_be_merged <- vector(mode="list")
    for (ord in c('A','B')){
      temp <- emptyconceptsetsataset
      for (conceptset in  concept_set_seccomp[[SECCOMP]][[ord]]){
        print(paste0(SECCOMP,' - conceptset ',ord, ' ',conceptset))
        load(paste0(dirtemp,conceptset,".RData"))
        
        toadd <- get(conceptset)
        # delete records whose meanings are not to be used for this concepset dataset in this datasource
        if (thisdatasource %in% datasources_with_specific_algorithms){
          toadd <- toadd[eval(parse(text = selection_meanings_from_OUTCOME[[thisdatasource]][[conceptset]])),]
        }
        # delete records whose meanings are not observed in this whole subpopulation
        if (this_datasource_has_subpopulations == TRUE){ 
          toadd <- toadd[eval(parse(text = select_in_subpopulationsEVENTS[[subpop]])),]
        }
        toadd <- toadd[,conceptsetname := conceptset]
        temp <- as.data.table(rbind(temp,toadd,fill = TRUE))
        rm(toadd)
      
        rm(conceptset,list = conceptset)
        }
      temp <- temp[,.(person_id,date,end_date_record,codvar, event_record_vocabulary,meaning_of_event,conceptsetname)]
      for (col in c('date','end_date_record','codvar', 'event_record_vocabulary','meaning_of_event','conceptsetname')){
        setnames(temp, col,paste0(col,ord) )
      }
      datasets_to_be_merged[[ord]] <- temp
      rm(temp)
    }
    
    # select according to the rule of the component
    selection <- paste0("!is.na(dateA) & !is.na(dateB) & ", selectionrule_direction_seccomp[[direction_seccomp[[SECCOMP]]]])
    
    # merge datasets A and B
    unique_A_AND_B_timeframe <- MergeFilterAndCollapse(
        listdatasetL = list(datasets_to_be_merged[['A']]),
        datasetS = datasets_to_be_merged[['B']],
        condition = selection,
        key = c("person_id"),
        typemerge = 2,
        sorting= c("person_id","dateA"),
        saveintermediatedataset = T,
        nameintermediatedataset = paste0(dirtemp,'tempfile'),
        strata = c("person_id"),
        summarystat = list(
          list(c('first'),'dateA','date_event')
        )
      )
    
    load(paste0(dirtemp,'tempfile.RData') )
    all_A_AND_B_timeframe <- tempfile
      
    if (rule_seccomp[[SECCOMP]] == "AND NOT"){
      if (nrow(datasets_to_be_merged[['A']]) > 0){
          all_component_A <- datasets_to_be_merged[['A']]
          all_component_A <- all_component_A[all_component_A[, .I[sample(.N,1)] , by = c("person_id","dateA")]$V1]
          if (nrow(all_A_AND_B_timeframe)>0){ 
            to_exclude <- all_A_AND_B_timeframe[,.(person_id,dateA,dateB)]
            listevents <-  merge(to_exclude,all_component_A, by = c("person_id","dateA"), all.y = T)[(eval(parse(text = 'is.na(dateB)'))),]
            rm(to_exclude)
          }
          else{
            listevents <-  all_component_A
          }
          rm(all_component_A)
        }
        else{
         listevents <- datasets_to_be_merged[['A']]
        }
      }
    if (rule_seccomp[[SECCOMP]] == "AND"){
        listevents <- all_A_AND_B_timeframe
      }
      
    listevents <- setnames(listevents, 'dateA','date') 

    components <- MergeFilterAndCollapse(
        listdatasetL = list(listevents),
        condition = "date >= study_entry_date - 365 ",
        key = c("person_id"),
        datasetS = COHORT_TMP,
        saveintermediatedataset = T,
        nameintermediatedataset = paste0(dirtemp,'tempfile'),
        strata = c("person_id"),
        summarystat =  list(
          list(c('min'),'date','date_event')
        )
      )
      
    load(paste0(dirtemp,'tempfile.RData') )
      
    componentsSECCOMP<- tempfile 

  nameobjectSECCOMP <- paste0('D3_eventsSecondary',"_",SECCOMP,suffix[[subpop]])
  assign(nameobjectSECCOMP,componentsSECCOMP)
  save(nameobjectSECCOMP,file=paste0(dirtemp,paste0(nameobjectSECCOMP,".RData")),list = nameobjectSECCOMP)
  rm(nameobjectSECCOMP,list = nameobjectSECCOMP)
  rm(datasets_to_be_merged,componentsSECCOMP,tempfile,COHORT_TMP,components,listevents,all_A_AND_B_timeframe,unique_A_AND_B_timeframe,study_population)
  
  rm(list=paste0("D4_study_population", suffix[[subpop]]))
  }
}
  

