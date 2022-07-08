#------------------------------------------------------------------
# CREATE RISK FACTORS

# input: D3_Vaccin_cohort_no_risk, concept set datasets in DRUGS_conceptssets ("CV","COVCANCER","COVCOPD","COVHIV","COVCKD","COVDIAB","COVOBES","COVSICKLE","IMMUNOSUPPRESSANTS")
# output: D3_Vaccin_cohort_DP


print('CREATE RISK FACTORS (drugs only) at date_vax_1')


# covariate : =1 if there are at least 2 records during 365 days of lookback
# 

D3_study_population_DP <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  load(paste0(dirtemp,"D3_Vaccin_cohort_no_risk",suffix[[subpop]],".RData"))
  study_population <- as.data.table(get(paste0("D3_Vaccin_cohort_no_risk",suffix[[subpop]])))
  
  COHORT_TMP <- study_population[,.(person_id, date_vax1)]
  study_population_DP <- COHORT_TMP
  for (conceptset in DRUGS_conceptssets) {
    load(paste0(dirtemp,conceptset,".RData"))
    output <- MergeFilterAndCollapse(list(get(conceptset)),
                                     condition= "date >= start_lookback & date<=date_vax1",
                                     key = c("person_id"),
                                     datasetS = COHORT_TMP,
                                     additionalvar = list(
                                       list(c("n"),"1","date <= date_vax1")
                                     ),
                                     saveintermediatedataset= F,
                                     strata=c("person_id"),
                                     summarystat = list(
                                       list(c("sum"),"n","howmanyrecords")
                                     )
    )
    output <- output[howmanyrecords > 1 ,namevar := 1]
    output <- output[,.(person_id,namevar)]
    setnames(output,"namevar",paste0(conceptset,"_at_date_vax_1"))
    study_population_DP <- merge(study_population_DP,output,all.x = T, by="person_id")
    study_population_DP[is.na(study_population_DP)] <- 0
    rm(list = conceptset)
    rm(output)
  }

  tempname<-paste0("D3_Vaccin_cohort_DP",suffix[[subpop]])
  assign(tempname,study_population_DP)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  rm(list=paste0("D3_Vaccin_cohort_DP",suffix[[subpop]]))
  rm(list=paste0("D3_Vaccin_cohort_no_risk",suffix[[subpop]]))
  
}


rm(COHORT_TMP, D3_study_population_DP, study_population_DP, study_population)



