# COUNT CODES PER ALL OUTCOMES IN THE STUDY POPULATION

# input: D4_study_population, D3_events_ALL_OUTCOMES 
# output: QC_code_counts_in_study_population_OUTCOME_YYYY (exported to csv)
# parameters: list_outcomes_observed_only_diagnosis, list_outcomes_observed_for_QC

#-----------------------------------------------

print("CREATE CODE COUNT FOR ALL OUTCOMES PER YEAR IN THE STUDY POPULATION")

start_persontime_studytime = as.character(paste0(study_years[1],"0101"))
end_persontime_studytime = as.character(paste0(study_years[length(study_years)],"1231"))

FirstJan <- list()
conditionYear <- list()
for (year in study_years) {
  FirstJan[[year]] <- as.Date(as.character(paste0(year,"0101")), date_format)
  conditionYear[[year]] <- paste0("date_event >= study_entry_date & date_event <= as.Date('",FirstJan[[year]],"')+365 & date_event >= as.Date('",FirstJan[[year]],"')")
}

  list_outcomes <- c()
  for (subpop in subpopulations_non_empty) { 
    
    load(paste0(dirtemp,"list_outcomes_observed_for_QC",suffix[[subpop]],".RData"))

    list_outcomes <- c(list_outcomes, get(paste0("list_outcomes_observed_for_QC",suffix[[subpop]])))
  }
  list_outcomes <- unique(list_outcomes)


for (outcome in list_outcomes) {
  print(outcome)
  outcome_code_counts <- vector(mode = 'list')
  for (subpop in subpopulations_non_empty) { 
    
    load(paste0(dirtemp,"D3_events_ALL_OUTCOMES",suffix[[subpop]],".RData"))
    load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 

    study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
      events_ALL_OUTCOMES <-get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]])) 
      admissible_outcomes <- get(paste0("list_outcomes_observed_for_QC", suffix[[subpop]])) 

    if (outcome %in% admissible_outcomes){
      for (year in study_years) {
        print(paste(outcome,year,subpop))
        nameobject <- paste0("QC_code_counts_in_study_population","_",outcome,"_",year)
    temp <- MergeFilterAndCollapse(list(events_ALL_OUTCOMES[name_event == outcome,]),
                                          key = 'person_id',
                                          condition = conditionYear[[year]],
                                          datasetS = study_population,
                                          additionalvar = list(list(c("n"),1)),
                                          strata=c( "meaning_of_first_event", "coding_system_of_code_first_event", "code_first_event"),
                                          summarystat = list(
                                            list(c("count"),"n")
                                            )
                                          )
    
    assign(nameobject, temp)
    thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
    fwrite(temp, file = paste0(thisdirexp,nameobject,".csv"))

    save(nameobject, file = paste0(diroutput,nameobject,".RData"),list = nameobject )
    rm(nameobject, list = nameobject )
      }
    }
      rm(list=paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
      rm(list=paste0("D4_study_population", suffix[[subpop]]))
  }
}

for (subpop in subpopulations_non_empty){
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  thisdirsmallcountsremoved <- ifelse(this_datasource_has_subpopulations == FALSE,dirsmallcountsremoved,dirsmallcountsremovedsubpop[[subpop]])
  suppressWarnings(
    DRE_Treshold(
      Inputfolder = thisdirexp,
      Outputfolder = thisdirsmallcountsremoved,
      Delimiter = ",",
      Varlist = c(count_n=5),
      FileContains = "QC_code_counts_in_"
    )
  )
}
rm(study_population,events_ALL_OUTCOMES,temp,outcome_code_counts)
