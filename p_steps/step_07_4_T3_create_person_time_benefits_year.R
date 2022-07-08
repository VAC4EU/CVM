# COUNT PERSON TIME PER COVID
#-----------------------------------------------
#To estimate the weekly incidence rates of COVID-19 (overall and by severity level) in 2020 by data source

# input: D3_vaxweeks_including_not_vaccinated, D3_outcomes_covid ,list_outcomes_observed_COVID
# output: D4_persontime_benefit_year (exported to csv)

print("COUNT PERSON TIME PER COVID by year benefits")

D4_persontime_benefit_year <- vector(mode = 'list')

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  start_persontime_studytime = "20200101"
  
  load(paste0(dirpargen,"list_outcomes_observed_COVID",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_outcomes_covid",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_vaxweeks_including_not_vaccinated",suffix[[subpop]],".RData"))
  
  list_outcomes<-get(paste0("list_outcomes_observed_COVID", suffix[[subpop]]))
  events_ALL_OUTCOMES<-unique(get(paste0("D3_outcomes_covid", suffix[[subpop]])))
  study_population<-get(paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))

  endyear<- substr(study_population[,max(end_date_of_period)], 1, 4)
  end_persontime_studytime<-as.character(paste0(endyear,"1231"))
  
  Output_file<-CountPersonTime(
    Dataset_events = events_ALL_OUTCOMES,
    Dataset = study_population,
    Person_id = "person_id",
    Start_study_time = start_persontime_studytime,
    End_study_time = end_persontime_studytime,
    Start_date = "start_date_of_period",
    End_date = "end_date_of_period",
    #Birth_date = "date_of_birth",
    Strata = c("Dose","sex", "ageband_at_study_entry","type_vax","week_fup", "CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD",
               "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
    Name_event = "name_event",
    Date_event = "date_event",
    #Age_bands = c(0,19,29,39,49,59,69,79),
    Increment="year",
    Outcomes_nrec = list_outcomes, 
    # Unit_of_age = "year",
    # include_remaning_ages = T,
    Aggregate = T
  )

persontime_benefit_year <- Output_file

thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
fwrite(persontime_benefit_year,file=paste0(thisdirexp,"D4_persontime_benefit_year.csv"))

nameoutput<-paste0("D4_persontime_benefit_year",suffix[[subpop]])
assign(nameoutput,persontime_benefit_year)
save(nameoutput,file=paste0(diroutput,nameoutput,".RData"),list=nameoutput)

rm(list=paste0("D4_persontime_benefit_year",suffix[[subpop]]))
rm(list=paste0("D3_vaxweeks_including_not_vaccinated", suffix[[subpop]]))
rm(list=paste0("list_outcomes_observed_COVID", suffix[[subpop]]))
rm(list=paste0("D3_outcomes_covid", suffix[[subpop]]))

}


for (subpop in subpopulations_non_empty){
  tempname<-paste0("D4_persontime_benefit_year")
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  assign(tempname,fread(paste0(thisdirexp,tempname,".csv")))
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
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
      FileContains = "D4_persontime_benefit_year"
    )
  )
  rm(list=tempname)
}

rm(persontime_benefit_year,events_ALL_OUTCOMES,study_population,list_outcomes)
