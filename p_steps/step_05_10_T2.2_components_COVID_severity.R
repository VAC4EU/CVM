# CREATE ALGORITHMS FOR COVID SEVERITY
#-----------------------------------------------
# input: D4_study_population, D3_events_COVID_narrow, D3_events_DEATH, covid_registry, COVID_symptoms
# output: D3_components_covid_severity

print("CREATE ALGORITHMS FOR COVID SEVERITY")

load(paste0(dirtemp,"covid_registry.RData")) #unique
#load(paste0(dirtemp,"D3_events_COVID_narrow.RData"))  
load(paste0(dirtemp,"D3_events_DEATH.RData")) #unique
load(paste0(dirtemp,"emptydataset"))



#-------------------------------------
# symptoms from covid registry
covid_registry_symptoms <- emptydataset

# if there is a covid registry, define the symptoms an integer 1-5, in a datasource-specific way
if (thisdatasource =='ARS'){
  load(paste0(dirtemp,"COVID_symptoms.RData"))
  covid_registry_symptoms <- COVID_symptoms[,.(person_id,survey_id,so_source_value)][so_source_value == 'Asintomatico' | so_source_value == 'Pauci-sintomatico', covid_registry_symptoms := 1][(so_source_value == 'Lieve' ), covid_registry_symptoms := 2][so_source_value == 'Severo', covid_registry_symptoms := 3][so_source_value == 'Critico', covid_registry_symptoms := 4][so_source_value == 'Deceduto', covid_registry_symptoms := 5]
  rm(COVID_symptoms)
}

if (thisdatasource =='BIFAP'){
  load(paste0(dirtemp,"COVID_ICU.RData"))
  load(paste0(dirtemp,"COVID_hospitalised.RData"))
  covid_registry_symptoms_icu <- COVID_ICU[,.(person_id,survey_id,so_source_column,so_source_value)][so_source_column == 'Ingreso_uci' & so_source_value == '1', covid_registry_symptoms := 4]
  covid_registry_symptoms <- covid_registry_symptoms_icu
  covid_registry_symptoms_hospitalised <- COVID_hospitalised[,.(person_id,survey_id,so_source_column,so_source_value)][so_source_column == 'Ingreso_hospital' & so_source_value == '1', covid_registry_symptoms := 3]
  covid_registry_symptoms <- rbind(covid_registry_symptoms,covid_registry_symptoms_hospitalised, fill = T)
  covid_registry_symptoms <- covid_registry_symptoms[!is.na(covid_registry_symptoms),]
  rm(covid_registry_symptoms_icu)
  rm(covid_registry_symptoms_hospitalised)
  rm(COVID_hospitalised)
  rm(COVID_ICU)
}

if (thisdatasource =='TEST'){
  load(paste0(dirtemp,"COVID_ICU.RData"))
  load(paste0(dirtemp,"COVID_hospitalised.RData"))
  load(paste0(dirtemp,"COVID_symptoms.RData"))
  covid_registry_symptoms <- COVID_symptoms[,.(person_id,survey_id,so_source_value)][so_source_value == 'Asintomatico' | so_source_value == 'Pauci-sintomatico', covid_registry_symptoms := 1][(so_source_value == 'Lieve' ), covid_registry_symptoms := 2][so_source_value == 'Severo', covid_registry_symptoms := 3][so_source_value == 'Critico', covid_registry_symptoms := 4][so_source_value == 'Deceduto', covid_registry_symptoms := 5]
  covid_registry_symptoms_icu <- COVID_ICU[,.(person_id,survey_id,so_source_column,so_source_value)][so_source_column == 'Ingreso_uci' & so_source_value == '1', covid_registry_symptoms := 4]
  covid_registry_symptoms <- rbind(covid_registry_symptoms,covid_registry_symptoms_icu, fill = T)
  covid_registry_symptoms_hospitalised <- COVID_hospitalised[,.(person_id,survey_id,so_source_column,so_source_value)][so_source_column == 'Ingreso_hospital' & so_source_value == '1', covid_registry_symptoms := 3]
  covid_registry_symptoms <- rbind(covid_registry_symptoms,covid_registry_symptoms_hospitalised, fill = T)
  covid_registry_symptoms <- covid_registry_symptoms[!is.na(covid_registry_symptoms),]
  rm(COVID_symptoms)
  rm(covid_registry_symptoms_icu)
  rm(covid_registry_symptoms_hospitalised)
  rm(COVID_hospitalised)
  rm(COVID_ICU)
}

#-----------------------------------
# piece together components of covid severity

D3_components_covid_severity <- vector(mode = 'list')
list_outcomes_observed_COVID <- c()
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D3_events_COVID_narrow",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData"))
  
  events_COVID_narrow<-get(paste0("D3_events_COVID_narrow", suffix[[subpop]])) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
  
  components_covid_severity <- study_population
  
  # define components first_date_covid_narrow
  
  if (nrow(events_COVID_narrow) > 0){
    # condmeaning[['HOSP']]
    covid_dates <- MergeFilterAndCollapse(listdatasetL = list(events_COVID_narrow),
                                          condition = "date == date",
                                          key = c("person_id"),
                                          datasetS = study_population,
                                          additionalvar = list(
                                            list(
                                              c('inhosp'),"1",condmeaning[['HOSP']])
                                          ),
                                          saveintermediatedataset = T,
                                          nameintermediatedataset = paste0(dirtemp,'covid_dates_HOSP'),
                                          strata = c("person_id",'meaning_of_event'),
                                          summarystat = list(                                                                          list(c("min"),"date","first_date_covid_narrow_")
                                          )
    )
    load(paste0(dirtemp,'covid_dates_HOSP.RData'))
    covid_dates_HOSP <- covid_dates_HOSP[inhosp == 1,]
    if (nrow(covid_dates_HOSP) > 0){
      suppressWarnings( covid_dates_HOSP <- unique(covid_dates_HOSP[,first_date_covid_narrow_hosp_discharge:=min(end_date_record),by="person_id"][,.(person_id,first_date_covid_narrow_hosp_discharge)]))
    }
    if (nrow(covid_dates) > 0){
      covid_dates_reshaped <- dcast(covid_dates, person_id ~ meaning_of_event, value.var = "first_date_covid_narrow_")
      covid_dates <- unique(covid_dates[,first_date_covid_narrow:=min(first_date_covid_narrow_),by="person_id"][,.(person_id,first_date_covid_narrow)])
      covid_dates <- merge(covid_dates,covid_dates_reshaped, all.x = T, by="person_id")
    }else{
      covid_dates <- unique(covid_dates[,first_date_covid_narrow:=min(first_date_covid_narrow_),by="person_id"][,.(person_id,first_date_covid_narrow)])
    }
    if (nrow(covid_dates_HOSP) > 0){
      covid_dates <- merge(covid_dates,covid_dates_HOSP, all.x = T, by="person_id")
    }else{ 
      covid_dates <- covid_dates[is.na(person_id),first_date_covid_narrow_hosp_discharge := lubridate::ymd(study_start)]
    }
    components_covid_severity <- merge(components_covid_severity,covid_dates, all.x = T, by="person_id")
    rm(covid_dates_HOSP)
  }else{ # if there is no narrow concepsets, define the dates as missing
    covid_dates <- emptydataset
    components_covid_severity[is.na(person_id), first_date_covid_narrow_hosp_discharge := lubridate::ymd(study_start)][is.na(person_id), first_date_covid_narrow := lubridate::ymd(study_start)]
  }
  
  #-------------------------------- 
  # define first_date_covid_registry
  
  if (nrow(covid_registry) > 0){
    covid_dates_registry <- MergeFilterAndCollapse(listdatasetL = list(covid_registry),
                                                   condition = "!is.na(date)",
                                                   key = c("person_id"),
                                                   datasetS = study_population,
                                                   strata = c("person_id"),
                                                   summarystat = list(                                                                          list(c("min"),"date","first_date_covid_registry")
                                                   )
    )
    components_covid_severity <- merge(components_covid_severity, covid_dates_registry, all.x = T, by="person_id")
  }else{# if there is no covid registry, define first_date_covid_registry as missing
    components_covid_severity[is.na(person_id), first_date_covid_registry := lubridate::ymd(study_start)]
  }
  
  #-------------------------------------------------------
  # symptoms from covid registry
  if (nrow(covid_registry_symptoms) > 0){
    components_covid_severity <- merge(components_covid_severity, covid_registry_symptoms, all.x = T, by="person_id")
  }else{# if there is no covid registry, define the symptoms as a missing integer
    components_covid_severity[is.na(person_id),covid_registry_symptoms:=1]
  }
  
  
  #-------------------------------------------------------
  # define components identifying symptoms within 30 days from covid inception, either from registry or from covid_narrow records
  
  symptoms_within_covid_narrow_date <- vector(mode = 'list')
  symptoms_within_registry_date <- vector(mode = 'list')
  listcomplications <- vector(mode = 'list')
  # listcomplications[['Infection']] <- c("Pneumonia","Bronchitis","LowerRespiratoryTractInfection","Hypoxemia")
  # listcomplications[['Respiratory']] <- c("AcuteRespiratoryFailure", "RespiratoryFailure", "AcuteRespiratoryDistress", "RespiratoryDistressSyndrome")
  # listcomplications[['MechanicalVentilation']] <- c('MechanicalVentilation')
  listcomplications[['Infection']] <- c("COVIDSYMPTOM")
  listcomplications[['Respiratory']] <- c("ARD_narrow","ARD_possible")
  listcomplications[['MechanicalVentilation']] <- c('MechanicalVent')
  
  for (symptoms in c("MechanicalVentilation","Infection","Respiratory")){
    conceptsetdatasets <- emptydataset
    for (complication in listcomplications[[symptoms]]){
      load(paste0(dirtemp,complication,'.RData'))
      temp <- get(complication)
      if (this_datasource_has_subpopulations == TRUE){
        temp <- temp[eval(parse(text = select_in_subpopulationsEVENTS[[subpop]]))]
      }
      if (nrow(temp) > 0){
        temp <- temp[,.(person_id,date)]
        conceptsetdatasets <- rbind(conceptsetdatasets,temp,fill = T)
      }
      rm(complication,list=complication)
    }
    conceptsetdatasets <- conceptsetdatasets[,n:=1]
    
    if (nrow(conceptsetdatasets) > 0 & nrow(covid_dates) > 0){
      selection <- "!is.na(first_date_covid_narrow) & date >= first_date_covid_narrow & date <= first_date_covid_narrow + 30 "
      componentconceptsetdatasets_within_covid_narrow_date <- MergeFilterAndCollapse(listdatasetL = list(conceptsetdatasets),                                                 condition = selection,
                                                                                     key = c("person_id"),
                                                                                     datasetS = covid_dates,
                                                                                     strata = c("person_id"),
                                                                                     summarystat = list(                                                                          list(c("max"),"n",paste0(symptoms,"_within_covid_narrow_date"))
                                                                                     )
      )
      if(nrow(componentconceptsetdatasets_within_covid_narrow_date) >0){ components_covid_severity <- merge(components_covid_severity,componentconceptsetdatasets_within_covid_narrow_date, all.x = T, by="person_id")
      components_covid_severity <- components_covid_severity[is.na(get(paste0(symptoms,"_within_covid_narrow_date"))),paste0(symptoms,"_within_covid_narrow_date") := 0]
      rm(componentconceptsetdatasets_within_covid_narrow_date)
      }else{# if there is no componentconceptsetdatasets_within_covid_narrow_date, define the symptoms as a missing integer
        components_covid_severity <- components_covid_severity[is.na(person_id),paste0(symptoms,"_within_covid_narrow_date") := 0]
      }
    }else{# if there is no componentconceptsetdatasets_within_registry_date, define the symptoms as a missing integer
      components_covid_severity <- components_covid_severity[is.na(person_id),paste0(symptoms,"_within_covid_narrow_date") := 0]
    }
    
    if (nrow(conceptsetdatasets) > 0 & nrow(covid_registry) > 0){
      selection <- "!is.na(first_date_covid_registry) & date >= first_date_covid_registry & date <= first_date_covid_registry + 30 "
      componentconceptsetdatasets_within_registry_date <- MergeFilterAndCollapse(listdatasetL = list(conceptsetdatasets),
                                                                                 condition = selection,
                                                                                 key = c("person_id"),
                                                                                 datasetS = covid_dates_registry,
                                                                                 strata = c("person_id"),
                                                                                 summarystat = list(                                                                          list(c("max"),"n",paste0(symptoms,"_within_registry_date"))
                                                                                 )
      )
      if(nrow(componentconceptsetdatasets_within_registry_date) >0){
        components_covid_severity <- merge(components_covid_severity,componentconceptsetdatasets_within_registry_date, all.x = T, by="person_id")
        components_covid_severity <- components_covid_severity[is.na(get(paste0(symptoms,"_within_registry_date"))),paste0(symptoms,"_within_registry_date") := 0]
        rm(componentconceptsetdatasets_within_registry_date)
      }else{# if there is no componentconceptsetdatasets_within_registry_date, define the symptoms as a missing integer
        components_covid_severity <- components_covid_severity[is.na(person_id),paste0(symptoms,"_within_registry_date") := 0]
      }
    }else{# if there is no componentconceptsetdatasets_within_registry_date, define the symptoms as a missing integer
      components_covid_severity <- components_covid_severity[is.na(person_id),paste0(symptoms,"_within_registry_date") := 0]
    }
  }
  
  
  
  #-------------------------------------
  # death
  
  if (nrow(conceptsetdatasets) > 0 & nrow(covid_dates) > 0){
    selection <- "!is.na(first_date_covid_narrow) & date <= first_date_covid_narrow + 60 "
    death_after_covid_narrow_date <- MergeFilterAndCollapse(listdatasetL =                                                  list(D3_events_DEATH),                                                 condition = selection,
                                                            key = c("person_id"),
                                                            datasetS = covid_dates,
                                                            additionalvar = list(
                                                              list(c('n'),"1")
                                                            ),
                                                            strata = c("person_id"),
                                                            summarystat = list(                                                                          list(c("max"),"n","death_after_covid_narrow_date")
                                                            )
    )
    components_covid_severity <- merge(components_covid_severity, death_after_covid_narrow_date, all.x = T, by="person_id")
    components_covid_severity <- components_covid_severity[is.na(death_after_covid_narrow_date),death_after_covid_narrow_date := 0]
    rm(death_after_covid_narrow_date)
  }else{# if there is no conceptsetdataset, define death_after_covid_narrow_date as a missing integer
    components_covid_severity <- components_covid_severity[is.na(person_id),death_after_covid_narrow_date := 0]
  }
  
  if (nrow(conceptsetdatasets) > 0 & nrow(covid_dates) > 0){
    selection <- "!is.na(first_date_covid_narrow_hosp_discharge) & date == first_date_covid_narrow_hosp_discharge "
    death_at_hosp_discharge <- MergeFilterAndCollapse(listdatasetL =                                                  list(D3_events_DEATH),                                                 condition = selection,
                                                      key = c("person_id"),
                                                      datasetS = covid_dates,
                                                      additionalvar = list(
                                                        list(c('n'),"1")
                                                      ),
                                                      strata = c("person_id"),
                                                      summarystat = list(                                                      list(c("max"),"n","covid_death_discharge")
                                                      )
    )
    components_covid_severity <- merge(components_covid_severity, death_at_hosp_discharge, all.x = T, by="person_id")
    components_covid_severity <- components_covid_severity[is.na(covid_death_discharge),covid_death_discharge := 0]
    rm(death_at_hosp_discharge)
  }else{# if there is no conceptsetdatasets, define covid_death_discharge as a missing integer
    components_covid_severity <- components_covid_severity[is.na(person_id),covid_death_discharge := 0]
  }
  
  if (nrow(covid_registry) > 0){
    selection <- "!is.na(first_date_covid_registry) & date <= first_date_covid_registry + 30 "
    death_after_covid_registry <- MergeFilterAndCollapse(listdatasetL =                                                  list(D3_events_DEATH),                                                 condition = selection,
                                                         key = c("person_id"),
                                                         datasetS = covid_dates_registry,
                                                         additionalvar = list(
                                                           list(c('n'),"1")
                                                         ),
                                                         strata = c("person_id"),
                                                         summarystat = list(                                                                          list(c("max"),"n","death_after_covid_registry")
                                                         )
    )
    components_covid_severity <- merge(components_covid_severity, death_after_covid_registry, all.x = T, by="person_id")
    components_covid_severity <- components_covid_severity[is.na(death_after_covid_registry),death_after_covid_registry := 0] 
    rm(death_after_covid_registry)
  }else{# if there is no covid registry, define covid_death_discharge as a missing integer
    components_covid_severity <- components_covid_severity[is.na(person_id),death_after_covid_registry := 0]
  }
  
  # assign components_covid_severity to the output
  tempname<-paste0("D3_components_covid_severity",suffix[[subpop]])
  assign(tempname,components_covid_severity)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  
  rm(list=paste0("D3_components_covid_severity",suffix[[subpop]]))
  rm(list=paste0("D4_study_population", suffix[[subpop]]))
  rm(list=paste0("D3_events_COVID_narrow", suffix[[subpop]]))
} 

rm(events_COVID_narrow,study_population,D3_events_DEATH, components_covid_severity,temp,covid_dates_registry,covid_dates,covid_registry_symptoms,covid_registry,conceptsetdatasets)
