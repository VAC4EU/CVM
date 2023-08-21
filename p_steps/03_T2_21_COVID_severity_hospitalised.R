# CREATE COMPONENTS FOR COVID SEVERITY - HOSPITALISED
#-----------------------------------------------
# input: D4_study_population, D3_covid_episodes, conceptsetdataset: I_COVID19DX_AESI_narrow, itemsetdatasets: COVID_symptoms, propmtsetdataset: hospitalisation_automatically_referred_to_PC
# output: D3_covid_severity_components_hospitalisation, QC_covid_severity_components_hospitalisation

# in this step the list of unique episodes of covid is labelled with the components indicating that the episode was hospitalised. The list of components is data source-specific, relevant parameters are described below and specified in 07_algorithms

# data sources including records of hospitalisations with a covid diagnosis
# datasources_hosp_due_to_covid 

# data sources including records of hospitalizations irrespective of diagnosis
# datasources_hosp_after_covid 

# data sources including hospitalization from covid registry (to be handled in a data source-tailored manner below)
# datasources_hosp_from_covid_registry


print("CREATE COMPONENTS FOR COVID SEVERITY - HOSPITALISED")
 
emptydataset <- data.table(
  person_id = "",
  date = NA,
  origin_component = ""
)
emptydataset <- emptydataset[, date := lubridate::ymd(date)]
emptydataset <- emptydataset[person_id!= "",]

# set selection of records whose meanings are not appropriate for a specific datasource 
selectionOUTCOME <- paste0('!(', select_meanings_AESI[[thisdatasource]], ')')

# load recordings of covid diagnoses and remove inappropriate meanings
load(paste0(dirconceptsets,"I_COVID19DX_AESI_narrow.RData"))
dia_COVID_narrow <- I_COVID19DX_AESI_narrow
dia_COVID_narrow <- dia_COVID_narrow[eval(parse(text = selectionOUTCOME)),]



# OVERALL STRATEGY 
# 1 rbind all files that imply hospitalisation
# 2 associate each record to the corresponding episode (from the date of the episode until the date of the next episode)
# 3 reshape the dataset to obtain one record per episode, labelled with each component indicating that the episode was hospitalised
# 4 create a descriptive of all the components that contributed to classify the episode as hospitalised

for (subpop in subpopulations_non_empty) { 

  #---------------------------------
  #---------------------------------
  # 1 rbind all files that imply hospitalisation
  #---------------------------------
  #---------------------------------
  

  load(paste0(dirtemp,"D3_covid_episodes",suffix[[subpop]],".RData"))
  COVID_episodes <- as.data.table(get(paste0("D3_covid_episodes",suffix[[subpop]])))
  
  # load(paste0(dirtemp,"D3_events_COVID_narrow",suffix[[subpop]],".RData"))
  # dia_COVID_narrow <- as.data.table(get(paste0("D3_events_COVID_narrow",suffix[[subpop]])))
  
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData"))
  
  
  #-------------------------
  # inizialise the dataset of components for severity level 'covid hospitalised'
  components_covid_hospitalisations <- emptydataset
  
  
  #-------------------------
  # hospitalisations due to covid
  if (thisdatasource %in% datasources_hosp_due_to_covid){
    dia_COVID_narrow_hosp <-  dia_COVID_narrow[eval(parse(text = condmeaning[["HOSP"]])),]
    dia_COVID_narrow_hosp <- dia_COVID_narrow_hosp[,origin_component := paste0("covid_narrow_hosp_m_",meaning_renamed)]
    components_covid_hospitalisations <- rbind(components_covid_hospitalisations, dia_COVID_narrow_hosp, fill = TRUE)[,.(person_id, date, origin_component)]
  }
  
  
  #-------------------------
  # hospitalisations after covid

  if (thisdatasource %in% datasources_hosp_after_covid){
    load(paste0(dirpromptsets,"hospitalisation_automatically_referred_to_PC.RData"))
    setnames(hospitalisation_automatically_referred_to_PC,c("date"),c("date_hosp"))
    hospitalisation_automatically_referred_to_PC <- hospitalisation_automatically_referred_to_PC[,.(person_id,date_hosp)]
    hosp_after_covid <-  merge(COVID_episodes,hospitalisation_automatically_referred_to_PC, by = "person_id",allow.cartesian = TRUE)[date_hosp >= date & date_hosp <= date + 14 , ]
    hosp_after_covid <- hosp_after_covid[,origin_component := "hospitalisation_after_COVID"]
    components_covid_hospitalisations <- rbind(components_covid_hospitalisations, hosp_after_covid, fill = TRUE)[,.(person_id, date, origin_component)]
    rm(hospitalisation_automatically_referred_to_PC)
  }
  
  
  #-------------------------
  # hospitalisations from covid registry (data source-tailored)
  
  if (thisdatasource %in% datasources_hosp_from_covid_registry){
    hosp_from_covid_registry <- emptydataset
    
    if (thisdatasource %in% c("TEST","ARS","CASERTA")){
      load(paste0(diritemsets,"COVID_symptoms.RData"))
      hosp_from_covid_registrythisdatasource <- COVID_symptoms[,.(person_id,date,so_source_value)][so_source_value == 'Severo' | so_source_value == 'SEVERO',]
      hosp_from_covid_registry = rbind(hosp_from_covid_registry,hosp_from_covid_registrythisdatasource, fill = TRUE)
      rm(COVID_symptoms, hosp_from_covid_registrythisdatasource)
    }
    
    if (thisdatasource %in% c("TEST","BIFAP")){
      load(paste0(diritemsets,"COVID_hospitalised.RData"))
      hosp_from_covid_registrythisdatasource <- COVID_hospitalised[,.(person_id,date,so_source_column,so_source_value)][so_source_column == 'Ingreso_hospital' & so_source_value == '1',]
      hosp_from_covid_registry = rbind(hosp_from_covid_registry,hosp_from_covid_registrythisdatasource, fill = TRUE)
      rm(COVID_hospitalised, hosp_from_covid_registrythisdatasource)
    }
    
    hosp_from_covid_registry <- hosp_from_covid_registry[,origin_component := "hospitalisation_from_covid_registry"]
    components_covid_hospitalisations <- rbind(components_covid_hospitalisations, hosp_from_covid_registry, fill = TRUE)[,.(person_id, date, origin_component)]
    
    rm(hosp_from_covid_registry)
    
  }
  
  #---------------------------------
  #---------------------------------
  # 2 MERGE ALL COMPONENTS TO LIST OF EPISODES
  #---------------------------------
  #---------------------------------

  COVID_episodes <- COVID_episodes[,date_next_record := shift(date, n = 1, fill = NA, type = c("lead")), by = "person_id"]
  components_covid_hospitalisations <- components_covid_hospitalisations[,.(person_id,date,origin_component)]
  setnames(components_covid_hospitalisations,c("date"),c("date_component"))
  components_covid_hospitalisations <- merge(COVID_episodes,components_covid_hospitalisations, all.x = TRUE, by = "person_id",allow.cartesian = TRUE)[date_component >= date & (date_component < date_next_record | is.na(date_next_record)), ]
  components_covid_hospitalisations <- unique(components_covid_hospitalisations[,.(person_id,date,origin_component)])

  #---------------------------------
  #---------------------------------
  # 3 reshape the dataset to obtain one record per episode
  #---------------------------------
  #---------------------------------
  
    
  components_covid_hospitalisations <- components_covid_hospitalisations[, component := 1]
  if (nrow(components_covid_hospitalisations) > 0 ){
    components_covid_hospitalisations <- dcast(components_covid_hospitalisations,person_id + date ~ origin_component, value.var = "component", fill = 0 )
  }else{
    components_covid_hospitalisations <- components_covid_hospitalisations[,.(person_id,date)]
  }
  
    
  tempname <- paste0("D3_covid_severity_components_hospitalisation",suffix[[subpop]])
  assign(tempname,components_covid_hospitalisations)
  save(list = tempname, file = paste0(dirtemp,tempname,".RData"))

  #---------------------------------
  #---------------------------------
  # 4 count occurrence of each combination of components
  #---------------------------------
  #---------------------------------

  components_covid_hospitalisations <- as.data.table(components_covid_hospitalisations)
  components_covid_hospitalisations <- components_covid_hospitalisations[, year := year(date)]
  columns_components_covid_hospitalisations <- colnames(components_covid_hospitalisations)[colnames(components_covid_hospitalisations) %not in% c("person_id","date","date_next_record","date_component")]
  components_covid_hospitalisations <- components_covid_hospitalisations[, .N, by = columns_components_covid_hospitalisations]
  
  tempname <- paste0("QC_covid_severity_components_hospitalisation",suffix[[subpop]])
  assign(tempname,components_covid_hospitalisations)
  save(list = tempname, file = paste0(diroutput,tempname,".RData"))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(components_covid_hospitalisations, file = paste0(thisdirexp,tempname,".csv"))
  
  
    
  rm(components_covid_hospitalisations)
  
}
rm(dia_COVID_narrow)