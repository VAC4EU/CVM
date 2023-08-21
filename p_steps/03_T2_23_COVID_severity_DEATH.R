# CREATE COMPONENTS FOR COVID SEVERITY - DEATH
#-----------------------------------------------
# input: D4_study_population, I_COVID19DX_AESI_narrow,  D3_covid_episodes, COVID_symptoms
# output: D3_covid_severity_components_DEATH, QC_covid_severity_components_DEATH

# in this step the list of unique episodes of covid is labelled with the components indicating that the episode ended in death. The list of components is data source-specific, relevant parameters are described below and specified in 07_algorithms

# data sources using death within 56 days
# datasources_death_after_covid 

# data sources including access to death from covid registry (to be handled in a data source-tailored manner below)
# datasources_death_from_covid_registry 
 

print("CREATE COMPONENTS FOR COVID SEVERITY - DEATH")

emptydataset <- emptydataset <- data.table(
  person_id = "",
  date = NA,
  origin_component = ""
)
emptydataset <- emptydataset[, date := lubridate::ymd(date)]
emptydataset <- emptydataset[person_id!= "",]

load(paste0(direvents,"D3_events_DEATH.RData"))
death <- as.data.table(D3_events_DEATH)

# "covid_severity_1" "covid_severity_2" "covid_severity_3" "covid_severity_4" "covid_severity_5"


# OVERALL STRATEGY 
# 1 rbind all files that imply death 
# 2 associate each record to the corresponding episode (from the date of the episode until the date of the next episode)
# 3 reshape the dataset to obtain one record per episode, labelled with each component indicating that the episode was lethal
# 4 create a descriptive of all the components that contributed to classify the episode as lethal

for (subpop in subpopulations_non_empty) { 
  #---------------------------------
  #---------------------------------
  # 1 rbind all files that imply death
  #---------------------------------
  #---------------------------------
  

  load(paste0(dirtemp,"D3_covid_episodes",suffix[[subpop]],".RData"))
  COVID_episodes <- as.data.table(get(paste0("D3_covid_episodes",suffix[[subpop]])))


  #-------------------------
  # inizialise the dataset of components for severity level 'covid death'
  components_covid_death <- emptydataset

  #-------------------------
  # death after covid

  if (thisdatasource %in% datasources_death_after_covid){
    setnames(death,c("date"),c("date_event"))
    death <- death[,.(person_id,date_event)]
    death <-  merge(COVID_episodes,death, by = "person_id",allow.cartesian = TRUE)[date_event >= date & date_event <= date + 56 , ]
    death <- death[,origin_component := "death_within_56_days"]
    death <- death[,.(person_id,date,origin_component)]
    components_covid_death <- rbind(components_covid_death, death, fill = TRUE)[,.(person_id, date, origin_component)]
  }  

  #-------------------------
  # death from covid registry (data source-tailored)
  
  if (thisdatasource %in% datasources_death_from_covid_registry){
    death_from_covid_registry <- emptydataset
    
    if (thisdatasource %in% c("TEST","ARS","CASERTA")){
      load(paste0(diritemsets,"COVID_symptoms.RData"))
      death_from_covid_registrythisdatasource <- COVID_symptoms[so_source_value == 'Deceduto' | so_source_value == 'DECEDUTO',]
      if (this_datasource_has_subpopulations == TRUE){  
        death_from_covid_registrythisdatasource <- death_from_covid_registrythisdatasource[eval(parse(text = select_in_subpopulationsSO[[subpop]])),]
      death_from_covid_registry = rbind(death_from_covid_registry,death_from_covid_registrythisdatasource, fill = TRUE)
      rm(COVID_symptoms, death_from_covid_registrythisdatasource)
      }
    }
    death_from_covid_registry <- death_from_covid_registry[,origin_component := "death_from_covid_registry"]
    components_covid_death <- rbind(components_covid_death, death_from_covid_registry, fill = TRUE)[,.(person_id, date, origin_component)]
    
    rm(death_from_covid_registry)
    
  }
  
  #---------------------------------
  #---------------------------------
  # 2 MERGE ALL COMPONENTS TO LIST OF EPISODES
  #---------------------------------
  #---------------------------------

  COVID_episodes <- COVID_episodes[,date_next_record := shift(date, n = 1, fill = NA, type = c("lead")), by = "person_id"]
  components_covid_death <- components_covid_death[,.(person_id,date,origin_component)]
  setnames(components_covid_death,c("date"),c("date_component"))
  components_covid_death <- merge(COVID_episodes,components_covid_death, all.x = TRUE, by = "person_id",allow.cartesian = TRUE)[date_component >= date & (date_component < date_next_record | is.na(date_next_record)), ]
  components_covid_death <- unique(components_covid_death[,.(person_id,date,origin_component)])

  #---------------------------------
  #---------------------------------
  # 3 reshape the dataset to obtain one record per episode
  #---------------------------------
  #---------------------------------
  
    
  components_covid_death <- components_covid_death[, component := 1]
  if (nrow(components_covid_death) > 0 ){
    components_covid_death <- dcast(components_covid_death,person_id + date ~ origin_component, value.var = "component", fill = 0 )
  }else{
    components_covid_death <- components_covid_death[,.(person_id,date)]
  }
    
  tempname <- paste0("D3_covid_severity_components_death",suffix[[subpop]])
  assign(tempname,components_covid_death)
  save(list = tempname, file = paste0(dirtemp,tempname,".RData"))

  #---------------------------------
  #---------------------------------
  # 4 count occurrence of each combination of components
  #---------------------------------
  #---------------------------------

  components_covid_death <- as.data.table(components_covid_death)
  components_covid_death <- components_covid_death[, year := year(date)]
  columns_components_covid_death <- colnames(components_covid_death)[colnames(components_covid_death) %not in% c("person_id","date","date_next_record","date_component")]
  components_covid_death <- components_covid_death[, .N, by = columns_components_covid_death]
  
  tempname <- paste0("QC_covid_severity_components_death",suffix[[subpop]])
  assign(tempname,components_covid_death)
  save(list = tempname, file = paste0(diroutput,tempname,".RData"))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(components_covid_death, file = paste0(thisdirexp,tempname,".csv"))
  
  rm(components_covid_death)
  
}
rm(D3_events_DEATH)
