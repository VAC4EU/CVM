###################################################################
# DESCRIBE THE PARAMETERS OF SUBPOPULATIONS RESTRICTING
###################################################################

# datasources_with_subpopulations lists the datasources where some meanings of events should be excluded during some observation periods, associated with some op_meanings

datasources_with_subpopulations <- c("BIFAP", "TEST")


this_datasource_has_subpopulations <- ifelse(thisdatasource %in% datasources_with_subpopulations,TRUE,FALSE) 


# subpopulations associates with each datasource the label of its subpopulations
subpopulations <- list()
 


# op_meaning_sets labels each group of op_meaning that are relevant
op_meaning_sets <- vector(mode="list")

# op_meaning_sets separates all the op_meanings available in OBSERVATION_PERIODS in op_meaning_sets (3-levels list: the datasource, and the op_meaning_sets)

op_meanings_list_per_set <- vector(mode="list")

# op_meaning_sets associates to each subpopulation the corresponding overlap of op_meaning_sets
op_meaning_sets_in_subpopulations <- vector(mode="list")

# exclude_meaning_renamed associates to each subpopulation the corresponding meaning of events that should not be processed (3-levels list: the datasource, and the subpopulation) 
exclude_meaning_renamed <- vector(mode="list") 
exclude_itemset_of_so <- vector(mode="list") 


# datasource TEST
subpopulations[["TEST"]] = c("HOSP","ER_HOSP")


op_meaning_sets[["TEST"]] <- c("meaningsHOSP","meaningsER")
op_meanings_list_per_set[["TEST"]][["meaningsHOSP"]] <- c("REGION_1_HOSP","REGION_2_HOSP","REGION_3_HOSP")
op_meanings_list_per_set[["TEST"]][["meaningsER"]] <- c("REGION_1_ER","REGION_2_ER","REGION_3_ER")

op_meaning_sets_in_subpopulations[["TEST"]][["HOSP"]] <- c("meaningsHOSP")
op_meaning_sets_in_subpopulations[["TEST"]][["ER_HOSP"]] <- c("meaningsHOSP","meaningsER")

exclude_meaning_renamed[["TEST"]][["ER_HOSP"]] <- c()
exclude_meaning_renamed[["TEST"]][["HOSP"]] <- c("emergency_room_diagnosis")

exclude_itemset_of_so[["TEST"]][["HOSP"]] <- list(list("Covid19_UCI","Ingreso_uci"),list("Covid19_UCI","Fecha_ingreso_uci"))


# # BIFAP
subpopulations[["BIFAP"]] = c("PC","PC_HOSP")

op_meaning_sets[["BIFAP"]] <- c("meaningsPC","meaningsHOSP")
#op_meaning_sets[["BIFAP"]] <- c("meaningsPC","meaningsHOSP","meaningsWITH_ICU")
op_meanings_list_per_set[["BIFAP"]][["meaningsPC"]] <- c("region2_PC","region3_PC","region7_PC","region13_PC","region14_PC")
op_meanings_list_per_set[["BIFAP"]][["meaningsHOSP"]] <- c("region3_HOSP","region7_HOSP","region13_HOSP")
op_meanings_list_per_set[["BIFAP"]][["meaningsWITH_ICU"]] <- c("region2_PC","region3_PC","region7_PC","region14_PC")
op_meanings_list_per_set[["BIFAP"]][["meaningsCOVID"]] <- c("region2_COVID","region3_COVID","region7_COVID","region14_COVID")

op_meaning_sets_in_subpopulations[["BIFAP"]][["PC"]] <- c("meaningsPC")
op_meaning_sets_in_subpopulations[["BIFAP"]][["PC_HOSP"]] <- c("meaningsPC","meaningsHOSP")
op_meaning_sets_in_subpopulations[["BIFAP"]][["WITH_ICU"]] <- c("meaningsWITH_ICU")
op_meaning_sets_in_subpopulations[["BIFAP"]][["PC_COVID"]] <- c("meaningsPC","meaningsCOVID")

exclude_meaning_renamed[["BIFAP"]][["PC"]] <- c("hopitalisation_diagnosis_unspecified","hospitalisation_primary","
hospitalisation_secondary","hospitalisation_secundary")
exclude_meaning_renamed[["BIFAP"]][["PC_COVID"]]<-c("hopitalisation_diagnosis_unspecified","hospitalisation_primary","
hospitalisation_secondary")
exclude_meaning_renamed[["BIFAP"]][["PC_HOSP"]]<-c()
exclude_meaning_renamed[["BIFAP"]][["WITH_ICU"]] <- c("hopitalisation_diagnosis_unspecified","hospitalisation_primary","
hospitalisation_secondary")

exclude_itemset_of_so[["BIFAP"]][["PC"]] <- list(list("Covid19_UCI","Ingreso_uci"),list("Covid19_UCI","Fecha_ingreso_uci"))
exclude_itemset_of_so[["BIFAP"]][["PC_HOSP"]] <- list(list("Covid19_UCI","Ingreso_uci"),list("Covid19_UCI","Fecha_ingreso_uci"))
exclude_itemset_of_so[["BIFAP"]][["WITH_ICU"]] <- c()


# # SIDIAP
# subpopulations[["SIDIAP"]] = c("PC","PC_HOSP")
# op_meaning_sets[["SIDIAP"]] <- c("meaningsPC","meaningsHOSP")
# 
# op_meanings_list_per_set[["SIDIAP"]][["meaningsPC"]] <- c("enlisted_with_GP") 
# op_meanings_list_per_set[["SIDIAP"]][["meaningsHOSP"]] <- c("observed_in_hospital") 
# 
# op_meaning_sets_in_subpopulations[["SIDIAP"]][["PC"]] <- c("meaningsPC")
# op_meaning_sets_in_subpopulations[["SIDIAP"]][["PC_HOSP"]] <- c("meaningsPC","meaningsHOSP")
# 
# exclude_meaning_renamed[["SIDIAP"]][["PC"]]<-c("hospitalisation_primary", "hospitalisation_secondary","hospitalisation_secondar")
# exclude_meaning_renamed[["SIDIAP"]][["PC_HOSP"]]<-c()

# # PHARMO
# subpopulations[["PHARMO"]] = c("PC","HOSP","PC_HOSP")
# 
# op_meaning_sets[["PHARMO"]] <- c("meaningsPC","meaningsHOSP","meaningsPHARMA")
# 
# op_meanings_list_per_set[["PHARMO"]][["meaningsPC"]] <- c("primary_care") 
# op_meanings_list_per_set[["PHARMO"]][["meaningsHOSP"]] <- c("hospitalisation") 
# op_meanings_list_per_set[["PHARMO"]][["meaningsPHARMA"]] <- c("outpatient_pharmacy") 
# 
# 
# op_meaning_sets_in_subpopulations[["PHARMO"]][["PC"]] <- c("meaningsPHARMA","meaningsPC")
# op_meaning_sets_in_subpopulations[["PHARMO"]][["HOSP"]] <- c("meaningsPHARMA","meaningsHOSP")
# op_meaning_sets_in_subpopulations[["PHARMO"]][["PC_HOSP"]] <- c("meaningsPHARMA","meaningsHOSP","meaningsPC")
# 
# 
# exclude_meaning_renamed[["PHARMO"]][["PC"]]<-c("hospital_diagnosis","amb_diagnosis")
# exclude_meaning_renamed[["PHARMO"]][["HOSP"]]<-c("primary_care_event")
# exclude_meaning_renamed[["PHARMO"]][["PC_HOSP"]]<-c()


if (this_datasource_has_subpopulations == TRUE){ 
  # define selection criterion for events
  select_in_subpopulationsEVENTS <- vector(mode="list")
  for (subpop in subpopulations[[thisdatasource]]){
    select <- "!is.na(person_id) "
    for (meaningevent in exclude_meaning_renamed[[thisdatasource]][[subpop]]){
      select <- paste0(select," & meaning_renamed!= '",meaningevent,"'")
    }
    select_in_subpopulationsEVENTS[[subpop]] <- select
  }
    # define selection criterion for Survey_OBSERVATIONS
  select_in_subpopulationsSO <- vector(mode="list")
  for (subpop in subpopulations[[thisdatasource]]){
      select <- "(!is.na(person_id) "
      for (itemsetSO in exclude_itemset_of_so[[thisdatasource]][[subpop]]){
        select <- paste0(select," & so_source_table != '",itemsetSO[1],"' & so_source_column != '",itemsetSO[2],"'")
      }
      select <- paste0(select,")")
      select_in_subpopulationsSO[[subpop]] <- select
  }
  
  # create multiple directories for export
  direxpsubpop <- vector(mode="list")
  for (subpop in subpopulations[[thisdatasource]]){
    direxpsubpop[[subpop]] <- paste0(thisdir,"/g_export_", subpop,'/')
    suppressWarnings(if (!file.exists(direxpsubpop[[subpop]])) dir.create(file.path(direxpsubpop[[subpop]])))
    file.copy(paste0(dirinput,'/METADATA.csv'), direxpsubpop[[subpop]], overwrite = T)
    file.copy(paste0(dirinput,'/CDM_SOURCE.csv'), direxpsubpop[[subpop]], overwrite = T)
    file.copy(paste0(dirinput,'/INSTANCE.csv'), direxpsubpop[[subpop]], overwrite = T)
    
    file.copy(paste0(thisdir,'/to_run.R'), direxpsubpop[[subpop]], overwrite = T)
  }
}

if (this_datasource_has_subpopulations==F) {
  dummytables <- paste0(direxp, "Dummy tables/")
  suppressWarnings(if (!file.exists(dummytables)) dir.create(file.path(dummytables)))
}


suffix <- vector(mode="list")

if (this_datasource_has_subpopulations == FALSE) {
  subpopulations_non_empty <- c('ALL')
  #subpopulations[[thisdatasource]] <- c('ALL')
  suffix[['ALL']] <- ''
  direxpsubpop <- vector(mode="list")
  direxpsubpop[['ALL']] <- paste0(thisdir, "/g_export/")
  
  dirtablesubpop <- vector(mode="list")
  dirtablesubpop[['ALL']] <- paste0(direxpsubpop[['ALL']], "Dummy tables/")
  suppressWarnings(if (!file.exists(dirtablesubpop[['ALL']])) dir.create(file.path(paste0(dirtablesubpop[['ALL']]))))
  
  dirD4D5subpop <- vector(mode="list")
  dirD4D5subpop[['ALL']] <- paste0(direxpsubpop[['ALL']], "D4-D5 tables/")
  suppressWarnings(if (!file.exists(dirD4D5subpop[['ALL']])) dir.create(file.path(paste0(dirD4D5subpop[['ALL']]))))
  
  dircomponents <- vector(mode="list")
  dircomponents[['ALL']] <- paste0(direxpsubpop[['ALL']], "components/")
  suppressWarnings(if (!file.exists(dircomponents[['ALL']])) dir.create(file.path(paste0(dircomponents[['ALL']]))))
}else{
  subpopulations_non_empty <- subpopulations[[thisdatasource]]
  dirtablesubpop <- vector(mode="list")
  dirD4D5subpop <- vector(mode="list")
  for (subpop in subpopulations_non_empty) {
    suffix[[subpop]] <- paste0('_', subpop)
    dirtablesubpop[[subpop]] <- paste0(direxpsubpop[[subpop]], "Dummy tables/")
    suppressWarnings(if(!file.exists(dirtablesubpop[[subpop]])) dir.create(file.path(paste0(dirtablesubpop[[subpop]]))))
    
    dirD4D5subpop[[subpop]] <- paste0(direxpsubpop[[subpop]], "D4-D5 tables/")
    suppressWarnings(if (!file.exists(dirD4D5subpop[[subpop]])) dir.create(file.path(paste0(dirD4D5subpop[[subpop]]))))
    
    dircomponents[[subpop]] <- paste0(direxpsubpop[[subpop]], "components/")
    suppressWarnings(if (!file.exists(dircomponents[[subpop]])) dir.create(file.path(paste0(dircomponents[[subpop]]))))
  }
}

for (subpop in subpopulations[[thisdatasource]]) {
  fileConn <- file(paste0(direxpsubpop[[subpop]], "subpop.txt"))
  writeLines(subpop, fileConn)
  close(fileConn)
}

rm(subpopulations)


