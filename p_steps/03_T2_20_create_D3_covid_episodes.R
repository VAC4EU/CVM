# CREATE COVID EPISODES REPEATED
#-----------------------------------------------
# input: D4_study_population, D3_events_COVID_narrow, covid_registry, COVID_test
# output: D3_covid_episodes, D3_covid_episodes_description and QC_covid_episodes

# in this step the list of unique episodes of covid is created, each person may have multiple episodes, separated by at least 60 days. The sources for the episodes are data source-specific.

print("CREATE COVID EPISODES REPEATED")

# data sources having registry
datasources_covid_registry <- c("TEST","ARS","BIFAP","CASERTA")

# data sources having positive tests
datasources_positive_tests <- c("TEST","SIDIAP","PEDIANET","UOSL", "FISABIO")

# data sources including all records with a covid diagnosis
datasources_covid_diagnosis_all <- c("FISABIO","SIDIAP", "UOSL", "CPRD", "PHARMO")

# data sources including only records of covid diagnosis from hospitals
datasources_covid_diagnosis_only_hosp <- c("TEST","ARS","CASERTA")

# load the subpopulation-independent input datasets
load(paste0(dirconceptsets,"I_COVID19DX_AESI_narrow.RData"))
# load(paste0(dirconceptsets,"I_COVID19DX_COV_narrow.RData"))
# # according to PHARMO, in their data source covid is ICPC code R83.03, which is stored in I_COVID19DX_COV_narrow and not in I_COVID19DX_AESI_narrow, therefore this is assigned as follows
# if (thisdatasource == 'PHARMO'){ 
#   dia_COVID_narrow <- I_COVID19DX_COV_narrow}else{
     dia_COVID_narrow <- I_COVID19DX_AESI_narrow
#     }
rm(I_COVID19DX_AESI_narrow)
#rm(I_COVID19DX_COV_narrow)

# OVERALL STRATEGY 
# 1 rbind all files that imply covid
# 2 remove recods closer than 60 days to a previous record
# 3 filter only episodes occurring during the study period and save
# 4 create a descriptive of all the components that contributed to identify an episode (from the date of the episode until the date of the next episode)

#---------------------------------
#---------------------------------
# 1 create list_all_covid_notifications
#---------------------------------
#---------------------------------

list_all_covid_notifications <- data.table()

#-------------------------------------
# add to list_all_covid_notifications positive results from covid test (after data source-specific processing)

if (thisdatasource %in% datasources_positive_tests){
  load(paste0(diritemsets,"COVID_test.RData")) 
  covid_test_positive <- COVID_test[,.(person_id,date, mo_source_value)][mo_source_value == "positive"]
  rm(COVID_test)
  covid_test_positive <- covid_test_positive[,.(person_id,date)]
  covid_test_positive <- covid_test_positive[,origin_case := "covid_test_positive" ]
  list_all_covid_notifications <- rbind(list_all_covid_notifications,covid_test_positive,fill = T)
}

#-------------------------------------
# add to list_all_covid_notifications records from covid registry

if (thisdatasource %in% datasources_covid_registry){
  load(paste0(dirpromptsets,"covid_registry.RData")) 
  covid_registry <- covid_registry[,.(person_id,date)]
  covid_registry <- covid_registry[,origin_case := "covid_registry" ]
  list_all_covid_notifications <- rbind(list_all_covid_notifications,covid_registry,fill = T)
}

#-------------------------------------
# add to list_all_covid_notifications diagnoses of covid (only in selected data sources, and sometimes only in hospital)

if (thisdatasource %in% datasources_covid_diagnosis_all) {
  dia_COVID_narrow <- dia_COVID_narrow[,origin_case := paste0("covid_narrow_m_",meaning_renamed)]
  dia_COVID_narrow <- dia_COVID_narrow[,.(person_id,date,origin_case)]
  list_all_covid_notifications <- rbind(list_all_covid_notifications,dia_COVID_narrow,fill = TRUE)
}
if (thisdatasource %in% datasources_covid_diagnosis_only_hosp) {
  dia_COVID_narrow <- dia_COVID_narrow[eval(parse(text = condmeaning[["HOSP"]])),]
  dia_COVID_narrow <- dia_COVID_narrow[,origin_case := paste0("covid_narrow_hosp_m_",meaning_renamed)]
  dia_COVID_narrow <- dia_COVID_narrow[,.(person_id,date,origin_case)]
  list_all_covid_notifications <- rbind(list_all_covid_notifications,dia_COVID_narrow,fill = TRUE)
}

#-------------------------------------
# CREATE SUBPOPULATION-SPECIFIC LIST OF EPISODES

for (subpop in subpopulations_non_empty) {
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData"))
  
  study_pop <- as.data.table(get(paste0("D4_study_population",suffix[[subpop]])))[,.(person_id,study_entry_date,study_exit_date)]
  persons_in_pop <- as.data.table(get(paste0("D4_study_population",suffix[[subpop]])))[,.(person_id)]

  list_all_covid_notificationssubpop <- list_all_covid_notifications
  
  #---------------------------------
  #---------------------------------
  # 2 remove notifications closer than 60 days to a previous notification
  #---------------------------------
  #---------------------------------
  
  # DT <- DT[, date_next_record := shift(date, n = 1, fill = NA, type=c("lead")), by = "pers_id"] 
  setorder(list_all_covid_notificationssubpop,"person_id","date")
  # list_all_covid_notificationssubpop <- list_all_covid_notificationssubpop[,date_next_record := shift(date, n = 1, fill = NA, type=c("lead")), by = "person_id"]
  # clean to keep only persons who are in fact in the study population (to avoid using too much memory during the next merge)
  list_all_covid_notificationssubpop <- merge(list_all_covid_notificationssubpop,persons_in_pop, all = FALSE, by = "person_id") 
  
  list_all_covid_notificationssubpop <- list_all_covid_notificationssubpop[,n:=seq_along(.I), by = "person_id"]
  
  copy <- copy(list_all_covid_notificationssubpop)
  copy <- copy[,.(person_id,date,n)]
  setnames(copy,c("date","n"),c("date_previous","n_previous"))
  listtodrop <- merge(list_all_covid_notificationssubpop,copy, all.x = TRUE, by = "person_id",allow.cartesian = TRUE)[date < date_previous + 60 & date > date_previous, ]
  listtodrop <- unique(listtodrop[,.(person_id,n)])
  listtodrop <- listtodrop[,todrop:= 1]
  listunique <- merge(list_all_covid_notificationssubpop,listtodrop, all.x = TRUE, by = c("person_id","n"))
  listunique <- listunique[ is.na(todrop), ]
  listunique <- unique(listunique[,.(person_id,date) ])
  
  
  #---------------------------------
  #---------------------------------
  # 3 filter episodes occurred during the study period
  #---------------------------------
  #---------------------------------
  listunique <- merge(study_pop,listunique, by = c("person_id"))[date >= study_entry_date & date <= study_exit_date, ]
  listunique <- unique(listunique[,.(person_id,date) ])
  
  tempname <- paste0("D3_covid_episodes",suffix[[subpop]])
  assign(tempname,listunique)
  save(list = tempname, file = paste0(dirtemp,tempname,".RData"))

  #---------------------------------
  #---------------------------------
  # 4 describe how each episode was detected
  #---------------------------------
  #---------------------------------
  
  setorder(listunique,"person_id","date")
  listdescr <- listunique[,n:=seq_along(.I), by = "person_id"]
  listdescr <- listdescr[,date_next_record := shift(date, n = 1, fill = NA, type = c("lead")), by = "person_id"]
  list_all_covid_notificationssubpop <- list_all_covid_notificationssubpop[,.(person_id,date,origin_case)]
  setnames(list_all_covid_notificationssubpop,c("date"),c("date_descr"))
  listdescr <- merge(listdescr,list_all_covid_notificationssubpop, all.x = TRUE, by = "person_id",allow.cartesian = TRUE)[date_descr >= date & (date_descr < date_next_record | is.na(date_next_record)), ]
  listdescr <- unique(listdescr[,.(person_id,date,n,origin_case)])
  listdescr <- listdescr[, component := 1]
  listdescr <- dcast(listdescr,person_id + date + n  ~ origin_case, value.var = "component", fill = 0 )
  listdescr <- as.data.table(listdescr)[, year := year(date)]
  
  tempname <- paste0("D3_covid_episodes_description",suffix[[subpop]])
  assign(tempname,listdescr)
  save(list = tempname, file = paste0(dirtemp,tempname,".RData"))
  
  columns_listdescr <- colnames(listdescr)[colnames(listdescr) %not in% c("person_id","date")]
  listdescr <- listdescr[, .N, by = columns_listdescr]

  tempname <- paste0("QC_covid_episodes",suffix[[subpop]])
  assign(tempname,listdescr)
  save(list = tempname, file = paste0(diroutput,tempname,".RData"))
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  fwrite(listdescr, file = paste0(thisdirexp,tempname,".csv"))
  
  
  # clean memory
  
  rm(list_all_covid_notificationssubpop,copy,listunique,listtodrop,persons_in_pop,columns_listdescr,listdescr)
  rm(list = paste0("D4_study_population", suffix[[subpop]]))
  rm(list = tempname)
}

rm(list_all_covid_notifications)
rm(dia_COVID_narrow)






