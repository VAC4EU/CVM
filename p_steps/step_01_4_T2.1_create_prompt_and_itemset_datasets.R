# -----------------------------------------------------------------------
# RETRIEVE PROMPT DATASETS and ITEMSET DATASETS

# input: SURVEY_ID, SURVEY_OBSERVATIONS
# output: prompt datasets and itemset datasets

print("RETRIEVE RECORDS FROM SURVEY")

# RETRIEVE FROM SURVEY_ID ALL prompt datasets corresponding to "covid_registry" 

# collect and rbind from all files whose name starts with 'SURVEY_ID'
SURVEY_ID_COVID <- data.table()
files<-sub('\\.csv$', '', list.files(dirinput))
for (i in 1:length(files)) {
  if (str_detect(files[i],"^SURVEY_ID")) {
    SURVEY_ID_COVID <-rbind(SURVEY_ID_COVID,fread(paste0(dirinput,files[i],".csv"), colClasses = list( character="person_id"))[survey_meaning =="covid_registry",])
  }
}

covid_registry <- SURVEY_ID_COVID[,date:=ymd(survey_date)]
covid_registry <- covid_registry[,-"survey_date"]


# RETRIEVE FROM SURVEY_OBSERVATIONS ALL itemset datasets corresponding to "COVID_symptoms" 
#-----------------------------------------------------


CreateItemsetDatasets(EAVtables = ConcePTION_CDM_EAV_tables_retrieve,
                      datevar= ConcePTION_CDM_datevar_retrieve,
                      dateformat= "YYYYmmdd",
                      rename_col = list(person_id=person_id,date=date),
                      study_variable_names = study_variables_of_our_study,
                      itemset = itemset_AVpair_our_study_this_datasource,
                      dirinput = dirinput,
                      diroutput = dirtemp,
                      extension = c("csv"))
# SURVEY_OBS_COVIDSYMPT <- data.table()
# files<-sub('\\.csv$', '', list.files(dirinput))
# for (i in 1:length(files)) {
#   if (str_detect(files[i],"^SURVEY_OBS")) {
#     SURVEY_OBS_COVIDSYMPT <-rbind(SURVEY_OBS_COVIDSYMPT,fread(paste0(dirinput,files[i],".csv"))[so_source_column =="STATOCLINICO_PIU_GRAVE",])
#   }
# }
# 
# COVID_symptoms <- SURVEY_OBS_COVIDSYMPT[,date:=ymd(so_date)]
# COVID_symptoms <- COVID_symptoms[,-"so_date"]

save(covid_registry,file = paste0(dirtemp,"covid_registry.RData"))
#save(COVID_symptoms,file = paste0(dirtemp,"COVID_symptoms.RData"))

#rm(SURVEY_ID_COVID, covid_registry, SURVEY_OBS_COVIDSYMPT,COVID_symptoms)
rm(SURVEY_ID_COVID, covid_registry, files)
rm(list=names(itemset_AVpair_our_study_this_datasource))
