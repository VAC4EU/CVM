###################################################################
# DESCRIBE THE ATTRIBUTE-VALUE PAIRS
###################################################################
CDM_SOURCE<- fread(paste0(dirinput,"CDM_SOURCE.csv"))
thidatasource <- as.character(CDM_SOURCE[1,3])

# -itemset_AVpair_our_study- is a nested list, with 3 levels: foreach study variable, for each coding system of its data domain, the list of AVpair is recorded

study_variables_of_our_study <- c("COVID_symptoms","COVID_hospitalised","COVID_hospitalised_date","COVID_ICU","COVID_ICU_date","COVID_test")

itemset_AVpair_our_study <- vector(mode="list")
datasources<-c("TEST","ARS","BIPS","BIFAP","FISABIO","SIDIAP","PEDIANET","PHARMO")


# specification COVID_symptoms
files<-sub('\\.csv$', '', list.files(dirinput))
for (i in 1:length(files)) {
  if (str_detect(files[i],"^SURVEY_OB")) {
    itemset_AVpair_our_study[["COVID_symptoms"]][[files[i]]][["ARS"]] <- list(list("COVIDDATASET","STATOCLINICO_PIU_GRAVE"))
    itemset_AVpair_our_study[["COVID_symptoms"]][[files[i]]][["TEST"]] <- list(list("COVIDDATASET","STATOCLINICO_PIU_GRAVE"))
    itemset_AVpair_our_study[["COVID_hospitalised"]][[files[i]]][["TEST"]] <- list(list("Covid19_Hospitalizacion","Ingreso_hospital"))
    itemset_AVpair_our_study[["COVID_hospitalised_date"]][[files[i]]][["TEST"]] <- list(list("Covid19_Hospitalizacion","Fecha_ingreso_hosp"))
    itemset_AVpair_our_study[["COVID_ICU"]][[files[i]]][["TEST"]] <- list(list("Covid19_UCI","Ingreso_uci"))
    itemset_AVpair_our_study[["COVID_ICU_date"]][[files[i]]][["TEST"]] <- list(list("Covid19_UCI","Fecha_ingreso_uci"))
    itemset_AVpair_our_study[["COVID_hospitalised"]][[files[i]]][["BIFAP"]] <- list(list("Covid19_Hospitalizacion","Ingreso_hospital"))
    itemset_AVpair_our_study[["COVID_hospitalised_date"]][[files[i]]][["BIFAP"]] <- list(list("Covid19_Hospitalizacion","Fecha_ingreso_hosp"))
    itemset_AVpair_our_study[["COVID_ICU"]][[files[i]]][["BIFAP"]] <- list(list("Covid19_UCI","Ingreso_uci"))
    itemset_AVpair_our_study[["COVID_ICU_date"]][[files[i]]][["BIFAP"]] <- list(list("Covid19_UCI","Fecha_ingreso_uci"))
    
  }
  if (str_detect(files[i],"^MEDICAL_OB")) {
    itemset_AVpair_our_study[["COVID_test"]][[files[i]]][["TEST"]] <- list(list("COVID-Test","covid19_pcr_test"),list("COVID-Test","covid19_antigen_test"))
    itemset_AVpair_our_study[["COVID_test"]][[files[i]]][["SIDIAP"]] <- list(list("COVID-Test","covid19_pcr_test"),list("COVID-Test","covid19_antigen_test"))
    
  }  
}


itemset_AVpair_our_study_this_datasource<-vector(mode="list")

for (t in  names(itemset_AVpair_our_study)) {
  for (f in names(itemset_AVpair_our_study[[t]])) {
    for (s in names(itemset_AVpair_our_study[[t]][[f]])) {
      if (s==thidatasource ){
        itemset_AVpair_our_study_this_datasource[[t]][[f]]<-itemset_AVpair_our_study[[t]][[f]][[s]]
      }
    }
  }
}
    

# 
# for (t in  names(person_id)) {
#   if (t=="Diagnosis")  person_id = person_id [[t]]
# }
# 
# 
# for (t in  names(date)) {
#   if (t=="Diagnosis")  date = date [[t]]
# }

