###################################################################
# DESCRIBE THE ATTRIBUTE-VALUE PAIRS
###################################################################

# -study_variables_of_our_study- is the list of itemsets that will contribute to study variables
study_variables_of_our_study <- c("COVID_symptoms","COVID_hospitalised","COVID_hospitalised_date","COVID_ICU","COVID_ICU_date","COVID_test","extracted_from_free_text")


# -itemset_AVpair_our_study- is a nested list, with 3 levels: foreach itemset, for each coding system of its data domain, the list of pairs of attributes is recorded

itemset_AVpair_our_study <- vector(mode="list")
itemset_AVpair_our_study_meaning <- vector(mode="list")
datasources <- c("TEST","ARS","BIPS","BIFAP","FISABIO","SIDIAP","PEDIANET","PHARMO","UOSL")


# specification COVID_symptoms
for (file in files_ConcePTION_CDM_tables[["SURVEY_OBSERVATIONS"]]){ 
    itemset_AVpair_our_study[["COVID_symptoms"]][[file]][["ARS"]] <- list(list("COVIDDATASET","STATOCLINICO_PIU_GRAVE"))
    itemset_AVpair_our_study[["COVID_symptoms"]][[file]][["TEST"]] <- list(list("COVIDDATASET","STATOCLINICO_PIU_GRAVE"))
    itemset_AVpair_our_study[["COVID_symptoms"]][[file]][["CASERTA"]] <- list(list("TAB_COVID","SINTOMATOLOGIA"))
    itemset_AVpair_our_study[["COVID_hospitalised"]][[file]][["TEST"]] <- list(list("Covid19_Hospitalizacion","Ingreso_hospital"))
    itemset_AVpair_our_study[["COVID_hospitalised_date"]][[file]][["TEST"]] <- list(list("Covid19_Hospitalizacion","Fecha_ingreso_hosp"))
    itemset_AVpair_our_study[["COVID_ICU"]][[file]][["TEST"]] <- list(list("Covid19_UCI","Ingreso_uci"))
    itemset_AVpair_our_study[["COVID_ICU_date"]][[file]][["TEST"]] <- list(list("Covid19_UCI","Fecha_ingreso_uci"))
    itemset_AVpair_our_study[["COVID_hospitalised"]][[file]][["BIFAP"]] <- list(list("Covid19_Hospitalizacion","Ingreso_hospital"))
    itemset_AVpair_our_study[["COVID_hospitalised_date"]][[file]][["BIFAP"]] <- list(list("Covid19_Hospitalizacion","Fecha_ingreso_hosp"))
    itemset_AVpair_our_study[["COVID_ICU"]][[file]][["BIFAP"]] <- list(list("Covid19_UCI","Ingreso_uci"))
    itemset_AVpair_our_study[["COVID_ICU_date"]][[file]][["BIFAP"]] <- list(list("Covid19_UCI","Fecha_ingreso_uci"))
}
for (file in files_ConcePTION_CDM_tables[["MEDICAL_OBSERVATIONS"]]){ 
  # itemset_AVpair_our_study[["COVID_test"]][[file]][["TEST"]] <- list(list("SIDIAP.Covid_tests","PCR"),list("SIDIAP.Covid_tests","ANTIGENS"))
  itemset_AVpair_our_study[["COVID_test"]][[file]][["SIDIAP"]] <- list(list("SIDIAP.Covid_tests","PCR"),list("SIDIAP.Covid_tests","ANTIGENS"))
  itemset_AVpair_our_study[["extracted_from_free_text"]][[file]][["TEST"]] <- list(list("RICOVERI_FSE","PROCEDURE_free_text"))
  itemset_AVpair_our_study[["extracted_from_free_text"]][[file]][["PEDIANET"]] <- list(list("RICOVERI_FSE","PROCEDURE_free_text"))
  
  # TEST with mo_origin and mo_meaning
  itemset_AVpair_our_study_meaning[["COVID_test"]][[file]][["TEST"]] <- list(list("RedMIVA","covid19_pcr_test"),list("RedMIVA","covid19_antigen_test"))
  # fisabio uses mo_origin and mo_meaning
  itemset_AVpair_our_study_meaning[["COVID_test"]][[file]][["FISABIO"]] <- list(list("RedMIVA","covid19_pcr_test"),list("RedMIVA","covid19_antigen_test"))
  # PEDIANET uses mo_origin and mo_meaning
  itemset_AVpair_our_study_meaning[["COVID_test"]][[file]][["PEDIANET"]] <- list(list("TAMPONI_COVID19","path_test"))
  # UOSL uses mo_origin and mo_meaning
  itemset_AVpair_our_study_meaning[["COVID_test"]][[file]][["UOSL"]] <- list(list("MSIS","covid19_positive_test"))
}

itemset_AVpair_our_study_this_datasource <- vector(mode="list")
itemset_AVpair_our_study_this_datasource_meaning <- vector(mode="list")
study_variables_this_datasource <- vector(mode="list")
study_variables_this_datasource_meaning <- vector(mode="list")


for (t in study_variables_of_our_study){
  for (f in c(files_ConcePTION_CDM_tables[["SURVEY_OBSERVATIONS"]],files_ConcePTION_CDM_tables[["MEDICAL_OBSERVATIONS"]])) {
    if (length(itemset_AVpair_our_study[[t]][[f]][[thisdatasource]]) > 0){
    itemset_AVpair_our_study_this_datasource[[t]][[f]] <- itemset_AVpair_our_study[[t]][[f]][[thisdatasource]]
    study_variables_this_datasource <- unique(c(study_variables_this_datasource,t))
    }
    if (length(itemset_AVpair_our_study_meaning[[t]][[f]][[thisdatasource]]) > 0){
  itemset_AVpair_our_study_this_datasource_meaning[[t]][[f]] <- itemset_AVpair_our_study_meaning[[t]][[f]][[thisdatasource]]
  study_variables_this_datasource_meaning <- unique(c(study_variables_this_datasource_meaning,t))
  
    }
  }
}

rm(t, f, file)
