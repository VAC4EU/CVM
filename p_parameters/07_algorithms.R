#------------------------------------------------------------------------
# associate with each variable their algorithm
# also create lists of variables with different types of algorithms

# input: partly from VAC4EU format

# OUTPUT

# variable_definition: is a list of lists, arguments in unique(c(OUTCOME_variables,CONTROL_variables,COV_variables_raw,concept_sets_of_our_study)); values:   variable_definition[[varname]] is the list of conceptsets used to compute varname

# elements_for_TD_variables: it a list of list, arguments in the list of TD variables DIAB CANCER PULMON OBES CKD HIV IMMUNOSUP SICKLE CVD, assigns to each TD variables their corresponding conceptsets

# SECCOMPONENTS: list of variable names that are computed using a complex algorithm (secondary components)

# meanings_of_this_study: this is a list of lists, that groups meanings thatshould be used for homogeneous purposes across data sources, names are HOSP (meanings that mean that the diagnostic code is observed during a hospitalization); PC (meanings that mean that the diagnostic code is observed during a primary care contact); HOSPsec (meanings that mean that the diagnostic code is observed during a hospitalization as a secondary diagnosis)

# select_meanings_AESI: this is the selection rule for some data sources to apply to AESI (e.g., discard secondary diagnoses)

# parameters for COVID severity

# parameters for secondary components (to be described): ...

#----------------------------------------------------
# variable_definition
#----------------------------------------------------


#----------------------------------------------------
# Diagnosis
NoAlgo <- VAR_codelist[!(Algorithm & !Algorithm_input), ]
OUTCOME_concepts_in_var <- sapply(NoAlgo[(AESI), Varname], paste0, c("_narrow"))
OUTCOME_concepts_in_comp <- lapply(Map(identity, NoAlgo[(AESI), Varname]), paste0, c("_narrow", "_possible"))
COV_concepts_in_var <- lapply(Map(identity, NoAlgo[(COV), Varname]), paste0, c("_narrow", "_possible"))
rm(VAR_codelist)

NEG_concepts_in_var <- sapply(NoAlgo[(NEG), Varname], paste0, c("_narrow"))
NEG_concepts_in_comp <- lapply(Map(identity, NoAlgo[(NEG), Varname]), paste0, c("_narrow", "_possible"))

# Divide the definition of algorithm input which are neither AESI, NEG or COV
NoAlgo <- NoAlgo[(Algorithm_input) & !((AESI) | (NEG) | (COV))]
NoAlgo <- NoAlgo[, c("system", "concept", "type") := tstrsplit(Varname, "_", fixed = T)][, .(Varname, type)]
OUTCOME_algo_input <- sapply(NoAlgo[type == "AESI", Varname], paste0, c("_narrow"))
OUTCOME_algo_input_possible <- lapply(Map(identity, NoAlgo[type == "AESI", Varname]), paste0, c("_narrow", "_possible"))
COV_algo_input <- lapply(Map(identity, NoAlgo[type == "COV", Varname]), paste0, c("_narrow", "_possible"))

OUTCOME_concepts_in_var <- c(OUTCOME_concepts_in_var, OUTCOME_algo_input)
OUTCOME_concepts_in_comp <- c(OUTCOME_concepts_in_comp, OUTCOME_algo_input_possible)
COV_concepts_in_var <- c(COV_concepts_in_var, COV_algo_input)
rm(COV_algo_input, OUTCOME_algo_input, OUTCOME_algo_input_possible)

#Drugs
NoAlgo <- DRUG_codelist[!(Algorithm), ]
DRUG_concepts_in_var <- sapply(NoAlgo[(AESI) | (COV) | (Algorithm_input), Varname], identity)
rm(DRUG_codelist, NoAlgo)

variable_definition <- c(OUTCOME_concepts_in_var, COV_concepts_in_var, DRUG_concepts_in_var, NEG_concepts_in_var)
component_definition <- c(OUTCOME_concepts_in_comp, NEG_concepts_in_comp)
outcome_definition <- OUTCOME_concepts_in_comp
not_outcome_definition <- c(COV_concepts_in_var, DRUG_concepts_in_var, NEG_concepts_in_var)
complete_definition <- c(OUTCOME_concepts_in_comp, COV_concepts_in_var, DRUG_concepts_in_var, NEG_concepts_in_var)
rm(OUTCOME_concepts_in_var, COV_concepts_in_var, DRUG_concepts_in_var, OUTCOME_concepts_in_comp)

for (var in names(variable_definition)) {
  if (var %in% CONCEPTSETS_to_be_split) {
    variable_definition[[var]] <- paste(var, seq_len(numbers_split[[which(CONCEPTSETS_to_be_split %in% var)]]), sep = "_")
  }
}

for (var in names(component_definition)) {
  if (var %in% CONCEPTSETS_to_be_split) {
    component_definition[[var]] <- paste(var, seq_len(numbers_split[[which(CONCEPTSETS_to_be_split %in% var)]]), sep = "_")
  }
}

ALGO_codelist <- readxl::read_excel(File_variables_ALG_DP_ROC20, sheet = "ALG")
ALGO_codelist <- as.data.table(ALGO_codelist)
rm(File_variables_ALG_DP_ROC20)

SECCOMPONENTS <- c("B_TTS_AESI","E_DM1_AESI")

# algortihms for outcomes

ALGO_link <- ALGO_codelist[Algorithm %in% unique(ALGO_codelist[Algorithm %not in% SECCOMPONENTS, Algorithm]),
                           .(Algorithm, VariableName)]

ALGO_link_components <- copy(ALGO_link)[, test := complete_definition[VariableName]][, VariableName := NULL]
ALGO_link <- ALGO_link[, test := variable_definition[VariableName]][, VariableName := NULL]

ALGO_link <- split(ALGO_link, by = "Algorithm", keep.by = F)
ALGO_link <- lapply(ALGO_link, unlist, use.names = F)

variable_definition <- variable_definition[names(variable_definition) %in% intersect(names(variable_definition), names(ALGO_link)) == FALSE]

variable_definition <- c(variable_definition, ALGO_link)
rm(ALGO_link)

ALGO_link_components <- split(ALGO_link_components, by = "Algorithm", keep.by = F)
ALGO_link_components <- lapply(ALGO_link_components, unlist, use.names = F)

component_definition <- component_definition[names(component_definition) %in% intersect(names(component_definition), names(ALGO_link_components)) == FALSE]

component_definition <- c(component_definition, ALGO_link_components)
rm(ALGO_link_components, complete_definition)

# variables_with_split_conceptsets <- c()
# for (var in names(variable_definition)) {
#   if (any(variable_definition[[var]] %in% CONCEPTSETS_to_be_split)) {
#     variables_with_split_conceptsets <- append(variables_with_split_conceptsets, var)
#   }
# }


# elements_for_TD_variables DIAB CANCER PULMON OBES CKD HIV IMMUNOSUP SICKLE CVD

list_of_covariates_for_cohort <- c("DIAB", "CANCER", "PULMON", "OBES", "CKD", "HIV", "IMMUNOSUP", "SICKLE", "CVD")
list_of_comedication_for_cohort <- c("DP_ANTIBIO", "DP_ANTITHROMBOTIC", "DP_ANTIVIR", "DP_SEXHORMONES", "DP_VACCINES",
                                     "DP_IMMUNOSUPPR", "DP_LIPIDLOWER")
list_of_prior_history_for_cohort <- c("V_PRIORVTEALGORITHM_COV", "Im_ANAPHYLAXIS_COV", "Im_ALLERGY_COV")

elements_for_TD_variables <- vector(mode = "list")
elements_for_TD_variables[['CVD']] <- "C_CARDIOCEREBROVASCULARDESE_COV"
elements_for_TD_variables[['CANCER']] <- "Onc_ANYMALIGNANCYALGORITHM_COV"
elements_for_TD_variables[['PULMON']] <- "R_RESPCHRONICALGORITHM_COV"
elements_for_TD_variables[['HIV']] <- c("I_HIV_COV","I_AIDS_CH")
elements_for_TD_variables[['CKD']] <- c("D_LIVERCHRONIC_COV","G_KDCHRONIC_COV")
elements_for_TD_variables[['DIAB']] <- c("E_DM12ALGORITHM_COV")
elements_for_TD_variables[['OBES']] <- c("L_OBESITYALGORITHM_COV")
elements_for_TD_variables[['SICKLE']] <- c("B_SICKLECELLALGORITHM_COV")
elements_for_TD_variables[['IMMUNOSUP']] <- c("Im_INMUNODEFALLALGORITHM_COV")
elements_for_TD_variables[['DP_ANTIBIO']] <- c("DP_ANTIBIO")
elements_for_TD_variables[['DP_ANTITHROMBOTIC']] <- c("DP_ANTITHROMBOTIC")
elements_for_TD_variables[['DP_ANTIVIR']] <- c("DP_ANTIVIR")
elements_for_TD_variables[['DP_SEXHORMONES']] <- c("DP_SEXHORMONES")
elements_for_TD_variables[['DP_VACCINES']] <- c("DP_VACCINES")
elements_for_TD_variables[['DP_IMMUNOSUPPR']] <- c("DP_IMMUNOSUPPR")
elements_for_TD_variables[['DP_LIPIDLOWER']] <- c("DP_LIPIDLOWER")
elements_for_TD_variables[['V_PRIORVTEALGORITHM_COV']] <- c("V_PRIORVTEALGORITHM_COV")
elements_for_TD_variables[['Im_ANAPHYLAXIS_COV']] <- c("Im_ANAPHYLAXIS_COV")
elements_for_TD_variables[['Im_ALLERGY_COV']] <- c("Im_ALLERGY_COV")

#----------------------------------------------------
# meanings_of_this_study
#----------------------------------------------------


# we need to create two groups of meanings: one referring to hospitals HOSP (excluding emergency care) and one referring to primary care PC

meanings_of_this_study <- vector(mode="list")
meanings_of_this_study[["HOSP"]]=c("hospitalisation_primary","hospitalisation_secondary","hospitalisation_secundary","hospital_diagnosis","hopitalisation_diagnosis_unspecified","episode_primary_diagnosis","episode_secondary_diagnosis","diagnosis_procedure","hospitalisation_associated","hospitalisation_linked","HH","NH","hospitalisation_ICU_primary","hospitalisation_ICU_secondary","hospitalisation_ICU_unspecified")
meanings_of_this_study[["PC"]]=c("primary_care_event","primary_care_diagnosis","primary_care_events_BIFAP","primary_care_antecedents_BIFAP","primary_care_condicionants_BIFAP")
meanings_of_this_study[["HOSPsec"]]=c("hospitalisation_secondary","hospitalisation_secundary")
meanings_of_this_study[["assessment_chronic"]] = c("exemption","evaluation_during_home_or_residential_care")

# create conditions on the meaning_renamed variable, associated to the lists above ("PC", "HOSP", ...)

condmeaning <- list()
for (level1 in c("HOSP","PC","HOSPsec","assessment_chronic")) {
  for (meaning in meanings_of_this_study[[level1]]) {
    if (length(condmeaning[[level1]])==0) {condmeaning[[level1]]=paste0("meaning_renamed=='",meanings_of_this_study[[level1]][[1]],"'")
    }else{
      condmeaning[[level1]]=paste0(condmeaning[[level1]], " | meaning_renamed=='",meaning,"'")
    }
  }
}
rm(meanings_of_this_study, level1, meaning)

#----------------------------------------------------
# select_meanings_AESI
#----------------------------------------------------


#--------------------------
# create datasource-specific algorithms for AESIs: the conditions below are used (recersed) to discard from the conceptset datasets records with the specified meanings

select_meanings_AESI <- list()
select_meanings_AESI[[thisdatasource]] <- "(is.na(person_id))" 
select_meanings_AESI[["BIFAP"]] <- condmeaning[["HOSPsec"]]
select_meanings_AESI[["ARS"]] <- condmeaning[["assessment_chronic"]]
select_meanings_AESI[["TEST"]] <- condmeaning[["HOSPsec"]]

#----------------------------
# VARIABLES FOR COVID SEVERITY
#----------------------------

# INCLUSION

# data sources having registry
datasources_covid_registry <- c("TEST","ARS","BIFAP","CASERTA")

# data sources having positive tests
datasources_positive_tests <- c("TEST","SIDIAP","PEDIANET","UOSL", "FISABIO")

# data sources including all records with a covid diagnosis
datasources_covid_diagnosis_all <- c("FISABIO","SIDIAP", "UOSL", "CPRD", "PHARMO")

# data sources including only records of covid diagnosis from hospitals
datasources_covid_diagnosis_only_hosp <- c("TEST","ARS","CASERTA")


# HOSPITALISATIONS

# data sources including records of hospitalisations with a covid diagnosis
datasources_hosp_due_to_covid <- c("TEST","ARS","CASERTA","FISABIO","SIDIAP","UOSL")

# data sources including records of hospitalizations irrespective of diagnosis
datasources_hosp_after_covid <- c("TEST","PEDIANET")

# data sources including hospitalization from covid registry (to be handled in a data source-tailored manner below)
datasources_hosp_from_covid_registry <- c("TEST","ARS","BIFAP","CASERTA")

# ICU

# data sources including records of ARDS diagnosis (D3_ARDS_narrow and D3_ARDS possible)
datasources_ICU_from_ARDS <- c("TEST","ARS","CASERTA","FISABIO","SIDIAP","UOSL","CPRD")
# R_ARDS_AESI_narrow

# data sources including records of procedures of mechanical ventilation (conceptset ICU_VENTILATION)
datasources_proc_mechanical_ventilation <- c("TEST","ARS","CASERTA")

# data sources including records of access to ICU from hospitalisation (meaning "hospitalisation_ICU_unspecified" in D3_covid_narrow)
datasources_access_ICU <- c("TEST","FISABIO")

# data sources including records of access to ICU from free text in MEDICAL_OBSERVATIONS (itemset extracted_from_free_text)
datasources_ICU_free_text <- c("TEST","PEDIANET")


# data sources including access to ICU from covid registry (from various itemsets, to be handled in a data source-tailored manner below)
datasources_ICU_from_covid_registry <- c("TEST","ARS","BIFAP","CASERTA")

# DEATH 

# data sources using death within 56 days
datasources_death_after_covid <- c("TEST","ARS","CASERTA","FISABIO","SIDIAP","BIFAP","PEDIANET","UOSL","CPRD")

# data sources including access to death from covid registry (to be handled in a data source-tailored manner below)
datasources_death_from_covid_registry <- c("TEST","ARS","CASERTA")



# 
#----------------------------
# SECONDARY COMPONENTS
#----------------------------

# SECCOMPONENTS <- c("ArterialNoTP", "ArterialTP", "VTENoTP", "VTETP", "ArterialVTENoTP", "ArterialVTETP", "CVSTNoTP", "CVSTTP")

concept_set_seccomp <- vector(mode="list")
rule_seccomp <- vector(mode="list")
distance_seccomp <- vector(mode="list")
direction_seccomp <- vector(mode="list")
upper_distance_seccomp <- vector(mode="list")
lower_distance_seccomp <- vector(mode="list")
selectionrule_direction_seccomp <- vector(mode="list")


for (SECCOMP in SECCOMPONENTS) {
  test <- ALGO_codelist[Algorithm == SECCOMP]
  
  upper <- as.numeric(test[, unique(upper)])
  lower <- as.numeric(test[, unique(lower)])
  
  if (SECCOMP == "E_DM1_AESI") {
    upper <- 99999
    lower <- 99999
  }
  
  upper_distance_seccomp[[SECCOMP]] = upper
  lower_distance_seccomp[[SECCOMP]] = abs(lower)
  
 # if (upper == abs(lower)) {
  #  distance_seccomp[[SECCOMP]] = upper
    direction_seccomp[[SECCOMP]] = "Either direction"
  # }
  
  # distance_seccomp[[SECCOMP]] = '10'
  # direction_seccomp[[SECCOMP]] = "Either direction"
  
  selectionrule_direction_seccomp[[SECCOMP]]["A before B"] <- paste0("dateA <= dateB  & dateB <= dateA + ", upper_distance_seccomp[[SECCOMP]])
  selectionrule_direction_seccomp[[SECCOMP]]["B before A"] <- paste0("dateB <= dateA  & dateA <= dateB + ", lower_distance_seccomp[[SECCOMP]])
  selectionrule_direction_seccomp[[SECCOMP]]["Either direction"] <- paste0('((',selectionrule_direction_seccomp[[SECCOMP]]["A before B"],') | (',selectionrule_direction_seccomp[[SECCOMP]]["B before A"],'))')
  
}
rm(distance_seccomp, SECCOMP, upper, lower)

for (SECCOMP in SECCOMPONENTS) {
  test <- ALGO_codelist[Algorithm == SECCOMP]
  for (i in seq_len(test[, max(group)])) {
    for (var in unlist(test[group == i, .(VariableName)])) {
      concept_set_seccomp[[SECCOMP]][[LETTERS[i]]] <- c(concept_set_seccomp[[SECCOMP]][[LETTERS[i]]],
                                                        variable_definition[[var]])
      component_definition[[SECCOMP]] <- c(component_definition[[SECCOMP]], variable_definition[[var]])
    }
  }
  
  rule_seccomp[[SECCOMP]] <- test[, unique(between)]
  
}

rm(ALGO_codelist, test, i, SECCOMP, var)

# concept sets specific for datasources

# if (thisdatasource == 'ARS'){
#   # TODO check which definition is correct
#   #concept_set_codes_our_study_pre[["COVID_narrow"]][["ICD9"]] <- c(concept_set_codes_our_study_pre[["COVID_narrow"]][["ICD9"]],'043','48041','51891','51971')
#   concept_set_codes_our_study_pre[["ARD_narrow"]][["ICD9"]] <- c(concept_set_codes_our_study_pre[["ARD_narrow"]][["ICD9"]],'5189')
# 
#   concept_set_codes_our_study_pre[["COVID_narrow"]][["ICD9"]] <- c(concept_set_codes_our_study_pre[["COVID_narrow"]][["ICD9"]],'043','48041','51891','51971')
#   concept_set_codes_our_study_pre[["ARDS_narrow"]][["ICD9"]] <- c(concept_set_codes_our_study_pre[["ARDS_narrow"]][["ICD9"]],'5189')
# 
# }

#-------------------------------------
# set concept sets

# # augment ICPC codes
# for (outcome in OUTCOME_variables){
#   outnarrow <- paste0(outcome,'_narrow')
#   outpossible <- paste0(outcome,'_possible')
#   if (length(concept_set_codes_our_study[[outnarrow]][["ICPC"]]) == 0 & length(concept_set_codes_our_study[[outnarrow]][["ICPC2P"]]) > 0){
#     concept_set_codes_our_study[[outpossible]][["ICPC"]] <- unique(c(concept_set_codes_our_study[[outpossible]][["ICPC"]],substr(concept_set_codes_our_study[[outnarrow]][["ICPC2P"]],1,3)))
#   }
# }
# rm(outcome, outnarrow, outpossible)

# for (conceptset in c(COV_variables, DP_variables, VACCINES_variable)){
#   if (length(concept_set_codes_our_study[[conceptset]][["ICPC2P"]]) > 0){
#     concept_set_codes_our_study[[conceptset]][["ICPC"]] <- unique(c(concept_set_codes_our_study[[conceptset]][["ICPC"]],substr(concept_set_codes_our_study[[conceptset]][["ICPC2P"]],1,3)))
#   }
# }

#-------------------------------------
# fix for ICD10GM

for (conceptset in concept_sets_of_our_study){
  if (concept_set_domains[[conceptset]] == "Diagnosis"){
    concept_set_codes_our_study[[conceptset]][["ICD10GM"]] <- concept_set_codes_our_study[[conceptset]][["ICD10CM"]]
  }
}

#-------------------------------------
# fix for ICD10CM
for (conceptset in concept_sets_of_our_study){
  if (concept_set_domains[[conceptset]] == "Diagnosis"){
    concept_set_codes_our_study[[conceptset]][["ICD10"]] <- concept_set_codes_our_study[[conceptset]][["ICD10CM"]]
  }
}
#-------------------------------------
# fix for ICD9CM
for (conceptset in concept_sets_of_our_study){
  if (concept_set_domains[[conceptset]] == "Diagnosis"){
    concept_set_codes_our_study[[conceptset]][["ICD9"]] <- concept_set_codes_our_study[[conceptset]][["ICD9CM"]]
  }
}

# ### Select only conceptset for SCRI
# variable_definition <- unlist(variable_definition[SCRI_list_variables])

# concept_set_codes_our_study <- concept_set_codes_our_study[names(concept_set_codes_our_study) %in% c(variable_definition, "COVID_VACCINES")]

# concept_set_codes_our_study <- concept_set_codes_our_study[names(concept_set_codes_our_study) %in% c(variable_definition, "COVID_VACCINES", "I_COVID19DX_AESI_narrow", "I_COVID19DX_COV_narrow")]


save(concept_set_codes_our_study,file=paste0(direxp,"concept_set_codes_our_study.RData"))

conceptsets_exact_matching <- intersect(concept_sets_of_our_study,
                                        c(unlist(outcome_definition, use.names = F),
                                          "Im_HYPERSENS_AESI_narrow", "Im_HYPERSENS_AESI_possible"))

conceptsets_children_matching <- setdiff(concept_sets_of_our_study,
                                         c(unlist(outcome_definition, use.names = F),
                                           "Im_HYPERSENS_AESI_narrow", "Im_HYPERSENS_AESI_possible"))

if (this_datasource_has_subpopulations == TRUE){
  for (subpop in subpopulations_non_empty){
    save(concept_set_codes_our_study, file = paste0(direxpsubpop[[subpop]], "concept_set_codes_our_study.RData"))
    
  }
  rm(subpop)
}

rm(conceptset, datasources_with_subpopulations)

