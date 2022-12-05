# associate with each variable their algorithm, partly from VAC4EU format

# also create lists of variables with different types of algorithms

# Diagnosis
NoAlgo <- VAR_codelist[!(Algorithm), ]
OUTCOME_concepts_in_var <- sapply(NoAlgo[(AESI) | (NEG), Varname], paste0, c("_narrow"))
COV_concepts_in_var <- lapply(Map(identity, NoAlgo[(COV) | (Algorithm_input), Varname]), paste0, c("_narrow", "_possible"))
rm(VAR_codelist)

# Divide the definition of algorithm input which are neither AESI, NEG or COV
NoAlgo <- NoAlgo[(Algorithm_input) & !((AESI) | (NEG) | (COV))]
NoAlgo <- NoAlgo[, c("system", "concept", "type") := tstrsplit(Varname, "_", fixed = T)][, .(Varname, type)]
OUTCOME_algo_input <- sapply(NoAlgo[type == "AESI", Varname], paste0, c("_narrow"))
COV_algo_input <- lapply(Map(identity, NoAlgo[type == "COV", Varname]), paste0, c("_narrow", "_possible"))

OUTCOME_concepts_in_var <- c(OUTCOME_concepts_in_var, OUTCOME_algo_input)
COV_concepts_in_var <- c(COV_concepts_in_var, COV_algo_input)
rm(COV_algo_input)

#Drugs
NoAlgo <- DRUG_codelist[!(Algorithm), ]
DRUG_concepts_in_var <- sapply(NoAlgo[(AESI) | (COV) | (Algorithm_input), Drug_proxie], identity)
rm(DRUG_codelist, NoAlgo)

variable_definition <- c(OUTCOME_concepts_in_var, COV_concepts_in_var, DRUG_concepts_in_var)
rm(OUTCOME_concepts_in_var, COV_concepts_in_var, DRUG_concepts_in_var)

for (var in names(variable_definition)) {
  if (var %in% CONCEPTSETS_to_be_split) {
    variable_definition[[var]] <- paste(var, seq_len(numbers_split[[which(CONCEPTSETS_to_be_split %in% var)]]), sep = "_")
  }
}

ALGO_codelist <- readxl::read_excel(paste0(thisdir,"/p_parameters/archive_parameters/Variables_ALG_DP_ROC20_July22.xlsx"),
                                    sheet = "ALG")
ALGO_codelist <- as.data.table(ALGO_codelist)

SECCOMPONENTS <- "B_TTS_AESI"

# algortihms for outcomes

ALGO_link <- ALGO_codelist[Algorithm %in% unique(ALGO_codelist[Algorithm != SECCOMPONENTS, Algorithm]),
                           .(Algorithm, VariableName)]

ALGO_link <- ALGO_link[, test := variable_definition[VariableName]][, VariableName := NULL]

ALGO_link <- split(ALGO_link, by = "Algorithm", keep.by = F)
ALGO_link <- lapply(ALGO_link, unlist, use.names = F)

variable_definition <- c(variable_definition, ALGO_link)
rm(ALGO_link)

# variables_with_split_conceptsets <- c()
# for (var in names(variable_definition)) {
#   if (any(variable_definition[[var]] %in% CONCEPTSETS_to_be_split)) {
#     variables_with_split_conceptsets <- append(variables_with_split_conceptsets, var)
#   }
# }

# we need to create two groups of meanings: one referring to hospitals HOSP (excluding emergency care) and one referring to primary care PC

meanings_of_this_study<-vector(mode="list")
meanings_of_this_study[["HOSP"]]=c("hospitalisation_primary","hospitalisation_secondary","hospital_diagnosis","hopitalisation_diagnosis_unspecified","episode_primary_diagnosis","episode_secondary_diagnosis","diagnosis_procedure","hospitalisation_associated","hospitalisation_linked","HH","NH","hospitalisation_ICU_primary","hospitalisation_ICU_secondary","hospitalisation_ICU_unspecified")
meanings_of_this_study[["PC"]]=c("primary_care_event","primary_care_diagnosis","primary_care_events_BIFAP","primary_care_antecedents_BIFAP","primary_care_condicionants_BIFAP")

# create two conditions on the meaning_renamed variable, associated to HOSP and to PC as listed above

condmeaning <- list()
for (level1 in c("HOSP","PC")) {
  for (meaning in meanings_of_this_study[[level1]]) {
    if (length(condmeaning[[level1]])==0) {condmeaning[[level1]]=paste0("meaning_renamed=='",meanings_of_this_study[[level1]][[1]],"'")
    }else{
      condmeaning[[level1]]=paste0(condmeaning[[level1]], " | meaning_renamed=='",meaning,"'")
    }
  }
}
rm(meanings_of_this_study, level1, meaning)

# 
#----------------------------
# SECONDARY COMPONENTS

# SECCOMPONENTS <- c("ArterialNoTP", "ArterialTP", "VTENoTP", "VTETP", "ArterialVTENoTP", "ArterialVTETP", "CVSTNoTP", "CVSTTP")

concept_set_seccomp <- vector(mode="list")
rule_seccomp <- vector(mode="list")
distance_seccomp <- vector(mode="list")
direction_seccomp <- vector(mode="list")

for (SECCOMP in SECCOMPONENTS) {
  distance_seccomp[[SECCOMP]] = '10'
  direction_seccomp[[SECCOMP]] = "Either direction"

  selectionrule_direction_seccomp <- vector(mode="list")
  selectionrule_direction_seccomp["A before B"] <- paste0("dateA <= dateB  & dateB <= dateA + ", distance_seccomp[[SECCOMP]])
  selectionrule_direction_seccomp["B before A"] <- paste0("dateB <= dateA  & dateA <= dateB + ", distance_seccomp[[SECCOMP]])
  selectionrule_direction_seccomp["Either direction"] <- paste0('((',selectionrule_direction_seccomp["A before B"],') | (',selectionrule_direction_seccomp["B before A"],'))')

}
rm(distance_seccomp, SECCOMP)

test <- ALGO_codelist[Algorithm == "B_TTS_AESI"]
rm(ALGO_codelist)

test_vect <- c()
for (a in unlist(test[group == 2, .(VariableName)])) {
  test_vect <- c(test_vect, variable_definition[[a]])
}

test_vect_1 <- c()
for (a in unlist(test[group == 1, .(VariableName)])) {
  test_vect_1 <- c(test_vect_1, variable_definition[[a]])
}

# ArterialNoTP
concept_set_seccomp[["B_TTS_AESI"]][['A']] <- test_vect
concept_set_seccomp[["B_TTS_AESI"]][['B']] <- test_vect_1
rule_seccomp[["B_TTS_AESI"]] <- "AND"

# ArterialNoTP
concept_set_seccomp[["B_TTS_AESI"]][['B']] <- test_vect
rule_seccomp[["B_TTS_AESI"]] <- "AND"

rm(test, test_vect, test_vect_1)

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

# augment ICPC codes
for (outcome in OUTCOME_variables){
  outnarrow <- paste0(outcome,'_narrow')
  outpossible <- paste0(outcome,'_possible')
  if (length(concept_set_codes_our_study[[outnarrow]][["ICPC"]]) == 0 & length(concept_set_codes_our_study[[outnarrow]][["ICPC2P"]]) > 0){
  concept_set_codes_our_study[[outpossible]][["ICPC"]] <- unique(c(concept_set_codes_our_study[[outpossible]][["ICPC"]],substr(concept_set_codes_our_study[[outnarrow]][["ICPC2P"]],1,3)))
  }
}
rm(outcome, outnarrow, outpossible)

for (conceptset in c(COV_variables, DP_variables, VACCINES_variable)){
  if (length(concept_set_codes_our_study[[conceptset]][["ICPC2P"]]) > 0){
    concept_set_codes_our_study[[conceptset]][["ICPC"]] <- unique(c(concept_set_codes_our_study[[conceptset]][["ICPC"]],substr(concept_set_codes_our_study[[conceptset]][["ICPC2P"]],1,3)))
  }
}

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


save(concept_set_codes_our_study,file=paste0(direxp,"concept_set_codes_our_study.RData"))

if (this_datasource_has_subpopulations == TRUE){
  for (subpop in subpopulations_non_empty){
    save(concept_set_codes_our_study, file = paste0(direxpsubpop[[subpop]], "concept_set_codes_our_study.RData"))
    
  }
  rm(subpop)
}

rm(a, conceptset, datasources_with_subpopulations)
