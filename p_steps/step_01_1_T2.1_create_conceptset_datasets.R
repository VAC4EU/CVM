
# APPLY THE FUNCTION CreateConceptSetDatasets TO CREATE ONE DATASET PER CONCEPT SET CONTAINING ONLY RECORDS WITH A CODE OF INTEREST

# input: EVENTS, MEDICINES, SURVEY_OBSERVATIONS, MEDICAL_OBSERVATIONS
# output: concept set datasets, one per concept set, named after the concept set itself


print('RETRIEVE FROM CDM RECORDS CORRESPONDING TO CONCEPT SETS')
CreateConceptSetDatasets(concept_set_names = c(vaccine__conceptssets),
                         dataset = ConcePTION_CDM_tables,
                         codvar = ConcePTION_CDM_codvar,
                         datevar= ConcePTION_CDM_datevar,
                         dateformat= "YYYYmmdd",
                         rename_col = list(person_id=person_id,date=date),
                         concept_set_domains = concept_set_domains,
                         concept_set_codes =	concept_set_codes_our_study,
                         discard_from_environment = T,
                         dirinput = dirinput,
                         diroutput = dirtemp,
                         extension = c("csv")
                         )

CreateConceptSetDatasets(concept_set_names = c(concept_sets_of_our_study),
                         dataset = ConcePTION_CDM_tables,
                         codvar = ConcePTION_CDM_codvar,
                         datevar = ConcePTION_CDM_datevar,
                         EAVtables = ConcePTION_CDM_EAV_tables,
                         EAVattributes = ConcePTION_CDM_EAV_attributes_this_datasource,
                         dateformat= "YYYYmmdd",
                         vocabulary = ConcePTION_CDM_coding_system_cols,
                         rename_col = list(person_id=person_id,date=date),
                         concept_set_domains = concept_set_domains,
                         concept_set_codes =	concept_set_codes_our_study,
                         concept_set_codes_excl = concept_set_codes_our_study_excl,
                         discard_from_environment = T,
                         dirinput = dirinput,
                         diroutput = dirtemp,
                         extension = c("csv"),
                         vocabularies_with_dot_wildcard=c("READ"),
                         vocabularies_with_exact_search = "Free_text")


# Works if concept has _free_text which contains two _ before the real concept
for (ft_concept in FREE_TEXT_conceptsets) {
  original_concept <- paste(head(strsplit(ft_concept, "_")[[1]],-2), collapse = "_")
  load(paste0(dirtemp,original_concept,".RData"))
  load(paste0(dirtemp,ft_concept,".RData"))
  final_concept <- unique(rbind(get(original_concept), get(ft_concept), fill=TRUE))
  assign(original_concept, final_concept)
  save(original_concept, file = paste0(dirtemp, "/", original_concept, ".RData"), list = original_concept)
  rm(final_concept, list = c(original_concept, ft_concept))
}

