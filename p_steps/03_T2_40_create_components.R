#------------------------------------------------------------------
# split outcomes in components

# input: D4_study_population, CONCEPT, named OUTCOME
# output: QC_all_components_OUTCOME.RData

print('BREAK DOWN OUTCOMES PER MEANINGS')

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  COHORT_TMP <- get(load(paste0(diroutput, "D4_study_population", suffix[[subpop]], ".RData"))[[1]])
  rm(list = paste0("D4_study_population", suffix[[subpop]]))
  
  COHORT_TMP <- as.data.table(COHORT_TMP)
  COHORT_TMP <- COHORT_TMP[, .(person_id, study_entry_date, study_exit_date)]
  
  firstyear = min(as.numeric(ComponentAnalysisYears))
  secondyear = max(as.numeric(ComponentAnalysisYears))
  
  COHORT_TMP <- COHORT_TMP[year(study_entry_date) <= secondyear & year(study_exit_date) >= firstyear,]
  
  for (OUTCOME in OUTCOME_variables) {
    print(OUTCOME)
    OUTCOME_components <- vector(mode = 'list')
    
    import_and_clean_components <- function(x) {
      if (!file.exists(paste0(dirconceptsets, x, ".RData"))) {
        return(NULL)
      }
      
      y <- get(load(paste0(dirconceptsets, x, ".RData"))[[1]])
      
      narrow_or_possible <- last(unlist(strsplit(x, "_")))
      if (narrow_or_possible %not in% c("narrow", "possible")) narrow_or_possible <- NA
      y[, type_concept_set := narrow_or_possible]
      
      if ('mo_meaning' %in% colnames(x)) y <- y[is.na(meaning_renamed), meaning_renamed := mo_meaning]
      y <- y[is.na(meaning_renamed), meaning_renamed := "survey"]
      y[, year := year(date)] 
      y[, n := 1]
      
      return(y)
    }
    
    nameconceptsetdatasetOUTCOMEtype <- component_definition[[OUTCOME]]
    conceptsets_list <- lapply(nameconceptsetdatasetOUTCOMEtype, import_and_clean_components)
    
    OUTCOME_detailed_components <- MergeFilterAndCollapse(listdatasetL = conceptsets_list,
                                                          datasetS = COHORT_TMP,
                                                          condition = paste0("date >= study_entry_date - 365 & year >= ",
                                                                             firstyear, " & year <= ", secondyear),
                                                          key = "person_id",
                                                          strata = c("person_id", "type_concept_set",
                                                                     "meaning_renamed", "year"),
                                                          summarystat = list(list("max", "n", "has_component")))
    
    OUTCOME_detailed_components[, todrop := fifelse(year == firstyear, 1, 0)]
    # OUTCOME_detailed_components <- merge(OUTCOME_detailed_components,
    #                                      OUTCOME_detailed_components[, .(todrop_min = min(todrop)), by = "person_id"],
    #                                      all.x = T, by = "person_id")
    # OUTCOME_detailed_components[, todrop := fifelse(todrop_min == 1, 0, 1)]
    
    OUTCOME_todrop <- OUTCOME_detailed_components[todrop == 1, ]
    OUTCOME_todrop <- unique(OUTCOME_todrop[, .(person_id, todrop)])
    
    OUTCOME_detailed_components <- OUTCOME_detailed_components[todrop == 0, ]
    OUTCOME_detailed_components <- OUTCOME_detailed_components[, component := paste0(type_concept_set, '_', meaning_renamed)]
    OUTCOME_detailed_components <- OUTCOME_detailed_components[, .(person_id, component, has_component)]
    
    # TODO test for empty datasets
    if (nrow(OUTCOME_detailed_components) > 0) {
      OUTCOME_reshaped <- data.table::dcast(OUTCOME_detailed_components, person_id ~ component, value.var = "has_component")
      
      OUTCOME_merged <- merge(COHORT_TMP, OUTCOME_reshaped, by = "person_id", all.x =T)
      rm(OUTCOME_reshaped)
    } else {
      OUTCOME_merged <- COHORT_TMP
    }
    
    OUTCOME_merged <- merge(OUTCOME_merged, OUTCOME_todrop, by = "person_id", all.x =T)
    OUTCOME_merged[, todrop := nafill(todrop, fill = 0)]
    
    cols_to_aggregate_by = setdiff(colnames(OUTCOME_merged), c("person_id", "study_entry_date", "study_exit_date"))
    OUTCOME_aggregated <- OUTCOME_merged[, .(.N), by = cols_to_aggregate_by]
    setnafill(OUTCOME_aggregated, fill = 0)
    
    nameobject <- paste0("QC_all_components_", OUTCOME)
    assign(nameobject, OUTCOME_aggregated)
    
    fwrite(get(nameobject), file = paste0(dircomponents[[subpop]], paste0(nameobject, ".csv")))
    save(nameobject, file = paste0(dirtemp, paste0(nameobject, ".RData")), list = nameobject)
    rm(nameobject, list = nameobject)
    
    rm(OUTCOME_detailed_components, OUTCOME_todrop, OUTCOME_merged, OUTCOME_aggregated)
  }
}
