#------------------------------------------------------------------
# create events and components of OUTCOMES 

# input: concept set datasets of outcomes (narrow and possible), D4_study_population
# output: for each outcome OUTCOME, D3_components_OUTCOME.RData and D3_events_OUTCOME_type.RData, for type = narrow, possible

print('create events and create components of OUTCOMES and CONTROLS')

# 
# load(paste0(dirpargen,"subpopulations_non_empty.RData"))

##for each var in OUTCOME and for each negative outcome create D3_var including all dates when that outcome is observed (use the corresponding conceptsets)

for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  COHORT_TMP <- get(load(paste0(diroutput, "D4_study_population", suffix[[subpop]], ".RData"))[[1]])
  rm(list = paste0("D4_study_population", suffix[[subpop]]))
  
  COHORT_TMP <- as.data.table(COHORT_TMP)
  COHORT_TMP <- COHORT_TMP[, .(person_id, study_entry_date)]
  
  for (OUTCOME in c(OUTCOME_variables[OUTCOME_variables %not in% SECCOMPONENTS], CONTROL_variables)) {
    print(OUTCOME)
    
    summarystatOUTCOME <- list()
    addvarOUTCOME <- list()
    
    # for (year in ComponentAnalysisYears) {
    #   FirstJan_1 <- ymd(paste0(year, "0101"))
    #   
    #   for (level1 in c("HOSP", "PC")) {
    #     
    #     namenewvar <- paste(OUTCOME, level1, year, sep = "_")
    #     summarystatOUTCOME <- append(summarystatOUTCOME, list(list(c("max"), namenewvar, namenewvar)))
    #     
    #     condition_str <- paste0("(", condmeaning[[level1]], ") & date <= as.Date('", FirstJan_1,
    #                             "') + 365 & date >= as.Date('", FirstJan_1, "')")
    #     addvarOUTCOME <- append(addvarOUTCOME, list(list(c(namenewvar), "1", condition_str)))
    #     addvarOUTCOME <- append(addvarOUTCOME, list(list(c(namenewvar), "0", paste0("is.na(", namenewvar, ")"))))
    #     
    #   }
    # }
    
    for (year in ComponentAnalysisYears) {
      
      for (level1 in c("HOSP", "PC")) {
        
        namenewvar <- paste(OUTCOME, level1, year, sep = "_")
        summarystatOUTCOME <- append(summarystatOUTCOME, list(list(c("max"), namenewvar, namenewvar)))
        
        condition_str <- paste0("(", condmeaning[[level1]], ") & year(date) == ", year)
        addvarOUTCOME <- append(addvarOUTCOME, list(list(c(namenewvar), "1", condition_str)))
        addvarOUTCOME <- append(addvarOUTCOME, list(list(c(namenewvar), "0", paste0("is.na(", namenewvar, ")"))))
        
      }
    }
    
    selectionOUTCOME <- "date >= study_entry_date - 365"
    
    # delete records that are not observed in this whole subpopulation
    if (this_datasource_has_subpopulations == TRUE){
      selectionOUTCOME <- paste(selectionOUTCOME, '&', select_in_subpopulationsEVENTS[[subpop]])
    }
    
    # delete records whose meanings are not appropriate for a specific datasource 
    selectionOUTCOME <- paste0(selectionOUTCOME, ' & !(', select_meanings_AESI[[thisdatasource]], ')')
    
    nameconceptsetdatasetOUTCOMEtype <- variable_definition[[OUTCOME]]
    conceptsets_list <- lapply(nameconceptsetdatasetOUTCOMEtype,
                               function(x) get(load(paste0(dirconceptsets, x,".RData"))[[1]]))
    
    components <- MergeFilterAndCollapse(
      listdatasetL = conceptsets_list,
      condition = selectionOUTCOME,
      key = "person_id",
      datasetS = COHORT_TMP,
      additionalvar = addvarOUTCOME,
      nameintermediatedataset = paste0(dirtemp,'tempfile'),
      strata = "person_id",
      saveintermediatedataset = T,
      summarystat = list(list("count", "person_id", "count"))
    )
    
    load(paste0(dirtemp, 'tempfile.RData'))
    
    nameobjectOUTCOMEtype <- paste0('D3_events_', OUTCOME, '_simple', suffix[[subpop]])
    assign(nameobjectOUTCOMEtype, tempfile)
    save(nameobjectOUTCOMEtype, file = paste0(direvents, nameobjectOUTCOMEtype, ".RData"), list = nameobjectOUTCOMEtype)
    rm(components, tempfile, nameobjectOUTCOMEtype, list = nameobjectOUTCOMEtype)
    
    # components <- MergeFilterAndCollapse(
    #   listdatasetL = conceptsets_list,
    #   condition = selectionOUTCOME,
    #   key = "person_id",
    #   datasetS = COHORT_TMP,
    #   additionalvar = addvarOUTCOME,
    #   strata = "person_id",
    #   summarystat = summarystatOUTCOME
    # )
    # 
    # componentsOUTCOMEfinal <- merge(COHORT_TMP, components, by = "person_id", all.x  = T)
    # setnafill(componentsOUTCOMEfinal, fill = 0)
    # 
    # nameobjectOUTCOME <- paste0("D3_components_", OUTCOME, suffix[[subpop]])
    # assign(nameobjectOUTCOME, componentsOUTCOMEfinal)
    # 
    # save(nameobjectOUTCOME, file = paste0(dircomponents, nameobjectOUTCOME, ".RData"), list = nameobjectOUTCOME)
    # rm(components, componentsOUTCOMEfinal, componentsOUTCOMEfinal, nameobjectOUTCOME, list = nameobjectOUTCOME)
    rm(conceptsets_list)
  }
}