#-----------------------------------------------
# set D3_events_ALL_OUTCOMES which contains the first outcome per person

# input: D4_study_population, D3_events_OUTCOME/CONTROL
# output: D3_events_ALL_OUTCOMES

for (subpop in subpopulations_non_empty) {
  
  variable_list <- lapply(c(OUTCOME_variables, CONTROL_variables, "DEATH"), function(x) {
    print(paste("Merging variable:", x))
    
    algo_suffix <- fifelse(x %in% SECCOMPONENTS, "_complex", "_simple")
    if (x == "DEATH") {
      df <- get(load(paste0(direvents, "D3_events_", x, ".RData"))[[1]])
    } else {
      df <- get(load(paste0(direvents, "D3_events_", x, algo_suffix, suffix[[subpop]], ".RData"))[[1]])
    }
    
    if (nrow(df) == 0) {
      df <- data.table(person_id = character(),
                       date = as.Date(as.POSIXct(character())),
                       type_outcome = character(),
                       meaning_renamed = character(),
                       codvar = character(),
                       event_record_vocabulary = character())
    } else {
      if (x %in% SECCOMPONENTS) {
        setnames(df, c("meaning_renamedA", "codvarA", "event_record_vocabularyA"),
                 c("meaning_renamed", "codvar", "event_record_vocabulary"))
      }
      if (x != "DEATH") {
        df <- df[, .(person_id, date, type_outcome = x, meaning_renamed,
                     codvar, event_record_vocabulary)]
      }
      if (x == "DEATH") {
        df <- df[, .(person_id, date, type_outcome = x)]
      }
    }
    return(df)
  })
  
  events_ALL_OUTCOMES <- rbindlist(variable_list, fill = T)
  
  nametemp <- paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]])
  assign(nametemp, events_ALL_OUTCOMES)
  save(nametemp, file = paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"), list = nametemp)
}
