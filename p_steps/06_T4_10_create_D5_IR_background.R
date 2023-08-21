for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  #D4_persontime_risk_year-----------------------------------------------------
  
  load(paste0(diroutput,"D4_persontime_background_aggregated",suffix[[subpop]],".RData"))
  persontime_windows<-get(paste0("D4_persontime_background_aggregated", suffix[[subpop]]))
  rm(list=paste0("D4_persontime_background_aggregated", suffix[[subpop]]))
  
  for (ev in c(OUTCOME_variables, CONTROL_variables, "DEATH")) {
    name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
    name_count <- paste0(ev,"_b")
    name_pt <- paste0("Persontime_",ev)
    persontime_windows[, (name_cols) := exactPoiCI(persontime_windows, name_count, name_pt)]
  }
  
  nameoutput <- paste0("D5_IR_background")
  assign(nameoutput, persontime_windows)
  save(nameoutput, file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".RData"), list = nameoutput)
  
  fwrite(get(nameoutput), file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".csv"))
}
