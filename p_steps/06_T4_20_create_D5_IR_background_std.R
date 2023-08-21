## Reference population
pop.eustat<- read.csv(paste0(dirmacro,"/ESP_ageband.csv"),sep = "")

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  #D4_persontime_risk_year-----------------------------------------------------
  
  load(paste0(diroutput,"D4_persontime_background_aggregated",suffix[[subpop]],".RData"))
  persontime_windows<-get(paste0("D4_persontime_background_aggregated", suffix[[subpop]]))
  rm(list=paste0("D4_persontime_background_aggregated", suffix[[subpop]]))
  persontime_windows <- persontime_windows[Ageband != "total" & sex == 'total' & year %in% c("2019","2020"),]
  persontime_windows <- persontime_windows[year == "2019",stratum_bkr := 1]
  persontime_windows <- persontime_windows[year == "2020" & COVID19 == 0 ,stratum_bkr := 2]
  persontime_windows <- persontime_windows[year == "2020"  & COVID19 == 1,stratum_bkr := 3]
  colstokeep <- c('COVID19','year')
  for (ev in c(OUTCOME_variables, CONTROL_variables, "DEATH")) {
    name_count <- paste0(ev,"_b")
    name_pt <- paste0("Persontime_",ev)
    suppressMessages(
    my_results_CVM <- dsr(data = persontime_windows,
                          event = get(name_count),
                          fu = get(name_pt),
                          subgroup = stratum_bkr,
                          Ageband,
                          refdata = pop.eustat,
                          method = "gamma",
                          sig = 0.95,
                          mp = 36525000, # 100,000 * 365.25
                          decimals = 2)
    )
    my_results_CVM <- my_results_CVM[,c(1,7,8,9) ]
    colnames(my_results_CVM) <- c('stratum_bkr',paste0(c("IR_std_", "lb_std_", "ub_std_"), ev))
    my_results_CVM <- as.data.table(my_results_CVM)
    my_results_CVM <- my_results_CVM[stratum_bkr == 1, COVID19 := 0   ]
    my_results_CVM <- my_results_CVM[stratum_bkr == 1, year := '2019'  ]
    my_results_CVM <- my_results_CVM[stratum_bkr == 2, COVID19 := 0   ]
    my_results_CVM <- my_results_CVM[stratum_bkr == 2, year := '2020'  ]
    my_results_CVM <- my_results_CVM[stratum_bkr == 3, COVID19 := 1   ]
    my_results_CVM <- my_results_CVM[stratum_bkr == 3, year := '2020'  ]
    persontime_windows <- merge(persontime_windows,my_results_CVM,by = c('COVID19','year','stratum_bkr'))
    colstokeep <- c(colstokeep,paste0(c("IR_std_", "lb_std_", "ub_std_"), ev))
  }
  persontime_windows <- persontime_windows[,..colstokeep]
  # persontime_windows <- persontime_windows[,-stratum_bkr]
  persontime_windows <- unique(persontime_windows)
  
  nameoutput <- paste0("D5_IR_background_std")
  assign(nameoutput, persontime_windows)
  save(nameoutput, file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".RData"), list = nameoutput)
  
  fwrite(get(nameoutput), file = paste0(dirD4D5subpop[[subpop]], nameoutput, ".csv"))
}

# data <- D4_persontime_background_aggregated_HOSP[Ageband != "total" & sex == 'total' & year == "2019/2020",]
# 
# my_results_CVM <- dsr(data = data,
#                       event = B_COAGDIS_AESI_b,
#                       fu = Persontime_B_COAGDIS_AESI,
#                       subgroup = COVID19,
#                       Ageband,
#                       refdata = pop.eustat,
#                       method = "gamma",
#                       sig = 0.95,
#                       mp = 36525000, # 100,000 * 365.25
#                       decimals = 2)