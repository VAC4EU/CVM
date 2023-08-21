# CREATE COVID SEVERITY 
#-----------------------------------------------
# input: D4_study_population,  D3_covid_episodes, D3_covid_severity_components_death, D3_covid_severity_components_ICU, D3_covid_severity_components_hospitalisation
# output: D3_covid_episodes_severity

# in this step the list of unique episodes of covid is labelled with levels of severity

print("CREATE COVID SEVERITY")

# "covid_severity_1" "covid_severity_2" "covid_severity_3" "covid_severity_4" "covid_severity_5"

for (subpop in subpopulations_non_empty) { 
  #---------------------------------
  #---------------------------------
  # 1 rbind all files that imply death
  #---------------------------------
  #---------------------------------
  

  load(paste0(dirtemp,"D3_covid_episodes",suffix[[subpop]],".RData"))
  COVID_episodes <- as.data.table(get(paste0("D3_covid_episodes",suffix[[subpop]])))
  
  covid_severity <- COVID_episodes

  for (lev in c("hospitalisation","ICU","death")){
    load( paste0(dirtemp,"D3_covid_severity_components_",lev,suffix[[subpop]],".RData"))
    levsev <- as.data.table(get(paste0("D3_covid_severity_components_",lev,suffix[[subpop]])))
    levsev <- levsev[,lev:= 1]
    levsev <-levsev[,.(person_id,date, lev)]
    setnames(levsev,c("lev"),c(lev))
    covid_severity <- merge(covid_severity,levsev,by = c("person_id","date"), all.x = TRUE)[]
    rm(levsev)
  }
  
  covid_severity <- covid_severity[death == 1, severity := 4]
  covid_severity <- covid_severity[is.na(severity) & ICU == 1 , severity := 3]
  covid_severity <- covid_severity[is.na(severity) & hospitalisation == 1 , severity := 2]
  covid_severity <- covid_severity[is.na(severity) , severity := 1]
  
  covid_severity <- covid_severity[severity == 4, value_of_variable := 'severity4']
  covid_severity <- covid_severity[severity == 3, value_of_variable := 'severity3']
  covid_severity <- covid_severity[severity == 2, value_of_variable := 'severity2']
  covid_severity <- covid_severity[severity == 1, value_of_variable := 'severity1']
  
  covid_severity <- covid_severity[,.(person_id,date,value_of_variable)]
  
  tempname<-paste0("D3_TD_variable_COVID",suffix[[subpop]])
  assign(tempname,covid_severity)
  save(list=tempname,file = paste0(dirTD,tempname,".RData"))
  
  # covid_severity_vert <- emptydataset
  # 
  # for (j in c(1,2,3,4) ){
  #   level <- paste0('covid_severity_',j,'_plus')
  #   covid_severity <- covid_severity[severity >= j , level := 1]
  #   setnames(covid_severity,c("level"),c(level))
  #   
  #   temp <- covid_severity[severity >= j]
  #   temp <- temp[, name_event := level]
  #   temp <- temp[, date_event := date]
  #   temp <- temp[,.(person_id,name_event,date_event)]
  #   covid_severity_vert <- rbind(covid_severity_vert, temp, fill = TRUE)
  #   rm(temp)
  #   temp <- covid_severity[severity == j]
  #   temp <- temp[, date_event := date]
  #   temp <- temp[,.(person_id,date_event)]
  #   temp <- temp[, name_event := paste0('covid_severity_',j)]
  #   covid_severity_vert <- rbind(covid_severity_vert, temp, fill = TRUE)
  #   rm(temp)
  # }  
  # tempname<-paste0("D3_outcomes_severity_episodes_covid",suffix[[subpop]])
  # assign(tempname,covid_severity_vert)
  # save(list=tempname,file = paste0(dirtemp,tempname,".RData"))
  # 
  # tempname<-paste0("D3_severity_episodes_covid",suffix[[subpop]])
  # assign(tempname,covid_severity)
  # save(list=tempname,file = paste0(dirtemp,tempname,".RData"))
  
  
  
}
