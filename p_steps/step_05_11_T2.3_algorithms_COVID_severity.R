# SUMMARIZE ALGORITHMS FOR COVID SEVERITY
#-----------------------------------------------
# input: D3_components_covid_severity, D4_study_population
# output: D3_algorithm_covid, table_QC_covid_diagnosis, D3_outcomes_covid_multiple, D3_outcomes_covid

print("SUMMARIZE ALGORITHMS FOR COVID SEVERITY")




D3_algorithm_covid <- vector(mode = 'list')
D3_outcomes_covid <- vector(mode = 'list')
list_outcomes_observed_COVID <- c()
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D3_components_covid_severity",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData"))
  
  algorithm_covid<-get(paste0("D3_components_covid_severity", suffix[[subpop]])) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]]))
  
   # date: covid registry
  algorithm_covid <- algorithm_covid[!is.na(first_date_covid_registry), date_covid := first_date_covid_registry]
  algorithm_covid <- algorithm_covid[!is.na(date_covid), origin_date_covid:= "covid_registry"]
  # date: covid narrow (except for BIFAP)
  if (thisdatasource != 'BIFAP'){
    algorithm_covid <- algorithm_covid[is.na(date_covid) & !is.na(first_date_covid_narrow),date_covid := first_date_covid_narrow]
    algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(origin_date_covid),origin_date_covid:="covid_narrow"]
  }
  # SEVERITY LEVEL
  
  # level 5: level from covid registry
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & covid_registry_symptoms == 5, severity_level_covid := 5][!is.na(date_covid) & !is.na(severity_level_covid), origin_severity_level_covid :='5-covid_registry']
  # level 5: death at discharge from hospital
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & covid_death_discharge != 0,severity_level_covid := 5][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='5-death_at_discharge']
  # level 5: death after covid registry
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & death_after_covid_registry != 0,severity_level_covid :=5][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='5-after_covid_registry']
  # level 5: death within days from covid narrow (except for BIFAP)
  if (thisdatasource != 'BIFAP'){
    algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & death_after_covid_narrow_date != 0, severity_level_covid := 5][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='5-death_after_covid_narrow_date']
  }
  
  # level 4: level from covid registry
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & covid_registry_symptoms == 4,severity_level_covid :=4][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='4-covid_registry']
  # level 4: conceptset from covid registry date
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & 
is.na(severity_level_covid) & MechanicalVentilation_within_registry_date != 0, severity_level_covid := 4][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='4-after_covid_registry']
  # level 4: conceptset from covid narrow date
  if (thisdatasource != 'BIFAP'){ 
    algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & MechanicalVentilation_within_covid_narrow_date != 0,severity_level_covid := 4][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='4-after_covid_narrow']
  }
  
  # level 3: level from covid registry
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & covid_registry_symptoms == 3,severity_level_covid :=3][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='3-covid_registry']
  # level 3: conceptset from covid registry 
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & Respiratory_within_registry_date != 0,severity_level_covid := 3][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='3-after_covid_registry']
  # level 3: conceptset from covid narrow  
  if (thisdatasource != 'BIFAP'){ 
    algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & Respiratory_within_covid_narrow_date != 0,severity_level_covid := 3][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='3-after_covid_narrow']
  }
  
  # level 2: level from covid registry
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & covid_registry_symptoms == 2,severity_level_covid :=2][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='2-covid_registry']
  # level 2: conceptset from covid registry 
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & Infection_within_registry_date != 0,severity_level_covid := 2][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='2-after_covid_registry']
  # level 2: conceptset from covid narrow  
  if (thisdatasource != 'BIFAP'){  
    algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) & Infection_within_covid_narrow_date != 0,severity_level_covid := 2][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='2-after_covid_narrow']
  }
  
  # level 1: level by exclusion
  algorithm_covid <- algorithm_covid[!is.na(date_covid) & is.na(severity_level_covid) ,severity_level_covid :=1][!is.na(date_covid) & !is.na(severity_level_covid) & is.na(origin_severity_level_covid), origin_severity_level_covid :='1-by_exclusion']
  
  
  # cumulative binary variables
  outcomes_covid <- data.table()
  list_outcomes_observed_COVID <- c()
  for (j in c(1,2,3,4,5) ){
    temp <- algorithm_covid[severity_level_covid >= j]
    level <- paste0('COVID_L',j,'plus')
    temp <- temp[, name_event := level]
    temp <- temp[, date_event := date_covid]
    temp <- temp[,.(person_id,name_event,date_event,origin_severity_level_covid)]
    outcomes_covid <- rbind(outcomes_covid, temp)
    list_outcomes_observed_COVID <- c(list_outcomes_observed_COVID, level)
      
  }
  
  outcomes_covid_wrong <- outcomes_covid[date_event < start_COVID_diagnosis_date, ][, covid_year := year(date_event)][, covid_month := month(date_event)]
  outcomes_covid_wrong <- outcomes_covid_wrong[, .N, by = c("covid_year", "covid_month")]
  setorder(outcomes_covid_wrong, covid_year, covid_month)
  fwrite(outcomes_covid_wrong, file = paste0(direxp, "table_QC_covid_diagnosis.csv"))
  rm(outcomes_covid_wrong)
  outcomes_covid <- outcomes_covid[date_event >= start_COVID_diagnosis_date, ]

  
  # save the COVID outcomes as a dataset and their list as a parameter
  tempname<-paste0("list_outcomes_observed_COVID",suffix[[subpop]])
  assign(tempname,list_outcomes_observed_COVID)
  save(list=tempname,file=paste0(dirpargen,tempname,".RData"))
  rm(tempname,list=tempname)
  
  tempname<-paste0("D3_algorithm_covid",suffix[[subpop]])
  assign(tempname,algorithm_covid)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  rm(tempname,list=tempname)

  outcomes_covid_multiple<-outcomes_covid
  tempname<-paste0("D3_outcomes_covid_multiple",suffix[[subpop]])
  assign(tempname,outcomes_covid_multiple)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  rm(tempname,list=tempname)
  
  outcomes_covid<-unique(outcomes_covid)
  tempname<-paste0("D3_outcomes_covid",suffix[[subpop]])
  assign(tempname,outcomes_covid)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  rm(tempname,list=tempname)
  
  temp2<-paste0("D3_components_covid_severity", suffix[[subpop]])
  rm(temp2,list=temp2)
  temp2<-paste0("D4_study_population", suffix[[subpop]])
  rm(temp2,list=temp2)
  
}

rm(study_population,algorithm_covid,temp,outcomes_covid,outcomes_covid_multiple)

