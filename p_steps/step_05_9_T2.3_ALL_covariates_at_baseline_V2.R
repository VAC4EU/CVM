# ----------------------------------
# for all covariates create binary variable drug proxy OR diagnosis; also create binary 'overall'

# input: D4_study_population_cov , D3_study_population_DP, D4_study_population
# output: D3_study_population_cov_ALL

print('create RISK FACTORS at baseline as either diagnosis or drugs')

COVnames<-c("CV","COVCANCER","COVCOPD","COVHIV","COVCKD","COVDIAB","COVOBES","COVSICKLE")

# create variable added to study population




D3_study_population_cov_ALL <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(diroutput,"D4_study_population",suffix[[subpop]],".RData")) 
  study_population<-get(paste0("D4_study_population", suffix[[subpop]])) 
  
  load(paste0(dirtemp,"D3_study_population_DP",suffix[[subpop]],".RData")) 
  study_population_DP<-get(paste0("D3_study_population_DP", suffix[[subpop]])) 
  
  load(paste0(diroutput,"D4_study_population_cov",suffix[[subpop]],".RData")) 
  study_population_cov<-get(paste0("D4_study_population_cov", suffix[[subpop]])) 
  

  study_population_cov_ALL <- merge(study_population[,-c("study_entry_date","sex")], study_population_cov, by=c("person_id", "date_of_death", "start_follow_up"), all.x = T)
  
  study_population_cov_ALL <- merge(study_population_cov_ALL, study_population_DP, by="person_id", all.x = T)
  study_population_cov_ALL <- study_population_cov_ALL[, all_covariates_non_CONTR := 0 ]
  
  # for (cov in COVnames ){
  #   if ( cov!="CV" ){
  #     nameDP =  paste0("DP_",cov,"_at_study_entry")
  #   }
  #   else{
  #     nameDP = "DP_CVD_at_study_entry"
  #   }
  #   study_population_cov_ALL <- study_population_cov_ALL[get(paste0(cov,"_at_study_entry")) == 1 | get(nameDP) == 1, namevar := 1]
  #   study_population_cov_ALL <- study_population_cov_ALL[namevar == 1 ,all_covariates_non_CONTR :=1]
  # 
  #   setnames(study_population_cov_ALL,"namevar",paste0(cov,"_either_DX_or_DP"))
  # 
  #   is.data.table(study_population_cov_ALL)
  #   for (i in names(study_population_cov_ALL)){
  #     study_population_cov_ALL[is.na(get(i)), (i):=0]
  #   }
  # }
  
  
  for (cov in COVnames ){
    if ( cov!="CV" ){
      nameDP =  paste0("DP_",cov,"_at_study_entry")
      study_population_cov_ALL <- study_population_cov_ALL[get(paste0(cov,"_at_study_entry")) == 1 | get(nameDP) == 1, namevar := 1]
    }
    else{
        nameDP1 = "DP_CVD_at_study_entry"
        nameDP2 = "DP_CONTRHYPERT_at_study_entry"
      study_population_cov_ALL <- study_population_cov_ALL[get(paste0(cov,"_at_study_entry")) == 1 | get(nameDP1) == 1  | get(nameDP2) == 1, namevar := 1]
    }
    
    study_population_cov_ALL <- study_population_cov_ALL[namevar == 1 ,all_covariates_non_CONTR :=1]
    setnames(study_population_cov_ALL,"namevar",paste0(cov,"_either_DX_or_DP"))

    is.data.table(study_population_cov_ALL)
    for (i in names(study_population_cov_ALL)){
      study_population_cov_ALL[is.na(get(i)), (i):=0]
    }
  }
  study_population_cov_ALL <- study_population_cov_ALL[IMMUNOSUPPR_at_study_entry == 1, all_covariates_non_CONTR :=1]
  
  tempname<-paste0("D3_study_population_cov_ALL",suffix[[subpop]])
  assign(tempname,study_population_cov_ALL)
  save(list=tempname,file=paste0(diroutput,tempname,".RData"))
  
  rm(list=paste0("D4_study_population", suffix[[subpop]]))
  rm(list=paste0("D3_study_population_cov_ALL",suffix[[subpop]]))
  rm(list=paste0("D3_study_population_DP", suffix[[subpop]]))
  rm(list=paste0("D4_study_population_cov",suffix[[subpop]]))
     
}



rm( study_population_DP, study_population,study_population_cov, study_population_cov_ALL)

