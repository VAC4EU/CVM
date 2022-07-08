# ----------------------------------
# for all covariates create binary variable drug proxy OR diagnosis; also create binary 'overall'

# input: D4_population_c_no_risk, D4_population_c_cov, D3_population_c_DP
# output: D3_population_c_cov_ALL

print('create RISK FACTORS at date_vax_1 as either diagnosis or drugs')

COVnames<-c("CV","COVCANCER","COVCOPD","COVHIV","COVCKD","COVDIAB","COVOBES","COVSICKLE")

# create variable added to study population

D3_Vaccin_cohort_cov_ALL <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D4_population_c_no_risk",suffix[[subpop]],".RData"))
  load(paste0(diroutput,"D4_population_c_cov",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_population_c_DP",suffix[[subpop]],".RData"))
  
  study_population<-get(paste0("D4_population_c_no_risk", suffix[[subpop]]))
  study_population_cov<-get(paste0("D4_population_c_cov", suffix[[subpop]]))
  study_population_DP<-get(paste0("D3_population_c_DP", suffix[[subpop]]))
  
  keep_col <- c("person_id", paste0(COVnames, "_at_covid"))
  study_population_cov_ALL <- merge(study_population, study_population_cov[, ..keep_col], by=c("person_id"), all.x = T)
  
  study_population_cov_ALL <- merge(study_population_cov_ALL, study_population_DP[, cohort_entry_date_MIS_c := NULL], by="person_id", all.x = T)
  
  study_population_cov_ALL <- study_population_cov_ALL[, all_covariates_non_CONTR := 0 ]
  
  for (cov in COVnames){
    if ( cov!="CV" ){
      nameDP =  paste0("DP_",cov,"_at_covid")
      study_population_cov_ALL <- study_population_cov_ALL[get(paste0(cov,"_at_covid")) == 1 | get(nameDP) == 1, namevar := 1]
      
    }
    else{
      nameDP1 = "DP_CVD_at_covid"
      nameDP2 = "DP_CONTRHYPERT_at_covid"
      study_population_cov_ALL <- study_population_cov_ALL[get(paste0(cov,"_at_covid")) == 1 | get(nameDP1) == 1  | get(nameDP2) == 1, namevar := 1]
    }

    study_population_cov_ALL <- study_population_cov_ALL[namevar == 1 ,all_covariates_non_CONTR :=1]
    
    setnames(study_population_cov_ALL,"namevar",paste0(cov,"_either_DX_or_DP"))
    
    for (i in names(study_population_cov_ALL)){
      study_population_cov_ALL[is.na(get(i)), (i):=0]
    }
  }
  
  study_population_cov_ALL <- study_population_cov_ALL[IMMUNOSUPPR_at_covid == 1, all_covariates_non_CONTR :=1]
  
  tempname<-paste0("D3_population_c_cov_ALL",suffix[[subpop]])
  assign(tempname,study_population_cov_ALL)
  save(list=tempname,file=paste0(diroutput,tempname,".RData"))
  
  rm(list=paste0("D3_population_c_cov_ALL", suffix[[subpop]]))
  rm(list=paste0("D4_population_c_no_risk", suffix[[subpop]]))
  rm(list=paste0("D3_population_c_DP", suffix[[subpop]]))
  rm(list=paste0("D4_population_c_cov", suffix[[subpop]]))
}


rm( study_population_DP, study_population,study_population_cov, study_population_cov_ALL)

