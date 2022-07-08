#------------------------------------------------------------------
# create age and other simple covariates at baseline

# input: D4_population_c_no_risk, D3_population_c_covariates
# output: D4_population_c_cov

print('create age and other simple covariates at date_vax_1')

D4_Vaccin_cohort_cov <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D4_population_c_no_risk",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_population_c_covariates",suffix[[subpop]],".RData"))
  
  study_population<-get(paste0("D4_population_c_no_risk", suffix[[subpop]]))
  study_population_covariates<-get(paste0("D3_population_c_covariates", suffix[[subpop]]))
  
  population_var<- study_population[, year_at_date_vax_1:=year(cohort_entry_date_MIS_c)]
  
  study_population_cov<-merge(population_var, study_population_covariates[,-"cohort_entry_date_MIS_c"],
                              by="person_id", all.x = T)
  
  study_population_cov <- study_population_cov[,-c("cohort_entry_date_MIS_c")]

  tempname<-paste0("D4_population_c_cov",suffix[[subpop]])
  assign(tempname,study_population_cov)
  save(list=tempname,file=paste0(diroutput,tempname,".RData"))
  
  rm(list=paste0("D4_population_c_cov", suffix[[subpop]]))
  rm(list=paste0("D4_population_c_no_risk", suffix[[subpop]]))
  rm(list=paste0("D3_population_c_covariates", suffix[[subpop]]))
  
} 

rm(population_var, study_population_cov,study_population,study_population_covariates)

