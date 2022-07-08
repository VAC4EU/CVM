#------------------------------------------------------------------
# create age and other simple covariates at date of first vaccination

# input: D3_Vaccin_cohort_no_risk, D3_Vaccin_cohort_covariates
# output: D4_Vaccin_cohort_cov

print('create age and other simple covariates at date_vax_1')

D4_Vaccin_cohort_cov <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D3_Vaccin_cohort_no_risk",suffix[[subpop]],".RData"))
  load(paste0(dirtemp,"D3_Vaccin_cohort_covariates",suffix[[subpop]],".RData"))
  
  study_population<-get(paste0("D3_Vaccin_cohort_no_risk", suffix[[subpop]]))
  study_population_covariates<-get(paste0("D3_Vaccin_cohort_covariates", suffix[[subpop]]))
  
  population_var<- study_population[, year_at_date_vax_1:=year(date_vax1)] 
  population_var<-population_var[, age_strata_at_date_vax_1 := cut(age_at_date_vax_1, breaks = Agebands, 
                                                                    labels = Agebands_labels)]
  
  study_population_cov<-merge(population_var, study_population_covariates[,-"date_vax1"],
                              by="person_id", all.x = T)
  
  study_population_cov <- study_population_cov[,-c("date_of_birth","date_vax1","study_exit_date")]

  tempname<-paste0("D4_Vaccin_cohort_cov",suffix[[subpop]])
  assign(tempname,study_population_cov)
  save(list=tempname,file=paste0(diroutput,tempname,".RData"))
  
  rm(list=paste0("D4_Vaccin_cohort_cov", suffix[[subpop]]))
  rm(list=paste0("D3_Vaccin_cohort_no_risk", suffix[[subpop]]))
  rm(list=paste0("D3_Vaccin_cohort_covariates", suffix[[subpop]]))
  
} 

rm(population_var, study_population_cov,study_population,study_population_covariates)

