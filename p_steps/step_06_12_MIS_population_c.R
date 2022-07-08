# ----------------------------------
# Calculate the population for cohort c

# input: D4_population_c_no_risk, D3_population_c_cov_ALL
# output: D4_population_c

for (subpop in subpopulations_non_empty) {
  
  load(paste0(dirtemp,"D4_population_c_no_risk",suffix[[subpop]],".RData"))
  D4_population_c_no_risk<-get(paste0("D4_population_c_no_risk", suffix[[subpop]]))
  load(paste0(diroutput,"D3_population_c_cov_ALL",suffix[[subpop]],".RData"))
  D3_population_c_cov_ALL<-get(paste0("D3_population_c_cov_ALL", suffix[[subpop]]))
  
  cols_risk_factors = c("CV_either_DX_or_DP", "COVCANCER_either_DX_or_DP", "COVCOPD_either_DX_or_DP",
                        "COVHIV_either_DX_or_DP", "COVCKD_either_DX_or_DP", "COVDIAB_either_DX_or_DP",
                        "COVOBES_either_DX_or_DP", "COVSICKLE_either_DX_or_DP", "all_covariates_non_CONTR")
  
  D3_population_c_cov_ALL <- D3_population_c_cov_ALL[, c("person_id", ..cols_risk_factors, "IMMUNOSUPPR_at_covid")]
  
  setnames(D3_population_c_cov_ALL, c(cols_risk_factors, "IMMUNOSUPPR_at_covid"),
           c("CV_at_covid", "COVCANCER_at_covid", "COVCOPD_at_covid", "COVHIV_at_covid", "COVCKD_at_covid",
             "COVDIAB_at_covid", "COVOBES_at_covid", "COVSICKLE_at_covid", "at_risk_at_covid",
             "immunosuppressants_at_covid"))
  
  D4_population_c <- merge(D4_population_c_no_risk, D3_population_c_cov_ALL, all.x = T, by="person_id")
  
  tempname<-paste0("D4_population_c",suffix[[subpop]])
  assign(tempname,D4_population_c)
  save(tempname, file = paste0(diroutput, tempname,".RData"),list=tempname)
  rm(list=tempname)
}
