#------------------------------------------------------------------
# CREATE RISK FACTORS

# input: D3_Vaccin_cohort_no_risk, concept set datasets of covariates (plus the six concept sets of the three outcomes CAD, MYOCARD and HF which form the covariate CV): "CV","COVCANCER","COVCOPD","COVHIV","COVCKD","COVDIAB","COVOBES","COVSICKLE"
# output: D3_Vaccin_cohort_covariates

print('CREATE RISK FACTORS (diagnosis only) at date_vax_1')

CV_string<-c("HF_narrow","HF_possible","MYOCARD_narrow","MYOCARD_possible","CAD_narrow","CAD_possible") 

lapply(paste0(dirtemp,CV_string,".RData"),load,.GlobalEnv)


CV<- rbind(HF_narrow, HF_possible,MYOCARD_narrow,MYOCARD_possible,CAD_narrow, CAD_possible, fill = T)
rm(HF_narrow,MYOCARD_narrow,HF_possible,MYOCARD_possible, CAD_narrow, CAD_possible)


COVnames<-c("CV","COVCANCER","COVCOPD","COVHIV","COVCKD","COVDIAB","COVOBES","COVSICKLE")

D3_study_population_covariates <- vector(mode = 'list')
for (subpop in subpopulations_non_empty) {
  print(subpop)
  load(paste0(dirtemp,"D3_Vaccin_cohort_no_risk", suffix[[subpop]],".RData")) 
  
  study_population<-get(paste0("D3_Vaccin_cohort_no_risk", suffix[[subpop]]))
  study_population_covariates <- study_population[, .(person_id, date_vax1)]
  
  
  #file<-"COVCANCER_narrow"
  for (file in COVnames) {
    if ( file!="CV" ){
      
      load(paste0(dirtemp,file,".RData"))
      filecovariate <- get(file)
      
      if (this_datasource_has_subpopulations == TRUE){
        selection <- select_in_subpopulationsEVENTS[[subpop]]
        filecovariate = filecovariate[eval(parse(text = selection)),]
      }
      
      temp<-merge(study_population,filecovariate, all.x = T, by="person_id")[,.(person_id,date_vax1,date)]
      temp<-temp[date>=start_lookback & date<date_vax1,file:=1][is.na(file),file:=0]
      suppressWarnings(temp<-unique(temp[,file1:=max(file),by="person_id"][,.(person_id,file1)]))
      setnames(temp,"file1",paste0(file,"_at_date_vax_1"))
      study_population_covariates<-merge(study_population_covariates,temp,all.x = T, by="person_id")
      
    } else {
      
      temp<-merge(study_population,get(file), all.x = T, by="person_id")[,.(person_id,date_vax1,date)]
      temp<-temp[date>=start_lookback & date<date_vax1,file:=1][is.na(file),file:=0]
      suppressWarnings(temp<-unique(temp[,file1:=max(file),by="person_id"][,.(person_id,file1)]))
      setnames(temp,"file1",paste0(file,"_at_date_vax_1"))
      study_population_covariates<-merge(study_population_covariates,temp,all.x = T, by="person_id")
    }
  }
  
  tempname<-paste0("D3_Vaccin_cohort_covariates",suffix[[subpop]])
  assign(tempname,study_population_covariates)
  save(list=tempname,file=paste0(dirtemp,tempname,".RData"))
  
  rm(list=paste0("D3_Vaccin_cohort_no_risk", suffix[[subpop]]))
  rm(tempname,list=tempname)
}    


rm(CV,COVCANCER,COVCOPD,COVHIV,COVCKD,COVDIAB,COVOBES,COVSICKLE) 
rm(temp,D3_study_population_covariates, study_population, study_population_covariates,filecovariate)


