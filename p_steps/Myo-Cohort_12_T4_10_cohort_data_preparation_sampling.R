# Program Information  ----------------------------------------------------
#
# Program:      Myo-cohort-01-T5_20_.R # in the 'steps' directory 
# Author:       Svetlana Belitser
# Description:  calls functions to run cohort analyses on a specified dataset 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival" 
#               input:   g_intermediate/data_vax_SCRI.RData
#               output:  g_output/scri/*  
#
#               parameters: in 08_cohort_parameters.R 
#  
#               function: p_macro/cohort_tools.R                    
#
#

# Housekeeping  -----------------------------------------------------------
# install and load packages

if(!any(ls()=="thisdir"))   thisdir   <- getwd()
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"/g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"/g_output/")

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp,   "cohort_causal")),            showWarnings = FALSE, recursive = TRUE)
#dir.create(file.path(paste0(thisdir,   "/log_files/cohort_causal")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "cohort_causal")),            showWarnings = FALSE, recursive = TRUE)


for (subpop in subpopulations_non_empty) {
  
  # cohort output directory for export:  
  sdr0 <- paste0(direxpsubpop[[subpop]], "cohort_causal/")
  dir.create(file.path(sdr0), showWarnings = FALSE, recursive = TRUE)
  
  sdr0_flowchart <- paste0(direxpsubpop[[subpop]], "cohort_causal/flowchart_for_vax/")
  dir.create(file.path(sdr0_flowchart), showWarnings = FALSE, recursive = TRUE)
  
  
  # cohort output_directory for models NOT for export:
  sdr_models0 <- paste0(diroutput, "cohort_causal/")
  if(length(subpopulations_non_empty)>1) sdr_models0 <- paste0(diroutput, "cohort_causal/",subpop,"/")
  dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
  
  # cohort intermediate directory: NOT for export:
  dirtemp_cohort <- paste0(dirtemp, "cohort_causal/")
  dir.create(dirtemp_cohort, showWarnings = FALSE, recursive = TRUE)
  

  catt(paste0('\n\t"',subpop,'":\n\n'))
  

  
  ############################################################################################################################################
  ############################################################################################################################################
  ############################################################################################################################################
  ############################################################################################################################################
  #
  #                              PREPARATION for SAMPLING APPROACH:
  #
  
  #########################
  # 1. dataset with person characteristics (one row per person): (1-2min)
  load(paste0(diroutput, "D3_study_population_cohort", suffix[[subpop]], ".RData"))  #  "person_id","sex","birthdate","age_at_study_entry","study_entry_date","study_exit_date","datasource","date_vax1","date_of_death"
  # D3_study_population_cohort <- D3_study_population_cohort[sample.int(nrow(D3_study_population_cohort),nrow(D3_study_population_cohort)/1000),]
  D3_study_population_cohort <- D3_study_population_cohort[,!(names(D3_study_population_cohort) %in% c(paste0("date_vax",2:4),paste0("type_vax",2:4)) )]
  D3_study_population_cohort <- D3_study_population_cohort[,!(names(D3_study_population_cohort) %in% 
                                                                c("covid19_date","gout_date","myocarditis_date","otitis_externa_date","pericarditis_date","valvular_heart_disease_date","gout","myocarditis","pericarditis","otitis_externa","valvular_heart_disease") )]
  D3_study_population_cohort$death_date <- D3_study_population_cohort$date_of_death
  
  ids_study_all_sampling <- D3_study_population_cohort[!duplicated(D3_study_population_cohort[,id_original]),id_original]
  D3_study_population_cohort[,id] <- match(D3_study_population_cohort[,id_original],ids_study_all_sampling)
  
  
  data_vax <- D3_study_population_cohort
  
  
  ##############################
  #
  #         VAX
  #
  # 2. dataset with vaccines (multiple rows per person): (till merge: 1-2min)
  load(paste0(dirtemp, "D3_vaccines_curated",        suffix[[subpop]], ".RData"))  #  "person_id","dose_curated","manufacturer_curated","date_curated" 
  names(D3_vaccines_curated)[names(D3_vaccines_curated)=="manufacturer_curated"] <- "vax_brand"
  names(D3_vaccines_curated)[names(D3_vaccines_curated)=="date_curated"]         <- "vax_date"
 
  
  D3_vaccines_curated <- id_add(D3_vaccines_curated,ids_study=ids_study_all_sampling, id=id_original, id_new=id,lprint=T)
 
  
  D3_vaccines_curated <- D3_vaccines_curated[order(D3_vaccines_curated[,id],as.numeric(D3_vaccines_curated$vax_date)),]
  
  D3_vaccines_curated$vax_n <- unlist(tapply(D3_vaccines_curated$dose_curated,D3_vaccines_curated[,id],function(x)1:length(x)))
  
  D3_vaccines_curated$not_same_id <- c( D3_vaccines_curated[-nrow(D3_vaccines_curated),id] != D3_vaccines_curated[-1,id], F)
  D3_vaccines_curated$next_vax_date <- c( D3_vaccines_curated$vax_date[-1], NA)
  D3_vaccines_curated$next_vax_date[ D3_vaccines_curated$not_same_id ] <- NA
  
  D3_vaccines_curated$next_vax_brand <- c( D3_vaccines_curated$vax_brand[-1], NA)
  D3_vaccines_curated$next_vax_brand[ D3_vaccines_curated$not_same_id] <- NA
  
  ######################
  # find end of the same brand:
  cond <- !is.na(D3_vaccines_curated$vax_date) & !is.na(D3_vaccines_curated$next_vax_date) & 
    D3_vaccines_curated$vax_brand != D3_vaccines_curated$next_vax_brand & !D3_vaccines_curated$not_same_id
  D3_vaccines_curated$not_same_vax_date[cond] <- D3_vaccines_curated$next_vax_date[cond]
  D3_vaccines_curated$not_same_vax_date  <- as.Date(D3_vaccines_curated$not_same_vax_date,origin="1970-01-01")
  
  tmp <- D3_vaccines_curated[,c(id,"not_same_vax_date")]
  tmp <- tmp[!is.na(tmp$not_same_vax_date),]
  tmp <- tmp[order(tmp[,id],tmp$not_same_vax_date),]
  tmp <- tmp[!duplicated(tmp[,id]),]
  names(tmp)[names(tmp)=="not_same_vax_date"] <- "end_same_vax_date"
  dim(D3_vaccines_curated)
  D3_vaccines_curated <- merge.data.frame(D3_vaccines_curated,tmp,all.x=T)
  dim(D3_vaccines_curated)
  
  
  ###############################
  # find all combinations in the dataset:
  tmp <- rep("",nrow(D3_vaccines_curated))
  cond <- !is.na(D3_vaccines_curated$vax_date) & !is.na(D3_vaccines_curated$next_vax_date) & !D3_vaccines_curated$not_same_id
  cond1 <- cond & (D3_vaccines_curated$next_vax_date - D3_vaccines_curated$vax_date)>90
  tmp[cond1] <- paste0( tmp[cond1], " ==> >90d novax" )
  tmp[cond] <- paste0( tmp[cond], " ==> ", c(D3_vaccines_curated$vax_brand[-1],"")[cond] )
  
  D3_vaccines_curated$vax_all <- substring(D3_vaccines_curated$vax_brand,1,5)
  i <- 1
  while(T){ 
    cond1 <- cond2 <- c(D3_vaccines_curated[-(1:i),id_original]==D3_vaccines_curated[-(nrow(D3_vaccines_curated)-0:(i-1)),id_original],rep(F,i))
    if(i>1) cond2 <- c(rep(F,i-1),cond1[-(nrow(D3_vaccines_curated)-0:(i-2))])
    if(!any(cond1)) break
    D3_vaccines_curated$vax_all[cond1] <- paste0( D3_vaccines_curated$vax_all[cond1], tmp[cond2] )
    i <- i+1
  }
  

  vax_combi_tab <- list(with_90d_interval=vax_combi_tab, without_interval=tb)
  #
  ################
  
  D3_vaccines_curated[,c("not_same_id","not_same_vax_date","next_vax_date","next_vax_brand")] <- NULL
  
  
  # only first vax or no vax ( One row per person.)
  D3_vaccines_curated <- D3_vaccines_curated[order(D3_vaccines_curated[,id],D3_vaccines_curated$vax_date),]
  D3_vaccines_curated <- D3_vaccines_curated[!duplicated(D3_vaccines_curated[,id]),]
  
  
  
  
  gc()
  
  dim(data_vax);dim(D3_vaccines_curated)  
  data_vax <- merge.data.frame(data_vax,D3_vaccines_curated[!duplicated(D3_vaccines_curated[,id]),names(D3_vaccines_curated)!=id_original], by=id,all.x=T)
  dim(data_vax) 
  
  
  data_vax$vax_n[is.na(data_vax$vax_n)] <- 0
  
 
  # calculate the 'days'-variables:
  names(data_vax)[names(data_vax)=="date_vax1"] <- "vax1_date" 
  for(idate_vars in substring(names(data_vax),6)[substring(names(data_vax),1,5)=="date_"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( data_vax[,paste0("date_",idate_vars)], as.Date("2020-08-31"), units="days"))
  for(idate_vars in substring(names(data_vax),1,nchar(names(data_vax))-5)[substring(names(data_vax),nchar(names(data_vax))-4,nchar(names(data_vax)))=="_date"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( as.Date(data_vax[,paste0(idate_vars,"_date")]), as.Date("2020-08-31"), units="days"))
  names(data_vax)[names(data_vax)=="of_death_days"] <- "death_days" 
  
  
  #############################################################
  # check:
  tb<-table((cond<-data_vax$study_entry_days < data_vax$study_exit_days ))
  
  if(length(tb)>1){
    warningg(paste("There are ",tb["FALSE"],"rows with 'study_entry_date' >= 'study_exit_date' !"))
    catt('\n"study_entry_days" < "study_exit_days":\n')
    printt(table(tb))
    catt(paste0('\n"study_entry_days" == "study_exit_days": ', sum(data_vax$study_entry_days == data_vax$study_exit_days,na.rm=T ),' rows\n\n'))
    data_vax_excluded_wrong_study_period <- data_vax[!cond,]
    save(data_vax_excluded_wrong_study_period,file=paste0(sdr_models0,"excluded_rows_wrong_study_period", suffix[[subpop]], ".RData"))
    sink(paste0(sdr_models0,"excluded_rows_wrong_study_period.txt")); old_sink = options (width=300, max.print=99999 );printt(data_vax_excluded_wrong_study_period);options(old_sink);sink()
    rm(data_vax_excluded_wrong_study_period)
    data_vax <- data_vax[cond,]
  }
  
  cond <- !is.na(data_vax$study_entry_days) & ( is.na(data_vax$vax_days) | ( !is.na(data_vax$vax_days) & data_vax$study_entry_days <= data_vax$vax_days  & data_vax$vax_days <= data_vax$study_exit_days ) )
  if(any(!cond)){
    warningg(paste( sum(!cond), "rows with the vaccination date not in ['study_entry_days';'study_exit_days']"))
    data_vax_excluded_wrong_vax_date <- data_vax[!cond,]
    save(data_vax_excluded_wrong_vax_date,file=paste0(sdr_models0,"excluded_rows_wrong_vax_date", suffix[[subpop]], ".RData"))
    sink(paste0(sdr_models0,"excluded_rows_wrong_vax_date.txt")); old_sink = options (width=300, max.print=99999 );printt(data_vax_excluded_wrong_vax_date);options(old_sink);sink()
    rm(data_vax_excluded_wrong_vax_date)
    data_vax <- data_vax[cond,]
  }
  
  cond <- !is.na(data_vax$study_entry_days) & ( is.na(data_vax$death_days) | ( !is.na(data_vax$death_days) & data_vax$study_entry_days <= data_vax$death_days & data_vax$death_days <= data_vax$study_exit_days ) )
  if(any(!cond)){
    warningg(paste( sum(!cond), "rows with the date of death not in ['study_entry_days';'study_exit_days']"))
    data_vax_excluded_wrong_death_date <- data_vax[!cond,]
    save(data_vax_excluded_wrong_death_date,file=paste0(sdr_models0,"excluded_rows_wrong_death_date", suffix[[subpop]], ".RData"))
    sink(paste0(sdr_models0,"excluded_rows_wrong_death_date.txt")); old_sink = options (width=300, max.print=99999 );printt(data_vax_excluded_wrong_death_date);options(old_sink);sink()
    rm(data_vax_excluded_wrong_death_date)
    data_vax <- data_vax[cond,]
  }
  #    end of check
  #################################################################################
  
  
  
  
  save(ids_study_all_sampling, file=paste0(dirtemp_cohort,"ids_study_all_sampling", suffix[[subpop]], ".RData"))
  save(data_vax, file=paste(dirtemp_cohort,"data_vax_tmp_sampling", suffix[[subpop]], ".RData"))
 
   
   
   
  ######################################################### 
  # 3. EVENTS
  
  #  3a. dataset with myocarditis (multiple rows per person):
  
  # load( file=paste(dirtemp_cohort,"data_vax_tmp", suffix[[subpop]], ".RData"))
  # load(file=paste0(dirtemp_cohort,"ids_study_all_sampling", suffix[[subpop]], ".RData"))
  
  # from parameters:
  #ae_events <-  c("myocarditis", "pericarditis", "myopericarditis", "otitis_externa", "valvular_heart_disease")
  
  
  #par(mfcol=c(4,2))
  for(event_name in ae_events){
    
    if(event_name =="myocarditis"           ) dataset_name = "D3_events_C_MYOCARD_AESI_simple"
    if(event_name =="pericarditis"          ) dataset_name = "D3_events_C_PERICARD_AESI_simple"
    if(event_name =="myopericarditis"       ) next #dataset_name = "???"
    if(event_name =="otitis_externa"        ) dataset_name = "D3_events_SO_OTITISEXT_COV_simple"
    if(event_name =="valvular_heart_disease") dataset_name = "D3_events_C_VALVULAR_COV_simple"
    
    add_td_cov_res <- add_td_covariate( data_vax, event_name , name=paste0(dataset_name,suffix[[subpop]]), dir=paste0(dirtemp,"events/"),  create_cov_var=T,
                                        methods = c("last_before","first_after"), ids_study = ids_study_all_sampling, save_all_events=T )     
    
    all_events <- add_td_cov_res$all_events
    data_vax <- add_td_cov_res$data
    rm(add_td_cov_res)
    
    data_vax[is.na(data_vax[,paste0(event_name, "_first_after")]), paste0(event_name, "_first_after")] <- 0
    data_vax[is.na(data_vax[,paste0(event_name, "_last_before")]), paste0(event_name, "_last_before")] <- 0
    
 
  }     
    
    # 4.  COVID variable
    
     # create minimal datasets with covid and event info:
    # covid dataset:
    data_vax <- add_td_covariate( data_vax, "COVID", methods = c("last_before","first_after"), na_value=0, ids_study = ids_study_all_sampling,
                                  name=paste0("D3_TD_variable_COVID",suffix[[subpop]]), dir=paste0(dirtemp,"td/") )
 
  
  
  # calculate the 'days'-variables:
  names(data_vax)[names(data_vax)=="date_vax1"] <- "vax1_date" 
  for(idate_vars in substring(names(data_vax),6)[substring(names(data_vax),1,5)=="date_"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( data_vax[,paste0("date_",idate_vars)], as.Date("2020-08-31"), units="days"))
  for(idate_vars in substring(names(data_vax),1,nchar(names(data_vax))-5)[substring(names(data_vax),nchar(names(data_vax))-4,nchar(names(data_vax)))=="_date"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( as.Date(data_vax[,paste0(idate_vars,"_date")]), as.Date("2020-08-31"), units="days"))
  names(data_vax)[names(data_vax)=="of_death_days"] <- "death_days" 
  
  
  
  gc()
  
  data_vax_for_sampling <- data_vax  
  save(data_vax_for_sampling, file=paste0(dirtemp_cohort, "data_vax_for_sampling", suffix[[subpop]], ".RData"))
  
  
  
  #
  #    end of preparation for sampling
  #
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  
  
  
  
}  # end of 'subpop' loop




