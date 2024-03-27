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
  #
  #                        PREPARATION for MATCHING APPROACH:
  #
  
  # 1. dataset with person characteristics (one row per person):
    load(paste0(diroutput, "D3_study_population_cohort", suffix[[subpop]], ".RData"))  #  "person_id","sex","birthdate","age_at_study_entry","study_entry_date","study_exit_date","datasource","date_vax1","date_of_death"
    D3_study_population_cohort <- D3_study_population_cohort[,!(names(D3_study_population_cohort) %in% c(paste0("date_vax",2:4),paste0("type_vax",2:4)) )]
    D3_study_population_cohort <- D3_study_population_cohort[,!(names(D3_study_population_cohort) %in% 
                                                                  c("covid19_date","gout_date","myocarditis_date","otitis_externa_date","pericarditis_date","valvular_heart_disease_date","gout","myocarditis","pericarditis","otitis_externa","valvular_heart_disease") )]
    D3_study_population_cohort$death_date <- D3_study_population_cohort$date_of_death
    
    # D3_study_population_cohort <- D3_study_population_cohort[sample.int(nrow(D3_study_population_cohort),  nrow(D3_study_population_cohort)/100  ),]
    
    ids_study_all_matching <- D3_study_population_cohort[!duplicated(D3_study_population_cohort[,id_original]),id_original]
    save(ids_study_all_matching,file=paste0(dirtemp_cohort,"ids_study_all_matching", suffix[[subpop]], ".RData"))
    
    D3_study_population_cohort[,id] <- match(D3_study_population_cohort[,id_original],ids_study_all_matching)
    
    # D3_study_population_cohort <- D3_study_population_cohort[ !is.na(D3_study_population_cohort[,id]), ]
    
    data_vax_small <- D3_study_population_cohort[,c(id,"study_entry_date","study_exit_date")]
    data_vax_small$start_date <- data_vax_small$study_entry_date
    data_vax_small$stop_date  <- data_vax_small$study_exit_date
    
    
    # 2.  create periods for vax-independent variables
    
    diseases <- c("IMMUNOSUP","NUMBER_CONDITIONS")
    #diseases <- c("COVID","CANCER","CKD","CVD","DIAB","HIV","IMMUNOSUP","OBES","PULMON","SICKLE","NUMBER_CONDITIONS")
    for(iziekte in diseases){
      if("data_tmp" %in% ls()) rm(data_tmp)
      load(paste0(dirtemp, "TD/D3_TD_variable_",iziekte, suffix[[subpop]], ".RData"))  #   "person_id","value_of_variable","date"
      data_tmp <- get(paste0("D3_TD_variable_",iziekte)); rm(list=paste0("D3_TD_variable_",iziekte))
      
      data_tmp[,"date"] <- as.Date(data_tmp[,"date"] )
      data_tmp <- id_add(data_tmp,ids_study=ids_study_all_matching, id=id_original, id_new=id,lprint=T)
      data_tmp <- data_tmp[order(data_tmp[,id],as.numeric(data_tmp$date)),]
      
      catt(paste("\n",iziekte,":\n\nDelete double rows:  from",(tmp <- nrow(data_tmp)), "rows to " ))
      data_tmp <- unique(data_tmp)
      catt(paste(nrow(data_tmp),"rows   ==>   ", tmp[1]-nrow(data_tmp), "rows deleted.\n\n"))
      
      catt(paste0("\nFrequency for 'value of variable' from dataset 'D3_TD_variable_",iziekte,"':\n"))
      printt(table1(data_tmp[,"value_of_variable"])); catt("\n")
      
      if(sum(is.na(data_tmp[,"date"]))>0) {
        warning(paste0(sum(is.na(data_tmp[,"date"]))," rows with missing 'date' deleted from 'D3_TD_variable_",iziekte,"'!"))
        data_tmp <- data_tmp[!is.na(data_tmp[,"date"]),]
      }
      
      
      names(data_tmp)[names(data_tmp)=="date"]               <- paste0(tolower(iziekte),"_start_date")
      names(data_tmp)[names(data_tmp)=="value_of_variable"]  <- tolower(iziekte)
      
      if(tolower(iziekte) == "number_conditions"){ 
        #browser()
        data_tmp$number_conditions_4_start_date <- data_tmp$number_conditions_start_date
        data_tmp$number_conditions_4            <- data_tmp$number_conditions
        data_tmp$number_conditions_4[data_tmp$number_conditions_4>4 & !is.na(data_tmp$number_conditions)] <- 4
        iziekte <- "number_conditions_4"
      }
      
      data_tmp[,paste0(tolower(iziekte),"_stop_date")] <- c(data_tmp[-1,paste0(tolower(iziekte),"_start_date")]-1,NA)
      data_tmp[ c( data_tmp[-nrow(data_tmp),id] != data_tmp[-1,id], F), paste0(tolower(iziekte),"_stop_date") ] <- NA
      data_vax_small <- add_periods(data_vax_small, data_new=data_tmp, start_state_date=paste0(tolower(iziekte),"_start_date"), stop_state_date=paste0(tolower(iziekte),"_stop_date"), state_var=iziekte, id=id, start_date="start_date", stop_date="stop_date", lsort=T)
      gc()
    }
    rm(data_tmp)
    
 

  # merge: 1-2 min
  dim(D3_study_population_cohort);dim(data_vax_small) # 15214165 x 11 ; 20571878 x 15
  data_vax <- merge.data.frame(D3_study_population_cohort[,!(names(D3_study_population_cohort) %in% c("study_entry_date","study_exit_date"))],data_vax_small, by=id,all.x=T)
  dim(data_vax) # 20571878 x 23

  rm(data_vax_small)
  rm(D3_study_population_cohort)  
 
  
  
  
  
  
  
  ##############################
  #
  #         VAX
  #
  # 2. dataset with vaccines (multiple rows per person): (till merge: 1-2min)
  load(paste0(dirtemp, "D3_vaccines_curated",        suffix[[subpop]], ".RData"))  #  "person_id","dose_curated","manufacturer_curated","date_curated" 
  names(D3_vaccines_curated)[names(D3_vaccines_curated)=="manufacturer_curated"] <- "vax_brand"
  names(D3_vaccines_curated)[names(D3_vaccines_curated)=="date_curated"]         <- "vax_date"
  
  # summary: from the original dataset:
  tb0_12   <- table1(D3_vaccines_curated[,c("dose_curated","vax_brand")])
  tb0_1    <- table1(D3_vaccines_curated[,c("dose_curated")])
  tb0_2    <- table1(D3_vaccines_curated[,c("vax_brand")])
  tb0_2v1  <- table1(D3_vaccines_curated[D3_vaccines_curated$dose_curated==1,c("dose_curated","vax_brand")])
  
  
  D3_vaccines_curated <- id_add(D3_vaccines_curated,ids_study=ids_study_all_matching, id=id_original, id_new=id,lprint=T)
  
  ##########
  # summary: in study:
  tb1_12   <- table1(D3_vaccines_curated[,c("dose_curated","vax_brand")])
  tb1_1    <- table1(D3_vaccines_curated[,c("dose_curated")])
  tb1_2    <- table1(D3_vaccines_curated[,c("vax_brand")])
  tb1_2v1  <- table1(D3_vaccines_curated[D3_vaccines_curated$dose_curated==1,c("dose_curated","vax_brand")])
  
  
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
  
  tb <- table1(D3_vaccines_curated$vax_all[D3_vaccines_curated$vax_n==1])
  tb <- tb[order(tb[,"n"],decreasing = T),]
  tb[,"cum_n"] <- cumsum(tb[,"n"]); tb[,"cum_percent"] <- cumsum(tb[,"percent"])
  #print(tb)
  vax_combi_tab <- tb
  
  tb <- table1(gsub(" ==> >90d novax","",D3_vaccines_curated$vax_all[D3_vaccines_curated$vax_n==1]))
  tb <- tb[order(tb[,"n"],decreasing = T),]
  tb[,"cum_n"] <- cumsum(tb[,"n"]); tb[,"cum_percent"] <- cumsum(tb[,"percent"])
  #print(tb)
  
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
  
  
  # one row per person: vax1 or no vax:
  tb2_12   <- table1(data_vax[,c("vax_n","vax_brand")])
  tb2_1    <- table1(data_vax[,c("vax_n")])
  tb2_2    <- table1(data_vax[,c("vax_brand")])
  tb2_2v1  <- table1(data_vax[data_vax$vax_n==1,c("vax_n","vax_brand")])
  
  
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
  
  
  
  
  # one row per person: vax1 or no vax:
  tb3_12   <- table1(data_vax[,c("vax_n","vax_brand")])
  tb3_1    <- table1(data_vax[,c("vax_n")])
  tb3_2    <- table1(data_vax[,c("vax_brand")])
  tb3_2v1  <- table1(data_vax[data_vax$vax_n==1,c("vax_n","vax_brand")])
  
  tabb12    <- tabb( list( in_dataset=tb0_12,     in_study=tb1_12 ,    vax1=tb2_12  , vax1_after_check=tb3_12    ))
  tabb1     <- tabb( list( in_dataset=tb0_1 ,     in_study=tb1_1  ,    vax1=tb2_1   , vax1_after_check=tb3_1     ))
  tabb2     <- tabb( list( in_dataset=tb0_2 ,     in_study=tb1_2  ,    vax1=tb2_2   , vax1_after_check=tb3_2     ))
  tabb2_v1  <- tabb( list( in_dataset=tb0_2v1 ,   in_study=tb1_2v1  ,  vax1=tb2_2v1 , vax1_after_check=tb3_2v1   ))
  
  
  
  
  # save for report:
  vax_tables_list1 <- list( vax_n_brand=tabb12, vax_n=tabb1, vax_brand=tabb2, vax1_brand=tabb2_v1 )
  vax_tables_list2 <- vax_combi_tab
  
  
  save(list=c("vax_tables_list1","vax_tables_list2"), file=paste(dirtemp_cohort,"vax_tables_list", suffix[[subpop]], ".RData"))
  
  
  
  
  # ids_study_all_matching <- ids_study_all_matching[sample.int(length(ids_study_all_matching), 1000000)]
  # data_vax <- data_vax[data_vax[,id_original] %in% ids_study_all_matching,]
  
  
  save(ids_study_all_matching, file=paste0(dirtemp_cohort,"ids_study_all_matching", suffix[[subpop]], ".RData"))
  save(data_vax, file=paste(dirtemp_cohort,"data_vax_tmp", suffix[[subpop]], ".RData"))
  
  
  
  
  
 
  
  
  
  
  ######################################################### 
  # 3. EVENTS
  
  #  3a. dataset with events (multiple rows per person):
  
  # load( file=paste(dirtemp_cohort,"data_vax_tmp", suffix[[subpop]], ".RData"))
  # load(file=paste0(dirtemp_cohort,"ids_study_all_matching", suffix[[subpop]], ".RData"))
  
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
                                        methods = c("last_before","first_after"), ids_study = ids_study_all_matching, save_all_events=T )     
    
    all_events <- add_td_cov_res$all_events
    data_vax <- add_td_cov_res$data
    rm(add_td_cov_res)
    
    data_vax[is.na(data_vax[,paste0(event_name, "_first_after")]), paste0(event_name, "_first_after")] <- 0
    data_vax[is.na(data_vax[,paste0(event_name, "_last_before")]), paste0(event_name, "_last_before")] <- 0
    
    hist_res <-  hist_vars( list(  all_events      = all_events[ all_events[,event_name]==1 , paste0(event_name,"_date")],  
                                   first_after_vax = data_vax[  (!is.na(data_vax[,paste0(event_name,"_date_first_after")]) & data_vax[,paste0(event_name,"_first_after")]==1 ) , paste0(event_name,"_date_first_after") ], 
                                   last_before_vax = data_vax[  (!is.na(data_vax[,paste0(event_name,"_date_last_before")]) & data_vax[,paste0(event_name,"_last_before")]==1 ) , paste0(event_name,"_date_last_before") ] ),
                            periods=c(7,30.4), tit1=event_name, lplot=F )
    
    # save for report:
    assign(paste0("hist_",event_name,"_data_prep"),hist_res)
    save(list=paste0("hist_",event_name,"_data_prep"), file=paste0(sdr0, "hist_",event_name,"_data_prep", suffix[[subpop]], ".RData"))
    
    load(file=paste0(dirtemp,"events/",dataset_name,suffix[[subpop]],".RData") )
    assign(paste0("event_dataset_",event_name),  id_add(get(dataset_name), ids_study=ids_study_all_matching, id=id_original, id_new=id,lprint=T))
    save(list=paste0("event_dataset_",event_name), file=paste0(dirtemp_cohort, "event_dataset_", event_name, suffix[[subpop]], ".RData"))
    
    
  } 
  
  
  ############################
  #
  #   4.
  
  # create minimal datasets with covid and event info:
  # covid dataset:
  data_vax <- add_td_covariate( data_vax, "COVID", methods = c("last_before","first_after"), na_value=0, 
                                name=paste0("D3_TD_variable_COVID",suffix[[subpop]]), dir=paste0(dirtemp,"td/") )
  
  load(file=paste0(dirtemp,"td/D3_TD_variable_COVID",suffix[[subpop]],".RData") )
  covid_dataset <-  id_add(D3_TD_variable_COVID, ids_study=ids_study_all_matching, id=id_original, id_new=id,lprint=T)
  save(covid_dataset, file=paste0(dirtemp_cohort, "covid_dataset", suffix[[subpop]], ".RData"))
  
 
  
  
  #
  #####################################################################################################
  
  
  # calculate the 'days'-variables:
  names(data_vax)[names(data_vax)=="date_vax1"] <- "vax1_date" 
  for(idate_vars in substring(names(data_vax),6)[substring(names(data_vax),1,5)=="date_"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( data_vax[,paste0("date_",idate_vars)], as.Date("2020-08-31"), units="days"))
  for(idate_vars in substring(names(data_vax),1,nchar(names(data_vax))-5)[substring(names(data_vax),nchar(names(data_vax))-4,nchar(names(data_vax)))=="_date"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( as.Date(data_vax[,paste0(idate_vars,"_date")]), as.Date("2020-08-31"), units="days"))
  names(data_vax)[names(data_vax)=="of_death_days"] <- "death_days" 
  
  
  
  gc()
  
  data_vax_for_matching <- data_vax  
  save(data_vax_for_matching, file=paste0(dirtemp_cohort, "data_vax_for_matching", suffix[[subpop]], ".RData"))
  
  #
  #    end of preparation for matching
  #
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  

}  # end of 'subpop' loop




