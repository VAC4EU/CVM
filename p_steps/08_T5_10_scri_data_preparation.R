# Program Information  ----------------------------------------------------
#
# Program:      step_08_..._scri_data_preparation.R 
# Author:       Svetlana Belitser
# Description:  calls functions which runs an SCRI on specified datasets 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival", "MASS", "qpdf", "metafor""
#               input:   D3_study_population_SCRI  
#               output:  data_vax_SCRI 
#                        excluded_rows.RData
#                        excluded_rows_2.RData
#
#               parameters: in 08_parameters_scri_inputs.R 
#  
#

# Housekeeping  -----------------------------------------------------------
# install and load packages

if(!any(ls()=="thisdir"))   thisdir   <- getwd()
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"/g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"/g_output/")

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp,   "scri")),            showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,   "/log_files/scri")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "scri")),            showWarnings = FALSE, recursive = TRUE)


for (subpop in subpopulations_non_empty) {
  
  # scri output_directory for export:  
  sdr0 <- paste0(direxpsubpop[[subpop]], "scri/")
  dir.create(file.path(sdr0), showWarnings = FALSE, recursive = TRUE)
  
  # scri output_directory for models not for export:
  sdr_models0 <- paste0(diroutput, "scri/")
  if(length(subpopulations_non_empty)>1) sdr_models0 <- paste0(diroutput, "scri/",subpop,"/")
  dir.create(sdr_models0, showWarnings = FALSE, recursive = TRUE)
  
  
  cat(paste0('\n\t"',subpop,'":\n\n'))
  
  
  memory.limit(10^13)
  #detach("package:data.table", unload = TRUE)
  
 #  max_number_doses <- 4
  
  
  
  # Import Data -------------------------------------------------------------
  
  # Load the D3_study_population_SCRI
  load(paste0(dirtemp, "D3_study_population_SCRI", suffix[[subpop]], ".RData"))
  scri_input <- as.data.frame(get(paste0("D3_study_population_SCRI", suffix[[subpop]])))
  rm(list = paste0("D3_study_population_SCRI", suffix[[subpop]]))
  
  # scri_input <- D3_study_population_SCRI  
  #################
  # create dataset 'data_vax' (with multiple rows per person) from 'scri_input' (with one row per person)
  #
  
  # change the id variable:
  scri_input$pat_n <- 1:nrow(scri_input)
  scri_input[,id] <- NULL
  id <- "pat_n"
  
  scri_input[,names(scri_input) %in% ae_events] <- NULL
  
  if(any(c("date_vax1","vax1_date") %in% names(scri_input))){
    date_var <- "date_vax"; type_var <- "type_vax"
    if(any(c("vax1_date") %in% names(scri_input))) date_var <- "vax_date"
    if(any(c("vax1_type") %in% names(scri_input))) type_var <- "vax_type"
    scri_n_vax_var <- as.integer(substring(names(scri_input),9)[substring(names(scri_input),1,8)==date_var][1:max_number_doses]) # ==> 1,2,3  (from "data_vax1", "date_vax2", ...)
    scri_n_vax_var <- scri_n_vax_var[!is.na(scri_n_vax_var)]
    
    for(ivax in scri_n_vax_var){
      print(Sys.time())  
      data_tmp <- scri_input[, !substring(names(scri_input),1,8) %in% c(type_var,date_var) | 
                               ( substring(names(scri_input),1,8) %in% c(type_var,date_var) & names(scri_input) %in% paste0(c(type_var,date_var),ivax) ) ]
      names(data_tmp)[names(data_tmp)==paste0(date_var,ivax)] <- "vax_date" 
      names(data_tmp)[names(data_tmp)==paste0(type_var,ivax)] <- "vax_brand"
      data_tmp$vax_n <- ivax
      data_tmp <- data_tmp[!is.na(data_tmp$vax_date),c( names(data_tmp)[!names(data_tmp) %in% c("vax_brand","vax_date")], "vax_brand", "vax_date" )]
      
      if(ivax==1) data_vax <- data_tmp
      else        data_vax <- rbind.data.frame(data_vax, data_tmp)
      gc()
    }
    print(Sys.time())
    data_tmp <- data_tmp[F,]
    if("vax1_date" %in% names(scri_input)) data_tmp <- scri_input[is.na(scri_input[,"vax1_date"]),]
    if("date_vax1" %in% names(scri_input)) data_tmp <- scri_input[is.na(scri_input[,"date_vax1"]),]
    data_tmp$vax_n <- 0; data_tmp$vax_date <- rep(NA,nrow(data_tmp)); data_tmp$vax_brand <- rep(NA,nrow(data_tmp))
    if(nrow(data_tmp)>0) {gc(); data_vax <- rbind.data.frame(data_vax, data_tmp[,names(data_vax)]) }
    data_vax <- data_vax[!(data_vax$vax_n>0 & is.na(data_vax$vax_date)),]
  } else {
    data_vax <- scri_input
    names(data_vax)[names(data_vax)  ==    "date_vax"            ]  <- "vax_date" 
    names(data_vax)[names(data_vax) %in% c("type_vax","vax_type")]  <- "vax_brand" 
  }
  rm(scri_input)
  gc()
  #
  #######
  Sys.time()  
  
  data_vax[,"vax_days"]  <- as.integer( difftime( data_vax[,"vax_date"], as.Date("2020-08-31"), units="days"))
  # sort per id, vax_days
  data_vax <- data_vax[order(data_vax$pat_n,data_vax[,"vax_days"]),]
  #data_vax <- data_vax[order(data_vax[,id],data_vax[,"vax_days"]),]
  
  # dap:
  if(!any(ls()=="dap")) dap <- ifelse( any(tolower(names(data_vax))=="dap"), data_vax[1,tolower(names(data_vax))=="dap"], "" )
  if(dap=="" & any(tolower(names(data_vax))=="datasource")) dap <- data_vax[1,tolower(names(data_vax))=="datasource"]
  
  # calculate the 'days'-variables:
  for(idate_vars in substring(names(data_vax),6)[substring(names(data_vax),1,5)=="date_"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( data_vax[,paste0("date_",idate_vars)], as.Date("2020-08-31"), units="days"))
  for(idate_vars in substring(names(data_vax),1,nchar(names(data_vax))-5)[substring(names(data_vax),nchar(names(data_vax))-4,nchar(names(data_vax)))=="_date"])
    data_vax[,paste0(idate_vars,"_days")]  <- as.integer( difftime( as.Date(data_vax[,paste0(idate_vars,"_date")]), as.Date("2020-08-31"), units="days"))
  names(data_vax)[names(data_vax)=="of_death_days"] <- "death_days" 
  
  
  #############   SCRI models ############################
  #
  #
  old_width = options(width=300)
  
  # check:
  tb<-table((cond<-data_vax$study_entry_days < data_vax$study_exit_days ))
  
  if(length(tb)>1){
    warning(paste("There are ",tb["FALSE"],"rows with 'study_entry_date' >= 'study_exit_date' !"))
    cat('\n"study_entry_days" < "study_exit_days":\n')
    print(table(tb))
    cat(paste0('\n"study_entry_days" == "study_exit_days": ', sum(data_vax$study_entry_days == data_vax$study_exit_days,na.rm=T ),' rows\n\n'))
    data_vax_excluded <- data_vax[!cond,]
    save(data_vax_excluded,file=paste0(sdr0,"excluded_rows.RData"))
    sink(paste0(sdr0,"excluded_rows.txt")); old_sink = options (width=300, max.print=99999 );print(data_vax_excluded);options(old_sink);sink()
    rm(data_vax_excluded)
    data_vax <- data_vax[cond,]
  }
  
  cond <- !is.na(data_vax$study_entry_days) & ( is.na(data_vax$vax_days) | ( !is.na(data_vax$vax_days) & data_vax$study_entry_days <= data_vax$vax_days ) )
  if(any(!cond)){
    warning(paste( sum(!cond), "rows with 'study_entry_days' > 'vax_days'"))
    data_vax_excluded_2 <- data_vax[!cond,]
    save(data_vax_excluded_2,file=paste0(sdr0,"excluded_rows_2.RData"))
    sink(paste0(sdr0,"excluded_rows_2.txt")); old_sink = options (width=300, max.print=99999 );print(data_vax_excluded_2);options(old_sink);sink()
    rm(data_vax_excluded_2)
    data_vax <- data_vax[cond,]
  }
  #
  ########################################################
  
  ##########################################
  #
  #   create variables: 
  #   dose number: 'vax_n'      (1,2,...) 
  #                'vax_number' ("dose 1" ,"dose 2", "dose 3", ...)
  #
  gc()
  #### 'vax_n':
  data_vax <- data_vax[order(data_vax[,id],data_vax[,"vax_days"]),]
  ivax <- 1; 
  data_vax[,"vax_n"] <- data_vax[,"vax_days"]
  data_vax[!is.na(data_vax[,"vax_n"]),"vax_n"] <- ivax 
  data_vax[ is.na(data_vax[,"vax_n"]),"vax_n"] <- 0 
  while(T){
    cond_vax_i <- data_vax$vax_n == ivax & !is.na(data_vax$vax_n) 
    cond2<-duplicated(data_vax[cond_vax_i,id])
    if(!any(cond2)) break
    ivax <- ivax + 1
    data_vax[cond_vax_i,"vax_n"][cond2] <- ivax
  }
  #### 'vax_number':
  data_vax$vax_number <- paste0("dose ",data_vax[,"vax_n"]," ")
  
  #### 'dist'
  data_vax$dist <- c(NA,  diff(data_vax[,"vax_days"])  )
  data_vax[ c(F, data_vax[-1,id]!=data_vax[-nrow(data_vax),id] ), "dist"  ] <- NA
  #### 'dist_gt_60'
  data_vax$dist_gt_60 <- c("<=60d"," >60d")[(data_vax$dist > 60) +1 ]
  data_vax[is.na(data_vax$dist_gt_60),"dist_gt_60"] <- ""
  
  #### 'vax_name'
  data_vax[ data_vax$vax_n==1, "vax_name"  ] <- "dose  1.1 " 
  data_vax[ data_vax$vax_n==2, "vax_name"  ] <- "dose  1.2 " 
  data_vax[ data_vax$vax_n==3, "vax_name"  ] <- "booster 1" 
  data_vax[ data_vax$vax_n==4, "vax_name"  ] <- "boost2or3" 
  
  data_vax[ data_vax$vax_n==2 & data_vax$dist>60, "vax_name" ]     <- "booster 1"
  data_vax[ c(F, data_vax[-1,"vax_n"]==3 & data_vax[-nrow(data_vax),"vax_name"]=="booster 1"), "vax_name" ] <- "booster 2" 
  data_vax[ c(F, data_vax[-1,"vax_n"]==4 & data_vax[-nrow(data_vax),"vax_name"]=="booster 2"), "vax_name" ] <- "booster 3" 
  
  #### create 'vax_brand_short':
  data_vax$vax_brand_short <- format(substring(data_vax[,"vax_brand"],1,5))
  
  ####################
  # add information about the first dose in separate variables:
  data_v1 <- data_vax[!duplicated(data_vax[,id]),c(id,"vax_days","vax_number","vax_name","vax_date","vax_brand","vax_brand_short")]  
  names(data_v1)[-1] <- paste0(names(data_v1)[-1],"_v1")
  gc()
  print(Sys.time())
  data_vax <- merge.data.frame(data_vax,data_v1,by=id,all.x=T )
  
  print(Sys.time()  )
  
  #######################################################################
  #
  #       strata analyse with brand for: age30, age30_50, sex, sex_age30 
  #
  
  data_vax$age30     <- paste0("age",as.character(cut(data_vax$age_at_study_entry, c(-1,30,   Inf)))); data_vax$age30[   is.na(data_vax$age_at_study_entry)] <- NA
  data_vax$age30_50  <- paste0("age",as.character(cut(data_vax$age_at_study_entry, c(-1,30,50,Inf)))); data_vax$age30_50[is.na(data_vax$age_at_study_entry)] <- NA
  
  if("o" %in% tolower(data_vax$sex))  data_vax$sex[tolower(data_vax$sex)=="o"]<-NA
  data_vax$sexc      <- paste0("sex:",data_vax$sex );                     data_vax$sexc[     is.na(data_vax$sex)                        ] <- NA
  data_vax$sex_age30 <- paste0("sex:",data_vax$sex, " ", data_vax$age30); data_vax$sex_age30[is.na(data_vax$sex) | is.na(data_vax$age30)] <- NA
  
  
  # Save dataset 'data_vax' in file 'data_vax_SCRI.RData':
  save(data_vax, file=paste0(dirtemp, "data_vax_SCRI", suffix[[subpop]], ".RData"))

  gc()
    
  ##########
  # restore options:
  options(old_width)
  
  
}  # end of 'subpop' loop

