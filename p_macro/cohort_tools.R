# Program Information  ----------------------------------------------------
#
#  functions for cohort analysis 
#
# Functions:  add_td_covariate
#             sample_controls
#             add_periods
#             hist_distribution
#             hist_vars
#             id_add
#             match_pop
#             get_matched_dataset
#             tabb
#             table1
#             
#
# Author:       Svetlana Belitser 
#               dec 2022 - apr 2023
#
#######################################################################




add_td_covariate <- function( dd,
                              cov_name="DIAB", # g_cov_name, # g_cov_name="if(x>4)x=4"
                              methods = "last_before",  # or c("laste_before","first_after", "diff")   may more than one
                              index_var="vax_date", 
                              diff_value, diff_unit="days",    # used if methods=="diff"
                              na_value,                        # if NA ==> change to na_value
                              ids_study = ids_study_all, 
                              cov_dataset,   # dataset with variables: id, "date","value_of_variable" 
                              name=paste0("D3_TD_variable_", cov_name,suffix[[subpop]]), dir=paste0(dirtemp,"TD/"),
                              date_var = "date", cov_var = "value_of_variable", create_cov_var=F, event_value=1, 
                              id_original="person_id", id="id_n", 
                              lprint=T, lfreq=F, lplot=F,
                              save_all_events = F){   
  
  if(lprint) catt(paste0("\n",cov_name,"     (",Sys.time(),") ","\n"))
  methods <- gsub(" ","_",methods) 
  
  if(missing(cov_dataset)){
    if(!file.exists(paste0(dir,name, ".RData"))) stop(paste0("File '",dir,name, ".RData","' is not found"))
    load(paste0(dir,name, ".RData"))  #   "person_id","value_of_variable","date"
    cov_dataset <- get(name); rm(list=name)
    cov_dataset <- id_add(cov_dataset,ids_study=ids_study, id=id_original, id_new=id,lprint=T)
  }

  
  if(date_var != "date" | !(date_var %in% names(cov_dataset))){
    if(date_var %in% names(cov_dataset)) cov_dataset[,"date"] <- cov_dataset[,date_var]  
    else stop(paste("date var is not found in",name) )
  }
  if(   cov_var %in% names(cov_dataset) & cov_var!="value_of_variable") cov_dataset[,"value_of_variable"] <- cov_dataset[,cov_var ]
  if( !(cov_var %in% names(cov_dataset))) {
      if(create_cov_var) cov_dataset[,"value_of_variable"] <- event_value
      else stop(paste("status var is not found for",cov_name) )
    }

  
  cov_dataset <- cov_dataset[cov_dataset[,id] %in% dd[!duplicated(dd[,id]),id],]
  cov_dataset[,"date"] <- as.numeric(as.Date(cov_dataset[,"date"] ))
  cov_dataset <- cov_dataset[order(cov_dataset[,id],cov_dataset$date),]
  
  if(any(cond <- duplicated(cov_dataset))){
    if(lprint) catt(paste("Delete double rows:  from",(tmp <- nrow(cov_dataset)), "rows to " ))
    cov_dataset <- cov_dataset[!cond,]
    if(lprint) catt(paste(nrow(cov_dataset),"rows   ==>   ", tmp[1]-nrow(cov_dataset), "rows deleted.\n"))
  }
  if(lfreq & lprint){
    catt(paste0("\nFrequency for 'value of variable' from dataset '",cov_name,"':\n"))
    print(table1(cov_dataset[,"value_of_variable"])); catt("\n")
  }
  cov_name <- tolower(cov_name)
  
  if(sum(is.na(cov_dataset[,"date"]))>0) {
    warning(paste0(sum(is.na(cov_dataset[,"date"]))," rows with missing 'date' deleted from '",cov_name,"'!"))
    cov_dataset <- cov_dataset[!is.na(cov_dataset[,"date"]),]
  }
  
  names(cov_dataset)[names(cov_dataset)=="value_of_variable"]  <- cov_name
  if(save_all_events) { 
    all_events <- cov_dataset[,c(id,"date",cov_name)]
    all_events$date <- as.Date(all_events$date,origin="1970-01-01")  # for histogram
    names(all_events)[names(all_events)=="date"] <- paste0(cov_name,"_date")
  }
  
  tmp0 <- merge.data.frame( cov_dataset ,dd[!is.na(dd[,index_var]),], by=id, all=F)
  #tmp <- merge.data.frame( cov_dataset[,c(id,"date",cov_name)] ,dd[, !(names(dd) %in% c("date_tmp",cov_name))], by=id, all.x=T)
  
  for(imethod in methods){
    
    if("last_before" %in% substring(imethod,1,11)) {
      # index_var =""vax_date; "date" - date of event or covariate
      tmp <- tmp0[ as.numeric(tmp0$date) < as.numeric(tmp0[,index_var]) ,  ]
      tmp <- tmp[order(tmp[,id],tmp$date,decreasing=c(F,T)),]
    }
    if("first_after" %in% substring(imethod,1,11)) {
      # index_var =""vax_date; "date" - date of event or covariate
      tmp <- tmp0[ as.numeric(tmp0$date) > as.numeric(tmp0[,index_var]) ,  ]
      tmp <- tmp[order(tmp[,id],tmp$date),]
    }
    names(tmp)[names(tmp)=="date"  ] <- paste0(cov_name,"_date_",imethod)
    names(tmp)[names(tmp)==cov_name] <- paste0(cov_name,     "_",imethod)
    tmp <- tmp[!duplicated(tmp[,id]),]
 
    if(lprint) catt(paste("dd_nrow=",nrow(dd), "cov_data_nrow=",nrow(tmp),"     ", Sys.time() ))
    dd <- merge.data.frame(dd[, !(names(dd) %in% c("date_tmp",cov_name))], tmp[,c(id, paste0(cov_name,c("_date_","_"),imethod) )], by=id, all.x=T)
    if(lprint) catt(paste("  ==> dd_nrow=",nrow(dd), "\n"))
    
    dd[,paste0(cov_name,"_date_",imethod)] <- as.Date(dd[,paste0(cov_name,"_date_",imethod)], origin="1970-01-01")
    
  }
  if(F){ 
    
    names(cov_dataset)[names(cov_dataset)=="date"]  <- "date_tmp"
    new_date_vars <- c()
    
    dd[,index_var] <- as.numeric(dd[,index_var])
    for(imethod in methods){
      dd[,paste0(cov_name,"_date_",imethod)] <- dd[,paste0(cov_name,"_", imethod)] <- rep(NA,nrow(dd))
      new_date_vars <- c(new_date_vars, paste0(cov_name,"_date_",imethod))
    }
    
    while(nrow(cov_dataset)>0){ 
      
      if(lprint) catt(paste("dd_nrow=",nrow(dd), "cov_data_nrow=",nrow(cov_dataset),"     ", Sys.time() ))
      dd <- merge.data.frame(dd[, !(names(dd) %in% c("date_tmp",cov_name))], cov_dataset[!duplicated(cov_dataset[,id]),c(id,"date_tmp",cov_name)], by=id, all.x=T)
      if(lprint) catt(paste("  ==> dd_nrow=",nrow(dd), "\n"))
      
      for(imethod in methods){
        
        if("last_before" %in% substring(imethod,1,11)) {
          cond <-  dd[,"date_tmp"] < dd[,index_var]   & !is.na(dd[,"date_tmp"]) & !is.na(dd[,index_var])
          cond <- cond & ( is.na(dd[,paste0(cov_name,"_date_",imethod)]) | ( !is.na(dd[,paste0(cov_name,"_date_",imethod)]) & dd[,paste0(cov_name,"_date_",imethod)] < dd[,"date_tmp"] )  )
        }
        if("first_after" %in% substring(imethod,1,11)) {
          cond <- dd[,index_var] <=  dd[,"date_tmp"]  & !is.na(dd[,"date_tmp"]) & !is.na(dd[,index_var])
          cond <- cond & ( is.na(dd[,paste0(cov_name,"_date_",imethod)]) | ( !is.na(dd[,paste0(cov_name,"_date_",imethod)]) & dd[,"date_tmp"] < dd[,paste0(cov_name,"_date_",imethod)] )  )
        }
        
        if(substring(imethod,1,11) %in% c("last_before","first_after")){
          dd[cond, paste0(cov_name,"_date_",imethod)] <- dd[cond,"date_tmp"]
          dd[cond, paste0(cov_name,"_",imethod)]     <- dd[cond, cov_name]
        }
      }
      
      dd[,c("date_tmp",cov_name)] <- NULL
      cov_dataset <- cov_dataset[duplicated(cov_dataset[,id]),]
      
    }
  }
  
  if(!missing(na_value)) for(imethod in methods) dd[ is.na(dd[, paste0(cov_name,"_", imethod)]), paste0(cov_name,"_", imethod) ]  <- na_value
  
  dd[,index_var] <- as.Date(dd[,index_var], origin="1970-01-01")
  
  
  if(F){
    par(mfrow=c(2,1)) 
    hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1) )
    par(mfrow=c(2,2))
    hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1,365/12) )
    
    par(mfrow=c(2,1))
    hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1) )
    par(mfrow=c(2,2))
    hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1,365/12) )
    
    
    
    h2<-hist(dd$diab_date_last_before[dd$diab_last_before==1],100,xlab="Time",main=paste0(cov_name," in vaccinated"),freq=T)
    hist(dd$diab_date_last_before[dd$diab_last_before==1],100,xlab="Time",main=paste0(cov_name," in vaccinated"),freq=T,
         ylab=paste0("Freq. till ",max(h2$counts[-1])), ylim=c(0,max(h2$counts[-1])) )
    
    if("case_control" %in% names(dd)){
      h3<-hist(dd$diab_date_last_before[dd$diab_last_before==1 & dd$case_control==1 ],100,xlab="Time",main=paste0(cov_name," in vaccinated"),freq=T)
      hist(dd$diab_date_last_before[dd$diab_last_before==1 & dd$case_control==1 ],100,xlab="Time",main=paste0(cov_name," in vaccinated"),freq=T,
           ylab=paste0("Freq. till ",max(h3$counts[-1])), ylim=c(0,max(h3$counts[-1])) )
      
      h4<-hist(dd$diab_date_last_before[dd$diab_last_before==1 & dd$case_control==1 ],100,xlab="Time",main=paste0(cov_name," in sampled"),freq=T)
      hist(dd$diab_date_last_before[dd$diab_last_before==1 & dd$case_control==1 ],100,xlab="Time",main=paste0(cov_name," in sampled"),freq=T,
           ylab=paste0("Freq. till ",max(h4$counts[-1])), ylim=c(0,max(h4$counts[-1])) )
    }
    
  }
  
  
  if(lprint) print(Sys.time())
  
  if(save_all_events) res <- list(data= dd, all_events = all_events )
  else res <- dd
  
  res
  
} # end function "add_td_covariate"





sample_controls <- function( data, 
                             cases,            # name of logical variable (example, T if vax_n=1 & vax_brand=="Pfizer" ==> sampling only for Pfizer vax1)
                             event_date,
                             cond_before_vax,  # list with elements: list( list( var="covid_date_last_before",       period= 30, unit="days", unit_short="d" ),
                             #                           list( var="myocarditis_date_last_before", period=365, unit="days", unit_short="d" )     )
                             vax_date="vax_date", id="id_n", begin_obs="study_entry_date", end_obs="study_exit_date",
                             ids_study
){
  
  if(missing(cases)){ cases <- "_tmp_cases"; data[,"_tmp_cases"] <- !is.na(data[,vax_date]) }
  cond <- !data[,cases] & !is.na(data[,vax_date]) &  data[,vax_date]  <= data[,end_obs]
  data[ cond, end_obs ] <- data[ cond, vax_date ] - 1

  # select id's with conditions before vax. (example: nocovid_before_vax_days during 30 days, or no myocarditits before vax during 1 year)
  if(!missing(cond_before_vax)){
    if(missing(ids_study)) stop(paste("'ids_study' must be specified for sampling"))
    for(icond in 1:length(cond_before_vax)){
      if( !("unit"           %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$unit           <- "days"
      if( !("unit_short"     %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$unit_short     <- "d"
      if( !("create_cov_var" %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$create_cov_var <- F
      if( !("method"         %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$method         <- "last_before"
      cond_before_vax[[icond]]$var <- paste0(cond_before_vax[[icond]]$cov_name,"_date_",cond_before_vax[[icond]]$method)
      
      if( !(cond_before_vax[[icond]]$var %in% names(data)) )
        data <- add_td_covariate( data, cond_before_vax[[icond]]$cov_name, methods=cond_before_vax[[icond]]$method, 
                                  name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                  create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                  ids_study = ids_study  )    
      
      diff_days <- as.numeric( difftime( data[,vax_date], data[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
      if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      
      data[!cond,end_obs] <- data[!cond,vax_date] - 1
      data[!cond,cases]       <- F 
      
      print(table1(cond))       
      print(tmp<-nrow(data))
      data <- data[ !data[,cases] | (data[,cases] & cond),  ]
      print(dim(data))
      print(table1(data[,cases]))
      
    }
  }
  
  # names(data) <- gsub("_last_before","",names(data) )  
  #
  
  # start date of sampling is the first vax date in the dataset:
  first_vax_date    <- min(data[,vax_date],na.rm=T)
  data$start_sample <- pmax(first_vax_date,  data[,begin_obs])
  data$stop_sample  <- pmin(data[,vax_date], data[,event_date], data[,end_obs], na.rm=T)
  
  data$days_before_vax <- as.numeric(difftime(data$stop_sample,data$start_sample,units="days"))
  data$days_before_vax[data$days_before_vax<0 | is.na(data$days_before_vax)] <- 0
  data$days_before_vax[duplicated(data[,id])] <- 0
  
  data <- data[order(data[,id],as.numeric(data[,vax_date])),]
  
  data$days_before_vax_cum_stop  <- cumsum(data$days_before_vax)
  data$days_before_vax_cum_start <- c(1,data$days_before_vax_cum_stop[-nrow(data)])
  
  
  #  ???????????????     
  #? hist(data$days_before_vax[!duplicated(data[,id])],1000)
  #  table(data$days_before_vax[!duplicated(data[,id])])
  
  # sum(data$days_before_vax[!duplicated(data[,id])]) # 4314594 days(age<30);   2808219103 (all)
  
  # table1(data$vax_n)
  
  if(sum(data$days_before_vax)==0) stop("no data for sampling")
  
  n_cases <- sum(data[,cases], na.rm=T)
  cat(paste0( "#days for sampling = ", sum(data$days_before_vax), ",  #cases = ", n_cases, "\n"))
  
  i<-1; n_controls<-0
  while(T) {
    # sample days:
    sampled_n <- sample.int(sum(data$days_before_vax),n_cases-n_controls)
    #sampled rows (category):
    table(table(sampled_cat <- findInterval(sampled_n, c(0,data$days_before_vax_cum_stop))))
    
    # controls:
    if(i>1) cdata0 <- cdata
    cdata <- data[sampled_cat,]
    cdata$sampled_days<- sampled_n 
    cdata$sampled_date <- as.Date(cdata$start_sample) + (cdata$sampled_days - cdata$days_before_vax_cum_start)
    cdata$index_date <- cdata$sampled_date
    cdata$case_control <- 0
    
    # select id's with conditions before sampled date (example: nocovid_before_vax_days during 30 days, or no myocarditits before vax during 1 year)
    if(!missing(cond_before_vax))
      for(icond in 1:length(cond_before_vax)){
        names_order <- names(cdata)
        cdata <- add_td_covariate( cdata[, !(names(cdata) %in% paste0(cond_before_vax[[icond]]$cov_name,c("_date_","_"),cond_before_vax[[icond]]$method)) ], 
                                   cond_before_vax[[icond]]$cov_name, methods=cond_before_vax[[icond]]$method,
                                   name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                   create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                   ids_study = ids_study  )    
        if(any(sort(names(cdata))!=sort(names_order))) {print(sort(names_order)); print(sort(names(cdata))); stop("problem: other variable names.")}
        cdata <- cdata[,names_order]
        
        diff_days <- as.numeric( difftime( data[,vax_date], data[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
        if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
        
        print(table1(cond))       
        print(dim(cdata))
        cdata <- cdata[cond,  ]
        print(dim(cdata))
      }
    if(i>1) cdata <- rbind(cdata0,cdata)
    
    
    if(nrow(cdata)>=n_cases){
      if(nrow(cdata)>n_cases) cdata <- cdata[1:n_cases,]
      break
    }   
    n_controls <- nrow(cdata)
    i <- i+1
  }  
  # cases:
  tmp <- data[ data[,cases], ];        
  tmp$sampled_days <- tmp$sampled_date <- NA
  tmp$index_date <- tmp[,vax_date]
  tmp$case_control <- 1
  
  # merge vaxed (cases) and sampled unvax dates (controls)
  dim(cdata); dim(tmp)
  data <- rbind.data.frame(tmp[,names(tmp)],cdata[,names(tmp)])
  dim(data)
  
  
  # data$index_days <- as.numeric(difftime(data$index_date,as.Date("2020-08-31")))
  
  data
}  # end function 'sample_controls'



add_periods <- function(data, data_new, start_state_date, stop_state_date, state_var, id="id_n", start_date="start_date", stop_date="stop_date", lsort=T, lprint=F){ 
  
  state_var <- tolower(state_var)
  data_new[,paste0(state_var,"_new")] <- ""
  data_new_rest <- data_new[,c(id,state_var,paste0(state_var,"_new"),start_state_date,stop_state_date)]
  npat <- max(table(data_new_rest[,id]))
  data_list <- vector("list",length=npat+1)
  data[,state_var] <- ""; data[,paste0(state_var,"_new")] <- "" ; data[,start_state_date] <- NA; data[,stop_state_date] <- NA
  for(i in 1:npat){ 
    catt(paste0("step ",i," (from ",npat,"): ", nrow(data_new_rest)," rows\n"))
    
    data_new_tmp <- data_new_rest[!duplicated(data_new_rest[,id]),]
    if(i>1) { data_new_tmp[,paste0(state_var,"_new")] <- data_new_tmp[,state_var]; data_new_tmp[,state_var]<-NULL }
    else data_new_tmp[,paste0(state_var,"_new")] <- ""
    
    data_list[[i]] <- data[!(data[,id] %in% data_new_tmp[!duplicated(data_new_tmp[,id]),id]),]
    data <- data[data[,id] %in% data_new_tmp[!duplicated(data_new_tmp[,id]),id],!(names(data) %in% c(start_state_date,stop_state_date,paste0(state_var,"_new")))]
    if(i==1) data <- data[,!(names(data) %in% state_var)] 
    if(lprint){ print(dim(data));print(dim(data_new_tmp)) }
    data <- merge.data.frame(data,data_new_tmp, by=id,all.x=T)
    if(lprint) print(dim(data))
    
    
    
    if(i==1) var0 <- var <- state_var
    if(i>1)  var <- paste0(state_var,"_new") 
    
    data[is.na(data[,var]),var] <- ""
    data[!is.na(data[,start_state_date]) &  data[,stop_date] < data[,start_state_date], var] <- ""
    data[!is.na(data[,stop_state_date])  &  data[,stop_state_date] < data[,start_date], var] <- ""
    
    
    # for start:
    cond_start <- !is.na(data[,start_state_date]) &  data[,start_date] < data[,start_state_date] & data[,start_state_date] <= data[,stop_date]
    if(any(cond_start)){
      
      data_start_extra <- data[cond_start,]
      data_start_extra[, start_date] <-  data_start_extra[, start_state_date] # 
      
      data[cond_start, stop_date] <- data[cond_start,start_state_date]-1 
      data[cond_start, var] <- ""  # ready       
      
      cond_start_stop <- !is.na(data_start_extra[,stop_state_date])  & data_start_extra[,start_date] <= data_start_extra[,stop_state_date] &  data_start_extra[,stop_state_date] < data_start_extra[,stop_date]
      if(any(cond_start_stop)){
        
        data_start_extra2 <- data_start_extra[cond_start_stop,]
        data_start_extra2[, start_date] <-  data_start_extra2[, stop_state_date]+1 
        data_start_extra2[, var] <-  "" 
        
        data_start_extra[cond_start_stop, stop_date] <- data_start_extra[cond_start_stop,stop_state_date]
      }
      data <- rbind.data.frame(data[,names(data)],data_start_extra[,names(data)])
      if(any(cond_start_stop)) data <- rbind.data.frame(data[,names(data)],data_start_extra2[,names(data)])
      
    }
    
    # for stop:
    cond_stop <- !is.na(data[,stop_state_date])  & data[,start_date] <= data[,stop_state_date] &  data[,stop_state_date] < data[,stop_date]
    if(any(cond_stop)){
      
      data_stop_extra <- data[cond_stop,]
      data_stop_extra[, start_date] <-  data_stop_extra[, stop_state_date] + 1 
      data_stop_extra[, var] <-  "" 
      
      data[cond_stop, stop_date] <- data[cond_stop,stop_state_date]  # ready 
      
      data <- rbind.data.frame(data[,names(data)],data_stop_extra[,names(data)])
    }
    
    
    if(i>1){ 
      for(j in 1:length(var0)) { cond <- !(as.character(data[,var0[j]]) %in% c("",NA,"NA")) & !(as.character(data[,var[j]]) %in% c("",NA,"NA")); data[cond , var0[j]] <- paste0( data[cond, var0[j]], " & ", data[cond, var[j]] ); data[cond , var[j]] <- "" }
      for(j in 1:length(var0)) { cond <-   as.character(data[,var0[j]]) %in% c("",NA,"NA")  & !(as.character(data[,var[j]]) %in% c("",NA,"NA")); data[cond, var0[j]]  <- data[cond, var[j]]; data[cond , var[j]] <- "" }
    }
    
    if(lsort) data <- data[order(data[,id],as.numeric(data[,start_date]),as.numeric(data[,stop_date])),]
    
    data_new_rest <- data_new_rest[duplicated(data_new_rest[,id]),] 
    
  } # end 'for'
  
  if(nrow(data)>0) data_list[[npat+1]] <- data
  
  colnames_in_list2 <- names(data_list[[2]])
  cat("create dataset...")
  data_list <- lapply(data_list, function(x)x[,colnames_in_list2])
  data <- do.call("rbind.data.frame",data_list)
  if(lsort) {cat("Sorting...");data <- data[order(data[,id],as.numeric(data[,start_date]),as.numeric(data[,stop_date])),]}
  cat("\n")
  
  data <- data[, !(names(data) %in% c(start_state_date, stop_state_date, paste0(state_var,"_new"))) ]
  print(table1(data[,state_var]))
  
  data
  
} # end of function 'add_periods'





hist_distribution <- function(dd, legend=c("vaccinated","sampled controls"),  
                              tit1="",
                              tit2=paste0("distribution of vaccinated and sampled per ",round(iperiod,1)," days"), 
                              periods=c(7,1,365/12),
                              smooth_line=T,
                              col=c("red","blue"),  alpha=0.7 ){
  
  vertical_lines <- as.numeric(difftime( as.Date(paste0("20",rep(18:24,each=4),"-",rep(c("01","04","07","10"),7),"-01")),as.Date("2020-08-31")),units="days")
  cal_time_range_days <- as.numeric(difftime( range( dd$index_date, na.rm=T), as.Date("2020-08-31"),units="days"))
  
  for(iperiod in periods){
    vax_int <- seq( cal_time_range_days[1],cal_time_range_days[2]+iperiod, by=iperiod) 
    h1 <- hist(dd$index_days[dd$case_control==1],breaks=vax_int,plot=F)
    h0 <- hist(dd$index_days[dd$case_control==0],breaks=vax_int,plot=F)
    xxlim <- range(dd$index_days)
    yylim <- c(0,max(h1$counts,h0$counts))
    plot(h1, xlim=xxlim, ylim=yylim, col=col[1], xlab="Time", ylab="Frequency",axes=F, main=paste(tit1,tit2))
    abline(v=vertical_lines,col="lightgray"); par(new=T)
    plot(h1, xlim=xxlim, ylim=yylim, col=col[1], xlab="", ylab="",axes=F, main=""); par(new=T)
    plot(h0, xlim=xxlim, ylim=yylim, col= rgb(t(col2rgb(col[2])/255),alpha=alpha), axes=F,xlab="", ylab="",main="")
    axis(2); axis(1,at=vertical_lines, labels =as.Date("2020-08-31")+vertical_lines, cex=0.8); box()
    if(length(legend)>0) legend("topright",legend=legend, fill=col, bty="n")
    if(smooth_line){
      try(lines( smooth.spline(h0$mids,h0$counts,df=15)  ,col=col[2]))
      try(lines( smooth.spline(h1$mids,h1$counts,df=15)  ,col=col[1]))
    }
  }
}

# par(mfrow=c(2,1))
# hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1) )
# par(mfrow=c(2,2))
# hist_distribution(sampled_data,legend=c("vaccinated","sampled controls"), tit1="Approach B: ", periods=c(7,1,365/12) )



hist_vars <- function(vars, legend=names(vars),  
                      tit1="",
                      tit2="", 
                      periods=c(7,1,365/12),
                      smooth_line=F,
                      col=c("yellow","deeppink","skyblue","aquamarine","cornsilk","azure"), alpha=0.6,
                      xlim, ylim, xlab="Time", ylab="Frequency",
                      lplot=T){  
  
  vertical_lines <- as.numeric(difftime( as.Date(paste0("20",rep(18:24,each=4),"-",rep(c("01","04","07","10"),7),"-01")),as.Date("2020-08-31")),units="days")
  cal_time_range_days <- as.numeric(difftime( as.Date(range( unlist(vars),na.rm=T),origin="1970-01-01"), as.Date("2020-08-31") ,units="days"))
  vars_days <- lapply(vars, function(x) as.numeric(difftime( x, as.Date("2020-08-31"), units="days")) )                                  
  lplot <- rep(lplot,2)
#browser()
  res<-vector("list",length=length(periods)); names(res) <- paste0("period_",periods,"days")
  for(iperiod in periods){
    
    vax_int <- seq( cal_time_range_days[1]-iperiod, cal_time_range_days[2]+iperiod, by=iperiod) 
    vars_hist <- lapply(vars_days, function(x) hist(x,breaks=vax_int,plot=F) )
    
    res[[paste0("period_",iperiod,"days")]] <- vars_hist
    
    if(missing(xlim)) xxlim <- range(unlist(vars_days),na.rm=T)  else  xxlim <- xlim
    if(missing(ylim)) yylim1 <- c(0,  max( sapply(vars_hist, function(x)max(x$counts,na.rm=T) )) )  else  yylim1<- ylim   
  
    yy_max_small <- max(sapply(vars_hist,function(x)if(any(x$counts!=max(x$counts,na.rm=T))) max(x$counts[x$counts!=max(x$counts,na.rm=T)]) else 0 ),na.rm=T)

    if(any(lplot)){
      for(i in 1:2){
        if(!lplot[i]) next
        if(i==2 & (yy_max_small > 0.5*yylim1[2] | yylim1[2]<=15)) next
        yylim <- yylim1; if(i==2) yylim[2] <- yy_max_small
        plot( vars_hist[[1]], xlim=xxlim, ylim=yylim, col=col[1], xlab=xlab, ylab=ylab, axes=F, main=paste(tit1,tit2))
        abline(v=vertical_lines,col="lightgray")
        axis(2); axis(1,at=vertical_lines, labels =as.Date("2020-08-31")+vertical_lines, cex=0.8); box(); par(new=T)
        plot(vars_hist[[1]], xlim=xxlim, ylim=yylim, col=col[1], xlab="", ylab="",axes=F, main=""); par(new=T)
        if(length(vars)>1) for(ivar in 2:length(vars)){
          plot(vars_hist[[ivar]], xlim=xxlim, ylim=yylim, col=rgb(t(col2rgb(col[ivar])/255),alpha= alpha) , xlab="", ylab="",axes=F, main=""); par(new=T) }
        if(missing(legend)) legend <- names(vars)
        if(length(legend)>0) legend("topright",legend=legend, fill=col[1:length(vars)], bty="n")
        if(smooth_line) for(ivar in 2:length(vars))
          try(lines( smooth.spline(vars_hist[[ivar]]$mids,vars_hist[[ivar]]$counts,df=20)  ,col=col[ivar]))
        if(i==2 & "vars_hist_small" %in% ls()) vars_hist <- vars_hist0
        par(new=F)
      }
    }
  }
  
  attributes(res) <- c(attributes(res),
                       list(vertical_lines=vertical_lines, cal_time_range_days=cal_time_range_days ))
  res
}




id_add <- function(data, ids_study, id, id_new="id_n",lprint=T){
  if( any(!(tmp<-(data[!duplicated(data[,id]),id] %in% ids_study))) ){ 
    if(lprint){
      cat(paste0(sum(!tmp,na.rm=T)," id's are not in the study population.\n"))
      print(tmp<-dim(data))
    }
    data <- data[data[,id] %in% ids_study,]
    if(lprint) {print(dim(data));cat(paste("==> ", tmp[1]-nrow(data)," rows deleted.\n\n"))}
  }
  data[,id_new] <- match(data[,id],ids_study)
  data
}



match_pop <- function(dd0, vax_time, vax1_time, vax_date="vax_date", #event_time, 
                      cases,            # name of logical variable (example, T if vax_n=1 & vax_brand=="Pfizer" ==> sampling only for Pfizer vax1)
                      matching_indep_factors=c(), 
                      cond_before_vax,  # list with elements: list( list( var="covid_date_last_before",       period= 30, unit="days", unit_short="d" ),
                      #                   list( var="myocarditis_date_last_before", period=365, unit="days", unit_short="d" )     )
                      #                   list( var="myocarditis_date_last_before", period=365, unit="days", unit_short="d" )     )
                      matching_dep_factors=c(), # list of vectors:  list( c( date="covid_date_last_before", event="covid_last_before", abs_diff=30,       units="days" ),
                                                #                         c( date="birth_date"                                       , abs_diff=2*365.25, units="days" )     )
                      ids_study,
                      vax_vars=c("vax_brand","vax_n"), 
                      start_interval="study_entry_days", stop_interval="study_exit_days", death_time="death_days",
                      result_name="", file_name="", dir="", nboot=0, lprint=F, lprint_str=T, lprint_warn=T, lprint_all=F, create_i=T, lparal=F){
  
  if(lprint_str | lprint) print(Sys.time())
  
  if(missing(cases)){ cases <- "_tmp_cases";dd0[,"_tmp_cases"] <- !is.na(dd0[,vax_date]) }
  cond <- !dd0[,cases] & !is.na(dd0[,vax_date]) &  dd0[,vax_date]  <= dd0[,stop_interval]
  dd0[  cond, stop_interval ] <- dd0[  cond, vax_date ] - 1
  
  # for vaccinated:
  if(!missing(cond_before_vax)){
    if(missing(ids_study)) stop(paste("'ids_study' must be specified for sampling"))
    for(icond in 1:length(cond_before_vax)){  
      if( !("unit"           %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$unit           <- "days"
      if( !("unit_short"     %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$unit_short     <- "d"
      if( !("create_cov_var" %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$create_cov_var <- F
      if( !("method"         %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$method         <- "last_before"
      if( !("condition"      %in% names(cond_before_vax[[icond]])) ) cond_before_vax[[icond]]$condition      <- ""
      cond_before_vax[[icond]]$var <- paste0(cond_before_vax[[icond]]$cov_name,"_date_",cond_before_vax[[icond]]$method)
      
      diff_days <- as.numeric( difftime( dd0[,vax_date], dd0[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
      if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
   
      if( !(cond_before_vax[[icond]]$var %in% names(dd0)) )
        dd0 <- add_td_covariate( dd0, cond_before_vax[[icond]]$cov_name, methods=cond_before_vax[[icond]]$method, 
                                  cov_dataset = cond_before_vax[[icond]]$dataset,
                                  #name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                  create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                  ids_study = ids_study, lprint=F  )    
      
      if(cond_before_vax[[icond]]$condition!=""){
        diff_days <- as.numeric( difftime( dd0[,vax_date], dd0[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
        if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
      }
      dd0[!cond,stop_interval] <- dd0[!cond,vax_date] - 1
      if(any(!cond)) dd0[!cond,cases]         <- F 
      
      #print(table1(cond))       
      #print(tmp<-nrow(dd0))
      dd0 <- dd0[ !dd0[,cases] | (dd0[,cases] & cond),  ]
      #print(dim(dd0))
      #print(table1(dd0[,cases]))
      
    }
  }
  
  nboot <- nboot+1
  start_interval_all <- dd0[,start_interval]
  end_obs_all <- dd0[,stop_interval]
  #end_novax_or_death <- pmin(dd0[,vax_time]-1,dd0[,stop_interval],dd0[,death_time],na.rm=T)
  #event_all <- dd0[,event_time]
  if(create_i) dd0$i <- 1:nrow(dd0)
  cond_novax <- is.na(dd0[,vax_time])
  vax_days_all  <- dd0[,vax_time]  
  vax1_days_all <- dd0[,vax1_time]  
  idata <- dd0$i
  cond_vax <- !cond_novax  & start_interval_all<=vax_days_all & vax_days_all<=end_obs_all
  #cond_event_after_vax <- !cond_novax &   !( !is.na(dd0[,event_time]) & dd0[,event_time] < vax_days_all ) & start_interval_all<=vax_days_all & vax_days_all<=end_obs_all
  icasedata <- dd0$i[cond_vax]
  vax_days  <- dd0[cond_vax,vax_time]  
  exposure_days_unique <- sort(vax_days[!duplicated(vax_days)])
  
  
  
  
  
 
  
  if(length(matching_indep_factors)==0) dd0$matching_factor <- ""
  else dd0$matching_factor <- format(dd0[,matching_indep_factors[1]])
  
  
  if(length(matching_indep_factors)>1) for(i in 2:length(matching_indep_factors)) dd0$matching_factor <- paste(dd0$matching_factor,format(dd0[,matching_indep_factors[i]]))
  matching_factor_levels <- levels(as.factor(dd0$matching_factor))
  dd0$matching_factor <- as.numeric(as.factor(dd0$matching_factor))
  
  
  # N3: 
  ####### with strata:                   ==>  31 min 33 sec
  (start_time <- Sys.time())
  
  if(sum(is.na(dd0$matching_factor))>0) cat(paste0("\n",sum(is.na(dd0$matching_factor))," rows with missing values in matching variables.\n\n"))
  
  if(!lparal){
    str_tab <- table(dd0[cond_vax,"matching_factor"]) #; str_tab <- str_tab[!is.na(str_tab)]
    str_tab_len <- length(str_tab)
    if(length(matching_indep_factors)>0){
      if(lprint){
        cat(paste0( 'matching factors: "',paste0(matching_factor_levels,collapse='","'),'"\n' ))
        print(c( str_tab, Sum=sum(str_tab) ))
      }
      str_case_nlevels <- names(str_tab)[length(str_tab)]
      str_tab_names <- as.numeric(names(str_tab)); names(str_tab_names) <- as.character(1:length(str_tab_names)) 
    }
    
    matched_cases    <- c()
    matched_controls <- vector("list",length=nboot); names(matched_controls) <- paste0("boot_",0:(nboot-1))
    less_controls_matched_cases <- less_controls_unmatched_cases <- list()
  
    for(istr in names(str_tab)){ 
      istr_i<-(1:str_tab_len)[names(str_tab)==istr]
      if(length(matching_indep_factors)>0){
        str_text <- paste0(" for stratum ", istr_i," ('", matching_factor_levels[as.numeric(istr)],"')")
        if( lprint_str & ( lprint | ( !lprint & istr_i %% 10 == 0 ) ) ) 
          cat(paste0("stratum ", istr_i," (",istr, ": '", matching_factor_levels[as.numeric(istr)],"') from ",str_tab_len," (",str_case_nlevels,");      "))
      } else str_text <- ""
      str_cond_all_data <- dd0$matching_factor==istr
      str_cond_cases    <- str_cond_all_data[cond_vax] 
      if(sum(str_cond_cases)==0) next
      
      cond_novax_str         <- cond_novax[str_cond_all_data]
      start_interval_str     <- start_interval_all[str_cond_all_data]
      end_obs_all_str        <- end_obs_all[str_cond_all_data]
      vax_days_str           <- vax_days[str_cond_cases]
      idata_str              <- idata[str_cond_all_data]
      icasedata_str          <- icasedata[str_cond_cases]
      vax_days_all_str       <- vax_days_all[str_cond_all_data]
      vax1_days_all_str      <- vax1_days_all[str_cond_all_data]
      #event_str              <- event_all[str_cond_all_data]
      

      exposure_days_unique_str <- sort(vax_days_str[!duplicated(vax_days_str)])
      matched_cases_list_str <- vector("list",length(exposure_days_unique_str))
      names(matched_cases_list_str) <- exposure_days_unique_str
      matched_controls_list_str <- matched_cases_list_str 
      
      
      if( lprint_str |  lprint | ( !lprint & istr_i %% 10 == 0) ) 
        cat(paste0(str_text, ":   #case-rows = ",length(icasedata_str), ";  number of days = ",length(exposure_days_unique_str), ";   ", Sys.time(), "\n"))
      
      for(iday in 1:length(exposure_days_unique_str)){ #print(iday) ; #if(iday==3)browser()
        ivaxday <- exposure_days_unique_str[iday]
        cond_contr <- start_interval_str<=ivaxday  & ivaxday<=end_obs_all_str # & !( !is.na(event_str) & event_str<ivaxday ) 
        cond_contr <- cond_contr &  ( is.na(vax1_days_all_str) | (!is.na(vax1_days_all_str) &  ivaxday  < vax1_days_all_str ) )  
        #cond_contr <- cond_contr &  !( !is.na(vax1_days_all_str) & vax1_days_all_str <= ivaxday )  

        # for unvaccinated periods or persons
        # select id's with conditions before sampled date (example: nocovid_before_vax_days during 30 days, or no myocarditits before vax during 1 year)
        if(!missing(cond_before_vax) & any(cond_contr)){
          for(icond in 1:length(cond_before_vax)){
            cdata <- dd0[ match( idata_str,dd0$i), c(id, vax_date)]
            cdata <- cdata[cond_contr,]
            cdata$ii <- 1:nrow(cdata); nrow0 <- nrow(cdata)
            if(nrow(cond_before_vax[[icond]]$dataset)==0) next
            cdata <- add_td_covariate( cdata[, !(names(cdata) %in% paste0(cond_before_vax[[icond]]$cov_name,c("_date_","_"),cond_before_vax[[icond]]$method)) ], 
                                       cond_before_vax[[icond]]$cov_name, methods=cond_before_vax[[icond]]$method,
                                       cov_dataset = cond_before_vax[[icond]]$dataset,
                                       #name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                       create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                       ids_study = ids_study, lprint=F  )    
            #if(any(sort(names(cdata))!=sort(names_order))) {print(sort(names_order)); print(sort(names(cdata))); stop("problem: other variable names.")}
            if(nrow(cdata)!= nrow0) {stop("problem: # rows should be the same.")} 
            cdata <- cdata[cdata$ii,]
            
            if(cond_before_vax[[icond]]$condition!=""){
              diff_days <- as.numeric( difftime( as.Date(ivaxday,"2020-08-31"), cdata[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
              if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])
            }
            #else
            cond_contr[cond_contr] <- cond
            #print(table1(cond))       
            #print(dim(cdata))
            #cdata <- cdata[cond,  ]
            #print(dim(cdata))
          }
        } # end if: !missing(cond_before_vax) & any(cond_contr) 





          
        if(!missing(matching_dep_factors) & any(cond_contr)){
          
         dep_matrix <- sapply( matching_dep_factors, function(x) {
            if(length(x)==4) if(is.null(names(x))) names(x) <- c("date","event","abs_diff","units") else x<-x[c("date","event","abs_diff","units")] 
            if(length(x)==3){ x <- c(x[1], gsub("_date","",x[1]), x[2:3]); names(x) <- c("date","event","abs_diff","units")}
            if(length(x)<3) stop("vector for 'matching_dep_factors' should have length 3 of 4 ('date','event','abs_diff','units') or ('date', 'abs_diff','units')")
            x
          } )
         
         # cases for day 'ivaxday' 
         csdata0 <- dd0[ icasedata_str[vax_days_str==ivaxday], names(dd0) %in%c("i",id, dep_matrix[c("date","event"),]), drop=F ]
         for(i in 1:length(cond_before_vax)) if( !(dep_matrix["event",i] %in% names(csdata0)) ) csdata0[,dep_matrix["event",i]] <- 1
      
         non_empty_case_events <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix[c("event","date"),, drop=F]]) )>0 , dep_matrix[c("event","date"),] ]
         #non_empty_case_events <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix["date",, drop=F]]) )>0 , dep_matrix["event",] ]
         non_empty_case_events <- non_empty_case_events[!duplicated(non_empty_case_events),]

       if(nrow(non_empty_case_events)>0){
       
           for(irow in 1:nrow(non_empty_case_events)){  # per each combi in cases
           # potential matched controls
            contr_available <-  dd0[ (match_contr <- match( idata_str[cond_contr] , dd0$i )) , 
                              c("i",id, vax_date, unlist(matching_dep_factors)[unlist(matching_dep_factors) %in% names(dd0)] ) ]
            #cdata <- dd0[ match( idata_str,dd0$i), c(id, vax_date, unlist(matching_dep_factors) ) ]
            #cdata <- cdata[cond_contr,]
            for(i in 1:length(cond_before_vax)) if( !(dep_matrix["event",i] %in% names(contr_available)) ) contr_available[,dep_matrix["event",i]] <- 1
            contr_available[, paste0( dep_matrix["event",],  "_case")] <- non_empty_case_events[ irow, dep_matrix["event",],drop=F ]
            contr_available[, paste0( dep_matrix["date" ,],  "_case")] <- non_empty_case_events[ irow, dep_matrix["date" ,],drop=F ]
            contr_available <- contr_available[ apply(contr_available[ , dep_matrix["event",] , drop=F],1,paste, collapse=" ") == paste( non_empty_case_events[ irow, dep_matrix["event",],drop=F ], collapse=" "), ]
            if(nrow(contr_available)>0) 
              for(ii in 1:length(dep_matrix["date",])){
                ivar <- dep_matrix["date",ii]
                if(!is.na(non_empty_case_events[irow,ivar]))
                  contr_available <- contr_available[ abs(as.numeric( difftime(non_empty_case_events[irow,ivar], contr_available[,ivar], units=dep_matrix["units",ii] ))) < as.numeric( dep_matrix["abs_diff",ii] ), ]
              }
           
            contr_available$case_combi <- apply( contr_available[,dep_matrix[c("event","date"),] ], 1, paste,collapse=" ")  
 
            ncontr  <- nrow(contr_available)  # ==> dd0[contr_available$i,] - controls in dd0
            csd <- csdata0[ apply(csdata0[ , dep_matrix[c("event","date"),] ],1,paste,collapse=" ") == apply( non_empty_case_events[ irow,,drop=F], 1, paste,collapse=" "),]
            ncases <- nrow(csd)   #   ==>  dd0[csd$i,] -  cases in dd0
            
            icases_str <- csd$i #csdata0$i[cond_cases]
            #non_empty_case_events <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix[c("event","date"),, drop=F]]) )>0 , dep_matrix[c("event","date"),] ]
            
            #icases_str <- icasedata_str[vax_days_str==ivaxday]
            #ncases  <- length(icases_str)
            
            
            
          #  if(lprint_str & iday%%10==0 & !lprint_all) cat(paste0(  iday, " "))
          #  if(lprint_str & iday%%50==0 & !lprint_all) cat(paste0( "\nday = ", iday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
          #  if(lprint_all) cat(paste0( "day=", iday, ";    time = ",ivaxday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
          #  
            if(ncases > ncontr) {
              if(lprint_warn) cat(paste0("Less potential controls for time = ",ivaxday," (day=", iday,");  #cases = ",ncases,";  #controls = ",ncontr, str_text," & ",apply( non_empty_case_events[ irow,,drop=F], 1, paste,collapse=" "), "\n"))  
              if(ncontr==0){ 
                less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                         matrix(rep(icases_str,nboot), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
                next
              }
              else {
                ii <- lapply(1:nboot, function(i,n,size) sample.int(n,size), n=ncases, size=ncontr)
                less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                         matrix(unlist(lapply(ii, function(x,rows)rows[!(c(1:length(rows)) %in% x)], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
                less_controls_matched_cases   <- c(less_controls_matched_cases,   list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                         matrix(unlist(lapply(ii, function(x,rows)rows[x], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1)))) )))
                
                
                matched_cases_list_str[[iday]] <- c( matched_cases_list_str[[iday]], icases_str[ii[[1]]] )
                ncases <- ncontr
              }
            } # enf if: ncases > ncontr 
            else matched_cases_list_str[[iday]] <- c(matched_cases_list_str[[iday]] , icases_str)
         
            if(length(matched_controls_list_str[[iday]])>0)
              matched_controls_list_str[[iday]] <- rbind.data.frame( matched_controls_list_str[[iday]],   
                                                                     matrix(unlist(lapply(1:nboot, function(i,n,size,idata) idata[sample.int(n,size)], n=ncontr, size=ncases, idata=idata_str[cond_contr] )),ncol=nboot) )
            else
              matched_controls_list_str[[iday]] <- matrix(unlist(lapply(1:nboot, function(i,n,size,idata) idata[sample.int(n,size)], n=ncontr, size=ncases, idata=idata_str[cond_contr] )),ncol=nboot)

            
            
          
          } # end for irow
          
          
          
          
          
    if(F) {     
          # variables with non-missing values for controls
          ii <- (1:length(matching_dep_factors))[colSums(!is.na(cntdata0[, dep_matrix["date",] ]))>0]
          # matching_dep_factors[colSums(!is.na(cdata[, dep_matrix["date",] ]))>0]
 
          # if there are non-missing date in cntdata dataset:
          if(length(ii)>0) for(idep in ii){
              
             # potential conrols 
            # cntdata <- cntdata0[, names(cntdata0) %in% c(id, dep_matrix[c("date","event"),ii]), drop=F ]
              
            # cases for day 'ivaxday' 
            csdata0 <- dd0[ icasedata_str[vax_days_str==ivaxday], names(dd0) %in%c(id, dep_matrix[c("date","event"),ii]), drop=F ]
            
            if( !(dep_matrix["event",ii] %in% names(csdata0)) ){ csdata0[,dep_matrix["event",ii]] <- 1; cntdata[,dep_matrix["event",ii]] <- 1 }
            
            
            for(ievent in csdata0[ !duplicated(csdata0[,dep_matrix["event",ii]]), dep_matrix["event",ii]] ){
              csdata <- csdata0[,dep_matrix["event",ii]]==ievent
              for(idate in csdata[!duplicated(csdata[,dep_matrix["date",ii]]),dep_matrix["date",ii]]){
                cond <- abs(as.numeric( difftime(idate, cntdata[,dep_matrix["date",ii]], units=dep_matrix["units",ii] ))) < as.numeric( dep_matrix["abs_diff",ii] )
                if(any(cond)) { cntdata[cond, idate_case] <-  idate; }
                
              }
              
              
            } #
               
            
          } # end 'for' idep
          
    }
          
       } # if nrow(non_empty_case_events)>0
         
      } # end if !missing(matching_dep_factors) & any(cond_contr)







      if(missing(matching_dep_factors) & any(cond_contr)){ # end if missing(matching_dep_factors) & any(cond_contr)
        
        ncontr  <- sum(cond_contr)
        icases_str <- icasedata_str[vax_days_str==ivaxday]
        ncases  <- length(icases_str)
        
        
        
        if(lprint_str & iday%%10==0 & !lprint_all) cat(paste0(  iday, " "))
        if(lprint_str & iday%%50==0 & !lprint_all) cat(paste0( "\nday = ", iday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
        if(lprint_all) cat(paste0( "day=", iday, ";    time = ",ivaxday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
        
        if(ncases > ncontr) {
          if(lprint_warn) cat(paste0("Less potential controls for time = ",ivaxday," (day=", iday,");  #cases = ",ncases,";  #controls = ",ncontr, str_text, "\n"))  
          if(ncontr==0){ 
            less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                     matrix(rep(icases_str,nboot), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
            next
          }
          else {
            ii <- lapply(1:nboot, function(i,n,size) sample.int(n,size), n=ncases, size=ncontr)
            less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                     matrix(unlist(lapply(ii, function(x,rows)rows[!(c(1:length(rows)) %in% x)], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
            less_controls_matched_cases   <- c(less_controls_matched_cases,   list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                     matrix(unlist(lapply(ii, function(x,rows)rows[x], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1)))) )))
            
            
            matched_cases_list_str[[iday]] <- icases_str[ii[[1]]]
            ncases <- ncontr
          }
        } # end if: ncases > ncontr 
        else matched_cases_list_str[[iday]] <- icases_str
        
        matched_controls_list_str[[iday]] <- matrix(unlist(lapply(1:nboot, function(i,n,size,idata) idata[sample.int(n,size)], n=ncontr, size=ncases, idata=idata_str[cond_contr] )),ncol=nboot)
       
        
      } # end of missing(matching_dep_factors) & any(cond_contr)



if(F){

#????          
          
          ncontr  <- sum(cond_contr)
          icases_str <- icasedata_str[vax_days_str==ivaxday]
          ncases  <- length(icases_str)
          
          
          
          if(lprint_str & iday%%10==0 & !lprint_all) cat(paste0(  iday, " "))
          if(lprint_str & iday%%50==0 & !lprint_all) cat(paste0( "\nday = ", iday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
          if(lprint_all) cat(paste0( "day=", iday, ";    time = ",ivaxday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"))
          
          if(ncases > ncontr) {
            if(lprint_warn) cat(paste0("Less potential controls for time = ",ivaxday," (day=", iday,");  #cases = ",ncases,";  #controls = ",ncontr, str_text, "\n"))  
            if(ncontr==0){ 
              less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                       matrix(rep(icases_str,nboot), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
              next
            }
            else {
              ii <- lapply(1:nboot, function(i,n,size) sample.int(n,size), n=ncases, size=ncontr)
              less_controls_unmatched_cases <- c(less_controls_unmatched_cases, list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                       matrix(unlist(lapply(ii, function(x,rows)rows[!(c(1:length(rows)) %in% x)], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1))) ) )))  
              less_controls_matched_cases   <- c(less_controls_matched_cases,   list(cbind.data.frame( strata_n=istr, ivaxday=ivaxday, 
                                                                                                       matrix(unlist(lapply(ii, function(x,rows)rows[x], rows=icases_str)), ncol=nboot, dimnames=list(c(),paste0("boot_",0:(nboot-1)))) )))
              
              
              matched_cases_list_str[[iday]] <- icases_str[ii[[1]]]
              ncases <- ncontr
            }
          } # enf if: ncases > ncontr 
          else matched_cases_list_str[[iday]] <- icases_str
          
          matched_controls_list_str[[iday]] <- matrix(unlist(lapply(1:nboot, function(i,n,size,idata) idata[sample.int(n,size)], n=ncontr, size=ncases, idata=idata_str[cond_contr] )),ncol=nboot)
          
          
          
}        
          
          
          
          
        }  # end "for" iday 
        
        n_day_matched_cases_str      <- sapply(matched_cases_list_str,length)
        n_day_matched_controls_str   <- unlist(lapply(matched_controls_list_str, nrow))
        #n_day_less_controls_list_str <- sapply(less_controls_list_str,function(x)length(x$unmatched_cases))
        
        if(!all(n_day_matched_cases_str[n_day_matched_cases_str>0]==n_day_matched_controls_str[n_day_matched_controls_str>0])) { print("problem 1.1");}  #browser() }
        
        matched_cases   <- c(matched_cases, unlist(matched_cases_list_str) )
        #matched_day    <- rep(as.numeric(names(n_day_matched_cases)),n_day_matched_cases)
        #matched_controls <- vector("list",length=nboot+1); names(matched_controls) <- paste0("boot_",0:nboot)
        for(iboot in 1:nboot)
          matched_controls[[paste0("boot_",iboot-1)]] <- c(matched_controls[[paste0("boot_",iboot-1)]], unlist(lapply(matched_controls_list_str,function(x,iboot)x[,iboot],iboot=iboot)) )

        gc()
        
      }# end 'for' istr
      
      if(length(less_controls_matched_cases)>0){
        less_controls_matched_cases   <- do.call("rbind.data.frame",less_controls_matched_cases)
        less_controls_matched_cases$strata <- matching_factor_levels[as.numeric(less_controls_matched_cases$strata_n)]
        less_controls_matched_cases   <- less_controls_matched_cases[,  c(ncol(less_controls_matched_cases),1:(ncol(less_controls_matched_cases)-1))]
      }
      if(length(less_controls_unmatched_cases)>0){
        less_controls_unmatched_cases <- do.call("rbind.data.frame",less_controls_unmatched_cases)
        less_controls_unmatched_cases$strata <- matching_factor_levels[as.numeric(less_controls_unmatched_cases$strata_n)]
      less_controls_unmatched_cases <- less_controls_unmatched_cases[,  c(ncol(less_controls_unmatched_cases),1:(ncol(less_controls_unmatched_cases)-1))]
    }
  }
  
  
  print(Sys.time()-start_time); print(Sys.time())
  
  
  if(length(matching_indep_factors)==0){matching_factor_levels <- c()}
  
  match_res <- list( matched_cases=matched_cases, matched_controls=matched_controls, 
                     less_controls_matched_cases=less_controls_matched_cases, less_controls_unmatched_cases=less_controls_unmatched_cases, 
                     cond_before_vax = lapply(cond_before_vax,function(x)x[names(x)!="dataset"]),
                     matching_indep_factors=matching_indep_factors, matching_factor_levels=matching_factor_levels, matching_dep_factors=matching_dep_factors,
                     parameters = list( vax_time=vax_time, vax1_time=vax1_time, #event_time=event_time, 
                                        matching_indep_factors=matching_indep_factors, matching_dep_factors=matching_dep_factors, cond_before_vax=lapply(cond_before_vax,function(x)x[names(x)!="dataset"]),
                                        vax_vars=vax_vars, 
                                        start_interval=start_interval, stop_interval=stop_interval, death_time=death_time,
                                        result_name=result_name, file_name=file_name, dir=dir )
  )
  
  if(file_name!=""){ 
    if(result_name!="") assign(result_name,match_res)
    else result_name <- "match_res"  
    if(dir!="") file_name <- paste0(dir,file_name)
    save( list=result_name, file=file_name )
  }
  
  if(lprint_str | lprint) print(Sys.time())
  
  match_res 
  
} # end of function 'match_pop'





get_matched_dataset <- function(matched_list, data_vax, update_event=T, stop_obs="study_exit_days", next_vax_time="", iboot=0, create_i=F, lprint=T){
  
  vax_time <- matched_list$parameters$vax_time
  event_time <- matched_list$parameters$event_time
  matching_indep_factors <- matched_list$parameters$matching_indep_factors
  vax_vars <- matched_list$parameters$vax_vars
  start_interval <- matched_list$parameters$start_interval
  stop_interval <- matched_list$parameters$stop_interval
  death_time <- matched_list$parameters$death_time
  result_name <- matched_list$parameters$result_name
  file_name <- matched_list$parameters$file_name
  
  
  if(create_i) data_vax$i <- 1:nrow(data_vax)
  data_vax$matching_factor <- ""
  if(length(matched_list$matching_indep_factors)>0){
    for(ivar in matched_list$matching_indep_factors)
      data_vax$matching_factor <- paste0(data_vax$matching_factor,"_",data_vax[,ivar])
    matching_factor_levels <- levels(as.factor(data_vax$matching_factor))
    data_vax$matching_factor <- as.numeric(as.factor(data_vax$matching_factor))
  }
  
  # cases:
  if(lprint) cat(paste0("#rows in dataset = ",nrow(data_vax),";   #rows with vax date = ",sum(!is.na(data_vax[,vax_time])), " (",round(100*sum(!is.na(data_vax[,vax_time]))/nrow(data_vax),2),"%);  " ))
  tmp <- data_vax[matched_list[["matched_cases"]],]
  if(lprint) cat(paste0( "#matched cases = ", nrow(tmp)," (",round(100*nrow(tmp)/sum(!is.na(data_vax[,vax_time])),2),"%)     and      "))
  
  tmp$i_vaxed   <- matched_list[["matched_cases"]]
  tmp$i_unvaxed <- matched_list[["matched_controls"]][[paste0("boot_",iboot)]]
  tmp[,paste0(vax_vars,"_case")] <- tmp[,vax_vars]
  tmp$case_control <- 1
  
  # controls:
  data_vax <- data_vax[matched_list[["matched_controls"]][[paste0("boot_",iboot)]],]
  if(lprint) cat(paste0("#matched controls = ",nrow( data_vax),"\n"))
  
  data_vax$i_vaxed   <- matched_list[["matched_cases"]]
  data_vax$i_unvaxed <- matched_list[["matched_controls"]][[paste0("boot_",iboot)]]
  data_vax[,vax_time] <- tmp[,vax_time]
  data_vax[,paste0(vax_vars,"_case")] <- tmp[,vax_vars]
  data_vax$case_control <- 0
  
  if(nrow(tmp)!=nrow(data_vax)) stop("#cases and #controls should be the same!")
  
  data_vax <- rbind.data.frame(tmp, data_vax)
  
  
  if(update_event & event_time!="" & event_time %in% names(data_vax)){    
    #data_vax
    data_vax$event_or_end_days <- with(data_vax, pmin(get(event_time), get(stop_obs), get(death_time), na.rm=T) ) 
    data_vax[,iae] <- as.numeric( !is.na(data_vax[,event_time]) &  data_vax[,event_time] == data_vax$event_or_end_days )
    data_vax$event_after_vax_days <- with(data_vax, event_or_end_days - get(vax_time)) 
    
    tmp <- pmin( data_vax$event_or_end_days[data_vax$case_control==1], data_vax$event_or_end_days[data_vax$case_control==0] )
    data_vax$event_or_end_days_matched_min <- c( tmp, tmp )
    data_vax[,paste0(iae,"_matched_min")] <- as.numeric( !is.na(data_vax[,event_time]) &  data_vax[,event_time] == data_vax$event_or_end_days_matched_min )
    data_vax$event_after_vax_days_matched_min <- with(data_vax, event_or_end_days_matched_min - get(vax_time)) 
    
    if(next_vax_time!="" & next_vax_time %in% names(data_vax)){
      data_vax$event_or_next_vax_or_end_days <- with(data_vax, pmin( get(next_vax_time), event_or_end_days, na.rm=T) ) 
      data_vax[,paste0(iae,"_before_next_vax")] <- as.numeric( !is.na(data_vax[,event_time]) &  data_vax[,event_time] == data_vax$event_or_next_vax_or_end_days )
      data_vax$event_or_next_after_vax_days <- with(data_vax, event_or_next_vax_or_end_days - get(vax_time)) 
      
      tmp <- pmin( data_vax$event_or_next_vax_or_end_days[data_vax$case_control==1], data_vax$event_or_next_vax_or_end_days[data_vax$case_control==0] )
      data_vax$event_or_next_vax_or_end_days_matched_min <- c( tmp, tmp )
      data_vax[,paste0(iae,"_before_next_vax_matched_min")] <- as.numeric( !is.na(data_vax[,event_time]) &  data_vax[,event_time] == data_vax$event_or_next_vax_or_end_days_matched_min )
      data_vax$event_or_next_after_vax_days_matched_min <- with(data_vax, event_or_next_vax_or_end_days_matched_min - get(vax_time)) 
    } 
  }
  
  data_vax
  
} # end of func 'get_matched_dataset'



tabb <-  function(tb, text_n="n", text_percent='%',lorder=T, decreasing=T){  
  if(missing(text_n) & !is.null(names(tb))) text_n <- names(tb)
  tb <- lapply(tb,function(x,all_names) if(!is.null(dim(x))) x[,1] )
  all_names <- unique(unlist(lapply(tb,names)))
  tb <- lapply(tb,function(x,all_names){  
    if(!is.null(dim(x))) x <- x[,1]
    if(length(all_names)!=length(x)){
      res<-c(x,rep(0,length(all_names)-length(x))) 
      names(res)<-c(names(x),all_names[ !(all_names %in% names(x))]) 
      x <- res[all_names] }; x
  }, all_names=all_names)
  if(length(text_n)      ==1) text_n       <- rep(text_n,length(tb))
  if(length(text_percent)==1) text_percent <- rep(text_percent,length(tb))
  
  for(i in 1:length(tb)){
    if(any(is.na(names(tb[[i]])))) names(tb[[i]])[is.na(names(tb[[i]]))] <- "NA"
    if(i==1) res <- cbind.data.frame( tb[[i]], round(100*tb[[i]]/sum(tb[[i]],na.rm=T),2) ) 
    else res <- cbind.data.frame(res, tb[[i]], round(100*tb[[i]]/sum(tb[[i]],na.rm=T),2) )
    names(res)[ncol(res)+c(-1,0)] <- c(text_n[[i]],text_percent[[i]])
  }
  for(i in 1:(length(tb)-1))
    for(j in (i+1):length(tb)){
      res <- cbind.data.frame( res, diff = ( res[,2*i-1] - res[,2*j-1] ), proc = round(100 * ( res[,2*i-1] - res[,2*j-1] ) / res[,2*i-1],2) )
      names(res)[ncol(res) +c(-1,0)] <- c(paste0("diff_",i,"_",j), "%")
    }
  if(lorder) res   <- res[order(res[,1],   decreasing = decreasing),]
  res
}



cov_balance <- function(cov_name, dd, categorical=T, tab,
                        ps="ps", swt="swt", wt="wt", treat="case_control",
                        breaks_all=c("Sturges","Scott"), breaks_ps="Sturges",strata_freq=T, strata_dens=T,
                        lprint=F){
  
  for(ibr in breaks_all){
    hist_contr <- hist(dd[dd[,treat]==0,cov_name], ibr,plot=F)
    hist_case  <- hist(dd[dd[,treat]==1,cov_name], ibr,plot=F)
    ylimm <- c(0,max(hist_contr$counts,hist_case$counts))
    plot(hist_contr, col="cyan", xlim=range(dd[,cov_name]), ylim=ylimm, xlab=cov_name,main=cov_name); par(new=T)
    plot(hist_case, col=rgb(t(col2rgb("red"))/255,alpha=0.5),xlim=range(dd[,cov_name]),ylim=ylimm,main="",xlab="",ylab="")
    legend("topright",legend=c("vaccinated","sampled"),fill=c(rgb(t(col2rgb("red"))/255,alpha=0.5),"cyan"))
  }
  
  dd$ps_str <- cut(dd[,ps], unique(quantile(dd[,ps], (0:5)/5 )) )
  for(ips in levels(dd$ps_str)){
    cond <- dd$ps_str ==ips
    hist_contr <- hist(dd[cond & dd[,treat]==0,cov_name],"Sturges", plot=F)   #100,
    hist_case  <- hist(dd[cond & dd[,treat]==1,cov_name],"Sturges", plot=F)   #100,
    
    if(strata_freq){
      ylimm <- c(0,max(hist_contr$counts,hist_case$counts))
      plot(hist_contr, col="cyan", xlim=range(dd[,cov_name]), ylim=ylimm, xlab=cov_name, main=paste0( cov_name," PS in ",ips) ); par(new=T)
      plot(hist_case, col=rgb(t(col2rgb("red"))/255,alpha=0.5),xlim=range(dd[,cov_name]), ylim=ylimm, main="", xlab="", ylab="")
      legend("topright",legend=c("vaccinated","sampled"),fill=c(rgb(t(col2rgb("red"))/255,alpha=0.5),"cyan"))
    }
    if(strata_dens){
      ylimm <- c(0,max(hist_contr$density,hist_case$density))
      plot(hist_contr,freq=F, col="skyblue", xlim=range(dd[,cov_name]), ylim=ylimm, xlab=cov_name, main=paste0( cov_name," PS in ",ips) ); par(new=T)
      plot(hist_case, freq=F,col=rgb(t(col2rgb("red"))/255,alpha=0.5),xlim=range(dd[,cov_name]), ylim=ylimm, main="", xlab="", ylab="")
      legend("topright",legend=c("vaccinated","sampled"),fill=c(rgb(t(col2rgb("red"))/255,alpha=0.5),"skyblue"))
    }
  }
  if(categorical){
    T
  }
  else{
    glm_balance <- lm(get(cov_name) ~ get(treat), data=dd  )
    glm_balance_swt <- lm(get(cov_name) ~ get(treat), data=dd, weights = swt)
    glm_balance_wt <- lm(get(cov_name) ~ get(treat), data=dd, weights =  wt)
    
    if(lprint){
      print(summary(glm_balance))
      print(summary(glm_balance_swt))
      print(summary(glm_balance_wt))
    }
    
    if(!missing(tab))
      for(iwt in c("","swt_","wt_")){
        if(iwt==""    ) glm_res <- glm_balance
        if(iwt=="swt_") glm_res <- glm_balance_swt
        if(iwt=="wt_" ) glm_res <- glm_balance_wt
        tab[cov_name,paste0("diff_",iwt,"value"   )] <- glm_res$coefficients["get(treat)"]
        tab[cov_name,paste0("diff_",iwt,"p_value" )] <- summary(glm_res)$coefficients["get(treat)","Pr(>|t|)"]
        tab[cov_name,paste0("diff_",iwt,"CI_",c("left","right"))] <- confint(glm_res)["get(treat)",]
      }
  }
  if(!missing(tab))
    for(itr in c(1,0)){
      vax_sampl_match_text <- switch(as.character(itr), "0"="sampled_", "1"="vax_")
      for(iwt in c("","swt_","wt_")){
        if(iwt==""    ) var <- dd[ dd[,treat]==itr, cov_name]
        if(iwt=="swt_") var <- dd[ dd[,treat]==itr, cov_name]*dd[ dd[,treat]==itr, swt] / sum(dd[ dd[,treat]==itr, swt]) * sum(dd[,treat]==itr)
        if(iwt=="wt_" ) var <- dd[ dd[,treat]==itr, cov_name]*dd[ dd[,treat]==itr,  wt] / sum(dd[ dd[,treat]==itr,  wt]) * sum(dd[,treat]==itr)
        # not weighted or swt or wt:
        tab[cov_name,paste0(vax_sampl_match_text,iwt,"n_nonmissing")] <- sum(!is.na(var))
        tab[cov_name,paste0(vax_sampl_match_text,iwt,"mean"        )] <- mean(      var, na.rm=T)
        tab[cov_name,paste0(vax_sampl_match_text,iwt,"sd"          )] <- sd(        var, na.rm=T)
        tab[cov_name,paste0(vax_sampl_match_text,iwt,c("min","Q1","median","Q3","max"))] <- quantile(var,c(0,0.25,0.5,0.75,1),na.rm=T)
      }
    }
  tab
} # end of function "cov_balance"





catt <- function(text, log_text=".log_text", start=F){
  if(exists(log_text,envir=.GlobalEnv) & !start) assign(log_text,list( get(log_text,envir=.GlobalEnv), text), envir=.GlobalEnv)
  else assign(log_text, text, envir=.GlobalEnv)
  cat(text)
  invisible()
}

printt <- function(text, log_text=".log_text", start=F){
  if(exists(log_text,envir=.GlobalEnv) & !start) assign(log_text,list( get(log_text,envir=.GlobalEnv), text), envir=.GlobalEnv)
  else assign(log_text, text, envir=.GlobalEnv)
  print(text)
  invisible()
}

warningg <- function(text, log_text=".log_text", start=F){
  if(exists(log_text,envir=.GlobalEnv) & !start) assign(log_text,list( get(log_text,envir=.GlobalEnv), text), envir=.GlobalEnv)
  else assign(log_text, text, envir=.GlobalEnv)
  warning(text)
  invisible()
}

stopp <- function(text, log_text=".log_text", start=F){
  if(exists(log_text,envir=.GlobalEnv) & !start) assign(log_text,list( get(log_text,envir=.GlobalEnv), text), envir=.GlobalEnv)
  else assign(log_text, text, envir=.GlobalEnv)
  stop(text)
  invisible()
}


# Function:     table1
# Description:  create a table with counts and percentages for one categorical variable. Similar to 'tabyl'.
#

table1 <- function(x, title="", digits=2, sep=" & ", print=c(T,T,T,T) ){
  if(!is.null(dim(x)) & length(dim(x))==2){
    for(icol in 2:ncol(x)) x[,1] <- paste(x[,1],x[,icol],sep=sep)
    x <- x[,1]
  }  
  else x <- as.factor(x)
  if(any(print) & title!="") cat(paste(title,"\n"))
  cbind( 
    n=(tb<-table( x, useNA="ifany" )), 
    cum_n=cumsum(tb), 
    percent=round(100*tb/sum(tb),digits), 
    cum_percent=round(100*cumsum( tb/sum(tb) ),digits), 
    percent2=c(round(100*(tb2<-table(x))/sum(tb2),digits),rep(NA,length(tb)-length(tb2))),
    cum_percent2 = c(round(100*cumsum( (tb2<-table(x))/sum(tb2) ),digits),rep(NA,length(tb)-length(tb2))) 
  )[,c( print, any(is.na(x)) & print[3] , any(is.na(x)) & print[4] ), drop=F]
}
#
#
###

