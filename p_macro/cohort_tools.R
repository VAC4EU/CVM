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
                              index_var,   # ="vax_date", 
                              diff_value, diff_unit="days",    # used if methods=="diff"
                              na_value,                        # if NA ==> change to na_value
                              ids_study = ids_study_all, 
                              cov_dataset,   # dataset with variables: id, "date","value_of_variable" 
                              name=paste0("D3_TD_variable_", cov_name,suffix[[subpop]]), dir=paste0(dirtemp,"TD/"),
                              date_var = "date", cov_var = "value_of_variable", create_cov_var=F, event_value=1, 
                              id_original="person_id", id="id_n", 
                              lprint=T, lfreq=F, lplot=F, log_text=".log_text",
                              save_all_events = F){   
  
  catt(paste0("\n",cov_name,"     (",Sys.time(),") ","\n"), log_text=log_text, lprint=lprint)
  methods <- gsub(" ","_",methods) 
  
  if(missing(cov_dataset)){
    if(!file.exists(paste0(dir,name, ".RData"))) stop(paste0("File '",dir,name, ".RData","' is not found"))
    load(paste0(dir,name, ".RData"))  #   "person_id","value_of_variable","date"
    cov_dataset <- as.data.frame(get(name)); rm(list=name)
    if(nrow(cov_dataset)==0) {
      cov_name <- tolower(cov_name)
      for(imethod in methods){
        dd[,paste0(cov_name      ,"_", imethod)] <- NA;
        dd[,paste0(cov_name ,"_date_", imethod)] <- NA;
        if(!missing(na_value)) dd[ is.na(dd[, paste0(cov_name,"_", imethod)]), paste0(cov_name,"_", imethod) ]  <- na_value
      }
      return(dd)
    }
    cov_dataset <- id_add(cov_dataset, data_name=cov_name, ids_study=ids_study, id=id_original, id_new=id, log_text=log_text, lprint=lprint)  
  }
  
  
  if(date_var != "date" | !(date_var %in% names(cov_dataset))){
    if(date_var %in% names(cov_dataset)) cov_dataset[,"date"] <- cov_dataset[,date_var]  
    else stopp(paste("date var is not found in",name), log_text=log_text, lprint=lprint)
  }
  if(   cov_var %in% names(cov_dataset) & cov_var!="value_of_variable") cov_dataset[,"value_of_variable"] <- cov_dataset[,cov_var ]
  if( !(cov_var %in% names(cov_dataset))) {
      if(create_cov_var) cov_dataset[,"value_of_variable"] <- event_value
      else stopp(paste("status var is not found for",cov_name), log_text=log_text, lprint=lprint )
  }
  
  dd <- as.data.frame(dd)
  cov_dataset <- cov_dataset[cov_dataset[,id] %in% dd[!duplicated(dd[,id]),id],]
  cov_dataset[,"date"] <- as.numeric(as.Date(cov_dataset[,"date"] ))
  cov_dataset <- cov_dataset[order(cov_dataset[,id],cov_dataset$date),]
  cov_dataset <- cov_dataset[,c(id,"value_of_variable","date")]
 
  if(any(cond <- duplicated(paste(cov_dataset[,id],cov_dataset$value_of_variable, as.numeric(cov_dataset$date)))) ){
    tmp <- nrow(cov_dataset)
    cov_dataset <- cov_dataset[!cond,]
    catt(paste("Delete double 'id'-'date'-'state' rows:  from", tmp, "rows to    ", nrow(cov_dataset),"rows   ==>   ", tmp[1]-nrow(cov_dataset), "rows deleted.\n"), log_text=log_text, lprint=lprint)
  }
  if(lfreq){
    catt(paste0("\nFrequency for 'value of variable' from dataset '",cov_name,"':\n"), log_text=log_text, lprint=lprint)
    print(table1(cov_dataset[,"value_of_variable"]), log_text=log_text, lprint=lprint); catt("\n", log_text=log_text, lprint=lprint)
  }
  cov_name <- tolower(cov_name)
  
  if(any(is.na(as.numeric(cov_dataset[,"date"]))) ) {
    warning(paste0(sum(is.na(as.numeric(cov_dataset[,"date"])))," rows with missing 'date' deleted from '",cov_name,"'!"), log_text=log_text, lprint=lprint)
    cov_dataset <- cov_dataset[!is.na(as.numeric(cov_dataset[,"date"])),]
  }

  names(cov_dataset)[names(cov_dataset)=="value_of_variable"]  <- cov_name
  if(save_all_events) { 
    all_events <- cov_dataset[,c(id,"date",cov_name)]
    all_events$date <- as.Date(all_events$date,origin="1970-01-01")  # for histogram
    names(all_events)[names(all_events)=="date"] <- paste0(cov_name,"_date")
  }

  dd_small <- dd[ !is.na(as.numeric(as.Date(dd[,index_var]))), c(id,index_var)]
  tmp0 <- merge.data.frame( cov_dataset ,dd_small, by=id, all=F)
 
  for(imethod in methods){
    
    if("last_before" %in% substring(imethod,1,11)) {
      # index_var =""vax_date; "date" - date of event or covariate
      tmp <- tmp0[ as.numeric(tmp0$date) < as.numeric(tmp0[,index_var]) ,  ]
      tmp <- tmp[order(tmp[,id],as.numeric(tmp$date),decreasing=c(F,T)),]
    }
    if("first_after" %in% substring(imethod,1,11)) {
      # index_var =""vax_date; "date" - date of event or covariate
      tmp <- tmp0[ as.numeric(tmp0$date) > as.numeric(tmp0[,index_var]) ,  ]
      tmp <- tmp[order(tmp[,id],as.numeric(tmp$date)),]
    }
    
    tmp <- tmp[!duplicated(tmp[,id]),]
    names(tmp)[names(tmp)=="date"  ] <- paste0(cov_name,"_date_",imethod)
    names(tmp)[names(tmp)==cov_name] <- paste0(cov_name,     "_",imethod)
    
    tmp2 <- nrow(dd)
    dd <- merge.data.frame(dd[, !(names(dd) %in% c("date_tmp",cov_name))], tmp[,c(id, paste0(cov_name,c("_date_","_"),imethod) )], by=id, all.x=T)
    catt(paste("data_nrow=",tmp2, " & cov_data_nrow=",nrow(tmp),"   ==> data_nrow=",nrow(dd), "(",Sys.time(),") \n"), log_text=log_text, lprint=lprint)
   
    dd[,paste0(cov_name,"_date_",imethod)] <- as.Date(dd[,paste0(cov_name,"_date_",imethod)], origin="1970-01-01")
    
  }
  
  
  if(!missing(na_value)) for(imethod in methods) dd[ is.na(dd[, paste0(cov_name,"_", imethod)]), paste0(cov_name,"_", imethod) ]  <- na_value
  
  dd[,index_var] <- as.Date(dd[,index_var], origin="1970-01-01")
  
  
    
  printt(Sys.time(), log_text=log_text, lprint=lprint)
  
  if(save_all_events) res <- list(data= dd, all_events = all_events )
  else res <- dd
  
  res
  
} # end function "add_td_covariate"


smds <- function(vars, vars_star, case_control_var="case_control", data, swt="swt", wt="wt", 
                 var_name, # used inside ==> do not use outside
                 check=F ){ 
  if(!missing(vars_star)) {
    vars0_star <- vars[vars %in% vars_star] 
    vars_star_extra <- vars_star[ !(vars_star %in% vars) ] 
  }
  else {vars0_star<-c(); vars_star_extra <- c()}
  ncol <- 5; lwt <- c(T,F,F)
  if(swt %in% names(data)) {ncol <- ncol + 5; lwt[2] <- T }
  if( wt %in% names(data)) {ncol <- ncol + 5; lwt[3] <- T }
  smd <- matrix(NA, nrow=length(vars), ncol=ncol+3, dimnames=
                  list(vars, c( paste0( rep(c("","swt_","wt_")[lwt],each=5) , rep(c("SMD","mean_case","mean_control", "sd_case","sd_control"),sum(lwt)) ), "n_cases","n_controls", "var_name") ))
  smd <- as.data.frame(smd)
  smd[,"n_cases"     ] <- sum(data[,case_control_var]==1,na.rm=T) 
  smd[,"n_controls"  ] <- sum(data[,case_control_var]==0,na.rm=T)
  if(!missing(var_name)) smd[,"var_name"    ] <- var_name
  else smd[,"var_name"    ] <- vars
  
  irow <- 1

  for(ivar in vars[vars %in% names(data)]){ 
    
    if(mode(data[,ivar])=="numeric"){
      
      for(i in 1:3){ 
        x0 <- data[data[,case_control_var]==0,ivar]
        x1 <- data[data[,case_control_var]==1,ivar]
        pref=""
        if(i==2) if(swt %in% names(data)) {pref <- "swt_"; wti <- swt} else next
        if(i==3) if( wt %in% names(data)) {pref <-  "wt_"; wti <-  wt} else next
        if(i==1){ 
          smd[irow,paste0(pref,"mean_case"   )] <- mean(x1) 
          smd[irow,paste0(pref,"mean_control")] <- mean(x0) 
          smd[irow,paste0(pref,"sd_case"     )] <- sd(x1) 
          smd[irow,paste0(pref,"sd_control"  )] <- sd(x0) 
          smd[irow,paste0(pref,"SMD")] <- ( smd[irow,paste0(pref,"mean_case")] - smd[irow,paste0(pref,"mean_control")] ) / sqrt( (var(x1) + var(x0) ) / 2 ) 
        }
        else{
          wt0 <- data[data[,case_control_var]==0,wti]
          wt1 <- data[data[,case_control_var]==1,wti]
          
          #x0 <- x0 * wt0 / sum(wt0)  * length(x0) 
          #x1 <- x1 * wt1 / sum(wt1) * length(x1) 
          
          smd[irow,paste0(pref,"mean_case"   )] <- sum(wt1*x1)/sum(wt1) 
          smd[irow,paste0(pref,"mean_control")] <- sum(wt0*x0)/sum(wt0)
          var_x1 <- sum( wt1  * (x1 - smd[irow,paste0(pref,"mean_case"   )] )^2 )  /  sum(wt1)
          var_x0 <- sum( wt0  * (x0 - smd[irow,paste0(pref,"mean_control")] )^2 )  /  sum(wt0)
          smd[irow,paste0(pref,"sd_case"     )] <- sqrt(var_x1)
          smd[irow,paste0(pref,"sd_control"  )] <- sqrt(var_x0)
          smd[irow,paste0(pref,"SMD")] <- ( smd[irow,paste0(pref,"mean_case")] - smd[irow,paste0(pref,"mean_control")] ) / sqrt( (var_x1 + var_x0)/2 )  
        }
        
        
        
        if(check & i>1){
          print(i)
          print(round(smd[irow,],5))
          
          x0 <- data[data[,case_control_var]==0,ivar]
          x1 <- data[data[,case_control_var]==1,ivar]
          
          wt0_stand <- data[data[,case_control_var]==0,wti] / sum(data[data[,case_control_var]==0,wti])
          wt1_stand <- data[data[,case_control_var]==1,wti] / sum(data[data[,case_control_var]==1,wti])
          
          print( c( x1_sum_stand  <- sum( x1*wt1_stand) ,
                    x0_sum_stand  <- sum( x0*wt0_stand) ,
                    # sd:
                    x1_sd_stand  = sqrt(sum( (x1 - sum(x1*wt1_stand) )^2 * wt1_stand  )), 
                    x0_sd_stand  = sqrt(sum( (x0 - sum(x0*wt0_stand) )^2 * wt0_stand  ))   ))
          
          #mean:
          print( c( x1_mean <- sum(x1*wt1)/sum(wt1) ,
                    x0_mean <- sum(x0*wt0)/sum(wt0) ,
                    # sd:
                    x1_sd  = sqrt( sum(wt0+wt1)  / (  sum(wt0+wt1)^2 - sum( wt0^2+wt1^2 )  ) * sum( wt1  *  (x1-x1_mean)^2 )  ), 
                    x0_sd  = sqrt( sum(wt0+wt1)  / (  sum(wt0+wt1)^2 - sum( wt0^2+wt1^2 )  ) * sum( wt0  *  (x0-x0_mean)^2 )  ),
                    x1_sd2 = sqrt( sum(wt1)  / (  sum(wt1)^2 - sum( wt1^2 )  ) * sum( wt1  *  (x1-x1_mean)^2 )  ), 
                    x0_sd2 = sqrt( sum(wt0)  / (  sum(wt0)^2 - sum( wt0^2 )  ) * sum( wt0  *  (x0-x0_mean)^2 )  ),
                    x1_sd3 = sqrt( sum( wt1  *  (x1-x1_mean)^2 )  / (((length(x1)-1)*sum(wt1))/length(x1)) ), 
                    x0_sd3 = sqrt( sum( wt0  *  (x0-x0_mean)^2 )  / (((length(x0)-1)*sum(wt0))/length(x0)) ) ))
          print("appr 2")
          #mean:
          print( c( x1_mean <- sum(x1*wt1)/sum(wt1) ,
                    x0_mean <- sum(x0*wt0)/sum(wt0) ,
                    # sd:
                    x1_sd=sqrt( sum( wt1^2 )  /   sum(wt1)^2  * var(x1) ), 
                    x0_sd=sqrt( sum( wt0^2 )  /   sum(wt0)^2  * var(x0) ),
                    x1_sd= sqrt( sum( wt1  * (x1 - sum(x1*wt1)/sum(wt1) )^2 )  /  sum(wt1) ), 
                    x0_sd= sqrt( sum( wt0  * (x0 - sum(x0*wt0)/sum(wt0) )^2 )  /  sum(wt0) )
          ))
          
          
          #mean:
          print("appr 3:")
          wt_3_0 <- data[data[,case_control_var]==0,wti]
          wt_3_1 <- data[data[,case_control_var]==1,wti]
          wt_3_0 <- wt_3_0/sum(wt_3_0)*length(wt_3_0)
          wt_3_1 <- wt_3_1/sum(wt_3_1)*length(wt_3_1)
          
          print( c( x1_mean_3 <- mean(x1*wt_3_1) ,
                    x0_mean_3 <- mean(x0*wt_3_0) ,
                    # sd:
                    x1_sd_3  = sd(x1*wt_3_1) , 
                    x0_sd_3  = sd(x0*wt_3_0) ))

 
          
   glm_balance <- lm(get(ivar) ~ -1 + as.factor(case_control), data=data  )
   glm_balance_swt <- lm(get(ivar) ~ -1 + as.factor(case_control), data=data, weights = data$swt)
   glm_balance_wt <- lm(get(ivar) ~ -1 + as.factor(case_control), data=data, weights =  data$wt)
   
   print(summary(glm_balance))
   print(summary(glm_balance_swt))
   print(summary(glm_balance_wt))
   
   
   glm_balance <- lm(get(ivar) ~ case_control, data=data  )
   glm_balance_swt <- lm(get(ivar) ~ case_control, data=data, weights = data$swt)
   glm_balance_wt <- lm(get(ivar) ~ case_control, data=data, weights =  data$wt)
   
   
            print(summary(glm_balance))
            print(summary(glm_balance_swt))
            print(summary(glm_balance_wt))

        
          
 
        }
      }
      if(ivar %in% vars0_star) dimnames(smd)[[1]][irow] <- paste0( dimnames(smd)[[1]][irow],"*")
      irow <- irow + 1
    }
    
    if(mode(data[,ivar])=="character"){ 
      data[,ivar] <- as.factor(data[,ivar]); cat_names <- levels(data[,ivar])
      if(nlevels(data[,ivar])==1){
        smd[irow,"mean_case"] <- smd[irow,"mean_control"] <- 1
        smd[irow,"sd_case"]   <- smd[irow,"sd_control"]   <- 0
        irow <- irow + 1
      }
      if(nlevels(data[,ivar])>=2){
        for(icat in cat_names[-1]) {
          data[,paste0(ivar,"_",icat)] <- as.numeric(data[,ivar]==icat);  
          if((1:length(cat_names))[icat==cat_names]<length(cat_names))  smd <- rbind( smd[1:irow,,drop=F], smd[irow:nrow(smd),,drop=F] )
          smd[irow,] <- smds(paste0(ivar,"_",icat), data=data, case_control_var=case_control_var, swt=swt, wt=wt, var_name=ivar)
          irow <- irow + 1
        } 
        dimnames(smd)[[1]][irow-(1:(length(cat_names)-1))] <- rev(paste0(ivar,"_",cat_names[-1]))
      }
    }
  }
  attributes(smd) <- c( attributes(smd), vars_in_model = list(vars0_star), extra_vars_in_model = list(vars_star_extra)  ) 
  if(any( !(vars %in% names(data)) )) attributes(smd) <- c( attributes(smd), not_in_data=list(vars[ !(vars %in% names(data)) ]) ) 
  smd
}


plot_smds <- function(x, col_names=c("SMD","wt_SMD","swt_SMD"),
                      col1=c("red","green4","blue"), col2=c("orange","green","skyblue"),pch=c(21,22,24),
                      legend=c("before weighting", "weighted", "weighted stabilised"), 
                      RR, tit_RR="RR = ", round=c(2,8), xlim_min_range=c(-0.1,0.1),
                      tit1="", tit2="in PS model: * ", tit_add_extra=T, tit1_tit2="\n", no_title=F
){ 
  if(!missing(RR)){
    if(length(round)==1) round <- rep(round,2)
    if(is.matrix(RR)) tmp <- cbind(round(RR[,c("exp(coef)","lower .95","upper .95"),drop=F], round[1]), round(RR[,c("Pr(>|z|)"),drop=F],round[2]) )
    else tmp <- cbind(round(summary(RR)$conf.int[,-2,drop=F], round[1]), round(summary(RR)$coefficients[,c("Pr(>|z|)"),drop=F],round[2]))
    dimnames(tmp)[[2]] <- c(""," (",";","); pval=")
    tit_RR <- paste0(tit_RR,paste0(paste0(dimnames(tmp)[[2]],round(tmp,round)),collapse=""))
    tit1 <- paste0(tit1,tit_RR)
  }
  if(length(attr(x,"vars_in_model"))==0) tit2 <- ""
  if(length(attr(x,"extra_vars_in_model"))>0) tit2 <- paste0(c(tit2, paste0(attr(x,"extra_vars_in_model"),collapse = "+")), collapse=" + " )
  
  x <- x[nrow(x):1,]  
  cond3 <- col_names %in% dimnames(x)[[2]]
  x <- x[, match(col_names[cond3],dimnames(x)[[2]]), drop=F]
  xlim1=c( min(x,xlim_min_range,na.rm=T) - 0.3*diff(range(x,xlim_min_range,na.rm=T)),max(x,xlim_min_range,na.rm=T) )
  for(i in 1:sum(cond3)){
    plot(x[,i], 1:nrow(x), axes=F, bg=col2[cond3][i],pch=pch[cond3][i], col=col1[cond3][i],cex=1.5, type="o", xlim=xlim1,xlab="SMD",ylab="")
    if(!no_title) title(paste0(tit1,tit1_tit2,tit2))
    axis(1); box(col="gray")
    text(xlim1[1], 1:nrow(x), dimnames(x)[[1]], adj=0, offset=-3) 
    abline(h=seq(5,100,by=5)+0.5,col="gray")
    abline(v=c(-0.1,0,0.1), col=c("red","gray","red"),lty=2); par(new=T)
  }
  legend("topright",legend=legend[cond3],col=col1[cond3],pt.bg=col2[cond3],pch=pch[cond3],bty="n")
}

summary_var <- function(x, digits=c(3,3), quantiles=T, cut_points){
  if(length(digits)==1) digigts <- rep(digits,2)
  
  res <- list( summary = c( round(summary(x),digits), sd= round(sd(x),digits[1]), n=length(x), n_nonmissing = sum(!is.na(x)) ) )
  if(quantiles){
    if(missing(cut_points)) cut_points <-  unique(c(min(x,na.rm=T)-0.0001, quantile(x,c( seq(0,0.99,by=0.05), seq(0.99,0.999,by=0.001), seq(0.999,1,by=0.0001))) ))
    res <- c(res, quantiles = list(table1(cut(x, cut_points), digits=digits[2])) )  
  }
  res
}

sample_controls <- function( data, 
                             cases,            # name of logical variable (example, T if vax_n=1 & vax_brand=="Pfizer" ==> sampling only for Pfizer vax1)
                             event_date,
                             cond_before_vax,  # list with elements: list( list( var="covid_date_last_before",       period= 30, unit="days", unit_short="d" ),
                             #                           list( var="myocarditis_date_last_before", period=365, unit="days", unit_short="d" )     )
                             vax_date="vax_date", id="id_n", begin_obs="study_entry_date", end_obs="study_exit_date",
                             ids_study, lprint=F, lprint_short=T, log_text=".log_text"
){ 
  sys_time<-Sys.time()
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
        data <- add_td_covariate( data, cond_before_vax[[icond]]$cov_name, index_var="vax_date",
                                  methods=cond_before_vax[[icond]]$method, 
                                  name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                  create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                  ids_study = ids_study, lprint=lprint, log_text=log_text  )    
      
      diff_days <- as.numeric( difftime( data[,vax_date], data[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
      if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(data[,vax_date]) | is.na(data[,cond_before_vax[[icond]]$var])   
      
      data[!cond,end_obs] <- data[!cond,vax_date] - 1
      data[!cond,cases]       <- F 
      
       data <- data[ !data[,cases] | (data[,cases] & cond),  ]
       
    }
  }
  
  
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
  
    
  if(sum(data$days_before_vax)==0) stop("no data for sampling")
  
  n_cases <- sum(data[,cases], na.rm=T)
  catt(paste0( "#days for sampling = ", sum(data$days_before_vax), ",  #cases = ", n_cases, ";   ", Sys.time()), log_text=log_text, lprint= (lprint|lprint_short) )
 
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
                                   index_var="sampled_date",
                                   name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                   create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                   ids_study = ids_study, lprint=lprint, log_text=log_text)    
        if(any(sort(names(cdata))!=sort(names_order))) {print(sort(names_order)); print(sort(names(cdata))); stop("problem: other variable names.")}
        cdata <- cdata[,names_order]
      
        diff_days <- as.numeric( difftime( cdata$sampled_date, cdata[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
        if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(cdata$sampled_date) | is.na(cdata[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(cdata$sampled_date) | is.na(cdata[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(cdata$sampled_date) | is.na(cdata[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(cdata$sampled_date) | is.na(cdata[,cond_before_vax[[icond]]$var])   
        
        cdata <- cdata[cond,  ]
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
  data <- rbind.data.frame(tmp[,names(tmp)],cdata[,names(tmp)])
  
  catt(paste0("==> duration = ",format(difftime(Sys.time(),sys_time))," (till ",Sys.time(),")\n"), log_text=log_text, lprint=lprint_short)  
  
  
  data
}  # end function 'sample_controls'



add_periods <- function(data, data_new, start_state_date, stop_state_date, state_var, id="id_n", start_date="start_date", stop_date="stop_date", lsort=T, log_text=".log_text", lprint=F){ 
  
  state_var <- tolower(state_var)
  data_new[,paste0(state_var,"_new")] <- ""
  data_new_rest <- data_new[,c(id,state_var,paste0(state_var,"_new"),start_state_date,stop_state_date)]
  npat <- max(table(data_new_rest[,id]))
  data_list <- vector("list",length=npat+1)
  data[,state_var] <- ""; data[,paste0(state_var,"_new")] <- "" ; data[,start_state_date] <- NA; data[,stop_state_date] <- NA
  for(i in 1:npat){ 
    catt(paste0("step ",i," (from ",npat,"): ", nrow(data_new_rest)," rows\n"), log_text=log_text, lprint=lprint )
    
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
  printt(table1(data[,state_var]),log_text=log_text,lprint=lprint)
  
  data
  
} # end of function 'add_periods'





hist_vars <- function(vars, legend=names(vars),  
                      tit1="",
                      tit2="", 
                      periods=c(35,100), #c("Sturges","Scott"),
                      smooth_line=F, spar=0.5, df, lwd=2,
                      col=c("yellow","deeppink","skyblue","aquamarine","cornsilk","azure"), alpha=rep(0.6,length(col)),
                      xlim, ylim, xlab="", ylab=ifelse(freq,"Frequency","Density"), freq=T,
                      lplot=T, add_dens=F, bw,
                      separate,   # first all on one plot; and then vars[[1]] with vars[[2]]; vars[[1]] with vars[[3]], ...
                      xlim_separate=F
){  
  
  if(add_dens) freq=F
  
  if(missing(separate)) {
    separate <- list( 1:length(vars))
    if(length(vars)>2) for(ivar in 2:length(vars)) separate <- c(separate, list(c(1,ivar)))
  }
  if(mode(vars[[1]])=="numeric"){
    xrange <- range(sapply(vars, function(x) range(x,na.rm=T) )    )
    hist_res <- vector("list",length=length(periods)); names(hist_res) <- periods
    hist_res$range <- xrange 
    for(iper in as.character(periods)){
      stepp <- diff(xrange)/as.numeric(iper)
      if( abs(stepp - round(stepp)) < 0.25*stepp) stepp <- round(stepp)
      breaks <- c( xrange[1]-0.001*diff(xrange), seq(xrange[1], xrange[2]+stepp, by=stepp) )
      hist_res[[iper]] <- vector("list",length=length(vars)); names(hist_res[[iper]]) <- legend
      for(ivar in 1:length(vars)) {
        hist_res[[iper]][[ivar]] <- hist(vars[[ivar]], breaks,  plot=F)
        if(freq) hist_res[[iper]][[ivar]]$density <- hist_res[[iper]][[ivar]]$counts
        else hist_res[[iper]][[ivar]]$density <- hist_res[[iper]][[ivar]]$counts/sum(hist_res[[iper]][[ivar]]$counts)
        if(add_dens){  
          if(missing(bw)) hist_res[[iper]][[ivar]]$add_dens <- density(vars[[ivar]])#[c("x","y")]
          else hist_res[[iper]][[ivar]]$add_dens <- density(vars[[ivar]],bw=bw)#[c("x","y")]
        }
        hist_res[[iper]][[ivar]] <- hist_res[[iper]][[ivar]][c("breaks","density")]
        
        zero_intervals <- (1:length(hist_res[[iper]][[ivar]]$density))[ hist_res[[iper]][[ivar]]$density==0 ]
        if(length(zero_intervals)>0){
          stop_points <- unique( c( (1:length(zero_intervals))[ diff(zero_intervals)>1 ],length(zero_intervals) ))
          start_points <- c(1, stop_points[-length(stop_points)] + 1 )
          if(length(stop_points)>1){
            start_points <- rev(start_points)
            stop_points  <- rev(stop_points)
            for(i in 1:length(stop_points)) {
              if(length(start_points[i]:stop_points[i])>1) {
                hist_res[[iper]][[ivar]]$breaks  <- hist_res[[iper]][[ivar]]$breaks[   -c( zero_intervals[start_points[i]+1]:zero_intervals[stop_points[i]])    ]
                hist_res[[iper]][[ivar]]$density <- hist_res[[iper]][[ivar]]$density[  -c( zero_intervals[start_points[i]]:zero_intervals[stop_points[i]-1])    ]
              }
            }
          }
        }
        
        
      }

      ymax <- max(sapply(hist_res[[iper]], function(x) max(x$density,na.rm=T) )    )
      for(isep in separate){
        isep_name <- paste0(isep,collapse="_")
        tit=paste0(tit1,tit2); xlab1=xlab; ylab1=ylab     
        
        if(xlim_separate & !identical(isep, 1:length(vars) )) { 
          xlim1 <- range(sapply(vars[isep], function(x) range(x,na.rm=T) )) 
          stepp <- diff(xlim1)/as.numeric(iper)
          if( abs(stepp - round(stepp)) < 0.25*stepp) stepp <- round(stepp)
          breaks1 <- c( xrange[1]-0.001*diff(xrange), seq(xrange[1], xrange[2]+stepp, by=stepp) )
          hist_res[[iper]][[isep_name]] <- vector("list",length=length(isep)); names(hist_res[[iper]][[isep_name]]) <- legend[isep]
          for(ivar in 1:length(isep)){
            hist_res[[iper]][[isep_name]][[ivar]]       <- hist(vars[[isep[ivar]]], breaks1,  plot=F)
            if(freq) hist_res[[iper]][[isep_name]][[ivar]]$density <- hist_res[[iper]][[isep_name]][[ivar]]$counts
            else hist_res[[iper]][[isep_name]][[ivar]]$density <- hist_res[[iper]][[isep_name]][[ivar]]$counts/sum(hist_res[[iper]][[isep_name]][[ivar]]$counts)
            if(add_dens){  
              if(missing(bw)) hist_res[[iper]][[isep_name]][[ivar]]$add_dens <- density(vars[[isep[ivar]]])#[c("x","y")]
              else hist_res[[iper]][[isep_name]][[ivar]]$add_dens <- density(vars[[isep[ivar]]],bw=bw)#[c("x","y")]
            }
            hist_res[[iper]][[isep_name]][[ivar]] <- hist_res[[iper]][[isep_name]][[ivar]][c("breaks","density")]
            

            zero_intervals <- (1:length(hist_res[[iper]][[isep_name]][[ivar]]$density))[ hist_res[[iper]][[isep_name]][[ivar]]$density==0 ]
            if(length(zero_intervals)>0){
              stop_points <- unique( c( (1:length(zero_intervals))[ diff(zero_intervals)>1 ],length(zero_intervals) ))
              start_points <- c(1, stop_points[-length(stop_points)] + 1 )
              if(length(stop_points)>1){
                start_points <- rev(start_points)
                stop_points  <- rev(stop_points)
                for(i in 1:length(stop_points)) {
                  if(length(start_points[i]:stop_points[i])>1) {
                    hist_res[[iper]][[isep_name]][[ivar]]$breaks  <- hist_res[[iper]][[isep_name]][[ivar]]$breaks[   -c( zero_intervals[start_points[i]+1]:zero_intervals[stop_points[i]])    ]
                    hist_res[[iper]][[isep_name]][[ivar]]$density <- hist_res[[iper]][[isep_name]][[ivar]]$density[  -c( zero_intervals[start_points[i]]:zero_intervals[stop_points[i]-1])    ]
                  }
                }
              }
            }
            
          }
          ymax1 <- max(sapply(hist_res[[iper]][[isep_name]], function(x) max(x$density,na.rm=T) )    )
          hist_res[[iper]][[isep_name]]$range <- xlim1 
        }
        else {xlim1 <- xrange ; ymax1 <- ymax}
        if(lplot) { 
          for(ivar in isep){
            if(xlim_separate & !identical(isep, 1:length(vars) )) hist_res1 <- hist_res[[iper]][[isep_name]][[match(ivar,isep)]] 
            else hist_res1 <- hist_res[[iper]][[ivar]] 
   
            attributes(hist_res1) <- c(attributes(hist_res1), class="histogram")
            plot(hist_res1, xlim=xlim1, ylim=c(0,ymax1), xlab=xlab1,ylab=ylab1,main=tit, col=rgb(t(col2rgb(col[ivar])/255),alpha= alpha[ivar])); par(new=T)
            if(add_dens   ) { try(lines(hist_res1$add_dens,col=col[ivar], lwd=lwd))  ; par(new=T)  }        
            if(smooth_line) { 
              if(!missing(spar)) {try(lines(smooth.spline(hist_res1$mids,hist_res1$density,spar=spar),col=col[ivar], lwd=lwd))  ; par(new=T)  }
              else               {try(lines(smooth.spline(hist_res1$mids,hist_res1$density,df  =df  ),col=col[ivar], lwd=lwd))  ; par(new=T)  }   
            }
            tit <- ""; xlab1=""; ylab1=""
          }
          legend("topright", legend=legend[isep], fill=rgb(t(col2rgb(col[isep])/255),alpha= alpha[isep]), bty="n")
        }
      }
    }
    res <-   hist_res
  }
  if(mode(vars[[1]])=="character"){
    
    hist_res <- vector("list",length=length(periods)); names(hist_res) <- periods
    
    hist_res <- vector("list",length=length(vars)); names(hist_res) <- legend
    for(ivar in 1:length(vars)) hist_res[[ivar]] <- table(vars[[ivar]])
    
    if(!freq) for(ivar in 1:length(vars)) hist_res[[ivar]] <- hist_res[[ivar]] / sum(hist_res[[ivar]])
    ymax <- max(sapply(hist_res, function(x) max(x,na.rm=T) )    )
    all_names <- unlist( lapply(hist_res, function(x)names(x) ))
    all_names <- all_names[!duplicated(all_names)]
    all_data <- matrix(NA, nrow=length(all_names), ncol=length(vars), dimnames=list(all_names,names(vars)))
    for(ivar in 1:length(vars)) all_data[ match( names(hist_res[[ivar]]), dimnames(all_data)[[1]] ), ivar] <- hist_res[[ivar]]
    if(lplot){ 
      barplot( all_data,col=rep( rgb(t( ( col2rgb(col[1:nrow(all_data)]) * rep((nrow(all_data):1)/nrow(all_data),each=3) )/255)),length(vars)),beside=T ); title(paste0(tit1,tit2)) 
      if(nrow(all_data)>2) legend("topright", legend=dimnames(all_data)[[1]], fill=rgb(t( ( col2rgb(col[1:nrow(all_data)]) * rep((nrow(all_data):1)/nrow(all_data),each=3) )/255)),bty="n")
      barplot( t(all_data),col=rep(rgb(t(col2rgb(col[1:length(vars)])/255)),nrow(all_data)),beside=T ); title(paste0(tit1,tit2)) 
      if(length(vars)>2) legend("topright", legend=legend[1:length(vars)], fill=rgb(t(col2rgb(col[1:length(vars)])/255)),bty="n")
    }
    res <- all_data
  }
  
  invisible(res)
} 





plot_hist_vars <- function( hist_res, legend=names(hist_res[[1]][ sapply(hist_res[[1]],function(x)"breaks" %in% names(x)) ]),  
                      separate,   # first all on one plot; and then vars[[1]] with vars[[2]]; vars[[1]] with vars[[3]], ...
                      tit1="",
                      tit2="", 
                      xlim, ylim, xlab="", ylab="", 
                      col=c("yellow","deeppink","skyblue","aquamarine","cornsilk","azure"), alpha=rep(0.6,length(col)),
                      smooth_line=F, spar=0.5, df, lwd=2,
                      add_dens=F, bw,
                      periods=names(hist_res[-length(hist_res)])
){  
  

  if(is.matrix(hist_res)){
    legend <- rownames(hist_res)
    barplot( hist_res,col=rep( rgb(t( ( col2rgb(col[1:nrow(hist_res)]) * rep((nrow(hist_res):1)/nrow(hist_res),each=3) )/255)),nrow(hist_res)),beside=T ); title(paste0(tit1,tit2)) 
    if(nrow(hist_res)>2 | ncol(hist_res)>2) 
      legend("topright", legend=dimnames(hist_res)[[1]], fill=rgb(t( ( col2rgb(col[1:nrow(hist_res)]) * rep((nrow(hist_res):1)/nrow(hist_res),each=3) )/255)),bty="n")
    barplot( t(hist_res),col=rep(rgb(t(col2rgb(col[1:ncol(hist_res)])/255)),ncol(hist_res)),beside=T ); title(paste0(tit1,tit2)) 
    if(nrow(hist_res)>2 | ncol(hist_res)>2) 
      legend("topright", legend=legend[1:ncol(hist_res)], fill=rgb(t(col2rgb(col[1:ncol(hist_res)])/255)),bty="n")
    return(invisible(NULL))
  }
  

  nvars    <- sum( sapply(hist_res[[1]],function(x)"breaks" %in% names(x)) )
  xlim_separate <- nvars < length(hist_res[[1]])
  
  if(missing(separate)){
    separate <- as.list(names(hist_res[[1]][ !sapply(hist_res[[1]],function(x)"breaks" %in% names(x)) ] ))
    if(length(separate)>0) for(i in 1:length(separate)) separate[[i]] <- as.numeric(strsplit(separate[[i]],"_")[[1]])
    separate <- c(list(c(1:nvars)), separate)
  }
  
  for(iper in as.character(periods)){
      xrange <- hist_res$range
      ymax   <- max(  sapply(hist_res[[iper]][1:nvars], function(x) max(  x$density,na.rm=T) )    )
      for(isep in separate){
        isep_name <- paste0(isep,collapse="_")
        tit=paste0(tit1,tit2); xlab1=xlab; ylab1=ylab     
        
        if(xlim_separate & !identical(isep, 1:nvars )) {
          xlim1 <- hist_res[[iper]][[isep_name]]$range
          ymax1 <- max(  sapply(hist_res[[iper]][[isep_name]][(1:length(isep))], function(x) max(  x$density,na.rm=T) )    )
        }
        else {xlim1 <- xrange ; ymax1 <- ymax}
        for(ivar in isep){
            if(xlim_separate & !identical(isep, 1:nvars )) hist_res1 <- hist_res[[iper]][[isep_name]][[match(ivar,isep)]] 
            else hist_res1 <- hist_res[[iper]][[ivar]] 
            
            attributes(hist_res1) <- c( attributes(hist_res1), class="histogram")
            plot(hist_res1, xlim=xlim1, ylim=c(0,ymax1), xlab=xlab1,ylab=ylab1,main=tit, col=rgb(t(col2rgb(col[ivar])/255),alpha= alpha[ivar])); par(new=T)
            if(add_dens   ) { try(lines(hist_res1$add_dens,col=col[ivar], lwd=lwd))  ; par(new=T)  }        
            if(smooth_line) { 
              if(!missing(spar)) {try(lines(smooth.spline(hist_res1$breaks[-1] - diff(hist_res1$breaks)/2,hist_res1$density,spar=spar),col=col[ivar], lwd=lwd))  ; par(new=T)  }
              else               {try(lines(smooth.spline(hist_res1$breaks[-1] - diff(hist_res1$breaks)/2,hist_res1$density,df  =df  ),col=col[ivar], lwd=lwd))  ; par(new=T)  }   
            }
            tit <- ""; xlab1=""; ylab1=""
        }
        legend("topright", legend=legend[isep], fill=rgb(t(col2rgb(col[isep])/255),alpha= alpha[isep]), bty="n")
      }
    }
  invisible(NULL)
} # end plot_hist_var





hist_dates <- function(vars, legend=names(vars),  
                      tit1="",
                      tit2="", 
                      periods=c(7,1,365/12),
                      smooth_line=F,
                      col=c("yellow","deeppink","skyblue","aquamarine","cornsilk","azure"), alpha=0.6,
                      xlim, ylim, xlab="Time", ylab="Frequency",
                      vline=T,
                      lplot=T){  
  
  vertical_lines <- as.numeric(difftime( as.Date(paste0("20",rep(18:24,each=4),"-",rep(c("01","04","07","10"),7),"-01")),as.Date("2020-08-31")),units="days")
  cal_time_range_days <- as.numeric(difftime( as.Date(range( unlist(vars),na.rm=T),origin="1970-01-01"), as.Date("2020-08-31") ,units="days"))
  vars_days <- lapply(vars, function(x) as.numeric(difftime( x, as.Date("2020-08-31"), units="days")) )                                  
  lplot <- rep(lplot,2)

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
        if(vline) abline(v=vertical_lines,col="lightgray")
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




plot_hist_dates <- function(hist_res, legend=names(hist_res[[1]][ sapply(hist_res[[1]],function(x)"breaks" %in% names(x)) ]),  
                       tit1="",
                       tit2="", 
                       periods=names(hist_res[-length(hist_res)]),
                       smooth_line=F,
                       col=c("yellow","deeppink","skyblue","aquamarine","cornsilk","azure"), alpha=0.6,
                       xlim, ylim, xlab="Time", ylab="Frequency",
                       vline=T,
                       lplot=T){  

  vertical_lines <- as.numeric(difftime( as.Date(paste0("20",rep(18:24,each=4),"-",rep(c("01","04","07","10"),7),"-01")),as.Date("2020-08-31")),units="days")
  cal_time_range_days <- as.numeric(difftime( as.Date(range( unlist(vars),na.rm=T),origin="1970-01-01"), as.Date("2020-08-31") ,units="days"))
  vars_days <- lapply(vars, function(x) as.numeric(difftime( x, as.Date("2020-08-31"), units="days")) )                                  
  lplot <- rep(lplot,2)

  res<-vector("list",length=length(periods)); names(res) <- paste0("period_",periods,"days")
  for(iper in as.character(periods)) {
    
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
        if(vline) abline(v=vertical_lines,col="lightgray")
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


id_add <- function(data, ids_study, id, id_new="id_n", data_name="", log_text=".log_text", lprint=T){
  if( any(!(tmp_ids<-(data[!duplicated(data[,id]),id] %in% ids_study))) ){ 
    tmp_ids <- sum(!tmp_ids,na.rm=T)
    tmp_dim<-dim(data)
    data <- data[data[,id] %in% ids_study,]
    n_rows_deleted <- tmp_dim[1]-nrow(data)
  }
  else tmp_ids  <- n_rows_deleted <- 0
  n_ids         <- sum(!duplicated(data[,id]),na.rm=T)
  data[,id_new] <- match(data[,id],ids_study)
  ids_info      <- c(n_ids=n_ids, n_ids_deleted=tmp_ids, n_rows_deleted=n_rows_deleted )
  
  catt(paste0("dataset:'",data_name,"':     #id's in the study population = ",n_ids," (",nrow(data)," rows)  + ", tmp_ids," id's are not in the study population (", n_rows_deleted," rows deleted).\n\n"),log_text=log_text,lprint=lprint)
  
  attributes(data) <- c(attributes(data)[ !(names(attributes(data)) %in% c("data_name","ids_info")) ], data_name=data_name, ids_info = list(ids_info ))
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
                      ids_study, id="id_n",
                      vax_vars=c("vax_brand","vax_n"), 
                      start_interval="study_entry_days", stop_interval="study_exit_days", death_time="death_days",
                      result_name="", file_name="", dir="", nboot=0, log_text=".log_text", log_dir="", lprint=F, lprint_str=T, lprint_warn=T, lprint_all=F, create_i=T, lparal=F){
 
  printt(Sys.time(),log_text=log_text, lprint= (lprint_str | lprint) )
  
  if(missing(cases)){ cases <- "_tmp_cases";dd0[,"_tmp_cases"] <- !is.na(dd0[,vax_date]) }
  cond <- !dd0[,cases] & !is.na(dd0[,vax_date]) &  dd0[,vax_date]  <= dd0[,stop_interval]
  dd0[  cond, stop_interval ] <- dd0[  cond, vax_date ] - 1
  
  # for vaccinated:
  if(!missing(cond_before_vax)){
    if(missing(ids_study)) stopp(paste("'ids_study' must be specified for sampling"),log_text=log_text, save=T, dir=log_dir)
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
        dd0 <- add_td_covariate( dd0, cond_before_vax[[icond]]$cov_name,  index_var="vax_date",
                                 methods=cond_before_vax[[icond]]$method, 
                                  cov_dataset = cond_before_vax[[icond]]$dataset,
                                  #name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                  create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                  ids_study = ids_study, lprint=F, log_text=log_text  )    
      
      if(cond_before_vax[[icond]]$condition!=""){
        diff_days <- as.numeric( difftime( dd0[,vax_date], dd0[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
        if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
        if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(dd0[,vax_date]) | is.na(dd0[,cond_before_vax[[icond]]$var])   
      }
      dd0[!cond,stop_interval] <- dd0[!cond,vax_date] - 1
      if(any(!cond)) dd0[!cond,cases]         <- F 
      
      dd0 <- dd0[ !dd0[,cases] | (dd0[,cases] & cond),  ]
       
    }
  }
  
  nboot <- nboot+1
  start_interval_all <- dd0[,start_interval]
  end_obs_all <- dd0[,stop_interval]
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
  printt(table1(dd0$matching_factor),log_text=log_text)
  matching_factor_levels <- levels(as.factor(dd0$matching_factor))
  dd0$matching_factor <- as.numeric(as.factor(dd0$matching_factor))
  
  
  # N3: 
  ####### with strata:                   ==>  31 min 33 sec
  (start_time <- Sys.time())
  
  if(sum(is.na(dd0$matching_factor))>0) catt(paste0("\n",sum(is.na(dd0$matching_factor))," rows with missing values in matching variables.\n\n"),log_text)
  
  if(!lparal){
    str_tab <- table(dd0[cond_vax,"matching_factor"]) #; str_tab <- str_tab[!is.na(str_tab)]
    str_tab_len <- length(str_tab)
    if(length(matching_indep_factors)>0){
      catt(paste0( 'matching factors: "',paste0(matching_factor_levels,collapse='","'),'"\n' ),log_text=log_text,lprint=lprint)
      printt(c( str_tab, Sum=sum(str_tab) ),log_text=log_text,lprint=lprint)
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
          catt(paste0("stratum ", istr_i," (",istr, ": '", matching_factor_levels[as.numeric(istr)],"') from ",str_tab_len," (",str_case_nlevels,");      "),
               log_text=log_text,lprint=(lprint_str & lprint ) ) #lprint=(lprint_str & ( lprint | ( !lprint & istr_i %% 10 == 0 ) )) ) 
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
      
      
      catt(paste0(str_text, ":   #case-rows = ",length(icasedata_str), ";  number of days = ",length(exposure_days_unique_str), ";   ", Sys.time(), "\n"),
            log_text=log_text, lprint = lprint_str |  lprint  ) # lprint = lprint_str |  lprint | ( !lprint & istr_i %% 10 == 0))
      
      
      
 
      
      
      
      
      
      if(!missing(matching_dep_factors)){
      #if(!missing(matching_dep_factors) & any(cond_contr)){
          
        dep_matrix <- sapply( matching_dep_factors, function(x) {
          if(length(x)==4) if(is.null(names(x))) names(x) <- c("date","event","abs_diff","units") else x<-x[c("date","event","abs_diff","units")] 
          if(length(x)==3){ x <- c(x[1], gsub("_date","",x[1]), x[2:3]); names(x) <- c("date","event","abs_diff","units")}
          if(length(x)<3) stopp("vector for 'matching_dep_factors' should have length 3 of 4 ('date','event','abs_diff','units') or ('date', 'abs_diff','units')", log_text=log_text, save=T, dir=log_dir)
          x
        } )
      iday=1; ivaxday=exposure_days_unique_str[iday]     
        # cases for day 'ivaxday' 
#?        csdata0 <- dd0[ icasedata_str[vax_days_str==ivaxday], names(dd0) %in%c("i",id, dep_matrix[c("date","event"),]), drop=F ]

      csdata0 <- dd0[ , names(dd0) %in%c("i",id, dep_matrix[c("event"),]), drop=F ]
      for(i in 1:length(cond_before_vax)) if( !(dep_matrix["event",i] %in% names(csdata0)) ) csdata0[,dep_matrix["event",i]] <- 1
      
      non_empty_case_events_cat <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix[c("event"),, drop=F]]) )>0 , dep_matrix[c("event"),] ]
      non_empty_case_events_cat <- non_empty_case_events_cat[!duplicated(non_empty_case_events_cat),]
      
      
      
  
      csdata0 <- dd0[ dd0[,start_interval] <= dd0[,vax_date] & dd0[,vax_date]<= dd0[,stop_interval], names(dd0) %in% c("i",id, vax_date, start_interval, stop_interval,  dep_matrix[c("date","event"),]), drop=F ]
      for(i in 1:length(cond_before_vax)) if( !(dep_matrix["event",i] %in% names(csdata0)) ) csdata0[,dep_matrix["event",i]] <- 1
      
      non_empty_case_events <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix[c("event","date"),, drop=F]]) )>0 , c("vax_date",start_interval, stop_interval, dep_matrix[c("event","date"),]) ]
      #non_empty_case_events <- csdata0[ rowSums( !is.na(csdata0[, dep_matrix["date",, drop=F]]) )>0 , dep_matrix["event",] ]
      non_empty_case_events <- non_empty_case_events[!duplicated(non_empty_case_events),]
      




        
        
        
        
        if(nrow(non_empty_case_events_cat)>0){ 
          
#          print(cbind.data.frame( day=iday, from=length(exposure_days_unique_str), non_empty_case_events, as.character(Sys.time())))
          
          for(irow in 1:nrow(non_empty_case_events_cat)){  #print(cbind.data.frame(irow,non_empty_case_events[irow,],as.character(Sys.time())) )# per each combi in cases
            
            cond_contr2 <- rep(T,nrow(dd0))
            for(idep in 1:ncol(non_empty_case_events_cat))
              if(dep_matrix[c("event"),idep] %in% names(dd0)) 
                cond_contr2 <- cond_contr2 & dd0[,dep_matrix[c("event"),idep],] == non_empty_case_events_cat[irow, dep_matrix[c("event"),idep]]
            
            cond_contr_short <- rep(T,nrow(non_empty_case_events))
            for(idep in 1:ncol(non_empty_case_events_cat))
                cond_contr_short <- cond_contr_short & non_empty_case_events[,dep_matrix[c("event"),idep]] == non_empty_case_events_cat[irow, dep_matrix[c("event"),idep]]
            non_empty_case_events_irow <- non_empty_case_events[cond_contr_short,]
            
            contr_available <-  dd0[cond_contr2, ]
            
            
            
            
            
            
            
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
        
        
        # for unvaccinated periods or persons
        # select id's with conditions before sampled date (example: nocovid_before_vax_days during 30 days, or no myocarditits before vax during 1 year)
        if(!missing(cond_before_vax) ){
        #if(!missing(cond_before_vax) & any(cond_contr)){
          for(icond in 1:length(cond_before_vax)){
            cdata <- dd0[ match( idata_str,dd0$i), c(id, vax_date)]
            cdata <- cdata[cond_contr,]
            cdata$ii <- 1:nrow(cdata); nrow0 <- nrow(cdata)
            if(nrow(cond_before_vax[[icond]]$dataset)==0) next
            cdata <- add_td_covariate( cdata[, !(names(cdata) %in% paste0(cond_before_vax[[icond]]$cov_name,c("_date_","_"),cond_before_vax[[icond]]$method)) ], 
                                       cond_before_vax[[icond]]$cov_name, index_var="???_date",
                                       methods=cond_before_vax[[icond]]$method,
                                       cov_dataset = cond_before_vax[[icond]]$dataset,
                                       #name=cond_before_vax[[icond]]$name, dir=cond_before_vax[[icond]]$dir,  
                                       create_cov_var = cond_before_vax[[icond]]$create_cov_var,
                                       ids_study = ids_study, lprint=F, log_text=log_text  )    
            #if(any(sort(names(cdata))!=sort(names_order))) {print(sort(names_order)); print(sort(names(cdata))); stop("problem: other variable names.")}
            if(nrow(cdata)!= nrow0) {stopp("problem: # rows should be the same.",log_text=log_text,save=T,dir=log_dir)} 
            cdata <- cdata[cdata$ii,]
            
            if(cond_before_vax[[icond]]$condition!=""){
              #???
              diff_days <- as.numeric( difftime( as.Date(ivaxday,"2020-08-31"), cdata[,cond_before_vax[[icond]]$var],  units=cond_before_vax[[icond]]$unit))
              if(cond_before_vax[[icond]]$condition==">" ) cond <- diff_days >  cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition==">=") cond <- diff_days >= cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition=="<" ) cond <- diff_days <  cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])   
              if(cond_before_vax[[icond]]$condition=="<=") cond <- diff_days <= cond_before_vax[[icond]]$period | is.na(cdata[,vax_date]) | is.na(cdata[,cond_before_vax[[icond]]$var])
            }
            #else
            cond_contr[cond_contr] <- cond
           }
        } # end if: !missing(cond_before_vax) & any(cond_contr) 


            if(ncases > ncontr) {
              catt(paste0("Less potential controls for time = ",ivaxday," (day=", iday,");  #cases = ",ncases,";  #controls = ",ncontr, str_text," & ",apply( non_empty_case_events[ irow,,drop=F], 1, paste,collapse=" "), "\n"),
                   log_text=log_text, lprint=lprint_warn)  
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
          
        } # if nrow(non_empty_case_events)>0
         
      } # end if !missing(matching_dep_factors) & any(cond_contr)







      if(missing(matching_dep_factors)){ # end if missing(matching_dep_factors) & any(cond_contr)
      #if(missing(matching_dep_factors) & any(cond_contr)){ # end if missing(matching_dep_factors) & any(cond_contr)
          
        
        
        for(iday in 1:length(exposure_days_unique_str)){
          ivaxday <- exposure_days_unique_str[iday]
          cond_contr <- start_interval_str<=ivaxday  & ivaxday<=end_obs_all_str # & !( !is.na(event_str) & event_str<ivaxday ) 
          cond_contr <- cond_contr &  !( !is.na(vax1_days_all_str) & vax1_days_all_str <= ivaxday )  
          
 
        ncontr  <- sum(cond_contr)
        icases_str <- icasedata_str[vax_days_str==ivaxday]
        ncases  <- length(icases_str)
        
        
        
        #catt(paste0(  iday, " "), log_text=log_text, lprint=(lprint_str & iday%%10==0 & !lprint_all) )
        #catt(paste0( "\nday = ", iday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"), log_text=log_text, lprint=(lprint_str & iday%%50==0 & !lprint_all) )
        catt(paste0( "day=", iday, ";    time = ",ivaxday, ";   n_cases = ",ncases, ";   available_controls = ",ncontr, ";        ", format(Sys.time()),"\n"), log_text=log_text, lprint=lprint_all)
        
        if(ncases > ncontr) {
          catt(paste0("Less potential controls for time = ",ivaxday," (day=", iday,");  #cases = ",ncases,";  #controls = ",ncontr, str_text, "\n"), log_text=log_text, lprint=lprint_warn)  
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
       
        }  # end "for" iday         
      } # end of missing(matching_dep_factors) & any(cond_contr)

          
        
        n_day_matched_cases_str      <- sapply(matched_cases_list_str,length)
        n_day_matched_controls_str   <- unlist(lapply(matched_controls_list_str, nrow))
        #n_day_less_controls_list_str <- sapply(less_controls_list_str,function(x)length(x$unmatched_cases))
        
        if(!all(n_day_matched_cases_str[n_day_matched_cases_str>0]==n_day_matched_controls_str[n_day_matched_controls_str>0])) { printt("problem 1.1", log_text=log_text);}  
        
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
  
  
  printt(Sys.time()-start_time,log_text=log_text); printt(Sys.time(),log_text=log_text)
  
  
  if(length(matching_indep_factors)==0){matching_factor_levels <- c()}
  if(!missing(cond_before_vax)) cond_before_vax <- lapply(cond_before_vax,function(x)x[names(x)!="dataset"])
  else cond_before_vax <- list()
  
  
  match_res <- list( matched_cases=matched_cases, matched_controls=matched_controls, 
                     less_controls_matched_cases=less_controls_matched_cases, less_controls_unmatched_cases=less_controls_unmatched_cases, 
                     cond_before_vax = cond_before_vax,
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
  
  printt(Sys.time(), log_text=log_text, lprint= (lprint_str | lprint) )
  
  match_res 
  
} # end of function 'match_pop'





get_matched_dataset <- function(matched_list, data_vax, update_event=F, stop_obs="study_exit_days", next_vax_time="", iboot=0, create_i=F, log_text=".log_text",dir="", lprint=T){
 
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
  catt(paste0("#rows in dataset = ",nrow(data_vax),";   #rows with vax date = ",sum(!is.na(data_vax[,vax_time])), " (",round(100*sum(!is.na(data_vax[,vax_time]))/nrow(data_vax),2),"%);  " ), log_text=log_text, lprint=lprint )
  tmp <- data_vax[matched_list[["matched_cases"]],]
  catt(paste0( "#matched cases = ", nrow(tmp)," (",round(100*nrow(tmp)/sum(!is.na(data_vax[,vax_time])),2),"%)     and      "), log_text=log_text, lprint=lprint)
  
  tmp$i_vaxed   <- matched_list[["matched_cases"]]
  tmp$i_unvaxed <- matched_list[["matched_controls"]][[paste0("boot_",iboot)]]
  tmp[,paste0(vax_vars,"_case")] <- tmp[,vax_vars]
  tmp$case_control <- 1
  
  # controls:
  data_vax <- data_vax[matched_list[["matched_controls"]][[paste0("boot_",iboot)]],]
  catt(paste0("#matched controls = ",nrow( data_vax),"\n"), log_text=log_text, lprint=lprint)
  
  data_vax$i_vaxed   <- matched_list[["matched_cases"]]
  data_vax$i_unvaxed <- matched_list[["matched_controls"]][[paste0("boot_",iboot)]]
  data_vax[,vax_time] <- tmp[,vax_time]
  data_vax[,paste0(vax_vars,"_case")] <- tmp[,vax_vars]
  data_vax$case_control <- 0
  
  if(nrow(tmp)!=nrow(data_vax)) stopp("#cases and #controls should be the same!", log_text=log_text, save=T, dir=dir)
  
  data_vax <- rbind.data.frame(tmp, data_vax)
  
  
  if(update_event) 
    if(event_time!="" & event_time %in% names(data_vax)){    
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



balance_one <- function(cov_name, dd, tab,
                        ps="ps", swt="swt", wt="wt", treat="case_control",
                        strata_freq=T, strata_dens=T, lprint=F ){
    
    glm_balance <- lm(get(cov_name) ~ get(treat), data=dd  )
    if(swt %in% names(dd)) glm_balance_swt <- lm(get(cov_name) ~ get(treat), data=dd, weights = swt)
    if( wt %in% names(dd)) glm_balance_wt <- lm(get(cov_name) ~ get(treat), data=dd, weights =  wt)
    
    if(lprint){
      print(summary(glm_balance))
      if(swt %in% names(dd)) print(summary(glm_balance_swt))
      if( wt %in% names(dd)) print(summary(glm_balance_wt))
    }
    
    if(!missing(tab))
      for(iwt in c("","swt_","wt_")){
        if(iwt=="swt_" & !(swt %in% names(dd)) ) next  
        if(iwt== "wt_" & !( wt %in% names(dd)) ) next  
        if(iwt==""    ) glm_res <- glm_balance
        if(iwt=="swt_") glm_res <- glm_balance_swt
        if(iwt=="wt_" ) glm_res <- glm_balance_wt
        tab[cov_name,paste0("diff_",iwt,"value"   )] <- glm_res$coefficients["get(treat)"]
        tab[cov_name,paste0("diff_",iwt,"p_value" )] <- summary(glm_res)$coefficients["get(treat)","Pr(>|t|)"]
        tab[cov_name,paste0("diff_",iwt,"CI_",c("left","right"))] <- confint(glm_res)["get(treat)",]
      }
    
    if(!missing(tab))
      for(itr in c(1,0)){
        vax_sampl_match_text <- switch(as.character(itr), "0"="sampled_", "1"="vax_")
        for(iwt in c("","swt_","wt_")){
          if(iwt=="swt_" & !(swt %in% names(dd)) ) next  
          if(iwt== "wt_" & !( wt %in% names(dd)) ) next  
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
} # end of function "balance_one"

balance <- function(vars, dd, 
                    ps="ps", swt="swt", wt="wt", treat="case_control",
                    strata_freq=T, strata_dens=T ){ 
  
  hist_res_all <- list()
  for(icov in vars[vars %in% names(dd)]){
    
    
    hist_res <- hist_vars(list( controls = dd[dd[,treat]==0,icov],
                                cases     = dd[dd[,treat]==1,icov]),  lplot=F)
    #if(lplot) plot_hist_vars(  hist_res  , xlab=icov,tit1="Approach: Sampling" )
    hist_res_one <- c(hist_all = list(hist_res))
    names(hist_res_one) <- icov
    
    if(ps %in% names(dd)){
      dd$ps_str <- cut(round(dd[,ps],2), unique(quantile(dd[,ps], (0:5)/5 )) )
      for(ips in 1:nlevels(dd$ps_str)){
        cond <- dd$ps_str ==levels(dd$ps_str)[ips]
        
        if(strata_freq){
          hist_res <- hist_vars(list( controls = dd[cond & dd[,treat]==0,icov],
                                      cases    = dd[cond & dd[,treat]==1,icov]),  lplot=F)
          
          #if(lplot) plot_hist_vars(  hist_res  , xlab=icov,tit2="\nApproach: Sampling" , tit1=paste0( icov,":  PS in ",levels(dd$ps_str)[ips]) )
          hist_res_one <- c(hist_res_one, list(hist_res)); 
          names(hist_res_one) <- c(icov,paste0( "PS_in_",levels(dd$ps_str)[1:ips]))
        }
      }
    }
    hist_res_all <- c(hist_res_all, icov=list(hist_res_one)); names(hist_res_all)[length(hist_res_all)] <- icov
    
    
    
    
    cat(paste0(icov," "))
    if(mode(dd[,icov])=="numeric"  ){
      tab_one <- matrix(NA,nrow=length(icov),ncol=9*6+4*3, 
                        dimnames=list(icov,
                                      c( paste0(rep(c("vax_","sampled_","vax_swt_","sampled_swt_","vax_wt_","sampled_wt_"),each=9), 
                                                rep( c("n","n_nonmissing","mean","sd","min","Q1","median","Q3","max"),  6)),
                                         paste0(rep(c("diff_","diff_swt_","diff_wt_"),each=4), rep(c("value","CI_left","CI_right","p_value"),3) )
                                      )  ) )
      tab_one <- balance_one(icov, dd,  tab=tab_one )
      if("tab_res" %in% ls()) tab_res <- rbind.data.frame(tab_res, tab_one)  else tab_res <- tab_one
      #hist_all <- c(hist_all,  list(tab_one$hist) ); names(hist_all)[length(hist_all)] <-  icov
    }
    if(mode(dd[,icov])=="character") if( nlevels(as.factor(dd[,icov]))>1)
      for(il in 2:nlevels(as.factor(dd[,icov])) ){
        new_name <- paste0(icov,"_",levels(as.factor(dd[,icov]))[il])
        dd[,new_name] <- as.numeric(as.factor(dd[,icov])==levels(as.factor(dd[,icov]))[il] )
        tab_one <- matrix(NA,nrow=length(new_name),ncol=9*6+4*3, 
                          dimnames=list(new_name,
                                        c( paste0(rep(c("vax_","sampled_","vax_swt_","sampled_swt_","vax_wt_","sampled_wt_"),each=9), 
                                                  rep( c("n","n_nonmissing","mean","sd","min","Q1","median","Q3","max"),  6)),
                                           paste0(rep(c("diff_","diff_swt_","diff_wt_"),each=4), rep(c("value","CI_left","CI_right","p_value"),3) )
                                        )  ) )
        tab_one <- balance_one(new_name, dd, tab=tab_one )
        if("tab_res" %in% ls()) tab_res <- rbind.data.frame(tab_res, tab_one)  else tab_res <- tab_one
        # hist_all <- c(hist_all,  list(tab_one$hist) ); names(hist_all)[length(hist_all)] <-  new_name
      }        
  }; cat("\n\n")
  list( balance_tab=tab_res, hist=hist_res_all )
}

plot_balance <- function(balance_res, mfrow,mfcol, elements=1:length(balance_res$hist) ){ 
  if(!missing(mfrow)) par(mfrow=mfrow)
  if(!missing(mfcol)) par(mfcol=mfcol)
  for(icov in elements)
    for(ih in 1:length(balance_res$hist[[icov]]))
      plot_hist_vars(  balance_res$hist[[icov]][[ih]],  xlab=names(balance_res$hist)[icov], 
                       tit2="\nApproach: Matching", tit1=names(balance_res$hist[[icov]])[ih]  )
}


catt <- function(text, log_text=".log_text", append=T, lprint=T, save=F, new_name=log_text, dir="", envir = .GlobalEnv){
  if(exists(log_text,envir=envir) & append) assign(log_text,c( get(log_text,envir=envir), list(text)), envir=envir)
  else assign(log_text, text, envir=envir)
  if(lprint) cat(text)
  if(save) save_log_file(log_text=log_text, new_name=new_name, dir=dir, envir=envir)
  invisible()
}

printt <- function(text, log_text=".log_text", append=T, lprint=T, save=F, new_name=log_text, dir="", envir = .GlobalEnv){
  if(exists(log_text,envir=envir) & append) assign(log_text,c( get(log_text,envir=envir), list(text)), envir=envir)
  else assign(log_text, text, envir=envir)
  if(lprint) print(text)
  if(save) save_log_file(log_text=log_text, new_name=new_name, dir=dir, envir=envir)
  invisible()
}

warningg <- function(text, log_text=".log_text", append=T, save=F, new_name=log_text, dir="", envir = .GlobalEnv){
  if(exists(log_text,envir=envir) & append) assign(log_text,c( get(log_text,envir=envir), list(text)), envir=envir)
  else assign(log_text, text, envir=envir)
  warning(text)
  if(save) save_log_file(log_text=log_text, new_name=new_name, dir=dir, envir=envir)
  invisible()
}

stopp <- function(text, log_text=".log_text", append=T, save=F, new_name=log_text, dir="", envir = .GlobalEnv){
  if(exists(log_text,envir=envir) & append) assign(log_text,c( get(log_text,envir=envir), list(text)), envir=envir)
  else assign(log_text, text, envir=envir)
  stop(text)
  if(save) save_log_file(log_text=log_text, new_name=new_name, dir=dir, envir=envir)
  invisible()
}

save_log_file <- function(log_text=".log_text", new_name=log_text, dir="", envir = .GlobalEnv){
  if(new_name!=log_text) assign(new_name, get(log_text, envir=envir), envir=envir)
  save(list=new_name, file=paste0(dir,new_name,".RData"), envir = envir)
  sink(file=paste0(dir,new_name,".txt")) 
  opt <- options(max.print=99999)
  if(exists(new_name,envir=envir)) print(get(new_name,envir=envir))
  options(opt)
  sink()
}

# Function:     table1
# Description:  create a table with counts and percentages for one categorical variable. Similar to 'tabyl'.
#

table1 <- function(x, title="", digits=2, sep=" & ", lprint=c(T,T,T,T) ){
  if(!is.null(dim(x)) & length(dim(x))==2){
    for(icol in 2:ncol(x)) x[,1] <- paste(x[,1],x[,icol],sep=sep)
    x <- x[,1]
  }  
  else x <- as.factor(x)
  res <- cbind( 
    n=(tb<-table( x, useNA="ifany" )), 
    cum_n=cumsum(tb), 
    percent=round(100*tb/sum(tb),digits), 
    cum_percent=round(100*cumsum( tb/sum(tb) ),digits), 
    percent2=c(round(100*(tb2<-table(x))/sum(tb2),digits),rep(NA,length(tb)-length(tb2))),
    cum_percent2 = c(round(100*cumsum( (tb2<-table(x))/sum(tb2) ),digits),rep(NA,length(tb)-length(tb2))) 
  )[,c( lprint, any(is.na(x)) & lprint[3] , any(is.na(x)) & lprint[4] ), drop=F]
  if(title!=""){ 
    attributes(res) <- c(attributes(res), title)
    if(any(lprint)) cat(paste(title,"\n") )
  }
  res
}
#
#
###

