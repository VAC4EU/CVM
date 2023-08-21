
# Program Information  ----------------------------------------------------
#
#  functions for SCRI analysis 

# Functions:    scri_fit
#               refresh_event_variable
#               create_rws
#               split_intervals
#               summary_tab
#               factor_ref
#               combine_vars
#               plot_res
#               add_to_models_list
#               more functions .....
#
# Author:       Svetlana Belitser 
#               nov 2021 - feb 2022
#


scri_fit <- function( formula="",  
                      vax_def,  
                      event_time, event,
                      id,
                      rws,                          # list of risk/control windows definitions
                      time_seq=c(), split_seq_name = "cal_time_cat", time_seq_ref="with events", # time_seq_ref="most events", 
                      combine_vars = c(),           # list of parameters to create new one variable from two other variables
                      start_obs, end_obs,
                      data, 
                      strata_var="", strata_value=NA, time_adj_events="strata",   # "strata"-all events in data[strata_var==T,];         #use_all_events=F,
                      # "all"-all events in 'data'; 
                      # "unadj"-only that in unadjusted model        
                      data_event,
                      data_source = "",
                      nvax,                          # ???if missing ==> maximum of doses
                      lab_orders = NA,
                      ref=1,
                      rw_observed_percentage=0,     # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                      censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                      #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                      event_in_rw=T,                # if event in rw ==> this rw should not be deleted even if not completely observed
                      delete_coef_no_events     = T,
                      delete_rows_start_with_no = T,
                      delete_no_ref_cols        = T,
                      delete_median_cols        = T,
                      lprint                    = T,       test = F,
                      lplot_hist=T, add_histplot=F, sdr_tabs="", width=14,
                      save_data = F,
                      ...
){  
  
														
  if(missing(rws       )) stop("'rws' is missing.")
  if(missing(event_time)) stop("'event_time' is missing.")
  if(missing(event     )) stop("'event' is missing.")
  if(missing(id        )) stop("'id' is missing.")
  if(missing(data      )) stop("'data' is missing.")
  if(missing(start_obs )) stop("'start_obs' is missing.")
  if(missing(end_obs   )) stop("'end_obs' is missing.")
  if(missing(nvax      )) stop("'nvax' is missing.")
  
  if(nrow(data)==0)     
    return(  list( tabs     = NULL, 
                   tab_full = NULL,
                   model    = NULL,
                   call     = list( match.call())  ))
  
  # stratum:
  if(strata_var!=""){ 
    
    if(mode(data[,strata_var])=="logical") data$strata_cond <- data[,strata_var] & !is.na(data[,strata_var])
    else {
      if(!is.na(strata_value)) data$strata_cond <- data[,strata_var]==strata_value & !is.na(data[,strata_var])
      else stop(paste0("'strata_var' should be logical or a variable with value 'strata_value'. (Now: 'strata_var'= '",strata_var,"'; 'strata_value'=",ifelse(is.na(strata_value),"NA",paste0("'",strata_value,"'")),")"))
    }
    if( sum(data$strata_cond)==0) return(  list( tab      = NULL, 
                                                 tab_full = NULL,
                                                 model    = NULL,
                                                 call     = list( match.call())  ))
													   
  }
  
  if(formula!=""){
    
    formula  <- formula(formula)
    tb       <- attributes(terms(formula))$factor
    tab_vars <- dimnames(tb)[[2]]
    
    if(missing(event)      & any(rowSums(tb)==0))  event      <- dimnames(tb)[[1]][rowSums(tb)==0]
    if(missing(event_time) & any(rowSums(tb)==0))  event_time <- dimnames(tb)[[1]][rowSums(tb)==0]
    
    no_formula <- F
  }
  else no_formula <- T
  
  
  ###########################################################################
  ################# create_rws (v3)  ###############
  data_rws  <- create_rws(
    obj = vax_def,
    rws = rws,
    data = data,
    strata_cond = ifelse(strata_var=="", F, T), strata_value = strata_value,
    start_obs = start_obs, end_obs = end_obs,
    event_time = event_time, event = event,
    id = id,
    lab_orders = lab_orders,
																								
    #ref="pre-",       #  ref=5 OR ref= "pre-exposure [-90;-30]"  or a part of a ref.category name
    rw_observed_percentage = rw_observed_percentage,   
    censored_vars = censored_vars,           
    event_in_rw = event_in_rw               
  ) 
  if(!is.data.frame(data_rws)) {
    sep_vars <- data_rws$sep_vars
    data_rws <- data_rws$data_rws
  }
  else sep_vars <- c()
  
  
  
  data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, event, event_time)
  
  data_rws$interval <- data_rws$rw_end - data_rws$rw_start + 1
  data_rws <- data_rws[data_rws$interval>0,]
  
  # delete id's without events in the windows:
  id_no_events <- names((tb<-tapply(data_rws[,event],data_rws[,id], sum, na.rm=T))[tb==0 & !is.na(tb)])
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  data_rws <- data_rws[ !(as.character(data_rws[,id]) %in% id_no_events),  ]
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  
  
  if(nrow(data_rws)==0) return(  list( tab      = NULL, 
                                       tab_full = NULL,
                                       model    = NULL,
                                       call     = list( match.call())  ))
  
  risk_time_tab <- list( event_tab=tapply(data_rws[,event],data_rws$lab,sum,na.rm=T), atrisk_tab=tapply(data_rws$interval,data_rws$lab,sum,na.rm=T) )
  
  split_seq_name0 <- split_seq_name
  sep_time_vars_cond <- time_seq_ref %in% c("with_events","with events","without_events","without events")

  #######
  #  create time intervals:
  #
  if(length(time_seq)>0 & nrow(data_rws)>0){ 

    data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, event, event_time)
    data_rws <- split_intervals( data =data_rws, 
                                 start_interval = "rw_start", end_interval = "rw_end", 
                                 splits_names = split_seq_name, # "cal_time_cat",
                                 splits       = time_seq,
                                 ref          = ifelse(sep_time_vars_cond,1,time_seq_ref), #"most events", 
                                 event        = event,        #   myopericarditis 
                                 event_time   = event_time 
    )
    
    # "standard" ==> i.e., if at risk only 'ref' or only non-"ref" ==> add this interval to a control category!!!
    #            ==>  ref_cat for time_intervals is union of all intervals with days at risk only 'ref' or only non-"ref".   
    rws <- eval(parse(text=rws))
    if(is.null(rws$ref)) {
      if( !is.null( lapply(rws,function(x) is.null(x$ref))) ) 
        rws$ref <- unlist(lapply(rws,function(x) x$ref))[1]
      else rws$ref <- data_rws$lab[order(data_rws$rw_start)][1]
    }
    
   
    data_rws[,paste0(split_seq_name,"0")] <- data_rws[,split_seq_name]
    
    risk_time_tab0 <-  list( event_tab0=tapply(data_rws[,event],list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T), 
                             at_risk_tab0=tapply(data_rws$interval,list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T) )
    
    if(sep_time_vars_cond){ 
      
      if( time_seq_ref %in% c("with_events","with events"))        ref_and_nonref_time_interval <- with( data_rws[data_rws[,event]>0,], tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),lab), get(split_seq_name), function(x) length(table(x))>1) ) 
      if( time_seq_ref %in% c("without_events","without events"))  ref_and_nonref_time_interval <- with( data_rws,                      tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),lab), get(split_seq_name), function(x) length(table(x))>1) ) 
      ref_and_nonref_time_interval[is.na( ref_and_nonref_time_interval)] <- F
      
      if(any(!ref_and_nonref_time_interval)){
        ref_time_names <- names(ref_and_nonref_time_interval)[!ref_and_nonref_time_interval]
        ref_time_name <- range(as.numeric(unlist(strsplit( substring(ref_time_names,2,nchar(ref_time_names)-1) ,";"))))
        ref_time_name <- paste0("[",ref_time_name[1],";...",ref_time_name[2],"]")
        data_rws[ data_rws[,split_seq_name] %in% ref_time_names, split_seq_name ] <- ref_time_names[1]
        levels(data_rws[,split_seq_name])[levels(data_rws[,split_seq_name])==ref_time_names[1]] <- ref_time_name
        
        data_rws[,split_seq_name] <- factor(  as.character(data_rws[,split_seq_name]), levels=unique(as.character(data_rws[,split_seq_name]))[match( levels(data_rws[,split_seq_name]) , unique(as.character(data_rws[,split_seq_name])))] )
        if(nlevels(data_rws[,split_seq_name])>1)
          contrasts(data_rws[,split_seq_name])[,] <- contr.treatment(nlevels(data_rws[,split_seq_name]), base=(1:nlevels(data_rws[,split_seq_name]))[levels(data_rws[,split_seq_name])==ref_time_name]) 
        
        risk_time_tab <-  c( list( event_tab=tapply(data_rws[,event],list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T), 
                                   at_risk_tab=tapply(data_rws$interval,list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T) ),
                             risk_time_tab0 )
        
        if(length(sep_vars)==0){  # ==> one "cal_time_cat" variable!
          
          for(iv in rev(names(rws)[names(rws)!="ref"][-1]))
            if(paste0(iv,"_lab") %in% names(data_rws)) {       
              risk_time_tab0 <- c( list( events  = with(data_rws[data_rws[,event]==1 & !(substring(tolower(data_rws[,paste0(iv,"_lab")]),1,3) %in% c("no ","no_")),], 
                                                        tapply(get(event),list(get(paste0(iv,"_lab")), cal_time_cat0),sum,na.rm=T) ), 
                                         at_risk = with(data_rws[!(substring(tolower(data_rws[,paste0(iv,"_lab")]),1,3) %in% c("no ","no_")),],
                                                        tapply(interval,list(get(paste0(iv,"_lab")), cal_time_cat0),sum,na.rm=T) )), 
                                   risk_time_tab0)
              names(risk_time_tab0)[1:2] <- paste0(c("events","at_risk"),"_tab0_",iv)
              
              
              if( time_seq_ref %in% c("with_events","with events"))        ref_and_nonref_time_interval <- with( data_rws[data_rws[,event]>0,], tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),get(paste0(iv,"_lab"))), cal_time_cat0, function(x) length(table(x))>1) ) 
              if( time_seq_ref %in% c("without_events","without events"))  ref_and_nonref_time_interval <- with( data_rws,                      tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),get(paste0(iv,"_lab"))), cal_time_cat0, function(x) length(table(x))>1) ) 
              ref_and_nonref_time_interval[is.na( ref_and_nonref_time_interval)] <- F
              if(any(!ref_and_nonref_time_interval)){
                ref_time_names <- names(ref_and_nonref_time_interval)[!ref_and_nonref_time_interval]
                ref_time_name <- range(as.numeric(unlist(strsplit( substring(ref_time_names,2,nchar(ref_time_names)-1) ,";"))))
                ref_time_name <- paste0("[",ref_time_name[1],";...",ref_time_name[2],"]")

                data_rws[ data_rws[,split_seq_name] %in% ref_time_names, split_seq_name] <- ref_time_names[1]
                levels(data_rws[,split_seq_name])[levels(data_rws[,split_seq_name])==ref_time_names[1]] <- ref_time_name
                
                data_rws[,split_seq_name] <- factor(  as.character(data_rws[,split_seq_name]), levels=unique(as.character(data_rws[,split_seq_name]))[match( levels(data_rws[,split_seq_name]) , unique(as.character(data_rws[,split_seq_name])))] )
                if(nlevels(data_rws[,split_seq_name])>1)
                  contrasts(data_rws[,split_seq_name])[,] <- contr.treatment(nlevels(data_rws[,split_seq_name]), base=(1:nlevels(data_rws[,split_seq_name]))[levels(data_rws[,split_seq_name])==ref_time_name]) 
                
                risk_time_tab <-  c( list( event_tab=tapply(data_rws[,event],list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T), 
                                           at_risk_tab=tapply(data_rws$interval,list(data_rws$lab, data_rws[,split_seq_name]),sum,na.rm=T) ),
                                     risk_time_tab0 )
              }
            }
        }
        else {  # length(sep_vars)>0  & sep_time_vars_cond  ==> example, per brand    ==>  "cal_time_cat_d1", "cal_time_cat_d2", "cal_time_cat_d3", ...
       
          for(iv in rev(names(rws)[names(rws)!="ref"][-1]))   # "d1", "d2", ...
            if(paste0(iv,"_lab") %in% names(data_rws)) {       
              risk_time_tab0 <- c( list( events  = with(data_rws[data_rws[,event]==1 & !(substring(tolower(data_rws[,paste0(iv,"_lab")]),1,3) %in% c("no ","no_")),], 
                                                        tapply(get(event),list(get(paste0(iv,"_lab")), cal_time_cat0),sum,na.rm=T) ), 
                                         at_risk = with(data_rws[!(substring(tolower(data_rws[,paste0(iv,"_lab")]),1,3) %in% c("no ","no_")),],
                                                        tapply(interval,list(get(paste0(iv,"_lab")), cal_time_cat0),sum,na.rm=T) )), 
                                   risk_time_tab0)
              names(risk_time_tab0)[1:2] <- paste0(c("events","at_risk"),"_tab0_",iv)
              
              
              if( time_seq_ref %in% c("with_events","with events"))        ref_and_nonref_time_interval <- with( data_rws[data_rws[,event]>0,], tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),get(paste0(iv,"_lab"))), cal_time_cat0, function(x) length(table(x))>1) ) 
              if( time_seq_ref %in% c("without_events","without events"))  ref_and_nonref_time_interval <- with( data_rws,                      tapply(grepl(paste0(rws[[1]]$lab,collapse="|"),get(paste0(iv,"_lab"))), cal_time_cat0, function(x) length(table(x))>1) ) 
              ref_and_nonref_time_interval[is.na( ref_and_nonref_time_interval)] <- F
              
              if(any(!ref_and_nonref_time_interval)){
                ref_time_names <- names(ref_and_nonref_time_interval)[!ref_and_nonref_time_interval]
                ref_time_name <- range(as.numeric(unlist(strsplit( substring(ref_time_names,2,nchar(ref_time_names)-1) ,";"))))
                ref_time_name <- paste0("[",ref_time_name[1],";...",ref_time_name[2],"]")
                data_rws[,paste0(split_seq_name,"_",iv)] <- data_rws[,paste0(split_seq_name,"0")]
                data_rws[ data_rws[,paste0(split_seq_name,"_",iv)] %in% ref_time_names, paste0(split_seq_name,"_",iv)] <- ref_time_names[1]
                levels(data_rws[,paste0(split_seq_name,"_",iv)])[levels(data_rws[,paste0(split_seq_name,"_",iv)])==ref_time_names[1]] <- ref_time_name
                
                data_rws[,paste0(split_seq_name,"_",iv)] <- factor(  as.character(data_rws[,paste0(split_seq_name,"_",iv)]), levels=unique(as.character(data_rws[,paste0(split_seq_name,"_",iv)]))[match( levels(data_rws[,paste0(split_seq_name,"_",iv)]) , unique(as.character(data_rws[,paste0(split_seq_name,"_",iv)])))] )
                if(nlevels(data_rws[,paste0(split_seq_name,"_",iv)])>1)
                  contrasts(data_rws[,paste0(split_seq_name,"_",iv)])[,] <- contr.treatment(nlevels(data_rws[,paste0(split_seq_name,"_",iv)]), base=(1:nlevels(data_rws[,paste0(split_seq_name,"_",iv)]))[levels(data_rws[,paste0(split_seq_name,"_",iv)])==ref_time_name]) 
                
                risk_time_tab <-  c( list( event_tab=tapply(data_rws[,event],list(data_rws$lab, data_rws[,paste0(split_seq_name,"_",iv)]),sum,na.rm=T), 
                                           at_risk_tab=tapply(data_rws$interval,list(data_rws$lab, data_rws[,paste0(split_seq_name,"_",iv)]),sum,na.rm=T) ),
                                     risk_time_tab0 )
              }
            }
          
        }
      }
    }
    
  }
  else {   # length(time_seq)==0
    
    if(lplot_hist){
 
      gc()
      if(file.exists(paste0(sdr_tabs,"histplots_tmp.pdf"))){
        if(add_histplot) file.copy(from=paste0(sdr_tabs,"histplots_tmp.pdf"), to=paste0(sdr_tabs,"histplots_tmp1.pdf"))
        suppressWarnings( file.remove(paste0(sdr_tabs,"histplots_tmp.pdf")) )
      }
      
      pdf(file=paste0(sdr_tabs,"histplots_tmp.pdf"), width=width,  ... )
      {
        par(mfrow=c(2,2))
        data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, event, event_time)
        ids_with_events <- names((tb<-tapply(data_rws[,event],data_rws[,id], sum, na.rm=T))[tb>0])
        
        suppressWarnings(try(hist_events_model(data,data_rws[data_rws[,id] %in% ids_with_events,], event=event, tit=paste(ifelse(is.na(strata_value),"",strata_value), substring(rws,1,regexpr(" = ",rws)-1), collapse="; ") )))
        
        
        data$death_date <- data$date_of_death
        data_rws$death_date <- data_rws$date_of_death
        data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, "death","death_days") 
        ids_with_deaths <- names((tb<-tapply(data_rws[,"death"],data_rws[,id], sum, na.rm=T))[tb>0])
        
        suppressWarnings(try(hist_events_model(data,data_rws[data_rws[,id] %in% ids_with_deaths,], event="death",tit=paste(ifelse(is.na(strata_value),"",strata_value), substring(rws,1,regexpr(" = ",rws)-1), collapse="; ") )))
      }
      dev.off()
      
      #if(file.exists(paste0(sdr_tabs,"histplots_tmp1.pdf")) & add_histplot ){ 
      #  qpdf::pdf_combine( c(paste0(sdr_tabs,"histplots_tmp1.pdf"),paste0(sdr_tabs,"histplots_tmp.pdf")) , 
      #                     paste0(sdr_tabs,"histplots_tmp2.pdf") )  
      #  file.rename(from=paste0(sdr_tabs,"histplots_tmp2.pdf"), to=paste0(sdr_tabs,"histplots_tmp.pdf"))
      #  file.remove(from=paste0(sdr_tabs,"histplots_tmp1.pdf"))
      #}
    }
  }
  
  
  if(length(combine_vars)>0 & nrow(data_rws)>0){
    data_rws$lab  <- combine_vars_func( data_rws[, c(combine_vars, "lab"), drop=F ], lab_orders = lab_orders, ref=ref, event = data_rws[,event] )
    data_rws$lab  <- factor_ref(  data_rws$lab, lab_orders=lab_orders, ref=ref, event_var=data_rws[,event] )  
  }
  
  
  #####
  #  create 'event'and 'interval' variables
  #
  data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, event, event_time)
  
  data_rws$interval <- data_rws$rw_end - data_rws$rw_start + 1
  data_rws <- data_rws[data_rws$interval>0,]
  
  # delete id's without events in the windows:
  id_no_events <- names((tb<-tapply(data_rws[,event],data_rws[,id], sum, na.rm=T))[tb==0 & !is.na(tb)])
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  data_rws <- data_rws[ !(as.character(data_rws[,id]) %in% id_no_events),  ]
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  
  
  if(nrow(data_rws)==0) return(  list( tab      = NULL, 
                                       tab_full = NULL,
                                       model    = NULL,
                                       call     = list( match.call())  ))
  #return(  list( res_tab = NULL, 
  #               model   = NULL,
  #               call    = list( match.call())  ))
  
  
										  
																						
																																											  
																 
																			 
  
  
																							 
				  
										   
																  
																							
   
  
  if(any(names(data_rws) %in% split_seq_name))
    for(iname in split_seq_name)
      data_rws[,iname]  <- factor_ref(  as.character(data_rws[,iname]), 
                                        lab=levels(data_rws[,iname])[levels(data_rws[,iname]) %in% unique(as.character(data_rws[,iname])) ],
                                        ref=time_seq_ref, 
                                        event_var=data_rws[,event] )  
  
  if(length(sep_vars)>0){  
      
    mods <- list()
    if(length(time_seq)>0 & sep_time_vars_cond) split_seq_name <- paste0(split_seq_name,"_", names(rws)[names(rws)!="ref"][-1])
								   
    
    res_tab <- NULL
    for(isep in 1:length(sep_vars)){
     
      # combine variables in formula with risk windows variable 'lab' (the begining of the code):
      if(!no_formula){
        tb <- attributes(terms(formula))$factor
        if(length(time_seq)>0 & sep_time_vars_cond & split_seq_name0 %in% dimnames(tb)[[2]]){
          dimnames(tb)[[1]][dimnames(tb)[[1]]==split_seq_name0] <- split_seq_name[isep]
          dimnames(tb)[[2]][dimnames(tb)[[2]]==split_seq_name0] <- split_seq_name[isep]
          tb <- tb[, c( dimnames(tb)[[2]][dimnames(tb)[[2]]!=split_seq_name[isep]], split_seq_name[isep] ) ]
          
        }
      }
      
      var_names0 <- sep_vars[[isep]]
      if(length(time_seq)>0){
        if(length(split_seq_name)>=isep) var_names0 <- c( var_names0, split_seq_name[isep] )   #& split_seq_name %in% dimnames(tb)[[2]])
        else var_names0 <- c( var_names0, split_seq_name )   #& split_seq_name %in% dimnames(tb)[[2]])
      }
	  
      # delete id's without events in the windows for this model:
      ids_with_events <- names((tb<-tapply(data_rws[,event],data_rws[,id], sum, na.rm=T))[tb>0])
      data_regr <- data_rws[data_rws[,id] %in% ids_with_events, names(data_rws) %in% c(event,var_names0,id,"interval")]
      
      # delete covariates if there is only one level with event[s]
      var_names <- var_names0 <- var_names0[ sapply(var_names0, function(x)sum(table( data_regr[ data_regr[,event]>0 & !(tolower(substring(data_regr[,x],1,3))%in%c("no_","no ")), x])>0,na.rm=T)>0) ]
      
      if(length(var_names0)==0) next
      if(length(var_names0)==1 & ( length(time_seq)==0 | !sep_time_vars_cond ) ) if(var_names0==split_seq_name) next
      if(length(var_names0)==1 &  length(time_seq)>0   &  sep_time_vars_cond )   if(var_names0==split_seq_name[isep]) next
      
      var_names  <- var_names[ sapply(var_names, function(x)sum(table( data_regr[ data_regr[,event]>0 & !(tolower(substring(data_regr[,x],1,3))%in%c("no_","no ")), x])>0,na.rm=T)>1) ]
      if(length(var_names)==1 & (length(time_seq)==0 | !sep_time_vars_cond) ) if(var_names[1]==split_seq_name) next 
      if(length(var_names)==1 & length(time_seq) >0  &  sep_time_vars_cond  ) if(var_names[1]==split_seq_name[isep]) next 
      
      #  data_rws$cal_time_mid <- I(unlist(lapply(strsplit(substring(as.character(data_rws[,split_seq_name]),2,nchar(as.character(data_rws[,split_seq_name]))-1),";",fixed=T),function(x)mean(as.numeric(x)))))
      Poisson_formula <- as.formula(paste(event,"~", 
                                          paste( var_names,  collapse=" + "), "+", 
                                          "strata(",id,")", "+", "offset(log(interval))")) 
      Poisson_formula_coxph <- as.formula(paste("Surv( rep(1,",nrow(data_regr), "),",event,") ~", 
                                                paste( var_names,  collapse=" + "), "+", 
                                                "strata(",id,")", "+", "offset(log(interval))")) 
      
	  
	  
	  
      if(!no_formula & sum(data_regr[,event],na.rm=T)>1 & length(var_names)>0){
        suppressWarnings(
          mod <- try( clogit(formula = Poisson_formula, data = data_regr, control=coxph.control(iter.max=1000) ), silent = T )
        )
        if(class(mod)[[1]]== "try-error")
          suppressWarnings(
            mod <- try( coxph(formula = Poisson_formula_coxph, data = data_regr, control=coxph.control(iter.max=1000) ), silent = T )
          )
		
		
		
		
      }
      else mod <- NA
      
      
      if( class(mod)[[1]]== "try-error" | no_formula | sum(data_regr[,event],na.rm=T)<=1 | length(var_names)==0 ){
        
        attributes(mod)$var_names_data      <- var_names0
        attributes(mod)$var_names_mod       <- var_names
        attributes(mod)$formula             <- deparse(Poisson_formula)
        attributes(mod)$formula_only_paste  <- paste(event,"~", paste( var_names,  collapse=" + "), "+", "strata(",id,")", "+", "offset(log(interval))")
        attributes(mod)$summary_with_events <- try(summary(data_regr[data_regr[,event]==1,var_names0]), silent = T )
        attributes(mod)$summary             <- try(summary(data_regr[,var_names0]), silent = T )
        
        
        if(F){
          if(nrow(table1(data_rws[  data_rws[,event]==1 ,var_names]))<=8)
            print( table1(data_rws[  data_rws[,event]==1 ,var_names]) ) 
          print(cbind.data.frame(var_names, nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names])) ))
          if(any(names(data_rws)==split_seq_name))   #"cal_time_cat"))
            print(cbind.data.frame(var_names[var_names!=split_seq_name], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!=split_seq_name]])) ))
          
          if(length(time_seq)>0 )
            warning("Error in Poisson regression. \ntime_seq=",median(diff(time_seq),na.rm=T),"; Formula: ", deparse(Poisson_formula))
          else warning("Error in Poisson regression. \nno time_seq; Formula: ", deparse(Poisson_formula))
        }        
        
        if(is.null(res_tab))
          res_tab <- summary_tab( var_names_data=var_names0,#sep_vars[[isep]], 
                                  var_names_mod = var_names,
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  lab_orders=lab_orders,
                                  model_number = isep )
        else  
          res_tab <- summary_tab( var_names_data=var_names0,#sep_vars[[isep]], 
                                  var_names_mod = var_names,
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  lab_orders=lab_orders,
                                  model_number = isep,
                                  res_tab = res_tab )
        
        mods <- c(mods, list(list(mod)) )
        names(mods)[length(mods)] <-  names(sep_vars)[isep]
        
        #mods[[isep]] <- list( mod ) 
        
      }
      else{
        
        if(is.null(res_tab))
          res_tab <- summary_tab( var_names_data=var_names0,#sep_vars[[isep]], 
                                  var_names_mod = var_names,
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  #coef_cond = coef_cond, 
                                  model_number = isep )
        else  
          res_tab <- summary_tab( var_names_data=var_names0,#sep_vars[[isep]], 
                                  var_names_mod = var_names,
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  #coef_cond=coef_cond, 
                                  model_number = isep,
                                  res_tab = res_tab )
        
        mods <- c(mods, list(list(mod)) )
        names(mods)[length(mods)] <-  names(sep_vars)[isep]
        
        # mods[[isep]] <- list(summary(mod)) 
      }
    }         
	
    mod <- mods
    
    if(lprint)  print(format( res_tab, digits=3, nsmall=2, justify="left" ))
    
  } else {  # ==> length(sep_vars)==0
    
    # combine variables in formula with risk windows variable 'lab' (the begining of the code):
    if(!no_formula){
      tb <- attributes(terms(formula))$factor
      if(length(time_seq)>0 & split_seq_name %in% dimnames(tb)[[2]])
        tb <- tb[, c( dimnames(tb)[[2]][dimnames(tb)[[2]]!=split_seq_name], split_seq_name ) ]
    }
    
    if(nrow(data_rws)==0) res_tab <- mod <- NULL
    else {
      
      if(any(cond<-(colSums(tb)>1)) & nrow(data_rws)>0)
        for(i in 1:sum(cond)){
          vars_tmp <- gsub(" ","",dimnames(tb)[[2]][cond][i])
          dimnames(tb)[[2]][cond][i] <- paste(strsplit( vars_tmp, ":" )[[1]],collapse="_")
          data_rws[,paste(strsplit( gsub(" ","",vars_tmp,fixed = TRUE) ,":")[[1]],collapse="_")]  <- combine_vars_func( 
            data_rws[, strsplit(gsub(" ","",vars_tmp,fixed = TRUE)  ,":")[[1]], drop=F ], lab_orders = lab_orders, event = data_rws[,event] )   #ref=ref, 
        }  
      
      var_names0 <- dimnames(tb)[[2]]
      
      if(length(time_seq)>0 & split_seq_name %in% dimnames(tb)[[2]])
        var_names0 <- unique(c( var_names0, split_seq_name ))
    
      # delete id's without events in the windows for this model:
      ids_with_events <- names((tb_events<-tapply(data_rws[,event],data_rws[,id], sum, na.rm=T))[tb_events>0])
      data_regr <- data_rws[data_rws[,id] %in% ids_with_events, names(data_rws) %in% c(event,var_names0,id,"interval")]
      
      # delete covariates if there is only one level with event[s]
      var_names <- var_names0 <- var_names0[ sapply(var_names0, function(x)sum(table( data_regr[ data_regr[,event]>0 & !(tolower(substring(data_regr[,x],1,3))%in%c("no_","no ")), x])>0,na.rm=T)>0) ]
									 
																	
																																													   
																								   
	  
			   
							   
															 
																	
						 
							  
																																										
	   
													 
																				   
																						   
	  
																								  
																						 
																								 
	  
	  
      
      run_model <- T
      if(length(var_names0)==0) run_model <- F
      if(length(var_names0)==1) if(var_names0==split_seq_name)  run_model <- F
      
      if(run_model){
        var_names  <- var_names[ sapply(var_names, function(x)sum(table( data_regr[ data_regr[,event]>0 & !(tolower(substring(data_regr[,x],1,3))%in%c("no_","no ")), x])>0,na.rm=T)>1) ]
        if(length(time_seq)>0 & length(var_names)==1 & var_names[1]==split_seq_name) var_names <- c()
		 
										 
						   
																																	 
		   
		
        
        if(test){
          print(paste(event,"~", 
                      paste( var_names,  collapse=" + "), "+", 
                      "strata(",id,")", "+", "offset(log(interval))"))
          print("ref_cat:")
          for(ivar in var_names)
            print( c( ivar, dimnames(contrasts(data_rws[,ivar]))[[1]][ !( dimnames(contrasts(data_rws[,ivar]))[[1]] %in% dimnames(contrasts(data_rws[,ivar]))[[2]] )  ] ))
        }
        Poisson_formula <- as.formula(paste(event,"~", 
                                            paste( var_names,  collapse=" + "), "+", 
                                            "strata(",id,")", "+", "offset(log(interval))")) 
        
        Poisson_formula_coxph <- as.formula(paste("Surv( rep(1,",nrow(data_regr), "),",event,") ~", 
                                                  paste( var_names,  collapse=" + "), "+", 
                                                  "strata(",id,")", "+", "offset(log(interval))")) 
        
		
					
	  
																												  
        
        if(!no_formula & sum(data_regr[,event],na.rm=T)>1 & length(var_names)>0){
          suppressWarnings(
            mod <- try( clogit(formula = Poisson_formula, data = data_regr, control=coxph.control(iter.max=1000) ), silent = T )
          )
          if(class(mod)[[1]]== "try-error")
            suppressWarnings(
              mod <- try( coxph(formula = Poisson_formula_coxph, data = data_regr, control=coxph.control(iter.max=1000) ), silent = T )
            )
          
        } 
        else mod <- NA
        
        if( class(mod)[[1]]== "try-error" | no_formula | sum(data_regr[,event],na.rm=T)<=1 | length(var_names)==0 ){
																		 
																		
																											
																	  
																																									
																																								   
          
          attributes(mod)$var_names_data      <- var_names0
          attributes(mod)$var_names_mod       <- var_names
          attributes(mod)$formula             <- deparse(Poisson_formula)
          attributes(mod)$formula_only_paste  <- paste(event,"~", paste( var_names,  collapse=" + "), "+", "strata(",id,")", "+", "offset(log(interval))")
          attributes(mod)$summary_with_events <- try(summary(data_regr[data_regr[,event]==1,var_names0]), silent = T )
          attributes(mod)$summary             <- try(summary(data_regr[,var_names0]), silent = T )
          
          if(F){
            if(nrow(table1(data_rws[  data_rws[,event]==1 ,var_names]))<=8)
              print( table1(data_rws[  data_rws[,event]==1 ,var_names]) ) 
            print(cbind.data.frame(var_names, nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names])) ))
            if(any(names(data_rws)==split_seq_name))   #"cal_time_cat"))
              print(cbind.data.frame(var_names[var_names!=split_seq_name], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!=split_seq_name]])) ))
            #print(cbind.data.frame(var_names[var_names!="cal_time_cat"], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!="cal_time_cat"]])) ))
            
            if(length(time_seq)>0)
              warning("Error in Poisson regression. \ntime_seq=",median(diff(time_seq),na.rm=T),"; Formula: ", deparse(Poisson_formula))
            else warning("Error in Poisson regression. \nno time_seq; Formula: ", deparse(Poisson_formula))
          }   
          
          res_tab <- summary_tab( var_names_data=var_names0, var_names_mod=var_names, event=event, data=data_rws, id_name=id, lab_orders=lab_orders )
        }
        else{
          res_tab <- summary_tab( var_names_data=var_names0, var_names_mod=var_names, event=event, data=data_rws, id_name=id,  mod=mod, lab_orders=lab_orders  )
          mod <- summary(mod)
        }
        
        if(lprint)  print(format( res_tab, digits=2, justify="left" ))
	   
		   
																																							  
						   
      }
	  
																	
    }
  }
  
  
  if(is.null(res_tab)) return( list( tab = NULL,
                                     tab_full = NULL,
                                     model = NULL,
                                     call = list( match.call()) ))
 
  
 
  
  if(data_source!="") res_tab <- cbind.data.frame( data = data_source, res_tab )
  rownames(res_tab) <- NULL
  
  res_tab0 <- res_tab
  if( delete_coef_no_events & "RR" %in% names(res_tab) ){ 
    res_tab[ is.na(res_tab$events_rw) | res_tab$events_rw==0 , c("events_rw","RR_data","RR") ] <- NA 
    res_tab[ is.infinite(res_tab$uci),                                                 "RR"  ] <- NA 
    res_tab[ is.na(res_tab$events_ref) & res_tab$uci>10^9 & !is.na(res_tab$uci),       "RR"  ] <- NA 
  }
  if(delete_rows_start_with_no)
    res_tab <- res_tab[ !(substring(tolower(res_tab$all_cat),1,3) %in% c("no ","no_")), ]
  #if(delete_no_brand_rows & any(ls()=="brands"))
  #res_tab <- res_tab[ !grepl(paste0("no ",brands,collapse="|"),res_tab$all_cat,perl=T), ]
  if(delete_no_ref_cols & any(  !grepl("noref",names(res_tab)) ))
    res_tab <- res_tab[ ,!grepl("noref",names(res_tab)) ]
  if(delete_median_cols & any(  !grepl("median",names(res_tab)) ))
    res_tab <- res_tab[ ,!grepl("median",names(res_tab)) ]
  
  if(strata_var!="" | !is.na(strata_value)) 
    res_tab <- cbind.data.frame( res_tab[,1:match("i",names(res_tab))], 
                                 strata = ifelse(!is.na(strata_value),strata_value,strata_var), 
                                 res_tab[,match("i",names(res_tab)):ncol(res_tab)]               )
  
  
  ret <- list( tab      = res_tab,
               tab_full = res_tab0,
               model    = mod,
               risk_time_tab = risk_time_tab,
               call     = list( match.call()) 
  )
  if(save_data) ret <- c( ret, data_rws=list(data_rws) )
  
  ret
  
}  # end of function 'scri_fit'



######################################
#######################################
######################################  
#  
#
#
refresh_event_variable <- function(  start_interval, end_interval, 
                                     data,
                                     event, event_time#, 
                                     #data_event
){
  data[ ,event] <- as.numeric( data[,start_interval]<=data[,event_time] & data[, event_time] <= data[,end_interval] & !is.na(data[,event_time]) )
  data[is.na(data[,event]),event] <- 0
  data
}  

refresh_event_variable_horiz <- function(  start_interval, end_interval, 
                                           data,
                                           event, event_time
){
  data[ ,event] <- as.numeric( data[,start_interval]<=data[,event_time] & data[, event_time] <= data[,end_interval] & data[,event]==1 )
  data[is.na(data[,event]),event] <- 0
  data
}  

######################################################################  
########################################################################  
########################################################################
#
#             create rws (version 3)  31.01.2022
#
#
create_rws <- function( rws,                          # list of risk/control windows definitions
                        data,                         # dataset
                        obj,
                        strata_cond = F, strata_value=NA,
                        start_obs, end_obs,           # start and end of observation variables
                        var_type,                     # if missing ==> vax1, vax2, vax3, ...
                        id,                           # id variable to get new start and end of observation  ('rw_min' and 'rw_max')
                        event_time, event,
                        lab_orders,
                        ref=1,                        # reference category for new variable: number OR [the beginning of ] category name OR "most events"
                        event_name="event",           # used if ref=="most events" to define reference category with most events
                        rw_observed_percentage=0,     # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                        censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                        #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                        event_in_rw=T                 # if event in rw ==> this rw shouldnot be deleted even if not completely observed
){ 
  
  # create and calculate variable from formula from list 'rws' and dataset 'data': 'lab' (or another name), 'rw_start' and 'rw_end':
  calc_wind <- function(x, data, obj){ 
    
    cond_rws <- rep(T, nrow(data))
    if(!is.null(x$cond)) cond_rws <- eval(x$cond,data)
    else {
      if(!is.null(x$cond_var)){
        if(!is.null(x$cond_values)) cond_rws <- data[,x$cond_var] %in% x$cond_values               
        if(!is.null(x$cond_value))  cond_rws <- data[,x$cond_var] %in% x$cond_value                
        if( is.null(x$cond_values) & is.null(x$cond_value) )  cond_rws <- data[, x$cond_var ]      
      }
    }
    
    data <- data[cond_rws & !is.na(data[, x$t0]),]
    
    if(nrow(data)==0) return(NULL)
    
    data$rws_cat <- x$name
    
    
    t0   <- data[ , x$t0] 
    #if( !x$no_last_interval ){
    if(           x$overlap_priority      == "next"     ) tend <- data[ , obj$data_parameters$next_vax_time ] - 1 
    #      if(          x$overlap_priority      == "overlap"  ) ??? tend <- eval(  data[ , "next"     ] )
    #???      if(       x$overlap_priority      == "previous" ) tend <- eval(  data[ , obj$data_parameters$prev_vax_time ] )
    #      if( substring(x$overlap_priority,1,2) == "no"       ) ?? nothing do   tend <- eval(  data[ , "no"       ] )
    #}      
    
    
    if(!any(names(x)=="cutpoints")) x$cutpoints <- 0
    
    #####
    # labels create
    if(!any(names(x)=="lab")) x$lab <- ""
    
    if( length(x$lab)<length(x$cutpoints) ) x$lab <- c( x$lab, rep(x$lab[length(x$lab)], length(x$cutpoints)-length(x$lab)) )
    
    
    if(x$lab_add_interval){
      x$lab[-length(x$lab)] <- paste0( x$lab[-length(x$lab)], "[", x$cutpoints[-length(x$cutpoints)],";",x$cutpoints[-1]-1,"]")
      if(!x$no_last_interval) x$lab[ length(x$lab)] <- paste0( x$lab[ length(x$lab)], ">", x$cutpoints[ length(x$lab)]-1)
    }
    
    if( any(names(x)=="no_last_interval")) x$no_last_interval <- eval(x$no_last_interval)
    else x$no_last_interval <- F
    
    
    
    
    
    for(i in 1:length(x$cutpoints)){
      
      if(any(ls()=="data_tmp")) rm(data_tmp)
      
      if(i==length(x$cutpoints) & x$no_last_interval) next
      
      if(any(ls()=="tend")) {
        
        if(i <length(x$cutpoints)) rw_enddd <- pmin(tend,t0 + x$cutpoints[i+1]-1,na.rm=T)
        if(i==length(x$cutpoints)) rw_enddd <- tend
        
        cond <- (t0+x$cutpoints[i]) <= tend  # start must be < the end of these intervals 'tend'
        cond[is.na(cond)] <- rep(T,sum(is.na(cond)))  # id's without next dose; later, outside of the small function, check and adding the end of observation if necessary
        
        if(i<length(x$cutpoints) & rw_observed_percentage>0){  #  ?? replace??? from the small function
          cond <- cond &   ( 100 -  100 * ( (t0 + x$cutpoints[i+1]-1) - rw_enddd ) / ( x$cutpoints[i+1] - x$cutpoints[i] ) >= rw_observed_percentage )  
          
          if(event_in_rw) cond <- cond | ( data[,event]==1 &  t0+x$cutpoints[i] <= data[,event_time] & data[,event_time] <= t0+x$cutpoints[i+1]-1  ) 
          
          if(length(censored_vars)>0) 
            for(icensor in censored_vars)
              cond <- cond | ( !is.na(data[,icensor]) & rw_enddd == tend & tend == data[,icensor] & (t0+x$cutpoints[i]) <= tend )
        }
      }     
      else {
        if(i <length(x$cutpoints)) rw_enddd <- t0 + x$cutpoints[i+1]-1
        if(i==length(x$cutpoints)) rw_enddd <- rep(NA,nrow(data))
        cond <- rep(T,nrow(data))
      }  
      
      if( any(cond)) 
        data_tmp <- cbind.data.frame( data,
                                      rw_start = t0 + x$cutpoints[i],
                                      rw_end   = rw_enddd,
                                      lab      = x$lab[i]    )[cond,]
      
      if(any(ls()=="data_tmp")){
        if(any(ls()=="data_tmp0")) data_tmp0 <- rbind.data.frame(data_tmp0, data_tmp)  
        else                       data_tmp0 <- data_tmp
      }
    } # end of 'for'
    
    if(!any(ls()=="data_tmp0")) return(NULL)
    
    rownames(data_tmp0) <- NULL
    
    data_tmp0
  }   # end of sub-function 'calc_wind'
  
  # create and calculate variable from formula from list 'rws' and dataset 'data': 'lab' (or another name), 'rw_start' and 'rw_end':
  calc_wind_horiz <- function(x, data){ 
    
    
    
    t0   <- unlist(lapply(as.expression(x$t0),  eval, data))
    tend <- unlist(lapply(as.expression(x$tend),eval, data))
    if(any(is.na(tend))) tend <- pmin(tend, data[,end_obs], na.rm=T)
    
    if(!any(names(x)=="cuts")) x$cuts <- 0
    else x$cuts <- eval(x$cuts)
    
    #####
    # labels create
    if(any(names(x)=="lab")) x$lab <- eval(x$lab)
    else x$lab <- ""
    
    if(length(x$lab) < length(x$cuts)) x$lab <- c(x$lab, rep( x$lab[length(x$lab)],length(x$cuts)-length(x$lab)) )
    if( any(names(x)=="lab_add_interval")) x$lab_add_interval <- eval(x$lab_add_interval)
    else x$lab_add_interval <- T
    if(x$lab_add_interval){
      x$lab[-length(x$lab)] <- paste0( x$lab[-length(x$lab)],"[",x$cuts[-length(x$cuts)],";",x$cuts[-1]-1,"]")
      if(any(names(x)=="tend")) x$lab[length(x$lab)] <- paste0( x$lab[length(x$lab)],">",x$cuts[length(x$lab)]-1)
    }
    if( any(names(x)=="no_last_interval")) x$no_last_interval <- eval(x$no_last_interval)
    else x$no_last_interval <- F
    
    for(i in 1:length(x$cuts)){
      
      if(any(ls()=="data_tmp")) rm(data_tmp)
      
      if(i<length(x$cuts)){
        
        if(!is.null(tend)) {
          rw_enddd <- pmin(tend,t0 + x$cuts[i+1]-1,na.rm=T)
          cond <- (t0+x$cuts[i]) <= tend  # start must be < the end of these intervals 'tend'
          if(rw_observed_percentage>0){
            cond <- cond &   ( 100 -  100 * ( (t0 + x$cuts[i+1]-1) - rw_enddd ) / ( x$cuts[i+1] - x$cuts[i] ) >= rw_observed_percentage )  
            
            if(event_in_rw) cond <- cond | ( data[,event]==1 &  t0+x$cuts[i] <= data[,event_time] & data[,event_time] <= t0+x$cuts[i+1]-1  ) 
            
            if(length(censored_vars)>0) 
              for(icensor in censored_vars)
                cond <- cond | ( !is.na(data[,icensor]) & rw_enddd == tend & tend == data[,icensor] & (t0+x$cuts[i]) <= tend )
          }
        }     
        else {
          rw_enddd <- t0 + x$cuts[i+1]-1
          cond <- rep(T,nrow(data))
        }  
        cond[is.na(cond)] <- rep(F,sum(is.na(cond)))
        
        if( any(cond)) 
          data_tmp <- cbind.data.frame( data,
                                        rw_start = t0 + x$cuts[i],
                                        rw_end   = rw_enddd,
                                        lab      = x$lab[i]    )[cond,]
      }  
      else {
        if(x$no_last_interval) 
          next
        else {
          if(!is.null(tend)) {
            cond <- (t0+x$cuts[i]) <= tend  # start must be < the end of these intervals 'tend'
            cond[is.na(cond)] <- rep(F,sum(is.na(cond)))
            
            if( any(cond) )     #if( any(names(x)=="tend"))
              data_tmp <- cbind.data.frame( data,
                                            rw_start = t0 + x$cuts[i],
                                            rw_end   = tend,
                                            lab = x$lab[i])[cond,]
          }
        }    
      }  
      if(any(ls()=="data_tmp"))
        if(!any(ls()=="data_tmp0")) data_tmp0 <- data_tmp
      else data_tmp0 <- rbind.data.frame(data_tmp0, data_tmp)
    }
    
    if(!any(ls()=="data_tmp0")) return(NULL)
    
    data_tmp0 <- data_tmp0[!is.na(data_tmp0$rw_start),]
    
    data_tmp0
  }   # end of sub-function 'calc_wind_horiz' for horizontal data
  
  
  if(!is.list(rws)) rws <- eval(parse(text = rws))
  rws <- lapply(rws, function(x) {
    names(x)[ substring(names(x),1,3)=="cut" ] <- "cutpoints"
    if(is.null(x$lab_add_interval)) x$lab_add_interval <- T
    if(is.null(x$no_last_interval)) x$no_last_interval <- F
    if(is.null(x$overlap_priority)) x$overlap_priority <- "next"
    x
  } )
  
  if(any(lapply(rws, function(x) x$overlap_priority)=="next")){
    if(is.null(obj$data_parameters$next_vax_time)){
      if(any(names(data)=="next_vax_time"))   obj$data_parameters$next_vax_time <- "next_vax_time"
      else {
        if(any(names(data)=="next_vax_days")) obj$data_parameters$next_vax_time <- "next_vax_days"
        else                                  obj$data_parameters$next_vax_time <- "next_vax_time"
      }
    }
    
    if(!any(names(data)==obj$data_parameters$next_vax_time)){
      data <- data[order(data[,obj$data_parameters$id],data[,obj$data_parameters$vax_time]),]
      data[,obj$data_parameters$next_vax_time] <- c(data[-1,obj$data_parameters$vax_time], NA)
      data[ c( data[-nrow(data),obj$data_parameters$id] != data[-1,obj$data_parameters$id], T)   ,obj$data_parameters$next_vax_time] <- NA
      obj$data_parameters$next_vax_time <- obj$data_parameters$next_vax_time
    }
  }
  
  rws <- lapply(rws,function(x){ if(any(names(x)=="vax_dep"))if(is.null(names(x$vax_dep))) names(x$vax_dep)<-"before"; x} )
  
  
  if( is.null(names(rws)) ) names(rws) <- lapply(rws, function(x)x$name)
  
  if( is.null(names(rws)) ) names(rws) <- lapply(rws, function(x) gsub( ".", "_", gsub( " |,|-|_", "", tolower(x$lab[1]) ), fixed=T ) )
  if( any(names(rws)=="") ) names(rws)[names(rws)==""] <- lapply(rws, function(x) gsub( ".", "_", gsub( " |,|-|_", "", tolower(x$lab[1]) ), fixed=T ) )[names(rws)==""] 
  
  
  if(!is.null( tmp<-unlist( lapply(rws,function(x)x$name)) )) names(rws) <- tmp
  
  
  
  data$id_n <- 1:nrow(data)
  
  rws_datalist <- lapply( rws, calc_wind, data,  obj )
  
  # create one dataset 'data_rws' from list of dataset 'rws_datalist':
  data_rws <- do.call("rbind.data.frame",rws_datalist)
  
  rownames(data_rws) <- NULL
  
  
  # check:
  if(!is.null(obj$data_parameters$start_obs) & !is.null(obj$data_parameters$end_obs)) 
    if(any(data_rws[,obj$data_parameters$start_obs] >= data_rws[,obj$data_parameters$end_obs])) stop("'start_obs' should be < 'end_obs'")
  
  
  # 1. fill in missing values in 'rw_end' 
  # 2. if necessary check 'rw_observed_percentage'
  # 3. delete [parts of] intervals if interval after 'end_obs'.
  if(!is.null(obj$data_parameters$end_obs)){
    
    # 1. fill in missing values in 'rw_end':
    if(any( (cond <- is.na(data_rws$rw_end) & data_rws$rw_start <= data_rws[,obj$data_parameters$end_obs]) )) 
      data_rws$rw_end[cond] <-  data_rws[cond,obj$data_parameters$end_obs]
    
    # 2. ??? if necessary check 'rw_observed_percentage'
    if(F){  # ???
      #if(i<length(x$cutpoints) & rw_observed_percentage>0){  # 
      cond <- cond &   ( 100 -  100 * ( (t0 + x$cutpoints[i+1]-1) - rw_enddd ) / ( x$cutpoints[i+1] - x$cutpoints[i] ) >= rw_observed_percentage )  
      
      if(event_in_rw) cond <- cond | ( data[,event]==1 &  t0+x$cutpoints[i] <= data[,event_time] & data[,event_time] <= t0+x$cutpoints[i+1]-1  ) 
      
      if(length(censored_vars)>0) 
        for(icensor in censored_vars)
          cond <- cond | ( !is.na(data[,icensor]) & rw_enddd == tend & tend == data[,icensor] & (t0+x$cutpoints[i]) <= tend )
    }
    
    # 3. delete [parts of] intervals if interval after 'end_obs':
    data_rws <- data_rws[ data_rws$rw_start <= data_rws[,obj$data_parameters$end_obs],  ]
    cond <-   data_rws$rw_start <= data_rws[,obj$data_parameters$end_obs] & data_rws[,obj$data_parameters$end_obs] < data_rws$rw_end & 
      !is.na(data_rws$rw_start) & !is.na(data_rws$rw_end) & !is.na(data_rws[,obj$data_parameters$end_obs])
    data_rws$rw_end[cond] <- data_rws[cond, obj$data_parameters$end_obs]
  }
  
  # delete [parts of] intervals if interval before 'start_obs'.
  if(!is.null(obj$data_parameters$start_obs)){ 
    data_rws <- data_rws[ data_rws[,obj$data_parameters$start_obs] <= data_rws$rw_end, ]
    cond <-     data_rws$rw_start < data_rws[,obj$data_parameters$start_obs] & data_rws[,obj$data_parameters$start_obs] <= data_rws$rw_end & 
      !is.na(data_rws$rw_start) & !is.na(data_rws$rw_end) & !is.na(data_rws[,obj$data_parameters$start_obs])
    data_rws$rw_start[cond] <- data_rws[cond, obj$data_parameters$start_obs]
  }
  
  if(strata_cond) {
    no_strata_name <- paste0("no ",ifelse(!is.na(strata_value),strata_value,"stratum"))
    data_rws[!data_rws$strata_cond,"lab"] <- no_strata_name
  }
  else  no_strata_name <- ""
  
  labels   <- data_rws$lab[!duplicated(data_rws$lab)]
  
  if(strata_cond) labels <- c(labels[labels!=no_strata_name],no_strata_name)
  
  # vax-dependend analysis (logical: T or F):    
  lvax_dep <- any( sapply(rws,function(x)!is.null(x$vax_dep)))
  
  if( !lvax_dep ){# no vax_dep variable ){
    
    if(is.null(rws$ref)) {
      if( !is.null( lapply(rws,function(x) is.null(x$ref))) ) 
        rws$ref <- unlist(lapply(rws,function(x) x$ref))[1]
      else rws$ref <- data_rws$lab[order(data_rws$rw_start)][1]
    }
    
    data_rws$lab <- factor_ref(  data_rws$lab, 
                                 #lab_orders = lab_orders,
                                 lab=labels, 
                                 lab_sort=F,
                                 ref=rws$ref 
                                 #event_var=data_rws[,event]
    ) 
    levels(data_rws$lab) <- format(levels(data_rws$lab))
  }
  else {  # with 'brand' or other vax_dependent_variable (e.g., vax_dep <- "vax_brand_short"  "type_vax_short")
    
    if(is.null(rws$ref)) {
      if( !is.null( lapply(rws,function(x) is.null(x$ref))) ) 
        rws$ref <- unlist(lapply(rws,function(x) x$ref))[1]
      else rws$ref <- data_rws$lab[order(data_rws$rw_start)][1]
    }
    
    
    ref_place <- data_rws$rws_cat[ grepl( unlist(lapply(rws[sapply(rws,length)>1],function(x) x$ref))[1]  ,  data_rws$lab ) ][1] # i.e., probably "v0" 
    #ref_place <- data_rws$rws_cat[ grepl(rws$ref,data_rws$lab) ][1] # i.e., probably "v0" 
    
    
    
    sep_vars <- vector("list",length=sum(sapply(rws,length)>1)); names(sep_vars) <- paste0("vars_for_",names(rws)[sapply(rws,length)>1])
    for(irw in (1:length(rws))[sapply(rws,length)>1]){
      
      x <- rws[[irw]]
      if(!any(names(data_rws) %in% x$vax_dep)) stop(paste0("The vaccine dependent variable '", x$vax_dep, "' not found. Check the risk windows definitions."))
      
      if( x$name != ref_place ){
        
        # add variable with the rws$vax_dep (e.g., variable with "Pfize","Moder",...) for x$name (e.g., "v1" or "v2" or "b1")
        tmp <- unique(data_rws[data_rws$rws_cat == x$name,c(obj$data_parameters$id, x$vax_dep)])
        names(tmp)[match(x$vax_dep,names(tmp))] <- paste0(x$name,"_",x$vax_dep) 
        if(any(duplicated(tmp[,obj$data_parameters$id]))) stop(paste0("There are more then one '", x$vax_dep, "' for dose name '", x$name, "'"))
        #dim(data_rws); sum(cond); dim(tmp)
        data_rws <- merge.data.frame( data_rws, tmp, by=obj$data_parameters$id, all.x= T )
        #dim(data_rws); sum(data_rws$rws_cat %in% c(ref_place,x$name))
        
        
        # create new variable for category x$name (e.g,"v1" or "v2" or "b1" ==> 'v1_lab' or ...)
        #     create new labels for all intervals corresponding to category with reference interval (probably, "v0" with "pre-exposure","buffer") 
        
        data_rws[, paste0(x$name,"_lab")] <- ""
        data_rws$tmp_cond <- rep(T,nrow(data_rws))
        
        # create new labels for all categories where ref_cat (e.g., for all categories before the first dose, because one of the categories is reference category)
        cond <- data_rws$rws_cat == ref_place
        
        tmp  <- ""
        for(iplace in x$vax_dep[names(x$vax_dep)=="before"]){
          data_rws$tmp_cond[cond] <- data_rws$tmp_cond[cond] & !is.na(data_rws[cond,paste0(x$name,"_",iplace)])
          tmp  <- paste0( tmp, data_rws[cond,paste0(x$name,"_",iplace)], " " )
        }  
        tmp  <- paste0( tmp ,ifelse( length(x$lab)==1, x$lab, x$name ), data_rws[cond,"lab"])
        
        for(iplace in x$vax_dep[names(x$vax_dep)=="after"]){
          data_rws$tmp_cond[cond] <- data_rws$tmp_cond[cond] & !is.na(data_rws[cond,paste0(x$name,"_",iplace)])
          tmp  <- paste0( tmp, " ", data_rws[cond,paste0(x$name,"_",iplace)] )
        }
        data_rws[cond & data_rws$tmp_cond, paste0(x$name,"_lab")] <- tmp[ data_rws$tmp_cond[cond] ]
        
        #     create new labels for intervals with x$name (e.g, "v2")
        cond <- data_rws$rws_cat ==x$name
        
        tmp  <- ""
        for(iplace in x$vax_dep[names(x$vax_dep)=="before"])
          tmp  <- paste0( tmp, data_rws[cond,paste0(x$name,"_",iplace)], " " )
        tmp  <- paste0( tmp, data_rws[cond,"lab"] )
        
        for(iplace in x$vax_dep[names(x$vax_dep)=="after"])
          tmp  <- paste0( tmp, " ", data_rws[cond,paste0(x$name,"_",iplace)] )
        data_rws[cond, paste0(x$name,"_lab")] <- tmp
        
        
        
        
        
        
        # create new labels (e.g., "no v2") for intervals without x$name (e.g, "v2")
        data_rws[ data_rws[, paste0(x$name,"_lab")]=="" , paste0(x$name,"_lab")] <- paste0("no ",x$name)
        #data_rws[!( data_rws$rws_cat %in% c(ref_place,x$name) ) | (is.na(data_rws[ ,paste0(x$name,"_",x$vax_dep)]) ) , paste0(x$name,"_lab")] <- paste0("no ",x$name)
        if(strata_cond) data_rws[!data_rws$strata_cond,paste0(x$name,"_lab")] <- no_strata_name
        
        
        # create new 'brand'-variables from this variable (e.g., 'v2_lab'  ==> 'vaxdep1_v2_lab', 'vaxdep2_v2_lab', ...)
        data_rws[,"tmp_var"] <- data_rws[,paste0(x$name,"_",x$vax_dep[1])]
        if(length(x$vax_dep)>1) 
          for(i in 2:length(x$vax_dep)) data_rws[,"tmp_var"] <- paste( data_rws[,"tmp_var"], "&", data_rws[,paste0(x$name,"_",x$vax_dep[i])])
        #data_rws[ data_rws[,paste0(x$name,"_lab")] == paste0("no ",x$name), "tmp_var"] <- NA
        
        vax_dep_cat    <-  unique(data_rws[ !(data_rws[,paste0(x$name,"_lab")] %in% c(paste0("no ",x$name),no_strata_name) ),"tmp_var"])
        if(length(vax_dep_cat)>0)
          for(ibr in 1:length(vax_dep_cat) ){
            if(any(data_rws[data_rws$rws_cat==x$name,"tmp_var"]==vax_dep_cat[ibr])){
              
              cond <- data_rws[ ,"tmp_var"]==vax_dep_cat[ibr]  & !(data_rws[ ,paste0(x$name,"_lab")] %in% c(paste0("no ",x$name),no_strata_name))
              data_rws[  cond, paste0("vaxdep",ibr,"_",x$name,"_lab")] <- data_rws[  cond, paste0(x$name,"_lab")]
              
              data_rws[ !cond, paste0("vaxdep",ibr,"_",x$name,"_lab")] <- paste0("no ",substring(vax_dep_cat[ibr],1,5),",",x$name, ifelse(no_strata_name=="","",paste0(",",substring(no_strata_name,4))) )  # e.g."no Pfize or v2"
              
              data_rws[, paste0("vaxdep",ibr,"_",x$name,"_lab")] <- factor_ref(  data_rws[, paste0("vaxdep",ibr,"_",x$name,"_lab")], lab_orders=lab_orders,  ref=x$ref )  
              labels <- levels(data_rws[, paste0("vaxdep",ibr,"_",x$name,"_lab")])
              labels <- c(labels[!(tolower(substring(labels,1,3)) %in% c("no ","no_"))], labels[tolower(substring(labels,1,3)) %in% c("no ","no_")])
              data_rws[, paste0("vaxdep",ibr,"_",x$name,"_lab")] <- factor_ref(  data_rws[, paste0("vaxdep",ibr,"_",x$name,"_lab")], lab=labels, lab_sort=F,  ref=x$ref )  
              
              sep_vars[[irw]] <- c( sep_vars[[irw]], paste0("vaxdep",ibr,"_",x$name,"_lab") )
            }
          }
        data_rws[,"tmp_var"] <- NULL
      }
    }
    sep_vars <- sep_vars[sapply(sep_vars,length)>0]
  }
  
  data_rws <- data_rws[ order(data_rws$id_n) , ]
  data_rws$id_n <- NULL
  
  if(lvax_dep) data_rws <- list( data_rws = data_rws,
                                 sep_vars = sep_vars       )
  
  data_rws
} # end of function 'create_rws'  (version 3) 




#########################################
#########################################
#########################################
#
#  add time dependent variable  => split some intervals
#


split_intervals <- function(  data,                            # dataset
                              start_interval, end_interval,    # two variables in 'data'
                              splits_names,                    # name of the new variable
                              splits,                          # 1. a number; or 2. vector with numbers; or 
                              # 3. a variable in 'data'; or 4. 'cbind' or 'cbind.data.frame' of variables in 'data'
                              lab = c("before","after"), # labels for split intervals. The length should be (#splits + 1)
                              # if  2 intervals ==> c("before","after")
                              # if  3 intervals ==> c("before","during","after")
                              # if >3 intervals ==> paste0("(",time_seq[-length(time_seq)],";",time_seq[-1],"]")
                              lab_add_interval = T,            # add intervals at the end of labels 'lab'
                              lab_orders       = NA,
                              ref=1,                           # reference category for new variable: number OR category name OR "most events"
                              event,                           # used if ref=="most events" to define reference category with most events
                              event_time                       # used if ref=="most events" to define reference category with most events
                              #event_var = substitute(event)    # used if ref=="most events" to define reference category with most events
){
  
  # splits can be a value   or   a vector of values    or  a variable     or a vector of variables:
  if( is.data.frame(splits)){   # splits are variables from 'data' 
    names(splits) <- paste0("_split_",1:ncol(splits))
    
    data <- cbind.data.frame(data, as.data.frame(splits))
    splits <- names(splits)
  }
  else splits <- splits[ min( data[,start_interval], na.rm=T) < splits & splits < max( data[,end_interval], na.rm=T) ]
  
  
  data$i_ <- 1:nrow(data)
  
  # split intervals:
  for(isplit in 1:length(splits)){
    if(mode(splits)=="character") 
      split <- data[,paste0("_split_",isplit)] # for 4. 'cbind' or 'cbind.data.frame' of variables in 'data'
    else  split <- splits[[isplit]]   #for:  # 1. a number; or 2. vector with numbers; or # 3. a variable in 'data';
    
    start_int <- data[,start_interval]   # eval(substitute(data$start_interval))
    end_int   <- data[,end_interval  ]   # eval(substitute(data$end_interval))
    
    # "before" the split: entire intervals to the left of the split
    if(isplit==1) {
      data[,splits_names] <- NA
      data[ end_int <  split & !is.na(end_int), splits_names] <- isplit
      data[ end_int <  split & !is.na(end_int), splits_names] <- isplit
    }
    
    # "after" the split: entire intervlas to the right of the split
    data[ split <= start_int & !is.na(start_int) & !is.na(split) , splits_names] <- isplit + 1
    
    #######
    # split interval into 2 :
    cond <-  start_int < split & split <= end_int & !is.na(split) & !is.na(start_int) & !is.na(end_int)
    if(sum(cond)>0){
      if(length(split)==1) split <- rep(split,sum(cond)) 
      else                 split <- split[cond]
      
      # 'after' the split: part of the interval to the right of the split
      data_tmp <- data[cond,]
      data_tmp[ , splits_names ]   <- isplit + 1
      data_tmp[ , start_interval ] <- split                     #  deparse(substitute(start_interval)) ] <- split
      
      # 'before' the split: part of the interval to the left of the split
      if(isplit==1) data[ cond, splits_names ] <- 1
      data[ cond, end_interval ] <- split - 1    #      deparse(substitute(end_interval)) 
      
      data <- rbind.data.frame(data, data_tmp)
    }
  } 
  
  if(mode(splits)=="character") data <- data[,!( names(data) %in% paste0("_split_",1:length(splits)) ) ] 
  ####
  # add labels
  data[,splits_names] <- as.factor(data[,splits_names])
  if(missing(lab)){ 
    if(mode(splits)=="character"){
      if(nlevels(data[,splits_names])==2) lab <- c("before","after") 
      if(nlevels(data[,splits_names])==3) lab <- c("before","during","after") 
    }
    if(mode(splits)=="numeric")
      lab <- paste0("[",c(min(start_int,na.rm=T), splits),";",c(splits-1, max(end_int,na.rm=T)),"]")
  }  
  
  
  
  data              <- refresh_event_variable( start_interval, end_interval, data, event, event_time)
  if(length(lab_orders)==1 & is.na(lab_orders[1]))
    data[,splits_names] <- factor_ref(  data[,splits_names], lab=lab, ref=ref, event_var=data[,event] ) 
  else
    data[,splits_names] <- factor_ref(  data[,splits_names], lab_orders = lab_orders, lab=lab, ref=ref, event_var=data[,event] )  
  
  
  data <- data[order( data$i_, data[,start_interval] ), ]   #   eval(substitute(data$start_interval))),] 
  data$i_ <- NULL
  data
}  # end of function 'split_intervals_horiz'

# data_rws$interval <- data_rws$rw_end - data_rws$start_rw + 1
# 
# data_rws2 <- split_intervals(data_rws, start_rw, rw_end, "tdvar", interval+50, ref=2 )
# data_rws2 <- split_intervals(data_rws, start_rw, rw_end, "tdvar", cbind.data.frame(interval, interval+50), c("befo","dur","afte"), ref=2 )
# 
# data_rws2 <- split_intervals(data_rws, start_rw, rw_end, "tdvar", 170, ref=2 )
# data_rws2 <- split_intervals(data_rws, start_rw, rw_end, "tdvar", c(150,170), c("b","d","after") )
# 
# data_rws2 <- data_rws2[order(data_rws2$person_id,data_rws2$start_rw),]  
# data_rws2 
# 








################################################################
################################################################
################################################################
#    summary_tab
# summary per risk/control windows and results from Poisson regression
# in function 'var_sum': 'var_name', in 'data' must be variables: 'event','enterval',id (can be changed)
summary_tab <- function(  var_names_data, # var_names <- c("lab", "cal_time_cat")
                          var_names_mod,
                          event,
                          data, id_name = "id",
                          mod,
                          lab_orders,
                          #coef_cond, 
                          model_number=1,
                          res_tab  #,
                          #print = F, 
                          #digits = 3
){
  add <- F
  if(!missing(res_tab)) { res_tab0 <- res_tab; add <- T}
  # check:
  if(missing(var_names_data)) stop("'var_names_data' is missing")
  if(missing(data) & missing(mod)) stop("'data' and/or 'mod' must be specified. ")
  
  ######
  # create summary table for variable: 'var_names_data' from 'data' (variable 'interval', 'event',id_name must be in 'data')
  if(!missing(data)){
    
    cond <- !is.na(data$interval) & data$interval>0
    res_tab <- vector("list",length(var_names_data))
    for(ivar in var_names_data){ 
      
      if(length(table(data[,ivar]))<=1) next
      if( is.factor(data[,ivar]) ){ 
        if(any(cond2<-rowSums(contrasts(data[,ivar]))==0)) ivar_ref_cat <- levels(data[,ivar])[cond2]
        else  ivar_ref_cat <- levels(data[,ivar])[1]   # here m.b. study other contrasts (for the future?)
      }
      else    ivar_ref_cat <- levels(data[,ivar])[1]
      
      n_events_per_id_ref <-  table( data[ cond &  data[,ivar]==ivar_ref_cat & data[,event]>0 , id_name] )
      ids_ref_names       <- unique( data[ cond &  data[,ivar]==ivar_ref_cat                  , id_name] )
      
      events_ref <- rep(NA,nlevels(data[,ivar])); names(events_ref) <- levels(data[,ivar])
      ids_rw  <- ids_ref <- observed_rw  <- observed_ref <- events_rw <- median_days_per_id_rw <- events_ref 
      ids_rw_noref <- observed_rw_noref  <- events_rw_noref <- median_days_per_id_rw_noref <- events_ref 
      for(ilab in levels(data[,ivar])){
        #  ids with events in the risk window:
        ids_event_rw_names  <- unique( data[ cond &  data[,ivar]==ilab & data[,event]>0 , id_name] )
        #  ids with events in the risk window but not observed in the ref.window:
        ids_event_rw_and_not_observed_ref <- ids_event_rw_names[ !(ids_event_rw_names %in% ids_ref_names) ]
        
        #  ids with events in the risk window AND observed in the ref.window:
        ids_event_rw_and_observed_ref <- ids_event_rw_names[ ids_event_rw_names %in% ids_ref_names ]
        #  ids with events in ref.window and observed in risk window:
        ids_event_ref_and_observed_rw <- names(n_events_per_id_ref)[ names(n_events_per_id_ref) %in% data[ cond & data[,ivar]==ilab, id_name]  ]
        
        #  both:
        if(ilab==ivar_ref_cat) ids_rw_names       <-    ids_event_rw_and_observed_ref                                     # ids_event_rw_names 
        else                   ids_rw_names       <- c( ids_event_rw_and_observed_ref, ids_event_ref_and_observed_rw )    # c( ids_event_rw_names, ids_event_ref_and_observed_rw )
        
        #  #events, #ids with events and #observed days in the risk window (and observed in the ref.window):
        if(length(ids_event_rw_and_not_observed_ref)>0){
          events_rw_noref[            ilab] <- sum(           data[ cond &  data[,ivar]==ilab & data[,id_name] %in% ids_event_rw_and_not_observed_ref                   ,   event    ] )
          ids_rw_noref[               ilab] <- length(unique( data[ cond &  data[,ivar]==ilab & data[,id_name] %in% ids_event_rw_and_not_observed_ref &  data[,event]>0 ,   id_name  ] ) )
          observed_rw_noref[          ilab] <- sum(           data[ cond &  data[,ivar]==ilab & data[,id_name] %in% ids_event_rw_and_not_observed_ref                   ,  "interval"] )
          median_days_per_id_rw_noref[ilab] <- median(   with(data[ cond &  data[,ivar]==ilab & data[,id_name] %in% ids_event_rw_and_not_observed_ref, c(id_name,"interval")], tapply(interval,get(id_name),sum)) )
        }
        #  #events, #ids with events and #observed days in the risk window (and observed in the ref.window):
        events_rw[            ilab] <- sum(           data[ cond &  data[,ivar]==ilab & data[,id_name] %in% c( ids_rw_names, ids_event_rw_and_not_observed_ref)                    ,   event    ] )
        ids_rw[               ilab] <- length(unique( data[ cond &  data[,ivar]==ilab & data[,id_name] %in% c( ids_rw_names, ids_event_rw_and_not_observed_ref)  &  data[,event]>0 ,   id_name  ] ) )
        observed_rw[          ilab] <- sum(           data[ cond &  data[,ivar]==ilab & data[,id_name] %in% c( ids_rw_names, ids_event_rw_and_not_observed_ref)                    ,  "interval"] )
        median_days_per_id_rw[ilab] <- median(   with(data[ cond &  data[,ivar]==ilab & data[,id_name] %in% c( ids_rw_names, ids_event_rw_and_not_observed_ref) , c(id_name,"interval")], tapply(interval,get(id_name),sum)) )
        #  #events, #ids with events and #observed days in ref.window (and observed in the risk windows):
        if(ilab!=ivar_ref_cat){ 
          events_ref[  ilab] <- sum(    n_events_per_id_ref[ names(n_events_per_id_ref) %in% ids_rw_names ] ) 
          ids_ref[     ilab] <- length( n_events_per_id_ref[ names(n_events_per_id_ref) %in% ids_rw_names ] ) 
          observed_ref[ilab] <- sum( data[ cond &  data[,ivar]==ivar_ref_cat &  data[,id_name] %in% ids_rw_names, "interval"] )
        }
        #  #events, #ids with events and #observed days in the risk window and NOT observed in the ref.window:
      }
      
      res_tab[ivar==var_names_data][[1]] <- cbind.data.frame(
        
        events_rw  = as.integer(events_rw),  # with(data[cond & data[,event] >0, ], tapply(get(id_name),get(ivar),function(x) sum( x %in% ids_pre_names ))),
        events_ref = events_ref,
        #events_ref = with(data[cond & data[,event]==0,], tapply(get(id_name),get(ivar),function(x) sum( n_events_per_id_pre[ names(n_events_per_id_pre) %in% unique(x) ])   )),
        
        days_rw    = observed_rw ,
        days_ref   = observed_ref, 
        
        ids_rw     = ids_rw,   #with(data[cond & data[,event] >0,], tapply(get(id_name),get(ivar),function(x)sum(unique(x) %in% ids_pre_names))),
        ids_ref    = ids_ref,  #with(data[cond & data[,event]==0,], tapply(get(id_name),get(ivar),function(x)sum(unique(x) %in% names(n_events_per_id_pre))))
        
        median_days_per_id_rw       = median_days_per_id_rw,
        median_days_per_id_rw_noref = median_days_per_id_rw_noref,  # in risk window but NOT observed in the reference window:
        
        # in risk window but NOT observed in the reference window:
        events_rw_noref = events_rw_noref,
        ids_rw_noref    = ids_rw_noref,
        days_rw_noref   = observed_rw_noref
      )
    }
    res_tab <- do.call("rbind.data.frame",res_tab)
    res_tab <- cbind.data.frame( i=1:nrow(res_tab), all_cat=rownames(res_tab) , res_tab)
    
    res_tab$events_ref[!is.na(res_tab$events_ref) & res_tab$events_ref==0] <- NA      
    res_tab$ids             <- rowSums(res_tab[,c("ids_rw", "ids_ref")], na.rm=T)
    
    res_tab$days_per_id_rw  <- round( res_tab$days_rw  / res_tab$ids , 2 )
    res_tab$days_per_id_ref <- round( res_tab$days_ref / res_tab$ids , 2 )
    
    # in risk window but NOT observed in the reference window:
    res_tab$days_per_id_rw_noref  <- round( res_tab$days_rw_noref  / res_tab$ids_rw_noref , 2 )
    
    res_tab$rw_portion           <- res_tab$events_rw  / res_tab$days_rw
    res_tab$ref_portion          <- res_tab$events_ref / res_tab$days_ref
    res_tab$RR_data              <- res_tab$rw_portion / res_tab$ref_portion
    
    res_tab <- res_tab[,c( "i",
                           "all_cat", "events_rw", "events_ref","days_per_id_rw", "days_per_id_ref", "median_days_per_id_rw", 
                           "days_rw", "days_ref", "ids_rw", "ids_ref", "ids", 
                           "events_rw_noref", "ids_rw_noref", "days_rw_noref", "days_per_id_rw_noref", "median_days_per_id_rw_noref",    # these are added later 
                           "rw_portion", "ref_portion", "RR_data" )]
  }
  
  ######
  # create summary table for result from Poisson regression with variables 'var_names_mod':
  #
  if(!missing(mod)){ 
    res_tab_model <- cbind.data.frame( summary(mod)$conf.int[,-2], 
                                       summary(mod)$coef[,match( c("Pr","co","se"), substring(dimnames(summary(mod)$coef)[[2]],1,2) ) ]  ) #  "Pr(>|z|)" "coef"   "se(coef)" 
    res_tab_model[ is.na(res_tab_model$coef), substring(dimnames(res_tab_model)[[2]],1,2)=="se" ] <- NA
    # delete variable names at the beginning of parameter names in results from Poisson regression:
    res_tab_model$all_cat <- dimnames(res_tab_model)[[1]]
    var_names_mod <- unique(unlist(strsplit(var_names_mod,":")))
    var_names_mod <- var_names_mod[order(nchar(var_names_mod),decreasing = T)]
    for(ivar in var_names_mod){
      res_tab_model$all_cat[cond] <- substring(res_tab_model$all_cat, nchar(ivar)+1)[(cond<-substring(res_tab_model$all_cat,1,nchar(ivar))==ivar)]
      res_tab_model$all_cat       <- gsub( paste0(":",ivar), ":", res_tab_model$all_cat)
    }
    names_matched <- match(c("exp","low","upp","Pr(","se("),substring(names(res_tab_model),1,3))
    
    if(any(!is.na(names_matched)))
      names(res_tab_model)[names_matched[!is.na(names_matched)]] <- c("RR","lci","uci","pval","se_coef")[!is.na(names_matched)]
    #names(res_tab_model)[match(c("exp","low","upp","Pr("),substring(names(res_tab_model),1,3))]] <- c("RR","2.5%","97.5%","pval")
    
    res_tab_model$model <- model_number
    
  }
  
  #####
  # combine (merge) these two tables:
  if(!missing(data) &  !missing(mod)){
    res_tab$all_cat_tmp       <- gsub(" ","",res_tab$all_cat)
    res_tab_model$all_cat_tmp <- gsub(" ","",res_tab_model$all_cat)
    res_tab <- merge.data.frame( res_tab, res_tab_model, by="all_cat_tmp", all=T, sort=F )
    res_tab <- res_tab[, names(res_tab)!="all_cat_tmp" ]
    
    if( any(names(res_tab)=="all_cat.x") ){
      if( any(names(res_tab)=="all_cat.y") ){ 
        res_tab$all_cat.x[ is.na(res_tab$all_cat.x) & !is.na(res_tab$all_cat.y) ] <- res_tab$all_cat.y[ is.na(res_tab$all_cat.x) & !is.na(res_tab$all_cat.y) ]
        res_tab <- res_tab[, names(res_tab)!="all_cat.y" ]
      }
      names(res_tab)[names(res_tab)=="all_cat.x"] <- "all_cat"
    }  
    res_tab <- res_tab[order(res_tab$i),]
  } 
  
  if(missing(data) & !missing(mod)) { 
    res_tab <- res_tab_model
    res_tab <- cbind.data.frame( i=1:nrow(res_tab), res_tab[, c(ncol(res_tab),1:(ncol(res_tab)-1))])
  }
  
  res_tab <- cbind.data.frame( event=substring(event,1,7), res_tab )
  
  res_tab$model <- model_number
  if(add){
    
    model_res_names <- c("RR","lci","uci","pval","coef","se_coef","model")
   
    res_tab_new <- res_tab0
    res_tab_new$all_cat_without_space <- gsub(" ","",res_tab_new$all_cat)
    res_tab$all_cat_without_space     <- gsub(" ","",res_tab$all_cat)
    res_tab_new <- merge.data.frame(res_tab_new,res_tab, by=c("all_cat_without_space"), all=T, sort=F )
    
    names(res_tab_new)[substring(names(res_tab_new),nchar(names(res_tab_new))-1,nchar(names(res_tab_new)))==".x"] <- substring(names(res_tab_new),1,nchar(names(res_tab_new))-2)[substring(names(res_tab_new),nchar(names(res_tab_new))-1,nchar(names(res_tab_new)))==".x"]
    for(ivar in c(names(res_tab)[ !(names(res_tab) %in% c("all_cat_without_space",model_res_names))]) ){  # c( "event","all_cat","n_events","atrisk_days","atrisk_ids","days_pp","relative_rate","relative_perc")){
      #for(ivar in c( "event","all_cat","n_events","atrisk_days","atrisk_ids","days_pp","relative_rate","relative_perc")){
      #names(res_tab_new)[names(res_tab_new)==paste0(ivar,".x")] <- ivar
      if(any( !is.na(res_tab_new[,"model"]) & !is.na(res_tab_new[,"model.y"]) & res_tab_new[,ivar]!=res_tab_new[,paste0(ivar,".y")] & 
              ( !is.na(res_tab_new[, ivar])  | !is.na(res_tab_new[,paste0(ivar,".y")]) ) & substring(res_tab_new$all_cat,1,1)!="[" )) next
      
      res_tab_new[is.na(res_tab_new[,"model"]) & !is.na(res_tab_new[,"model.y"]) ,ivar] <- res_tab_new[ is.na(res_tab_new[,"model"]) & !is.na(res_tab_new[,"model.y"]), paste0(ivar,".y")]
      res_tab_new[,paste0(ivar,".y")] <- NULL
    }
    if( "model.y" %in% names(res_tab_new) & !("RR.y" %in% names(res_tab_new))) res_tab_new <- res_tab_new[,names(res_tab_new)!="model.y"]
    
    if(!missing(mod)){
      model_res_names <- c("RR","lci","uci","pval","coef","se_coef","model")
      #model_res_names <- c("i","RR","2.5%","97.5%","pval","coef","se(coef)","model")
      
      var_numbers <- match(model_res_names,names(res_tab_new))
      var_numbers <- var_numbers[!is.na(var_numbers)]
      names(res_tab_new)[var_numbers] <-  model_res_names[!is.na(var_numbers)]
      
      if(  any( names(res_tab_new)=="RR" ) & any(names(res_tab_new)=="RR.y" ) ){ 
        # duplicated sets of columns with different values
        if(any(  cond<-!is.na(res_tab_new$model.y) & substring(res_tab_new$all_cat,1,1)=="[" )){
          #if(sum( cond<- !is.na( res_tab_new[,"pval"] ) & !is.na( res_tab_new[,paste0("pval.y")]) & res_tab_new[,paste0("pval")]!=res_tab_new[,paste0("pval.y")] )>0){
          res_tab_new_dupl <- res_tab_new[cond, , drop=F]
          res_tab_new_dupl[, model_res_names ] <- NULL
          
          match_res <- match(paste0(model_res_names,".y"),names(res_tab_new_dupl))
          names(res_tab_new_dupl)[match_res[!is.na(match_res)]] <-  model_res_names[!is.na(match_res)]
          if(any(names(res_tab_new_dupl)=="i.y")){
            res_tab_new_dupl <- res_tab_new_dupl[order(res_tab_new_dupl$i.y),]
            res_tab_new_dupl$i.y <- NULL  
          }
          res_tab_new_dupl <- res_tab_new_dupl[!duplicated(res_tab_new_dupl[,names(res_tab_new_dupl)!="i"]),]
        }
        
        if(sum( (cond <- is.na( res_tab_new[,"RR"] ) & !is.na( res_tab_new[,paste0("RR.y")])) )>0) 
          res_tab_new[cond, model_res_names ] <- res_tab_new[cond, paste0(model_res_names,".y") ]
        
        match_res <- match(paste0(model_res_names,".y"),names(res_tab_new))
        res_tab_new[,match_res[!is.na(match_res)]] <-  NULL
      }
    }
    
    if(any(names(res_tab_new)=="i.y")) res_tab_new <- res_tab_new[,names(res_tab_new)!="i.y"]
    if( any(ls()=="res_tab_new_dupl") ){
      match_res <- match(names(res_tab_new),names(res_tab_new_dupl))
      if(any(is.na(match_res))) res_tab_new_dupl[, names(res_tab_new)[is.na(match_res)] ] <- NA
      res_tab_new <- rbind.data.frame(res_tab_new, res_tab_new_dupl[,match(names(res_tab_new),names(res_tab_new_dupl))]  )
    }
    res_tab_new$all_cat_without_space <- NULL
    
    
    if(all(names(res_tab_new)!="model")){
      if(any(names(res_tab_new)=="i.x"))     names(res_tab_new)[names(res_tab_new)=="i.x"    ] <-  "i"
      if(any(names(res_tab_new)=="model.x")) names(res_tab_new)[names(res_tab_new)=="model.x"] <-  "model"
      if(any(cond<-!is.na(match(names(res_tab_new),c("i.y","model.y")))))   res_tab_new[,cond] <-  NULL
      res_tab_new$model <- 1
    }
    
    # order table on 'all_cat'
    res_tab_new$tablab_order <- match( res_tab_new$all_cat, levels(factor_ref(res_tab_new$all_cat[ substring(res_tab_new$all_cat,1,1)!="[" ], lab_orders=lab_orders)) )
    res_tab_new$tablab_order[is.na(res_tab_new$tablab_order)] <- 100000
    res_tab_new <- res_tab_new[order(res_tab_new$tablab_order, res_tab_new$model),]
    #res_tab_new <- res_tab_new[order(res_tab_new$tablab_order, res_tab_new$model, res_tab_new$i),]
    res_tab_new$tablab_order <- NULL
    res_tab_new[,"ids"]   <- as.integer(res_tab_new[,  "ids"])
    res_tab_new[,"model"] <- as.integer(res_tab_new[,"model"])
    
    res_tab <- res_tab_new
    
    res_tab$i <- 1:nrow(res_tab)
    
    
  }
  
  res_tab
}  # end of function 'summary_tab'

#summary_tab(var_names=c("lab", "cal_time_cat"), data=data_rws, id_name="person_id",  print=T, digits=2)
#summary_tab(var_names=c("lab", "cal_time_cat"), data=data_rws, id_name="person_id", mod=mod, print=T,digits=1)
#summary_tab(var_names=c("lab", "cal_time_cat"),  mod=mod, print=T,digits=1)


#########################################
#########################################
#########################################
#
#  add a level in factor variable and choose reference category
#

factor_ref <- function(  var,
                         lab,
                         lab_sort   = F,
                         lab_orders = NA,
                         ref=1,          # reference category for new variable: number OR [beginning of ] category name OR "most events"
                         event_var,      # used if ref=="most events" to define reference category with most events
                         keep_ref = T,  # if 'var' is factor with reference category and 'ref' is missing, then keep this reference category in the updated variable 
                         na_replace = ""
){ 
  
  if( length(lab_orders)>1 | !is.na(lab_orders[1]) ){
    
    # lab_orders <- list( c("pre-","buf", "dose 1",     "dose 2" ),
    #                c("[0;0]","[1;28]",">28"),
    #                c("Pfizer","Moderna","AstraZeneca", ... ) )
    
    if(!is.list(lab_orders)) lab_orders <- list(lab_orders)
    
    var_unique <- unique(var)
    big_number <- max(sapply(lab_orders,length))+10
    
    orders <- as.data.frame( matrix( big_number, nrow=length(var_unique), ncol=length(lab_orders),
                                     dimnames=list(var_unique, 1:length(lab_orders)) ))
    
    for(itype in 1:length(lab_orders))
      for(icat in 1:length(lab_orders[[itype]]))
        if(any( cond<-(grepl( lab_orders[[itype]][icat], var_unique, fixed=T)  &  orders[,itype]==big_number) ))
          orders[ cond, itype] <- icat
    
    
    orders <- eval(parse(text=paste0("orders[ order(", paste( paste0("orders[,", 1:length(lab_orders)),collapse="],") ,"]),]")) )
    orders$order <- 1:nrow(orders)
    
    var <- factor(var, levels=rownames(orders), labels=rownames(orders) )
  }
  else {
    if(is.factor(var)){
      if(nlevels(var)>1) ref0 <- levels(var)[rowSums(contrasts(var))==0]
      levels0 <- levels(var)
      var <- as.character(var)
      if(all(!grepl("[a-z,;:$@|]",tolower(levels0))))
        if( all(!is.na(as.numeric(levels0))) ) var <- as.numeric(var)
    }  
    if(mode(var)=="character") var[is.na(var)] <- ""
    if(missing(lab)) lab <- var[!duplicated(var)] 
    if(lab_sort)  lab <- sort(lab)
    
    if(is.numeric(var)) {
      if(length(lab)== length(unique(var)))
        var <- factor(var, labels=lab)
      else {
        if(max(var,na.rm=T)<=length(lab) & all(!is.na(var)))
          var <- factor(var, labels=lab[sort(unique(var))])
        else
          var <- factor(var)
      }    
    }  
    else var <- factor(var, levels=lab, labels=lab)
  }
  
  # choose reference category ref
  if(!missing(ref)){
    if(is.character(ref)){
      if( any( cond<-grepl(ref,levels(var),fixed=T) ) ) 
        ref <- (1:nlevels(var))[cond][1]
      if(ref=="most events"){
        tb <- table( var[ !is.na(event_var) & event_var>0 ])
        ref <- match(names(tb)[  tb == max(tb) ][1], levels(var))
      }
    }
  }
  else{  # 'ref' is missing
    if(keep_ref & any(ls()=="ref0")) 
      if( ref0 %in% levels(var)) ref <- (1:nlevels(var))[levels(var)==ref0]
  }
  if(is.numeric(ref) & nlevels(var)>1){
    contrasts(var) <- contr.treatment(nlevels(var),base=ref)
    dimnames(contrasts(var))[[2]] <- levels(var)[-ref]
  } 
  var
} # end of function 'factor_ref'

#
###########################################
###########################################
###########################################



###########################################
###########################################
#
#
combine_vars_func <- function(var1, var2, sep=" & ", 
                              lab_orders,
                              ref, event          # reference category for new variable: number OR [beginning of ] category name OR "most events"
){
  
  if(!missing(var2)) var0 <- cbind.data.frame(var1, var2)
  else var0 <- var1
  
  for(icol in 2:ncol(var0)){
    var2 <- var0[,icol]
    var1 <- var0[,icol-1]
    
    ref_cat <- ref_cat1 <- ref_cat2 <- c()  
    if(is.factor(var1)) ref_cat1 <- levels(var1)[rowSums(contrasts(var1))==0]
    if(is.factor(var2)) ref_cat2 <- levels(var2)[rowSums(contrasts(var2))==0]
    
    if(!is.factor(var1)) var1 <- as.factor(var1)
    if(!is.factor(var2)) var2 <- as.factor(var2)
    
    if(length(unique(var1))==1) ref_cat1 <- unique(var1)
    if(length(unique(var2))==1) ref_cat2 <- unique(var2)
    
    all_combi <- paste( rep( levels(var1),  each=nlevels(var2)),  
                        rep( levels(var2),       nlevels(var1)), sep=sep)
    var12 <- paste(var1, var2, sep=sep)
    var_levels <- all_combi[ all_combi %in% unique(var12) ]
    
    # search for ref category:
    if(!missing(ref)){
      if(length(ref)>1){
        if(!is.na(ref[[1]])){
          if(is.numeric(ref[1])) ref_cat1 <- var1[ref[1]]
          if(is.character(ref[1])) {
            if(ref[1]=="most events") ref_cat2 <- names(table(var2))[table(var2)==max(table(var2))]
            else
              if(any(grepl(ref[1],var1)))
                ref_cat1 <- grep(ref[1],var1, value=T)
          }  
        }  
        if(!is.na(ref[[2]])){
          if(is.numeric(ref[2])) ref_cat2 <- var2[ref[2]]
          if(is.character(ref[2])) {
            if(ref[2]=="most events") ref_cat2 <- names(table(var2))[table(var2)==max(table(var2))]
            else  
              if(any(grepl(ref[2],var2))) ref_cat2 <- grep(ref[2],var2, value=T) 
          }    
        }  
      }  
      else {
        if(is.numeric( ref )) ref_cat <- var_levels[ref]
        if(is.character(ref)) {
          if(ref=="most events") ref_cat <- ref
          else  
            if(any(grepl(ref,var_levels)))    ref_cat <- grep(ref, var_levels, value=T)
        }    
      }
    }
    
    if(length(ref_cat1)>0) ref_cat1 <- grep( ref_cat1, var_levels, value=T, fixed=T)
    if(length(ref_cat2)>0) ref_cat2 <- grep( ref_cat2, var_levels, value=T, fixed=T)
    
    if(length(ref_cat)==0){
      if(length(ref_cat1)>0  & length(ref_cat2) >0) ref_cat <- var_levels[ var_levels %in% ref_cat1 & var_levels %in% ref_cat2 ]
      if(length(ref_cat1)>0  & length(ref_cat2)==0) ref_cat <- var_levels[ var_levels %in% ref_cat1                            ]
      if(length(ref_cat1)==0 & length(ref_cat2) >0) ref_cat <- var_levels[                            var_levels %in% ref_cat2 ]
    }  
    
    if(length(ref_cat)>1 & length(ref_cat1)>0){
      tmp <- c()
      for(icat in ref_cat1)
        tmp <- c( tmp, grep(icat, ref_cat, value=T, fixed=T) )
      if(length(tmp)>0) ref_cat <- tmp
    }  
    if(length(ref_cat)>1 & length(ref_cat2)>0){
      tmp <- c()
      for(icat in ref_cat2)
        tmp <- c( tmp, grep(icat, ref_cat, value=T, fixed=T) )
      if(length(tmp)>0) ref_cat <- tmp
    }  
    
    #if(length(ref_cat)>0 & !is.na(ref_cat))
    #  ref <- (1:length(var_levels))[var_levels==ref_cat[1]]
    
    #ref <- (1:nlevels(var_levels))[var_levels==ref_cat[1]]
    if(length(ref_cat)==0) ref_cat <- 1
    
    var12   <- factor_ref( var12, lab_orders = lab_orders, lab=var_levels, ref=ref_cat[1], event_var = event )
    var0[,icol] <- var12
  }
  var12
}   # end of function 'combine_vars'
#table1(data_rws$brand_lab)
#contrasts(data_rws$brand_lab)


###############################
#
plot_res <- function(res, main="", 
                     col=c("red", "green3", "orange",  "deepskyblue", "magenta2", "cyan2", "chocolate1", "gray" ), 
                     time_cat_i=length(strata), ylim, CI = T, 
                     correct_max_time_adj =  T){
  
  plotted <- F
  if(all(unlist(lapply(res, function(x)sum(!is.na(x$RR))))==0)) return(plotted) 
  
  res <- lapply(res,function(x){if(length(dim(x))==2)x<-as.data.frame(x);x})
  
  ncoef     <- nrow(res[[1]])
  ncoef_max <- max(unlist( lapply(res,nrow) ))
  if(ncoef_max>ncoef) # +1 because for ref.category of cal_time_var 'model' is NA
    ncoef_max <- ncoef + 1 + max(unlist( lapply(res,function(x) {
      cond <- (1:nrow(x))[(ncoef+1):nrow(x)]
      if(all(dimnames(x)[[2]]!="model")) return(0)
      if(all(is.na(x[cond,"model"])))    return(0)
      tapply( x[cond,"i"], x[cond,"model"], length ) 
    } )))
  
  if(missing(ylim)){
    ymax   <- unlist(lapply(res,function(x) {
      if(all(dimnames(x)[[2]]!="RR")) return(NA) 
      x_tmp   <- x$RR[x$RR<1000] 
      x_tmp_1 <- x$RR[1:ncoef][x$RR[1:ncoef]<1000] 
      if(any(!is.na(x_tmp))) max(c(x_tmp_1*1.2,x_tmp),na.rm=T) else NA 
    } ))
    
    if(any(!is.na(ymax))) ymax <- max(ymax, na.rm=T)
    else                  ymax <- 10
    ylim <- c( 0, ymax )
  }
  
  
  text_col_cond <- 1 + as.numeric(!is.na(res[[1]]$RR[1:ncoef]))
  tmp_1         <-    tmp_05    <-    rep(F,length(text_col_cond)) 
  for(imod in 1:length(res)){
    if(any(dimnames(res[[imod]])[[2]]=="RR")){
      tmp_1  <-  tmp_1  | ( !is.na(res[[imod]]$RR[1:ncoef]) & !is.na(res[[imod]]$pval[1:ncoef]) & res[[imod]]$pval[1:ncoef]<=0.1  )
      tmp_05 <-  tmp_05 | ( !is.na(res[[imod]]$RR[1:ncoef]) & !is.na(res[[imod]]$pval[1:ncoef]) & res[[imod]]$pval[1:ncoef]<=0.05 )
    }
  }
  text_col_cond <- text_col_cond + as.numeric(tmp_1 )
  text_col_cond <- text_col_cond + as.numeric(tmp_05)
  
  x_deltas <- 0
  if(length(res)>1) x_deltas <- 0.4* (1/(length(res)-1) * (0:(length(res)-1)) - 0.5)
  
  # function for colors:
  col_alpha <- function(col,alpha=0) rgb(t(col2rgb(col)/255),alpha=alpha)
  
  ###########
  #
  # plot 1: all coefficients:
  #
  if(ncoef_max > ncoef){ 

    plot( c(0,ncoef_max), ylim, type="n", main=main, xlab="effect number  |     time adjustment effects", ylab="RR", axes=F)
    axis(2); box()
    axis(1, at=1:ncoef )
    plotted <- T
    
    if(ncoef_max > ncoef+1) { #res[[i]]$all_cat[substring(res[[i]]$all_cat,1,1)=="["]
      suppressWarnings(
        numbers <- as.numeric(unique(unlist(lapply(res, function(x, ncoef) 
          strsplit( x$all_cat[ substring(x$all_cat,1,1)=="[" ], "\\[|;|,|\\]" ), ncoef=ncoef ))))
      )  
      axis(1, at=c(ncoef+0.5, ncoef_max+0.5), labels = paste0(c(min(numbers, na.rm=T), max(numbers, na.rm=T)),"days"), padj=1, tck=-0.04  )
      axis(1, at=c(ncoef+0.5, ncoef_max+0.5), labels = paste0(c(min(numbers, na.rm=T), max(numbers, na.rm=T)),"days"), padj=1, tck=1  )
    }  
    
    grid();abline(h=1, col="darkgray",lty=1)
    #abline(v=ncoef+0.5, col="darkgray", lty=1)   # "orange"
    abline(v=10*(1:(ncoef%/%10)), col="lightgray", lty=1)
    abline(v= 5*(1:(ncoef%/% 5)), col="lightgray", lty=2)
    
    #  CI's for unadjusted and adjusted RR's:
    if(CI)
      for(imod in 1:length(res)){
        if(sum(c("RR","lci","uci") %in% dimnames(res[[imod]])[[2]])==3 & nrow(res[[imod]])>0 )
          if(any(!is.na( res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("lci","uci")] )))
            matlines( rbind( (1:ncoef+x_deltas[imod]),(1:ncoef+x_deltas[imod]))[,!is.na(res[[imod]]$RR[1:ncoef])],
                      t(res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("lci","uci")]),
                      lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
      }
    
    # RR's:
    for(imod in 1:length(res)){
      if(all(dimnames(res[[imod]])[[2]]!="RR") | nrow(res[[imod]])==0) next
      lines( 1:ncoef+x_deltas[imod],res[[imod]]$RR[1:ncoef], type="o", col=col[imod],lwd=ifelse(imod==1,2,1)); 
      if(imod==1) text( 1:ncoef,res[[1]]$RR, labels=as.character(res[[1]]$i), pos=3, col=col[1], cex= ifelse(ncoef<=50,1,0.7) ) 
      if(any( (cond <- !is.na(res[[imod]]$pval[1:ncoef]) & res[[imod]]$pval[1:ncoef]<=0.05) ))                # check for significant p-values
        points( (1:ncoef+x_deltas[imod])[cond], res[[imod]]$RR[1:ncoef][cond], pch="*",cex=3, col=col[imod]) 
    }
    
    # calendar time adjusted
    if(length(res)>1){
      for(imod in 2:length(res)){
        
        if(all(dimnames(res[[imod]])[[2]]!="RR") | nrow(res[[imod]])==0) next
        
        cond_after_ncoef <- (1:nrow(res[[imod]]))>ncoef
        
        if(ncoef_max>ncoef){
          if(!correct_max_time_adj){
            # CI's for adjusted RR's:
            if(CI & sum(c("RR","lci","uci") %in% dimnames(res[[imod]])[[2]])==3) 
              if(any(!is.na(  res[[imod]][cond_after_ncoef,][!is.na(res[[imod]]$RR[ cond_after_ncoef ]), c("lci","uci")] )))
                matlines( rbind( ncoef +  1:(nrow(res[[imod]])-ncoef), ncoef +  1:(nrow(res[[imod]])-ncoef)  )[,!is.na(res[[imod]]$RR[ cond_after_ncoef ])],
                          t(res[[imod]][cond_after_ncoef,][ !is.na(res[[imod]]$RR[ cond_after_ncoef ])  ,c("lci","uci")]),
                          lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
            # RR's:
            lines( ncoef +  1:(length(res[[imod]]$RR)-ncoef), res[[imod]]$RR[ cond_after_ncoef ], type="o", col=col[imod],cex=0.5)
          }  
          else {
            if(correct_max_time_adj & sum(c("RR","lci","uci") %in% dimnames(res[[imod]])[[2]])==3 ){
              
              cond_time_adj <- substring(res[[imod]][cond_after_ncoef,"all_cat"],1,1)=="["
              
              suppressWarnings(
                mids <- unlist(   lapply( strsplit(  res[[imod]]$all_cat[cond_after_ncoef][ cond_time_adj ], "\\[|;|,|\\]") , 
                                          function(x) mean(as.numeric(x),na.rm=T) ) )
              )
              
              xx_time_adj    <-  ( (mids - min(numbers, na.rm=T) ) * ncoef_max  +( max(numbers,na.rm=T) - mids) * (ncoef+1) ) / ( max(numbers,na.rm=T) -  min(numbers, na.rm=T) )              
              #xx_time_adj    <- ncoef +  1   +   (mids - min(numbers, na.rm=T) ) * ( ncoef_max - ncoef ) / ( max(numbers,na.rm=T) -  min(numbers, na.rm=T) )              
              RR_time_adj    <- res[[imod]]$RR[   cond_after_ncoef][cond_time_adj] 
              CI_time_adj    <- res[[imod]][cond_after_ncoef,][cond_time_adj,c("lci","uci")] 
              pval_time_adj  <- res[[imod]]$pval[ cond_after_ncoef][cond_time_adj] 
              model_time_adj <- res[[imod]]$model[cond_after_ncoef][cond_time_adj] 
              
              if(any(!is.na(model_time_adj))){
                ref_models <- sort(unique( model_time_adj[!is.na(model_time_adj)] ))
                x_ref_deltas <- 0 
                if(length(ref_models)>1) x_ref_deltas <- 0.2* (1/(length(ref_models)-1) * (0:(length(ref_models)-1)) - 0.5)
                
                for(ii in 1:length(ref_models) ){
                  cond_ref_model <- !is.na(model_time_adj) & model_time_adj==ref_models[ii]
                  # CI's:
                  if(CI & any(!is.na(CI_time_adj[cond_ref_model & !is.na(RR_time_adj),])))
                    matlines( rbind( xx_time_adj, xx_time_adj )[, cond_ref_model & !is.na(RR_time_adj) ]+x_ref_deltas[ii], 
                              t(CI_time_adj)[                   , cond_ref_model & !is.na(RR_time_adj) ],
                              lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
                  
                  # RR's:
                  if(any(!is.na(RR_time_adj[ cond_ref_model ])))
                    lines( xx_time_adj[ cond_ref_model ] + x_ref_deltas[ii], RR_time_adj[ cond_ref_model ], type="o", col=col[imod],
                           pch=ifelse( length(unique( model_time_adj[!is.na(model_time_adj)] ))==1, 1, as.character(ii)), cex=0.5)
                  
                  if(any( (cond <- !is.na(pval_time_adj[cond_ref_model]) & pval_time_adj[cond_ref_model]<=0.05) )){  # check for significant p-values
                    if( !correct_max_time_adj)
                      points( ((ncoef+1):nrow(res[[imod]]))[cond_ref_model][cond], RR_time_adj[cond_ref_model][cond], pch="*",cex=3, col=col[imod]) 
                    else
                      points( xx_time_adj[cond_ref_model][cond], RR_time_adj[cond_ref_model][cond], pch="*",cex=3, col=col[imod])
                  }
                }
              }  
            }
          }
        }  
      }  
    } 
    
    legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), 
           text.col=c("gray65","black","blue","red")[text_col_cond], cex=0.6,
           pt.cex=0.6, box.lty=0,bg="transparent", ncol= (ncoef %/% 38) + as.numeric((ncoef %% 38)>0) )
    #legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), cex=0.6, box.lty=0, bg="transparent", ncol= (ncoef %/% 38) + as.numeric((ncoef %% 38)>0) )
    legend("topleft", legend= c("significant","",names(res)), cex=0.7, 
           pt.cex=c( 3,0,rep(1,length(res)) ),
           pch=c( 42,32,rep(1,length(res)) ), lwd=2, lty=c( 0,0,rep(1,length(res)) ), 
           col=c("black","black",col[1:length(res)]), box.lty=0,bg="transparent")  
  }
  
  
  
  ############
  #
  # plot 2 en 3 (if ymax>30): coefficients without coefficients for time adjusting:
  #
  for(i in 2:3) {   # plot 2    and   plot 3
    
    if(i==2){ 
      ymax   <- unlist(lapply(res,function(x) {
        if(is.null(x)) return(NA) 
        if(all(dimnames(x)[[2]]!="RR") | nrow(x)==0 ) return(NA) 
        x_tmp <- x$RR[1:ncoef][x$RR[1:ncoef]<1000]
        if(any(!is.na(x_tmp))) max(x_tmp*1.2,na.rm=T) else NA 
      } ))
      
      if(any(!is.na(ymax))) ymax <- max(ymax, na.rm=T)
      else                  ymax <- 10
      ylim <- c( 0, ymax )
    }  
    
    if(i==3 & max(ylim)<=30) next
    if(i==3 & max(ylim)>30 ) ylim=c(0,20)

    plot( c(0,ncoef), ylim, type="n", main=main, xlab="effect number", ylab=ifelse(i==2, "RR", "RR under 20") )
    plotted <- T
    
    grid();abline(h=1, col="darkgray",lty=1)
    abline(v=length(res[[1]]$RR[1:ncoef])+0.5, col="orange", lty=3)
    
    #  CI's for unadjusted and adjusted RR's:
    if(CI)
      for(imod in 1:length(res)){ 
        if(is.null(res[[imod]])) next
        if(sum(c("RR","lci","uci") %in% dimnames(res[[imod]])[[2]])!=3 | nrow(res[[imod]])==0) next
        if(any(!is.na( t(res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("lci","uci")]) )))
          matlines( rbind( (1:ncoef+x_deltas[imod]),(1:ncoef+x_deltas[imod]))[,!is.na(res[[imod]]$RR[1:ncoef])],
                    t(res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("lci","uci")]),
                    lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
      }
    # RR's:
    for(imod in 1:length(res)){
      if(is.null(res[[imod]])) next
      if(sum(c("RR","lci","uci") %in% dimnames(res[[imod]])[[2]])!=3 | nrow(res[[imod]])==0) next
      lines( 1:ncoef+x_deltas[imod],res[[imod]]$RR[1:ncoef], type="o", col=col[imod],lwd=ifelse(imod==1,2,1)); 
      if(imod==1) text( 1:ncoef,res[[1]]$RR, labels=as.character(res[[1]]$i), pos=3, col=col[1], cex= ifelse(ncoef<=50,1,0.7) ) 
      if(any( (cond <- !is.na(res[[imod]]$pval[1:ncoef]) & res[[imod]]$pval[1:ncoef]<=0.05) ))                # check for significant p-values
        points( (1:ncoef+x_deltas[imod])[cond], res[[imod]]$RR[1:ncoef][cond], pch="*",cex=3, col=col[imod]) 
    }
    
    
    legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), 
           text.col=c("gray65","black","blue","red")[text_col_cond], cex=0.6,
           pt.cex=0.6, box.lty=0,bg="transparent", ncol= (ncoef %/% 38) + as.numeric((ncoef %% 38)>0) )
    #legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), cex=0.6, box.lty=0, bg="transparent", ncol= (ncoef %/% 38) + as.numeric((ncoef %% 38)>0) )
    legend("topleft", legend= c("significant","",names(res)), cex=0.7, 
           pt.cex=c( 3,0,rep(1,length(res)) ),
           pch=c( 42,32,rep(1,length(res)) ), lwd=2, lty=c( 0,0,rep(1,length(res)) ), 
           col=c("black","black",col[1:length(res)]), box.lty=0,bg="transparent")
  }
  invisible(plotted)
}  # end of plot_res

###############################
#
plot_res2 <- function(res, tit="", col=col_list ){
  par(mfrow=c(2,1))
  ymax <- max(unlist(lapply(res,function(x)max(x$RR[x$RR<500],na.rm=T))))
  
  
  plot(res[[2]]$RR, type="o", col="blue",ylim=c(0,  ymax ), 
       main=tit,xlab="coefficients number",ylab="RR"); 
  grid();abline(h=1, col="gray",lty=3)
  abline(v=length(res[[1]]$RR)+0.5, col="orange", lty=3)
  lines(1:length(res[[1]]$RR)+0.15,res[[1]]$RR, type="o", col=col[1],lwd=2); 
  for(i in 2:length(res)){
    time_adj <- strsplit( substring(dimnames(res[[i]])[[1]],1,(length(dimnames(res[[i]])[[1]])-1)), ":")
    lapply(time_adj,function(x)x)
    
    lines( res[[i]]$RR, type="o", col=col[i])
  }  
  
  text_col_cond <- 1 + as.numeric(!is.na(res[[1]]$RR))
  text_col_cond <- text_col_cond + as.numeric(!is.na(res[[1]]$RR) & !is.na(res[[1]]$pval<=0.1))
  text_col_cond <- text_col_cond + as.numeric(!is.na(res[[1]]$RR) & !is.na(res[[1]]$pval<=0.05))
  legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), 
         text.col=c("gray65","black","blue","red")[text_col_cond], 
         cex=0.7, box.lty=0,bg="transparent")
  #legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), cex=0.7, box.lty=0,bg="transparent")
  legend("topleft",legend= names(res), cex=0.7, pch=1, lwd=1, col=col[1:length(res)], box.lty=0,bg="transparent")
}



###############################
#
#  functions to add results to report list and model list and save them.  ==>  delete ???
#
add_to_report_list <- function(x, name, list=report_list, add=T){ 
  if(!add) list <- list()
  else if(!is.list(list)) stop("'list' must be a list") 
  list <- c( list, list(x) )
  if(!missing(name)) names(list)[length(list)] <- name
  list
}
add_to_models_list <- function(x, name, list=models_list, add=T){ 
  if(!add) list <- list()
  else if(!is.list(list)) stop("'list' must be a list") 
  list <- c( list, list(x) )
  if(!missing(name)) names(list)[length(list)] <- name
  list
}

###########
#  save report_list and model_list  ==>  ?? delete
#
save_results <- function(name, report_list, models_list, sep="" ){
  
  name <- gsub(":","",name)
  
  begin_report <- paste0(dap, sep, name, "_report")
  begin_report <- trimws( begin_report,"left","[ _:;,.]" )
  
  
  # copy report [and models] lists to lists with other names:
  assign(begin_report,  report_list )
  
  # save report [and models] list as .RData:
  save(list=begin_report, file = paste0(sdr, begin_report,".RData" ))
  
  if( !missing(models_list))
    if(length(models_list)>0 ){
      begin_models <- paste0(dap, sep, name, "_models")
      assign(begin_models,  models_list )
      save(list=begin_models, file = paste0(sdr_models, begin_models,".RData" ))
    }
  
  # print tables from the 'report' list in .txt file:
  sink(paste0(sdr, begin_report, ".txt" ))
  
  old_width = options (width=200, max.print=99999 )
  cat(paste0("\n\nNumber of rows in the dataset = ", nrow(scri_input),".\n\n\n"))
  
  for(i in 1:length(report_list)){
    cat(paste0("\n\n",names(report_list[i]),":"))
    for(j in 1:length(report_list[[i]])){
      cat(paste0("\n\n",names(report_list[[i]][j]),":\n\n"))
      print(lapply(report_list[[i]][[j]], format, justify="left", digits=3))
    }
  }  
  options(old_width)
  
  sink()
  
}

#
#####################


#################################################################
#################################################################
#################################################################
#
scri <- function(vax_def,
                 extra_name = "", 
                 formula    = "", 
                 time_seq   = c(), time_seq_ref = "with_events", # c("with_events","without_events"),
                 event_time, event, event_date, event_info,
                 id,
                 rws ,
                 ref,
                 combine_vars = c(),
                 start_obs, end_obs,
                 lab_orders   = NA,
                 data,  inclusion="", inclusion_name="",
                 strata_var   = "", strata_value = NA, use_all_events = T,
                 data_event   = NA,
                 data_source = "",
                 censored_vars = c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                 rw_observed_percentage = 0,   # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                 #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                 event_in_rw = T,                # if event in rw ==> this rw should not be deleted even if not completely observed
                 lplots      = T,
                 lplot_hist  = T, add_histplot=F,
                 leventplot  = T, max_n_points = NA,  ngrid = c(500,300),  only_add_plot_to_file=F, eventplot_file_separate = F, path_cohort="", image_performance = F, warn_image_plots =-1, 
                 lplot       = T, CI = T,
                 lforest     = T, forest_nrows = 50 ,forest_cex = 0.5, forest_cex_head = 0.5,
                 nvax,
                 delete_coef_no_events      = T,
                 delete_rows_start_with_no  = T,
                 delete_no_ref_cols         = T,
                 delete_median_cols         = T,
                 lparal          = T,  n_cores = NA,          # paral_vars = c(), parallel programming (using more cores)
                 lprint          = T, 
                 save_data       = F,
                 width           = 14, 
                 extra_parameters,
                 add_to          = "",
                 add_to_itself   = F,
                 ndigits         = 2,
                 sdr_tabs, sdr_models, 
                 cut_points_name = "",
                 col             = c("red", "green3", "orange",  "deepskyblue", "magenta2", "cyan2", "chocolate1", "gray" ),
                 performance     = F, 
                 ...
){  
  
  if(missing(data)) stop("Dataset 'data' is missing.")
  
  start_sys_time <- Sys.time()
  
  if(!missing(vax_def))
    if(class(vax_def)=="scri_parameters"){
      # vax variables:
      if(missing(nvax))          nvax          <- vax_def$nvax    # delete?
      
      # variable names:
      if(missing(id))            id            <- vax_def$data_parameters$id
      if(missing(start_obs))     start_obs     <- vax_def$data_parameters$start_obs
      if(missing(end_obs))       end_obs       <- vax_def$data_parameters$end_obs
      if(missing(censored_vars)) 
        if(any(names(vax_def$data_parameters)=="censored_vars")) 
          censored_vars <- vax_def$data_parameters$censored_vars
      
      # model parameters:
      if(missing(rws))              rws             <- vax_def$rws$rws_def
      if(missing(cut_points_name))  cut_points_name <- vax_def$rws$cut_points_name
      
      # lab order:
      if(missing(lab_orders))    lab_orders    <- vax_def$lab_orders
      
    }
  
  if(!missing(lplots))
    if(!lplots){
      if( missing(leventplot)) leventplot <- F  
      if( missing(lplot     )) lplot      <- F  
      if( missing(lforest   )) lforest    <- F  
    }
  
  
  if(!missing(extra_parameters)){
    
    if( missing(extra_name  ) & "extra_name"   %in% names(extra_parameters) ) extra_name  <- extra_parameters[["extra_name"  ]] 
    if( missing(data_source ) & "data_source"  %in% names(extra_parameters) ) data_source <- extra_parameters[["data_source" ]] 
    
    if( missing(inclusion      ) & "inclusion"       %in% names(extra_parameters) ) inclusion      <- extra_parameters[["inclusion"      ]] 
    if( missing(inclusion_name ) & "inclusion_name"  %in% names(extra_parameters) ) inclusion_name <- extra_parameters[["inclusion_name" ]] 
    
    
    
    
    #    if(missing(lplots)){  # ???
    #      
    #      if( missing(lplots                 ) & "lplots"                  %in% names(extra_parameters) ){
    #        lplots                  <- extra_parameters[["lplots"            ]] 
    
    if( missing(lplot_hist             ) & "lplot_hist"              %in% names(extra_parameters) ) lplot_hist              <- extra_parameters[["lplot_hist"              ]]   
    if( missing(add_histplot           ) & "add_histplot"            %in% names(extra_parameters) ) add_histplot            <- extra_parameters[["add_histplot"            ]]   
    if( missing(leventplot             ) & "leventplot"              %in% names(extra_parameters) ) leventplot              <- extra_parameters[["leventplot"              ]]   
    if( missing(only_add_plot_to_file )  & "only_add_plot_to_file"   %in% names(extra_parameters) ) only_add_plot_to_file   <- extra_parameters[["only_add_plot_to_file"   ]]   
    if( missing(max_n_points           ) & "max_n_points"            %in% names(extra_parameters) ) max_n_points            <- extra_parameters[["max_n_points"            ]]   
    if( missing(eventplot_file_separate) & "eventplot_file_separate" %in% names(extra_parameters) ) eventplot_file_separate <- extra_parameters[["eventplot_file_separate" ]] 
    
    if( missing(lplot           ) & "lplot"             %in% names(extra_parameters) ) lplot            <- extra_parameters[["lplot"        ]]   
    if( missing(CI              ) & "CI"                %in% names(extra_parameters) ) CI               <- extra_parameters[["CI"           ]]
    if( missing(CI              ) & "CI_draw"           %in% names(extra_parameters) ) CI               <- extra_parameters[["CI_draw"      ]]
    
    if( missing(lforest         ) & "lforest"         %in% names(extra_parameters) ) lforest          <- extra_parameters[["lforest"        ]]   
    if( missing(forest_nrows    ) & "forest_nrows"    %in% names(extra_parameters) ) forest_nrows     <- extra_parameters[["forest_nrows"   ]]   
    if( missing(forest_cex      ) & "forest_cex"      %in% names(extra_parameters) ) forest_cex       <- extra_parameters[["forest_cex"     ]]   
    if( missing(forest_cex_head ) & "forest_cex_head" %in% names(extra_parameters) ) forest_cex_head  <- extra_parameters[["forest_cex_head"]]   
    
    #      }
    #    }
    #    else {
    #      if( missing(leventplot)) leventplot <- lplots  
    #      if( missing(lplot     )) lplot      <- lplots  
    #      if( missing(lforest   )) lforest    <- lplots  
    #    }
    
    if( missing(delete_coef_no_events    ) & "delete_coef_no_events"     %in% names(extra_parameters) ) delete_coef_no_events     <- extra_parameters[[ "delete_coef_no_events"     ]]   
    if( missing(delete_coef_no_events    ) & "delete_coef_no_events"     %in% names(extra_parameters) ) delete_coef_no_events     <- extra_parameters[[ "delete_coef_no_events"     ]]   
    if( missing(delete_rows_start_with_no) & "delete_rows_start_with_no" %in% names(extra_parameters) ) delete_rows_start_with_no <- extra_parameters[[ "delete_rows_start_with_no" ]]   
    if( missing(delete_no_ref_cols       ) & "delete_no_ref_cols"        %in% names(extra_parameters) ) delete_no_ref_cols        <- extra_parameters[[ "delete_no_ref_cols"        ]]   
    if( missing(delete_median_cols       ) & "delete_median_cols"        %in% names(extra_parameters) ) delete_median_cols        <- extra_parameters[[ "delete_median_cols"        ]]   
    
    
    if( missing(lparal      ) & "lparal"      %in% names(extra_parameters) ) lparal       <- extra_parameters[["lparal"      ]]   
    if( missing(n_cores     ) & "n_cores"     %in% names(extra_parameters) ) n_cores      <- extra_parameters[["n_cores"     ]] 
    
    if( missing(lprint      ) & "lprint"      %in% names(extra_parameters) ) lprint       <- extra_parameters[["lprint"      ]]   
    if( missing(save_data   ) & "save_data"   %in% names(extra_parameters) ) save_data    <- extra_parameters[["save_data"   ]]   
    if( missing(width       ) & "width"       %in% names(extra_parameters) ) width        <- extra_parameters[["width"       ]]
    
    if( missing(time_seq     ) & "time_seq"     %in% names(extra_parameters) ) time_seq      <- extra_parameters[["time_seq"     ]] 
    if( missing(time_seq_ref ) & "time_seq_ref" %in% names(extra_parameters) ) time_seq_ref  <- extra_parameters[["time_seq_ref" ]] 
    
    if( missing(sdr_tabs    ) & "sdr_tabs"    %in% names(extra_parameters) ) sdr_tabs     <- extra_parameters[["sdr_tabs"    ]]   
    if( missing(sdr_models  ) & "sdr_models"  %in% names(extra_parameters) ) sdr_models   <- extra_parameters[["sdr_models"  ]] 
    if( missing(path_cohort ) & "path_cohort" %in% names(extra_parameters) ) path_cohort  <- extra_parameters[["path_cohort" ]]   
    
    if( missing(performance       ) & "performance"       %in% names(extra_parameters) ) performance       <- extra_parameters[["performance"       ]]    
    if( missing(image_performance ) & "image_performance" %in% names(extra_parameters) ) image_performance <- extra_parameters[["image_performance" ]]    
    if( missing(warn_image_plots  ) & "warn_image_plots"  %in% names(extra_parameters) ) warn_image_plots  <- extra_parameters[["warn_image_plots"  ]]    
    
  }
  
  
  
  
  
  
  
  
  
  if(!missing(event_info)){
    
    if(missing(event_time)) event_time <- event_info[["event_time"]]
    if(missing(event_date)) event_date <- event_info[["event_date"]]
    if(missing(event     )) event      <- event_info[["event"     ]]
    
    if(any(!is.na(data[,event_date])) & all( is.na(data[,event_time]))) data[,event_time] <- as.numeric(difftime(data[,event_date],as.Date("2020-08-31"),units="days")) 
    if(all( is.na(data[,event_date])) & any(!is.na(data[,event_time]))) data[,event_date] <- as.Date("2020-08-31") + data[,event_time]
  }
  data[,event] <- as.numeric(!is.na(data[,event_date]))
  
  ldata <- T
  if(!any(data[,event]==1)) ldata <- F
  
  #######################################################
  #
  #             output structure, content of output files, directories:
  #
  ################
  # directories:
  #
  #     directory 1: g_export / scri / { data_source }_{ event }  /   ?  /{output_name}.txt; {output_name}.RData; {output_name}_{stratum ?or all strata? }.pdf
  #
  #     directory 2: g_?local? / scri / { event }  / ?    /{output_name}_models.RData
  #
  ###############
  # 'output_name'
  #
  #     output_name <- { data_source }_{ event }_[ extra_name ]_[ nothing or covid_selection ]_{ no_split; brand  via vax_dep 'before' }_{ dist  via vax_dep 'after'}_{ nothing or stratum }
  #
  ###############
  # content of {output_name}.RData, {output_name}.txt, {output_name}_models.RData files; {output_name}.pdf
  #   
  #    without strata:      |    with strata:                                                                                                                                           
  #    _____________________|________________________________________________________________________                                                                           
  #                         |                                                                                                                                                                              
  #                         |                                                                                                                                                         
  #     - cut_points_28d    |      -- cut_points_28d                                                                
  #       -- $no_adj        |         -- $no_adj                    ( <- if only stratum events: use_all_events=F )                                                                                              
  #       -- $adj_30ds_1    |         -- $adj_30ds_1                          --  ||  --                                                                                                    
  #                         |         -- $adj_20ds_10                         --  ||  --                                                                                
  #                         |         -- ...                                                                                                                          
  #       -- $adj_20ds_10   |         -- $no_adj_all_events         ( <- all events from are used: use_all_events=T )                                                                                                                   
  #                         |         -- $adj_30ds_1_all_events     (    they are place in on separate category for risk windows )                                                                                      
  #                         |         -- $adj_20ds_10_all_events    (    but they are used for calendar time intervals           )                                                                                                        
  #       -- ...            |         -- ...                                  --  ||  --                                                                                                                      
  #                         |                                                                                                                                            
  #     - cut_points_7d     |      -- cut_points_7d                                                                                                                          
  #       -- ...                      -- ...                                                                                                                             
  #
  
  output_name <- data_source
  if(inclusion_name!="")                          output_name <- c( output_name, inclusion_name                                                       )
  if(any(ls()=="event"))                          output_name <- c( output_name, event                                                                )
  if(extra_name != "")                            output_name <- c( output_name, extra_name                                                           )
  if(any(ls()=="covid_selection"))                output_name <- c( output_name, "covid_selection"                                                    )
  if(any(names(vax_def$rws)=="vax_dep")){                                                                  
    if(any(names(vax_def$rws$vax_dep)=="before")) output_name <- c( output_name, paste0( gsub("_","",vax_def$rws$vax_dep["before"]), collapse = "_")  )
    if(any(names(vax_def$rws$vax_dep)=="after" )) output_name <- c( output_name, paste0( gsub("_","",vax_def$rws$vax_dep["after"] ), collapse = "_")  )         
  } else                                          output_name <- c( output_name, "nosplit"                                                            )
  if(!is.na(strata_value))                        output_name <- c( output_name, strata_value                                                         )
  
  output_name <- paste0(output_name, collapse="_")
  output_name <- trimws(output_name, "left", "_")
  
  cat(output_name)
  
  output_file_name <- gsub("[:;,.-]","_",output_name)
  
  
  id        <- vax_def$data_parameters$id
  vax_time  <- vax_def$data_parameters$vax_time
  vax_date  <- vax_def$data_parameters$vax_date
  vax_name  <- vax_def$data_parameters$vax_name
  vax_dep   <- vax_def$rws$vax_dep
  
  plotted <- F
 
  # create leventplot (image) plots:
  if(leventplot & ldata){ 
    
    if(any( !( names(list(...)) %in% names(formals(pdf)) ) )) 
      stop(paste0("Argument[s] '", paste(names(list(...))[ !( names(list(...)) %in% names(formals(pdf)) )], collapse="', '" ),
                  "' is not an argumnent of function 'pdf' from library 'qpdf'."))
    
    gc()
    pdf(file=paste0(sdr_tabs,"eventplot_tmp.pdf"), width=width,  ... )
    
    image_plots( vax_def=vax_def, file=ifelse(path_cohort=="","",paste0(path_cohort,output_file_name)), data=data, event_info=event_info,  strata_var=strata_var, strata_value=strata_value, 
                 tit=output_name, max_n_points=max_n_points, warn=warn_image_plots, performance=image_performance )
    plotted <- T
    dev.off()
    
    if(only_add_plot_to_file){
      
      if(plotted & file.exist( paste0(sdr_tabs,output_file_name,".pdf") )) {         
        files_to_copy <- paste0(sdr_tabs,c(output_file_name,"eventplot_tmp"),".pdf")
        for(iname in files_to_copy) if(!file.exists(iname)) files_to_copy <- files_to_copy[files_to_copy!=files_to_copy]
        if(length(files_to_copy)>0) qpdf::pdf_combine( files_to_copy , paste0(sdr_tabs,output_file_name,".pdf")  )  
      }
      
      if(file.exists(paste0(sdr_tabs,"eventplot_tmp.pdf" )))  suppressWarnings( file.remove(paste0(sdr_tabs, "eventplot_tmp.pdf"   )) )
      if(!is.null(dev.list())) dev.off()
      
      if(performance) cat(paste0(": duration = ",format(difftime(Sys.time(),start_sys_time))," (till ",Sys.time(),")\n"))  
      else            cat("\n")  
      
      return(NULL)
    }
    
  }  # end if leventplot
 
  if(ldata){
	
    if(inclusion!="") {
      if(inclusion %in% names(data)) data <- data[data[,inclusion],]
      else stop(paste0("inclusion variable '",inclusion,"' is not in the dataset."))
    }
    # stratum:
    # if only for subset (i.e., use_all_events==F) ==> delete all other rows
    if(strata_var!="" & !use_all_events & ldata) {
      if(mode(data[,strata_var])=="logical") data <- data[ data[,strata_var] & !is.na(data[,strata_var]), ]
      else {
        if(!is.na(strata_value)) data <- data[ data[,strata_var]== strata_value & !is.na(data[,strata_var]), ]
        else stop(paste0("'strata_var' should be logical or a variable with value 'strata_value'. (Now: 'strata_var'= '",strata_var,"'; 'strata_value'=",ifelse(is.na(strata_value),"NA",paste0("'",strata_value,"'")),")"))
      }
      strata_var <- ""
    }
    
    data <- data[data[,event]==1,] 
    if(nrow(data)==0) ldata <- F
  }

  # time_seq:
  # add an empty time sequence for a model without time category:
  if(length(time_seq)>=1){
    time_seq_ref <- c("", rep( time_seq_ref,  each=length(time_seq)) )
    time_seq     <-       rep( time_seq, (length(time_seq_ref)-1)/length(time_seq) )
  }
  time_seq <- c( no_adj = list(c()), time_seq )
																									  
  
  
  if(any(ls()=="res"      )) rm(res)
  if(any(ls()=="res_paral")) rm(res_paral)
  
  res <- vector("list", length=length(time_seq))
  if(length(res)==1) names(res) <- names(time_seq)
  else names(res) <- paste0(names(time_seq),"_",time_seq_ref)
  names(res)[substring(names(res),nchar(names(res)))=="_"] <- substring(names(res),1,nchar(names(res))-1)[substring(names(res),nchar(names(res)))=="_"]
  
  if(strata_var!="" & use_all_events) names(res) <- paste0(names(res), "_all_events")
  
  #############################
  #
  # parallel use of cores:
  #
  if(length(time_seq)>1 & lparal & ldata){
    
    library(parallel)
    
    if(is.na(n_cores)) n_cores <- detectCores() - 2 
    n_cores <- min( length(time_seq), n_cores, na.rm=T )
    
    cl      <- makeCluster( n_cores )    # outfile = paste0(sdr,"log_parallel.txt")
    
    clusterEvalQ( cl, { library("survival") } )
    clusterExport(cl, c("scri_fit", "refresh_event_variable", "create_rws", "split_intervals", 
                        "summary_tab", "factor_ref", "combine_vars_func" ), 
                  envir=environment()  )
    
    for(ipar in names(as.list(match.call()))[-1])
      clusterExport(cl, ipar, envir=environment() )        
    #if(!missing(paral_vars)) clusterExport(cl, paral_vars )
    
    res_paral <- parLapply(cl,
                           time_seq,  # list with different sets of intervals
                           function(i_time_seq) 
                             scri_fit(formula =  paste(formula, ifelse(is.null(i_time_seq),"", paste( " +  cal_time_cat") ) ),
                                      #formula = ifelse(formula=="","", paste( formula, " +  cal_time_cat") ),
                                      vax_def=vax_def,
                                      event_time = event_time, event = event, id=id,
                                      rws = rws, 
                                      ref=ref,
                                      time_seq = i_time_seq,  split_seq_name = "cal_time_cat", #time_seq_ref="most events",time_seq_ref=time_seq_ref[i_time_seq],
                                      combine_vars = combine_vars,
                                      start_obs = start_obs, end_obs = end_obs,
                                      lab_orders = lab_orders,
                                      data = data, strata_var=strata_var, strata_value=strata_value, #use_all_events=use_all_events,
                                      data_event  = data_event,
                                      data_source = data_source,
                                      nvax       = nvax,
                                      delete_coef_no_events      = delete_coef_no_events,    
                                      delete_rows_start_with_no  = delete_rows_start_with_no,
                                      delete_no_ref_cols         = delete_no_ref_cols,       
                                      delete_median_cols         = delete_median_cols,       
                                      rw_observed_percentage     = rw_observed_percentage,   
                                      censored_vars = censored_vars,           
                                      event_in_rw = event_in_rw,               
                                      lprint=lprint,
                                      lplot_hist = lplot_hist, width=width, sdr_tabs = sdr_tabs, add_histplot=add_histplot, 
                                      save_data = save_data, ...)
    )
    
    stopCluster(cl)
    
    names_overlap <- match(names(res),names(res_paral))
    res[ names(res_paral)[ names_overlap[!is.na(names_overlap)] ] ] <- res_paral[!is.na(match(names(res_paral),names(res)))]
    
  }  # end of parallel

  if( (!lparal | length(time_seq)==1) & ldata )
    for(i in 1:length(res)) 
      res[[i]] <-  scri_fit(formula =  paste(formula, ifelse(is.null(time_seq[[i]]),"", paste( " +  cal_time_cat") ) ),
                            vax_def=vax_def,
                            event_time = event_time, event = event, id=id,
                            rws = rws, 
                            ref=ref,
                            time_seq = time_seq[[i]] , time_seq_ref=time_seq_ref[i], split_seq_name = "cal_time_cat", #time_seq_ref="most events", 
                            combine_vars = combine_vars,
                            start_obs = start_obs, end_obs = end_obs,
                            lab_orders = lab_orders,
                            data = data, strata_var=strata_var, strata_value=strata_value, #use_all_events=use_all_events,
                            data_event   = data_event,
                            data_source = data_source,
                            nvax       = nvax,
                            delete_coef_no_events      = delete_coef_no_events,    
                            delete_rows_start_with_no  = delete_rows_start_with_no,
                            delete_no_ref_cols         = delete_no_ref_cols,       
                            delete_median_cols         = delete_median_cols,       
                            rw_observed_percentage     = rw_observed_percentage,     
                            censored_vars = censored_vars,           
                            event_in_rw = event_in_rw,               
                            lprint=lprint,
                            lplot_hist = lplot_hist, width=width, sdr_tabs = sdr_tabs,  add_histplot=add_histplot,
                            save_data = save_data, ...)

  if(ldata){
    
    # for output file with 'tabs'  
    tabs             <- lapply(res,function(x) x$tab )
    attributes(tabs) <- c( attributes(tabs), lapply(res,function(x)x[names(x)==c("tab_full")]) )
    class(tabs)      <- "scri_tabs"
    attributes(tabs) <- c( attributes(tabs), 
                           name      = output_name,
                           vax_def = list(vax_def),                  
                           event   = list( c(event=event, event_time=event_time, event_date=event_date) ) )
  } 
  else tabs <- res
  
  res <- c (res, 
            name    = output_name,
            vax_def = list(vax_def),                  
            event   = list( c( event=event, event_time=event_time, event_date=event_date))   )
  
  res  <- list(res) ;  if(cut_points_name!="") names( res) <- cut_points_name
  tabs <- list(tabs);  if(cut_points_name!="") names(tabs) <- cut_points_name
  
  if(add_to_itself) add_to <- output_name
  
  if(add_to!="") {
    
    load(paste0(sdr_tabs,   gsub("[:;,.-]","_",add_to),        ".RData"))
    load(paste0(sdr_models, gsub("[:;,.-]","_",add_to), "_models.RData"))
    
    tabs <- c(        get(add_to)           , tabs )
    res  <- c( get(paste0(add_to,"_models")), res  ) 
    
  }
  
  
  ##################  
  # save the results in .RData and .txt files:  
  
  assign(       output_name           , tabs)
  assign(paste0(output_name,"_models"), res )
  
  save( list =        output_name           , file = paste0(sdr_tabs  , output_file_name,       ".RData" ))
  save( list = paste0(output_name,"_models"), file = paste0(sdr_models, output_file_name,"_models.RData" ))
  
  # print tables from the 'report' list in .txt file:
  sink(paste0(sdr_tabs, output_file_name, ".txt" ))
  {
    old_width = options (width=300, max.print=99999 )
    #cat(paste0("\n\nNumber of ids in the dataset = ", length(data_vax[,?id]),".\n\n\n"))   #???
    #cat(paste0("\n\nNumber of rows in the dataset = ", nrow(scri_input),".\n\n\n"))   #???
    
    cat(paste("\n\n Name:\t",output_name,"\n"))
    
    if(length(tabs)==0) cat("no data\n\n")
    else if(all(sapply(tabs,function(x)is.null(x[[1]])))) cat("no data\n\n")
    
    for(i1 in 1:length(tabs)){
      
      if(length(tabs[[i1]])==0){ cat("no data\n\n"); next}
      
      txt_i1 <- paste0("\n\n",names(tabs)[i1])
      if(length(dim(tabs[[i1]]))==2){
        cat(paste(txt_i1,":\n\n"))
        print(format(tabs[[i1]],              justify="left", digits=ndigits))
      }
      
      for(i2 in 1:length(tabs[[i1]])){
        txt_i2 <- paste0(txt_i1," $ ",names(tabs[[i1]])[i2])
        
        if(length(tabs[[i1]][[i2]])==0){ cat("no data\n\n"); next}
        
        if(length(dim(tabs[[i1]][[i2]]))==2){
          cat(paste(txt_i2,":\n\n"))
          print(format(tabs[[i1]][[i2]],        justify="left", digits=ndigits, nsmall=ndigits))
        }
        if(length(tabs[[i1]][[i2]])>0)
          for(i3 in 1:length(tabs[[i1]][[i2]])){
            txt_i3 <- paste0(txt_i2," $ ",names(tabs[[i1]][[i2]])[i3])
            if(length(dim(tabs[[i1]][[i2]][[i3]]))==2){
              cat(paste(txt_i3,":\n\n"))
              print(format(tabs[[i1]][[i2]][[i3]], justify="left", digits=ndigits, nsmall=ndigits))
            }
          } # for i3
      } # for i2
    } # for i1  
    
    
    cat("\n\n\nSpecified parameters:\n\n")
    
    for(i1 in 1:length(tabs)){
      cat(paste(" ", names(tabs)[i1],":\n\n"))
      cat("   Event info:\t")
      cat( paste0( "c( ", paste0( paste0( names(attributes(tabs[[i1]])$event), ' = "', attributes(tabs[[i1]])$event, '"' ), collapse=', ' ), " )" ))
      
      cat("\n\n\n   Vaccine info:\n")
      print(attributes(tabs[[i1]])$vax_def)
      cat("\n\n")
    }
    
    options(old_width)
  }  
  sink()
  
  # print tables from the 'report' list in .txt file:
  sink(paste0(sdr_models, output_file_name, "_models.txt" ))
  {
    old_width = options (width=300, max.print=99999 )
    #cat(paste0("\n\nNumber of ids in the dataset = ", length(data_vax[,?id]),".\n\n\n"))   #???
    #cat(paste0("\n\nNumber of rows in the dataset = ", nrow(scri_input),".\n\n\n"))   #???
    
    cat(paste("\n\n Name:\t",output_name,"\n"))
    
    if(length(res)==0) cat("no data\n\n")
    else if(all(sapply(res,function(x)is.null(x[[1]])))) cat("no data\n\n")
    
    print(res)
 
    options(old_width)
  }  
  sink()
  
  ret <- list( tabs   = tabs,
               models = res
  )
  class(ret) <- "scri_output"
  
  
  
  if( (lplot | lforest | leventplot)){ 
    if(lplot | lforest ){
      pdf(file=paste0(sdr_tabs,"tmp.pdf"), width=width,  ...)
      if(!is.null(tabs[[1]])){
        # create forest plot:
        if(lforest) {
          if(missing(forest_nrows)) 
            forest_nrows <-  pmax(1, 58 - 4*length(tabs[[i1]]) + pmin(1,length(tabs[[i1]])%/%11)*(-2*10 + 2*length(tabs[[i1]]) ) + pmin(1,length(tabs[[i1]])%/%16)*( -15 + length(tabs[[i1]]) ) ) 
          par(mfrow=c(1,2))
          for(i1 in 1:length(tabs)){
            
            forest_res <- try( forest_plots_tab( tabs[[i1]], nrows_forest_plot=forest_nrows,cex=forest_cex, cex_head=forest_cex_head, ltable=F, col=col ) )
            
            if(class(forest_res)[[1]] != "try-error")  plotted <- plotted | forest_res
          }
          
          
        }
        # create plots with coefficients:
        if(lplot){
          par(mfrow=c(1,1))
          for(i1 in 1:length(tabs))
            plotted <- plotted | plot_res(tabs[[i1]], main=paste( event, formula,"; \n",output_name,"; ", names(tabs)[i1]), col=col, CI=CI) 
        }
      }
      dev.off()  # end pdf
    }
    
    if(!plotted) {         
      pdf(file=paste0(sdr_tabs,output_file_name,".pdf"), width=width,  ...)
      plot(1,1,axes=F,type="n",xlab="",ylab="")
      text(1,1,"no data", cex=1.5) 
      dev.off()
    } 
    else { 
      files_to_copy <- c()
      if(lplot | lforest                      ) files_to_copy <- c( files_to_copy, paste0(sdr_tabs,"tmp.pdf"          ) ) 
      if(lplot_hist & add_histplot            ) files_to_copy <- c( files_to_copy, paste0(sdr_tabs,"histplots_tmp1.pdf") )
      if(lplot_hist                           ) files_to_copy <- c( files_to_copy, paste0(sdr_tabs,"histplots_tmp.pdf") )
      if(leventplot & !eventplot_file_separate) files_to_copy <- c( files_to_copy, paste0(sdr_tabs,"eventplot_tmp.pdf") )
     
      for(iname in files_to_copy) if(!file.exists(iname)) files_to_copy <- files_to_copy[files_to_copy!=files_to_copy]
      if(length(files_to_copy)>0) qpdf::pdf_combine( files_to_copy , paste0(sdr_tabs,output_file_name,".pdf")  )  
    }

    if(file.exists(paste0(sdr_tabs,"histplots_tmp1.pdf")))  suppressWarnings( file.remove(paste0(sdr_tabs,"histplots_tmp1.pdf"   )) )
    if(file.exists(paste0(sdr_tabs,"histplots_tmp.pdf" )))  suppressWarnings( file.remove(paste0(sdr_tabs,"histplots_tmp.pdf"    )) )
    if(file.exists(paste0(sdr_tabs,"eventplot_tmp.pdf" )))  suppressWarnings( file.remove(paste0(sdr_tabs, "eventplot_tmp.pdf"   )) )
    if(file.exists(paste0(sdr_tabs,          "tmp.pdf" )))  suppressWarnings( file.remove(paste0(sdr_tabs,           "tmp.pdf"   )) )
    if(!is.null(dev.list())) dev.off()
  }  # end plots
  
  if(performance) cat(paste0(": duration = ",format(difftime(Sys.time(),start_sys_time))," (till ",Sys.time(),")\n"))  
  else            cat("\n")  
  
  ret
  
}  # the end of function 'scri'


#
#################################################################
#################################################################
#################################################################






#################################################################
#################################################################
#################################################################
#
#    functions for images and 3D plots:
# 
################


image_plots <- function(vax_def, file="", event_info, vax1_time="vax_days_v1", event_time, death_time="death_days", data, date_axes=T, small_range=c(-90,90),
                        tit="", strata_var="", strata_value=NA, max_n_points=NA, ngrid=c(500,300), cex=0.5, warn=-1, performance=F){  
  
  cat(paste(" [ image plot:",(ttt<-Sys.time())," " ))
  id        <- vax_def$data_parameters$id
  vax_name  <- vax_def$data_parameters$vax_name
  vax_date  <- vax_def$data_parameters$vax_date
  vax_days  <- vax_def$data_parameters$vax_time
  vax_time  <- vax_def$data_parameters$vax_time
  start_obs <- vax_def$data_parameters$start_obs
  stop_obs  <- vax_def$data_parameters$end_obs
  
  event      <- event_info$event           
  if(missing(event_time)) event_time <- event_info$event_time   
  event_date <- event_info$event_date  
  
  if(strata_var!="") {
    if(!is.na(strata_value)){
      if(mode(data[,strata_var])=="logical") data <- data[data[,strata_var],]
      else  data <- data[data[,strata_var]==strata_value & !is.na(data[,strata_var]), ]
    }
    else data <- data[data[!is.na(data[,strata_var]),strata_var], ]
  }
  if(nrow(data)==0) {
    plot(1,1,axes=F,type="n",xlab="",ylab="")
    text(1,1,"no data", cex=1.5) 
    return()
  }
  
  # define unique values of vaccine-dependent variables (e.g., brand, distance between doses)
  vaxdep_values <- c()
  if(!is.null(vax_def$rws$vax_dep)){  
    vax_dep   <- vax_def$rws$vax_dep
    vaxdep_vars <- c( vax_dep[names(vax_dep)=="before"],vax_dep[names(vax_dep)=="after"] )
    data$vaxdep_all <- data[,vaxdep_vars[1]]
    if(length(vaxdep_vars)>1)
      for(ivaxdep_var in vaxdep_vars[-1] )
        data$vaxdep_all <-  paste0(data$vaxdep_all," & ",data[,ivaxdep_var])
    vaxdep_values <- data$vaxdep_all[!duplicated(data$vaxdep_all)]
    vaxdep_values <- levels(factor_ref(vaxdep_values, lab_orders = vax_def$lab_orders))
  }
  if(length(vaxdep_values)==0) vaxdep_values <- "all"
  
  
  # time from vax for  event or death
  data$time_event_min_vax   <- data[,event_time] - data[,vax_time ]
  data$time_event_min_vax1  <- data[,event_time] - data[,vax1_time]
  data$time_death_min_vax   <- data[,death_time] - data[,vax_time ]
  data$time_death_min_vax1  <- data[,death_time] - data[,vax1_time]
  
  
  # create dataset with death rows (vaxed and unvaxed):
  data_with_deaths <- as.data.frame(data[!is.na(data[, death_time]), ])
  #if(any(class(data_with_deaths[,death_time]) %in% c("POSIXct","POSIXt","Date"))) data_with_deaths$death_date <- data_with_deaths[,death_time]  
  if(any(class(data_with_deaths[,death_time]) %in% c("numeric","integer"))) data_with_deaths$death_days <- data_with_deaths[,death_time]  
  data_with_deaths$time_death_min_vax  <- data_with_deaths[,death_time] - data_with_deaths[,vax_time ]
  
  total_deaths <- sum(!duplicated(data_with_deaths[,id]))
  
  # the number of vaccinated and unvaccinated id's:
  n_unvaccinated <- sum(!duplicated(data[data$vax_n==0,id]))
  n_vaccinated   <- sum(!duplicated(data[data$vax_n==1,id]))
  
  cat("*")
  
  # prepar.for cohort plots
  obs_per_day_all    <- obs_per_day_calc(data[!duplicated(data[,"pat_n"]),], performance=performance)
  all_obs_days       <- attributes(obs_per_day_all)$all_obs_days
  obs_per_day_all_t0 <- obs_per_day_calc(data[data$vax_n!=0,], start=paste0(vax_days,"-",vax_days),  stop=paste0(stop_obs,"-",vax_days), performance=performance)
  all_obs_days_t0    <- attributes(obs_per_day_all_t0)$all_obs_days
  
  cat("*")  
  
  # condition: with events and vaccinated
  data$cond_event_vax <- data[,event]==1 & !is.na(data[,vax_date]) & !is.na(data[,vax_name])
  
  # dataset unvaccinated deaths 
  data_unvax_deaths <- as.data.frame(data_with_deaths[is.na(data_with_deaths[,vax_date]), ])
  data_unvax_deaths <- data_unvax_deaths[!duplicated(data_unvax_deaths[,id]),]
  
  # dataset vaccinated deaths
  data_vax_deaths <- as.data.frame(data_with_deaths[!is.na(data_with_deaths[,vax_date]), ])
  rm(data_with_deaths)
  
  if(sum(data$cond_event_vax)==0) {
    plot(1,1,axes=F,type="n",xlab="",ylab="")
    text(1,1,"no data", cex=1.5) 
    return()
  }
  
  all_vax_names <- levels(factor_ref(data[data$cond_event_vax,vax_name][!duplicated(data[data$cond_event_vax,vax_name])], lab_orders = vax_def$lab_orders))
  
  vax_names_col <- c("red","magenta2","skyblue","violet","green3","darkorchid")
  if(length(vax_names_col)<length(all_vax_names)) vax_names_col <- rep(vax_names_col,length(all_vax_names)/5)
  vax_names_col <- vax_names_col[1:length(all_vax_names)]
  
  
  par(mfcol=c(2,6))
  
  at_date     <- seq( 2015,  as.numeric(substring(Sys.Date(),1,4))+1, by=0.25 ) 
  at_date_lab <- month.abb[c(1,4,7,10)][1+4*at_date%%1] 
  at_date_lab[at_date_lab==month.abb[1]] <- paste0( at_date_lab, "\n", at_date%/%1 )[at_date_lab==month.abb[1]]
  
  to_plot_date <- function(x){ (1970 + as.numeric(as.Date("2020-08-31"))/365.25) + x/365.25 }
  
  cat("*") 

  cohort_vax_plots <- list()
  # create set of plots for each value of 'vaxdep_values' 
  for(idep in 1:length(vaxdep_values)){
    
    cat(" ") 
    
    if( !(length(vaxdep_values)==1 & vaxdep_values[1]=="all") ){
      data_idep            <- data[ data$cond_event_vax & data$vaxdep_all == vaxdep_values[idep],]
      data_vax_deaths_idep <- data_vax_deaths[ data_vax_deaths$vaxdep_all == vaxdep_values[idep],]
      tit_dep <- paste( ifelse(tit=="","",paste0(tit," & ")),vaxdep_values[idep] )
    } 
    else{  # no vaccine-dependent variables (==> no split analysis)
      data_idep            <- data[data$cond_event_vax,]
      data_vax_deaths_idep <- data_vax_deaths
      tit_dep <- tit
    }
    
    limits_time_event_min_vax       <- limits_time_event_date <- limits_time_vax_date <- c(NA,NA)
    
    # per 'vax_number' or 'vax_name'  
    for(ivax in 1:length(all_vax_names)){
      
      if(nrow(data_idep)>0){
        
        cat(ivax)
        
        # create temporary variables with the current dose info:
        data_ivax <- data_idep[ data_idep[,vax_name]==all_vax_names[ivax], c(id, vax_time, vax_date, vax_name, event_time, "time_event_min_vax")]
        if(nrow(data_ivax)>0){
          data_ivax$time_event_date    <- data_ivax[,event_time];  data_ivax[,event_time]    <- NULL
          data_ivax$time_vax_date      <- data_ivax[,  vax_time]
          
          names(data_ivax)[match(c(vax_time, vax_date, vax_name, "time_event_min_vax"),names(data_ivax))] <- paste0(c(vax_time, vax_date, vax_name, "time_event_min_vax"),"_ivax")
          data_plot <- merge(data_idep[, !(names(data_idep) %in% names(data_ivax)[-1]) ], data_ivax, by=id, all.x=T )
          
          data_plot <- data_plot[data_plot[,paste0(vax_name,"_ivax")]==all_vax_names[ivax] & !is.na(data_plot[,paste0(vax_name,"_ivax")]),] 
          
          limits_time_event_min_vax <- range(data_plot$time_event_min_vax ,     na.rm=T)
          limits_time_event_date    <- range(data_plot$time_event_date, na.rm=T)
          limits_time_vax_date      <- range(data_plot$time_vax_date, na.rm=T)
        }
        else data_plot <- data_idep[0,] 
      } 
      else data_plot <- data_idep[0,] 
      
      # create temporary variables with the current dose info in 'data_vax_deaths_idep' dataset:
      if(nrow(data_vax_deaths_idep)>0){
        
        data_ivax <- data_vax_deaths_idep[ data_vax_deaths_idep[,vax_name]==all_vax_names[ivax], c(id, vax_time, vax_date, vax_name, death_time, "time_death_min_vax") ]
        if(nrow(data_ivax)>0){
          #data_ivax$time_death_min_vax <- data_ivax$death_days - data_ivax[,vax_time]
          data_ivax$time_death_date    <- data_ivax$death_days
          data_ivax$time_vax_date      <- data_ivax[,vax_time]
          data_ivax$death_days         <- NULL
          names(data_ivax)[match(c(vax_time, vax_date, vax_name, "time_death_min_vax"),names(data_ivax))] <- paste0(c(vax_time, vax_date, vax_name, "time_death_min_vax"),"_ivax")
          data_vax_deaths_plot  <- merge(data_vax_deaths_idep[, !(names(data_vax_deaths_idep) %in% names(data_ivax)[-1]) ], data_ivax, by=id, all.x=T )
          
          data_vax_deaths_plot <- data_vax_deaths_plot[ data_vax_deaths_plot[,paste0(vax_name,"_ivax")]==all_vax_names[ivax] & !is.na(data_vax_deaths_plot[,paste0(vax_name,"_ivax")]), ]
          
          limits_time_event_min_vax <- range( data_vax_deaths_plot$time_death_min_vax, limits_time_event_min_vax , na.rm=T)
          limits_time_event_date    <- range( data_vax_deaths_plot$time_death_date,    limits_time_event_date,     na.rm=T)
          limits_time_vax_date      <- range( data_vax_deaths_plot$time_vax_date,      limits_time_vax_date,       na.rm=T)
        }
        else data_vax_deaths_plot <- data_idep[0,] 
      }
      else data_vax_deaths_plot <- data_idep[0,]
      
      if( nrow(data_plot)>0 | nrow(data_vax_deaths_plot)>0 ){
        limits_time_event_min_vax <- limits_time_event_min_vax + c(-1,1)*0.1*diff(limits_time_event_min_vax)
        limits_time_event_date    <- limits_time_event_date    + c(-1,1)*0.1*diff(limits_time_event_date);   limits_time_event_date[2] <- min(as.numeric(difftime(Sys.Date(),as.Date("2020-08-31"))),limits_time_event_date[2],         na.rm=T)
        limits_time_vax_date      <- limits_time_vax_date      + c(-1,1)*0.1*diff(limits_time_vax_date);     limits_time_vax_date[2]   <- min(as.numeric(difftime(Sys.Date(),as.Date("2020-08-31"))),limits_time_vax_date[2],           na.rm=T)
        
        tit_ivax_events        <- paste0(tit_dep,  "\n", event, "; ",          trimws(all_vax_names[ivax]), " ", vaxdep_values[idep])
        tit_ivax_deaths        <- paste0(tit_dep,  "\ndeaths; ",               trimws(all_vax_names[ivax]), " ", vaxdep_values[idep])
        tit_ivax_events_deaths <- paste0(tit_dep,  "\n", event, " or death; ", trimws(all_vax_names[ivax]), " ", vaxdep_values[idep])
        
        #####################
        
        
        #####################
        #  start of plots:
        image_start <- function(time_x, time_y, distr_2d, xlab, all_vax_names, ivax, ylab_extra, at_date, at_date_lab, no_image=F){
          xlim <- range(time_x,na.rm=T)
          ylim <- range(time_y,na.rm=T); ylim[1] <- min(c(0,ylim[1])); ylim[2] <- max(c(small_range[2],ylim[2]))
          ylab <- paste0('# days from "',trimws(all_vax_names)[ivax],'"', ylab_extra)
          plot(time_x, time_y, type="n", ylab=ylab, xlab=xlab, cex.lab=0.9, axes=F, xlim=xlim, ylim=ylim )
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="#FFFFC8",border=NA)
          
          if(!is.null(distr_2d) & length(unique(distr_2d$x))>1 & !no_image ) image(distr_2d, axes=F, add=T )
          
          axis(2); axis(1, at=at_date, labels=at_date_lab, las=1); box(col="#FFFFA0")
          
          # ablines:
          abline(h=c(0,7,14,21,28, 60),lwd=0.2); abline(h=0,lwd=0.5)
          if(ivax==1) abline(h=c(-30,small_range[1]),lwd=0.2)
          abline(v=c(2015 + 0.25*((-1):100)), col="pink", lwd=0.2)
          abline(v=c(2020,2021,2022), col="pink",lwd=0.5)
          
        } # end function 'image_start'
        
        
        for(ievent_death in c("events","deaths")){
          ylab_extra <- ""
          if(ievent_death == "events") 
            tit_ivax0 <- tit_ivax_events
          else tit_ivax0 <- tit_ivax_events_deaths
          
          
          #events:
          if(nrow(data_plot)>0 | nrow(data_vax_deaths_plot)>0){ 
            
            for(i_event_vax in c("event","vax")){
              
              time_x <- time_y <- time_event_min_vax1 <- time_event_min_vax  <- vax_name_col <- pat_ns <- not_deaths <- c()
              
              if(nrow(data_plot)>0){ 
                time_x              <- c(time_x,              to_plot_date(data_plot[,paste0("time_",i_event_vax,"_date")]) )
                time_y              <- c(time_y,              data_plot$time_event_min_vax_ivax )
                time_event_min_vax1 <- c(time_event_min_vax1, data_plot$time_event_min_vax1)
                time_event_min_vax  <- c(time_event_min_vax,  data_plot$time_event_min_vax)
                vax_name_col        <- c(vax_name_col,        data_plot[,vax_name])
                pat_ns              <- c(pat_ns,              data_plot[,id])
                not_deaths          <- c(not_deaths,          rep(T,nrow(data_plot)))
              }
              
              if(ievent_death == "deaths" & nrow(data_vax_deaths_plot)>0 ){ 
                time_x              <- c(time_x,              to_plot_date(data_vax_deaths_plot[,paste0("time_",ifelse(i_event_vax=="vax","vax","death"),"_date")]) )
                time_y              <- c(time_y,              data_vax_deaths_plot$time_death_min_vax_ivax )
                time_event_min_vax1 <- c(time_event_min_vax1, rep(small_range[1],nrow(data_vax_deaths_plot)))
                time_event_min_vax  <- c(time_event_min_vax,  data_vax_deaths_plot$time_death_min_vax)
                vax_name_col        <- c(vax_name_col,        data_vax_deaths_plot[,vax_name])
                pat_ns              <- c(pat_ns,              data_vax_deaths_plot[,id])
                not_deaths          <- c(not_deaths,          rep(F,nrow(data_vax_deaths_plot)))
              }
              
              #if(is.null(time_x))  for(i in 1:6){plot(1,1,type="n",axes=F,main=tit_ivax0,cex.main=0.8); text(1,1,"no data"); next}
              
              for(iall in 1:2){ # 2 ==> smaller limits
                data_this_plot <- cbind.data.frame(time_x              = time_x              ,
                                                   time_y              = time_y              ,
                                                   time_event_min_vax1 = time_event_min_vax1 ,
                                                   time_event_min_vax  = time_event_min_vax  ,
                                                   vax_name_col        = vax_name_col        ,
                                                   pat_n               = pat_ns              ,
                                                   events              = not_deaths            )
                
                cond <- !is.na(data_this_plot$time_x) 
                if(i_event_vax=="event" & iall==2) next
                if(iall==2) cond <- cond & small_range[1] <= data_this_plot$time_event_min_vax1   &  data_this_plot$time_y <= small_range[2] # smaller window
                data_this_plot <- data_this_plot[cond,]
                
                xlab <- paste0("date of ",ifelse(i_event_vax=="event", paste0(event,ifelse(ievent_death=="events",""," or death") ), trimws(all_vax_names[ivax])))
                if(iall==2) xlab <- paste0(xlab,"\n(",small_range[1]," <= y <= ",small_range[2],")")
                
                if(nrow(data_this_plot)==0) {
                  for(i in 1:2){
                    plot(1,1,type="n",axes=F,xlab=xlab,ylab="",main=tit_ivax0,cex.main=0.8,cex.lab=0.8)
                    text(1,1,ifelse(ievent_death=="events","no events","no events and deaths"))
                  }
                  next
                }
                
                data_this_plot_uniq <- data_this_plot[!duplicated(data_this_plot[,c("pat_n","events")]),]
                #data_this_plot_uniq <- data_this_plot[!duplicated(data_this_plot$pat_n),]
                if(nrow(data_this_plot_uniq)>0){ 
                  if(length(unique(data_this_plot_uniq$time_x))>1 & length(unique(data_this_plot_uniq$time_y))>1 ) {
                    distr_2d <- try(with(data_this_plot_uniq, kde2d(x=time_x, y=time_y,  n=ngrid )),silent=T)
                    if(class(distr_2d)[[1]]=="try-error") 
                      distr_2d <- try(with(data_this_plot_uniq, kde2d(x=time_x, y=time_y,  n=ngrid, h= 4 * 1.06 * sqrt(var(time_x)) * length(time_x)^(-1/5)  )),silent=T)
                    if(class(distr_2d)[[1]]=="try-error") distr_2d <- NULL
                  }
                  else distr_2d <- NULL
                  # add points:
                  for(with_points in c(0,2)){ # c(0,1,2)
                    image_start(data_this_plot_uniq$time_x, data_this_plot_uniq$time_y, distr_2d=distr_2d, 
                                xlab, all_vax_names, ivax, ylab_extra, at_date, at_date_lab, no_image=(with_points>=2))
                    if(ievent_death=="events") tit_ivax <- paste0(tit_ivax0,"\n[",sum(data_this_plot_uniq$events), " events]")  
                    else                       tit_ivax <- paste0(tit_ivax0,"\n[",sum(data_this_plot_uniq$events), " events; ", sum(!data_this_plot_uniq$events)," deaths]") 
                    
                    if(with_points>0){ 
                      ###########
                      # sample points if they are more than 'max_n_points'
                      if(!is.na(max_n_points) & max_n_points<nrow(data_this_plot_uniq)) { 
                        sum_events   <- sum( data_this_plot_uniq$events)
                        sum_deaths   <- sum(!data_this_plot_uniq$events)
                        data_this_plot_uniq <- data_this_plot_uniq[sample.int(nrow(data_this_plot_uniq),max_n_points),]
                        sampled_ids <- data_this_plot_uniq$pat_n
                        data_this_plot <- data_this_plot[data_this_plot$pat_n %in% sampled_ids,]
                        if(ievent_death=="events") 
                          tit_ivax <- paste0(tit_ivax0,"\n[", round(100*sum( data_this_plot_uniq$events)/sum_events,0),"% of ",sum_events," events]")
                        else tit_ivax <- paste0(tit_ivax0,"\n[", round(100*sum( data_this_plot_uniq$events)/sum_events,0),"% of ",sum_events," events; ", 
                                                round(100*sum(!data_this_plot_uniq$events)/sum_deaths,0),"% of ",sum_deaths," deaths]")
                      }
                      if(any( data_this_plot_uniq$events)) with(data_this_plot_uniq[ data_this_plot_uniq$events,], points(time_x, time_y, cex=cex, pch=1) ) 
                      if(any(!data_this_plot_uniq$events)) with(data_this_plot_uniq[!data_this_plot_uniq$events,], points(time_x, time_y, cex=cex, pch=3) ) 
                      for(ivax2 in 1:length(all_vax_names)){
                        cond2 <- data_this_plot$time_event_min_vax>0 & data_this_plot$vax_name_col==all_vax_names[ivax2] & !is.na(all_vax_names[ivax2])
                        if(any(cond2 &  data_this_plot$events)) with(data_this_plot[cond2 &  data_this_plot$events,],points(time_x, time_y, cex=cex,     col=vax_names_col[ivax2], pch=1 ) )
                        if(any(cond2 & !data_this_plot$events)) with(data_this_plot[cond2 & !data_this_plot$events,],points(time_x, time_y, cex=1.2*cex, col=vax_names_col[ivax2] ,pch=3 ) )
                      }                    }
                    title(tit_ivax,cex.main=0.8)
                  }
                }  # nrow > 0
              }
            }
          } # if nrow(data_plot)>0 | nrow(data_vax_deaths_plot)>0)
          
        } # for event or deaths
      } # 2 nrows >0
    
      cat("-")
      
      #################################################################
      #
      #                   cohort plots for vaccinated with ivax
      #    
      if(!(length(vaxdep_values)==1 & vaxdep_values[1]=="all") ) cond_vax_dep <- data$vaxdep_all == vaxdep_values[idep]  else cond_vax_dep <- T
      if(any( cond_cohort<-( cond_vax_dep & data[,vax_name]==all_vax_names[ivax] )  )){  rm(cond_vax_dep)
        data_cohort <- data[cond_cohort,c(vax_days, event, event_time, death_time, stop_obs, "next_vax_days","time_event_min_vax","time_death_min_vax")]
													 
		  
																																																							
																																																							
																																																							
		  
																																																							   
																																																																				 
																																																																			   
		  
																																											   
																																																																					
																																																																					
		  
																																															   
																																																																									
																																																																									
		  
																																														  
																																																																							  
																																																																							  
		  
																																																		 
																																																																											  
																																																																											 
		  
		  
		   
																					   
        
        cohort_res <- vector("list",length=27)
        names(cohort_res) <- c("name", "vax_name", "vax_dep_name",
                               "obs_vax",                        "obs_vax_till_next",                       "obs_vax_28d",
                               "obs_vax_t0",                     "obs_vax_till_next_t0",                    "obs_vax_28d_t0",
                               "events_per_day_vax",             "events_per_day_between_vax",              "events_per_day_vax_28d",
                               "events_per_day_vax_t0",          "events_per_day_between_vax_t0",           "events_per_day_vax_28d_t0",
                               "events_or_deaths_per_day_vax",   "events_or_deaths_per_day_between_vax",    "events_or_deaths_per_day_vax_28d",
                               "events_or_deaths_per_day_vax_t0","events_or_deaths_per_day_between_vax_t0", "events_or_deaths_per_day_vax_28d_t0",
                               "deaths_per_day_vax",             "deaths_per_day_between_vax",              "deaths_per_day_vax_28d",
                               "deaths_per_day_vax_t0",          "deaths_per_day_between_vax_t0",           "deaths_per_day_vax_28d_t0"
        )     
 
        cohort_res$name                 = paste0(ifelse(length(vaxdep_values)==1 & vaxdep_values[1]=="all", "", paste0(vaxdep_values[idep],"; ")), all_vax_names[ivax])
        cohort_res$vax_name             = all_vax_names[ivax]
        cohort_res$vax_dep_name         = ifelse(length(vaxdep_values)==1 & vaxdep_values[1]=="all", "", vaxdep_values[idep])
       
        cohort_res$obs_vax              = obs_per_day_calc( data_cohort, start=vax_days,  stop=stop_obs,                                                             all_obs_days=all_obs_days, performance=performance)
        cohort_res$obs_vax_till_next    = obs_per_day_calc( data_cohort, start=vax_days,  stop=paste0('pmin(next_vax_days-1,',stop_obs,',na.rm=T)'),                 all_obs_days=all_obs_days, performance=performance)
        cohort_res$obs_vax_28d          = obs_per_day_calc( data_cohort, start=vax_days,  stop=paste0('pmin(next_vax_days-1,',stop_obs,",",vax_days,'+28,na.rm=T)'), all_obs_days=all_obs_days, performance=performance)
       
        cohort_res$obs_vax_t0           = obs_per_day_calc( data_cohort, start=paste0(vax_days,'-',vax_days),  stop=paste(stop_obs,'-',vax_days),          all_obs_days=all_obs_days_t0, performance=performance)
        cohort_res$obs_vax_till_next_t0 = obs_per_day_calc( data_cohort, start=paste0(vax_days,'-',vax_days),  stop=paste0('pmin((next_vax_days-',vax_days,'),',stop_obs,'-',vax_days,',na.rm=T)'),    all_obs_days=all_obs_days_t0, performance=performance)
        cohort_res$obs_vax_28d_t0       = obs_per_day_calc( data_cohort, start=paste0(vax_days,'-',vax_days),  stop=paste0('pmin((next_vax_days-',vax_days,'),',stop_obs,'-',vax_days,',28,na.rm=T)'), all_obs_days=all_obs_days_t0, performance=performance)
        
        # for event:
        data_cohort_event <- data_cohort[ data_cohort[,event]==1 & data_cohort[,vax_days]<=data_cohort[, event_time],]
        cohort_res$events_per_day_vax            = event_calc(data_cohort_event[, event_time]           ,  all_obs_days=all_obs_days)
        cohort_res$events_per_day_vax_t0         = event_calc(data_cohort_event[,"time_event_min_vax"],  all_obs_days=all_obs_days_t0)
        
        data_cohort_event <- data_cohort_event[ data_cohort_event[, event_time] <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],na.rm=T),]
        cohort_res$events_per_day_between_vax    = event_calc( data_cohort_event[ ,  event_time         ]  ,  all_obs_days=all_obs_days)
        cohort_res$events_per_day_between_vax_t0 = event_calc( data_cohort_event[ , "time_event_min_vax"],  all_obs_days=all_obs_days_t0)

        data_cohort_event <- data_cohort_event[ data_cohort_event[, event_time] <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],(data_cohort_event[,vax_days]+28),na.rm=T),]
        cohort_res$events_per_day_vax_28d        = event_calc( data_cohort_event[,  event_time         ]  ,  all_obs_days=all_obs_days)
        cohort_res$events_per_day_vax_28d_t0     = event_calc( data_cohort_event[, "time_event_min_vax"],  all_obs_days=all_obs_days_t0)
        
        # for event or 'death':
        data_cohort_event <- data_cohort[ (!is.na(data_cohort[, event_time]) & data_cohort[,vax_days]<=data_cohort[, event_time]) | 
                                          (!is.na(data_cohort[, death_time]) & data_cohort[,vax_days]<=data_cohort[, death_time]),]
        data_cohort_event[!is.na(data_cohort_event[, event_time]) & data_cohort_event[,vax_days]>data_cohort_event[, event_time], event_time] <- NA
        data_cohort_event$event_death_days       <- pmin(data_cohort_event[,event_time],data_cohort_event[,death_time],na.rm=T)
        data_cohort_event$event_death_min_vax    <- data_cohort_event$event_death_days - data_cohort_event[,vax_days]
        
        cohort_res$events_or_deaths_per_day_vax            = event_calc( data_cohort_event[, "event_death_days"],       all_obs_days=all_obs_days)
        cohort_res$events_or_deaths_per_day_vax_t0         = event_calc( data_cohort_event[, "event_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        data_cohort_event <- data_cohort_event[ data_cohort_event$event_death_days <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],na.rm=T),]
        cohort_res$events_or_deaths_per_day_between_vax    = event_calc( data_cohort_event[, "event_death_days"   ]  ,  all_obs_days=all_obs_days)
        cohort_res$events_or_deaths_per_day_between_vax_t0 = event_calc( data_cohort_event[, "event_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        data_cohort_event <- data_cohort_event[ data_cohort_event$event_death_days <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],(data_cohort_event[,vax_days]+28),na.rm=T),]
        cohort_res$events_or_deaths_per_day_vax_28d        = event_calc( data_cohort_event[, "event_death_days"   ]  ,  all_obs_days=all_obs_days)
        cohort_res$events_or_deaths_per_day_vax_28d_t0     = event_calc( data_cohort_event[, "event_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        
        # for 'death':
        data_cohort_event <- data_cohort[ !is.na(data_cohort[, death_time]) & data_cohort[,vax_days]<=data_cohort[, death_time],]
        cohort_res$deaths_per_day_vax            = event_calc( data_cohort_event[ , death_time],  all_obs_days=all_obs_days)
        cohort_res$deaths_per_day_vax_t0         = event_calc( data_cohort_event[, "time_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        data_cohort_event <- data_cohort_event[ data_cohort_event[, death_time] <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],na.rm=T),]
        cohort_res$deaths_per_day_between_vax    = event_calc( data_cohort_event[,  death_time         ]  ,  all_obs_days=all_obs_days)
        cohort_res$deaths_per_day_between_vax_t0 = event_calc( data_cohort_event[, "time_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        data_cohort_event <- data_cohort_event[ data_cohort_event[, death_time] <= pmin(data_cohort_event[,"next_vax_days"]-1,data_cohort_event[,stop_obs],(data_cohort_event[,vax_days]+28),na.rm=T),]
        cohort_res$deaths_per_day_vax_28d        = event_calc( data_cohort_event[,  death_time         ]  ,  all_obs_days=all_obs_days)
        cohort_res$deaths_per_day_vax_28d_t0     = event_calc( data_cohort_event[, "time_death_min_vax"],  all_obs_days=all_obs_days_t0)
        
        
        rm(data_cohort_event); rm(data_cohort)
        
        cohort_res <- list(cohort_res); names(cohort_res) <- cohort_res[[1]]$name
        cohort_vax_plots <- c( cohort_vax_plots, cohort_res )
        
      }
      
    } # for ivax
  }  # for idep
  
  cat("| *")  
  
  
  #################################################################
  #
  #
  #                   cohort plots for unvaccinated
  #
  
  obs_per_day_all <- obs_per_day_calc(data[!duplicated(data[,"pat_n"]),], performance=performance)
  all_obs_days    <- attributes(obs_per_day_all)$all_obs_days
  
  obs_per_day_unvax_periods    <- obs_per_day_calc(data[data$vax_n%in%c(0,1),], start=start_obs,  stop=paste0('pmin(',vax_days,'-1,',stop_obs,',na.rm=T)'), all_obs_days=all_obs_days, performance=performance)
  events_per_day_unvax_periods <- event_calc( data[ data[,event]==1            & ( data$vax_n==0 | ( data$vax_n==1 & data[, event_time] < pmin(data[,vax_days]-1,data[,stop_obs],na.rm=T ) )), event_time],  all_obs_days=all_obs_days)
  deaths_per_day_unvax_periods <- event_calc( data[ !is.na(data[, death_time]) & ( data$vax_n==0 | ( data$vax_n==1 & data[, death_time] < pmin(data[,vax_days]-1,data[,stop_obs],na.rm=T ) )), death_time],  all_obs_days=all_obs_days)
  
  if(file!=""){
    cohort_per_day <- c( unvaccinated=list(list(obs=obs_per_day_unvax_periods, events_per_day=events_per_day_unvax_periods, deaths_per_day=deaths_per_day_unvax_periods), cohort_vax_plots) )
    sink(paste0(file,"_cohort.txt"))
    print(cohort_per_day)
    sink()
    save(cohort_per_day, file=paste0(file,"_cohort.RData"))
  }  
  
  #### 'next_vax_days'
  if(!("next_vax_days" %in% names(data))){
    data$next_vax_days <- data[,"vax_days"] + c(diff(data[,"vax_days"]), NA)
    data[ c(data[-nrow(data),id]!=data[-1,id], F), "next_vax_days"  ] <- NA
  }
  
  cat("**")
  
  xx    <- as.Date("2020-08-31") + all_obs_days
  xx_t0 <- all_obs_days_t0
  
  tit_plot <- paste0("unvac.periods\n")
  if(tit!="") tit_plot <- paste0(tit,"\n",tit_plot)
  
  warn_opt <- options(warn=warn) 
  
  par(mfcol=c(3,3),cex.main=0.8,cex.lab=0.8)  
  
  ##############
  #  events
  if(any(data[,event]>0)){
    # events:
    plot(xx, 100000*events_per_day_unvax_periods/obs_per_day_unvax_periods,type="s", xlab=paste("Date of ",event),ylab="10^5 #events/#observed",main=paste0(tit, "; unvaccinated\n100000 * #events/ #observed"))
    try(lines(smooth.spline( xx, 100000*events_per_day_unvax_periods/obs_per_day_unvax_periods,df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, 100000*events_per_day_unvax_periods/obs_per_day_unvax_periods,df=10),col="red",lwd=3),silent=T)
    plot(xx, events_per_day_unvax_periods, type="s", xlab=paste0("Date of ",event),ylab=paste0("# ",event),main=paste0(tit, "; unvaccinated\n# ",event))
    try(lines(smooth.spline( xx, events_per_day_unvax_periods,df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, events_per_day_unvax_periods,df=10),col="red",lwd=3),silent=T)
    plot(xx, obs_per_day_unvax_periods,    type="s", xlab="Date",ylab="number of persons",main=paste0(tit,"; unvaccinated\nobserved unvaccinated days"))
    
  } else  for(i in 1:3){plot(1,1,axes=F,type="n",xlab="",ylab=""); text(1,1,paste0("no ",event), cex=1.5)} 
  
  
  ###########
  #  deaths:
  if(any(!is.na(death_time))){  
    # events and deaths:
    plot(xx, 100000*(events_per_day_unvax_periods+deaths_per_day_unvax_periods)/obs_per_day_unvax_periods,type="s", xlab=paste0("Date of ",event," or deaths"),ylab="10^5 #(events,deaths)/#observed",main=paste(tit, " unvaccinated\n100000 * #(events,deaths) / #observed"))
    try(lines(smooth.spline( xx, 100000*(events_per_day_unvax_periods+deaths_per_day_unvax_periods)/obs_per_day_unvax_periods,df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, 100000*(events_per_day_unvax_periods+deaths_per_day_unvax_periods)/obs_per_day_unvax_periods,df=10),col="red",lwd=3),silent=T)
    plot(xx, events_per_day_unvax_periods+deaths_per_day_unvax_periods, type="s", xlab=paste0("Date of ",event," or deaths"),ylab=paste0("#(events,deaths)"),main=paste(tit, "; unvaccinated\n#(",event,",deaths)") )
    try(lines(smooth.spline( xx, (events_per_day_unvax_periods+deaths_per_day_unvax_periods),df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, (events_per_day_unvax_periods+deaths_per_day_unvax_periods),df=10),col="red",lwd=3),silent=T)
    plot(xx, obs_per_day_unvax_periods,    type="s", xlab="Date",ylab="number of persons",main=paste(tit,"; unvaccinated\nobserved unvaccinated days"))
    
    # deaths:
    plot(xx, 100000*deaths_per_day_unvax_periods/obs_per_day_unvax_periods,type="s", xlab="Date of death",ylab="10^5 #deaths/#observed",main=paste0(tit, "; unvaccinated\n100000 * #deaths / #observed"))
    try(lines(smooth.spline( xx, 100000*deaths_per_day_unvax_periods/obs_per_day_unvax_periods,df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, 100000*deaths_per_day_unvax_periods/obs_per_day_unvax_periods,df=10),col="red",lwd=3),silent=T)
    plot(xx, deaths_per_day_unvax_periods, type="s", xlab="Date of death",ylab="#deaths",main=paste0(tit, "; unvaccinated\n# deaths"))
    try(lines(smooth.spline( xx, deaths_per_day_unvax_periods,df=30),col="skyblue",lwd=1),silent=T)
    try(lines(smooth.spline( xx, deaths_per_day_unvax_periods,df=10),col="red",lwd=3),silent=T)
    plot(xx, obs_per_day_unvax_periods,    type="s", xlab="Date",ylab="number of persons",main=paste0(tit,"; unvaccinated\nobserved unvaccinated days"))
  }  else for(i in 1:3){plot(1,1,axes=F,type="n",xlab="",ylab=""); text(1,1,"no deaths", cex=1.5)}  
  
  
  #######################
  #                   cohort plots for vaccinated with ivax
  #
  xx0    <- xx     
  xx0_t0 <- xx_t0     
  
  for(icoh_vax in 1:length(cohort_vax_plots)){
    
    par(mfcol=c(4,6))
    
    ###################
    # cal.time      
    for(iii in 1:3){ 
      if(iii==1){
        obs_vax                 <- cohort_vax_plots[[icoh_vax]]$obs_vax
        events_vax              <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax
        events_or_deaths_vax    <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax
        deaths_vax              <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_vax
      }
      if(iii==2){
        obs_vax              <- cohort_vax_plots[[icoh_vax]]$obs_vax_till_next
        events_vax           <- cohort_vax_plots[[icoh_vax]]$events_per_day_between_vax
        events_or_deaths_vax <- cohort_vax_plots[[icoh_vax]]$events_or_deaths_per_day_between_vax
        deaths_vax           <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_between_vax
      }
      if(iii==3){
        obs_vax              <- cohort_vax_plots[[icoh_vax]]$obs_vax_28d
        events_vax           <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax_28d
        events_or_deaths_vax <- cohort_vax_plots[[icoh_vax]]$events_or_deaths_per_day_vax_28d
        deaths_vax           <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_vax_28d
      }
      tit_plot   <- paste0(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],switch(iii,"(after)\n","(before next)\n","(during 28d)\n"))
      
      cond_tmp             <- obs_vax>=min(obs_vax[obs_vax>0])
      events_vax           <- events_vax[cond_tmp]
      events_or_deaths_vax <- events_vax[cond_tmp]
      deaths_vax           <- deaths_vax[cond_tmp]
      obs_vax              <- obs_vax[   cond_tmp]
      xx                   <- xx0[       cond_tmp]
      
      if(length(xx)==0) for(i in 1:4) {plot(1,1,type="n",axes=F,xlab="",ylab="",main=tit_plot); text(1,1,"no data")}
    
      plot(xx, 100000*events_vax/obs_vax,type="s", xlab=paste("Date of ",event),ylab="10^5 #events/#observed",main=paste0(tit_plot, "100000 * #events/ #observed"))
      try(lines(smooth.spline( xx, 100000*events_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, 100000*events_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
      plot(xx, events_vax, type="s", xlab=paste0("Date of ",event),ylab=paste0("# ",event),main=paste(tit_plot, "# ",event," (n=",sum(events_vax,na.rm=T),")"))
      try(lines(smooth.spline( xx, events_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, events_vax,df=10),col="red",lwd=3),silent=T)
      
      # cal.time; events and deaths:
      if(F){
        plot(xx, 100000*events_or_deaths_vax/obs_vax,type="s", xlab=paste0("Date of ",event," or deaths"),ylab="10^5 #(events,deaths)/#observed",main=paste(tit_plot, "100000 * #(events,deaths) / #observed"))
        try(lines(smooth.spline( xx, 100000*events_or_deaths_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
        try(lines(smooth.spline( xx, 100000*events_or_deaths_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
        plot(xx, events_or_deaths_vax, type="s", xlab=paste0("Date of ",event," or deaths"),ylab=paste0("#(events,deaths)"),main=paste(tit_plot,  "#(",event,",deaths (n=",sum(events_vax+deaths_vax,na.rm=T),"))"))
        try(lines(smooth.spline( xx, events_or_deaths_vax,df=30),col="skyblue",lwd=1),silent=T)
        try(lines(smooth.spline( xx, events_or_deaths_vax,df=10),col="red",lwd=3),silent=T)
      }
      
      # cal.time; deaths:
      plot(xx, 100000*deaths_vax/obs_vax,type="s", xlab="Date of death",ylab="10^5 #deaths/#observed",main=paste0(tit_plot, "100000 * #deaths / #observed"))
      try(lines(smooth.spline( xx, 100000*deaths_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, 100000*deaths_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
      plot(xx, deaths_vax, type="s", xlab="Date of death",ylab="#deaths",main=paste0(tit_plot,  "# deaths (n=",sum(deaths_vax,na.rm=T),")"))
      try(lines(smooth.spline( xx, deaths_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, deaths_vax,df=10),col="red",lwd=3),silent=T)
    }
    #
    for(iii in 1:3){ 
      if(iii==1){
        obs_vax_t0              <- cohort_vax_plots[[icoh_vax]]$obs_vax_t0
        events_vax_t0           <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax_t0
        events_or_deaths_vax_t0 <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax_t0
        deaths_vax_t0           <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_vax_t0
      }
      if(iii==2){
        obs_vax_t0              <- cohort_vax_plots[[icoh_vax]]$obs_vax_till_next_t0
        events_vax_t0           <- cohort_vax_plots[[icoh_vax]]$events_per_day_between_vax_t0
        events_or_deaths_vax_t0 <- cohort_vax_plots[[icoh_vax]]$events_per_day_between_vax_t0
        deaths_vax_t0           <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_between_vax_t0
      }
      if(iii==3){
        obs_vax_t0              <- cohort_vax_plots[[icoh_vax]]$obs_vax_28d_t0
        events_vax_t0           <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax_28d_t0
        events_or_deaths_vax_t0 <- cohort_vax_plots[[icoh_vax]]$events_per_day_vax_28d_t0
        deaths_vax_t0           <- cohort_vax_plots[[icoh_vax]]$deaths_per_day_vax_28d_t0
      }
      tit_plot   <- paste0(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],switch(iii,"(after)\n","(before next)\n","(during 28d)\n"))
      
      #####################
      #  t0 is vax date;
      cond_tmp             <- obs_vax_t0>=min(obs_vax_t0[obs_vax_t0>0])
      events_vax           <- events_vax_t0[cond_tmp]
      events_or_deaths_vax <- events_vax_t0[cond_tmp]
      deaths_vax           <- deaths_vax_t0[cond_tmp]
      obs_vax              <- obs_vax_t0[   cond_tmp]
      xx                   <- xx0_t0[       cond_tmp]
      
      if(length(xx)==0) for(i in 1:4) {plot(1,1,type="n",axes=F,main=tit_plot); text(1,1,"no data")}
      
      #  t0 is vax date;  events
      plot(xx, 100000*events_vax/obs_vax,type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab="10^5 #events/#observed",main=paste0(tit_plot, "100000 * #events/ #observed"))
      try(lines(smooth.spline( xx, 100000*events_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, 100000*events_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
      plot(xx, events_vax, type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab=paste0("# ",event),main=paste(tit_plot, "# ",event," (n=",sum(events_vax,na.rm=T),")"))
      try(lines(smooth.spline( xx, events_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, events_vax,df=10),col="red",lwd=3),silent=T)
      
      #  t0 is vax date;  events and deaths:
      if(F){
        plot(xx, 100000*events_or_deaths_vax/obs_vax,type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab="10^5 #(events,deaths)/#observed",main=paste0(tit_plot, "100000 * #(events,deaths) / #observed"))
        try(lines(smooth.spline( xx, 100000*events_or_deaths_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
        try(lines(smooth.spline( xx, 100000*events_or_deaths_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
        plot(xx, events_or_deaths_vax, type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab=paste0("#(events,deaths)"),main=paste0(tit_plot, "#(",event,",deaths (n=",sum(events_vax+deaths_vax,na.rm=T),"))"))
        try(lines(smooth.spline( xx, events_or_deaths_vax,df=30),col="skyblue",lwd=1),silent=T)
        try(lines(smooth.spline( xx, events_or_deaths_vax,df=10),col="red",lwd=3),silent=T)
      }
      
      #  t0 is vax date;  deaths:  
      plot(xx, 100000*deaths_vax/obs_vax,type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab="10^5 #deaths/#observed",main=paste0(tit_plot, "100000 * #deaths / #observed"))
      try(lines(smooth.spline( xx, 100000*deaths_vax/obs_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, 100000*deaths_vax/obs_vax,df=10),col="red",lwd=3),silent=T)
      plot(xx, deaths_vax, type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name),ylab="#deaths",main=paste0(tit_plot, "# deaths (n=",sum(deaths_vax,na.rm=T),")"))
      try(lines(smooth.spline( xx, deaths_vax,df=30),col="skyblue",lwd=1),silent=T)
      try(lines(smooth.spline( xx, deaths_vax,df=10),col="red",lwd=3),silent=T)
    }
  }
  options(warn_opt) 
  
  ###########################
  # plots for #observed:
  par(mfrow=c(4,6))
  for(icoh_vax in 1:length(cohort_vax_plots)){
    tit_plot   <- paste0(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],switch(iii,"(after)\n","(before next)\n","(during 28d)\n"))
    # cal.time
    obs_vax    <- cohort_vax_plots[[icoh_vax]]$obs_vax
    cond_tmp <- obs_vax>=min(obs_vax[obs_vax>0])
    plot(xx0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab="Date",ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(after)\nobserved days"))
    plot(xx0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab="Date",ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(before next)\nobserved days"))
    plot(xx0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab="Date",ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(during 28d)\nobserved days"))
    
    # t0 is vax date
    obs_vax    <- cohort_vax_plots[[icoh_vax]]$obs_vax_t0
    cond_tmp <- obs_vax>=min(obs_vax[obs_vax>0])
    plot(xx0_t0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name), ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(after)\nobserved days"))
    plot(xx0_t0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name), ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(before next)\nobserved days"))
    plot(xx0_t0[cond_tmp], obs_vax[cond_tmp],    type="s", xlab=paste0("Days from ",cohort_vax_plots[[icoh_vax]]$vax_name), ylab="number of persons",main=paste(tit,"; vaccinated\n",names(cohort_vax_plots)[icoh_vax],"(during 28d)\nobserved days"))
  }
  if(F){   
    if(length(cohort_vax_plots)>0){
      # vaccinated deaths
      if(iset=="incidence"){
        plot(xx, 100000*n_vax_deaths_per_day/n_vax_per_day, type="s", xlab="Date of death",ylab="10^5 #deaths/#observed",main=paste(tit, "; vaccinated\n100000 * #deaths / #observed"))
        try(lines(smooth.spline( xx, 100000*n_deaths_per_vax1_day/n_per_vax1_day,df=10),col="red",lwd=3))
      }
      if(iset=="events"){
        plot(xx, n_vax_deaths_per_day,type="s", xlab="Date of death",ylab="#deaths",main=paste(tit, "; vaccinated\ndeaths"))
        plot(xx, n_vax_per_day,type="s", xlab="Observation date",ylab="number of persons",main=paste(tit,"; vaccinated\nobserved per day after dose 1"))
      }
      
      n_deaths_per_vax1_day <- rep(0,length(obs_vax1_days)); names(n_deaths_per_vax1_day) <- obs_vax1_days
      tb <- table(data_deaths[data_deaths$vax_n==1,"death_days"]-data_deaths[data_deaths$vax_n==1,"vax_days"])
      n_deaths_per_vax1_day[match(names(tb),names(n_deaths_per_vax1_day))] <- tb
      xx <- as.numeric(names(n_deaths_per_vax1_day))
      if(iset=="incidence"){
        plot(xx, 100000*n_deaths_per_vax1_day/n_per_vax1_day, type="s", xlab="Days from dose 1",ylab="10^5 #deaths/#observed",main=paste(tit, "; vaccinated\n100000 * #deaths / #observed"))
        try(lines(smooth.spline( xx, 100000*n_deaths_per_vax1_day/n_per_vax1_day,df=10),col="red",lwd=3))
      }
      if(iset=="events"){
        plot(xx, n_deaths_per_vax1_day,type="s", xlab="Days from dose 1",ylab="#deaths",main=paste(tit, "; vaccinated\ndeaths"))
        plot(xx, n_per_vax1_day,type="s", xlab="Days from dose 1", ylab="number of persons",main=paste(tit,"; vaccinated\nobserved per day after dose 1"))
      }
    }
  } 
  par(cex.main=1,cex.lab=1)
  cat(paste0(" ==> duration = ",format(Sys.time()-ttt,digits=2), " (till ", Sys.time()),") ]\n")
  
}  # the end of function 'brand_images'

# obs calculation


# obs calculation
obs_per_day_calc <- function(data, all_obs_days, start="study_entry_days", stop="study_exit_days", performance=F){  
  tmp<-Sys.time() 
  if( !(is.numeric(start) & length(start)==1) ) start <- with(data, eval(parse(text=start)))
  if( !(is.numeric(stop)  & length(stop )==1) ) stop  <- with(data, eval(parse(text=stop)))
  nrow_diff <- nrow(data)
  data <- data[!is.na(start) & !is.na(stop) & start<=stop,]
  nrow_diff <- nrow_diff - nrow(data)
  #print( start); print(stop)
  if(missing(all_obs_days)) all_obs_days <- min(start,na.rm=T):max(stop,na.rm=T)
  obs_per_day <- rep(0,length(all_obs_days))
  
  if(length(start)==1 & length(start)<length(stop)) obs_per_day[1] <- nrow(data)
  else obs_per_day[1] <- sum( start==all_obs_days[1], na.rm=T)
  for(iday in all_obs_days[-1] ){ #cat(iday)
    obs_per_day[all_obs_days==iday] <- obs_per_day[all_obs_days==iday-1] + sum( start==iday, na.rm=T) -  
      sum( stop==iday-1, na.rm=T)
  }
  attributes(obs_per_day) <- list(all_obs_days=all_obs_days, start_gt_stop_or_missing_rows=nrow_diff)
  if(performance) print(Sys.time()-tmp)
  obs_per_day
}

obs_per_day_calc_old <- function(data, all_obs_days, start="study_entry_days", stop="study_exit_days", performance=F){
  tmp<-Sys.time()
  start <- with(data, eval(parse(text=start)))
  stop  <- with(data, eval(parse(text=stop)))
  if(missing(all_obs_days)) all_obs_days <- min(start,na.rm=T):max(stop,na.rm=T)
  obs_per_day <- rep(0,length(all_obs_days))
  #names(obs_per_day) <- all_obs_days 
  #print(length(all_obs_days))
  #days <- sort(unique(start,stop))
  obs_per_day[1] <- sum( start==all_obs_days[1], na.rm=T)
  for(iday in all_obs_days[-1] ){ #cat(iday)
    obs_per_day[all_obs_days==iday] <- obs_per_day[all_obs_days==(iday-1)] + sum( start==iday, na.rm=T) -  
      sum( stop==(iday-1), na.rm=T)
    #if(iday %in% days) obs_per_day[all_obs_days==iday] <- obs_per_day[all_obs_days==iday] + sum( start==iday, na.rm=T) - 
    #    sum( stop==iday-1, na.rm=T)
  }
  #obs_per_day <- obs_per_day + sum( start==all_obs_days[1], na.rm=T)
  attributes(obs_per_day) <- list(all_obs_days=all_obs_days)
  if(performance) print(Sys.time()-tmp)
  obs_per_day
}

# events calculation  
event_calc <- function(event_time_var,all_obs_days){
  events_tab <- rep(0,length(all_obs_days))
  tb <- table(event_time_var)
  if(length(tb)>0){ 
    tb_part <- match(as.numeric(names(tb)),all_obs_days)
    if(length(tb_part)>0) {
      tb <- tb[!is.na(tb_part)]
      tb_part <- tb_part[!is.na(tb_part)]
    }
    if(length(tb_part)>0) events_tab[tb_part] <- tb[!is.na(tb_part)]
  }
  events_tab
}


if(F){
  Sys.time()
  res <- scri( formula = "~ lab", vax_def = vax_def,  data = data_vax, event_info=event_info,  extra_parameters = extra_options, add_to_itself=F, lplots=F, leventplot = T, max_n_points = 500  )  
  Sys.time()
}
#
#########################################################################################


###############################################
# create forest plots 
#   and 
# print markdown tables

forest_plots_tab <- function(tab_list, ndigits=2, nrows_forest_plot=40, cex=0.8, cex_head=0.8, 
                             col=c("red", "green3", "orange",  "deepskyblue", "magenta2", "cyan2", "chocolate1", "gray" ), 
                             lplot=T, ltable=T){
  
  plotted <- F
  if(is.null(tab_list)   ) return(plotted)
  if(length( tab_list)==0) return(plotted)
  
  if(!is.list(tab_list)) tab_list <- list(tab_list)
  
  if(all(sapply(tab_list,is.null)                    )) return(plotted)
  if(all(sapply(tab_list,function(x)length(dim(x)))<2)) return(plotted)
  if(all(sapply(tab_list,nrow)==0                    )) return(plotted)
  
  if(is.null(names(tab_list))) names(tab_list) <- 1:length(tab_list)
  if(length(col)<length(tab_list))col <- rep(col,length(tab_list))[1:length(tab_list)]
  
  i_without_time <- 1
  if(any(substring(names(tab_list),1,6)=="no_adj"))
    i_without_time <- (1:length(tab_list))[substring(names(tab_list),1,6)=="no_adj"][1]
  nrow_without_time <- nrow(tab_list[[i_without_time]])
  
  
  var_RR <- rep("RR_data",length(tab_list))
  RR_max <- NA
  
  for(i in 1:length(tab_list)){
    
    if(    is.null(tab_list[[i]])    ) next
    if( length(dim(tab_list[[i]]))<2 ) next
    if(       nrow(tab_list[[i]])==0 ) next
    
    tt0 <- tab_list[[i]][1:nrow_without_time, ]
    row.names(tt0) <- NULL
    tt <- tt0[, !(names(tt0) %in%  c("event","i"))]
    
    tt <- as.data.frame(lapply(tt,function(x,digits){if(mode(x)=="numeric")x<-round(x,digits=digits);x},digits=ndigits))
    tt$group_start <- 1:nrow(tt)
    tt$group_start[!grepl("reference",tt$all_cat) & substring(tt$all_cat,1,1)!="[" & c(F,substring(tt$all_cat[-1],1,8) == substring(tt$all_cat[-nrow(tt)],1,8))] <- NA
    #tt$group_start[!grepl("reference",tt$all_cat) & c(F,substring(tt$all_cat[-1],1,3) == substring(tt$all_cat[-nrow(tt)],1,3))] <- NA
    tt$group_start[ !c( T, tt$group_start[-1] - tt$group_start[-nrow(tt)] >1 ) ] <- NA
    
    tt$group <- 1
    tt$group <- rep( tt$group_start[!is.na(tt$group_start)], c(tt$group_start[!is.na(tt$group_start)][-1],nrow(tt)+1) - tt$group_start[!is.na(tt$group_start)] )
    
    no_RR <- F
    if(!any("lci" %in% names(tt))) no_RR <- T
    
    if(!no_RR){
      var_RR[i]   <- "RR"  
      
      tt$RR_CI <- paste0(format(tt$RR,nsmall=ndigits)," (",format(tt$lci,nsmall=ndigits),";",format(tt$uci,nsmall=ndigits),")")
      tt$RR_CI[is.na(tt$RR)] <- ""
    }
    tab_list[[i]] <- tt
    if(var_RR[i] %in% names(tt)) 
      if(any(!is.na(tt[,var_RR[i]]))) RR_max <- max(tt[,var_RR[i]], RR_max, na.rm=T)
  }
  
  if(lplot){
    
    tt <- tab_list[[i_without_time]]
    
    add_text <- function(){  
      # headings:
      op <- par(cex=cex_head, font=2)
      text(-( c(4,2.2)*10/(40/max(1,RR_max)) ),  1 + nrows_forest_plot/40, c("           # event in",  "# observed days per person in"), offset=0, pos=4, xpd=T )
      text(-(  ( 4:1 )*10/(40/max(1,RR_max)) ),  1                       , c("risk wind", "ref.wind",  "risk wind",    "ref.wind"     ), offset=0, pos=4 )
      par(font=4)  # bold italic font
      par(op) 
    }
    
    groups <- tt$group[!duplicated(tab_list[[1]]$group)]
    #groups <- tt$group[!duplicated(tt$group)]
    
    start_y <- 1
    #if(no_RR) print("??????????????????????????CI's are not correct!!! There are no estimates from Poisson regression!!!")
    
    for(igr in 1:length(groups)){
      
      if(sum(tt$group == groups[igr])==0) next
      
      cond_group <- !is.na(tt[,var_RR[i_without_time]]) & tt$group == groups[igr]
      if(sum(cond_group)==0) next
      
      group_part_start_row <- 1
      group_part_end_row   <- sum(cond_group)
      
      while(T){
        if(   start_y + (group_part_end_row-group_part_start_row+1)   > 0.9*nrows_forest_plot | start_y > nrows_forest_plot ){ start_y <- 1; par(new=F) }
        if( ( start_y + (group_part_end_row-group_part_start_row+1) ) > nrows_forest_plot                                   )  group_part_end_row <- nrows_forest_plot - start_y + group_part_start_row
        
        IRR_tit <- format("IRR (CI)",width=15)
        if(no_RR){
          IRR_tit <- format("IRR (CI*)",width=15)
          tt$lci <- pmax(0, tt[,var_RR[i_without_time]] - 0.1 )
          tt$uci <- (tt[,var_RR[i_without_time]] + 0.1)
          IRR_tit <- "IRR from data"
        }
        if(any(!is.na(tt[cond_group,][group_part_start_row:group_part_end_row,var_RR[i_without_time]]))){
          with(tt[cond_group,][group_part_start_row:group_part_end_row,], { 
            forest.default(slab=all_cat, x=get(var_RR[i_without_time]), ci.lb=lci, ci.ub=uci,
                           refline=1, cex=cex, 
                           header = c("Intervals",IRR_tit),
                           xlim=c(-(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), max(1,RR_max)*1.6), 
                           ylim=c(2, -nrows_forest_plot) ,
                           alim=c(  0, max(1,RR_max)*1.1),
                           ilab=cbind(events_rw,events_ref, days_per_id_rw, days_per_id_ref), 
                           ilab.xpos= -(  (4:1)*10/(40/max(1,RR_max)) ), rows= -start_y
            )} )
          plotted <- T
        }
        par(new=T) 
        
        legend_names <- names(tab_list)
        for(itab in 1:length(tab_list) ){
          
          cond_group <- !is.na(tt[,var_RR[i_without_time]]) & tt$group == groups[igr]
          
          names_as_unadj <- tab_list[[itab]]$all_cat[ match( tt[cond_group,][group_part_start_row:group_part_end_row,"all_cat"] , tab_list[[itab]]$all_cat ) ]
       
          if(itab!=i_without_time & !identical( names_as_unadj , tt[cond_group,][group_part_start_row:group_part_end_row,"all_cat"] )  )  { 
            warning(paste0("The main effect names of tables '",tt[cond_group,][group_part_start_row:group_part_end_row,"all_cat"]," and '",names_as_unadj,"' are not identical"))
            #legend_names <- legend_names[-itab]
            next
          }
    if(F) if(itab!=i_without_time & !identical( tab_list[[itab]]$all_cat[1:nrow_without_time] , tab_list[[i_without_time]]$all_cat ) ) { 
            warning(paste0("The main effect names of tables '",names(tab_list)[i_without_time]," and '",names(tab_list)[itab],"' are not identical"))
            legend_names <- legend_names[-itab]
            next
          }
          if(!no_RR & any(!is.na(tab_list[[itab]][match(names_as_unadj,tab_list[[itab]][,"all_cat"]),var_RR[itab]]))){
            with(tab_list[[itab]][match(names_as_unadj,tab_list[[itab]][,"all_cat"]),], { 
              forest.default(slab=rep("",length(all_cat)), x=get(var_RR[itab]), ci.lb=lci, ci.ub=uci,
                             header = F, annotate=F, cex=cex, refline=1, 
                             xlim=c(-(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), max(1,RR_max)*1.6), 
                             ylim=c(2, -nrows_forest_plot) ,
                             alim=c(  0, max(1,RR_max)*1.1),
                             rows= -start_y + 0.8*(itab-1)/length(tab_list), col=col[itab]
              )} )
            plotted <- T
          }
          par(new=T) 
        }
        
        if(start_y==1) add_text()
        legend( -(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), -1.04*nrows_forest_plot,legend=legend_names,
                box.col="lightgray", col=col[match(legend_names,names(tab_list))],xpd=T,lty=1,pch=15, cex=cex, xjust=0,yjust=1,ncol=max(1,round(length(legend_names)/3+0.1)))
        #legend( "bottomleft", inset=c(0,-0.05),legend=rep(legend_names,4),col=rep(col[match(legend_names,names(tab_list))],4),xpd=T,lty=1,pch=15, cex=cex, xjust=1,yjust=0,ncol=round(length(legend_names)/4+0.1))
        if(no_RR) text( -(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), -nrows_forest_plot, "* CI's are not correct!!! There are no estimates from Poisson regression!!!", cex=0.5,pos=4)
        par(new=T) 
        
        start_y <- start_y + group_part_end_row - group_part_start_row + 2
        group_part_start_row <- group_part_end_row + 1
        if(group_part_start_row > sum(cond_group)) break
        group_part_end_row <- group_part_start_row + ( sum(cond_group)-group_part_start_row+1) - 1
      }
    }
    par(new=F) 
  }
  
  if(ltable){
    for(i in 1:length(tab_list)){
      tt<-tab_list[[i]]
      tmp <- tt[!is.na(tt$group_start),]
      tmp$group <- tmp$group - 0.5
      tmp[,names(tmp)!="group"] <- ""
      tt_print <- rbind.data.frame(tt,tmp)
      tt_print <- tt_print[order(tt_print$group),]
      row.names(tt_print) <- NULL
      
      #tt_print$days_per_id_rw <- format(tt_print$days_per_id_rw,digits=2,nsmall=2)
      #tt_print$pval <- sprintf(".2%s",tt_print$pval )
      #tt <- as.data.frame(lapply(tt,function(x,digits){if(mode(x)=="numeric")x<-format(x,digits=ndigits);x},digits=ndigits))
      
      
      if(no_RR) print("CI's are not correct!!! There are no estimates from Poisson regression!!!")
      if(i==1){
        # columns description:
        cat('Column names:\n')
        cat('   - #ev.risk"      - the number of events in the risk window                                      \n'  )
        cat('   - #ev.ref."      - the number of events in the corresponding reference window                   \n'  )
        cat('   - #days_id_risk" - the number of observed days in the risk window per person                    \n'  )
        cat('   - #days_id_ref." - the number of observed days in the corresponding reference window per person \n\n')
      }
      # print table:
      if(!no_RR)
        knitr::kable( tt_print[,c("all_cat", "events_rw","events_ref", "days_per_id_rw", "days_per_id_ref", "RR_CI", "pval" )],
                      col.names = c("Intervals", "#ev.risk","#ev.ref.","#days_id_rw", "#days_id_ref", "IRR(CI)", "pval" ),
                      align="lrrrrrr"  )
      else
        knitr::kable( tt_print[,c("all_cat", "events_rw","events_ref", "days_per_id_rw", "days_per_id_ref", "RR_data")],
                      col.names = c("Intervals", "#ev.risk","#ev.ref.","#days_id_rw", "#days_id_ref", "IRR_from_data" ),
                      align="lrrrrrr"  )
    }
  }
  invisible(plotted)
} # end of the function


#################################################################
###############################    ####################
scri_data_parameters <- function(...){ 
  
  res   <- list()
  param <- list(...)
  if(any(names(param)=="data")) data <- param[["data"]]
  
  if(F){              
    # specify 'nvax':
    if(all(names(param)!="nvax")){
      if(any(names(param)=="data")) 
        nvax <- max(as.numeric(gsub("[^0-^9]+","",names(data)[substring(names(data),1,3)=="vax"])),na.rm=T)
      else nvax <- 3
    }
    else nvax <- param[["nvax"]]
    res$nvax <- nvax
  }       
  
  # parameters for dataset:
  if(any(names(param)=="data")){
    res$data_parameters <- list()
    if(F){        
      # specify vaccination time variable 'vax_time':
      tmp <- tmp0 <- c("vax%n_time",  "vax%n_t",  "days_vax%n","vax%n_days","time_vax%n","t_vax%n")
      if(any(names(param)=="vax_time")) tmp <- param[["vax_time"]]
      if(any( unlist(lapply(1:res$nvax,function(ivax) sum(tolower(names(data)) %in% gsub("%n",ivax,tmp))))!=1 )) 
        stop(paste0("Specify the time of vaccination variable 'vax_time' (default,'",paste0(tmp0,collapse="' or '"),"'). Use %n instead of the numbers in the vaccine's names."))
      else res$data_parameters$vax_time <- unlist(lapply(1:res$nvax,function(ivax) names(data)[tolower(names(data)) %in% gsub("%n",ivax,tmp)])) 
    }          
    
											  
    res$data_parameters$vax_name <- "vax_name"
    res$data_parameters$vax_time <- "vax_time"
    res$data_parameters$vax_date <- "vax_date"
    
    if("vax1" %in% names(param)) res$data_parameters$vax1 <- param[["vax1"]]
    else stop("parameter 'vax1' is missing.")
    
    res$data_parameters$vax1     <- res$data_parameters$vax1

    # specify the vax_name variable, 'vax_name':
    tmp <- tmp0 <- c("vax_name","name_vax")
    if(any(names(param)=="vax_name")) tmp <- param[["vax_name"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_name' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_name <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_name_v1 variable, 'vax_name_v1':
    if(any(names(param)=="vax_name_v1")){ 
      if(param[["vax_name_v1"]] %in% names(data)) res$data_parameters$vax_name_v1 <- names(data)[tolower(names(data)) %in% param[["vax_name_v1"]]]
      else stop(paste0("variable '",param[["vax_name_v1"]],"' is not found."))  
    }
    else res$data_parameters$vax_name_v1 <- res$data_parameters$vax_name
    
    if(!(res$data_parameters$vax1 %in% (tmp<-unique(data[,res$data_parameters$vax_name_v1])))){ 
      if(trimws(res$data_parameters$vax1) %in% trimws(tmp)) 
        res$data_parameters$vax1 <- tmp[trimws(tmp)==trimws(res$data_parameters$vax1)]
    }

    
    # specify the vax_time variable, 'vax_time':
    tmp <- tmp0 <- c("vax_time","time_vax",  "vax_t","t_vax", "days_vax", "vax_days")
    if(any(names(param)=="vax_time")) tmp <- param[["vax_time"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_time' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_time <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_time_v1 variable, 'vax_time_v1':
    if(any(names(param)=="vax_time_v1")){
      if( param[["vax_time_v1"]] %in% names(data)) res$data_parameters$vax_time_v1 <- names(data)[tolower(names(data)) %in% param[["vax_time_v1"]]]
      else stop(paste0("variable '",param[["vax_time_v1"]],"' is not found."))
    }
    else res$data_parameters$vax_time_v1 <- res$data_parameters$vax_time
    
    # specify the vax_date variable, 'vax_date':
    tmp <- tmp0 <- c("vax_date","date_vax", "vax_time","time_vax",   "vax_t","t_vax", "days_vax")
    if(any(names(param)=="vax_date")) tmp <- param[["vax_date"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_date' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_date <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_date_v1 variable, 'vax_date_v1':
    if(any(names(param)=="vax_date_v1")){ 
      if(param[["vax_date_v1"]] %in% names(data)) res$data_parameters$vax_date_v1 <- names(data)[tolower(names(data)) %in% param[["vax_date_v1"]]]
      else stop(paste0("variable '",param[["vax_date_v1"]],"' is not found."))  
    }
    else res$data_parameters$vax_date_v1 <- res$data_parameters$vax_date
    
    # specify the id variable, 'id':
    tmp <- tmp0 <- c("pat_n","id","person_id","pat","pid" )
    if(any(names(param)=="id")) tmp <- param[["id"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'id' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$id <- names(data)[tolower(names(data)) %in% tmp]
    
    
    # specify the start of observation 'start_obs':
    tmp <-  tmp0 <- c("start_obs","start_obs_time","t_start_obs","start_study","start_study_time","t_start_study","study_entry_days","study_entry_time","t_study_entry")
    if(any(names(param)=="start_obs")) tmp <- param[["start_obs"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the start of observation variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$start_obs <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the end of observation 'end_obs':
    tmp <- tmp0 <- c("end_obs","end_obs_time","t_end_obs","end_study","end_study_time","t_end_study","study_exit_days","study_exit_time","t_study_exit")
    if(any(names(param)=="end_obs")) tmp <- param[["end_obs"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the end of observation variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$end_obs <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify censored variables, for example, death_days:
    tmp <- c("death_time","death_days","days_death")
    if(any(names(param)=="censored_vars")) tmp <- param[["censored_vars"]]
    if(sum(tolower(names(data)) %in% tmp)>0) res$data_parameters$censored_vars <- tmp[tmp %in% names(data)]
    else res$data_parameters$censored_vars <- c()
    
  }
  else{
    if(any(names(param)=="vax_time"))      res$data_parameters$vax_time      <- param[["vax_time"]]      else res$data_parameters$vax_time      <- paste0("vax",1:res$nvax,"_days") 
    if(any(names(param)=="id"))            res$data_parameters$id            <- param[["id"]]            else res$data_parameters$id            <- "id" 
    if(any(names(param)=="start_obs"))     res$data_parameters$start_obs     <- param[["start_obs"]]     else res$data_parameters$start_obs     <- "start_obs"  
    if(any(names(param)=="end_obs"))       res$data_parameters$end_obs       <- param[["end_obs"]]       else res$data_parameters$end_obs       <- "end_obs"  
    if(any(names(param)=="censored_vars")) res$data_parameters$censored_vars <- param[["censored_vars"]] else res$data_parameters$censored_vars <- "death_time"   
  }
  
  class(res) <- "scri_parameters"
  res
}  # end of function 'scri_data_parameters'

define_rws <- function( obj,
                        cond_vax_var,  cond_vax_values,  vax_names_short, 
                        cond_vax1_var, cond_vax1_values, vax1_names_short, # for periods before vax1
                        before_vax_lab=c( "reference", "buffer" ), after_vax_lab,
                        #before_vax_lab=c( "pre-exposure", "buffer" ), after_vax_lab,
                        ref = "reference",
                        #ref = "pre-",
                        data,
                        lab_order = T,
                        ...                        #  for names as: 'begin_cut_points', 'after_overlap_priority', ...
){ 
  
  names_before_after <- c( "name","t0","cond","cond_var","cond_value","vax_dep","cut_points", "lab", "ref", "lab_add_interval", "no_last_interval", "overlap_priority_after" )
  names_rws_short    <- c("before","after", "rws_def")
  names_rws          <- c( names_rws_short, names_before_after, paste0("before_",names_before_after), paste0(names_before_after,"_before"), paste0("after_",names_before_after) , paste0(names_before_after,"_after")  )
  extra_pars         <- list(...)
  
  if(!missing(obj)){
    if(!is.list(obj)){ 
      if(mode(obj)=="character"){ if(substring(obj,1,4)=="list") obj<-list(rws=list(rws_def=obj)) }
      else stop("'obj' should be a list.")
    }  
    if(!("rws" %in% names(obj))) 
      if(any( names_rws_short %in% names(obj) )) {obj <- list(rws=obj); class(obj) <- "scri_parameters"}
    else if(any(names_rws %in% names(obj))) {}  # not ready yet
  }
 
  if(missing(cond_vax_var)){
    if(!is.null(obj$data_parameters$vax_name))  cond_vax_var <-  obj$data_parameters$vax_name
    else {
      if(missing(cond_vax_var) & missing(data)) stop("Define 'cond_vax_var' or/and 'data'.")
      # count vax  1,2,3,... !!!
      nvax <- 1
      ids <- data[,obj$data_parameters$id]
      while(T){
        ids <- ids[duplicated(ids)]
        if(length(ids)==0) break
        nvax <- nvax+1
        
      }
    }
  }

  
  if(missing(cond_vax_values))
    if(!missing(data)){
      if( !any(cond_vax_var %in% names(data))) stop(paste0("variable '",cond_vax_var,"' not found in dataset 'data'."))
      cond_vax_values_all <- data[data$vax_n>0,cond_vax_var]
      if(!is.null(obj$data_parameters$vax_time)) cond_vax_values_all <- cond_vax_values_all[order(data[data$vax_n>0,obj$data_parameters$vax_time])]
      cond_vax_values <- cond_vax_values_all[!duplicated(cond_vax_values_all)]
      cond_vax_var    <- rep(cond_vax_var, length(cond_vax_values))
      t0    <- rep(obj$data_parameters$vax_time, length(cond_vax_values))
      #cond_vax_values <- c( cond_vax_values[1], cond_vax_values )  # the first value for reference period ==> condition tov the start of the first vax

      cond_vax_values <- c( obj$data_parameters$vax1,        cond_vax_values )
      cond_vax_var    <- c( obj$data_parameters$vax_name_v1, cond_vax_var    )   
      t0              <- c( obj$data_parameters$vax_time_v1, t0              )
      
      
      if(missing(vax_names_short)){ 
        vax_names_short <- c()
        for(ichar in  substring(cond_vax_values[-1],1,1)[!duplicated(substring(cond_vax_values[-1],1,1))] ) 
          vax_names_short <- c( vax_names_short, paste0( ichar, 1:sum(substring(cond_vax_values[-1],1,1)==ichar) )  )
      }
      else {
        if(mode(vax_names_short)=="character" & length(vax_names_short)==1){
          if( !any(vax_names_short %in% names(data))) stop(paste0("variable '",vax_names_short,"' not found in dataset 'data'."))
          vax_names_short <- data[ !duplicated(cond_vax_values_all), vax_names_short]
        }
      }      
      if(missing(vax1_names_short)) vax_names_short <- c( paste0(substring(cond_vax_values[1],1,1),"0"), vax_names_short)
      else vax_names_short <- c(vax1_names_short, vax_names_short)
      
      if(missing(after_vax_lab)) after_vax_lab <- cond_vax_values[-1]
      else {
        if(mode(after_vax_lab)=="character" & length(after_vax_lab)==1){
          if( !any(after_vax_lab %in% names(data))) stop(paste0("variable '",after_vax_lab,"' not found in dataset 'data'."))
          after_vax_lab <- data[ !duplicated(cond_vax_values_all), after_vax_lab]
        }
      }
	  
    }

  
  
  if(length(vax_names_short)!=length(cond_vax_values)) stop(paste0("The length of 'vax_names_short' (",length(vax_names_short),") should be the same as the length of 'cond_vax_value' (",length(cond_vax_values),")"))
  if(length(vax_names_short)!=length(after_vax_lab)+1) stop(paste0("The length of 'vax_names_short' (",length(vax_names_short),") should be the same as the length of 'after_vax_lab' (",length(after_vax_lab),") + 1."))
  
  obj$rws$name       <- vax_names_short
  obj$rws$t0         <- t0  # obj$data_parameters$vax_time
  obj$rws$cond_var   <- cond_vax_var
  obj$rws$cond_value <- cond_vax_values
  
  # vaccine dependent variables. Adding "before"- or "after"-  risk windows 'lab'
  if("vax_dep" %in% names(extra_pars)) {
    obj$rws$vax_dep <- extra_pars$vax_dep
    extra_pars <- extra_pars[names(extra_pars)!="vax_dep"]  
    if(is.null(names(obj$rws$vax_dep))) names(obj$rws$vax_dep) <- rep("before",length(obj$rws$vax_dep))
  }
  
  # cut points name:
  if("cut_points_name" %in% names(extra_pars)) {
    obj$rws$cut_points_name <- paste0("cut_points_",extra_pars$cut_points_name)
    extra_pars <- extra_pars[names(extra_pars)!="cut_points_name"]  
  }
  else obj$rws$cut_points_name <- ""
  
  
  
  if(F){  
    if(!missing(obj)){
      if("rws" %in% names(extra_pars)){
        if("rws" %in% names(obj)){
          if(!is.list(rws)) stop("'rws' should be a list.")
          for(iname in names(rws)) obj[["rws"]][[iname]] <- rws[[iname]]
          rws <- obj$rws
        } 
        else obj$rws <- rws
      }
      else if("rws" %in% names(obj)) rws <- obj$rws
    }
    else{
      if("rws" %in% names(extra_pars)) if(!is.list(rws)) stop("'rws' should be a list.") 
      scri_obj <- list(rws=rws)  
    }
  }
  
  
  if( length(extra_pars)>0){
    for(ibefore_after in c("before","after")){
      if( any( (cond1<-substring(names(extra_pars),1,nchar(ibefore_after))==ibefore_after) | (cond2<-substring(names(extra_pars),nchar(names(extra_pars))-nchar(ibefore_after)+1,nchar(names(extra_pars)))==ibefore_after) )){
        
        if(any(cond1)) for(iname in names(extra_pars)[cond1]) {
          obj$rws[[ibefore_after]][[(tmp<-trimws(substring(iname,nchar(ibefore_after)+1),"left","[ _.$]"))]]  <- extra_pars[[iname]]
          if(!(tmp %in% names_before_after)) 
            stop("strange name[s] for '", paste0(tmp[tmp %in% names_before_after],collapse="','"), "'. (Possible names in this format: '", paste0( paste0(ibefore_after,"_",names_before_after),collapse="','"),"')" )
        }
        if(any(cond2)) for(iname in names(extra_pars)[cond2]) {
          obj$rws[[ibefore_after]][[(tmp<-trimws(substring(iname,1,nchar(iname)-nchar(ibefore_after)),"right","[ _.]"))]] <- extra_pars[[iname]] 
          if(!(tmp %in% names_before_after)) 
            stop("strange name[s] for '", paste0(tmp[tmp %in% names_before_after],collapse="','"), "'. (Possible names in this format: '", paste0( paste0("_",names_before_after,"_",ibefore_after),collapse="','"),"')" )
        }        
        extra_pars <- extra_pars[ !cond1 & !cond2  ]
      } 
    }
    if("rws_def" %in% names(extra_pars)) {
      obj$rws$rws_def <- rws_def
      extra_pars <- extra_pars[ names(extra_pars)!= "rws_def"  ]
    }
  }
  
  
  #### cut_points:
  
  ## before vax part:
  # cut_points for the part before vax: before$cut_points
  if("cut_points_before" %in% names(extra_pars)) obj$rws$before$cut_points <- extra_pars[[       "cut_points_before"]]
  if("before_cut_points" %in% names(extra_pars)) obj$rws$before$cut_points <- extra_pars[["before_cut_points"       ]]
  if(is.null(obj$rws$before$cut_points)) stop("'cut_points_before' must be defined." )
  if(is.list(obj$rws$before$cut_points) & length(obj$rws$before$cut_points)==1) obj$rws$before$cut_points <- obj$rws$before$cut_points[[1]] 
  
  obj$rws$before$lab <- before_vax_lab
  if(!missing(ref) | is.null(obj$rws$before$ref)) obj$rws$before$ref <- ref[1]
  
  # extra parameters for period before vax:
  if(is.null(obj$rws$before$lab_add_interval)){ if("lab_add_interval_before" %in% names(extra_pars)) obj$rws$before$lab_add_interval <- extra_pars[["lab_add_interval_before"]] else obj$rws$before$lab_add_interval <- T     } 
  if(is.null(obj$rws$before$no_last_interval)){ if("no_last_interval_before" %in% names(extra_pars)) obj$rws$before$no_last_interval <- extra_pars[["no_last_interval_before"]] else obj$rws$before$no_last_interval <- T     } 
  if(is.null(obj$rws$before$overlap_priority)){ if("overlap_priority_before" %in% names(extra_pars)) obj$rws$before$overlap_priority <- extra_pars[["overlap_priority_before"]] else obj$rws$before$overlap_priority <- "next"} 
  extra_pars <- extra_pars[!(names(extra_pars) %in% c("lab_add_interval_before","no_last_interval_before","overlap_priority_before"))]
  
  
  
  ## after vax part:
  # cut_points for the part after vax: after$cut_points
  if(!is.null(obj$rws$after$cut_points)) {
    if(!is.list(obj$rws$after$cut_points)) {
      obj$rws$after$cut_points <- list(obj$rws$after$cut_points)
      names(obj$rws$after$cut_points) <- vax_names_short[1]                        
    }
    if(length(vax_names_short)>2 & length(obj$rws$after$cut_points)==1) {          
      obj$rws$after$cut_points <- rep(obj$rws$after$cut_points,length(vax_names_short)-1)
      names(obj$rws$after$cut_points) <- vax_names_short[-1]
    }  
  }
  # after_vax_lab
  if(is.null(after_vax_lab)) after_vax_lab <- as.list(paste0("dose ",1:(length(vax_names_short)-1)," "))
  else{
    if( !is.list(after_vax_lab) ) {
      if(length(after_vax_lab)==length(vax_names_short)-1) after_vax_lab <- as.list(after_vax_lab)
      else after_vax_lab <- list(after_vax_lab)
    }
    if(length(vax_names_short)>2 & length(after_vax_lab)==1) after_vax_lab <- rep(after_vax_lab,length(vax_names_short)-1)
  }
  obj$rws$after$lab <- after_vax_lab
  names(obj$rws$after$lab) <- vax_names_short[-1]
  if(!missing(ref) | is.null(obj$rws$after$ref)){ 
    if(length(ref)>1) obj$rws$after$ref <- ref[2:length(vax_names_short)]
    else obj$rws$after$ref <- rep(ref,length(vax_names_short)-1)
  }
  
  # extra parameters for period after vax:
  if(is.null(obj$rws$after$lab_add_interval)){ if("lab_add_interval_after" %in% names(extra_pars)) obj$rws$after$lab_add_interval <- extra_pars[["lab_add_interval_after"]] else obj$rws$after$lab_add_interval <- T     } 
  if(is.null(obj$rws$after$no_last_interval)){ if("lab_add_interval_after" %in% names(extra_pars)) obj$rws$after$no_last_interval <- extra_pars[["no_last_interval_after"]] else obj$rws$after$no_last_interval <- F     } 
  if(is.null(obj$rws$after$overlap_priority)){ if("lab_add_interval_after" %in% names(extra_pars)) obj$rws$after$overlap_priority <- extra_pars[["overlap_priority_after"]] else obj$rws$after$overlap_priority <- "next"} 
  extra_pars <- extra_pars[!(names(extra_pars) %in% c("lab_add_interval_after","no_last_interval_after","overlap_priority_after"))]
  
  if(length(vax_names_short)>2){
    if(length(obj$rws$after$lab_add_interval) ==1) obj$rws$after$lab_add_interval <- rep(obj$rws$after$lab_add_interval,length(vax_names_short)-1)
    if(length(obj$rws$after$no_last_interval) ==1) obj$rws$after$no_last_interval <- rep(obj$rws$after$no_last_interval,length(vax_names_short)-1)
    if(length(obj$rws$after$overlap_priority) ==1) obj$rws$after$overlap_priority <- rep(obj$rws$after$overlap_priority,length(vax_names_short)-1)
  }  
  
  
  # before the first vaccine part:
  if(obj$rws$cut_points_name=="")
    rws_text <- paste0(                                                                               ' list( '    ) 
  else rws_text <- paste0(                             obj$rws$cut_points_name,                          ' = list( '  )  
  rws_text <- paste0(rws_text, '\n\t',                 obj$rws$name[1],                                  ' = list( '  )
  rws_text <- paste0(rws_text, 'name="',               obj$rws$name[1],                                  '"'          )
  rws_text <- paste0(rws_text, ', t0="',               obj$rws$t0[1],                                    '"'          )
  rws_text <- paste0(rws_text, ', cond_var="',         obj$rws$cond_var[1],                              '"'          )
  rws_text <- paste0(rws_text, ', cond_value="',       obj$rws$cond_value[1],                            '"'          )
  
  if(!is.null(obj$rws$vax_dep))
    rws_text <- paste0(rws_text, ', vax_dep=c(', paste0(paste0(names(obj$rws$vax_dep),'="',obj$rws$vax_dep,'"'), collapse=','),  ')' )
  
  rws_text <- paste0(rws_text, ', cut_points=c(',      paste0(obj$rws$before$cut_points,collapse=', '),  ')'          )
  rws_text <- paste0(rws_text, ', lab=c("',            paste0(obj$rws$before$lab,collapse='","'),        '")'         )
  rws_text <- paste0(rws_text, ', ref="',              obj$rws$before$ref,                               '"'          )
  rws_text <- paste0(rws_text, ', lab_add_interval=',  obj$rws$before$lab_add_interval                                )
  rws_text <- paste0(rws_text, ', no_last_interval=',  obj$rws$before$no_last_interval                                )
  rws_text <- paste0(rws_text, ', overlap_priority="', obj$rws$before$overlap_priority,                  '" )'        )
  
  # after each dose parts:
  for(ivax in 2:length(obj$rws$name)) {
    rws_text <- paste0(rws_text, ', \n\t',               obj$rws$name[ivax],                                       ' = list( '  )
    rws_text <- paste0(rws_text, 'name="',               obj$rws$name[ivax],                                       '"'          )
    rws_text <- paste0(rws_text, ', t0="',               obj$rws$t0[ivax],                                         '"'          )
    rws_text <- paste0(rws_text, ', cond_var="',         obj$rws$cond_var[ivax],                                   '"'          )
    rws_text <- paste0(rws_text, ', cond_value="',       obj$rws$cond_value[ivax],                                 '"'          )
    
    if(!is.null(obj$rws$vax_dep))
      rws_text <- paste0(rws_text, ', vax_dep=c(', paste0(paste0(names(obj$rws$vax_dep),'="',obj$rws$vax_dep,'"'), collapse=','),  ')'         )
    
    rws_text <- paste0(rws_text, ', cut_points=c(',      paste0(obj$rws$after$cut_points[[ivax-1]],collapse=', '),  ')'         )
    rws_text <- paste0(rws_text, ', lab=c("',            paste0(obj$rws$after$lab[[ivax-1]],collapse='","'),        '")'        )
    rws_text <- paste0(rws_text, ', ref="',              obj$rws$after$ref[ivax-1],                                 '"'         )
    rws_text <- paste0(rws_text, ', lab_add_interval=',  obj$rws$after$lab_add_interval[ivax-1]                                 )
    rws_text <- paste0(rws_text, ', no_last_interval=',  obj$rws$after$no_last_interval[ivax-1]                                 )
    rws_text <- paste0(rws_text, ', overlap_priority="', obj$rws$after$overlap_priority[ivax-1],                     '" )'      )
  }
  
  rws_text <- paste0(rws_text, ' )')
  
  obj$rws$rws_def <- rws_text
  
  # add lab order
  if(lab_order) obj <- create_lab_orders(obj, data=data)
  
  obj
}  # end of function 'define_rws'


#define_rws(vax_def,  cut_points_before = c(-91,-29), cut_points_after = c(0,1,29,62), cond_vax_var="vax_name", data=data_vax)
#define_rws(vax_def,  cut_points_before = c(-91,-29), cut_points_after = c(0,1,29,62), cond_vax_var="vax_name", data=data_vax, vax_dep = c( before="type_vax_short", after="dist_gt_60" ))


create_lab_orders <- function(obj, 
                              data, 
                              vax_dep_order = c("pfizer","pfize","moderna","moder","astrazeneca","astra","az","j&j","jj","janssen","janss")
){
  
  #if(!missing(obj) & any(is.na(rws))) if(!is.null(obj$rws)) 
  rws <- eval(parse(text=obj$rws$rws_def) )
  
  lab_orders <- list(  
    #c("F", "M" ),
    c("sex0","sex:0","sexF","sex:F", "sex1", "sex:1" , "sexM", "sex:M" ),
    paste0("age",c("(-1,30]","(30,40]","(30,50]","(30,60]","(30,Inf]",">=30",">30","(40,50]","(50,60]","(50,65]","(50,Inf]",">=50",">50",">=60","(60,Inf]", ">65" ))#,
    #c("(-1,30]","(30,40]","(30,50]","(30,60]","(30,Inf]",">=30",">30","(40,50]","(50,60]","(50,65]","(50,Inf]",">=50",">50",">=60","(60,Inf]", ">65" ),
    #paste0("d",1:nvax,":")
  )
  
  if(any(!is.na(rws))){
    
    if(!is.list(rws)) rws <- eval(parse(text=rws))
    rws <- rws[sapply(rws,length)>1]
    
    refs_list_name <- unique(sapply(rws, function(x) x$ref))
    if(length(refs_list_name)==0) warning(paste0("There are no reference category!"))
    if(length(refs_list_name) >1) warning(paste0('It should be just one reference category. They are ',length(refs_list_name),'now:"',paste0(refs_list_name,collapse='", "'),'".'))
    
    refs_list_name <- sapply(rws,function(x) ifelse( any(grepl(refs_list_name, x$lab)), x$name, NA) ) 
    refs_list_name <- refs_list_name[!is.na(refs_list_name)]
    
    
    vax_dep_before <- unique(unlist(lapply(rws,function(x)x$vax_dep[names(x$vax_dep)=="before"])))
    vax_dep_after  <- unique(unlist(lapply(rws,function(x)x$vax_dep[names(x$vax_dep)=="after"])))
    
    if(any( sapply(rws,function(x)!is.null(x$vax_dep)) )){  
      
      rws <- lapply(rws,function(x){ if(is.null(names(x$vax_dep))) names(x$vax_dep)<-"before"; x} )
      #rws <- unique(rws$after$ref)[1]   #lapply(rws,function(x){ if(is.null(names(x$vax_dep))) names(x$vax_dep)<-"before"; x} )
      
      
      if(!missing(data)){
        if(!is.list(vax_dep_order)) vax_dep_order <- list(vax_dep_order)
        for(iorder in  1:length(vax_dep_order) )
          for(ivar in  vax_dep_before  ){
            vax_dep_values <- unique(data[,ivar])
            tmp <-  match(gsub(" ","",tolower(vax_dep_order[[iorder]])),gsub(" ","",tolower(vax_dep_values)))
            vax_dep_values <- unique(c( vax_dep_values[tmp[!is.na(tmp)]], vax_dep_values) )
            vax_dep_values <- vax_dep_values[!is.na(vax_dep_values)]
            lab_orders <- c(lab_orders,  list(  c(rbind( vax_dep_values, gsub(" ","",vax_dep_values),
                                                         paste0("no ",vax_dep_values), paste0("no ",gsub(" ","",vax_dep_values))
            ), "no ","no_"   )) )
          }
      }
    }
    
    if(any( sapply(rws,function(x)!is.null(x$cond_value)) )){  
      dose_names <- unlist(lapply(rws,function(x)x$cond_value))
      lab_orders <- c(lab_orders,  list(  dose_names[!duplicated(dose_names)]   ))
    }
    
    if(any( sapply(rws,function(x)!is.null(x$lab)) )){ 
      labs <- unlist(lapply(rws[sapply(rws,function(x) !(x$name %in% refs_list_name) )],function(x)x$lab))
      names(labs) <- NULL
      lab_orders <- c(lab_orders,  list(  labs[!duplicated(labs)]   ))
    }
    
    if(any( sapply(rws,function(x)!is.null(x$name)) )){  
      names <- unlist(lapply(rws,function(x)x$name))
      names(names) <- NULL
      lab_orders <- c(lab_orders,  list(  names[!duplicated(names)]   ))
    }
    
  }
  
  lab_orders <- c( lab_orders, list(
    ##c("Pfi", "Mod", "Ast",  "JJ","J&J" ),
    ##c("Pfizer", "Moderna", "AstraZeneca", "JJ","J&J" ),
    substring( unlist( ifelse(any(ls()=="brands"), list( c(rbind( brands, paste0("no ",brands), paste0("no_",brands)), "no" ) ), "" )) , 1, 5) ,
    unlist( ifelse(any(ls()=="brands"), list( c(rbind( brands, paste0("no ",brands), paste0("no_",brands)), "no" ) ),  "" )) ,
    c("reference","buf", "dose", "vax", "boos", "booster"),
    c("ref",      "buf", "dose", "vax", "boos", "booster"),
    ##c("pre-","buf", "dose", "vax", "boos", "booster"),
    #paste0( rep(paste0("dose ",1:nvax),each=12), rep(c("ref"," ref","buf"," buf","<"," <","("," (","["," [",">"," >"),nvax) ),
    ##paste0( rep(paste0("dose ",1:nvax),each=12), rep(c("pre-"," pre-","buf"," buf","<"," <","("," (","["," [",">"," >"),nvax) ),
    #c("ref","buf", paste0("dose ",1:nvax), "only_for_time_adj" ),
    ##c("pre-","buf", paste0("dose ",1:nvax), "only_for_time_adj" ),
    c("[-91;-29]","[-90;-30]","[-29;-1]","[-28;-1]","[0;0]","[1;7]","[1;14]","[1;21]","[1;28]","[8;14]","[15;21]","[15;28]","[22;28]","[22;42]",">28",
      "[29;60]","[29;61]","[29;63]",">42",">60","[61;180]","[61;182]",">61","[63;180]","[63;182]",">63",">180",">182")   ))
  
  
  
  
  if(any(!is.na(rws)))
    if(length(vax_dep_after)>0)
      if(!missing(data)){
        if(!is.list(vax_dep_order)) vax_dep_order <- list(vax_dep_order)
        for(iorder in  1:length(vax_dep_order) )
          for(ivar in  vax_dep_after ){
            vax_dep_values <- unique(data[,ivar])
            tmp <-  match(gsub(" ","",tolower(vax_dep_order[[iorder]])),gsub(" ","",tolower(vax_dep_values)))
            vax_dep_values <- unique(c( vax_dep_values[tmp[!is.na(tmp)]], vax_dep_values) )
            vax_dep_values <- vax_dep_values[!is.na(vax_dep_values)]
            lab_orders <- c(lab_orders,  list(  vax_dep_values ))
          }
      }
  
  lab_orders <- c(lab_orders, list(
    paste0(":",c("(-1,30]","(30,60]","(60,Inf]","(60, Inf]")),
    paste0(":",c("(-Inf,-1]","(-1,70]","(70,Inf]")),
    paste0(":",c("(-\U221E,-1]","(-Inf,-1]","(-1,21]","(-1,30]","(21,35]","(30,60]","(35,56]","(56,84]","(60,Inf]","(60, Inf]","(84, Inf]","(84, \U221E]")),
    paste0(":",c("(-\U221E,-1]","(-Inf,-1]","(-1,175]","(175,Inf]","(175, \U221E]"))
  ))
  
  obj$lab_orders <- lab_orders
  
  obj
}  # end of function 'create_lab_orders'  


print.scri_parameters <- function(x, show_order=F){  
  
  cat("\n*****  SCRI parameters  *****\n")
  
  print_scri_help_format <- function( x, i, text, width=40){ 
    if(is.null(names(x[[i]]))) {
      main_part <- paste0(x[i],collapse="','")
      if(length(x[[i]])==1 & is.character(x[[i]])) main_part <- paste0('"', main_part, '"')
    }
    else  main_part <- paste0('c( ', paste0( paste0(names(x[[i]]), '="',x[[i]],'"'), collapse=','), ' )' )
    cat(paste0(format( text, width=width ),
               format(names(x)[i],width=max(nchar(names(x))), justify = "right")," = ",
               main_part, "\n"))
  }
  print_scri_help_noformat <- function( x, i, text, width=40){
    cat(paste0(format( text, width=width ),
               names(x)[i]," = ",
               ifelse(length(x[[i]])==1 & is.character(x[[i]]), "'",'') ,
               paste0(x[i],collapse="','"),
               ifelse(length(x[[i]])==1 & is.character(x[[i]]), "'", '' ) , "\n"))
  }
  
  xx <- unclass(x)
  # nvax
  if(any(names(xx)=="nvax"))
    cat(paste0("\nNumber of doses: nvax = ",xx$nvax,"\n"))
  
  # dataset parameters
  if(any(names(xx)=="data_parameters")){
    cat("\nThe data parameters:\n\n")
    variable_word <- " variable"
    for(ipar in 1:length(xx$data_parameters)){
      if(names(xx$data_parameters)[ipar]=="vax1") {cat("\n");variable_word=" value"}
      if(names(xx$data_parameters)[ipar]=="id") cat("\n")
      print_scri_help_format(xx$data_parameters, ipar, text=paste0(" - the ",names(xx$data_parameters)[ipar],variable_word,": ") )
      if(names(xx$data_parameters)[ipar]=="vax1") variable_word=" variable"
    }
    cat("\n")
  }
  # risk intervals parameters
  if(any(names(xx)=="rws")){
    cat("The risk intervals parameters:\n\n")
    for(irw in 1:length(xx$rws)){
      if(names(xx$rws)[irw]=="rws_def"){
        cat("\n 1. The risk window parameters: \n\n" )
        print_scri_help_noformat(xx$rws, irw, text=paste0("   a) the risk window definition:\n\n      ") )
      }
      
      # 'before'  and 'after' parameters
      if(names(xx$rws)[irw] %in% c("before")){
        cat("\n   b) windows before the first dose:\n\t")
        if(xx$rws$before$no_last_interval)
          cat(paste0(' risk windows = { ' , paste0( paste0('"',xx$rws$before$lab,'" = [',xx$rws$before$cut_points[-length(xx$rws$before$cut_points)],";",xx$rws$before$cut_points[-1]-1,"]"), collapse=", ")," }\n\t "))
        other_before_par <- list(lapply(xx$rws$before[!(names(xx$rws$before) %in% c("cut_points","lab"))],function(x,i)x[[i]],1))
        print_scri_help_noformat(other_before_par, 1, text="other_param ",width=0 )
        cat("\n")
      }
      if(names(xx$rws)[irw] %in% "after"){
        cat("   c) windows after vaccination:\n")
        for(ivax in 1:length(xx$rws$after$cut_points)){
          cat("       after dose ",names(xx$rws$after$cut_points)[ivax],":\t")
          if(xx$rws$after$no_last_interval[ivax])
            cat(paste0("risk windows = { [" , paste0( paste0(xx$rws$after$cut_points[[ivax]][-length(xx$rws$after$cut_points[[ivax]])],";",
                                                             xx$rws$after$cut_points[[ivax]][-1]-1), collapse="], ["),"] }; "))
          else 
            cat(paste0("risk windows = { [" , paste0( paste0(xx$rws$after$cut_points[[ivax]],";",c(xx$rws$after$cut_points[[ivax]][-1]-1,"max")), collapse="], ["),"] }; "))
          other_after_par <- list(sapply(xx$rws$after[names(xx$rws$after)!="cut_points"],function(x,i)x[[i]],ivax))
          print_scri_help_noformat(other_after_par, 1, text="other_param",width=11 )
        }
        cat("\n")
      }
      
      #if(names(xx)[irw] %in% "lab_orders" & show_order)
      #  for(i1 in 1:length(xx$rws$lab_orders)){
      #    print_scri_help_format(xx$rws$lab_orders, i1, text=paste0("      the ",names(xx$rws$lab_orders)[i1]," variable: ") )
      #  }  
      
      # other risk windows parameters
      if( !( names(xx$rws)[irw] %in% c("rws_def","before","after") )) 
        print_scri_help_format(xx$rws, irw, text=paste0("   - the ",names(xx$rws)[irw]," variable: ") )
    }
    cat("\n")
  }
  
  # lab_orders
  if("lab_orders" %in% names(xx) & show_order ){
    cat("\n\nParameter 'lab_orders:'\n")
    for(i1 in 1:length(xx$lab_orders)){
      print_scri_help_noformat(xx$lab_orders, i1, text=paste0("   ",names(xx$lab_orders)[i1]), width=0  )
    }  
  }
  
  # other parameters
  if(any(!(names(xx) %in% c("data_parameters","rws","lab_orders"))))
    for(ipar in 1:length(tmp<-xx[ !(names(xx) %in% c("data_parameters","rws","lab_orders")) ]) )
      print_scri_help_format(tmp, ipar, text=paste0("The ",names(tmp)[ipar]," variable: ") )
  invisible(x)
}   # end of function 'print.scri_input'


print.scri_tabs <- function(x,digits){  
  attributes(x)[names(attributes(x))!="names"] <- NULL
  if(missing(digits)) print(x)
  else lapply(1:length(x),function(i){ print(names(x)[i]);print(format(x[[i]],digits=digits, nsmall=digits)) })
  NULL
}  # end print.scri.tabs


###############################
###############################


# 
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
#table1 <- function(x, digits=2){
#  x <- as.factor(x)
#  cbind( n=(tb<-table(x,useNA="ifany")), 
#         percent=round(100*tb/sum(tb),digits), 
#         percent2=c(round(100*(tb2<-table(x))/sum(tb2),digits),rep(NA,length(tb)-length(tb2))) )[,c(T,T,any(is.na(x)))]
#}




###################################################
#
#  baseline tables
# 
characteristics <- function(data, event, path, condition_value="", data_source_var="datasource",
                            vax_name="vax_number", id="pat_n", start_obs="study_entry_days",
                            vax_time="vax_days", vax_date="vax_date", vax_dist="dist", vax_dist_60="dist_gt_60",
                            death_date="death_date", death_days="death_days", # death_date="date_of_death", 
                            age="age_at_study_entry", age30_50="age30_50", age30="age30",
                            sex="sex", sexc="sexc", sex_age30="sex_age30", 
                            subpopulations=c(0:6), #c(0,10) ==> all subpopulations
                            vax_part=F, event_vax_part=T,
                            lparal=T, n_cores=NA, lprint=F, flowchart_print=T, flowchart_type=c("table","all")
){  
  cat(" flowchart ")
  sys_time <- Sys.time()
  
  extra_attributes <- c( vax_name=vax_name, id=id, start_obs=start_obs, 
                         vax_time=vax_time, vax_date=vax_date, vax_dist=vax_dist, vax_dist_60=vax_dist_60,
                         death_date=death_date, death_days=death_days, 
                         age=age, age30_50=age30_50,  age30=age30,     
                         sex=sex, sexc=sexc, sex_age30=sex_age30,
                         condition_value=condition_value  )
  if(!missing(event)) extra_attributes <- c( event=event, extra_attributes )
  
  if(       id!="id"       ) names(data)[names(data)==id       ] <- "id"
  if(start_obs!="start_obs") names(data)[names(data)==start_obs] <- "start_obs"
  
  if(vax_name!="vax_name") names(data)[names(data)==vax_name] <- "vax_name"
  if(vax_time!="vax_time") names(data)[names(data)==vax_time] <- "vax_time"
  if(vax_date!="vax_date") names(data)[names(data)==vax_date] <- "vax_date"
  
  if(any(names(data)==vax_dist)) {
    if(vax_dist!="vax_dist") names(data)[names(data)==vax_dist] <- "vax_dist"
  }  else data$vax_dist <- rep(NA,nrow(data))
  
  if(any(names(data)==vax_dist_60)){
    if(vax_dist!="vax_dist_60") names(data)[names(data)==vax_dist_60] <- "vax_dist_60"
  } else data$vax_dist_60 <- rep(NA,nrow(data))
  
  
  if(event_vax_part){
    names(data)[names(data)==paste0(event,"_days")] <- "event_days"
    names(data)[names(data)==paste0(event,"_date")] <- "event_date"
    
    data$event  <- as.numeric(!is.na(data$event_days))
    data$eventc <- paste0("event:",data$event)
  }
  
  names(data)[names(data)==death_days] <- "death_days"
  names(data)[names(data)==death_date] <- "death_date"
  
  # create some strata variables:
  
  if(age     !="age"       ) names(data)[names(data)==age      ] <- "age"
  if(age30_50!="age30_50"  ) names(data)[names(data)==age30_50 ] <- "age30_50"
  if(age30   !="age30"     ) names(data)[names(data)==age30    ] <- "age30"
  
  if(sex      !="sex"      ) names(data)[names(data)==sex      ] <- "sex"
  if(sexc     !="sexc"     ) names(data)[names(data)==sexc     ] <- "sexc"
  if(sex_age30!="sex_age30") names(data)[names(data)==sex_age30] <- "sex_age30"
  
  # delete some variables:
  names_to_delete <- c(id,vax_name, vax_time, vax_date, vax_dist, vax_dist_60, death_days, death_date, age, sex, age30_50, age30, sexc, sex_age30 )
  if(event_vax_part) names_to_delete <- c(names_to_delete, event, paste0(event,"_days"), paste0(event,"_date"))
  good_names <- c("id",  "start_obs",
                  "vax_name","vax_brand", "vax_time", "vax_date", "vax_n", "vax_days_v1", "vax_dist","vax_dist_60",
                  "death_days","death_date","age","age30_50", "age30", "sexc", "sex_age30")
  if(event_vax_part)  good_names <- c(good_names, "event", "event_days", "event_date", "eventc")
  names_to_delete <- names_to_delete[ !(names_to_delete %in% good_names) ]
  data[,names_to_delete] <- NULL
  
  data0 <- data[ !duplicated(data[,c("id","vax_time")]), good_names ]
  gc()
  
  summary_with_NA <- function(x)c(mean=mean(x,na.rm=T),quantile(x,c(0,0.25,0.5,0.75,1),na.rm=T),NAs=sum(is.na(x)),n=length(x))
  
  characteristics_kern <- function( istrata_tab, strata_info, data, vax_part, event_vax_part, lprint=lprint ) { 
    
    if(lprint) print(strata_info[istrata_tab,])
    i <- strata_info[istrata_tab, "dataset"]
    
    if(!event_vax_part & i>0) return("only for vaccines")  
    
    days_1jan2020 <- as.numeric(difftime(as.Date("2020-08-31"),as.Date("2020-01-01"),units="days"))
    
    ####
    # all id's:
    if(i==0){ 
      data        <- data0
      data_deaths <- data0[!is.na(data0$death_days),]
    }
    ####
    # only for vaccinated with event:
    if(i==1) {
      data        <- data0[data0$event==1  & !is.na(data0$vax_date),]
      data_deaths <- data0[!is.na(data0$death_days) ,]
    }
    
    ####
    # only with event: from (vax1-90d) and observed: at vax1-90-365 or  1jan2020 or age<=2
    if(i==2) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365,days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1-90 <= data0$event_days , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365, days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1-90 <= data0$death_days, ]
    }
    
    ####################
    # event: from vax1
    ####
    # only with event: from vax1 and observed: at vax1
    if(i==3) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$event_days , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$death_days , ]
    }
    ####
    # only with event: from vax1 and observed: at vax1-90
    if(i==4) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$event_days , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$death_days , ]
    }
    ####
    # only with event: from vax1 and observed: at 1sep2020
    if(i==5) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             (data0$start_obs <= 0  | (data0$age<=1 & data0$start_obs<=data0$vax_days_v1-90) )& 
                             data0$vax_days_v1 <= data0$event_days , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             (data0$start_obs <= 0 | (data0$age<=1 & data0$start_obs<=data0$vax_days_v1-90) ) &
                             data0$vax_days_v1 <= data0$death_days , ]
    }
    ####
    # only with event: from vax1 and observed: at vax1-90-365 or  1jan2020 or age<=2
    if(i==6) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365,days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$event_days, ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365,days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$death_days , ]
    }
    ####
    # only with event: from vax1 and observed: at 1jan2020
    if(i==7) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             (data0$start_obs <= days_1jan2020 | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$event_days, ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             (data0$start_obs <= days_1jan2020 | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$death_days , ]
    }
    
    
    ####################
    # event: in [vax; vax+28d]
    
    # only with event: in [vax1; vax1+28d] and observed: at vax1
    if(i==8) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$event_days &
                             data0$event_days <= data0$vax_time + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$death_days &
                             data0$death_days <= data0$vax_time + 28 , ]
    }
    
    
    # only with event: in [vax1; vax1+28d] and observed: at vax1-90
    if(i==9) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$event_days &
                             data0$event_days <= data0$vax_time + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$death_days &
                             data0$death_days <= data0$vax_time + 28 , ]
    }
    
    
    ####
    # only with event: in [vax1; vax1+28d] and observed: at vax1-90-365 or  1jan2020 or age<=2
    if(i==10) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365,days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$event_days &
                             data0$event_days <= data0$vax_time + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             ( data0$start_obs <= pmax(data0$vax_days_v1-90-365,days_1jan2020) | (data0$age<=2 & data0$start_obs<=data0$vax_days_v1-90) ) & 
                             data0$vax_days_v1 <= data0$death_days &
                             data0$death_days <= data0$vax_time + 28 , ]
    }
    
    nrow_data        <- nrow(data)
    nrow_data_deaths <- nrow(data_deaths)
    
    
    if(nrow_data>0 | nrow_data_deaths>0){
      
      istrata_var <-  as.character(strata_info[istrata_tab,"strata"])
      
      if(istrata_var =="no_strata"){
        if(nrow_data>0       ) data$strata_variable        <- "no_strata"
        if(nrow_data_deaths>0) data_deaths$strata_variable <- "no_strata"
      }
      else {
        if(nrow_data>0       ) data$strata_variable        <- data[       ,istrata_var]
        if(nrow_data_deaths>0) data_deaths$strata_variable <- data_deaths[,istrata_var]
      } 
      
      if(nrow_data>0       ) data_strata        <- data[       !is.na(data$strata_variable       ) & !is.na(data$id       ),]
      if(nrow_data_deaths>0) data_deaths_strata <- data_deaths[!is.na(data_deaths$strata_variable) & !is.na(data_deaths$id),]
      
      gc()
      
      flowchart <- list()
      flowchart <- c(flowchart, info = list(strata_info[istrata_tab,]))
      
      if(nrow_data>0){ 
        if(nrow(data_strata)>0){
          
          
          
          flowchart <- c(flowchart, n_ids = list(table1(data_strata[!duplicated(data_strata$id),c("strata_variable","id")][,"strata_variable"] )))
          
          if(vax_part){ 
            flowchart <- c(flowchart, n_ids_per_vax_name_brand = list(table1( data_strata[, c("strata_variable","vax_name","vax_brand")])) )
          }
          
          if(event_vax_part){ 
            flowchart <- c(flowchart, n_ids_per_event_vax_name_brand = list(table1( data_strata[, c("strata_variable","eventc","vax_name","vax_brand")])) )
          }
          
          # age as continuous or integer:
          if(any(!is.na(data_strata$age))){
            
            data_strata_age <- data_strata[!is.na(data_strata$age),]
            
            if(       vax_part) flowchart <- c(flowchart, summary_id_age                      = list(do.call("rbind",with(data_strata_age[!duplicated(data_strata_age[,"id"]), c("strata_variable","id","age")], tapply(age,    strata_variable                           , summary_with_NA) )) ))
            if(       vax_part) flowchart <- c(flowchart, summary_id_age_vax_name             = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age",         "vax_name"            )], tapply(age, paste(strata_variable,       vax_name          ), summary_with_NA) )) ))
            if(       vax_part) flowchart <- c(flowchart, summary_id_age_vax_brand            = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age",                    "vax_brand")], tapply(age, paste(strata_variable,                vax_brand), summary_with_NA) )) ))
            if( event_vax_part) flowchart <- c(flowchart, summary_id_age_event                = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age","eventc"                       )], tapply(age, paste(strata_variable,eventc                   ), summary_with_NA) )) ))
            if( event_vax_part) flowchart <- c(flowchart, summary_id_age_event_vax_name_brand = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age","eventc","vax_name","vax_brand")], tapply(age, paste(strata_variable,eventc,vax_name,vax_brand), summary_with_NA) )) ))
            
            rm(data_strata_age)
          } 
          
          gc()
          
          # summaries for vaccination time and date:
          if(vax_part & any(!is.na(data_strata$vax_time))){
            
            data_strata_vax <- data_strata[!is.na(data_strata$vax_time),]
            
            flowchart <- c( flowchart, summary_id_vax_time = list( do.call("rbind",with(
              data_strata_vax[, c("strata_variable","vax_n","id", "vax_time")], tapply( vax_time, paste0(strata_variable," vax_n:",vax_n), summary_with_NA) )) ))
            
            flowchart <- c(flowchart, summary_id_vax_date = list(do.call("rbind",with(
              data_strata_vax[, c("strata_variable","vax_n","id", "vax_date")], tapply( vax_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
            
          } 
          gc()
          
          # vax_dist as continuous or integer:
          if(any(!is.na(data_strata$vax_dist))){
            
            data_strata_dist <- data_strata[!is.na(data_strata$vax_dist),]
            
            if(       vax_part) flowchart <- c(flowchart, summary_id_dist                      = list(do.call("rbind",with(data_strata_dist[!duplicated(data_strata_dist[,"id"]), c("strata_variable","id","vax_dist")], tapply(vax_dist,    strata_variable                          , summary_with_NA) )) ))
            if(       vax_part) flowchart <- c(flowchart, summary_id_dist_vax_name             = list(do.call("rbind",with(data_strata_dist[, c("strata_variable","id","vax_dist",         "vax_name"            )], tapply(vax_dist, paste(strata_variable,       vax_name          ), summary_with_NA) )) ))
            if(       vax_part) flowchart <- c(flowchart, summary_id_dist_vax_brand            = list(do.call("rbind",with(data_strata_dist[, c("strata_variable","id","vax_dist",                    "vax_brand")], tapply(vax_dist, paste(strata_variable,                vax_brand), summary_with_NA) )) ))
            if( event_vax_part) flowchart <- c(flowchart, summary_id_dist_event                = list(do.call("rbind",with(data_strata_dist[, c("strata_variable","id","vax_dist","eventc"                       )], tapply(vax_dist, paste(strata_variable,eventc                   ), summary_with_NA) )) ))
            if( event_vax_part) flowchart <- c(flowchart, summary_id_dist_event_vax_name_brand = list(do.call("rbind",with(data_strata_dist[, c("strata_variable","id","vax_dist","eventc","vax_name","vax_brand")], tapply(vax_dist, paste(strata_variable,eventc,vax_name,vax_brand), summary_with_NA) )) ))
            
            rm(data_strata_dist)
          } 
          
          gc()
          
          ######################
          # event_days:
          if( event_vax_part & any(!is.na(data_strata$event_days))){
            
            data_strata_vax <- data_strata[ !is.na(data_strata$event_days) & !is.na(data_strata$vax_time),]
            
            if( nrow(data_strata_vax)>0){
              if( any(( cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)<0 )) ){
                flowchart <- c(  flowchart, summary_id_event_min_vax_before_vax = list( do.call("rbind",with(
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time","vax_name")]), tapply( (event_days-vax_time), paste(strata_variable,vax_name), summary_with_NA) )) ))
              }
              if( any(( cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)>=0 )) ){
                flowchart <- c(  flowchart, summary_id_event_min_vax_after_vax = list( do.call("rbind",with(
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time","vax_name")]), tapply( (event_days-vax_time), paste(strata_variable,vax_name), summary_with_NA) )) ))
              }
            }
            
            flowchart <- c( flowchart, summary_id_event_time_vax_n = list( do.call("rbind",with(
              data_strata[!is.na(data_strata$event_days), c("strata_variable","vax_n","id", "event_days")], tapply( event_days, paste0(strata_variable," vax_n:",vax_n), summary_with_NA) )) ))
            
            flowchart <- c(flowchart, summary_id_event_date_vax_n = list(do.call("rbind",with(
              data_strata[!is.na(data_strata$event_days), c("strata_variable","vax_n","id", "event_date")], tapply( event_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
            
            if( nrow(data_strata_vax)){
              if( any( (cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)<0 )) ){
                flowchart <- c(flowchart, summary_id_event_min_vax_before_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time", "eventc","vax_name","vax_brand")]), tapply( (event_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), summary_with_NA ) ))  ) )
              }
              if( any( (cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)>=0 )) ){
                flowchart <- c(flowchart, summary_id_event_min_vax_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time", "eventc","vax_name","vax_brand")]), tapply( (event_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), summary_with_NA ) ))  ) )
              }
            }
            rm(data_strata_vax); rm(data_strata)
          } 
        } # if nrow(data_strata) > 0
      }  # end if nrow_data>0
      
      gc()
      
      ######################
      # deaths:
      if(nrow_data_deaths>0){
        
        if(vax_part & nrow(data_deaths_strata)>0){
          
          flowchart$deaths <- c( summary_id_death_date = list(do.call("rbind",with(
            data_deaths_strata[, c("strata_variable","vax_n","id", "death_date")], 
            tapply( death_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
          
          flowchart$deaths <- c( flowchart$deaths, summary_id_death_time = list( do.call("rbind",with(
            data_deaths_strata[, c("strata_variable","vax_n","id", "death_days")], 
            tapply( death_days, paste0(strata_variable," vax_n:",vax_n), summary_with_NA) )) ))
          
          data_deaths_strata <- data_deaths_strata[!is.na(data_deaths_strata$vax_time),]
          
          if(nrow(data_deaths_strata)>0){
            flowchart$deaths <- c( flowchart$deaths, summary_id_death_after_vax = list( do.call("rbind",with(
              data_deaths_strata[, c("strata_variable","id", "death_days","vax_time","vax_name")], 
              tapply( (death_days-vax_time), paste(strata_variable,vax_name), summary_with_NA) )) ))
            
            flowchart$deaths <- c(flowchart$deaths, summary_id_death_after_vax_per_vax_name_brand = list( do.call("rbind",with( 
              data_deaths_strata[, c("strata_variable","id", "death_days","vax_time", "vax_name","vax_brand")], 
              tapply( (death_days-vax_time), paste(strata_variable,vax_name,vax_brand), summary_with_NA ) ))  ) )
          }
          
        } # end if vax_part & nrow(data_deaths_strata)>0
        
        if(event_vax_part & nrow(data_deaths_strata)>0){
          flowchart$deaths <- c(flowchart$deaths, summary_id_death_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
            data_deaths_strata[, c("strata_variable","id", "death_days","vax_time", "eventc","vax_name","vax_brand")], 
            tapply( (death_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), summary_with_NA ) ))  ) )
        }
      } # end nrow_data_deaths>0
      
      
      
    }  # end if nrow(data)>0 | nrow(data_deaths)>0
    else return(NULL)
    
    flowchart
  } # end func "characteristics_kern"
  
  
  
  strata_tab <- expand.grid( dataset=subpopulations, strata=c("no_strata", "age30_50","age30","sexc","sex_age30", "vax_dist_60"))
  #strata_tab <- expand.grid( dataset=c(0:6), strata=c("no_strata", "age30_50","age30","sexc","sex_age30", "vax_dist_60"))
  strata_tab$output_number <- 1:nrow(strata_tab)
  
  strata_tab$selection_name <- c("0" = "all_data",
                                 
                                 "1" = "with_vax_event",
                                 
                                 # all with events from start of the control period:
                                 "2 "= "events_after_vax1-90_and_observed_at_vax1-90-365_or_1jan2020_or_baby",
                                 
                                 # all with events after vax1:
                                 "3" = "events_after_vax1_and_observed_at_vax1",
                                 "4" = "events_after_vax1_and_observed_at_vax1-90",
                                 "5" = "events_after_vax1_and_observed_at_1sept2020_or_baby",
                                 "6" = "events_after_vax1_and_observed_at_vax1-90-365_or_1jan2020_or_baby",
                                 "7" = "events_after_vax1_and_observed_at_1jan2020_or_baby",
                                 
                                 # all with events during 28 days after vax1:
                                 "8" = "events_during_28_days_after_vax1_and_observed_at_vax1",
                                 "9" = "events_during_28_days_after_vax1_and_observed_at_vax1-90",
                                 "10"= "events_during_28_days_after_vax1_and_observed_at_vax1-90-365_or_1jan2020_or_baby"
  )[ strata_tab$dataset + 1]
  strata_tab <- strata_tab[,c("dataset","selection_name","strata","output_number")] 
  
  if(lparal){
    library(parallel)
    if(is.na(n_cores)) n_cores <- detectCores() - 2 
    n_cores <- min( nrow(strata_tab), n_cores, na.rm=T )
    cl      <- makeCluster( n_cores ) 
    clusterExport(cl, c("table1","data0","strata_tab","vax_part","event_vax_part","characteristics_kern","lprint"), envir = environment() )
    #sink(file=paste0(path,"characteristics_log.txt"))
    #print(paste("Total:",nrow(strata_tab),"rows"))
    res <- parLapply(cl, 1:nrow(strata_tab), function(istrata_tab)
      characteristics_kern( istrata_tab, strata_info=strata_tab, data=data0, vax_part=vax_part, event_vax_part=event_vax_part, lprint=lprint ) )
    #sink()
    stopCluster(cl)
    
    #names_overlap <- match(names(res),names(res_paral))
    #res[ names(res_paral)[ names_overlap[!is.na(names_overlap)] ] ] <- res_paral[!is.na(match(names(res_paral),names(res)))]
    
    
  }
  else{
    res <- vector("list",length=nrow(strata_tab))
    for(istrata_tab in 1:nrow(strata_tab)){
      res[[istrata_tab]] <- characteristics_kern( istrata_tab, strata_info=strata_tab, data=data0, vax_part=vax_part, event_vax_part=event_vax_part, lprint=lprint )
      cat(paste0(istrata_tab,"\b"))
    }
  }
  
  names(res) <- strata_tab$strata
  
  data_sel <- sapply(res,function(x)ifelse(any(names(x[[1]])=="selection_name"),paste0(x[[1]][c("dataset","selection_name")],collapse = ": "),""))
  strata_tab <- strata_tab[data_sel!="",]
  data_sel_unique <- data_sel[data_sel!="" & !duplicated(data_sel)]
  res_per_selection <- vector("list", length(data_sel_unique))
  names(res_per_selection) <- data_sel_unique
  for(isel in data_sel_unique)
    res_per_selection[[isel]] <- res[data_sel==isel]
  
  flowchart <- list(data_selection   = res_per_selection, 
                    all_descriptions = strata_tab         )
  attributes(flowchart) <- c( attributes(flowchart), variables = list(extra_attributes) )
  
  if(!missing(path)) {
    if(substring(path,nchar(path),nchar(path))!="/") path <- paste0(path,"/")
    
    data_source <- ""
    if(data_source_var!="") {if(data_source_var %in% names(data)) data_source <- as.character(data[1,data_source_var])}
    if(data_source=="") path <- paste0(path,"characteristics")
    else                path <- paste0(path,data_source,"_characteristics")
    
    if(event_vax_part)      path <- paste0(path,"_",event,"_",vax_name)
    else if(vax_part)       path <- paste0(path,"_",vax_name   )
    if(condition_value!="") path <- paste0(path,"_",condition_value)
    
    save( flowchart, file=paste0(path,".RData") )
  }
  
  if(flowchart_print) 
    for(itype in flowchart_type) characteristics_print(flowchart, path=path, type=itype)
  
  cat(paste0(": duration = ",format(difftime(Sys.time(),sys_time))," (till ",Sys.time(),")\n"))  
  
  invisible(flowchart)
  
}  # end of function 'characteristics'


if(F){
  Sys.time()
  tmp_paral <- characteristics(data=data_vax, vax_name="vax_number", vax_part=T, event_vax_part=F, path=sdr0, id="pat_n", condition_value="", age="age_at_study_entry", lparal=T )
  Sys.time()  
  
  
  Sys.time()
  tmp <- characteristics(data=data_vax, vax_name="vax_number", vax_part=T, event_vax_part=F, path=sdr0, id="pat_n", condition_value="", age="age_at_study_entry", lparal=F )
  Sys.time()  
  
  
  Sys.time()
  tmp2<-characteristics(data=data_vax, vax_name="vax_number", vax_part=F, event_vax_part=T, event=iae, path=sdr, id="pat_n", condition_value="", age="age_at_study_entry", lparal=T )
  Sys.time()
  
  
  sink("tmp"); old<-options(max.print = 99999); print(tmp);options(old);sink()
  
  sink("tmp_paral"); old<-options(max.print = 99999); print(tmp_paral);options(old);sink()
  
  
  sink("tmp2"); old<-options(max.print = 99999); print(tmp_paral);options(old);sink()
}


# characteristics_print(tmp_paral)    

characteristics_print <- function(flowchart, path,
                                  type="all"    # "all", "table"
){
  flowchart_all <- flowchart
  
  if(!missing(path)) sink(paste0(path,".txt"))
  
  for(iname in names(attributes(flowchart_all)$variables))
    cat( paste0('\n',iname,'\t',ifelse(nchar(iname)<=5,'\t',''),'\t=\t"',  attributes(flowchart_all)$variables[iname],'"') )
  cat("\n\n")
  
  if(any(names(attributes(flowchart_all)$variables)=="event")) event <- attributes(flowchart_all)$variables["event"] 
  
  for(i in 0:(length(flowchart_all[[1]])-1) ){
    
    flowchart <- flowchart_all[[1]][[i+1]]
    
    cat("\n\n\n\n*************************************************************************************\n")
    cat("*************************************************************************************\n")
    cat("*************************************************************************************\n\n")
    
    ####
    # all id's:
    if(flowchart_all$all_descriptions[i+1,"dataset"]==0) cat("\t\t\tThe whole dataset\n\n")
    ####
    # only for vaccinated with event:
    if(flowchart_all$all_descriptions[i+1,"dataset"]==1) cat("\t\t\tvaccinated persons with ",event,"\n\n")
    ############### 
    # only with event: from vax1-90 days observed: at (vax_day1 - 90 - 365 days) or at 1jan2020 or age<=2y 
    if(flowchart_all$all_descriptions[i+1,"dataset"]==2) cat("\t\t\tvaccinated persons observed at (first dose - 90 - 365 days) or 1jan2020 or age<=2y with ",event," after (vax_day1 - 90 days)\n\n")
    ############### 
    # only with event: from vax1 and observed: at vax1
    if(flowchart_all$all_descriptions[i+1,"dataset"]==3) cat("\t\t\tvaccinated persons observed at first dose with ",event," after first dose\n\n") 
    ############### 
    # only with event: from vax1 and observed: at vax1-90
    if(flowchart_all$all_descriptions[i+1,"dataset"]==4) cat("\t\t\tvaccinated persons observed at (first dose-90 days) with ",event," after first dose\n\n")
    ############### 
    # only with event: from vax1 and observed: at 1sept2020 or age<=1y 
    if(flowchart_all$all_descriptions[i+1,"dataset"]==5) cat("\t\t\tvaccinated persons observed at the 1sept2020 or age<=1y with ",event," after first dose\n\n")
    ############### 
    # only with event: from vax1 and observed: at (vax_day1 - 90 - 365 days) or at 1jan2020 or age<=2y
    if(flowchart_all$all_descriptions[i+1,"dataset"]==6) cat("\t\t\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose\n\n")
    ############### 
    # only with event: from vax1 and observed: at 1jan2020 or age<=2y
    if(flowchart_all$all_descriptions[i+1,"dataset"]==7) cat("\t\t\tvaccinated persons observed at 1jan2020 or age<=2y with ",event," after first dose\n\n")
    ############### 
    # only with event: in [vax_date; vax_date + 28 days] and observed: after the first dos
    if(flowchart_all$all_descriptions[i+1,"dataset"]==8) cat("\t\t\tvaccinated persons observed at first dose with ",event," in [vax_date; vax_date + 28 days]\n\n")
    ############### 
    # only with event: in [vax_date; vax_date + 28 days] and observed:  after (first dose - 90 days)
    if(flowchart_all$all_descriptions[i+1,"dataset"]==9) cat("\t\t\tvaccinated persons observed at (first dose-90 days) with ",event," in [vax_date; vax_date + 28 days]\n\n")
    ############### 
    # only with event: in [vax_date; vax_date + 28 days] and observed: at (vax_day1 - 90 - 365 days) or at 1jan2020 or age<=2y 
    if(flowchart_all$all_descriptions[i+1,"dataset"]==10) cat("\t\t\tvaccinated persons observed at (first dose - 90 - 365 days) or 1jan2020 or age<=2y with ",event," in [vax_date; vax_date + 28 days]\n\n")
    
    cat("\n\n")
    print(flowchart_all$all_descriptions[i+1,])
    
    
    for(istrata_var in names(flowchart) ){ #c("no_strata", "age30_50","age30","sexc","sex_age30")){
      
      flowchart_strata <- flowchart[[istrata_var]]
      cat("\n\n*******************************************************\n")
      cat("*******************************************************\n")
      cat("*******************************************************\n")
      
      if(istrata_var =="no_strata")
        cat(paste0("\n\tSummary not stratified:\n\n\n"))
      else cat(paste0("\n\tSummary stratified by variable '",istrata_var,"':\n\n\n"))
      
      ####
      if(flowchart_all$all_descriptions[i+1,"dataset"]==0) cat(paste0(" ***\tall persons in the dataset:***\n\nthe number of persons in the dataset: \n" ))
      
      ####
      if(flowchart_all$all_descriptions[i+1,"dataset"]==1) cat(paste0("***\tonly for vaccinated persons with ",event,":***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==2) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or 1jan2020 or age<=2y with ",event," after (vax_day1 - 90 days):***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==3) cat(paste0("***\tvaccinated persons observed at first dose with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==4) cat(paste0("***\tvaccinated persons observed at (first dose-90 days) with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==5) cat(paste0("***\tvaccinated persons observed at the 1sept2020 or age<=1y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==6) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==7) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==8) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==9) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      ############### 
      if(flowchart_all$all_descriptions[i+1,"dataset"]==10) cat(paste0("***\tvaccinated persons observed at (first dose - 90 - 365 days) or at 1jan2020 or age<=2y with ",event," after first dose:***\n\nthe number of persons: \n" ))
      
      print(flowchart_strata$n_ids)
      
      if(!is.null(flowchart_strata$n_ids_per_vax_name_brand)){ 
        cat(paste0("\nthe numbers of persons  per vaccine name and brand:\n"))
        print(flowchart_strata$n_ids_per_vax_name_brand)
      }
      
      if(!is.null(flowchart_strata$n_ids_per_event_vax_name_brand)){ 
        cat(paste0("\nthe numbers for persons  per ",event,", vaccine name and brand:\n"))
        print(flowchart_strata$n_ids_per_event_vax_name_brand)
      }
      
      # age as continuous or integer:
      if(!is.null(flowchart_strata$summary_id_age)){
        
        cat(paste0("\n\nthe distribution of variable '",attributes(flowchart_all)$variables["age"] ,"':\n\n"))
        
        if(!is.null(flowchart_strata$summary_id_age                     )) print( flowchart_strata$summary_id_age                      ); cat("\n")
        if(!is.null(flowchart_strata$summary_id_age_vax_name            )) print( flowchart_strata$summary_id_age_vax_name             ); cat("\n")
        if(!is.null(flowchart_strata$summary_id_age_vax_brand           )) print( flowchart_strata$summary_id_age_vax_brand            ); cat("\n")
        if(!is.null(flowchart_strata$summary_id_age_event               )) print( flowchart_strata$summary_id_age_event                ); cat("\n")
        if(!is.null(flowchart_strata$summary_id_age_event_vax_name_brand)) print( flowchart_strata$summary_id_age_event_vax_name_brand )  
      }
      
      
      if(any(names(flowchart_strata) %in% c("summary_id_vax_time","summary_id_vax_date","summary_id_dist","summary_id_dist_vax_name","summary_id_dist_vax_brand","summary_id_dist_event","summary_id_dist_event_vax_name_brand"))){ 
        
        if(!is.null(flowchart_strata$summary_id_vax_time)){ 
          cat(paste0("\n\nthe distribution of the vaccination variable: '",attributes(flowchart_all)$variables["vax_time"],"':\n"))
          print(flowchart_strata$summary_id_vax_time)
        }
        if(!is.null(flowchart_strata$summary_id_vax_date)){ 
          cat(paste0("\nthe distribution of the vaccination variable: '",attributes(flowchart_all)$variables["vax_date"],"':\n"))
          print(flowchart_strata$summary_id_vax_date)
        }
        
        if(!is.null(flowchart_strata$summary_id_dist                     )){ cat("\nthe distribution of the distance between vaccines:\n"); print(flowchart_strata$summary_id_dist                     ) }                                    
        if(!is.null(flowchart_strata$summary_id_dist_vax_name            )){ cat("\nthe distribution of the distance between vaccines:\n"); print(flowchart_strata$summary_id_dist_vax_name            ) }                          
        if(!is.null(flowchart_strata$summary_id_dist_vax_brand           )){ cat("\nthe distribution of the distance between vaccines:\n"); print(flowchart_strata$summary_id_dist_vax_brand           ) }                         
        if(!is.null(flowchart_strata$summary_id_dist_event               )){ cat("\nthe distribution of the distance between vaccines:\n"); print(flowchart_strata$summary_id_dist_event               ) }                            
        if(!is.null(flowchart_strata$summary_id_dist_event_vax_name_brand)){ cat("\nthe distribution of the distance between vaccines:\n"); print(flowchart_strata$summary_id_dist_event_vax_name_brand) }               
        
      } else cat("\n\nno vaccinated.\n\n")
      
      ######################
      # event_days:
      
      if(any(names(flowchart_strata) %in% c("summary_id_event_min_vax_before_vax","summary_id_event_min_vax_after_vax",
                                            "summary_id_event_time_vax_n", "summary_id_event_date_vax_n",
                                            "summary_id_event_min_vax_before_vax_per_event_vax_name_brand",
                                            "summary_id_event_min_vax_after_vax_per_event_vax_name_brand"    ))){ 
        
        if(!is.null(flowchart_strata$summary_id_event_min_vax_before_vax)){ 
          cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) before vaccination:\n"))
          print(flowchart_strata$summary_id_event_min_vax_before_vax)
        }
        if(!is.null(flowchart_strata$summary_id_event_min_vax_after_vax)){ 
          cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) after vaccination:\n"))
          print(flowchart_strata$summary_id_event_min_vax_after_vax)
        }
        
        if(!is.null(flowchart_strata$summary_id_event_time_vax_n)){ 
          cat(paste0("\n\nthe distribution of the '",event,"_days' variable per vaccine number:\n"))
          print(flowchart_strata$summary_id_event_time_vax_n)
        }           
        if(!is.null(flowchart_strata$summary_id_event_date_vax_n)){ 
          cat(paste0("\n\nthe distribution of the '",event,"_date' variable per vaccine number:\n"))
          print(flowchart_strata$summary_id_event_date_vax_n)
        }           
        if(!is.null(flowchart_strata$summary_id_event_min_vax_before_vax_per_event_vax_name_brand)){ 
          cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) before vaccination per ",event,", vaccine name and brand:\n"))
          print(  flowchart_strata$summary_id_event_min_vax_before_vax_per_event_vax_name_brand)
        }
        if(!is.null(flowchart_strata$summary_id_event_min_vax_after_vax_per_event_vax_name_brand)){ 
          cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) after vaccination per ",event,", vaccine name and brand:\n"))
          print(  flowchart_strata$summary_id_event_min_vax_after_vax_per_event_vax_name_brand)
        }
      } else cat("\n\nno events.\n\n")
      
      
      
      ######################
      # deaths:
      if(!is.null(flowchart_strata$deaths$summary_id_death_date)){ 
        cat(paste0("\n\nthe distribution of the 'death_date' variable:\n"))
        print(flowchart_strata$deaths$summary_id_death_date)
      }            
      if(!is.null(flowchart_strata$deaths$summary_id_death_time)){ 
        cat(paste0("\n\nthe distribution of the 'death_days' variable:\n"))
        print(flowchart_strata$deaths$summary_id_death_time)
      }              
      
      if(!is.null(flowchart_strata$deaths$summary_id_death_after_vax)){ 
        cat(paste0("\n\nthe distribution of 'death_days': days after vaccination:\n"))
        print(flowchart_strata$deaths$summary_id_death_after_vax)
      }                
      if(!is.null(flowchart_strata$deaths$summary_id_death_after_vax_per_vax_name_brand)){ 
        cat(paste0("\nthe distribution of the 'death_days' (days after vaccination) per, vaccine name and brand:\n"))
        print(  flowchart_strata$deaths$summary_id_death_after_vax_per_vax_name_brand)
        
      }
      
      if(!is.null(flowchart_strata$deaths$summary_id_death_after_vax_per_event_vax_name_brand)){ 
        cat(paste0("\nthe distribution of the 'death_days' (days after vaccination) per ",event,", vaccine name and brand:\n"))
        print(  flowchart_strata$deaths$summary_id_death_after_vax_per_event_vax_name_brand)
      } 
      
    }  #end for istrata_var
    
  } # end 'i'
  
  cat("\n\n\nflowchart parts:\n")
  print(flowchart_all$all_descriptions)
  
  cat("\n\nattributes:\n\n")
  print(attributes(flowchart_all))
  
  if(!missing(path)) sink()
  
  invisible(NULL)
  
} # the end of function 'characteristics_print'

#characteristics_print(tmp) 


hist_events_model <- function(data_all, data_model, event, range=c(), id="pat_n",tit=""){
  #par(mfrow=c(1,2))
  
  if( ( nrow(data_all)==0 | all(is.na(data_all[,paste0(event,"_date")])) ) & (nrow(data_model)==0 | all(is.na(data_all[,paste0(event,"_date")])) ) ) return()

  if(length(range)==0) { 
    range <- range(c(data_all[,paste0(event,"_date")],data_model[,paste0(event,"_date")]),na.rm=T) 
    if("vax_date" %in% names(data_all)) range <- range( c(range,data_all[,"vax_date"]), na.rm=T)
    #range <- range(data_model[,paste0(event,"_date")])
    if(substring(as.character(range[1]),9,10)!="01") range[1] <- as.Date(paste0(substring(as.character(range[1]),1,8),"01"))
    if(!(substring(as.character(range[2]),9,10) %in% "01")) {  
      if(substring(as.character(range[2]),6,7)!="12") range[2] <- as.Date(paste0(substring(as.character(range[2]),1,5), format(as.numeric(substring(as.character(range[2]),6,7))+1,2),"-01"))
      else range[2] <- as.Date(paste0(format(as.numeric(substring(as.character(range[2]),1,4))+1,4), "-01-01"))
    }
  }
  # breaks:
  vertical_lines <- as.Date(paste0("20",rep(19:24,each=4),"-",rep(c("01","04","07","10"),6),"-01"))
  #vertical_lines <- as.numeric(difftime( as.Date(paste0("20",rep(20:24,each=4),"-",rep(c("01","04","07","10"),5),"-01")),as.Date("2020-08-31")),units="days")
  
  yy <- as.numeric(substring(as.character(range[1]),1,4)):as.numeric(substring(as.character(range[2]),1,4))
  yy_n <- break_months_all <- c()
  for(i in yy) {
    break_months <- 1
    if( as.numeric(substring(as.character(range[1]),1,4))==i ) break_months <- as.numeric(substring(as.character(range[1]),6,7))
    if( as.numeric(substring(as.character(range[2]),1,4))==i ) break_months <- break_months:as.numeric(substring(as.character(range[2]),6,7))
    else break_months <- break_months:12
    yy_n <- c(yy_n, length(break_months))
    break_months_all <- c(break_months_all, break_months)
  }
  
  breaks <- as.Date(paste0(rep(yy,yy_n),"-",format(break_months_all,2),"-01"))
  all_event_dates <- data_all[data_all[,paste0(event,"_date")]>=range[1] & data_all[,paste0(event,"_date")]<range[2],]
  all_vax_dates   <- all_event_dates[!duplicated(all_event_dates[,id],all_event_dates$vax_date),]
  all_atrisk      <- all_event_dates[,]
  all_event_dates <- all_event_dates[!duplicated(all_event_dates[,id],all_event_dates[,paste0(event,"_date")]), paste0(event,"_date")]
  
  
  # hist 3:
  if(length(all_event_dates)>0){
    
   # cal_time_range_days <- as.numeric(difftime( c(min(as.Date("2019-01-01"),as.Date(all_event_dates[!is.na(all_event_dates)])),max(as.Date(all_event_dates[!is.na(all_event_dates)]))+7),as.Date("2020-08-31"),units="days"))
   # breaks <- as.Date("2020-08-31") + seq( cal_time_range_days[1],cal_time_range_days[2], by=7) 
  
    hist_events_res <- hist(all_event_dates,breaks,plot=F)
    hist(all_event_dates,breaks,xlim=range,ylim=c(0,max(hist_events_res$counts)),axes=F, main=paste(event,";  (", tit,")"), xlab="Time", freq=T,col="skyblue");box();par(new=T)
    axis(2); axis(1,at=vertical_lines, labels =vertical_lines, cex=0.8); box()
    abline(v=vertical_lines,col="lightgray");par(new=T)
    hist(all_event_dates,breaks,xlim=range,ylim=c(0,max(hist_events_res$counts)),axes=F, main="", xlab="", ylab="", freq=T,col="skyblue");box();par(new=T)
    #print(table(all_vax_dates$vax_n))  
    #print(table(data_all[data_all[,paste0(event,"_date")]>=range[1] & data_all[,paste0(event,"_date")]<range[2],"vax_n"]))  
    legend("topleft", legend=c("in dataset","in SCRI (before and after)","in SCRI after vax"), fill=c("skyblue","blue","red"),bty="n")
  }
  
  model_event_dates <- data_model[data_model[,paste0(event,"_date")]>=range[1] & data_model[,paste0(event,"_date")]<range[2],]
  model_event_dates <- model_event_dates[model_event_dates[,id] %in% unique(model_event_dates[model_event_dates[event]==1,id]),]
  model_events_only_dates <-  model_event_dates[model_event_dates[,event]>0,]
  model_events_only_dates <-  model_events_only_dates[!duplicated(model_events_only_dates[,id],model_events_only_dates[,paste0(event,"_date")]),]
  # hist 3:
  if(nrow(model_events_only_dates)>0){ par(new=T)
    hist( model_events_only_dates[,paste0(event,"_date")], breaks, axes=F,main="",xlab="",xlim=range, ylim=c(0,max(hist_events_res$counts)), col="blue",freq=T);par(new=T)
    hist( model_events_only_dates[model_events_only_dates[,"vax_date_v1"]<=model_events_only_dates[,paste0(event,"_date")],paste0(event,"_date")], breaks, axes=F,main="",xlab="",xlim=range, ylim=c(0,max(hist_events_res$counts)), col="red",freq=T)
  }
  
  if(nrow(data_all)>0){
    obs_per_day_all_data <- obs_per_day_calc(data_all[!duplicated(data_all[,id]),], start="study_entry_date",  stop="study_exit_date", performance=F)
    obs_per_day_all_data <- tapply(obs_per_day_all_data,cut(as.Date(attributes(obs_per_day_all_data)$all_obs_days,origin="1970-01-01"),breaks),sum)
  }
  
  if(length(all_event_dates)>0){
    hist_res <- hist_events_res
    #print(hist_res$counts) 
    #print(obs_per_day_all_data) 
    #print(length(hist_res$counts) );  print(length(obs_per_day_all_data) )
    #print(as.Date(hist_res$breaks,origin="1970-01-01"))
    #print(hist_res)
    
    #plot(hist_res, xlim=range, ylim=c(0,max(hist_res$counts)), main=paste0("#",event,"; all data"))
    
    hist_res$counts <- obs_per_day_all_data
    plot(hist_res, xlim=range, freq=T,ylim=c(0,max(hist_res$counts)),axes=F,main=paste0("at risk; all data ( ",tit," )") )
    axis(2); axis(1,at=vertical_lines, labels =vertical_lines, cex=0.8); box()
    abline(v=vertical_lines,col="lightgray"); par(new=T)
    plot(hist_res, xlim=range, freq=T,ylim=c(0,max(hist_res$counts)),axes=F,main="",xlab="",ylab="")
    
    
    if(F){
      hist_res$counts <- 1000000*hist_events_res$counts/hist_res$counts
      print(hist_res$counts)
      plot(hist_res, xlim=range, freq=T, ylim=c(0,max(hist_res$counts)),main=paste0(1000000, " * #",event," / at risk;  all data"))
    }
  }  
  
  if(F){
    
    #  model data:
    if(nrow(model_events_only_dates)>0){
      start_obs_id <- tapply(model_event_dates$rw_start, model_event_dates[,id],min,na.rm=T)
      end_obs_id   <- tapply(model_event_dates$rw_end,   model_event_dates[,id],max,na.rm=T)
      obs_per_day_model_data <- obs_per_day_calc(cbind.data.frame(study_entry_date=as.Date("2020-08-31")+start_obs_id,
                                                                  study_exit_date =as.Date("2020-08-31")+end_obs_id),
                                                 start="study_entry_date",  stop="study_exit_date", performance=F)
      
      obs_per_day_model_data <- tapply(obs_per_day_model_data,cut(as.Date(attributes(obs_per_day_model_data)$all_obs_days,origin="1970-01-01"),breaks),sum)
      print(obs_per_day_model_data)
      obs_per_day_model_data[is.na(obs_per_day_model_data)] <- 0
      
      
      hist_res$counts <- obs_per_day_model_data
      plot(hist_res, xlim=range, ylim=c(0,max(hist_res$counts)),main="at risk; model data")
      
      hist_res$counts <- 1000000*hist(model_events_only_dates[,paste0(event,"_date")],breaks,xlim=range,main=paste0("#",event,"; model data"),freq=T)$counts/obs_per_day_model_data
      
      print(hist(model_events_only_dates[,paste0(event,"_date")],breaks,xlim=range,freq=T,plot=F)$counts)
      hist_res$counts[is.nan(hist_res$counts)] <- 0
      plot(hist_res, xlim=range, freq=T, ylim=c(0,max(hist_res$counts)),main=paste0("#",event," / at risk;  model data"))
    }
    
    
    
    if(length(all_vax_dates)>0){
      print(summary(all_vax_dates[,"vax_date"]))
      print(breaks)
      
      # hist 4,5:
      hist_res<-hist( all_vax_dates[,"vax_date"], breaks,freq=T,plot=F)
      hist( all_vax_dates[,"vax_date"], breaks,xlim=range,ylim=c(0,max(hist_res$counts)),freq=T,col="green3", main="all data",xlab="Date of vaccination")
      for(ivax in max(all_vax_dates$vax_n,na.rm=T):1)
        if(any(all_vax_dates$vax_n<=ivax)){par(new=T)
          hist( all_vax_dates[all_vax_dates$vax_n<=ivax,"vax_date"], breaks,freq=T,xlim=range,ylim=c(0,max(hist_res$counts)),col=rgb(0,ivax/max(all_vax_dates$vax_n,na.rm=T),0,0.5), main="all data",xlab="Date of vaccination")
        }
    }
    if(nrow(model_event_dates)>0) hist( model_event_dates[,"vax_date"], breaks,xlim=range, col=rgb(0,0,0.8,0.5))
    # hist 6:
    if(length(all_vax_dates)>0)
      hist_res <- hist( all_vax_dates[,"vax_date"], breaks,xlim=range,freq=T,col="green", main="all data",xlab="Date of vaccination"); box();par(new=T)
    if(nrow(model_event_dates)>0) 
      hist( model_event_dates[model_event_dates[,event]==1,"vax_date"], breaks,xlim=range, ylim=c(0,max(hist_res$counts)), axes=F,main="",xlab="",freq=T,col="red")
  }
  
  
} # end func "hist_events_model"




#
#   the end of functions :
#  
###############################
