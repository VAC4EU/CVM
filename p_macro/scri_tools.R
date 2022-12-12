
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
                      time_seq, split_seq_name = "cal_time_cat", time_seq_ref="most events", 
                      time_dep = NA,                # list of risk/control windows definitions
                      combine_vars = c(),           # list of parameters to create new one variable from two other variables
                      start_obs, end_obs,
                      data, strata_var="", strata_value=NA, use_all_events=F,
                      data_event,
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
                      save_data = F
){  
  
  #if(missing(formula   )) stop("'formula' is missing.")
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
  
  if(strata_var!=""){ 
    if(!is.na(strata_value)) data$strata_cond <- data[,strata_var]==strata_value & !is.na(data[,strata_var])
    else                     data$strata_cond <- data[,strata_var]               & !is.na(data[,strata_var])
    if( sum(data$strata_cond)==0) return(  list( tab      = NULL, 
                                                 tab_full = NULL,
                                                 model    = NULL,
                                                 call     = list( match.call())  ))
    if(!use_all_events) data <- data[data$strata_cond,]
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
    obj=vax_def,
    rws       = rws,
    data =data,
    strata_cond = ifelse(strata_var=="", F, T), strata_value = strata_value,
    start_obs = start_obs, end_obs = end_obs,
    event_time = event_time, event = event,
    id = id,
    lab_orders = lab_orders,
    ref="reference ",   #  ref=5 OR ref= "reference [-90;-30]"  or a part of a ref.category name
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
  
  
  
  #######
  #  create time intervals:
  #
  if(!missing(time_seq) & nrow(data_rws)>0){   
    data_rws <- refresh_event_variable( "rw_start", "rw_end", data_rws, event, event_time)
    data_rws <- split_intervals( data =data_rws, 
                                 start_interval = "rw_start", end_interval = "rw_end", 
                                 splits_names = split_seq_name, # "cal_time_cat",
                                 splits       = time_seq,
                                 ref          = time_seq_ref, #"most events", 
                                 event        = event,        #   myopericarditis 
                                 event_time   = event_time 
    )
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
  
  
  if(nrow(data_rws)==0) return(  list( res_tab = NULL, 
                                       model   = NULL,
                                       call    = list( match.call())  ))
  
  
  if(any(names(data_rws)==split_seq_name))
    data_rws[,split_seq_name]  <- factor_ref(  as.character(data_rws[,split_seq_name]), 
                                               lab=levels(data_rws[,split_seq_name])[levels(data_rws[,split_seq_name]) %in% unique(as.character(data_rws[,split_seq_name])) ],
                                               ref=time_seq_ref, 
                                               event_var=data_rws[,event] )  
  
  
  # combine variables in formula with risk windows variable 'lab' (the begining of the code):
  if(!no_formula){
    tb <- attributes(terms(formula))$factor
    if(!missing(time_seq) & split_seq_name %in% dimnames(tb)[[2]])
      tb <- tb[, c( dimnames(tb)[[2]][dimnames(tb)[[2]]!=split_seq_name], split_seq_name ) ]
  }
  
  
  if(length(sep_vars)>0){
    
    mods <- vector("list",length=length(sep_vars))
    names(mods) <- names(sep_vars)
    
    for(isep in 1:length(sep_vars)){
      
      var_names <- sep_vars[[isep]]
      if(!missing(time_seq)) #& split_seq_name %in% dimnames(tb)[[2]])
        var_names <- c( var_names, split_seq_name )
      var_names <- var_names[ lapply(var_names, function(x)nlevels(data_rws[,x]))>1 ]
      
      
      #  data_rws$cal_time_mid <- I(unlist(lapply(strsplit(substring(as.character(data_rws$cal_time_cat),2,nchar(as.character(data_rws$cal_time_cat))-1),";",fixed=T),function(x)mean(as.numeric(x)))))
      Poisson_formula <- as.formula(paste(event,"~", 
                                          paste( var_names,  collapse=" + "), "+", 
                                          "strata(",id,")", "+", "offset(log(interval))")) 
      #cat(paste("before ",Poisson_formula,":",Sys.time(),"\n"))  
      
      if(!no_formula){
        suppressWarnings(
          mod <- try( clogit(formula = Poisson_formula, data = data_rws, control=coxph.control(iter.max=1000) ), silent = T )
        )
      }
      else mod <- NA
      #cat(paste("after ",Poisson_formula,":",Sys.time(),"\n"))  
      
      
      if(class(mod)[[1]]== "try-error" | no_formula){
        
        if(F){
          if(nrow(table1(data_rws[  data_rws[,event]==1 ,var_names]))<=8)
            print( table1(data_rws[  data_rws[,event]==1 ,var_names]) ) 
          print(cbind.data.frame(var_names, nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names])) ))
          if(any(names(data_rws)==split_seq_name))   #"cal_time_cat"))
            print(cbind.data.frame(var_names[var_names!=split_seq_name], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!=split_seq_name]])) ))
          
          if(!missing(time_seq) )
            warning("Error in Poisson regression. \ntime_seq=",median(diff(time_seq),na.rm=T),"; Formula: ", deparse(Poisson_formula))
          else warning("Error in Poisson regression. \nno time_seq; Formula: ", deparse(Poisson_formula))
        }        
        
        if(isep==1)
          res_tab <- summary_tab( var_names=var_names,#sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  lab_orders=lab_orders,
                                  model_number = isep )
        else  
          res_tab <- summary_tab( var_names=var_names,#sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  lab_orders=lab_orders,
                                  model_number = isep,
                                  res_tab = res_tab )
        
        mods[[isep]] <- list( mod ) 
        
      }
      else{
        
        
        if(isep==1)
          res_tab <- summary_tab( var_names=var_names,#sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  #coef_cond = coef_cond, 
                                  model_number = isep )
        else  
          res_tab <- summary_tab( var_names=var_names,#sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  #coef_cond=coef_cond, 
                                  model_number = isep,
                                  res_tab = res_tab )
        
        mods[[isep]] <- list(summary(mod)) 
      }
    }
    
    mod <- mods
    
    if(lprint)  print(format( res_tab, digits=3, nsmall=2, justify="left" ))
    
  } else {
    
    if(nrow(data_rws)==0) res_tab <- mod <- NULL
    else {
      
      if(any(cond<-(colSums(tb)>1)) & nrow(data_rws)>0)
        for(i in 1:sum(cond)){
          vars_tmp <- gsub(" ","",dimnames(tb)[[2]][cond][i])
          dimnames(tb)[[2]][cond][i] <- paste(strsplit( vars_tmp, ":" )[[1]],collapse="_")
          data_rws[,paste(strsplit( gsub(" ","",vars_tmp,fixed = TRUE) ,":")[[1]],collapse="_")]  <- combine_vars_func( 
            data_rws[, strsplit(gsub(" ","",vars_tmp,fixed = TRUE)  ,":")[[1]], drop=F ], lab_orders = lab_orders, event = data_rws[,event] )   #ref=ref, 
        }  
      
      var_names <- dimnames(tb)[[2]]
      if(!missing(time_seq) & split_seq_name %in% dimnames(tb)[[2]])
        var_names <- unique(c( var_names, split_seq_name ))
      var_names <- var_names[ lapply(var_names, function(x)nlevels(data_rws[,x]))>1 ]
      
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
      #cat(paste("before ",Poisson_formula,":",Sys.time(),"\n"))  
      if(!no_formula){
        suppressWarnings(
          mod <- try( clogit(formula = Poisson_formula, data = data_rws, control=coxph.control(iter.max=1000) ), silent = T )
        )
      }
      #cat(paste("after ",Poisson_formula,":",Sys.time(),"\n"))  
      
      if(class(mod)[[1]]== "try-error" | no_formula) {
        
        if(F){
          if(nrow(table1(data_rws[  data_rws[,event]==1 ,var_names]))<=8)
            print( table1(data_rws[  data_rws[,event]==1 ,var_names]) ) 
          print(cbind.data.frame(var_names, nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names])) ))
          if(any(names(data_rws)==split_seq_name))   #"cal_time_cat"))
            print(cbind.data.frame(var_names[var_names!=split_seq_name], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!=split_seq_name]])) ))
          #print(cbind.data.frame(var_names[var_names!="cal_time_cat"], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!="cal_time_cat"]])) ))
          
          if(!missing(time_seq))
            warning("Error in Poisson regression. \ntime_seq=",median(diff(time_seq),na.rm=T),"; Formula: ", deparse(Poisson_formula))
          else warning("Error in Poisson regression. \nno time_seq; Formula: ", deparse(Poisson_formula))
        }   
        
        res_tab <- summary_tab( var_names=dimnames(tb)[[2]], event=event, data=data_rws, id_name=id, lab_orders=lab_orders )
      }
      else{
        res_tab <- summary_tab( var_names=dimnames(tb)[[2]], event=event, data=data_rws, id_name=id,  mod=mod, lab_orders=lab_orders  )
        mod <- summary(mod)
      }
      
      if(lprint)  print(format( res_tab, digits=2, justify="left" ))
    }
  }
  
  
  res_tab0 <- res_tab
  if(delete_coef_no_events){ 
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
    res_tab <- cbind.data.frame(res_tab[,1:match("i",names(res_tab))], 
                                strata=ifelse(!is.na(strata_value),strata_value,strata_var), 
                                res_tab[,match("i",names(res_tab)):ncol(res_tab)]               )
  
  
  ret <- list( tab      = res_tab,
               tab_full = res_tab0,
               model    = mod,
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
#   !!!!!!!!!!!!!!!!!!!!
#
#
#             create rws (version 3)  31.01.2022
#
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
    
    if(nrow(data)==0) return("no data")
    
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
summary_tab <- function(  var_names, # var_names <- c("lab", "cal_time_cat")
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
  if(missing(var_names)) stop("'var_names' is missing")
  if(missing(data) & missing(mod)) stop("'data' and/or 'mod' must be specified. ")
  
  ######
  # create summary table for variable: 'var_names' from 'data' (variable 'interval', 'event',id_name must be in 'data')
  if(!missing(data)){
    
    cond <- !is.na(data$interval) & data$interval>0
    res_tab <- vector("list",length(var_names))
    for(ivar in var_names){
      
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
      
      res_tab[ivar==var_names][[1]] <- cbind.data.frame(
        
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
  # create summary table for result from Poisson regression with variables 'var_names':
  #
  if(!missing(mod)){  
    res_tab_model <- cbind.data.frame( summary(mod)$conf.int[,-2], 
                                       summary(mod)$coef[,match( c("Pr","co","se"), substring(dimnames(summary(mod)$coef)[[2]],1,2) ) ]  ) #  "Pr(>|z|)" "coef"   "se(coef)" 
    res_tab_model[ is.na(res_tab_model$coef), substring(dimnames(res_tab_model)[[2]],1,2)=="se" ] <- NA
    # delete variable names at the beginning of parameter names in results from Poisson regression:
    res_tab_model$all_cat <- dimnames(res_tab_model)[[1]]
    var_names <- unique(unlist(strsplit(var_names,":")))
    var_names <- var_names[order(nchar(var_names),decreasing = T)]
    for(ivar in var_names){
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
              ( !is.na(res_tab_new[, ivar])  | !is.na(res_tab_new[,paste0(ivar,".y")]) ) )) next
      
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
  
  if(all(unlist(lapply(res, function(x)sum(!is.na(x$RR))))==0)) return(NULL) 
  
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
                
                for(ii in ref_models ){
                  cond_ref_model <- !is.na(model_time_adj) & model_time_adj==ii
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
                 formula="", time_seq=c(), time_seq_ref="most events",
                 event_time, event, event_date, event_info,
                 id,
                 rws ,
                 ref,
                 time_dep = NA,
                 combine_vars = c(),
                 start_obs, end_obs,
                 lab_orders = NA,
                 data,  strata_var="", strata_value=NA, use_all_events=F,
                 data_event  = NA,
                 rw_observed_percentage=0,   # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                 censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                 #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                 event_in_rw=T,                # if event in rw ==> this rw should not be deleted even if not completely observed
                 lplots     = T,
                 leventplot = T, max_n_points=NA,  eventplot_file_separate=F, 
                 lplot = T, CI = T,
                 lforest=T, forest_nrows=50 ,forest_cex=0.5, forest_cex_head=0.5,
                 nvax,
                 delete_coef_no_events      = T,
                 delete_rows_start_with_no  = T,
                 delete_no_ref_cols         = T,
                 delete_median_cols         = T,
                 paral = T,  n_cores=NA,          # paral_vars = c(), parallel programming (using more cores)
                 lprint=T, 
                 save_data = F,
                 width=14, 
                 extra_parameters,
                 add_to="",
                 add_to_itself = F,
                 ndigits = 2,
                 sdr_tabs, sdr_models, 
                 cut_points_name = "",
                 col = c("red", "green3", "orange",  "deepskyblue", "magenta2", "cyan2", "chocolate1", "gray" ),
                 performance=F,
                 ...
){  

  if(missing(data)) stop("Dataset 'data' is missing.")
  
  start_sys_time <- Sys.time()

  # stratum:
  # if only for subset (i.e., use_all_events==F) ==> delete all other rows
  if(strata_var!="" & !use_all_events) {
    data <- data[ data[,strata_var]== strata_value & !is.na(data[,strata_var]), ]
    strata_var <- ""
  }
 
  if(nrow(data)==0) return("no data") 
  
  
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
  
  if(!missing(extra_parameters)){
    
    if( missing(extra_name ) & "extra_name"  %in% names(extra_parameters) ) extra_name <- extra_parameters[["extra_name"              ]] 
    
    if( missing(lplots                 ) & "lplots"                  %in% names(extra_parameters) ) lplots                  <- extra_parameters[["lplots"            ]]   
    
    if( missing(leventplot             ) & "leventplot"              %in% names(extra_parameters) ) leventplot              <- extra_parameters[["leventplot"              ]]   
    if( missing(max_n_points           ) & "max_n_points"            %in% names(extra_parameters) ) max_n_points            <- extra_parameters[["max_n_points"            ]]   
    if( missing(eventplot_file_separate) & "eventplot_file_separate" %in% names(extra_parameters) ) eventplot_file_separate <- extra_parameters[["eventplot_file_separate" ]] 
    
    if( missing(lplot           ) & "lplot"             %in% names(extra_parameters) ) lplot            <- extra_parameters[["lplot"            ]]   
    if( missing(CI              ) & "CI"                %in% names(extra_parameters) ) CI               <- extra_parameters[["CI"               ]]
    if( missing(CI              ) & "CI_draw"           %in% names(extra_parameters) ) CI               <- extra_parameters[["CI_draw"          ]]
    
    if( missing(lforest         ) & "lforest"         %in% names(extra_parameters) ) lforest          <- extra_parameters[["lforest"        ]]   
    if( missing(forest_nrows    ) & "forest_nrows"    %in% names(extra_parameters) ) forest_nrows     <- extra_parameters[["forest_nrows"   ]]   
    if( missing(forest_cex      ) & "forest_cex"      %in% names(extra_parameters) ) forest_cex       <- extra_parameters[["forest_cex"     ]]   
    if( missing(forest_cex_head ) & "forest_cex_head" %in% names(extra_parameters) ) forest_cex_head  <- extra_parameters[["forest_cex_head"]]   
    
    if( missing(delete_coef_no_events    ) & "delete_coef_no_events"     %in% names(extra_parameters) ) delete_coef_no_events     <- extra_parameters[[ "delete_coef_no_events"     ]]   
    if( missing(delete_coef_no_events    ) & "delete_coef_no_events"     %in% names(extra_parameters) ) delete_coef_no_events     <- extra_parameters[[ "delete_coef_no_events"     ]]   
    if( missing(delete_rows_start_with_no) & "delete_rows_start_with_no" %in% names(extra_parameters) ) delete_rows_start_with_no <- extra_parameters[[ "delete_rows_start_with_no" ]]   
    if( missing(delete_no_ref_cols       ) & "delete_no_ref_cols"        %in% names(extra_parameters) ) delete_no_ref_cols        <- extra_parameters[[ "delete_no_ref_cols"        ]]   
    if( missing(delete_median_cols       ) & "delete_median_cols"        %in% names(extra_parameters) ) delete_median_cols        <- extra_parameters[[ "delete_median_cols"        ]]   
    
    
    if( missing(paral           ) & "paral"             %in% names(extra_parameters) ) paral            <- extra_parameters[["paral"            ]]   
    if( missing(n_cores         ) & "n_cores"           %in% names(extra_parameters) ) n_cores          <- extra_parameters[["n_cores"          ]] 
    
    if( missing(lprint          ) & "lprint"            %in% names(extra_parameters) ) lprint           <- extra_parameters[["lprint"           ]]   
    if( missing(save_data       ) & "save_data"         %in% names(extra_parameters) ) save_data        <- extra_parameters[["save_data"        ]]   
    if( missing(width           ) & "width"             %in% names(extra_parameters) ) width            <- extra_parameters[["width"            ]]
    
    if( missing(time_seq        ) & "time_seq"          %in% names(extra_parameters) ) time_seq         <- extra_parameters[["time_seq"         ]]   
    if( missing(sdr_tabs        ) & "sdr_tabs"          %in% names(extra_parameters) ) sdr_tabs         <- extra_parameters[["sdr_tabs"         ]]   
    if( missing(sdr_models      ) & "sdr_models"        %in% names(extra_parameters) ) sdr_models       <- extra_parameters[["sdr_models"       ]] 
    
    if( missing(performance     ) & "performance"       %in% names(extra_parameters) ) performance      <- extra_parameters[["performance"      ]]   
  }
  
  if(!missing(event_info)){
    
    if(missing(event_time)) event_time <- event_info[["event_time"]]
    if(missing(event_date)) event_date <- event_info[["event_date"]]
    if(missing(event     )) event      <- event_info[["event"     ]]
    
    if(any(!is.na(data[,event_date])) & all( is.na(data[,event_time]))) data[,event_time] <- as.numeric(difftime(data[,event_date],as.Date("2020-08-31"),units="days")) 
    if(all( is.na(data[,event_date])) & any(!is.na(data[,event_time]))) data[,event_date] <- as.Date("2020-08-31") + data[,event_time]
  }
  data[,event] <- as.numeric(!is.na(data[,event_date]))
  
  if(!any(data[,event]==1)) return("no data")
  
  #######################################################
  #
  #             output structure, content of output files, directories:
  #
  ################
  # directories:
  #
  #     directory 1: g_export / scri / { ?dap }_{ event }  /   ?  /{output_name}.txt; {output_name}.RData; {output_name}_{stratum ?or all strata? }.pdf
  #
  #     directory 2: g_?local? / scri / { event }  / ?    /{output_name}_models.RData
  #
  ###############
  # 'output_name'
  #
  #     output_name <- { dap }_{ event }_[ extra_name ]_[ nothing or covid_selection ]_{ no_split; brand  via vax_dep 'before' }_{ dist  via vax_dep 'after'}_{ nothing or stratum }
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
  
  output_name <- ""
  if(any(ls()=="dap"))                            output_name <- c( output_name, dap                                                                  )
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

#browser()  
#cat(paste("before image plot:",Sys.time(),"\n"))  

  
  
  # create leventplot (image) plots:
  if(leventplot & lplots){  
    
    if(any( !( names(list(...)) %in% names(formals(pdf)) ) )) 
      stop(paste0("Argument[s] '", paste(names(list(...))[ !( names(list(...)) %in% names(formals(pdf)) )], collapse="', '" ),
                  "' is not an argumnent of function 'pdf' from library 'qpdf'."))
    
    gc()
    pdf(file=paste0(sdr_tabs,"eventplot_tmp.pdf"), width=width,  ... )

    try(image_plots( vax_def=vax_def, data=data, event_info=event_info,  strata_var=strata_var, strata_value=strata_value, use_all_events=use_all_events, tit=output_name, max_n_points=max_n_points ),silent=T)

    dev.off()
    #cat(paste("after image plot:",Sys.time(),"\n"))  

  }
  
  
  data <- data[data[,event]==1 & !is.na(data[,vax_def$data_parameters$vax_time]),]
  
  
  if( length(time_seq)>1 & length(time_seq_ref)==1 ) time_seq_ref <- rep(time_seq_ref,length(time_seq))
  
  if(any(ls()=="res")) rm(res)
  
  res <- vector("list", length=length(time_seq)+1)
  names(res) <- c("no_adj", names(time_seq) )
  
  if(strata_var!="" & use_all_events) names(res) <- paste0(names(res), "_all_events")

  #cat(paste("before scri_fit:",Sys.time(),"\n"))  
  
  res[[1]] <-  scri_fit(formula = formula,
                        vax_def=vax_def,
                        event_time = event_time, event = event, id=id,
                        rws        = rws, 
                        ref=ref,
                        time_dep   = time_dep,
                        combine_vars = combine_vars,
                        start_obs  = start_obs, end_obs = end_obs,
                        lab_orders = lab_orders,
                        data       = data, strata_var=strata_var, strata_value=strata_value, use_all_events=use_all_events,
                        data_event = data_event,
                        nvax       = nvax,
                        delete_coef_no_events      = delete_coef_no_events,     
                        delete_rows_start_with_no  = delete_rows_start_with_no, 
                        delete_no_ref_cols         = delete_no_ref_cols,        
                        delete_median_cols         = delete_median_cols,        
                        rw_observed_percentage     = rw_observed_percentage,          
                        censored_vars = censored_vars,           
                        event_in_rw = event_in_rw,               
                        lprint=lprint,
                        save_data = save_data)
  #cat(paste("after scri_fit:",Sys.time(),"\n"))  
  
  #############################
  #
  # parallel use of cores:
  #
  if(length(time_seq)>0 & paral){
    
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
                             scri_fit(formula = ifelse(formula=="","", paste( formula, " +  cal_time_cat") ),
                                      vax_def=vax_def,
                                      event_time = event_time, event = event, id=id,
                                      rws = rws, 
                                      ref=ref,
                                      time_seq = i_time_seq,  split_seq_name = "cal_time_cat", #time_seq_ref="most events",time_seq_ref=time_seq_ref[i_time_seq],
                                      time_dep   = time_dep, 
                                      combine_vars = combine_vars,
                                      start_obs = start_obs, end_obs = end_obs,
                                      lab_orders = lab_orders,
                                      data = data, strata_var=strata_var, strata_value=strata_value, use_all_events=use_all_events,
                                      data_event  = data_event,
                                      nvax       = nvax,
                                      delete_coef_no_events      = delete_coef_no_events,    
                                      delete_rows_start_with_no  = delete_rows_start_with_no,
                                      delete_no_ref_cols         = delete_no_ref_cols,       
                                      delete_median_cols         = delete_median_cols,       
                                      rw_observed_percentage     = rw_observed_percentage,   
                                      censored_vars = censored_vars,           
                                      event_in_rw = event_in_rw,               
                                      lprint=lprint,
                                      save_data = save_data)
    )
    
    stopCluster(cl)
    
    names_overlap <- match(names(res),names(res_paral))
    res[ names(res_paral)[ names_overlap[!is.na(names_overlap)] ] ] <- res_paral[!is.na(match(names(res_paral),names(res)))]
    
  }  # end of parallel
  #cat(paste("after parallel scri_fit:",Sys.time(),"\n"))  
  
  
  if(length(time_seq)>0 & !paral)
    for(i in 2:length(res)){
      res[[i]] <-  scri_fit(formula = ifelse( formula=="","", paste( formula, " +  cal_time_cat") ),
                            vax_def=vax_def,
                            event_time = event_time, event = event, id=id,
                            rws = rws, 
                            ref=ref,
                            time_seq = time_seq[[i-1]] , time_seq_ref=time_seq_ref[i-1], split_seq_name = "cal_time_cat", #time_seq_ref="most events", 
                            time_dep   = time_dep,
                            combine_vars = combine_vars,
                            start_obs = start_obs, end_obs = end_obs,
                            lab_orders = lab_orders,
                            data = data, strata_var=strata_var, strata_value=strata_value, use_all_events=use_all_events,
                            data_event   = data_event,
                            nvax       = nvax,
                            delete_coef_no_events      = delete_coef_no_events,    
                            delete_rows_start_with_no  = delete_rows_start_with_no,
                            delete_no_ref_cols         = delete_no_ref_cols,       
                            delete_median_cols         = delete_median_cols,       
                            rw_observed_percentage     = rw_observed_percentage,     
                            censored_vars = censored_vars,           
                            event_in_rw = event_in_rw,               
                            lprint=lprint,
                            save_data = save_data)#[[1]]
      
    }  
  #cat(paste("after non-parallel scri_fit:",Sys.time(),"\n"))  
  
  
  # for output file with 'tabs'  
  tabs             <- lapply(res,function(x) x$tab )
  attributes(tabs) <- c( attributes(tabs), lapply(res,function(x)x[names(x)==c("tab_full")]) )
  class(tabs)      <- "scri_tabs"
  attributes(tabs) <- c( attributes(tabs), 
                         name      = output_name,
                         vax_def = list(vax_def),                  
                         event   = list( c(event=event, event_time=event_time, event_date=event_date) ) )
  
  res <- c (res, name      = output_name,
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
    
    if(length(tabs)>0)
      for(i1 in 1:length(tabs)){
        txt_i1 <- paste0("\n\n",names(tabs)[i1])
        if(length(dim(tabs[[i1]]))==2){
          cat(paste(txt_i1,":\n\n"))
          print(format(tabs[[i1]],              justify="left", digits=ndigits))
        }
        
        if(length(tabs[[i1]])>0)
          for(i2 in 1:length(tabs[[i1]])){
            txt_i2 <- paste0(txt_i1," $ ",names(tabs[[i1]])[i2])
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
  
  
  ret <- list( tabs   = tabs,
               models = res
  )
  class(ret) <- "scri_output"
  
  
  if((lplot | lforest) & lplots){  
    
    pdf(file=paste0(sdr_tabs,"tmp.pdf"), width=width,  ...)
    if(!is.null(tabs[[1]])){
      # create forest plot:
      if(lforest) {
        if(missing(forest_nrows)) 
          forest_nrows <-  pmax(1, 58 - 4*length(tabs[[i1]]) + pmin(1,length(tabs[[i1]])%/%11)*(-2*10 + 2*length(tabs[[i1]]) ) + pmin(1,length(tabs[[i1]])%/%16)*( -15 + length(tabs[[i1]]) ) ) 
        par(mfrow=c(1,2))
        for(i1 in 1:length(tabs))
          forest_plots_tab( tabs[[i1]], nrows_forest_plot=forest_nrows,cex=forest_cex, cex_head=forest_cex_head, ltable=F, col=col )  
      }
      # create plots with coefficients:
      if(lplot){
        par(mfrow=c(1,1))
        for(i1 in 1:length(tabs))
          plot_res(tabs[[i1]], main=paste( event, formula,"; \n",output_name,"; ", names(tabs)[i1]), col=col, CI=CI)
      }
    }
    dev.off()  # end pdf
    
    files_to_copy <- paste0(sdr_tabs,"tmp.pdf")
    if(leventplot & !eventplot_file_separate) files_to_copy <- c( files_to_copy, paste0(sdr_tabs,"eventplot_tmp.pdf") )
    
    for(iname in files_to_copy) if(!file.exists(iname)) files_to_copy <- files_to_copy[files_to_copy!=files_to_copy]
    if(length(files_to_copy)>0) qpdf::pdf_combine( files_to_copy , paste0(sdr_tabs,output_file_name,".pdf")  )  
    
    if(file.exists(paste0(sdr_tabs,"eventplot_tmp.pdf"  )))  suppressWarnings( file.remove(paste0(sdr_tabs, "eventplot_tmp.pdf"   )) )
    if(file.exists(paste0(sdr_tabs,          "tmp.pdf"  )))  suppressWarnings( file.remove(paste0(sdr_tabs,           "tmp.pdf"   )) )
    if(!is.null(dev.list())) dev.off()
  }  # end plots
  
  #cat(paste("end :",Sys.time(),"\n"))  
  if(performance) cat(paste0(": duration = ",round(Sys.time()-start_sys_time),"s (till ",Sys.time(),")\n"))  
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


image_plots <- function(vax_def, event_info, data, tit="", strata_var="", strata_value=NA, use_all_events=T, max_n_points=NA){  

  id        <- vax_def$data_parameters$id
  vax_name  <- vax_def$data_parameters$vax_name
  vax_date  <- vax_def$data_parameters$vax_date
  vax_time  <- vax_def$data_parameters$vax_time
  
  event      <- event_info$event           
  event_time <- event_info$event_time   
  event_date <- event_info$event_date  
  
  if(strata_var!="") {
    if(!is.na(strata_value)) data <- data[data[,strata_var]==strata_value & !is.na(data[,strata_var]), ]
    else                     data <- data[data[!is.na(data[,strata_var]),strata_var], ]
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
  
  # create dataset with death rows (vaxed and unvaxed):
  cond_death <- F
  death_date_possible_names <- c("date_death","death_date","date_of_death","dateofdate","datedeath","deathdate","dod")
  if(!is.null(vax_def$data_parameters$censored_vars))
    death_date_possible_names <- unique(c( death_date_possible_names, 
                                           grep("death",vax_def$data_parameters$censored_vars, value=T ),
                                           grep("dead" ,vax_def$data_parameters$censored_vars, value=T ) ))
  if(any( (tmp_cond<-names(data) %in% death_date_possible_names) ))
    for(iname in names(data)[tmp_cond]) cond_death <- cond_death | !is.na(data[, iname])
  data_with_deaths <- as.data.frame(data[cond_death, ])
  for(iname in names(data)[tmp_cond]){
    if(any(class(data_with_deaths[,iname]) %in% c("POSIXct","POSIXt","Date"))) data_with_deaths$death_date <- data_with_deaths[,iname]  
    if(any(class(data_with_deaths[,iname]) %in% c("numeric","integer")))       data_with_deaths$death_days <- data_with_deaths[,iname]  
  }
  total_deaths <- length(unique(data_with_deaths[,id]))
  
  # the number of vaccinated and unvaccinated id's:
  n_unvaccinated <- length(unique(data[ is.na(data[,vax_date]),id]))
  n_vaccinated   <- length(unique(data[!is.na(data[,vax_date]),id]))
  
  
  # dataset unvaccinated and event 
  data_unvax_events <- as.data.frame(data[ data[,event]==1 & is.na(data[,vax_date]), ])
  
  # dataset with events and vaccinated
  data <- as.data.frame(data[ data[,event]==1 & !is.na(data[,vax_date]) & !is.na(data[,vax_name]) , ])
  all_vax_names <- levels(factor_ref(data[!duplicated(data[,vax_name]),vax_name], lab_orders = vax_def$lab_orders))
  
  vax_names_col <- c("red","magenta2","skyblue","violet","green3","darkorchid")
  if(length(vax_names_col)<length(all_vax_names)) vax_names_col <- rep(vax_names_col,length(all_vax_names)/5)
  vax_names_col <- vax_names_col[1:length(all_vax_names)]
  
  
  par(mfcol=c(3,3))
  
  at_date     <- seq( 2015,  as.numeric(substring(Sys.Date(),1,4))+1, by=0.25 ) 
  at_date_lab <- month.abb[c(1,4,7,10)][1+4*at_date%%1] 
  at_date_lab[at_date_lab==month.abb[1]] <- paste0( at_date_lab, "\n", at_date%/%1 )[at_date_lab==month.abb[1]]
  
  
  # create set of plots for each value of 'vaxdep_values' 
  for(idep in 1:length(vaxdep_values)){
    
    if( !(length(vaxdep_values)==1 & vaxdep_values[1]=="all") ){
      if(sum(data$vaxdep_all==vaxdep_values[idep])==0) next
      data$idep_cond <- data$vaxdep_all==vaxdep_values[idep]
      data$idep_value[data$idep_cond] <- vaxdep_values[idep]
      tit_dep <- vaxdep_values[idep]
    } else{  # no vaccine-dependent variables (==> no split analysis)
      data$idep_cond  <- rep(T,nrow(data))
      data$idep_value <- ""
      tit_dep <- ""
    }
    
    # per 'vax_number' or 'vax_name'  
    for(ivax in 1:length(all_vax_names)){
      
      if(sum(data[data$idep_cond,vax_name]==all_vax_names[ivax])==0) next
      data$ivax_cond  <- data[,vax_name]==all_vax_names[ivax]
      data$ivax_value <- all_vax_names[ivax]
      
      # create temporary variables with the current dose info:
      data_ivax <- data[data$idep_cond & data$ivax_cond,c(id,vax_time,vax_date,vax_name, event_time)]
      data_ivax$time_x  <- data_ivax[,event_time] - data_ivax[,vax_time]
      data_ivax$time_y1 <- data_ivax[,event_time]
      data_ivax$time_y2 <- data_ivax[,  vax_time]
      data_ivax <- data_ivax[,names(data_ivax)!=event_time]
      names(data_ivax)[2:4] <- paste0(c("vax_time","vax_date","vax_name"),"_ivax")
      data_ivax$idep_ivax_id <- rep(T,nrow(data_ivax))
      data <- merge(data[, !(names(data) %in% names(data_ivax)[-1]) ], data_ivax, by=id, all.x=T )
      data$ivax_cond  <- data[,"vax_name_ivax"]==all_vax_names[ivax]
      
      data_plot <- data[data$idep_ivax_id & !is.na(data$idep_ivax_id),]
      
      limits_x  <- range(data_plot$time_x,  na.rm=T)
      limits_y1 <- range(data_plot$time_y1, na.rm=T)
      limits_y2 <- range(data_plot$time_y2, na.rm=T)
      
      # create temporary variables with the current dose info in 'data_with_deaths' dataset:
      if(nrow(data_with_deaths)>0){
        data_ivax <- data_with_deaths[data_with_deaths[,vax_name]==all_vax_names[ivax],c(id,vax_time,vax_date,vax_name, "death_days")]
        data_ivax$time_x  <- data_ivax$death_days - data_ivax[,vax_time]
        data_ivax$time_y1 <- data_ivax$death_days
        data_ivax$time_y2 <- data_ivax[,vax_time]
        data_ivax <- data_ivax[,names(data_ivax)!=event_time]
        names(data_ivax)[2:4] <- paste0(c("vax_time","vax_date","vax_name"),"_ivax")
        data_ivax$ivax_id <- rep(T,nrow(data_ivax))
        data_ivax$death_days <- NULL
        data_with_deaths <- merge(data_with_deaths[, !(names(data_with_deaths) %in% names(data_ivax)[-1]) ], data_ivax, by=id, all.x=T )
        data_with_deaths$ivax_cond  <- data_with_deaths[,"vax_name_ivax"]==all_vax_names[ivax] & !is.na(data_with_deaths[,"vax_name_ivax"])
        
        limits_x  <- range(data_plot$time_x,  data_with_deaths$time_x,  na.rm=T)
        limits_y1 <- range(data_plot$time_y1, data_with_deaths$time_y1, na.rm=T)
        limits_y2 <- range(data_plot$time_y2, data_with_deaths$time_y2, na.rm=T)
  
      }
      
      limits_x  <- limits_x  + c(-1,1)*0.02*diff(limits_x)
      limits_y1 <- limits_y1 + c(-1,1)*0.06*diff(limits_y1)
      limits_y2 <- limits_y2 + c(-1,1)*0.06*diff(limits_y2)
      
      tit_ivax <- tit_ivax_events <- tit_ivax_deaths <- paste0(tit_dep, ifelse(tit_dep=="",""," & "), trimws(all_vax_names[ivax]))
      
      ###########
      # sample points if they are more than 'max_n_points'
      #
      data_plot$cond_sampled_ids <- rep(T,nrow(data_plot))
      if(nrow(data_with_deaths)>0)
        data_with_deaths$cond_sampled_ids <- rep(T,nrow(data_with_deaths))
      
      if(!is.na(max_n_points)){ 
        # for events:
        ivax_ids <- unique(data_plot[data_plot$ivax_cond & data_plot$time_x>=0,id])
        if( max_n_points < length(ivax_ids) ) {
          ivax_all_ids <- data_plot[data_plot$ivax_cond,id]
          sampled_ids <- sample(ivax_all_ids, max_n_points/length(ivax_ids) * length(ivax_all_ids))
          data_plot$cond_sampled_ids <- data_plot[,id] %in% sampled_ids & data_plot$ivax_cond
          tit_ivax_events <- paste0(tit_ivax, ifelse(tit_ivax=="",""," & "),"\n[", round(100*max_n_points/length(ivax_ids),1),"% of all ",length(ivax_all_ids)," events after ",all_vax_names[ivax],"]")
        }
        # for deaths:
        if(nrow(data_with_deaths)>0){
          ivax_ids <- unique(data_with_deaths[data_with_deaths$ivax_cond & data_with_deaths$time_x>=0,id])
          if( max_n_points < length(ivax_ids) ) {
            ivax_all_ids <- unique(data_with_deaths[data_with_deaths$ivax_cond,id])
            sampled_ids <- sample(ivax_all_ids, max_n_points/length(ivax_ids) * length(ivax_all_ids))
            data_with_deaths$cond_sampled_ids <- data_with_deaths[,id] %in% sampled_ids & data_with_deaths$ivax_cond
            tit_ivax_deaths <- paste0(tit_ivax, ifelse(tit_ivax=="",""," & "),"\n[", round(100*max_n_points/length(ivax_ids),1),"% of all ",length(ivax_all_ids)," deaths after ",all_vax_names[ivax],"]")
          }
        }
      }
      #
      #####################

      
      cex=0.8
      for(ievent_death in c("events","deaths"))
        for(iplot_type in c(1,2,4)){
          ###################### the first and second types of plot:  ########################
          #
          xlab_extra <- ""
          if(iplot_type%in%c(1,2)) { data_plot$time_y_date <- data_plot$time_y1; limits<-c(limits_x,limits_y1); ylab <- ""                            }
          if(iplot_type%in%c(3,4)) { data_plot$time_y_date <- data_plot$time_y2; limits<-c(limits_x,limits_y2); ylab <- paste0('date of  "', trimws(all_vax_names[ivax]),'"') }
          data_plot$time_y        <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + data_plot$time_y_date        / 365.25
          
          if( nrow(data_with_deaths)>0){
            if(iplot_type%in%c(1,2)) { data_with_deaths$time_y_date <- data_with_deaths$time_y1 }
            if(iplot_type%in%c(3,4)) { data_with_deaths$time_y_date <- data_with_deaths$time_y2 }
            data_with_deaths$time_y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + data_with_deaths$time_y_date / 365.25
          }
          
          if(ievent_death=="events" & iplot_type%in%c(2,4)){ 
            limits[1:2] <- c( min( -90, data_plot$time_x[  data_plot[,paste0(vax_time,"_v1")]-90 <= data_plot[,event_time]]), 90 )
            xlab_extra <- paste0("\n[between vax1 - 90days and ", trimws(all_vax_names[ivax])," + 60 days]")
          }
          if(ievent_death=="deaths") {
            limits[1] <- -2
            if(iplot_type%in%c(2,4)) { limits[2] <- min(180,limits[2]); xlab_extra <- paste0("\n[after ", trimws(all_vax_names[ivax])," + 180 days]") }
          }
          
          
          if(ievent_death == "events"){
            cond <- data_plot$cond_sampled_ids & data_plot$ivax_cond &  ( data_plot$time_x>=0 |  data_plot[,event_time] < data_plot[,paste0(vax_time,"_v1")] ) 
            if(!any(cond)) next
            distr_2d <- with(data_plot[cond,], kde2d(x=time_x, y=time_y,  h=20, n=100, lims=limits ) )
            distr_2d$y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + distr_2d$y/365.25
          }
          if(ievent_death == "deaths" & nrow(data_with_deaths)>0){ 
            cond <- data_with_deaths$cond_sampled_ids & data_with_deaths$ivax_cond &  data_with_deaths$time_x>=0 
            if(!any(cond)) next
            distr_2d <- with(data_with_deaths[cond,], kde2d(x=time_x, y=time_y,  h=20, n=100, lims=limits ) )
            distr_2d$y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + distr_2d$y/365.25
          }
          
          image_start <- function(data_for_plot){
            if(length(unique(distr_2d$y))>1) 
              image(distr_2d, ylab=ylab, xlab=paste0(event,': number of days from "',trimws(all_vax_names)[ivax],'"', xlab_extra),axes=F )
            #image(distr_2d, ylab=ylab, xlab=paste0(event,': number of days from "',trimws(all_vax_names)[ivax],'"', xlab_extra),col=hcl.colors(100, "YlOrRd", rev = TRUE),axes=F )
            else with(data_for_plot[data_for_plot$cond_sampled_ids,],
                      plot(time_x,time_y, type="n", xlim=c(min(0,time_x-abs(0.02*time_x),na.rm=T),max(60,time_x*1.02,na.rm=T)), axes=F, 
                           xlab=paste0(event,': number of days from "',trimws(all_vax_names)[ivax],'"', xlab_extra), ylab=ylab )  )
            
            axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
            #points_legend("event_vax1_num", "event_num", data, place=place, cex=cex )
            #title(paste0("id's with ",event," & ",tit_ivax,'\nrelative to "', all_vax_names[ivax],'"'))
            
            ################
            # ablines:
            abline(v=c(0,7,14,21,28, 60)); abline(v=0,lwd=2)
            if(ivax==1) abline(v=c(-30,-90))
            abline(h=c(2021 + 0.25*((-1):10)), col="pink")
            abline(h=c(2021,2022), col="pink",lwd=2)
          } # end function 'image_start'
          
          # add points:
          if(ievent_death=="events"){
            if( iplot_type%in%c(1,2) ) ylab <- paste0("date of ", event)
            image_start(data_plot[data_plot$cond_sampled_ids,])
            if(tit!="") tit_plot <- paste0(tit,"\n",tit_ivax_events)
            else        tit_plot <- tit_ivax_events
            title(paste0("id's with ",event," & ",tit_plot))
            if(any(data_plot$cond_sampled_ids)){
              points(data_plot$time_x[data_plot$cond_sampled_ids],  data_plot$time_y[data_plot$cond_sampled_ids], cex=cex )
              for(ivax2 in 1:length(all_vax_names)){
                cond_ivax2 <- data_plot$cond_sampled_ids & data_plot[,vax_name]==all_vax_names[ivax2]
                with(data_plot[cond_ivax2 & data_plot[,vax_time]<=data_plot[,event_time],], points(time_x, time_y, cex=cex, col=vax_names_col[ivax2] ) )
              }
            }
          }
          if(ievent_death=="deaths"){  
            if( iplot_type%in%c(1,2) ) ylab <- "date of death"
            image_start(data_with_deaths[data_with_deaths$cond_sampled_ids,])
            if(tit!="") tit_plot <- paste0(tit,"\n",tit_ivax_deaths)
            else        tit_plot <- tit_ivax_deaths
            title(paste0("id's with deaths & ",tit_plot))
            for(ivax2 in ivax:length(all_vax_names)){
              # deaths:
              if(nrow(data_with_deaths)>0){ 
                if( length(vaxdep_values)==1 & vaxdep_values[1]=="all" ) data_with_deaths$vaxdep_all_cond <- rep(T,nrow(data_with_deaths))
                else                                                     data_with_deaths$vaxdep_all_cond <- data_with_deaths$ivax_cond 
                #data_with_deaths$vaxdep_all_cond <- data_with_deaths$vaxdep_all_cond & data_with_deaths$ivax_id
                
                if(!any(data_with_deaths$vaxdep_all_cond & data_with_deaths$cond_sampled_ids)) next
                
                data_tmp <- data_with_deaths[data_with_deaths$vaxdep_all_cond & data_with_deaths$cond_sampled_ids,]
                if(nrow(data_tmp)>0){
                  if(iplot_type %in% c(1,2)) data_tmp$time_y_date <- data_tmp$time_y1
                  if(iplot_type %in% c(3,4)) data_tmp$time_y_date <- data_tmp$time_y2
                  data_tmp$time_y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + data_tmp$time_y_date/365.25
                  with(data_tmp[data_tmp[,vax_time]<=data_tmp$death_days &  data_tmp[,vax_name]==all_vax_names[ivax2],], 
                       points(time_x, time_y, cex=cex, col=vax_names_col[ivax2], pch=3 ) )
                }
              }
            } 
          } # end if ievent_death=="deaths"
        } # for iplot_type
    } # for ivax
  }  # for idep
  
  if(is.na(strata_value))   strata_value <- ""
  else if(strata_value!="") strata_value <- paste0(strata_value,"; ") 
  
  if(nrow(data_unvax_events)>0){
    ###########
    # sample points if they are more than 'max_n_points'
    #
    tit_plot <- ""
    data_unvax_events$cond_sampled_ids <- rep(T,nrow(data_unvax_events))
    if(!is.na(max_n_points)){ 
      # for events:
      data_unvax_events_all_ids <- unique(data_unvax_events[,id])
      if( max_n_points < length(data_unvax_events_all_ids) ) {
        sampled_ids <- sample(data_unvax_events_all_ids, max_n_points )
        data_unvax_events$cond_sampled_ids <- data_unvax_events[,id] %in% sampled_ids 
        tit_ivax_events <- paste0("[",round(100*max_n_points/length(data_unvax_events_all_ids),1),"% of all unvac.events (",length(data_unvax_events_all_ids),")]")
        if(tit!="") tit_plot <- paste0(tit,"\n",tit_ivax_events)
      }
    }
    data_unvax_events <- data_unvax_events[data_unvax_events$cond_sampled_ids, ]
    data_unvax_events$time_y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + data_unvax_events[,event_time]/365.25
    plot(1:nrow(data_unvax_events), data_unvax_events$time_y, type="n", axes=F, xlab="just observation numbers", ylab=paste0("date of ", event ) )
    axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
    abline(h=c(2021 + 0.25*((-1):10)), col="pink")
    abline(h=c(2021,2022), col="pink",lwd=2)
    points(1:nrow(data_unvax_events), data_unvax_events$time_y, cex=cex, col=vax_names_col[1], pch=1 )
  } else { plot(1,1,axes=F,type="n",xlab="",ylab=""); text(1,1,paste0("no ",event), cex=1.5) }
  title(paste0(strata_value,"unvaccinated persons (n_unvax=",n_unvaccinated,", n_vax=",n_vaccinated,")\n",tit_plot))
  
  no_plot_yet <- T; tit_plot <- ""
  if(nrow(data_with_deaths)>0){
    if(any(is.na(data_with_deaths[,vax_time]) & !is.na(data_with_deaths$death_days) )){
      ###########
      # sample points if they are more than 'max_n_points'
      #
      data_with_deaths$cond_sampled_ids <- rep(T,nrow(data_with_deaths))
      if(!is.na(max_n_points)){  
        # for deaths:
        data_with_deaths$ivax_cond  <- is.na(data_with_deaths[,vax_time]) & !is.na(data_with_deaths$death_days)
        ivax_ids <- length(unique(data_with_deaths[data_with_deaths$ivax_cond,id]))
        if( max_n_points < ivax_ids  ) {
          ivax_total_ids <- unique(data_with_deaths[,id])
          sampled_ids <- sample(ivax_total_ids, max_n_points/ivax_ids * total_deaths )
          data_with_deaths$cond_sampled_ids <- data_with_deaths[,id] %in% sampled_ids
          tit_ivax_deaths <- paste0("[",round(100*max_n_points/ivax_ids,1),"% of total unvac.deaths (",total_deaths,")]")
          if(tit!="") tit_plot <- paste0(tit,"\n",tit_ivax_deaths)
        }
      }
      data_with_deaths <- data_with_deaths[ data_with_deaths$cond_sampled_ids, ]
      data_with_deaths$time_y <- (1970 + as.numeric(as.Date("2020-09-01"))/365.25) + data_with_deaths$death_days/365.25
      plot( 1:sum(is.na(data_with_deaths[,vax_time])), data_with_deaths[is.na(data_with_deaths[,vax_time]),"time_y"], type="n", axes=F, xlab="just observation numbers",  ylab=paste0("date of death")  )
      axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
      abline(h=c(2021 + 0.25*((-1):10)), col="pink")
      abline(h=c(2021,2022), col="pink",lwd=2)
      points( 1:sum(is.na(data_with_deaths[,vax_time])), data_with_deaths[is.na(data_with_deaths[,vax_time]),"time_y"], cex=cex, col=vax_names_col[1], pch=3 )
      no_plot_yet <- F
    }
  }
  if(no_plot_yet){ plot(1,1,axes=F,type="n",xlab="",ylab=""); text(1,1,"no deaths",cex=1.5) }
  title(paste0(strata_value,"unvaccinated persons (n_unvax=",n_unvaccinated,", n_vax=",n_vaccinated,")\n", tit_plot))
  
  
}  # the end of function 'brand_images'


###################################

# create forest plots 
#   and 
# print markdown tables

forest_plots_tab <- function(tab_list, ndigits=2, nrows_forest_plot=40, cex=0.8, cex_head=0.8, 
                             col=c("red", "green3", "orange",  "deepskyblue", "magenta2", "cyan2", "chocolate1", "gray" ), 
                             lplot=T, ltable=T){
  
  if(is.null(tab_list)) return()
  if(!is.list(tab_list)) tab_list <- list(tab_list)
  if(all(sapply(tab_list,is.null))) return()
  if(is.null(names(tab_list))) names(tab_list) <- 1:length(tab_list)
  if(length(col)<length(tab_list))col <- rep(col,length(tab_list))[1:length(tab_list)]
  
  i_without_time <- 1
  if(any(substring(names(tab_list),1,6)=="no_adj"))
    i_without_time <- (1:length(tab_list))[substring(names(tab_list),1,6)=="no_adj"][1]
  nrow_without_time <- nrow(tab_list[[i_without_time]])
  
  
  RR_max <- NA
  
  for(i in 1:length(tab_list)){
    
    if(is.null(tab_list[[i]])) next
    
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
    
    if(no_RR) var_RR <- "RR_data"
    else {
      var_RR   <- "RR"  
      tt$RR_CI <- paste0(format(tt$RR,nsmall=ndigits)," (",format(tt$lci,nsmall=ndigits),";",format(tt$uci,nsmall=ndigits),")")
      tt$RR_CI[is.na(tt$RR)] <- ""
    }
    tab_list[[i]] <- tt
    if(any(!is.na(tt[,var_RR]))) RR_max <- max(tt[,var_RR], RR_max, na.rm=T)
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
      
      cond_group <- !is.na(tt[,var_RR]) & tt$group == groups[igr]
      if(sum(cond_group)==0) next
      
      group_part_start_row <- 1
      group_part_end_row   <- sum(cond_group)
      
      while(T){
        if(   start_y + (group_part_end_row-group_part_start_row+1)   > 0.9*nrows_forest_plot | start_y > nrows_forest_plot ){ start_y <- 1; par(new=F) }
        if( ( start_y + (group_part_end_row-group_part_start_row+1) ) > nrows_forest_plot                                   )  group_part_end_row <- nrows_forest_plot - start_y + group_part_start_row
        
        IRR_tit <- format("IRR (CI)",width=15)
        if(no_RR){
          IRR_tit <- format("IRR (CI*)",width=15)
          tt$lci <- pmax(0, tt[,var_RR] - 0.1 )
          tt$uci <- (tt[,var_RR] + 0.1)
          IRR_tit <- "IRR from data"
        }
        if(any(!is.na(tt[cond_group,][group_part_start_row:group_part_end_row,var_RR])))
          with(tt[cond_group,][group_part_start_row:group_part_end_row,], { 
            forest.default(slab=all_cat, x=get(var_RR), ci.lb=lci, ci.ub=uci,
                           refline=1, cex=cex, 
                           header = c("Intervals",IRR_tit),
                           xlim=c(-(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), max(1,RR_max)*1.6), 
                           ylim=c(2, -nrows_forest_plot) ,
                           alim=c(  0, max(1,RR_max)*1.1),
                           ilab=cbind(events_rw,events_ref, days_per_id_rw, days_per_id_ref), 
                           ilab.xpos= -(  (4:1)*10/(40/max(1,RR_max)) ), rows= -start_y
            )} )
        par(new=T) 
        
        legend_names <- names(tab_list)
        for(itab in 1:length(tab_list) ){
          if(itab!=i_without_time & !identical( tab_list[[itab]]$all_cat[1:nrow_without_time] , tab_list[[i_without_time]]$all_cat ) ) {
            warning(paste0("The main effect names of tables '",names(tab_list)[i_without_time]," and '",names(tab_list)[itab],"' are not identical"))
            legend_names <- legend_names[-itab]
            next
          }
          if(!no_RR & any(!is.na(tab_list[[itab]][1:nrow_without_time,][cond_group,][group_part_start_row:group_part_end_row,var_RR])))
            with(tab_list[[itab]][1:nrow_without_time,][cond_group,][group_part_start_row:group_part_end_row,], { 
              forest.default(slab=rep("",length(all_cat)), x=get(var_RR), ci.lb=lci, ci.ub=uci,
                             header = F, annotate=F, cex=cex, refline=1, 
                             xlim=c(-(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), max(1,RR_max)*1.6), 
                             ylim=c(2, -nrows_forest_plot) ,
                             alim=c(  0, max(1,RR_max)*1.1),
                             rows= -start_y + 0.8*(itab-1)/length(tab_list), col=col[itab]
              )} )
          par(new=T) 
        }
        
        if(start_y==1) add_text()
        legend( -(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), -1.04*nrows_forest_plot,legend=legend_names,
                box.col="lightgray", col=col[match(legend_names,names(tab_list))],xpd=T,lty=1,pch=15, cex=cex, xjust=0,yjust=1,ncol=max(1,round(length(legend_names)/3+0.1)))
        #legend( "bottomleft", inset=c(0,-0.05),legend=rep(legend_names,4),col=rep(col[match(legend_names,names(tab_list))],4),xpd=T,lty=1,pch=15, cex=cex, xjust=1,yjust=0,ncol=round(length(legend_names)/4+0.1))
        if(no_RR) text( -(max(nchar(tt$all_cat))/40*max(1,RR_max)+max(1,RR_max)), -nrows_forest_plot, "* CI's are not correct!!! There are no estimates from Poisson regression!!!")
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
    
    res$data_parameters$vax_time <- "vax_time"
    res$data_parameters$vax_name <- "vax_name"
    res$data_parameters$vax_date <- "vax_date"
    
    # specify the id variable, 'id':
    tmp <- tmp0 <- c("id","person_id","pat","pid" )
    if(any(names(param)=="id")) tmp <- param[["id"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'id' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$id <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_time variable, 'vax_time':
    tmp <- tmp0 <- c("vax_time","time_vax",  "vax_t","t_vax", "days_vax")
    if(any(names(param)=="vax_time")) tmp <- param[["vax_time"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_time' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_time <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_date variable, 'vax_date':
    tmp <- tmp0 <- c("vax_date","date_vax", "vax_time","time_vax",   "vax_t","t_vax", "days_vax")
    if(any(names(param)=="vax_date")) tmp <- param[["vax_date"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_date' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_date <- names(data)[tolower(names(data)) %in% tmp]
    
    # specify the vax_name variable, 'vax_name':
    tmp <- tmp0 <- c("vax_name","name_vax")
    if(any(names(param)=="vax_name")) tmp <- param[["vax_name"]]
    if(sum(tolower(names(data)) %in% tmp)!=1) stop(paste0("Specify the 'vax_name' variable (default,'",paste0(tmp0,collapse="' or '"),"')."))
    else res$data_parameters$vax_name <- names(data)[tolower(names(data)) %in% tmp]
    
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
                        cond_vax_var, cond_vax_values, vax_names_short, 
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
      cond_vax_values <- c( cond_vax_values[1], cond_vax_values )  # the first value for reference period ==> condition tov the start of the first vax
      
      if(missing(vax_names_short)){ 
        vax_names_short <- paste0(substring(cond_vax_values[1],1,1),"0") 
        for(ichar in  substring(cond_vax_values[-1],1,1)[!duplicated(substring(cond_vax_values[-1],1,1))] ) 
          vax_names_short <- c( vax_names_short, paste0( ichar, 1:sum(substring(cond_vax_values[-1],1,1)==ichar) )  )
      }
      else {
        if(mode(vax_names_short)=="character" & length(vax_names_short)==1){
          if( !any(vax_names_short %in% names(data))) stop(paste0("variable '",vax_names_short,"' not found in dataset 'data'."))
          vax_names_short <- data[ !duplicated(cond_vax_values_all), vax_names_short]
        }
      }      
      
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
  obj$rws$t0         <- obj$data_parameters$vax_time
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
  rws_text <- paste0(rws_text, ', t0="',               obj$rws$t0,                                       '"'          )
  rws_text <- paste0(rws_text, ', cond_var="',         obj$rws$cond_var,                                 '"'          )
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
    rws_text <- paste0(rws_text, ', t0="',               obj$rws$t0,                                               '"'          )
    rws_text <- paste0(rws_text, ', cond_var="',         obj$rws$cond_var,                                         '"'          )
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
    cat("\nThe data parameters:\n")
    for(ipar in 1:length(xx$data_parameters))
      print_scri_help_format(xx$data_parameters, ipar, text=paste0(" - the ",names(xx$data_parameters)[ipar]," variable: ") )
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
          cat(paste0(' risk windows = { ' , paste0( paste0('"',xx$rws$before$lab,'" = [',xx$rws$before$cut_points,";",c(xx$rws$before$cut_points[-1],0)-1,"]"), collapse=", ")," }\n\t "))
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
                                                             c(xx$rws$after$cut_points[[ivax]][-1]-1)), collapse="], ["),"] }; "))
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

characteristics <- function(data, event, path_file_name, condition_value="", vax_name="vax_number", id="person_id", start_obs="study_entry_days",
                            death_date="date_of_death", vax_time="vax_days", vax_date="vax_date", age="age_at_study_entry", lab_orders){
  
  if(!missing(path_file_name)) {
    if( substring(path_file_name,nchar(path_file_name)-3) != ".txt") path_file_name <- paste0(path_file_name,".txt")
    sink(paste0(path_file_name))
  }
  
  cat( paste("\n",paste( paste0(c('vax_name', 'id\t', 'start_obs'),'\t=\t"',c(vax_name, id, start_obs),'"') ,collapse = '\n')))
  if(condition_value!="") cat(paste0("\n\n\n\tConditoin name: ", condition_value, ":\n\n\n"))
  
  data$vax_name <- data[,vax_name]
  data$vax_time <- data[,vax_time]
  data$vax_date <- data[,vax_date]
  
  data$event_days <- data[,paste0(event,"_days")]
  data$event_date <- data[,paste0(event,"_date")]
  data$eventc <- paste0("event:",data[,event])
  
  # create some strata variables:
  data$age           <- data[,age]
  data$age_cat_30_50 <- factor_ref( paste0("age",as.character(cut(data[,age], c(-1,30,50, Inf)))), lab_orders = lab_orders )
  data$age_cat_30    <- factor_ref( paste0("age",as.character(cut(data[,age], c(-1,30,    Inf)))), lab_orders = lab_orders )
  data$sexc          <- factor_ref( paste0("sex:",data$sex),                                       lab_orders = lab_orders )
  data$sex_age       <- factor_ref( paste0("age",as.character(cut(data[,age], c(-1,30,    Inf))), " sex:",data$sex), lab_orders = lab_orders )
  
  flowchart_all <- vector("list",length=7)
  names(flowchart_all) <- c("all_data",
                            "with_vax_event",
                            "observed_before_vax1-90_and_events_after_vax1-90",
                            "observed_before_vax1_and_events_after_vax1-90",
                            "observed_before_vax1_and_events_after_vax1",
                            "observed_before_vax1-90_and_events_during_28_days_after_vax",
                            "observed_before_vax1_and_events_during_28_days_after_vax")
  # create the 'vax1_days' variable:
  vax1_days <- data[ data$vax_n==1, c(id, vax_time) ]
  names(vax1_days) <- c(id,"vax1_days")
  data <- merge.data.frame(data, vax1_days, by=id, all.x=T)
  
  data0 <- data
  for( i in 0:6 ){
    
    cat("\n\n\n\n*************************************************************************************\n")
    cat("*************************************************************************************\n")
    cat("*************************************************************************************\n\n")
    
    ####
    # from here only for vaccinated with event:
    if(i==0){ 
      data        <- data0
      data_deaths <- data0[!is.na(data0$death_days),]
      cat("\t\t\tThe whole dataset\n\n")
    }
    ####
    # from here only for vaccinated with event:
    if(i==1) {
      data        <- data0[data0[,event]==1         & !is.na(data0$vax_date),]
      data_deaths <- data0[!is.na(data0$death_days) ,]
      cat("\t\t\tvaccinated persons with ",event,"\n\n")
    }
    
    ############### 
    # from here only with event at -91 days before vax1:
    if(i==2){
      data        <- data0[data0[,event]==1 & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days-90 & 
                             data0$vax1_days-90 <= data0[,"event_days"], ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days-90 & 
                             data0$vax1_days-90 <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at (vax_day1 - 90 days) with ",event," after (vax_day1 - 90 days)\n\n")
    }
    ############### 
    # from here only with event at -91 days before vax1: vaccinated persons observed at vax_day1 with ",event," after (vax_day1 - 90 days)
    if(i==3){
      data        <- data0[data0[,event]==1 & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days-90 <= data0[,"event_days"], ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days-90 <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at vax_day1 with ",event," after (vax_day1 - 90 days)\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," after (vax_day1 - 0 days)
    if(i==4) {
      data        <- data0[data0[,event]==1 & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days <= data0[,"event_days"], ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at the first dose with ",event," after (vax_day1 - 0 days)\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," in [vax_date; vax_date + 28 days]
    if(i==5) {
      data        <- data0[data0[,event]==1 & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days-90 & 
                             data0$vax1_days <= data0[,"event_days"] &
                             data0[,"event_days"] <= data0[,vax_time] + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days-90 & 
                             data0$vax1_days <= data0$death_days &
                             data0$death_days <= data0[,vax_time] + 28 , ]
      cat("\t\t\tvaccinated persons observed at (first dose - 90 days) with ",event," in [vax_date; vax_date + 28 days]\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," in [vax_date; vax_date + 28 days]
    if(i==6) {
      data        <- data0[data0[,event]==1 & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days <= data0[,"event_days"] &
                             data0[,"event_days"] <= data0[,vax_time] + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax1_days) &   
                             data0[,start_obs] <= data0$vax1_days & 
                             data0$vax1_days <= data0$death_days &
                             data0$death_days <= data0[,vax_time] + 28 , ]
      cat("\t\t\tvaccinated persons observed at the first dose with ",event," in [vax_date; vax_date + 28 days]\n\n")
    }
    
    for(istrata_var in c("all", "age_cat_30_50","age_cat_30","sexc","sex_age")){
      cat("\n\n*******************************************************\n")
      cat("*******************************************************\n")
      cat("*******************************************************\n")
      
      if(istrata_var =="all"){
        cat(paste0("\n\tSummary not stratified:\n\n\n"))
        if(nrow(data       )>0) data$strata_variable        <- "all"
        if(nrow(data_deaths)>0) data_deaths$strata_variable <- "all"
      }
      else {
        cat(paste0("\n\tSummary stratified by variable '",istrata_var,"':\n\n\n"))
        if(nrow(data       )>0) data$strata_variable        <- data[       ,istrata_var]
        if(nrow(data_deaths)>0) data_deaths$strata_variable <- data_deaths[,istrata_var]
      }
      
      data$strata_cond <- !is.na(data$strata_variable) & !is.na(data[,id])
      
      if(nrow(data)>0 & any(data$strata_cond)){
        
        flowchart <- list()
        flowchart <- c(flowchart, n_ids = list(table1( unique(data[data$strata_cond,c("strata_variable",id)])[,"strata_variable"] ) ))
        
        ####
        # the whole dataset
        ##########  for all persons in the dataset:
        if(i==0) cat(paste0(" ***\tall persons in the dataset:***\n\nthe number of persons in the dataset: \n" ))
        
        ####
        # from here only for vaccinated with event:
        if(i==1) cat(paste0("***\tonly for vaccinated persons with ",event,":***\n\nthe number of persons: \n" ))
        
        ############### 
        # from here only with event at -91 days before vax1:
        if(i==2) cat(paste0("***\tonly for vaccinated persons observed at (vax_day1 - 90 days) with ",event," after (vax_day1 - 90 days):***\n\nthe number of persons: \n" ))
        
        ############### 
        # from here only with event at -91 days before vax1:
        if(i==3) cat(paste0("***\tonly for vaccinated persons observed at vax_day1 with ",event," after (vax_day1 - 90 days):***\n\nthe number of persons: \n" ))
        
        ############### 
        # from here only 
        if(i==4) cat(paste0("***\tonly for vaccinated persons observed at the first dose with ",event," after (vax_day1 - 0 days):***\n\nthe number of persons: \n" ))
        
        ############### 
        # from here only 
        if(i==5) cat(paste0("***\tonly for vaccinated persons observed at (first dose - 90 days) with ",event," in [vax_date; vax_date + 28 days]:***\n\nthe number of persons: \n" ))
        
        ############### 
        # from here only 
        if(i==6) cat(paste0("***\tonly for vaccinated persons observed at the first dose with ",event," in [vax_date; vax_date + 28 days]:***\n\nthe number of persons: \n" ))
        
        print(flowchart$n_ids)
        
        cat(paste0("\nthe numbers for persons  per ",event,", vaccine name and brand:\n"))
        flowchart <- c(flowchart, n_ids_per_event_vax_name_brand = list(table1( unique(data[,c("strata_variable",id,"eventc","vax_name","vax_brand")]) [ , c("strata_variable","eventc","vax_name","vax_brand")])) )
        print(flowchart$n_ids_per_event_vax_name_brand)
        
        
        # age as continuous or integer:
        if(any(!is.na(data$age))){
          cat(paste0("\n\nthe distribution of variable '",age,"':\n\n"))
          flowchart <- c(flowchart, summary_id_age                      = list(do.call("rbind",with(unique(data[!is.na(data$age) & data$strata_cond, c("strata_variable",id,"age")]), tapply( age, strata_variable, function(x)c(summary(x),n=length(x))) )) ))
          flowchart <- c(flowchart, summary_id_age_vax_name             = list(do.call("rbind",with(unique(data[!is.na(data$age) & data$strata_cond, c("strata_variable",id,"age",         "vax_name"            )]), tapply(age, paste(strata_variable,       vax_name          ), function(x)c(summary(x),n=length(x))) )) ))
          flowchart <- c(flowchart, summary_id_age_vax_brand            = list(do.call("rbind",with(unique(data[!is.na(data$age) & data$strata_cond, c("strata_variable",id,"age",                    "vax_brand")]), tapply(age, paste(strata_variable,                vax_brand), function(x)c(summary(x),n=length(x))) )) ))
          flowchart <- c(flowchart, summary_id_age_event                = list(do.call("rbind",with(unique(data[!is.na(data$age) & data$strata_cond, c("strata_variable",id,"age","eventc"                       )]), tapply(age, paste(strata_variable,eventc                   ), function(x)c(summary(x),n=length(x))) )) ))
          flowchart <- c(flowchart, summary_id_age_event_vax_name_brand = list(do.call("rbind",with(unique(data[!is.na(data$age) & data$strata_cond, c("strata_variable",id,"age","eventc","vax_name","vax_brand")]), tapply(age, paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x))) )) ))
          
          print( flowchart$summary_id_age                      ); cat("\n")
          print( flowchart$summary_id_age_vax_name             ); cat("\n")
          print( flowchart$summary_id_age_vax_brand            ); cat("\n")
          print( flowchart$summary_id_age_event                ); cat("\n")
          print( flowchart$summary_id_age_event_vax_name_brand )  
        } else cat("\n\nno age.\n\n")
        
        
        # summaries for vaccination time and date:
        if(any(!is.na(data[,vax_time]))){
          
          cat(paste0("\n\nthe distribution of the vaccinatioen variable: '",vax_time,"':\n"))
          flowchart <- c( flowchart, summary_id_vax_time = list( do.call("rbind",with(
            unique(data[!is.na(data[,vax_time]) & data$strata_cond, c("strata_variable","vax_n",id, vax_time)]), 
            tapply( get(vax_time), paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
          print(flowchart$summary_id_vax_time)
          
          cat(paste0("\nthe distribution of the vaccination variable: '",vax_date,"':\n"))
          flowchart <- c(flowchart, summary_id_vax_date = list(do.call("rbind",with(
            unique(data[!is.na(data[,"vax_date"]) & data$strata_cond, c("strata_variable","vax_n",id, "vax_date")]), 
            tapply( vax_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
          print(flowchart$summary_id_vax_date)
          
        } else cat("\n\nno vax.\n\n")
        
        
        
        # event_days:
        if(any(!is.na(data$event_days))){
          
          if(any( (cond<-!is.na(data[,vax_time]) & !is.na(data$event_days) &  (data$event_days-data[,vax_time])<0) )){
            cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) before vaccination:\n"))
            flowchart <- c( summary_id_event_min_vax_before_vax = list( do.call("rbind",with(
              unique(data[cond & data$strata_cond, c("strata_variable",id, "event_days",vax_time,"vax_name")]), 
              tapply( (event_days-get(vax_time)), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
            print(flowchart$summary_id_event_min_vax_before_vax)
          }
          
          if(any( (cond<-!is.na(data[,vax_time]) & !is.na(data$event_days) &  (data$event_days-data[,vax_time])>=0) )){
            cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) after vaccination:\n"))
            flowchart <- c( summary_id_event_min_vax_after_vax = list( do.call("rbind",with(
              unique(data[cond & data$strata_cond, c("strata_variable",id, "event_days",vax_time,"vax_name")]), 
              tapply( (event_days-get(vax_time)), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
            print(flowchart$summary_id_event_min_vax_after_vax)
          }
          cat(paste0("\n\nthe distribution of the '",event,"_days' variable:\n"))
          flowchart <- c( flowchart, summary_id_event_time = list( do.call("rbind",with(
            unique(data[!is.na(data$event_days) & data$strata_cond, c("strata_variable","vax_n",id, "event_days")]), 
            tapply( event_days, paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
          print(flowchart$summary_id_event_time)
          
          cat(paste0("\n\nthe distribution of the '",event,"_date' variable:\n"))
          flowchart <- c(flowchart, summary_id_event_date = list(do.call("rbind",with(
            unique(data[!is.na(data$event_days) & data$strata_cond, c("strata_variable","vax_n",id, "event_date")]), 
            tapply( event_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
          print(flowchart$summary_id_event_date)
          
          if(any( (cond<-!is.na(data[,vax_time]) & !is.na(data$event_days) &  (data$event_days-data[,vax_time])<0) )){
            cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) before vaccination per ",event,", vaccine name and brand:\n"))
            flowchart <- c(flowchart, summary_id_event_min_vax_before_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
              unique(data[cond & data$strata_cond, c("strata_variable",id, "event_days",vax_time, "eventc","vax_name","vax_brand")]), 
              tapply( (event_days-get(vax_time)), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
            print(  flowchart$summary_id_event_min_vax_before_vax_per_event_vax_name_brand)
          }
          
          if(any( (cond<-!is.na(data[,vax_time]) & !is.na(data$event_days) &  (data$event_days-data[,vax_time])>=0) )){
            cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) after vaccination per ",event,", vaccine name and brand:\n"))
            flowchart <- c(flowchart, summary_id_event_min_vax_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
              unique(data[cond & data$strata_cond, c("strata_variable",id, "event_days",vax_time, "eventc","vax_name","vax_brand")]), 
              tapply( (event_days-get(vax_time)), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
            print(  flowchart$summary_id_event_min_vax_after_vax_per_event_vax_name_brand)
          }
        } else cat("\n\nno events.\n\n")
      }
      
      
      # death:
      data_deaths$strata_cond <- !is.na(data_deaths$strata_variable) & !is.na(data_deaths[,id])
      
      if(nrow(data_deaths)>0 & any(data_deaths$strata_cond)){
        
        data$strata_cond <- !is.na(data$strata_variable) & !is.na(data[,id])
        
        cat(paste0("\n\nthe distribution of 'death_days': days after vaccination:\n"))
        flowchart$deaths <- c( summary_id_death_after_vax = list( do.call("rbind",with(
          unique(data_deaths[!is.na(data_deaths[,vax_time]) & data_deaths$strata_cond, c("strata_variable",id, "death_days",vax_time,"vax_name")]), 
          tapply( (death_days-get(vax_time)), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
        print(flowchart$deaths$summary_id_death_after_vax)
        
        cat(paste0("\n\nthe distribution of the 'death_days' variable:\n"))
        flowchart$deaths <- c( flowchart$deaths, summary_id_death_time = list( do.call("rbind",with(
          unique(data_deaths[data_deaths$strata_cond, c("strata_variable","vax_n",id, "death_days")]), 
          tapply( death_days, paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
        print(flowchart$deaths$summary_id_death_time)
        
        cat(paste0("\n\nthe distribution of the 'death_date' variable after vaccination:\n"))
        flowchart$deaths <- c(flowchart$deaths, summary_id_death_date = list(do.call("rbind",with(
          unique(data_deaths[data_deaths$strata_cond, c("strata_variable","vax_n",id, death_date)]), 
          tapply( get(death_date), paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
        print(flowchart$deaths$summary_id_death_date)
        
        cat(paste0("\nthe distribution of the 'death_days' (days after vaccination) per ",event,", vaccine name and brand:\n"))
        flowchart$deaths <- c(flowchart$deaths, summary_id_death_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
          unique(data_deaths[!is.na(data_deaths[,vax_time]) & data_deaths$strata_cond, c("strata_variable",id, "death_days",vax_time, "eventc","vax_name","vax_brand")]), 
          tapply( (death_days-get(vax_time)), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
        print(  flowchart$deaths$summary_id_death_after_vax_per_event_vax_name_brand)
        
      } else cat("\n\nno deaths.\n\n")
      
      flowchart_all[[i+1]] <- c( flowchart_all[[i+1]], list(flowchart))
      names(flowchart_all[[i+1]])[length(flowchart_all[[i+1]])] <- istrata_var
      
    }  # end 'istrata_var'
  } # end 'i'
  
  attributes(flowchart_all) <- c( attributes(flowchart_all), variables=list( vax_name=vax_name, id=id, start_obs=start_obs, condition_value=condition_value )  )
  
  
  if(!missing(path_file_name))  sink()
  
  path_file_name <- paste0( substring(path_file_name,1,nchar(path_file_name)-3), "RData")
  if(!missing(path_file_name))  save( flowchart_all, file=path_file_name )
  
  invisible(flowchart_all)
  
} # the end of function 'characteristics'



#
#   the end of functions :
#  
###############################
