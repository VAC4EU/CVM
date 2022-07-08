
# Program Information  ----------------------------------------------------
#
#  functions for SCRI analysis 

# Functions:    scri_sv
#               refresh_event_variable
#               create_rws
#               split_intervals
#               summary_tab
#               factor_ref
#               combine_vars
#               plot_res
#               add_to_report_list
#               add_to_models_list
#               more functions .....
#
# Author:       Svetlana Belitser 
#               nov 2021 - feb 2022
#

  
scri_sv <- function(formula,    
                    event_time, event,
                    id,
                    rws,                          # list of risk/control windows definitions
                    time_seq, split_seq_name = "cal_time_cat", time_seq_ref="most events", 
                    time_dep = NA,                # list of risk/control windows definitions
                    combine_vars = c(),           # list of parameters to create new one variable from two other variables
                    start_obs, end_obs,
                    data, cond_var = "",
                    nvax = 2,
                    lab_orders = NA,
                    ref=1,
                    rw_observed_percentage=100,   # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                    censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                                                  #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                    event_in_rw=T,                # if event in rw ==> this rw should not be deleted even if not completely observed
                    delete_coef_no_events = T,
                    lprint    = T, test = F,
                    save_data = F
                    ){

  if(missing(formula   )) stop("'formula' is missing.")
  if(missing(rws       )) stop("'rws' is missing.")
  if(missing(event_time)) stop("'event_time' is missing.")
  if(missing(event     )) stop("'event' is missing.")
  if(missing(id        )) stop("'id' is missing.")
  if(missing(data      )) stop("'data' is missing.")
  if(missing(start_obs )) stop("'start_obs' is missing.")
  if(missing(end_obs   )) stop("'end_obs' is missing.")
  

  tb       <- attributes(terms(formula))$factor
  tab_vars <- dimnames(tb)[[2]]

  
  if(missing(event)      & any(rowSums(tb)==0))  event      <- dimnames(tb)[[1]][rowSums(tb)==0]
  if(missing(event_time) & any(rowSums(tb)==0))  event_time <- dimnames(tb)[[1]][rowSums(tb)==0]

 
  ###########################################################################
  ################# create_rws (v3)  ###############
  data_rws  <- create_rws(
    deparse(rws),
    start_obs = start_obs, end_obs = end_obs,
    event_time = event_time, event = event,
    id = id,
    lab_orders = lab_orders,
    data =data,
    ref="pre-",   #  ref=5 OR ref= "pre-exposure [-90;-30]"  or a part of a ref.category name
    rw_observed_percentage = rw_observed_percentage,   
    censored_vars = censored_vars,           
    event_in_rw = event_in_rw               
  ) 

  if(cond_var!="") {
    tmp <- levels(data_rws[,"lab"])
    data_rws[,"lab"] <- as.character(data_rws[,"lab"])
    data_rws[ !data_rws[,cond_var], "lab"] <- "only_for_time_adj"
    #data_rws[,"lab"] <- factor(data_rws[,"lab"],levels=c(tmp,"only_for_time_adj"))
    data_rws[,"lab"] <- factor_ref(  data_rws[,"lab"], lab=c(tmp,"only_for_time_adj"),  ref="pre-")  
  }  
  data_rws <- refresh_event_variable( "rw_start", "rw_end", event, event_time, data_rws)


  #######
  #  create time intervals:
  #
  if(!missing(time_seq) & split_seq_name %in% tab_vars & nrow(data_rws)>0)
    data_rws <- split_intervals( data =data_rws, 
                                 start_interval = "rw_start", end_interval = "rw_end", 
                                 splits_names = split_seq_name, # "cal_time_cat",
                                 splits       = time_seq,
                                 ref          = time_seq_ref, #"most events", 
                                 event        = event,        #   myopericarditis 
                                 event_time   = event_time 
                                 )

  ####################
  #
  #  create tim-dependent variables.
  #
  #   For example, variable 'vax123' with values: "pre-" before vax1 ,"vax1" after vax1; "vax2" after vax2, ...
  #
  sep_vars <- list()
  if(!is.na(time_dep) & nrow(data_rws)>0){
    time_dep <- eval(time_dep,data_rws)
    #time_dep <- eval(parse(text=deparse(substitute(time_dep))),data_rws)
    for(itd in 1:length(time_dep)){
      ###
      # create time-dep variable (ex,'vax123'). For example,  "vax1" after vax1; vax2 after vax2, ...
      data_rws <- split_intervals( data =data_rws, 
                                   start_interval = "rw_start", end_interval = "rw_end", 
                                   splits_names = eval(time_dep[[itd]]$splits_names),    #    "cal_time_cat"   or "dose12"
                                   splits       = eval(time_dep[[itd]]$splits,data_rws), #  c("days_vax1","days_vax2"),
                                   lab          = eval(time_dep[[itd]]$lab),             #  c("pre-vax","dose 1","dose 2"),
                                   ref          = eval(time_dep[[itd]]$ref),             #    "pre-"  or  "most events", 
                                   event        = event,
                                   event_time   = event_time
      )
  
      if(!is.null(time_dep[[itd]]$change)){
        
        changes <- eval(time_dep[[itd]]$change)
        if( is.list(changes) & length(changes)>1 )
          if( !is.list(changes[[1]]) | !is.list(changes[[2]])  ) changes <- list(changes)
        
        for(ich in 1:length(changes)){
          
          #
          # copy the just created time-dependent variable into varibale with name 'data_rws[, changes$name ]'
          #             ex of changes$name: brand123
          data_rws[,changes[[ich]]$new_name] <- as.character(data_rws[,changes[[ich]]$from])
          #data_rws[,changes$name] <- as.character(data_rws[,time_dep[[itd]]$splits_names])
          for(ireplace in 1:length(changes[[ich]]$replace)){
            #
            # cond: if the copied variable (from just created time-dep variable) == specific value 'changes$replace[[ireplace]][["value"]]' ==> 
            #       ==> this value in this copied variable is changed to the value from variable  'changes$replace[[ireplace]][["var_name"]]' 
            #       [example: cond for "dose 2"  ==>  brand of vax2]
            #
            cond <- data_rws[,changes[[ich]]$condition_var] == changes[[ich]]$replace[[ireplace]][["value"]] 
            if(any(cond)){
              if( any(names(changes[[ich]])=="add_begin_sep") )
                data_rws[ cond, changes[[ich]]$new_name ] <- paste( as.character( data_rws[ cond,  changes[[ich]]$replace[[ireplace]][["var_name"]] ] ), 
                                                                    data_rws[ cond, changes[[ich]]$new_name      ], 
                                                                    sep = changes[[ich]]$add_begin_sep )
              if( any(names(changes[[ich]])=="add_end_sep") )
                data_rws[ cond, changes[[ich]]$new_name ] <- paste( data_rws[ cond, changes[[ich]]$new_name ],
                                                                    as.character( data_rws[ cond,  changes[[ich]]$replace[[ireplace]][["var_name"]] ] ), 
                                                                    sep = changes[[ich]]$add_end_sep )
              if( !any(names(changes[[ich]]) %in% c("add_begin_sep","add_end_sep")) )
                data_rws[ cond, changes[[ich]]$new_name ] <- as.character( data_rws[ cond,  changes[[ich]]$replace[[ireplace]][["var_name"]] ] )
            }  
          }  
          data_rws[,changes[[ich]]$new_name] <- factor_ref(  data_rws[,changes[[ich]]$new_name], lab_orders=lab_orders,  ref=ref, event_var=data_rws[,event] )  
          
          
          if(any(names(changes[[ich]])=="separate_variables")){
        
            sep_origin         <- paste0(changes[[ich]]$new_name,"_lab")
            sep_origin_formula <- paste0(changes[[ich]]$new_name,":lab")
            
            sep_vars <- vector("list", length=nvax)
            names(sep_vars) <- paste0("ref_dose",1:nvax)
            
            brands <- unique(as.character(data_rws[,changes[[ich]]$new_name]))
            brands <- brands[brands != "no_vax"]
            for(ibr in 1:length(brands)){
              
              # create vars: paste0(changes[[ich]]$new_name,ibr,"_123")
              var_name <- paste0(changes[[ich]]$new_name,ibr)
              data_rws[,var_name] <- as.character(data_rws[,changes[[ich]]$new_name])
              data_rws[ !(data_rws[,var_name] %in% c("no_vax",brands[ibr]) ), var_name] <- paste0("no_",brands[ibr])
              # select only patients with minimaal one Moderna or Pfizer or ...:
              data_rws[ data_rws[,var_name] %in% "no_vax"  & data_rws[,id] %in% unique(data_rws[data_rws[,var_name]==brands[ibr],id]) , var_name] <- brands[ibr]
              
              # create vars: paste0(changes[[ich]]$new_name,ibr,"_123_lab")
              var_name <- paste0(changes[[ich]]$new_name,ibr,"_lab")
              data_rws[,var_name] <- paste(data_rws[,paste0(changes[[ich]]$new_name,ibr)], data_rws$lab)
              data_rws[ substring(data_rws[,var_name],1,nchar(brands[ibr])) != brands[ibr], var_name] <- paste0("no_",brands[ibr])
              
              if(any(names(changes[[ich]])=="add_dist")){
                
                data_rws$tmp <- as.character(data_rws[,var_name])
                
                for(ireplace in 1:length(changes[[ich]]$add_dist)){
                  #
                  # cond: if the copied variable (from just created time-dep variable) == specific value 'changes$replace[[ireplace]][["value"]]' ==> 
                  #       ==> this value in this copied variable is changed to the value from variable  'changes$replace[[ireplace]][["var_name"]]' 
                  #       [example: cond for "dose 2"  ==>  brand of vax2]
                  #
                  cond <- data_rws[,changes[[ich]]$condition_var] == changes[[ich]]$add_dist[[ireplace]][["value"]] 
                  if(any(cond))
                    data_rws[ cond, "tmp" ] <- paste( data_rws[ cond, "tmp" ],
                                                      as.character( data_rws[ cond,  changes[[ich]]$add_dist[[ireplace]][["var_name"]] ] ), 
                                                      sep = " " )
                } 
                data_rws$tmp[ substring(data_rws$tmp, 1, nchar(paste0("no_",brands[ibr])) ) == paste0("no_",brands[ibr]) ] <- paste0("no_",brands[ibr])
                data_rws[,var_name] <- data_rws$tmp
              }
              
              data_rws[, var_name] <- factor_ref(  data_rws[, var_name], lab_orders=lab_orders,  ref="pre-" )  
              
              #######
              #
              # models with reference categories: pre-exposure dose 1
              #
              var_name <- paste0(changes[[ich]]$new_name,ibr,"_lab_refd1")
              data_rws[,var_name] <- as.character( data_rws[,paste0(changes[[ich]]$new_name,ibr,"_lab")] )
              data_rws[ substring(data_rws[,var_name],1,nchar(brands[ibr])) != brands[ibr], var_name ] <- paste0("no ",brands[ibr]," or dose 1")
              
              ids_tmp <- data_rws$person_id[  grepl("dose1",gsub(" ","",data_rws[, var_name])) &  
                                                grepl(brands[ibr],data_rws[, var_name], fixed=T) &  
                                                substring(data_rws[, var_name],1,3) != "no " ]
              cond <- data_rws$person_id %in% ids_tmp  
              data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & cond , var_name] <- paste0(brands[ibr]," dose 1 pre-exposure[-90;-30]")
              data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & !cond , var_name] <- paste0("no ",brands[ibr]," or dose 1")
              
              ref_tmp <- ""
              if( length(grep("dose 1 pre-",unique(data_rws[,var_name])))>0) ref_tmp <- "dose 1 pre-"
              else if( length(grep("buffer",unique(data_rws[,var_name])))>0) ref_tmp <- "buffer"
              
              if(ref_tmp=="") ref_tmp <- 1
              data_rws[, var_name] <- factor_ref(  data_rws[, var_name], lab_orders=lab_orders,  ref=ref_tmp )  
              
              if( ref_tmp != 1 & nlevels(data_rws[, var_name])>1 ) sep_vars$ref_dose1 <- c( sep_vars$ref_dose1 , var_name)
              
              
              #######
              #
              # models with reference categories: pre-exposure dose 2
              #
              var_name <- paste0(changes[[ich]]$new_name,ibr,"_lab_refd2")
              data_rws[,var_name] <- as.character( data_rws[,paste0(changes[[ich]]$new_name,ibr,"_lab")] )
              data_rws[ substring(data_rws[,var_name],1,nchar(brands[ibr])) != brands[ibr], var_name ] <- paste0("no ",brands[ibr]," or dose 2")
              
              ids_tmp <- data_rws$person_id[  grepl("dose2",gsub(" ","",data_rws[, var_name])) &   
                                                grepl(brands[ibr],data_rws[, var_name], fixed=T) &  
                                                substring(data_rws[, var_name],1,3) != "no " ]
              cond <- data_rws$person_id %in% ids_tmp      
              data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & cond , var_name] <- paste0(brands[ibr]," dose 2 pre-exposure[-90;-30]")
              data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & !cond , var_name] <- paste0("no ",brands[ibr]," or dose 2")
              
              ref_tmp <- ""
              if( length(grep("dose 2 pre-",unique(data_rws[,var_name])))>0) ref_tmp <- "dose 2 pre-"
              else if( length(grep("buffer",unique(data_rws[,var_name])))>0) ref_tmp <- "buffer"
              
              if(ref_tmp=="") ref_tmp <- 1
              data_rws[, var_name] <- factor_ref(  data_rws[, var_name], lab_orders=lab_orders,  ref=ref_tmp )  
              
              if( ref_tmp != 1 & nlevels(data_rws[, var_name])>1 ) sep_vars$ref_dose2 <- c( sep_vars$ref_dose2 , var_name)
              
              
              
              if(nvax>2){
                #######
                #
                # models with reference categories: pre-exposure dose 3
                #
                var_name <- paste0(changes[[ich]]$new_name,ibr,"_lab_refd3")
                data_rws[,var_name] <- as.character( data_rws[,paste0(changes[[ich]]$new_name,ibr,"_lab")] )
                data_rws[ substring(data_rws[,var_name],1,nchar(brands[ibr])) != brands[ibr], var_name ] <-paste0("no ",brands[ibr]," or dose 3")
                
                ids_tmp <- data_rws$person_id[  grepl("dose3",gsub(" ","",data_rws[, var_name])) &   
                                                  grepl(brands[ibr],data_rws[, var_name], fixed=T) &  
                                                  substring(data_rws[, var_name],1,3) != "no " ]
                cond <- data_rws$person_id %in% ids_tmp      
                data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & cond , var_name] <- paste0(brands[ibr]," dose 3 pre-exposure[-90;-30]")
                data_rws[   data_rws[, var_name] ==paste0(brands[ibr]," pre-exposure[-90;-30]")  & !cond , var_name] <- paste0("no ",brands[ibr]," or dose 3")
                
                ref_tmp <- ""
                if( length(grep("dose 3 pre-",unique(data_rws[,var_name])))>0) ref_tmp <- "dose 3 pre-"
                else if( length(grep("buffer",unique(data_rws[,var_name])))>0) ref_tmp <- "buffer"
                
                if(ref_tmp=="") ref_tmp <- 1
                data_rws[, var_name] <- factor_ref(  data_rws[, var_name], lab_orders=lab_orders,  ref=ref_tmp )  
                
                if( ref_tmp != 1 & nlevels(data_rws[, var_name])>1 ) sep_vars$ref_dose3 <- c( sep_vars$ref_dose3 , var_name)
              }
              
            }
          } 
        }  
      }
    }  
  }  
  

  if(length(combine_vars)>0 & nrow(data_rws)>0){
    data_rws$lab  <- combine_vars_func( data_rws[, c(combine_vars, "lab"), drop=F ], lab_orders = lab_orders, ref=ref, event = data_rws[,event] )
    data_rws$lab  <- factor_ref(  data_rws$lab, lab_orders=lab_orders, ref=ref, event_var=data_rws[,event] )  
  }
  
  
 
  #####
  #  create 'event'and 'interval' variables
  #
  data_rws <- refresh_event_variable( "rw_start", "rw_end", event, event_time, data_rws)
  
  data_rws$interval <- data_rws$rw_end - data_rws$rw_start + 1
  data_rws <- data_rws[data_rws$interval>0,]
  
  # delete id's without events in the windows:
  id_no_events <- names((tb<-tapply(data_rws[,event],data_rws[,id], sum))[tb==0])
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  data_rws <- data_rws[ !(as.character(data_rws[,id]) %in% id_no_events),  ]
  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  
  
#  # delete id's without events in the windows:
#  id_no_events <- as.numeric(names((tb<-tapply(data_rws[,event],data_rws[,id], sum))[tb==0]))
#  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
#  data_rws <- data_rws[ !(data_rws[,id] %in% id_no_events),  ]
#  sum(data_rws[,event]); length(unique(data_rws[,id])); nrow(data_rws)
  
  
  
  
  
    
  if(nrow(data_rws)==0) return(  list( res_tab = NULL, 
                                       model   = NULL,
                                       call    = list( match.call())  ))
  
 
  if(any(names(data_rws)==split_seq_name))
    data_rws[,split_seq_name]  <- factor_ref(  as.character(data_rws[,split_seq_name]), 
                                               lab=levels(data_rws[,split_seq_name])[levels(data_rws[,split_seq_name]) %in% unique(as.character(data_rws[,split_seq_name])) ],
                                               ref=time_seq_ref, 
                                               event_var=data_rws[,event] )  

  
  # combine variables in formula with risk windows variable 'lab' (the begining of the code):
  tb <- attributes(terms(formula))$factor
  if(!missing(time_seq) & split_seq_name %in% dimnames(tb)[[2]])
    tb <- tb[, c( dimnames(tb)[[2]][dimnames(tb)[[2]]!=split_seq_name], split_seq_name ) ]
 
  
  
  if(length(sep_vars)>0){

    sep_vars0 <- sep_vars
    tb0 <- tb
    if(any(cond<-(colSums(tb0)>1) & grepl(":",dimnames(tb0)[[2]]))  & nrow(data_rws)>0)
      dimnames(tb0)[[2]] <- gsub(" ","",dimnames(tb0)[[2]])
      tb <- tb0
      for(i in 1:sum(cond)){
        for(idose in 1:length(sep_vars0))
          for(ibr in 1:length(sep_vars0[[idose]])){
            vars_tmp <- gsub( sep_origin_formula, sep_vars0[[idose]][ibr], dimnames(tb0)[[2]][i]) 
            for(ivar in 1:length(vars_tmp)){
              if(length(strsplit( vars_tmp[ivar]  ,":")[[1]])>1)
                data_rws[,paste(strsplit( vars_tmp[ivar] ,":")[[1]],collapse="_")]  <- combine_vars_func( 
                  data_rws[, strsplit( vars_tmp[ivar]  ,":")[[1]], drop=F ], lab_orders = lab_orders, event = data_rws[,event] )  #ref=ref, 
              if( grepl(":",vars_tmp[ivar],fixed = TRUE) ){
                sep_vars[[idose]] <- c( sep_vars[[ idose ]], paste(strsplit( vars_tmp[ivar] ,":")[[1]],collapse="_") )   
                sep_vars[[idose]] <- sep_vars[[idose]][ !( sep_vars[[idose]] %in% strsplit( vars_tmp[ivar]  ,":")[[1]] ) ]
              }  
            }
        }
      }  
    
    extra_vars <- dimnames(tb0)[[2]][ !grepl( sep_origin_formula, dimnames(tb0)[[2]] ) ]
    if(length(extra_vars)) sep_vars <- lapply(sep_vars, function(x)c(x,extra_vars))  
    
    mods <- vector("list",length=length(sep_vars))
    names(mods) <- names(sep_vars)
    
    for(isep in 1:length(sep_vars)){
      
      var_names <- sep_vars[[isep]]
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
      suppressWarnings(
        mod <- try( clogit(formula = Poisson_formula, data = data_rws, control=coxph.control(iter.max=1000) ), silent = T )
      )  

      if(class(mod)[[1]]== "try-error"){

        if(F){
          if(nrow(table1(data_rws[  data_rws[,event]==1 ,var_names]))<=8)
            print( table1(data_rws[  data_rws[,event]==1 ,var_names]) ) 
          print(cbind.data.frame(var_names, nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names])) ))
          if(any(names(data_rws)==split_seq_name))   #"cal_time_cat"))
            print(cbind.data.frame(var_names[var_names!=split_seq_name], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!=split_seq_name]])) ))
            #print(cbind.data.frame(var_names[var_names!="cal_time_cat"], nrow=nrow(table1(data_rws[  data_rws[,event]==1 ,var_names[var_names!="cal_time_cat"]])) ))
          
          if(!missing(time_seq) )
               warning("Error in Poisson regression. \ntime_seq=",median(diff(time_seq),na.rm=T),"; Formula: ", deparse(Poisson_formula))
          else warning("Error in Poisson regression. \nno time_seq; Formula: ", deparse(Poisson_formula))
        }        
        
        if(isep==1)
          res_tab <- summary_tab( var_names=sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  #mod=mod, 
                                  lab_orders=lab_orders,
                                  delete_coef_no_events = delete_coef_no_events, 
                                  coef_cond = coef_cond, model_number = isep )
        else  
          res_tab <- summary_tab( var_names=sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  #mod=mod, 
                                  lab_orders=lab_orders,
                                  delete_coef_no_events=delete_coef_no_events, 
                                  coef_cond=coef_cond, model_number = isep,
                                  res_tab = res_tab )
        
        mods[[isep]] <- list( mod ) 
        
      }
      else{

        coef_cond <-  grepl( paste0( ibr, " buffer"), names(mod$coefficients) )
        for(ibr in brands)
          coef_cond <- coef_cond | grepl( paste0( ibr, " dose ", isep), names(mod$coefficients) )
        for(iextra in extra_vars)
          coef_cond <- coef_cond | ( substring(names(mod$coefficients),1,nchar(iextra))==iextra )
        
        if(isep==1)
          res_tab <- summary_tab( var_names=sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  delete_coef_no_events = delete_coef_no_events, 
                                  coef_cond = coef_cond, model_number = isep )
        else  
          res_tab <- summary_tab( var_names=sep_vars[[isep]], 
                                  event=event, 
                                  data=data_rws, 
                                  id_name=id,  
                                  mod=mod, 
                                  lab_orders=lab_orders,
                                  delete_coef_no_events=delete_coef_no_events, 
                                  coef_cond=coef_cond, model_number = isep,
                                  res_tab = res_tab )
        
        mods[[isep]] <- list(summary(mod)) 
      }
    }
    
    mod <- mods
    
    if(lprint)  print(format( res_tab, digits=2, justify="left" ))
    
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
      suppressWarnings(
        mod <- try( clogit(formula = Poisson_formula, data = data_rws, control=coxph.control(iter.max=1000) ), silent = T )
      )  

      if(class(mod)[[1]]== "try-error"){

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
        res_tab <- summary_tab( var_names=dimnames(tb)[[2]], event=event, data=data_rws, id_name=id,  mod=mod, lab_orders=lab_orders, delete_coef_no_events=delete_coef_no_events )
        mod <- summary(mod)
      }
  
      if(lprint)  print(format( res_tab, digits=2, justify="left" ))
    }
  }
  
  ret <- list( res_tab = res_tab, 
               model   = mod,
               call    = list( match.call()) 
  )
  if(save_data) ret <- c( ret, data_rws=list(data_rws) )
 
  
  ret
  
}  # end of function 'scri'

  
  
######################################
#######################################
######################################  
#  
#
#
refresh_event_variable <- function(  start_interval, end_interval, 
                                     event, event_time, 
                                     data
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
                        start_obs, end_obs,           # start and end of observation variables
                        event_time, event,
                        id,                           # id variable to get new start and end of observation  ('rw_min' and 'rw_max')
                        data,                         # dataset
                        lab_orders,
                        ref=1,                        # reference category for new variable: number OR [the beginning of ] category name OR "most events"
                        event_name="event",           # used if ref=="most events" to define reference category with most events
                        rw_observed_percentage=100,   # 100% - the whole interval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                        censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for variables 'censored_vars'. 
                                                      #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                        event_in_rw=T                 # if event in rw ==> this rw shouldnot be deleted even if not completely observed
){ 
  
    # create and calculate variable from formula from list 'rws' and dataset 'data': 'lab' (or another name), 'rw_start' and 'rw_end':
    calc_wind <- function(x, data){

      t0   <- unlist(lapply(as.expression(x$t0),  eval,data))
      tend <- unlist(lapply(as.expression(x$tend),eval,data))

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
          cond[is.na(cond)] <- F
      
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
              cond[is.na(cond)] <- F

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
    }   # end of sub-function 'calc_wind'
    

    
    if(substitute(rws)[[1]]=="list") 
      rws_datalist <- lapply( substitute(rws)[-1], calc_wind ,data)
    else
      rws_datalist <- lapply( eval(parse(text = rws), data), calc_wind, data )

    # create one dataset 'data_rws' from list of dataset 'rws_datalist':
    data_rws <- do.call("rbind.data.frame",rws_datalist)
    
    # check:
    if(!missing(start_obs) & !missing(end_obs)) 
      if(any(data_rws[,start_obs] >= data_rws[,end_obs])) stop("'start_obs' should be < 'end_obs'")

    # delete [parts of] intervals if interval before 'start_obs'.
    if(!missing(start_obs)){ 
      data_rws$rw_end[ data_rws$rw_end < data_rws[,start_obs] ] <- data_rws$rw_start[ data_rws$rw_end < data_rws[,start_obs] ] <- NA
      cond <-      data_rws$rw_start < data_rws[,start_obs] & data_rws[,start_obs] <= data_rws$rw_end & !is.na(data_rws$rw_start) & !is.na(data_rws$rw_end)
      data_rws$rw_start[cond] <- data_rws[cond,start_obs]
    }
    # delete [parts of] intervals if interval after 'end_obs'.
    if(!missing(end_obs)){
      data_rws$rw_start[ data_rws[,end_obs] < data_rws$rw_start ] <- data_rws$rw_end[ data_rws[,end_obs] < data_rws$rw_start ] <- NA
      cond <- data_rws[,end_obs] < data_rws$rw_end & data_rws$rw_start <= data_rws[,end_obs] & !is.na(data_rws$rw_start) & !is.na(data_rws$rw_end)
      data_rws$rw_end[cond] <- data_rws[cond,end_obs]
    }
    
    
    
    labels   <- data_rws$lab[!duplicated(data_rws$lab)]
    data_rws <- data_rws[!is.na(data_rws$rw_start),] 
 
    data_rws <- refresh_event_variable( "rw_start", "rw_end", event, event_time, data_rws)
    data_rws$lab <- factor_ref(  data_rws$lab, 
                                 lab_orders = lab_orders,
                                 lab=labels, 
                                 lab_sort=F,
                                 ref=ref,    event_var=data_rws[,event] ) 
    
    levels(data_rws$lab) <- format(levels(data_rws$lab))
    
    data_rws <- data_rws[order( data_rws[,id] ), ]

    
     # add new variables for start and end of observation ( 'rw_min' and 'rw_max' )
    rw_min <- tapply(data_rws$rw_start, data_rws[,id], min,na.rm=T)    # eval(substitute(id),data_rws),
    rw_max <- tapply(data_rws$rw_end,   data_rws[,id], max,na.rm=T)    # eval(substitute(id),data_rws),
    tmp <- cbind.data.frame(id=names(rw_min),rw_min=rw_min,rw_max=rw_max) 
    dimnames(tmp)[[2]][1] <- id
 
    data_rws <- merge.data.frame(data_rws,tmp, all=T) 
    
    data_rws
} # end of function 'create_rws'  (version 3) 
  
#  data_rws  <- create_rws( 
#    list(   
#      prevax  = list( t0 = days_vax1, cuts = c( -90, -29, 0),  lab = c( "pre-exposure",  "buffer" ), lab_add_interval=T, no_last_interval=T   ), # also can be added:    add_interval = T 
#      v1      = list( t0 = days_vax1, cuts = c( 0, 1, 29),     tend = pmin(days_vax2-1, study_exit_days, na.rm=T),  lab = "v1 "   ), 
#      v2      = list( t0 = days_vax2, cuts = c( 0, 1, 29),     lab = "v2 "   )
#      #all      = list( t0 = study_entry_days, tend=study_exit_days,     lab = "all", lab_add_interval=F   )
#    ), 
#    start_obs = study_entry_days, end_obs = study_exit_days,
#    id = person_id,
#    data = scri_input[scri_input$age_at_study_entry<=30,][1:5,],
#    ref="pre-exp"    #  ref=5 OR ref= "pre-exposure [-90;-30]"  or a part of a ref.category name
#  ) 
  
  
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

  
  
  data              <- refresh_event_variable( start_interval, end_interval, event, event_time, data)
  if(length(lab_orders)==1 & is.na(lab_orders[1]))
    data[,splits_names] <- factor_ref(  data[,splits_names], lab=lab, ref=ref, event_var=data[,event] ) 
  else
    data[,splits_names] <- factor_ref(  data[,splits_names], lab_orders = lab_orders, lab=lab, ref=ref, event_var=data[,event] )  
  
  
  data <- data[order( data$i_, data[,start_interval] ), ]   #   eval(substitute(data$start_interval))),] 
  data$i_ <- NULL
  data
}  # end of function 'split_intervals'

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
                          delete_coef_no_events = T,
                          coef_cond, model_number=1,
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
    for(ivar in var_names)
      res_tab[ivar==var_names][[1]] <- cbind.data.frame(
        n_events      = with(data[cond,],         tapply(get(event), get(ivar), sum, na.rm=T)),                   #   event
        cum_ev        = cumsum( with(data[cond,], tapply(get(event), get(ivar), sum, na.rm=T)) ),                 #   event
        atrisk_days   = with(data[cond,], tapply(interval,get(ivar), sum   , na.rm=T)),
        atrisk_ids    = with(data[cond,], tapply(get(id_name),get(ivar),function(x)length(unique(x))))
      )
    res_tab <- do.call("rbind.data.frame",res_tab)
    res_tab <- cbind.data.frame( i=1:nrow(res_tab), all_cat=rownames(res_tab) , res_tab)
    res_tab$days_pp        <- round(res_tab$atrisk_days / res_tab$atrisk_ids,2)
    res_tab$relative_rate  <- res_tab$n_events / res_tab$atrisk_days
    res_tab$relative_perc  <- res_tab$n_events / res_tab$atrisk_ids
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
    names(res_tab_model)[match(c("exp","low","upp","Pr("),substring(names(res_tab_model),1,3))] <- c("RR","2.5%","97.5%","pval")
  
    res_tab_model$model <- model_number
    
    if(!missing(coef_cond)) 
      res_tab_model       <- res_tab_model[coef_cond,]
    
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
    
    if(delete_coef_no_events) 
      res_tab[ res_tab$n_events==0 & !is.na(res_tab$n_events), "RR" ] <- NA  
  } 
  
  if(missing(data) & !missing(mod)) { 
    res_tab <- res_tab_model
    res_tab <- cbind.data.frame( i=1:nrow(res_tab), res_tab[, c(ncol(res_tab),1:(ncol(res_tab)-1))])
  }
  
  if(!missing(data) & !missing(mod) & delete_coef_no_events) 
    res_tab[ res_tab$n_events==0 & !is.na(res_tab$n_events), "RR" ] <- NA  
  
  res_tab <- cbind.data.frame( event=substring(event,1,7), res_tab )
  
  res_tab$model <- model_number
 
  if(add){
  
    res_tab_new <- res_tab0
    res_tab_new$all_cat_without_space <- gsub(" ","",res_tab_new$all_cat)
    res_tab$all_cat_without_space     <- gsub(" ","",res_tab$all_cat)
    res_tab_new <- merge.data.frame(res_tab_new,res_tab, by=c("all_cat_without_space"), all=T )
    names(res_tab_new)

    for(ivar in c( "event","all_cat","n_events","atrisk_days","atrisk_ids","days_pp","relative_rate","relative_perc")){
      names(res_tab_new)[names(res_tab_new)==paste0(ivar,".x")] <- ivar
      if( any( (cond<-is.na(res_tab_new[,ivar]) & !is.na(res_tab_new[,paste0(ivar,".y")])) ) )
        res_tab_new[cond,ivar] <- res_tab_new[cond,paste0(ivar,".y")]
      if( all( res_tab_new[,ivar] == res_tab_new[,paste0(ivar,".y")], na.rm=T ) )
        res_tab_new[,paste0(ivar,".y")] <- NULL
    }
    
    if(!missing(mod)){
      model_res_names <- c("i","RR","2.5%","97.5%","pval","coef","se(coef)","model")
      
      var_numbers <- match(paste0(model_res_names,".x"),names(res_tab_new))
      var_numbers <- var_numbers[!is.na(var_numbers)]
      names(res_tab_new)[var_numbers] <-  model_res_names[!is.na(var_numbers)]
      
      
      if( any( names(res_tab_new)=="pval") & any(names(res_tab_new)=="pval.y" ) ){
        
        # duplicated sets of columns with different values
        if(sum( cond<- !is.na( res_tab_new[,"pval"] ) & !is.na( res_tab_new[,paste0("pval.y")]) )>0){
          res_tab_new_dupl <- res_tab_new[cond, , drop=F]
          res_tab_new_dupl[, model_res_names ] <- NULL
          names(res_tab_new_dupl)[match(paste0(model_res_names,".y"),names(res_tab_new_dupl))] <-  model_res_names
        }
        
        cond <- is.na( res_tab_new[,"pval"] ) & !is.na( res_tab_new[,paste0("pval.y")]  )
        res_tab_new[cond, model_res_names ] <- res_tab_new[cond, paste0(model_res_names,".y") ]
        res_tab_new[,match(paste0(model_res_names,".y"),names(res_tab_new))] <-  NULL
        
      }
    }

    if( any(ls()=="res_tab_new_dupl") )
      res_tab_new <- rbind.data.frame(res_tab_new, res_tab_new_dupl[,match(names(res_tab_new),names(res_tab_new_dupl))]  )
    
    res_tab_new$all_cat_without_space <- NULL
    
   
    if(all(names(res_tab_new)!="model")){
      if(any(names(res_tab_new)=="i.x")) names(res_tab_new)[names(res_tab_new)=="i.x"    ] <-  "i"
      if(any(names(res_tab_new)=="i.x")) names(res_tab_new)[names(res_tab_new)=="i.model"] <-  "model"
      if(any(cond<-!is.na(match(names(res_tab_new),c("i.y","model.y")))))   res_tab_new[,cond] <-  NULL
      res_tab_new$model <- 1
    }
    
    # order table on 'all_cat'
    res_tab_new$tablab_order <- match( res_tab_new$all_cat, levels(factor_ref(res_tab_new$all_cat[ substring(res_tab_new$all_cat,1,1)!="[" ], lab_orders=lab_orders)) )
    res_tab_new$tablab_order[is.na(res_tab_new$tablab_order)] <- 100000
    res_tab_new <- res_tab_new[order(res_tab_new$tablab_order, res_tab_new$model, res_tab_new$i),]
    res_tab_new$tablab_order <- NULL
    
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

  if( !( length(lab_orders)==1 & is.na(lab_orders[1]) ) ){
    
    # lab_orders <- list( c("pre-","buf", "dose 1",     "dose 2" ),
    #                c("[0;0]","[1;28]",">28"),
    #                c("Pfizer","Moderna","AstraZeneca", ... ) )
    
      if(!is.list(lab_orders)) lab_orders <- list(lab_orders)
      
      var_unique <- unique(var)
      orders <- as.data.frame( matrix(0, nrow=length(var_unique), ncol=length(lab_orders),
                               dimnames=list(var_unique, 1:length(lab_orders)) ))
      
      for(itype in 1:length(lab_orders))
        for(icat in 1:length(lab_orders[[itype]]))
          orders[ grep( lab_orders[[itype]][icat], var_unique, fixed=T), itype] <- icat

      
      orders <- eval(parse(text=paste0("orders[ order(", paste( paste0("orders[,", 1:length(lab_orders)),collapse="],") ,"]),]")) )
      orders$order <- 1:nrow(orders)
      
      var <- factor(var, levels=rownames(orders), labels=rownames(orders) )
  }
  else {
    if(is.factor(var)){
      if(nlevels(var)>1) ref0 <- dimnames(contrasts(var))[[1]][rowSums(contrasts(var))==0]
      levels0 <- levels(var)
      var <- as.character(var)
      if( all(levels0==as.character((1:length(levels0)))) ) var <- as.numeric(var)
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
    if(is.factor(var1)) ref_cat1 <- dimnames(contrasts(var1))[[1]][rowSums(contrasts(var1))==0]
    if(is.factor(var2)) ref_cat2 <- dimnames(contrasts(var2))[[1]][rowSums(contrasts(var2))==0]
    
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
              if(any(grepl(ref[2],var2)))
                ref_cat2 <- grep(ref[2],var2, value=T) 
          }    
        }  
      }  
      else {
        if(is.numeric( ref )) ref_cat <- var_levels[ref]
        if(is.character(ref)) {
          if(ref=="most events") ref_cat <- ref
          else  
            if(any(grepl(ref,var_levels)))
              ref_cat <- grep(ref, var_levels, value=T)
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
plot_res <- function(res, main="", col=col_list, time_cat_i=length(strata), ylim, CI = T, 
                     correct_max_time_adj =  T){
 
  if(all(unlist(lapply(res, function(x)sum(!is.na(x$RR))))==0)) return(NULL) 
  
  ncoef     <- nrow(res[[1]])
  ncoef_max <- max(unlist( lapply(res,nrow) ))
  if(ncoef_max>ncoef) # +1 because for ref.category of cal_time_var 'model' is NA
    ncoef_max <- ncoef + 1 + max(unlist( lapply(res,function(x) {
      cond <- (1:nrow(x))[(ncoef+1):nrow(x)]; 
      if(all(dimnames(x)[[2]]!="model")) return(0)
      if(all(is.na(x[cond,"model"])))    return(0)
      tapply( x[cond,"i"], x[cond,"model"], length ) 
      } )))

  if(missing(ylim)){
    ymax   <- unlist(lapply(res,function(x) {
      if(all(dimnames(x)[[2]]!="RR")) return(NA) 
      x_tmp <- x$RR[x$RR<1000]; 
      if(any(!is.na(x_tmp))) max(x_tmp,na.rm=T) else NA 
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
        if(any(dimnames(res[[imod]])[[2]]=="RR"))
          matlines( rbind( (1:ncoef+x_deltas[imod]),(1:ncoef+x_deltas[imod]))[,!is.na(res[[imod]]$RR[1:ncoef])],
                    t(res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("2.5%","97.5%")]),
                    lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
      }
    
    # RR's:
    for(imod in 1:length(res)){
      if(all(dimnames(res[[imod]])[[2]]!="RR")) next
      lines( 1:ncoef+x_deltas[imod],res[[imod]]$RR[1:ncoef], type="o", col=col[imod],lwd=ifelse(imod==1,2,1)); 
      if(imod==1) text( 1:ncoef,res[[1]]$RR, labels=as.character(res[[1]]$i), pos=3, col=col[1], cex= ifelse(ncoef<=50,1,0.7) ) 
      if(any( (cond <- !is.na(res[[imod]]$pval[1:ncoef]) & res[[imod]]$pval[1:ncoef]<=0.05) ))                # check for significant p-values
        points( (1:ncoef+x_deltas[imod])[cond], res[[imod]]$RR[1:ncoef][cond], pch="*",cex=3, col=col[imod]) 
    }
    
    # calendar time adjusted
    if(length(res)>1){
      for(imod in 2:length(res)){
        
        if(all(dimnames(res[[imod]])[[2]]!="RR")) next

        cond_after_ncoef <- (1:nrow(res[[imod]]))>ncoef

        if(ncoef_max>ncoef){
          if(!correct_max_time_adj){
            # CI's for adjusted RR's:
            if(CI)
              matlines( rbind( ncoef +  1:(nrow(res[[imod]])-ncoef), ncoef +  1:(nrow(res[[imod]])-ncoef)  )[,!is.na(res[[imod]]$RR[ cond_after_ncoef ])],
                        t(res[[imod]][cond_after_ncoef,][ !is.na(res[[imod]]$RR[ cond_after_ncoef ])  ,c("2.5%","97.5%")]),
                        lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
            # RR's:
            lines( ncoef +  1:(length(res[[imod]]$RR)-ncoef), res[[imod]]$RR[ cond_after_ncoef ], type="o", col=col[imod],cex=0.5)
          }  
          else
            if(correct_max_time_adj){
              
              cond_time_adj <- substring(res[[imod]][cond_after_ncoef,"all_cat"],1,1)=="["
              
              suppressWarnings(
                mids <- unlist(   lapply( strsplit(  res[[imod]]$all_cat[cond_after_ncoef][ cond_time_adj ], "\\[|;|,|\\]") , 
                                          function(x) mean(as.numeric(x),na.rm=T) ) )
              )
              
              xx_time_adj    <-  ( (mids - min(numbers, na.rm=T) ) * ncoef_max  +( max(numbers,na.rm=T) - mids) * (ncoef+1) ) / ( max(numbers,na.rm=T) -  min(numbers, na.rm=T) )              
              #xx_time_adj    <- ncoef +  1   +   (mids - min(numbers, na.rm=T) ) * ( ncoef_max - ncoef ) / ( max(numbers,na.rm=T) -  min(numbers, na.rm=T) )              
              RR_time_adj    <- res[[imod]]$RR[   cond_after_ncoef][cond_time_adj] 
              CI_time_adj    <- res[[imod]][cond_after_ncoef,][cond_time_adj,c("2.5%","97.5%")] 
              pval_time_adj  <- res[[imod]]$pval[ cond_after_ncoef][cond_time_adj] 
              model_time_adj <- res[[imod]]$model[cond_after_ncoef][cond_time_adj] 
              
              if(any(!is.na(model_time_adj))){
                ref_models <- sort(unique( model_time_adj[!is.na(model_time_adj)] ))
                x_ref_deltas <- 0
                if(length(ref_models)>1) x_ref_deltas <- 0.2* (1/(length(ref_models)-1) * (0:(length(ref_models)-1)) - 0.5)
                
                for(ii in ref_models ){
                  cond_ref_model <- !is.na(model_time_adj) & model_time_adj==ii
                  # CI's:
                  if(CI)
                    matlines( rbind( xx_time_adj, xx_time_adj )[, cond_ref_model & !is.na(RR_time_adj) ]+x_ref_deltas[ii], 
                              t(CI_time_adj)[                   , cond_ref_model & !is.na(RR_time_adj) ],
                              lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
                  # RR's:
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
  # plot 2 en 3(if ymax>30): coefficients without coefficients for time adjusting:
  #
  for(i in 2:3) {   # plot 2    and   plot 3
    
    if(i==3 & max(ylim)<=30) next
    if(i==3 & max(ylim)>30 ) ylim=c(0,20)
      
    if(i==2 & missing(ylim)){
      ymax   <- unlist(lapply(res,function(x) {
        if(all(dimnames(x)[[2]]!="RR")) return(NA) 
        x_tmp <- x$RR[x$RR<1000]; 
        if(any(!is.na(x_tmp))) max(x_tmp,na.rm=T) else NA 
      } ))
      
      if(any(!is.na(ymax))) ymax <- max(ymax, na.rm=T)
      else                  ymax <- 10
      ylim <- c( 0, ymax )
    }  
    
    
    plot( c(0,ncoef), ylim, type="n", main=main, xlab="effect number", ylab=ifelse(i==2, "RR", "RR under 20") )
    
    grid();abline(h=1, col="darkgray",lty=1)
    abline(v=length(res[[1]]$RR[1:ncoef])+0.5, col="orange", lty=3)

    #  CI's for unadjusted and adjusted RR's:
    if(CI)
      for(imod in 1:length(res)){
        if(all(dimnames(res[[imod]])[[2]]!="RR")) next
          matlines( rbind( (1:ncoef+x_deltas[imod]),(1:ncoef+x_deltas[imod]))[,!is.na(res[[imod]]$RR[1:ncoef])],
                    t(res[[imod]][1:ncoef,][ !is.na(res[[imod]]$RR[1:ncoef])  ,c("2.5%","97.5%")]),
                    lty=1, lwd=1, col=col_alpha(col[imod],0.15), type="o", pch="-", cex=2 )
      }
    # RR's:
    for(imod in 1:length(res)){
      if(all(dimnames(res[[imod]])[[2]]!="RR")) next
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

  
#   ######
#   # plot coefficients without coefficients for time adjusting:
# 
#   rest_strata <- res[[1]]$all_cat
#   strata <- list()
#   while(!all(regexpr("&", rest_strata, fixed=T )==-1)){
#     strata      <- c( strata, list(substring( rest_strata, 1, regexpr("&", rest_strata, fixed=T ) - 1)) )
#     rest_strata <- substring( rest_strata,    regexpr("&", rest_strata, fixed=T ) + 1) 
#   }
#   strata <- c(strata, list(rest_strata))
# 
#   # get start and stop of each interval:    
#   t_start <- substring( strata[[time_cat_i]], regexpr("[", strata[[time_cat_i]],  fixed=T ) + 1 )
#   t_stop  <- substring( t_start,              regexpr(";", t_start,         fixed=T ) + 1 )
#   t_stop  <- substring( t_stop,       1,      regexpr("]", t_stop,          fixed=T ) - 1 )
#   t_start <- substring( t_start,      1,      regexpr(";", t_start,         fixed=T ) - 1 )
#   
#   if(sum(cond <- regexpr(">",strata[[time_cat_i]][t_start==""],fixed=T)!=-1)>0)
#     t_start[t_start==""] <- substring( strata[[time_cat_i]][t_start==""], regexpr(">",strata[[time_cat_i]][t_start==""],fixed=T)+1)
#   if(sum(cond <- regexpr(">=",strata[[time_cat_i]][t_start==""],fixed=T)!=-1)>0)
#     t_start[t_start==""] <- substring( strata[[time_cat_i]][t_start==""], regexpr(">",strata[[time_cat_i]][t_start==""],fixed=T)+2)
#   if(sum(cond <- regexpr("<",strata[[time_cat_i]][t_stop==""],fixed=T)!=-1)>0)
#     t_stop[t_stop==""] <- substring( strata[[time_cat_i]][t_stop==""], regexpr("<",strata[[time_cat_i]][t_stop==""],fixed=T)+1)
#   if(sum(cond <- regexpr("<=",strata[[time_cat_i]][t_stop==""],fixed=T)!=-1)>0)
#     t_stop[t_stop==""] <- substring( strata[[time_cat_i]][t_stop==""], regexpr("<=",strata[[time_cat_i]][t_stop==""],fixed=T)+2)
#   
#   t_start <- as.numeric(t_start)
#   t_stop  <- as.numeric(t_stop )
#   
#   t_start[is.na(t_stop)] <- t_start[is.na(t_stop)] + 3
#   
#   
#   ymax <- max(unlist(lapply(res,function(x)max(x$RR[1:ncoef][x$RR[1:ncoef]<500],na.rm=T))))
#   matplot( rbind(t_start, t_stop), rbind(res[[1]]$RR,res[[1]]$RR),  type="b", col="blue", ylim=c(0,  max(ymax ,2)),
#            main=main,xlab="coefficients number",ylab="RR"); 
#   
#   
#   ymax <- max(unlist(lapply(res,function(x)max(x$RR[1:ncoef][x$RR[1:ncoef]<500],na.rm=T))))
#   plot( 0.5*(t_start+ t_stop), res[[1]]$RR,  type="o", col="blue", ylim=c(0,  max(ymax ,2)),
#            main=main,xlab="coefficients number",ylab="RR"); 
#   
#   
#   grid();abline(h=1, col="darkgray",lty=1)
#   abline(v=length(res[[1]]$RR[1:ncoef])+0.5, col="orange", lty=3)
#   lines( 1:length(res[[1]]$RR[1:ncoef])+0.15,res[[1]]$RR[1:ncoef], type="b", col=col[1],lwd=2); 
#   if(length(res)>1){
#     for(i in 2:length(res))
#       lines( res[[i]]$RR[1:ncoef], type="b", col=col[i]); 
#   }  
#   legend("topright",legend= paste(res[[1]]$i,res[[1]]$all_cat), cex=0.7, box.lty=0,bg="transparent")
#   legend("topleft",legend= names(res), cex=0.7, pch=1, lwd=1, col=col[1:length(res)], box.lty=0,bg="transparent")
  
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
#  functions to add results to report list and model list and save them.
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
#  save report_list and model_list
#
save_results <- function(name, report_list, models_list, sep="" ){
  
  name <- gsub(":","",name)

  begin_report <- paste0(dap, sep, name, "_report")
  
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
scri_strata <- function(strata_var = "all data", 
                        output_name, 
                        formula_text, time_seq=c(), time_seq_ref="most events",
                        event_time, event, id,
                        rws ,
                        time_dep = NA,
                        combine_vars = c(),
                        start_obs, end_obs,
                        lab_orders = NA,
                        data, cond_var = "",
                        rw_observed_percentage=100,   # 100% - the whole inteval must be observed; 0% - even one day is enough to include this 'id' in the current risk window
                        censored_vars=c(),            # The rule 'rw_observed_percentage' does not work for varuables 'censored_vars'. 
                                                      #  (for example, "death_days" ==> 'id' is included in the corresponding risk window till death date.)
                        event_in_rw=T,                # if event in rw ==> this rw shouldnot be deleted even if not completely observed
                        image_plots = T, image_file_separate=F, image_strata, image_tit="", image_brand=F, 
                        nvax = 2,
                        delete_coef_no_events = T,
                        paral = T,  paral_vars = c(), n_cores=NA,               # parallel programming (using more cores)
                        lprint=T, 
                        global_plot_name, add_global_plot = T, CI = T,
                        extra_plots=F, width=14, ...
){   

  if( length(time_seq)>1 & length(time_seq_ref)==1 ) time_seq_ref <- rep(time_seq_ref,length(time_seq))
  
  if(!missing(strata_var)){
      if(is.factor(data[,strata_var[1]])) strata <- levels(       data[,strata_var[1] ])
      else                                strata <- unique(unlist(data[,strata_var    ]))  
  }
  else strata <- strata_var

  res_strata <- vector("list", length=length(strata))
  names(res_strata) <- strata 
  tabs_strata <- res_strata
  
  for(istr in strata){
    if(any(ls()=="res")) rm(res)
    
    if(is.na(istr)){ 
      cond_stratum <- is.na(data[,strata_var])
      if(length(strata_var)>1)
        cond_stratum <- rowSums(cond_stratum)>0
    }  
    else { 
      if(istr %in% c("all data","all_data","alldata")) cond_stratum <- rep(T,nrow(data))
      else{
        cond_stratum <- data[,strata_var] == istr   & !is.na(data[,strata_var])
        if(length(strata_var)>1)
          cond_stratum <- rowSums(cond_stratum, na.rm=T)>0
      }
    }  
  
    # create image plots:
    if(image_plots){  
      
      pdf(file=paste0(sdr,"image_",gsub(" |:","",istr),"_tmp.pdf"), width=width,  ...)
        
        if(!missing(image_strata)){
          for(im_str in unique(data[ cond_stratum, image_strata ]) ){
            if(image_brand){
              for(im_br in unique(data[ cond_stratum & data[,image_strata ]==im_str, paste0("type_vax",1) ]) )
                brand_images( data[ cond_stratum & data[,image_strata ]==im_str, ], ae_event=event, brand=im_br, tit=paste0(image_tit,im_str) )
            }  
            else 
              brand_images( data[ cond_stratum & data[,image_strata ]==im_str, ], ae_event=event,   tit=paste0(image_tit,im_str) )
          }  
        }
        else {
          if(image_brand){
            for(im_br in unique(c(data[ cond_stratum, paste0("type_vax",1) ])) )
              brand_images( data[cond_stratum,], ae_event=event, brand=im_br, tit=image_tit )
          }  
          else brand_images( data[cond_stratum,], ae_event=event, tit=image_tit  )
        }  
      
      dev.off()
      
      if(image_file_separate) {

        if( add_global_plot &  file.exists(paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),"_images.pdf")) ){
          
          file.copy( from=paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),"_images.pdf"), to=paste0(sdr,"_temp_all_images.pdf"), overwrite=T)   
          files_to_copy <- c( paste0(sdr,"_temp_all_images.pdf"), paste0(sdr,"image_",gsub(" |:","",istr),"_tmp.pdf"))
          qpdf::pdf_combine( files_to_copy , paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),"_images.pdf")  )  
          if(file.exists(paste0(sdr,"_temp_all_images.pdf")))     suppressWarnings( file.remove(paste0(sdr,"_temp_all_images.pdf")) )
          
        } 
        else file.copy( from=paste0(sdr,"image_",gsub(" |:","",istr),"_tmp.pdf"), to=paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),"_images.pdf"), overwrite=T)
      }
    }
   
    res <- vector("list", length=length(time_seq)+1)
    names(res) <- c("no_adj", names(time_seq) )
    tabs <- res
    
    res[[1]] <-  scri_sv(formula = formula(formula_text),
                         event_time = event_time, event = event, id=id,
                         rws        = rws, 
                         time_dep   = time_dep,
                         combine_vars = combine_vars,
                         start_obs  = start_obs, end_obs = end_obs,
                         lab_orders = lab_orders,
                         data       = data[cond_stratum,], cond_var = cond_var,
                         nvax       = nvax,
                         delete_coef_no_events = delete_coef_no_events,
                         rw_observed_percentage = rw_observed_percentage,   
                         censored_vars = censored_vars,           
                         event_in_rw = event_in_rw,               
                         lprint=lprint )#[[1]]

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
      clusterExport(cl, c("scri_sv", "refresh_event_variable", "create_rws", "split_intervals", 
                          "summary_tab", "factor_ref", "combine_vars_func" ), 
                    envir=environment()  )
      for(ipar in names(as.list(match.call()))[-1])
        clusterExport(cl, ipar, envir=environment() )        
      if(!missing(paral_vars)) clusterExport(cl, paral_vars )
      
      res_paral <- parLapply(cl,
                             time_seq,  # list with different sets of intervals
                             function(i_time_seq) 
                               scri_sv(formula = formula(paste( formula_text, " +  cal_time_cat")),
                                       event_time = event_time, event = event, id=id,
                                       rws = rws, 
                                       time_seq = i_time_seq,  split_seq_name = "cal_time_cat", #time_seq_ref="most events",time_seq_ref=time_seq_ref[i_time_seq],
                                       time_dep   = time_dep, 
                                       combine_vars = combine_vars,
                                       start_obs = start_obs, end_obs = end_obs,
                                       lab_orders = lab_orders,
                                       data = data[cond_stratum,], cond_var = cond_var,
                                       nvax       = nvax,
                                       delete_coef_no_events = delete_coef_no_events,
                                       rw_observed_percentage = rw_observed_percentage,   
                                       censored_vars = censored_vars,           
                                       event_in_rw = event_in_rw,               
                                       lprint=lprint    )
      )
      
      stopCluster(cl)
     
      names_overlap <- match(names(res),names(res_paral))
      res[ names(res_paral)[ names_overlap[!is.na(names_overlap)] ] ] <- res_paral[!is.na(match(names(res_paral),names(res)))]

    }  # end of parallel
    
    if(length(time_seq)>0 & !paral)
      for(i in 2:length(res)){
        res[[i]] <-  scri_sv(formula = formula(paste( formula_text, " +  cal_time_cat")),
                             event_time = event_time, event = event, id=id,
                             rws = rws, 
                             time_seq = time_seq[[i-1]] , time_seq_ref=time_seq_ref[i-1], split_seq_name = "cal_time_cat", #time_seq_ref="most events", 
                             time_dep   = time_dep,
                             combine_vars = combine_vars,
                             start_obs = start_obs, end_obs = end_obs,
                             lab_orders = lab_orders,
                             data = data[cond_stratum,], cond_var = cond_var,
                             nvax       = nvax,
                             delete_coef_no_events = delete_coef_no_events,
                             rw_observed_percentage = rw_observed_percentage,   
                             censored_vars = censored_vars,           
                             event_in_rw = event_in_rw,               
                             lprint=lprint    )#[[1]]

      }  
    
    for(i in 1:length(res)) tabs[[i]] <- res[[i]][[1]]  
    
    res_strata[[istr]]  <- res
    tabs_strata[[istr]] <- tabs
  } # end of loop 'for' with 'istr'

  # create plots with coefficients:
  if(extra_plots) plot_name <- paste0(sdr,"_",gsub(" |:","", output_name),".pdf")
  else            plot_name <- paste0(sdr,"_temp_plot.pdf")
 
  pdf(file=plot_name, width=width,  ...)
    for(istr in names(tabs_strata) )
      if(!is.null(tabs_strata[[istr]][[1]]))
        plot_res(tabs_strata[[istr]], main=paste( event, formula_text,"; ",image_tit,"\n",output_name), col=col_list, CI=CI)
  dev.off()
  
  # add plot[s] into a global file
  if( add_global_plot &  file.exists(paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),".pdf")) ){
    file.copy( from=paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),".pdf"), to=paste0(sdr,"_temp_plot_2.pdf"), overwrite=T)   
    files_to_copy <- c( paste0(sdr,"_temp_plot_2.pdf"), plot_name)
  } 
  else files_to_copy <- plot_name

  if(image_plots & !image_file_separate) files_to_copy <- c( paste0(sdr,"image_",gsub(" |:","", strata),"_tmp.pdf"), files_to_copy )
  
  qpdf::pdf_combine( files_to_copy , paste0(sdr,dap, "_", gsub(" |:","",global_plot_name),".pdf")  )    

  if(!extra_plots & file.exists(plot_name))                                suppressWarnings( file.remove( plot_name                                           ) )
  if(file.exists(paste0(sdr,"_temp_plot_2.pdf")))                          suppressWarnings( file.remove(paste0(sdr,                       "_temp_plot_2.pdf")) )
  if(file.exists(paste0(sdr,"image_",gsub(" |:","",strata),"_tmp.pdf")))   suppressWarnings( file.remove(paste0(sdr,"image_",gsub(" |:","",strata),"_tmp.pdf")) )
     
  c( tabs      = list(tabs_strata),
     scri_all  = list(res_strata ))
}  # the end of function 'scri_strata'


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
#######
#  data preparation for plots:
data_plot_prepair <- function(plot_data, event, vax2, vax3, vax4){
  
  # delete rows with missing events or missing vax1:
  plot_data <- plot_data[!is.na(plot_data[,paste0(event,"_date")]) & !is.na(plot_data$date_vax1), ]
  
  if(missing(vax2)){  # TRUE if variable 'date_vax2' exists and is not empty.
    vax2 <- any(names(plot_data)=="date_vax2")
    if(vax2) vax2 <- any(!is.na(plot_data$date_vax2))
  }  
  if(missing(vax3)){  # TRUE if variable 'date_vax3' exists and is not empty.
    vax3 <- any(names(plot_data)=="date_vax3")
    if(vax3) vax3 <- any(!is.na(plot_data$date_vax3))
  }  
  if(missing(vax4)){  # TRUE if variable 'date_vax3' exists and is not empty.
    vax4 <- any(names(plot_data)=="date_vax4")
    if(vax4) vax4 <- any(!is.na(plot_data$date_vax4))
  }  
  
  
  plot_data$event_date  <-  as.Date(plot_data[,paste0(event,"_date")]) 
  plot_data$date_vax1   <-  as.Date(plot_data$date_vax1) 
  plot_data$date_vax2   <-  as.Date(plot_data$date_vax2) 
  if(vax3) 
    plot_data$date_vax3 <-  as.Date(plot_data$date_vax3) 
  if(vax4) 
    plot_data$date_vax4 <-  as.Date(plot_data$date_vax4) 
  
  #??? mb to use ???
  plot_data$event_weight <- rep(0, nrow(plot_data))
  for(i in 1:length(plot_data$event_date))
    plot_data$event_weight[i] <-  sum(  plot_data$event_date[i] >= as.Date(plot_data[,"study_entry_date"])  & 
                                          plot_data$event_date[i] <= as.Date(plot_data[,"study_exit_date" ])     )   
  
  
  plot_data$event_num   <- as.numeric(difftime( plot_data$event_date, as.Date("2020-09-01"), units="days"))
  plot_data$vax1_num   <- as.numeric(difftime( plot_data$date_vax1,  as.Date("2020-09-01"), units="days"))
  plot_data$vax2_num   <- as.numeric(difftime( plot_data$date_vax2,  as.Date("2020-09-01"), units="days"))
  if(vax3) 
    plot_data$vax3_num <- as.numeric(difftime( plot_data$date_vax3,  as.Date("2020-09-01"), units="days"))
  if(vax4) 
    plot_data$vax4_num <- as.numeric(difftime( plot_data$date_vax4,  as.Date("2020-09-01"), units="days"))
  
  plot_data$event_vax1_num    <- as.numeric(difftime( plot_data$event_date, plot_data$date_vax1,     units="days"))
  plot_data$event_vax2_num    <- as.numeric(difftime( plot_data$event_date, plot_data$date_vax2,     units="days"))
  if(vax3) 
    plot_data$event_vax3_num  <- as.numeric(difftime( plot_data$event_date, plot_data$date_vax3,     units="days"))
  if(vax4) 
    plot_data$event_vax4_num  <- as.numeric(difftime( plot_data$event_date, plot_data$date_vax4,     units="days"))
  
  
  #####
  # create list with conditions for points before vax1, between vax1 and vax2, between vax2 and vax3, after vax3:
  plot_data$cond_before_vax1 <- plot_data$event_vax1_num<0
  plot_data$cond_after_vax1  <- plot_data$event_vax1_num>=0
  
  plot_data$cond_between_vax1_vax2       <- plot_data$event_vax1_num>=0 
  cond <- !is.na(plot_data$event_vax2_num)
  plot_data$cond_between_vax1_vax2[cond] <-  plot_data$event_vax2_num[cond]<0 & plot_data$cond_between_vax1_vax2[cond]
  
  plot_data$cond_after_vax2        = plot_data$event_vax2_num>=0 
  plot_data$cond_after_vax2[is.na(plot_data$cond_after_vax2)] = F 
  
  
  plot_data$cond_between_vax2_vax3       <-  plot_data$cond_after_vax2 
  if(vax3){
    cond <- !is.na(plot_data$event_vax3_num)
    plot_data$cond_between_vax2_vax3[cond] <- plot_data$event_vax3_num[cond]<0 & plot_data$cond_between_vax2_vax3[cond]
    
    plot_data$cond_after_vax3        <- plot_data$event_vax3_num>=0 
    plot_data$cond_after_vax3[is.na(plot_data$cond_after_vax3)] <- F
    
    plot_data$cond_between_vax3_vax4       <-  plot_data$cond_after_vax3 
    if(vax4){
      cond <- !is.na(plot_data$event_vax4_num)
      plot_data$cond_between_vax3_vax4[cond] <- plot_data$event_vax3_num[cond]<0 & plot_data$cond_between_vax3_vax4[cond]
      
      plot_data$cond_after_vax4        <- plot_data$event_vax4_num>=0 
      plot_data$cond_after_vax4[is.na(plot_data$cond_after_vax4)] <- F
    }  
  }
  
  plot_data
}   # end of function 'data_plot_prepair'


#######
# plot with 3 images per vaccine :
images_plots_doses <- function(data, event, tit="", 
                               plots=1:100, 
                               mfrow, mfcol, add_plots=F,
                               xmin, xmax, ymin, ymax, 
                               xlab_extra = "",
                               cex=0.8,
                               place = "bottomright"  # for legend
){
  
  data <- data_plot_prepair(data, event)
  if(nrow(data)==0)  return()
  
  if(!add_plots){
    if(missing(mfrow) & missing(mfrow)){
      if(any(names(data)=="cond_after_vax3"))
        par(mfrow=c(3,3))
      else
        par(mfrow=c(2,3))
    }  
    else{
      if(!missing(mfrow)) par(mfrow=mfrow)
      else par(mfcol=mfcol)
    }
  }  
  
  points_legend <- function(xvar, yvar, data, points=c(0,1,2,3), place, pre_lines=T, cex=cex, cex.legend=0.6, pt.cex.legend=1){
    
    ################
    # ablines:
    abline(v=c(0,7,14,21,28, 60)); abline(v=0,lwd=2)
    if(pre_lines) abline(v=c(-30,-90))
    abline(h=c(2021 + 0.25*((-1):10)), col="pink")
    abline(h=c(2021,2022), col="pink",lwd=2)
    
    ##############
    # points:
    if(any(points==0))
      points(data[ data$cond_before_vax1 & !is.na(data[,xvar]),xvar], 
             1970 + as.numeric(as.Date("2020-09-01"))/365.25 + data[ data$cond_before_vax1     & !is.na(data[,xvar]), yvar]/365.25, 
             cex=cex )
    if(any(points==1) & any(names(data)=="cond_between_vax1_vax2"))
      points(data[ data$cond_between_vax1_vax2 & !is.na(data[,xvar]), xvar], 
             1970 + as.numeric(as.Date("2020-09-01"))/365.25 + data[ data$cond_between_vax1_vax2 & !is.na(data[,xvar]), yvar]/365.25,
             col="red", cex=cex)
    if(any(points==2) & any(names(data)=="cond_between_vax2_vax3"))
      points(data[ data$cond_between_vax2_vax3 & !is.na(data[,xvar]),xvar], 
             1970 + as.numeric(as.Date("2020-09-01"))/365.25 + data[ data$cond_between_vax2_vax3 & !is.na(data[,xvar]), yvar]/365.25,
             col="magenta", cex=cex)
    if(any(points==3) & any(names(data)=="cond_after_vax3"))
      points(data[ data$cond_after_vax3 & !is.na(data[,xvar]),xvar], 
             1970 + as.numeric(as.Date("2020-09-01"))/365.25 + data[ data$cond_after_vax3        & !is.na(data[,xvar]), yvar]/365.25,
             col="purple3", cex=cex)
    
    ##############
    # legend:
    if(!missing(place)){
      legend <- c( "vertical lines: -90,-30,0,7,14,21,28,60 days" , "horiz.lines: every 3 months", "",
                   paste(event,c("before dose 1", "between dose 1 and dose 2","between dose 2 and dose 3", "after dose 3")) )
      pch    <- c(32,32,32, 1,1,1,1)
      lty    <- c(1,1,0,    0,0,0,0)
      col    <- c("black","pink","black",   "black","red","magenta","blue")
      if(F)
        if( all(names(data)!="cond_after_vax3") ) {
          legend <- legend[ -( length(legend) - c(0:1) ) ]
          pch    <- pch[    -( length(pch)    - c(0:1) ) ]
          lty    <- lty[    -( length(lty)    - c(0:1) ) ]
          col    <- col[    -( length(col)    - c(0:1) ) ]
        }  
      legend(place,legend=legend, pch=pch, lty=lty, col=col, bg="transparent", box.lty=0, cex=cex.legend, pt.cex=pt.cex.legend)
    }
    
  } # end of sub-function 'points_legend' 
  
  at_date     <- seq( 2015,  as.numeric(substring(Sys.Date(),1,4))+1, by=0.25 ) 
  at_date_lab <- month.abb[c(1,4,7,10)][1+4*at_date%%1] 
  at_date_lab[at_date_lab==month.abb[1]] <- paste0( at_date_lab, "\n", at_date%/%1 )[at_date_lab==month.abb[1]]
  
  
  #####
  #  vax 1
  ylim <- range(data$event_num, na.rm=T)
  xlim <- range(data$event_vax1_num, na.rm=T); xlim[2] <- max(28,xlim[2])
  
  if(!missing(xmin)) { xlim[1] <- xmin; bound_x <- T }
  if(!missing(xmax)) { xlim[2] <- xmax; bound_x <- T }
  if(!missing(ymin)) ylim[1] <- ymin
  if(!missing(ymax)) ylim[2] <- ymax
  
  # 1    relative to dose 1
  if(any(plots==1)){
    distr_2d <- kde2d(x=data$event_vax1_num, y=data$event_num,  h=20, n=100, lims=c(xlim, ylim))
    distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
    
    try({
      image(distr_2d, ylab=paste("date of", event), xlab=paste(event,"days from dose 1", xlab_extra), axes=F )
      axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
      points_legend("event_vax1_num", "event_num", data, place=place, cex=cex )
      title(paste0(event,"; ",tit,"\nrelative to dose 1"))
    }, silent=T)
  }
  
  # 2    relative to dose 1
  if(any(plots==2)){
    distr_2d <- kde2d(x=data$event_vax1_num, y=data$vax1_num, h=20, n=100, lims=c(xlim, ylim))
    distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
    
    
    try({
      image(distr_2d, ylab="date of dose 1",  xlab=paste(event,"days from dose 1", xlab_extra), axes=F )
      axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
      points_legend( "event_vax1_num", "vax1_num", data,  cex=cex)  # place="topright",
      title(paste0(event,"; ",tit,"\nrelative to dose 1"))
    }, silent=T)
  }
  
  # 3   relative to dose 1
  if(any(plots==3)){
    cond <- data$cond_before_vax1 | data$cond_between_vax1_vax2
    if(any(cond)){
      distr_2d <- kde2d(x=data$event_vax1_num[cond], y=data$vax1_num[cond], h=20, n=100, lims=c(xlim, ylim))
      distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
      
      
      try({
        image(distr_2d, ylab="date of dose 1",  xlab=paste(event,"days from dose 1", xlab_extra), axes=F  )
        axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
        points_legend( "event_vax1_num", "vax1_num", data, points=c(0,1), cex=cex)
        title(paste0(event,"; ",tit,"\nrelative to dose 1, before dose 2"))
      }, silent=T)
    }
  }
  
  
  #####
  #  vax 2:
  
  if(any(names(data)=="cond_after_vax2")){
    
    xlim <- range(data$event_vax2_num, na.rm=T); xlim[2] <- max(28,xlim[2])
    
    if(!missing(xmin)) { xlim[1] <- xmin; bound_x <- T }
    if(!missing(xmax)) { xlim[2] <- xmax; bound_x <- T }
    
    # 4   relative to dose 2
    if(any(plots==4)){
      if(any(!is.na(data$event_vax2_num))){
        distr_2d <- kde2d(x=data$event_vax2_num[!is.na(data$event_vax2_num)], y=data$event_num[!is.na(data$event_vax2_num)],  h=20, n=100, lims=c(xlim, ylim))
        distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
        
        
        try({
          image(distr_2d, ylab=paste("date of", event),  xlab=paste(event,"days from dose 2", xlab_extra), axes=F  )
          axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
          points_legend( "event_vax2_num", "event_num", data, pre_lines=F, cex=cex )
          title(paste0(event,"; ",tit," \nrelative to dose 2; only persons with dose 2!"))
        }, silent=T)
      }
    }
    
    # 5   relative to dose 2
    if(any(plots==5)){
      if(any(!is.na(data$event_vax2_num))){
        distr_2d <- kde2d(x=data$event_vax2_num[!is.na(data$event_vax2_num)], y=data$vax2_num[!is.na(data$event_vax2_num)], h=20, n=100, lims=c(xlim, ylim))
        distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
        
        
        try({
          image(distr_2d, ylab="date of dose 2",  xlab=paste(event,"days from dose 2", xlab_extra), axes=F  )
          axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
          points_legend( "event_vax2_num", "vax2_num", data, pre_lines=F, cex=cex)
          title(paste0(event,"; ",tit,"\nrelative to dose 2;  only persons with dose 2!"))
        }, silent=T)
      }
    }
    
    # 6   relative to dose 2
    if(any(plots==6)){
      cond <- ( !data$cond_after_vax2 | data$cond_between_vax2_vax3 )  & !is.na(data$event_vax2_num)
      if(any(cond)){
        distr_2d <- kde2d(x=data$event_vax2_num[cond], y=data$vax2_num[cond], h=20, n=100, lims=c(xlim, ylim))
        distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
        
        
        try({
          image(distr_2d, ylab="date of dose 2",  xlab=paste(event,"days from dose 2", xlab_extra), axes=F  )
          axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
          points_legend( "event_vax2_num", "vax2_num", data, points=c(0,1,2), pre_lines=F, cex=cex )
          title(paste0(event,"; ",tit,"\nrelative to dose 2, before dose 3; only persons with dose 2!"))
        }, silent=T)
      }
    }
  }
  
  #####
  #  vax 3:
  if(any(names(data)=="cond_after_vax3")){
    xlim <- range(data$event_vax3_num, na.rm=T); xlim[2] <- max(28,xlim[2])
    
    if(!missing(xmin)) { xlim[1] <- xmin; bound_x <- T }
    if(!missing(xmax)) { xlim[2] <- xmax; bound_x <- T }
    
    # 7   relative to dose 3
    if(any(plots==7)){
      if(any(!is.na(data$event_vax3_num))){
        distr_2d <- kde2d(x=data$event_vax3_num[!is.na(data$event_vax3_num)], y=data$event_num[!is.na(data$event_vax3_num)],  h=20, n=100, lims=c(xlim, ylim))
        distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
        
        
        try({
          image(distr_2d, ylab=paste("date of", event),  xlab=paste(event,"days from dose 3", xlab_extra), axes=F  )
          axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
          points_legend( "event_vax3_num", "event_num", data, pre_lines=F, cex=cex )
          title(paste0(event,"; ",tit,"\nrelative to dose 3; only persons with dose 3!"))
        }, silent=T)
      }
    }
    
    # 8   relative to dose 3
    if(any(plots==8)){
      if(any(!is.na(data$event_vax3_num))){
        distr_2d <- kde2d(x=data$event_vax3_num[!is.na(data$event_vax3_num)], 
                          y=data$vax3_num[!is.na(data$event_vax3_num)], h=20, n=100, lims=c(xlim, ylim))
        distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
        
        
        try({
          image(distr_2d, ylab="date of dose 3",  xlab=paste(event,"days from dose 3", xlab_extra), axes=F  )
          axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
          points_legend("event_vax3_num", "vax3_num", data, pre_lines=F, cex=cex)
          title(paste0(event,"; ",tit,"\nrelative to dose 3; only persons with dose 3!"))
        }, silent=T)
      }
    }  
    
    # 9   relative to dose 3
    if(any(plots==9)){
      if(any(names(data)=="cond_between_vax3_vax4")){
        cond <- ( !data$cond_after_vax3 | data$cond_between_vax3_vax4 ) & !is.na(data$event_vax3_num)
        if(any(cond)){
          distr_2d <- kde2d(x=data$event_vax3_num[cond], y=data$vax3_num[cond], h=20, n=100, lims=c(xlim, ylim))
          distr_2d$y <- 1970 + as.numeric(as.Date("2020-09-01"))/365.25 + distr_2d$y/365.25
          
          
          try({
            image(distr_2d, ylab="date of dose 3", xlab=paste(event,"days from dose 3", xlab_extra), axes=F  )
            axis(1); axis(2, at=at_date, labels=at_date_lab, las=1); box(col="yellow3")
            points_legend("event_vax3_num", "vax3_num", data, points=c(0,1,2,3), pre_lines=F, cex=cex )
            title(paste0(event,"; ",tit,"\nrelative to dose 3, before dose 4; only persons with dose 3!"))
          }, silent=T)
        }
        # else {
        #   plot(1,1,type="n", axes=F, xlab="", ylab="" );  text(1,1,"no data")
        #   title(paste0(event,"; ",tit,"relative to dose 2, before dose 3 \n only persons with dose 2!"))
        # }
      }
    } 
  }  
}  #  end of function for max 6 images 'images_plots_doses'





brand_images <- function(plot_data, ae_event, brand="", tit=""){  

  cond <- plot_data[,ae_event]==1
  if(brand!=""){
    cond <- cond & ( plot_data$type_vax1==brand | ( !is.na(plot_data$type_vax2) & plot_data$type_vax2==brand ) )
    if(any(names(plot_data)=="type_vax3"))
      cond <- cond | ( plot_data$type_vax3==brand & !is.na(plot_data$type_vax3) & plot_data[,ae_event]==1  ) 
  }
  plot_data <- plot_data[ cond, ]
  
  ylim <- range( as.numeric(difftime( as.Date(plot_data[,paste0(ae_event,"_date")]), as.Date("2020-09-01"), units="days")), na.rm=T)
  ylim[1] <- min(0, ylim[1])
  ylim[2] <- max(as.numeric(difftime( as.Date("2022-01-01"), as.Date("2020-09-01"), units="days")), ylim[2])
  
  if(brand !="") brand_tit <- paste0(" is ",brand)
  else           brand_tit <- ""
  
  if(tit!="") tit <- paste0(tit,"; ")
  
  # dose 1 is brand 'brand': 
  cond <- cond_event <- plot_data[,ae_event]==1
  if(brand!="") cond <- cond_event & plot_data$type_vax1==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 1",brand_tit),  plots=c(1:3), mfrow=c(3,3),cex=1, ymin=ylim[1], ymax= ylim[2] )
  # dose 2 is brand 'brand': 
  if(brand!="" & any(names(plot_data)=="type_vax2")) 
    cond <- cond_event &  !is.na(plot_data$type_vax2) & plot_data$type_vax2==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 2",brand_tit), plots=c(4:6), add_plots = T, cex=1, ymin=ylim[1], ymax= ylim[2] )
  # dose 3 is brand 'brand': 
  if(brand!="" & any(names(plot_data)=="type_vax3")) 
    cond <- cond_event &  !is.na(plot_data$type_vax3) & plot_data$type_vax3==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 3",brand_tit), plots=c(7:9), add_plots = T, cex=1, ymin=ylim[1], ymax= ylim[2] )
  
  
  ###############    from -120 till 60 days: #################
  # dose 1 is brand 'brand': 
  xlab_extra <- "(between -120 and 60 days)"
  cond <- cond_event <- plot_data[,paste0(ae_event,"_date")]>=plot_data$date_vax1-90 & plot_data[,ae_event]==1
  if(brand!="") cond <- cond_event & plot_data$type_vax1==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 1",brand_tit), plots=c(1:3), mfrow=c(3,3),cex=1, ymin=ylim[1], ymax= ylim[2], xmin=-120, xmax=60, xlab_extra=xlab_extra, place="topleft" )
  # dose 2 is brand 'brand': 
  if(brand!="" & any(names(plot_data)=="type_vax2")) 
    cond <- cond_event &  !is.na(plot_data$type_vax2) & plot_data$type_vax2==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 2",brand_tit), plots=c(4:6), add_plots = T, cex=1, ymin=ylim[1], ymax= ylim[2], xmin=-120, xmax=60, xlab_extra=xlab_extra )
  # dose 3 is brand 'brand': 
  if(brand!="" & any(names(plot_data)=="type_vax3")) 
    cond <- cond_event &  !is.na(plot_data$type_vax3) & plot_data$type_vax3==brand
  if(sum(cond)>=3) images_plots_doses( plot_data[cond,], ae_event, tit=paste(tit,"dose 3",brand_tit), plots=c(7:9), add_plots = T, cex=1, ymin=ylim[1], ymax= ylim[2], xmin=-120, xmax=60, xlab_extra=xlab_extra )
  
}  # the end of function 'brand_images'

# brand_images(plot_data, ae_event=ae_event, brand="Moderna", tit="")
# brand_images(plot_data, ae_event=ae_event,  tit="")



#
#   the end of functions :
#  
###############################

###############################
###############################




























