

###################################################
#
#  baseline tables
# 

characteristics <- function(data, event, path_file_name, condition_value="", 
                            vax_name="vax_number", id="pat_n", start_obs="study_entry_days",
                            vax_time="vax_days", vax_date="vax_date", 
                            death_date="date_of_death", death_days="death_days", 
                            age="age_at_study_entry", age30_50="age30_50", age30="age30",
                            sex="sex", sexc="sexc", sex_age30="sex_age30", 
                            lab_orders){
  
  cat("flowchart")
  sys_time <- Sys.time()
  
  if(!missing(path_file_name)) {
    if( substring(path_file_name,nchar(path_file_name)-3) != ".txt") path_file_name <- paste0(path_file_name,".txt")
    sink(paste0(path_file_name))
  }
  
  cat( paste("\n",paste0(paste0(c('event\t', 'vax_name', 'id\t\t', 'start_obs', 'vax_time', 'vax_date', 'death_date', 'death_days', 'age\t\t', 'age30_50', 'age30\t', 'sex\t\t', 'sexc\t', 'sex_age30' ),'\t=\t"',
                                c( event ,  vax_name,   id,     start_obs ,  vax_time ,  vax_date ,  death_date ,  death_days ,  age   ,  age30_50 ,  age30,     sex,     sexc,     sex_age30  ),'"') ,collapse = '\n')))
  if(condition_value!="") cat(paste0("\n\n\n\tConditoin name: ", condition_value, ":\n\n\n"))
  
  if(       id!="id"       ) names(data)[names(data)==id       ] <- "id"
  if(start_obs!="start_obs") names(data)[names(data)==start_obs] <- "start_obs"
  
  if(vax_name!="vax_name") names(data)[names(data)==vax_name] <- "vax_name"
  if(vax_time!="vax_time") names(data)[names(data)==vax_time] <- "vax_time"
  if(vax_date!="vax_date") names(data)[names(data)==vax_date] <- "vax_date"
  
  names(data)[names(data)==paste0(event,"_days")] <- "event_days"
  names(data)[names(data)==paste0(event,"_date")] <- "event_date"
  
  data$event  <- as.numeric(!is.na(data$event_days))
  data$eventc <- paste0("event:",data$event)
  
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
  names_to_delete <- c(id,vax_name, vax_time, vax_date, event, paste0(event,"_days"), paste0(event,"_date"), death_days, death_date, age, sex, age30_50, age30, sexc, sex_age30 )
  names_to_delete <- names_to_delete[ !(names_to_delete %in% c("id","vax_name", "vax_time", "vax_date",
                                                               "event", "event_days", "event_date", "eventc", 
                                                               "death_days","death_date",
                                                               "age", "sex", "age30_50", "age30", "sexc", "sex_age30")) ]
  data[,names_to_delete] <- NULL
  
  flowchart_all <- vector("list",length=7)
  names(flowchart_all) <- c("all_data",
                            "with_vax_event",
                            "observed_before_vax1-90_and_events_after_vax1-90",
                            "observed_before_vax1_and_events_after_vax1-90",
                            "observed_before_vax1_and_events_after_vax1",
                            "observed_before_vax1-90_and_events_during_28_days_after_vax",
                            "observed_before_vax1_and_events_during_28_days_after_vax")
  
  data0 <- data[ !duplicated(data[,c("id","vax_time")]),c("id","event", "eventc","event_days","event_date",  "start_obs",
                                                          "vax_name","vax_brand", "vax_time", "vax_date", "vax_n", "vax_days_v1", 
                                                          "death_days","death_date","age","age30_50", "age30", "sexc", "sex_age30")]
  gc()
  
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
      data        <- data0[data0$event==1  & !is.na(data0$vax_date),]
      data_deaths <- data0[!is.na(data0$death_days) ,]
      cat("\t\t\tvaccinated persons with ",event,"\n\n")
    }
    
    ############### 
    # from here only with event at -91 days before vax1:
    if(i==2){
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1-90 <= data0$event_days, ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1-90 <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at (vax_day1 - 90 days) with ",event," after (vax_day1 - 90 days)\n\n")
    }
    ############### 
    # from here only with event at -91 days before vax1: vaccinated persons observed at vax_day1 with ",event," after (vax_day1 - 90 days)
    if(i==3){
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1-90 <= data0$event_days, ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1-90 <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at vax_day1 with ",event," after (vax_day1 - 90 days)\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," after (vax_day1 - 0 days)
    if(i==4) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$event_days, ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$death_days, ]
      cat("\t\t\tvaccinated persons observed at the first dose with ",event," after (vax_day1 - 0 days)\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," in [vax_date; vax_date + 28 days]
    if(i==5) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$event_days &
                             data0$event_days <= data0$vax_time + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1-90 & 
                             data0$vax_days_v1 <= data0$death_days &
                             data0$death_days <= data0$vax_time + 28 , ]
      cat("\t\t\tvaccinated persons observed at (first dose - 90 days) with ",event," in [vax_date; vax_date + 28 days]\n\n")
    }
    ############### 
    # from here only with event after vax1: vaccinated persons observed after the first dose with ",event," in [vax_date; vax_date + 28 days]
    if(i==6) {
      data        <- data0[data0$event==1 & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$event_days &
                             data0$event_days <= data0$vax_time + 28 , ]
      data_deaths <- data0[!is.na(data0$death_days) & !is.na(data0$vax_days_v1) &   
                             data0$start_obs <= data0$vax_days_v1 & 
                             data0$vax_days_v1 <= data0$death_days &
                             data0$death_days <= data0$vax_time + 28 , ]
      cat("\t\t\tvaccinated persons observed at the first dose with ",event," in [vax_date; vax_date + 28 days]\n\n")
    }
    
    nrow_data        <- nrow(data)
    nrow_data_deaths <- nrow(data_deaths)
    
    if(nrow_data==0       ) print("no events")
    if(nrow_data_deaths==0) print("no deaths")
    
    if(nrow_data>0 | nrow_data_deaths>0){
      
      flowchart_all[[i+1]] <- vector("list",length=5)
      names(flowchart_all[[i+1]]) <- c("all", "age30_50","age30","sexc","sex_age30")
      
      for(istrata_var in c("all", "age30_50","age30","sexc","sex_age30")){
        cat("\n\n*******************************************************\n")
        cat("*******************************************************\n")
        cat("*******************************************************\n")
        
        if(istrata_var =="all"){
          cat(paste0("\n\tSummary not stratified:\n\n\n"))
          if(nrow_data>0       ) data$strata_variable        <- "all"
          if(nrow_data_deaths>0) data_deaths$strata_variable <- "all"
        }
        else {
          cat(paste0("\n\tSummary stratified by variable '",istrata_var,"':\n\n\n"))
          if(nrow_data>0       ) data$strata_variable        <- data[       ,istrata_var]
          if(nrow_data_deaths>0) data_deaths$strata_variable <- data_deaths[,istrata_var]
        } 
        
        if(nrow_data>0       ) data_strata        <- data[       !is.na(data$strata_variable       ) & !is.na(data$id       ),]
        if(nrow_data_deaths>0) data_deaths_strata <- data_deaths[!is.na(data_deaths$strata_variable) & !is.na(data_deaths$id),]
        
        gc()
        
        if(nrow_data>0 & nrow(data_strata)>0){
          
          flowchart <- list()
          flowchart <- c(flowchart, n_ids = list(table1(data_strata[!duplicated(data_strata$id),c("strata_variable","id")][,"strata_variable"] )))
          
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
          flowchart <- c(flowchart, n_ids_per_event_vax_name_brand = list(table1( data_strata[, c("strata_variable","eventc","vax_name","vax_brand")])) )
          print(flowchart$n_ids_per_event_vax_name_brand)
          
          # age as continuous or integer:
          if(any(!is.na(data_strata$age))){
            
            cat(paste0("\n\nthe distribution of variable '",age,"':\n\n"))
            
            data_strata_age <- data_strata[!is.na(data_strata$age),]
            
            flowchart <- c(flowchart, summary_id_age                      = list(do.call("rbind",with(data_strata_age[!duplicated(data_strata_age[,"id"]), c("strata_variable","id","age")], tapply(age,    strata_variable                           , function(x)c(summary(x),n=length(x))) )) ))
            flowchart <- c(flowchart, summary_id_age_vax_name             = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age",         "vax_name"            )], tapply(age, paste(strata_variable,       vax_name          ), function(x)c(summary(x),n=length(x))) )) ))
            flowchart <- c(flowchart, summary_id_age_vax_brand            = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age",                    "vax_brand")], tapply(age, paste(strata_variable,                vax_brand), function(x)c(summary(x),n=length(x))) )) ))
            flowchart <- c(flowchart, summary_id_age_event                = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age","eventc"                       )], tapply(age, paste(strata_variable,eventc                   ), function(x)c(summary(x),n=length(x))) )) ))
            flowchart <- c(flowchart, summary_id_age_event_vax_name_brand = list(do.call("rbind",with(data_strata_age[, c("strata_variable","id","age","eventc","vax_name","vax_brand")], tapply(age, paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x))) )) ))
            
            rm(data_strata_age)
            
            print( flowchart$summary_id_age                      ); cat("\n")
            print( flowchart$summary_id_age_vax_name             ); cat("\n")
            print( flowchart$summary_id_age_vax_brand            ); cat("\n")
            print( flowchart$summary_id_age_event                ); cat("\n")
            print( flowchart$summary_id_age_event_vax_name_brand )  
          } else cat("\n\nno nonmissing age.\n\n")
          
          gc()
          
          # summaries for vaccination time and date:
          if(any(!is.na(data_strata$vax_time))){
            
            data_strata_vax <- data_strata[!is.na(data_strata$vax_time),]
            
            cat(paste0("\n\nthe distribution of the vaccination variable: '",vax_time,"':\n"))
            flowchart <- c( flowchart, summary_id_vax_time = list( do.call("rbind",with(
              data_strata_vax[, c("strata_variable","vax_n","id", "vax_time")], tapply( vax_time, paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
            print(flowchart$summary_id_vax_time)
            
            cat(paste0("\nthe distribution of the vaccination variable: '",vax_date,"':\n"))
            flowchart <- c(flowchart, summary_id_vax_date = list(do.call("rbind",with(
              data_strata_vax[, c("strata_variable","vax_n","id", "vax_date")], tapply( vax_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
            print(flowchart$summary_id_vax_date)
            
          } else cat("\n\nno vaccinated.\n\n")
          gc()
          
          ######################
          # event_days:
          if(any(!is.na(data_strata$event_days))){
            
            data_strata_vax <- data_strata[ !is.na(data_strata$event_days) & !is.na(data_strata$vax_time),]
            
            if( nrow(data_strata_vax)>0){
              if( any(( cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)<0 )) ){
                cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) before vaccination:\n"))
                flowchart <- c(  flowchart, summary_id_event_min_vax_before_vax = list( do.call("rbind",with(
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time","vax_name")]), tapply( (event_days-vax_time), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
                print(flowchart$summary_id_event_min_vax_before_vax)
              }
              if( any(( cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)>=0 )) ){
                cat(paste0("\n\nthe distribution of '",event,"_days' (days after vaccination) after vaccination:\n"))
                flowchart <- c(  flowchart, summary_id_event_min_vax_after_vax = list( do.call("rbind",with(
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time","vax_name")]), tapply( (event_days-vax_time), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
                print(flowchart$summary_id_event_min_vax_after_vax)
              }
            }
            
            cat(paste0("\n\nthe distribution of the '",event,"_days' variable:\n"))
            flowchart <- c( flowchart, summary_id_event_time = list( do.call("rbind",with(
              data_strata[!is.na(data_strata$event_days), c("strata_variable","vax_n","id", "event_days")], tapply( event_days, paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
            print(flowchart$summary_id_event_time)
            
            cat(paste0("\n\nthe distribution of the '",event,"_date' variable:\n"))
            flowchart <- c(flowchart, summary_id_event_date = list(do.call("rbind",with(
              data_strata[!is.na(data_strata$event_days), c("strata_variable","vax_n","id", "event_date")], tapply( event_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
            print(flowchart$summary_id_event_date)
            
            if( nrow(data_strata_vax)){
              if( any( (cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)<0 )) ){
                cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) before vaccination per ",event,", vaccine name and brand:\n"))
                flowchart <- c(flowchart, summary_id_event_min_vax_before_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time", "eventc","vax_name","vax_brand")]), tapply( (event_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
                print(  flowchart$summary_id_event_min_vax_before_vax_per_event_vax_name_brand)
              }
              if( any( (cond <- (data_strata_vax$event_days-data_strata_vax$vax_time)>=0 )) ){
                cat(paste0("\nthe distribution of the '",event,"_days' (days after vaccination) after vaccination per ",event,", vaccine name and brand:\n"))
                flowchart <- c(flowchart, summary_id_event_min_vax_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
                  unique(data_strata_vax[cond, c("strata_variable","id", "event_days","vax_time", "eventc","vax_name","vax_brand")]), tapply( (event_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
                print(  flowchart$summary_id_event_min_vax_after_vax_per_event_vax_name_brand)
              }
            }
            rm(data_strata_vax); rm(data_strata)
          } else cat("\n\nno events.\n\n")
        } # if nrow(data_strata) > 0
        else cat("\nno events.\n\n")
        
        gc()
        
        ######################
        # deaths:
        if(nrow(data_deaths)>0 & nrow(data_deaths_strata)>0){
          
          cat(paste0("\n\nthe distribution of the 'death_date' variable:\n"))
          flowchart$deaths <- c( summary_id_death_date = list(do.call("rbind",with(
            data_deaths_strata[, c("strata_variable","vax_n","id", "death_date")], 
            tapply( death_date, paste0(strata_variable," vax_n:",vax_n), function(x)c(as.character(summary(x)),n=as.character(length(x)))) )) ))
          print(flowchart$deaths$summary_id_death_date)
          
          cat(paste0("\n\nthe distribution of the 'death_days' variable:\n"))
          flowchart$deaths <- c( flowchart$deaths, summary_id_death_time = list( do.call("rbind",with(
            data_deaths_strata[, c("strata_variable","vax_n","id", "death_days")], 
            tapply( death_days, paste0(strata_variable," vax_n:",vax_n), function(x)c(summary(x),n=length(x))) )) ))
          print(flowchart$deaths$summary_id_death_time)
          
          data_deaths_strata <- data_deaths_strata[!is.na(data_deaths_strata$vax_time),]
          
          if(nrow(data_deaths_strata)>0){
            cat(paste0("\n\nthe distribution of 'death_days': days after vaccination:\n"))
            flowchart$deaths <- c( flowchart$deaths, summary_id_death_after_vax = list( do.call("rbind",with(
              data_deaths_strata[, c("strata_variable","id", "death_days","vax_time","vax_name")], 
              tapply( (death_days-vax_time), paste(strata_variable,vax_name), function(x)c(summary(x),n=length(x))) )) ))
            print(flowchart$deaths$summary_id_death_after_vax)
            
            cat(paste0("\nthe distribution of the 'death_days' (days after vaccination) per ",event,", vaccine name and brand:\n"))
            flowchart$deaths <- c(flowchart$deaths, summary_id_death_after_vax_per_event_vax_name_brand = list( do.call("rbind",with( 
              data_deaths_strata[, c("strata_variable","id", "death_days","vax_time", "eventc","vax_name","vax_brand")], 
              tapply( (death_days-vax_time), paste(strata_variable,eventc,vax_name,vax_brand), function(x)c(summary(x),n=length(x)) ) ))  ) )
            print(  flowchart$deaths$summary_id_death_after_vax_per_event_vax_name_brand)
          }
        } # end if nrow(data_deaths_strata)>0
        else cat("\n\nno deaths.\n\n")   
        
        gc()
        flowchart_all[[i+1]][[istrata_var]] <- list(flowchart)
        
      }  #end for istrata_var
      
      
    }  # end if nrow(data)>0 | nrow(data_deaths)>0
    
  } # end 'i'
  
  attributes(flowchart_all) <- c( attributes(flowchart_all), variables=list( event=event, vax_name=vax_name, id=id, start_obs=start_obs, 
                                                                             vax_time=vax_time, vax_date=vax_date, 
                                                                             death_date=death_date, death_days=death_days, 
                                                                             age=age, age30_50=age30_50,  age30=age30,     
                                                                             sex=sex, sexc=sexc, sex_age30=sex_age30,
                                                                             condition_value=condition_value )  )
  cat("\n\nattributes:\n\n")
  print(attributes(flowchart_all))
  
  if(!missing(path_file_name)){  
    sink()
    path_file_name <- paste0( substring(path_file_name,1,nchar(path_file_name)-3), "RData")
    save( flowchart_all, file=path_file_name )
  }
  
  cat(paste0(": duration = ",format(difftime(Sys.time(),sys_time))," (till ",Sys.time(),")\n"))  
  
  invisible(flowchart_all)
  
} # the end of function 'characteristics'


