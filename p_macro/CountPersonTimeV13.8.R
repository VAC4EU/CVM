
#Author: Roel Elbers MSc.
#email: r.j.h.elbers@umcutrecht.nl
#Organisation: UMC Utrecht, Utrecht, The Netherlands
#Date: 15/07/2021

CountPersonTime <- function(Dataset_events = NULL, Dataset, Person_id, Start_study_time, End_study_time, Start_date,
                            End_date, Birth_date = NULL,Rec_period = NULL, Strata = NULL,Outcomes_nrec = NULL,
                            Outcomes_rec = NULL, Name_event = NULL, Date_event = NULL, Age_bands = NULL,
                            Unit_of_age = "year" , Increment = "year", include_remaning_ages = T, Aggregate = T,
                            print = F, check_overlap = T, save_intermediate = NULL, load_intermediate = F){
  
  if(print) print("Version 13.8")
  # Check if demanded R packages are installed, install if not,  and activate
  ################################################################################################################################
  if(print) print("Check packages data.table and lubridate")
  if (!require("data.table")) install.packages("data.table")
  library(data.table)
  
  if (!require("lubridate")) install.packages("lubridate")
  library(lubridate)
  ################################################################################################################################
  
  
  if(any(Outcomes_rec %in% Outcomes_nrec)){stop("Overlapping event names for Outcomes_rec and Outcomes_nrec")}
  
  if (!is.logical(load_intermediate)) {stop("Parameter 'load_intermediate' accepts only logical constants")}
  if (!missing(save_intermediate)) {
    tryCatch(dirname(save_intermediate), error = function(e)
      stop("Parameter 'save_intermediate' accepts only valid filepaths with already existing folders"))
  }
  
  if (load_intermediate) {
    if (!missing(save_intermediate) && file.exists(save_intermediate)) {
      load(save_intermediate)
    } else {
      Dataset <- as.data.table(Dataset)
    }
  } else {
    
    #Set character input for study dates to date format
    ################################################################################################################################
    if(print) print("Assign date format to Start_study_time and End_study_time")
    Start_study_time<-as.IDate(as.character(Start_study_time),"%Y%m%d")
    End_study_time<-as.IDate(as.character(End_study_time),"%Y%m%d")
    ################################################################################################################################
    
    #
    ################################################################################################################################
    end_date_new <- as.Date(ifelse(End_study_time <= Sys.Date(), Sys.Date(), End_study_time + 1),origin = "1970-01-01")
    ################################################################################################################################
    
    #check if study start and stop dates are valid
    ################################################################################################################################
    
    if(print) print("Check if Start_study_time and End_study_time are valid")
    
    if(!sum(Start_study_time==seq.Date(as.Date("19000101","%Y%m%d"),Sys.Date(),by = Increment))==1){
      
      if(Increment == "year"){stop("Change the start date to the first of january. Wrong study start date can produce invalid results.")}
      if(Increment == "month"){stop("Change the start date to the first of month. Wrong study start date can produce invalid results.")}
      if(Increment == "week"){stop("Change the start date to a monday. Wrong study start date can produce invalid results.")}
      
    }
    
    if(!sum(End_study_time==seq.Date(as.Date("19000101","%Y%m%d"),end_date_new ,by = Increment)-1)==1){
      
      if(Increment == "year"){stop("Change the end date to the 31th of december. Wrong study start date can produce invalid results.")}
      if(Increment == "month"){stop("Change the end date to the last day of the month. Wrong study start date can produce invalid results.")}
      if(Increment == "week"){stop("Change the end date to a sunday. Wrong study start date can produce invalid results.")}
      
    }
    
    gc()
    
    ################################################################################################################################
    
    date_cols<-c(Start_date,End_date,Birth_date)
    
    # Reduce memory size using integers
    Dataset[, c(date_cols) := lapply(.SD, as.IDate), .SDcols = date_cols]
    #Dataset[, c(Strata) := lapply(.SD, as.integer), .SDcols = Strata]
    gc()
    
    
    #Check if start, end and birth date are all filled. If end date is not filled it will be replaced by the and study date
    ################################################################################################################################
    
    
    if(print) print("Check if date columns in input data are valid and in correct order")
    if(sum(is.na(Dataset[,.(get(Start_date))]))>0){stop("Empty start dates")}
    if(!is.null(Age_bands)){if(sum(is.na(Dataset[,.(get(Birth_date))]))>0){stop("Empty birth dates")}}
    if(sum(is.na(Dataset[,.(get(End_date))]))>0){print(paste0(sum(is.na(Dataset[,.(get(End_date))]))," empty end dates will be filled with the end study date. This may cause overlapping intervals"))}
    Dataset[is.na(get(End_date)),eval(End_date) := End_study_time]
    
    gc()
    
    #Check the order of dates
    ################################################################################################################################
    wrong_End_date<-nrow(Dataset[get(Start_date)>get(End_date),])
    if (wrong_End_date>0){warning(paste0(wrong_End_date," end date(s) prior to start date"))}
    wrong_Start_date<-nrow(Dataset[get(Start_date)>Sys.Date(),])
    if (wrong_Start_date>0){warning(paste0(wrong_Start_date," start date(s) in future"))}
    
    if(!is.null(Age_bands)){
      wrong_Birth_date<-nrow(Dataset[get(Start_date)<get(Birth_date),])
      if (wrong_Birth_date>0){warning(paste0(wrong_Start_date," start date(s) before birth date"))}}
    ################################################################################################################################
    
    #Check if birthdays are unique
    ################################################################################################################################
    #if(!is.null(Age_bands)){
    # if(nrow(Dataset[,uniqueN(get(Birth_date)), by = Person_id][V1>1,])!=0){stop("Persons with several birth dates") }
    #}
    ################################################################################################################################
    
    #Check if the subjects have overlap in the time intervals (within strata???), defined by end-start date.
    ################################################################################################################################
    
    if(check_overlap){
      
      if(print) print("Check if observation periods do not have overlap")
      test_overlap<-Dataset[!is.na(get(End_date))&!is.na(get(End_date))&get(End_date)>get(Start_date),][,.(get(Person_id), as.integer(get(Start_date)), as.integer(get(End_date)))]
      setkey(test_overlap,V1,V2,V3)
      test_overlap2<-as.data.table(foverlaps(test_overlap, test_overlap, type="any", which=TRUE))
      test_overlap2<-test_overlap2[xid!=yid,]
      test_overlap[,id:=as.integer(rownames(test_overlap))]
      overlap_subjects<-unlist(unique(test_overlap2[test_overlap, on = .(xid = id), nomatch=NULL][,.(V1)]))
      
      if(length(overlap_subjects) > 0){
        warning("Subjects have overlapping person time: ")
        warning(paste0(overlap_subjects," "))
        stop("Invalid results could be generated so the script is stopped")
      }
      
      rm(test_overlap,test_overlap2,overlap_subjects)
      gc()
      
    }
    ################################################################################################################################
    
    
    #Select relevant data
    ################################################################################################################################
    intv <- as.IDate(c(Start_study_time, End_study_time))
    Dataset <- Dataset[get(Start_date) %between% intv|get(End_date) %between% intv|(get(Start_date) < Start_study_time & get(End_date) > End_study_time)] 
    
    if(nrow(Dataset) == 0){
      Dataset <- NULL
      if(print) print("No subjects with any observation time within studyperiod. NULL is returned")
      return(Dataset)
    }
    
    Dataset[get(Start_date) < Start_study_time,eval(Start_date) := Start_study_time]
    Dataset[get(End_date) > End_study_time,eval(End_date) := End_study_time]
    
    start <- Start_study_time
    end <- End_study_time
    Time_Increment <- as.IDate(seq.Date(start, end, Increment))
    
    if(length(Time_Increment) > 1){
      Time_Increment_end <- as.IDate(seq.Date(Time_Increment[2], by = Increment, length = length(Time_Increment)))-1
    }else{
      Time_Increment_end <- end
    }
    
    
    #Determine the ages at the beginning and end of all observation periods. Output is a starting point for calculation and splitting of
    # age bands
    ################################################################################################################################
    if(!is.null(Age_bands)){
      
      if(print) print(paste0("Calculate ages at the start and end of every observation period by ", Increment))
      
      if (nrow(Dataset) > 0){
        
        Dataset[, age_start := floor(time_length(interval(get(Birth_date), get(Start_date)), Unit_of_age)) ]
        Dataset[, age_end := floor(time_length(interval(get(Birth_date), get(End_date)), Unit_of_age)) ]
        
      } else{
        Dataset[,age_start := NA]
        Dataset[,age_end := NA]
      }   
      
    }
    ################################################################################################################################
    
    #Calculate agebands in 2 steps ((1)split/recalculate start/end ages and (assign row to ageband) )
    ################################################################################################################################
    if(!is.null(Age_bands)){
      if(print) print("Create agebands")
      if(nrow(Dataset) > 0){
        #### New code from version 13.6
        
        #Produce a dataset with Agebands and the start and end age of that ageband. This can be used to merge top all cases in the Dataset that overlap with the start and end age. 
        Agebands_list <- list()
        
        for (k in 1:length(Age_bands)){
          
          if( k == 1) Agebands_list[[k]] <- paste0(Age_bands[k],"-",Age_bands[k+1])
          if( k > 1 &k!= length(Age_bands)) Agebands_list[[k]] <- paste0(Age_bands[k]+1,"-",Age_bands[k+1])
          if( k== length(Age_bands) & include_remaning_ages == T) Agebands_list[[k]] <- paste0(Age_bands[k]+1,"+")
          
        }
        
        Agebands_list <- as.data.table(do.call(rbind, Agebands_list))
        colnames(Agebands_list)<- "Ageband"
        
        Agebands_list[,row := row.names(Agebands_list) ]
        Agebands_list[,ST := as.numeric(gsub("[^[:digit:].]", "\\1",strsplit(as.character(Ageband),"-")[[1]][1])),by = row ]
        Agebands_list[,EN := as.numeric(gsub("[^[:digit:].]", "\\1",strsplit(as.character(Ageband),"-")[[1]][2])),by = row ]
        Agebands_list[is.na(EN),EN := 4000 ]
        
        #Merge the overlapping
        setkeyv(Dataset, c("age_start","age_end"))
        Dataset <- foverlaps(Agebands_list, Dataset, by.x = c("ST","EN"), nomatch = 0L, type = "any")
        
        # select the rows that doubled by the merge. In these, multiple agebands occur witing the obeservation period. So start and end dated need to be adapted
        Dataset <- Dataset[, row := row.names(Dataset)]
        Dataset <- Dataset[age_start < ST  ,eval(Start_date) := as.IDate(add_with_rollback(get(Birth_date), period(ST,units = Unit_of_age), roll_to_first = T, preserve_hms = T)), by = row]
        Dataset <- Dataset[age_end > EN  ,eval(End_date) := as.IDate(add_with_rollback(get(Birth_date), period(EN + 1,units = Unit_of_age), roll_to_first = T, preserve_hms = T)) - 1, by = row]
        Dataset <- Dataset[,':=' (age_start = NULL, age_end = NULL,ST = NULL, EN = NULL, row = NULL)]
        }
    }
    
    
  
  if(!is.null(Age_bands)){Age_band_coln<-"Ageband"} else Age_band_coln<-"Ageband"<-NULL
  
  
  #Enlarge table by time increment. 
  ################################################################################################################################
  
  
  if(print) print(paste0("Transform input date to a dataset per ", Increment, ". This step increases the size of the file with respect to the choosen increment" ))
  
  Dummy <- as.data.table(cbind(Increment = Time_Increment,End = Time_Increment_end))
  Dummy <- Dummy[, Increment := as.IDate(Increment)]
  Dummy <- Dummy[, End := as.IDate(End)]
  colnames(Dummy) <- c(Increment,"End")
  
  
  setkeyv(Dataset, c(Start_date, End_date))
  Dataset <- foverlaps(Dummy, Dataset, by.x=c(Increment, "End"), nomatch = 0L, type = "any")
  
  Dataset <- Dataset[get(Start_date) <= get(Increment) & get(End_date) >= get(Increment),eval(Start_date) := get(Increment)]
  Dataset <-Dataset[get(End_date) >= End & get(Start_date) <= End, eval(End_date) := End][, End := NULL]
  
  rm(Dummy)
  gc()
  
  ################################################################################################################################
  
  if (!missing(save_intermediate)) {
    save(Dataset, file = save_intermediate)
  } 
  
  }
  
  
  #NEW CODE V11 If recurrent events is true. This is a whole different approach compared to the situation where only the first
  # event is used. When joining multiple doubling will occur. I choose to do not do this and choose a method only allowing joins 
  # that have unique combinations
  ################################################################################################################################
  
  #Make boolean for reccurent and not recurrent events that are needed to decide which codes need to be applied
  
  if (is.null(Outcomes_rec)) {Rec_events <- F} else{Rec_events <- T}
  if (is.null(Outcomes_nrec)) {nRec_events <- F} else{nRec_events <- T}
  
  #Combine two input datasets to one dataset if only the first occurrence of an event is evaluated or exclude events that are following 
  #the previous event within a time period. 
  ################################################################################################################################
  
  #Create a boolean for if there are any events at all so we can remove Dataset_events and spare RAM
  if(!(is.null(Dataset_events) | (is.null(Outcomes_rec) & is.null(Outcomes_nrec)))){Dataset_events_b <- T} else{Dataset_events_b <- F}

  
  if (Dataset_events_b) {
  
  #Make a vector for all events that can be used to processes that are independent of recurrent or not  
  Outcomes <- c(Outcomes_nrec, Outcomes_rec)
  
  Dataset_events <- copy(Dataset_events)    
  
  #NEW CODE FROM V11
  
  if(nrow(Dataset_events) > 0){
    Dataset_events <- Dataset_events[get(Name_event) %in% Outcomes,]
    
    # set date to integer to spare RAM
    Dataset_events[, c(Date_event) := lapply(.SD, as.IDate), .SDcols = Date_event]
    
    #Add a column that represents the occurence of the event. This is needed to select the first 1 if not recurrent and to see which events are removed if recurrent
    setorderv(Dataset_events,c(Person_id,Name_event,Date_event))
    Dataset_events[,Recurrent := cumsum(!is.na(get(Date_event))),by=c(Person_id,Name_event)]
    
    #Split events data in recurrent and not recurrent and select only the outcomes that are actually in the data to prevent processes that are not needed
    if(print) print("Split outcomes in event table in recurrent and not recurrent")
    Dataset_events_nrec <- copy(Dataset_events)[get(Name_event) %in% Outcomes_nrec,]
    colls_outcomes_nrec <- Outcomes_nrec[Outcomes_nrec %in% Dataset_events_nrec[,get(Name_event)]]
    Dataset_events_rec <- copy(Dataset_events)[get(Name_event) %in% Outcomes_rec ,]
    colls_outcomes_rec <- Outcomes_rec[Outcomes_rec %in% Dataset_events_rec[,get(Name_event)]]
   
    
    rm(Dataset_events)
    gc()
    
    if(Rec_events & nrow(Dataset_events_rec) > 0){
      if(print) print("If Rec_events = T then determine the censoring periods. This is a looped procedure that may need some improvement still")
      it=1
      if(!length(Rec_period)==length(Outcomes_rec)) stop("the vectors Outcomes and rec period have different lengths")
      
      events_rec_list <- copy(Dataset_events_rec[0])
      
      #Split events on if they are in data and have a rec_period > 0. If they are not in the data or have a Rec_period of 0 they do not need to be in all the processes.
      #Only the count of the events is mandatory if a Rec_period = 0. The while loop to remove events and the procedure to subtract persontime are not needed. 
      Outcomes_rec1 <- colls_outcomes_rec[colls_outcomes_rec  %in% Outcomes_rec[!Rec_period == 0]]  
      Outcomes_rec2 <- colls_outcomes_rec[!colls_outcomes_rec %in% Outcomes_rec1] 
      #Make a vector with only the rec_period for the outcomes in Outcomes_rec1. This in needed along with Outcomes_rec1 in the row 313 untill 352 
      Rec_period1 <- Rec_period[Outcomes_rec %in% Outcomes_rec1]
      
      #Delete Rec_period to be sure it is not used anymore further in the script. It should be matched with Outcomes_rec1  
      rm(Rec_period)
      
      
      #Split the events in 2 datasets 
      Dataset_events_rec2  <- copy(Dataset_events_rec)[get(Name_event) %in% Outcomes_rec1]
      Dataset_events_rec3  <- copy(Dataset_events_rec)[get(Name_event) %in% Outcomes_rec2]
      
      if(length(Outcomes_rec1) > 0){
            for (i in 1:length(Outcomes_rec1)){
              
              if(print) print(paste("Remove ",Outcomes_rec1[i], "outcomes within a censoring of ",Rec_period1[i]," days"))
              events_rec  <- copy(Dataset_events_rec2)[get(Name_event) == Outcomes_rec1[i],]
              
              
              while(nrow(events_rec) > 0){ 
                
                setorderv(events_rec,c(Person_id,Name_event,Date_event))
                events_rec <- events_rec[,D := shift(get(Date_event)),by = c(Person_id,Name_event) ]
                events_rec <- events_rec[,dif := get(Date_event) - D]
                events_rec <- events_rec[is.na(dif), dif := 0 ][,dif := as.numeric(dif)]
                events_rec <- events_rec[,cumdif := cumsum(dif),by = c(Person_id,Name_event)]
                
                events_rec2 <-events_rec[ cumdif <= Rec_period1[i],]
                setorderv(events_rec2,c(Person_id,Name_event,Date_event))
                events_rec2 <-events_rec2[,.SD[1],by = c(Person_id,Name_event)][,Iteration := it]
                
                events_rec_list <- rbindlist(list(events_rec_list, events_rec2),fill=T)
                rm(events_rec2)
                
                events_rec <- events_rec[cumdif > Rec_period1[i],]
                
                lapply(c("dif","cumdif","D"), function(x){events_rec <-events_rec[,eval(x) := NULL]})
                if(print) print(paste0("Cycle ",it))
                it=it+1
                gc()
              }
              rm(events_rec)
              gc()
            }
      
      
      Dataset_events_rec <- events_rec_list
      rm(events_rec_list,Dataset_events_rec2)
      gc()
      
      for(i in 1:length(Rec_period1)){
        
        if(print) print(paste("Set censoring start and end dates for outcome ",Outcomes_rec1[i], " with a duration of ",Rec_period1[i]," days"))
        Dataset_events_rec <- Dataset_events_rec[dif!=0 & get(Name_event) == Outcomes_rec1[i] & dif < Rec_period1[i],Delete :=T ][is.na(Delete),]
        Dataset_events_rec <- Dataset_events_rec[get(Name_event) == Outcomes_rec1[i], ":=" (RecSTDT=get(Date_event),RecENDT=get(Date_event) + Rec_period1[i])]  
        gc()
        
      }
      
      
      
      lapply(c("dif","cumdif","D","Iteration","Delete"), function(x){Dataset_events_rec <- Dataset_events_rec[,eval(x) := NULL]})
      }else{
        Dataset_events_rec <- Dataset_events_rec3[0]
      }
      
      setnames(Dataset_events_rec, eval(Person_id), "Person_id")
      
      setnames(Dataset_events_rec3, eval(Person_id), "Person_id")
      
    }  
      
    if(nRec_events & nrow(Dataset_events_nrec) > 0){ 
      
      if(print) print("If Rec_events = F then selecting only the first event")
      Dataset_events_nrec <- Dataset_events_nrec[Recurrent==1,]
      Dataset_events_nrec<-dcast(Dataset_events_nrec, get(Person_id) + Recurrent ~ get(Name_event), value.var = eval(Date_event))
      setcolorder(Dataset_events_nrec,neworder = c('Person_id','Recurrent',colls_outcomes_nrec))
      
      setkeyv(Dataset,Person_id)
      setkey(Dataset_events_nrec,Person_id)
      Dataset <- Dataset_events_nrec[Dataset,][,Recurrent := NULL]
      setnames(Dataset, "Person_id", eval(Person_id))
      
    }
    rm(Dataset_events_nrec)
    gc()
    
    
    lapply(Outcomes_nrec, function(x) if (!x %in% colnames(Dataset)) Dataset <- Dataset[, eval(x) := as.IDate(NA, format = "%d%m%Y")]) 
  } else{
    invisible(lapply(Outcomes, function(x) Dataset <- Dataset[, eval(x) := as.IDate(NA, format = "%d%m%Y")]))
    Dataset_events_nrec <- copy(Dataset_events)
    Dataset_events_rec <- copy(Dataset_events)
    Dataset_events_rec3 <- copy(Dataset_events)
  }
  
  ################################################################################################################################
  
  if(Rec_events){
    if(nrow(Dataset_events_rec) > 0 | nrow(Dataset_events_rec3) > 0){
    
    #Dataset[,row := row.names(Dataset)]
    setnames(Dataset, eval(Person_id), "Person_id")
      
    #Calculate for every observation period the amount of persontime that needs to subtracted for every event
    ################################################################################################################################
    
    Dataset_events_rec2 <- copy(Dataset_events_rec)
    if(nrow(Dataset_events_rec2) > 0){
    
      if(print) print("Calculate days to subtract from persontime for recurrent events")
      
      Dataset_temp <- copy(Dataset)[,c("Person_id",Start_date,End_date), with = F]
      
      setkeyv(Dataset_temp, c("Person_id",Start_date, End_date))
      Dataset_events_rec2 <- foverlaps(Dataset_events_rec2, Dataset_temp, by.x = c("Person_id","RecSTDT","RecENDT"), nomatch = 0L, type = "any")
      
      ###Bug fix, it can happen that recurrent events are without the study period minus the rec_period. It would also be an option to get rid of them on the event
      ###File. However, I advise to do that before the function if wanted
      
      if(nrow(Dataset_events_rec2) > 0){
      
      Dataset_events_rec2 <- Dataset_events_rec2[, row2 := row.names(Dataset_events_rec2)]
      Dataset_events_rec2 <- Dataset_events_rec2[,':=' (start_date2 = max(get(Start_date),RecSTDT,na.rm=T), end_date2 = min(get(End_date)+1,RecENDT,na.rm=T)), keyby = row2]     
      Dataset_events_rec2 <- Dataset_events_rec2[,SUBTR := (as.numeric(end_date2)-as.numeric(start_date2)),by = row2]
      
      Dataset_events_rec2 <- Dataset_events_rec2[, var := paste0("SUBTRCUM_", get(Name_event))]
      Dataset_events_rec2 <- dcast(Dataset_events_rec2, Person_id + get(Start_date) + get(End_date) ~  var, value.var = "SUBTR", fun.aggregate = sum)
      setnames(Dataset_events_rec2, "Start_date",eval(Start_date))
      setnames(Dataset_events_rec2, "End_date",eval(End_date))
      
      Dataset <- merge(x = Dataset, y = Dataset_events_rec2 , by = c("Person_id", eval(Start_date), eval(End_date)), all.x = T , all.y = F, allow.cartesian = F)
      }else{
        
        print("All recurrent events fall witout the period of interest")
      }
      
      rm(Dataset_temp)  
      gc()
    }
    
    # else{
    #   add <- paste0("SUBTRCUM_",Outcomes_rec)                 
    #   lapply(add, function(x){Dataset <- Dataset[,eval(x) := 0]})
    #   rm(add)
    # }
    # 
    rm(Dataset_events_rec2)
    gc()
    #Check if all outcomes are present. If not add column with value 0
    
    #if(length(Outcomes[!Outcomes %in% unique(p[["OUTC"]])]) > 0)lapply(paste0("SUBTRCUM_",Outcomes[!Outcomes %in% unique(p[["OUTC"]])]), function(x){Dataset <-Dataset[,eval(x) := 0]})
    SUBTRCUM_MISSING <- Outcomes_rec[!paste0("SUBTRCUM_",Outcomes_rec) %in% unique(colnames(Dataset))] 
    if(length(SUBTRCUM_MISSING) > 0)lapply(paste0("SUBTRCUM_",SUBTRCUM_MISSING), function(x){Dataset <- Dataset[,eval(x) := 0]})
  
    ################################################################################################################################
    
    #Calculate the number of events per observation period
    ################################################################################################################################
    #k <- c("Person_id",Start_date,End_date)
    
    if(print) print("Count the number of events per subject for recurrent events")
    Dataset_temp <- copy(Dataset)[,c("Person_id",Start_date,End_date), with = F]
    Event_temp <- rbind(copy(Dataset_events_rec)[,c("Person_id",Date_event,Name_event), with = F],copy(Dataset_events_rec3)[,c("Person_id",Date_event,Name_event), with = F])[, Dummy := get(Date_event)]
    
    setkeyv(Dataset_temp, c("Person_id",Start_date, End_date))
    Dataset_temp <- foverlaps(Event_temp,Dataset_temp, by.x=c("Person_id",Date_event, "Dummy"), nomatch = 0L, type = "any")
    
    
    if(nrow(Dataset_temp) > 0){
    Dataset_temp <- Dataset_temp[,.("b" = .N),by = c("Person_id",Start_date,End_date, Name_event)][,var := paste0(get(Name_event),"_b")]
    Dataset_temp <- dcast(Dataset_temp, Person_id + get(Start_date) + get(End_date) ~  var, value.var = "b", fill = 0)
    setnames(Dataset_temp, c("Start_date","End_date"), c(eval(Start_date), eval(End_date) ))
    
    Dataset <- merge(x = Dataset, y = Dataset_temp, by = c("Person_id", Start_date, End_date), allow.cartesian = F, all.x = T, all.y = F)
    }
    
    rm(Event_temp,Dataset_temp)
    gc()
    

    rm(Dataset_events_rec)
    gc()
    
    B_MISSING <- Outcomes_rec[!paste0(Outcomes_rec, "_b") %in% unique(colnames(Dataset))]
    if(length(B_MISSING) > 0)lapply(paste0(B_MISSING,"_b"), function(x){Dataset <-Dataset[,eval(x) := 0]})
    
    
    ################################################################################################################################   
    }else{
      lapply(Outcomes_rec, function(x){Dataset <- Dataset[,eval(paste0(x, "_b")) := 0]})
      lapply(Outcomes_rec, function(x){Dataset <- Dataset[,eval(paste0("SUBTRCUM_",Outcomes_rec)) := 0]})
      
    }
  }
  }
  ################################################################################################################################
  
  # If aggregate is TRUE create columns for aggregation. I do this earlier in the program, so sparated from the aggregation so that I can delete Dataset_events earlier for 
  #RAM optimalistation
  if (Aggregate) {
    
    if (Dataset_events_b){
      PT_colls <- c("Persontime",paste0("Persontime_",Outcomes),paste0(Outcomes,"_b"))
    }else{
      PT_colls <- "Persontime"
    }
    
    if (!is.null(Age_bands)){
      by_colls <- c(Strata, Increment, "Ageband")
    }else{
      by_colls <- c(Strata, Increment)
    }
  }
  ################################################################################################################################  
  
  #Calculate persontimes
  ################################################################################################################################
  
  if(print) print("Calculate persontimes")
  
  if(Dataset_events_b){
  
  #rm(Dataset_events)  
  #gc()    
    
  Outcomes_b <- paste0(Outcomes, "_b")
  Persontime_Outcomes <- paste0("Persontime_", Outcomes)
  #sort_order <- c(eval(Person_id), eval(Start_date))
  sort_order <- c(eval(Person_id), eval(Start_date),Age_band_coln,eval(Strata))
  coln <- c(eval(Person_id), eval(Strata), Age_band_coln, eval(Increment), "Persontime", eval(Persontime_Outcomes), eval(Outcomes_b))
  Dataset[,Persontime := .(get(End_date)-get(Start_date) + 1)]
        
          if(nRec_events){
          lapply(Outcomes_nrec,function(x)Dataset[,paste0(eval(x),"_b") := fifelse(!is.na(get(x)) & get(x) %between% list(get(Start_date),get(End_date)),1,0)])
          lapply(Outcomes_nrec, function(x) Dataset[,paste0("Persontime_",x) := fifelse(!is.na(get(x)) & get(x) < get(Start_date), 0, Persontime)])
          #lapply(Outcomes,function(x)Dataset[get(x) %between% list(get(Start_date),get(End_date)),`:=`(paste0(eval(x),"_b")= 1, paste0("Persontime_",x) = .(get(x)-get(Start_date)+1))])
          lapply(Outcomes_nrec,function(x)Dataset[get(x) %between% list(get(Start_date),get(End_date)),paste0("Persontime_",x) := .(get(x)-get(Start_date)+1)])
          }
        
        
        if(Rec_events){
          if("Person_id" %in% colnames(Dataset)) setnames(Dataset,"Person_id",eval(Person_id))
          
          #By merge NA's are introduced. Removed but not found a good method yet to do this in 1 step within data.table
          x <- c(paste0("SUBTRCUM_",Outcomes_rec),paste0(Outcomes_rec, "_b"))
          Dataset <- as.data.frame(Dataset)
          Dataset[,x][is.na(Dataset[,x])] <- 0
          Dataset <- as.data.table(Dataset)
          rm(x)
          
          lapply(Outcomes_rec,function(x)Dataset[,paste0("Persontime_",x) := .(get(End_date) - get(Start_date) + 1 - get(paste0("SUBTRCUM_",x)))])
        }
        
  } else {
    
    
    Dataset[,Persontime := .(get(End_date)-get(Start_date) + 1)]
    sort_order <- c(eval(Person_id), eval(Start_date),Age_band_coln,eval(Strata))
    coln <- c(eval(Person_id), eval(Strata), Age_band_coln, eval(Increment),"Persontime")
    
  }
          
  
  ################################################################################################################################
  
  #Create output table
  ################################################################################################################################
  
  if(print) print("Create output table")
  if(Increment=="month"){Dataset[,eval(Increment) := substr(get(Increment),1,7)]}
  if(Increment=="year"){Dataset[,eval(Increment) := substr(get(Increment),1,4)]}
  
  setorderv(Dataset, sort_order)
  Dataset <- Dataset[,coln,with=FALSE]
  
  
  #Aggregate based on PT_colls and by_colls created earlier  
  if (Aggregate) {    
    Dataset <- Dataset[, lapply(.SD, sum), .SDcols = PT_colls, by = by_colls]
    rm(PT_colls,by_colls)

  }
  return(Dataset)
  ################################################################################################################################
  rm(Dataset)
  gc()
  
  
  
}






