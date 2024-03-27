# Program Information  ----------------------------------------------------
#
# Program:      Myo-SCRI-01-T5_20_.R # in the 'steps' directory 
# Author:       Svetlana Belitser
# Description:  calls functions to run SCRI analyses on a specified dataset 
#               runs on all datasets in g_output/scri                  
# Requirements: 
#               dependencies: preceding steps, package "survival" 
#               input:   g_intermediate/data_vax_cohort.RData  ==> 'vax_data' dataset
#               output:  g_output/scri/*  
#
#               parameters: in 08_SCRI_parameters.R 
#  
#               function: p_macro/scri_tools.R                    
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
  
  
  cat(paste0('\n\t"',subpop,'":\n\n'))
  
  
  # Import Data -------------------------------------------------------------
  
  # Load dataset 'data_vax' by loading file 'data_vax_cohort.RData'
  load(file=paste0(dirtemp_cohort, "data_vax_for_matching", suffix[[subpop]], ".RData"))
  load(file=paste0(dirtemp_cohort,"ids_study_all_matching", suffix[[subpop]], ".RData"))
  
  data_vax <- data_vax_for_matching
  
  
  
  #############   cohort models ############################
  #
  #
  old_width = options(width=300)
  

 
  
  for(iae in ae_events){
    
    cat(paste("\n",iae,":\n\n"))
    
    
    
       print("main part")
      
      # cohort output_directory for the event:  EXPORT 
      sdr <- paste0(sdr0, iae,"/")
      dir.create(sdr, showWarnings = FALSE, recursive = TRUE)
      
      # cohort output_directory for the event:  LOCAL 
      sdr_models <- paste0(sdr_models0, iae,"/")
      dir.create(sdr_models, showWarnings = FALSE, recursive = TRUE)
      
      # cohort output_directory for the event:  EXPORT 
      sdr_cohort <- paste0(sdr,"cohort_causal/")
      dir.create(sdr_cohort, showWarnings = FALSE, recursive = TRUE)
      
      # cohort output_directory for the event:  LOCAL 
      dirtemp_event <- paste0(dirtemp_cohort,iae,"/")
      dir.create(dirtemp_event, showWarnings = FALSE, recursive = TRUE)
      
      
  
      
  
      #########################################
      #
      #   matching vax-independent factors:   "sex","immunosup","number_conditions"  
      #
      #   matching vax-dependent factors:     "age":    {abs(age_vax - age_unvax)[index_dat]<=2years}
      #                                       "covid":  {abs(last_covid_date_before_index_date_vaxxed - last_covid_date_before_index_date_unvax) <= 30days}
      #
      #   matching vax-dependent conditions:  "myocarditis" (or other event): {(index_date - last_myo_before_index_date) > 1 year}
      #                                       "covid":                        {(last_covid_date_before_vax - index_date) > 30days}
      #
      #
      
      match_name          <- paste0("match_info_sex_age_covid_2comorb_",iae)
      match_indep_factors <- c( "sex","immunosup","number_conditions_4" )
      match_dep_factors   <- c( "age","covid",iae)
      
      # create overlapping age categories, because abs(age_dif)<=2years:
      age_groups_start <- c( -0.01, 20, 30, 40, 50, 60, 70  )
      age_groups_stop  <- c(    22, 32, 42, 52, 62, 72, 150 )
      
      # event_dataset
      load( paste0(dirtemp_cohort, "event_dataset_", iae, suffix[[subpop]], ".RData") )
      event_dataset <- get(paste0("event_dataset_", iae)); rm(list=paste0("event_dataset_", iae))
      # covid_dataset
      load( file=paste0(dirtemp_cohort, "covid_dataset", suffix[[subpop]], ".RData"))
      
      for(isex in unique(data_vax$sex))
        for(iage in 1:length(age_groups_start)){
          cat(paste("\n\n",isex," age: (",age_groups_start[iage],";",age_groups_stop[iage],"]:\n\n"))
          data_vax_ij <- data_vax[ data_vax$sex==isex & age_groups_start[iage] < data_vax$age_at_study_entry & data_vax$age_at_study_entry <= age_groups_stop[iage],]
          event_dataset_ij <- event_dataset[event_dataset$person_id %in% data_vax_ij$person_id,]
          covid_dataset_ij <- covid_dataset[covid_dataset$person_id %in% data_vax_ij$person_id,]
          match_pop_res <- match_pop(  data_vax_ij , 
                                       matching_indep_factors=match_indep_factors,
                                       cond_before_vax = list( list( cov_name="covid", method="last_before", period = 30, unit="days", unit_short="d", condition=">", dataset=covid_dataset_ij ) , 
                                                             list( cov_name=iae, method="last_before", period=365, unit="days", unit_short="d", condition=">", dataset=event_dataset_ij, create_cov_var=T ) 
                                       ),
                                       matching_dep_factors = list( c( date="covid_date_last_before", event="covid_last_before", abs_diff=30   , units="days" ),
                                                                    c( date="birth_date"                                       , abs_diff=2*366, units="days" )     ),
                                       ids_study = ids_study_all_matching,
                                       vax_time="vax_days", vax1_time="vax1_days", vax_date="vax_date", #event_time=paste0(iae,"_days"), 
                                       vax_vars=c("vax_brand","vax_n"), 
                                       start_interval="start_days", stop_interval="stop_days", death_time="death_days",
                                       result_name = paste0(match_name,"_",isex,"_",round(age_groups_start[iage]),"_",age_groups_stop[iage]),
                                       file_name    = paste0(match_name,"_",isex,"_",round(age_groups_start[iage]),"_",age_groups_stop[iage],".RData"), 
                                       dir=dirtemp_event, lprint=F, lprint_str=T, lprint_warn=T, lparal=F, nboot=0) # nboot=200)
          
          gc()
        }
      
      save(match_pop_res,file=paste0(dirtemp_event,"match_pop_res.RData"))
      
    
      
      print(Sys.time())
      if("match_info" %in% ls()) rm(match_info)
      for(isex in unique(data_vax$sex))
        for(iage in 2:length(age_groups)){
          cat(paste0(iae,"_",isex,"_",round(age_groups[iage-1]),"_",age_groups[iage],": "))
          load(file=paste0(dirtemp_event,match_name,"_",isex,"_",round(age_groups[iage-1]),"_",age_groups[iage], suffix[[subpop]], ".RData"))
          tmp <- get(paste0(match_name,"_",isex,"_",round(age_groups[iage-1]),"_",age_groups[iage]))
          if("match_info" %in% ls()) { 
            match_info$matched_cases <- c( match_info$matched_cases, tmp$matched_cases )
            match_info$matched_controls$boot_0 <- c( match_info$matched_controls$boot_0, tmp$matched_controls$boot_0 )
            match_info$less_controls_matched_cases <- rbind.data.frame( match_info$less_controls_matched_cases, tmp$less_controls_matched_cases )
            match_info$less_controls_unmatched_cases <- rbind.data.frame( match_info$less_controls_unmatched_cases, tmp$less_controls_unmatched_cases )
            match_info$matching_factor_levels<- c( match_info$matching_factor_levels, tmp$matching_factor_levels )
          }
          else {
            match_info <- tmp
            match_info$parameters$result_name <- match_name
            match_info$parameters$file_name   <- paste0(match_name, suffix[[subpop]], ".RData")
          }
          rm(list=paste0(match_name,"_",isex,"_",round(age_groups[iage-1]),"_",age_groups[iage]))
        }
      assign(match_info$parameters$result_name, match_info)
      save(list=match_info$parameters$result_name,file=paste0(match_info$parameters$dir,match_info$parameters$file_name))
      print(Sys.time())
      #
      ############################
      
      
      
      
      
      
      
      
        #    match_info_sex_age_2comorb
        load(file=paste0(dirtemp_event,"match_info_sex_age_2comorb_",iae, suffix[[subpop]], ".RData"))
        matched_data <- get_matched_dataset(get(paste0("match_info_sex_age_2comorb_", iae)), data_vax, next_vax_time="next_vax_days", iboot=0)
        
        
        
        
        
        #    match_info_sex_age
        load(file=paste0(dirtemp_event,"match_info_sex_age_",iae, suffix[[subpop]], ".RData"))
        matched_data <- get_matched_dataset(get(paste0("match_info_sex_age_", iae)), data_vax, next_vax_time="next_vax_days", iboot=0)
        
        
         
            
      # match vaxed and unvaxed only on calendar day:  (about 15 minutes in CPRD)  
      match_pop_res <- match_pop(  data_vax, 
                                   vax_time="vax_days", vax1_time="vax1_days", event_time=paste0(iae,"_days"), vax_vars=c("vax_brand","vax_n"), 
                                   start_interval="start_days", stop_interval="stop_days", death_time="death_days",
                                   result_name = paste0("match_info_",iae), file_name = paste0("match_info_",iae,".RData"), dir=dirtemp_event, 
                                   lprint=T, lparal=F, nboot=0) # nboot=200)

      
      #  match only on datum:
      load(file=paste0(dirtemp_event,"match_info_",iae, suffix[[subpop]], ".RData"))
      matched_data <- get_matched_dataset(get(paste0("match_info_", iae)), data_vax, next_vax_time="next_vax_days", iboot=0)

     
      

      
      
      
      
      
      
      
      
            
      
      
      
      
      #####################################
      #
      #   PS -approach:
      #
      
      
      
   
      
      mdata <- matched_data[matched_data$vax_n_case==1 & tolower(matched_data$vax_brand_case)=="pfizer",]
      
      dd <- mdata
      dd <- sampled_data
      
      names(dd) <- gsub("_last_before","",names(dd) )  
      

      print(summary(log_regr_res2 <- glm( case_control ~ sex + age_at_study_entry + cancer + ckd + cvd + diab + hiv + 
                                            immunosup + obes + pulmon + sickle + as.numeric(number_conditions) , 
                                          data=dd, family="binomial")))
      print(summary(log_regr_res3 <- glm( case_control ~ sex + age_at_study_entry + I(age_at_study_entry^2) + I(age_at_study_entry^3) + cancer + ckd + cvd + diab + hiv + 
                                            immunosup + obes + pulmon + sickle + as.numeric(number_conditions) , 
                                          data=dd, family="binomial")))
      
      dd$ps2 <- predict(log_regr_res2, type="response")
      dd$swt2 <- NA
      dd$swt2[dd$case_control==1] <- ( (sum(dd$case_control==1)/nrow(dd)) / dd$ps2 )[dd$case_control==1] 
      dd$swt2[dd$case_control==0] <- ( (sum(dd$case_control==0)/nrow(dd)) / (1-dd$ps2) )[dd$case_control==0] 
      
      
      
      ps_contr <- hist(dd$ps2[dd$case_control==0], 100,plot=F)
      ps_case  <- hist(dd$ps2[dd$case_control==1], 100,plot=F)
      ylimm <- c(0,max(ps_contr$counts,ps_case$counts))
      plot(ps_contr, col="blue", xlim=range(dd$ps2), ylim=ylimm); par(new=T)
      plot(ps_case, col=rgb(1,0,0,alpha=0.8),xlim=range(dd$ps2),ylim=ylimm)
      
      summary(dd$ps2[dd$case_control==0])
      summary(dd$ps2[dd$case_control==1])
      
      
      swt_contr <- hist(dd$swt2[dd$case_control==0], 100,plot=F)
      swt_case  <- hist(dd$swt2[dd$case_control==1], 100,plot=F)
      ylimm <- c(0,max(swt_contr$counts,swt_case$counts))
      plot(swt_contr, col="blue", xlim=range(dd$swt2), ylim=ylimm); par(new=T)
      plot(swt_case, col=rgb(1,0,0,alpha=0.6),xlim=range(dd$swt2),ylim=ylimm)
      
      summary(dd$swt2[dd$case_control==0])
      summary(dd$swt2[dd$case_control==1])
      
      
      
      
      
      ##############################
      #  
      #   Cox regression:
      
      dd$survival_time <- pmin( dd[,paste0(iae,"_days")], dd$end_same_vax_days, dd$death_days, dd$study_exit_days, na.rm=T )
      dd$event_after_vax_days <- dd$survival_time - dd$index_days
      
      summary(dd$survival_time)
      summary(pmin( dd[,paste0(iae,"_days")],          dd$death_days, dd$study_exit_days, na.rm=T ))
      dd[,iae] <- as.numeric( dd$survival_time == dd[,paste0(iae,"_days")]  )
      dd[is.na(dd[,iae]),iae] <- 0
      
      summary(cox_res     <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control, data=dd))
      summary(cox_res_wt  <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control, data=dd, weight=swt2))
      summary(cox_res_cov <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control + sex + age_at_study_entry + cancer + ckd + cvd + diab + hiv + 
                                    immunosup + obes + pulmon + sickle + as.numeric(number_conditions) , data=dd, weight=swt2))
      
     
      
       
      
      
      #######################################
      #
      #         Survival
      # 
      #    only baseline information used!!! 
      #
      
      
      #matched_data <- match_pop_all
    
      ylimm <- c(0.99975,1)
      ylimm <- c(0.998,1)
      
      
      survfit_all <- survfit( Surv( event_after_vax_days, get(iae) ) ~ case_control, data=dd)

      plot(survfit_all, ylim=ylimm,col=c("blue","red"),lwd=2, main= paste0(iae," unweighted"));     abline(v=0,col="gray",lty=33);grid()
      legend("bottomleft",legend=names(survfit_all$strata),col=c("blue","red"),lty=1,lwd=3,bty="n")
      
      # weighted:
      survfit_all_wt <- survfit( Surv( event_after_vax_days, get(iae) ) ~ case_control, data=dd, weight=swt2)

      plot(survfit_all_wt, ylim=ylimm,col=c("blue","red"),lwd=2, main=paste0(iae," weighted"));       abline(v=0,col="gray",lty=33);grid()
      legend("bottomleft",legend=names(survfit_all$strata),col=c("blue","red"),lty=1,lwd=3,bty="n")
      
      
      
      #par(mfrow=c(2,2))
      for(iv in 1:max(dd$vax_n_case)){
        survfit_all <- survfit( Surv( event_after_vax_days, get(iae) ) ~ vax_n_case + vax_brand_case + case_control, data=dd[dd$vax_n_case==iv,])
        coll <- rainbow(length(survfit_all$strata))
        #coll <- colors()[sample.int(length(colors()),length(survfit_all$strata))]
        plot(survfit_all, ylim=ylimm, col=coll,lwd=2, lty=c(2,1), main=paste0(iae,"; vax n=",iv));      abline(v=0,col="gray",lty=33);grid()
        legend("bottomleft",legend=names(survfit_all$strata),col=coll,lty=c(3,1),lwd=4, box.lty=0, cex=0.8)
      }
      
      par(mfrow=c(2,2))
      for(iv in 1:max(dd$vax_n_case)){
        survfit_all <- survfit( Surv( event_after_vax_days, get(paste0(iae,"_before_next_vax")) ) ~ vax_n_case + vax_brand_case + case_control, data=dd[dd$vax_n_case==iv,])
        coll <- rainbow(length(survfit_all$strata))
        #coll <- colors()[sample.int(length(colors()),length(survfit_all$strata))]
        plot(survfit_all, ylim=ylimm, col=coll,lwd=2, lty=c(2,1), main=paste0(iae,"; vax n=",iv, " before ", iv+1));      abline(v=0,col="gray",lty=33);grid()
        legend("bottomleft",legend=names(survfit_all$strata),col=coll,lty=c(3,1),lwd=4, box.lty=0, cex=0.8)
      }
     
 
  } # end iae      
  ##########
  # restore options:
  options(old_width)
  
  
}  # end of 'subpop' loop

