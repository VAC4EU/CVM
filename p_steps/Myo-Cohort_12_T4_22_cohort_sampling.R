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
  
  # Load dataset 'data_vax_for_sampling' by loading file 'data_vax_for_sampling.RData'
  load(file=paste0(dirtemp_cohort, "data_vax_for_sampling", suffix[[subpop]], ".RData"))
  load(file=paste0(dirtemp_cohort,"ids_study_all_sampling", suffix[[subpop]], ".RData"))
  
  
  
  
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
      
      
  
      ###############################################################################################################################################      
      ###############################################################################################################################################      
      ###############################################################################################################################################      
      ###############################################################################################################################################      
      ###############################################################################################################################################      
      ###############################################################################################################################################      
      ###############################################################################################################################################
      #
      #
      #
      #          design B: random sampling from before vax:
      #
      #
      
      data_vax <- data_vax_for_sampling

      
      sampled_data <- sample_controls(data_vax_for_sampling, event_date = paste0(iae,"_date_first_after"), ids_study = ids_study_all_sampling,
                                      cond_before_vax=list( list( cov_name="covid", method="last_before", period= 30, unit="days", unit_short="d", condition=">", 
                                                                  name="D3_TD_variable_COVID", dir=paste0(dirtemp,"td/") ),
                                                            list( cov_name="myocarditis",method="last_before", period=365, unit="days", unit_short="d", condition=">", 
                                                                  name="D3_events_C_MYOCARD_AESI_simple", dir=paste0(dirtemp,"events/"), create_cov_var=T ) ) )
      
       
      
      # 4.  add covariates:
      #
      diseases <- c("CANCER","CKD","CVD","DIAB","HIV","IMMUNOSUP","OBES","PULMON","SICKLE","NUMBER_CONDITIONS")
      for(iziekte in diseases)
        sampled_data <- add_td_covariate( sampled_data, iziekte , methods = "last_before", na_value=0, 
                                          name=paste0("D3_TD_variable_", iziekte, suffix[[subpop]]), ids_study = ids_study_all_sampling)     

       assign(paste0("sampled_data_",iae),sampled_data)
       save(list=paste0("sampled_data_",iae), file=paste0(dirtemp_event,"sampled_data_",iae,suffix[[subpop]], ".RData"))
       
 
      
      #####################################
      #
      #   PS -approach:
      #
      
      dd <- sampled_data
      
      names(dd) <- gsub("_last_before","",names(dd) )  
      
      print(summary(log_regr_res <- glm( case_control ~ sex + age_at_study_entry + cancer + ckd + cvd + diab + hiv + 
                                            immunosup + obes + pulmon + sickle + as.numeric(number_conditions) , 
                                          data=dd, family="binomial")))
      print(summary(log_regr_res <- glm( case_control ~ sex + age_at_study_entry + cancer + ckd + cvd + diab + hiv + 
                                            immunosup + obes + pulmon + sickle  , 
                                          data=dd, family="binomial")))
      
      dd$ps <- predict(log_regr_res, type="response")
      dd$swt <- NA
      dd$swt[dd$case_control==1] <- ( (sum(dd$case_control==1)/nrow(dd)) / dd$ps )[dd$case_control==1] 
      dd$swt[dd$case_control==0] <- ( (sum(dd$case_control==0)/nrow(dd)) / (1-dd$ps) )[dd$case_control==0] 
      
      
      dd$wt <- NA
      dd$wt[dd$case_control==1] <- ( 1 / dd$ps )[dd$case_control==1] 
      dd$wt[dd$case_control==0] <- ( 1/ (1-dd$ps) )[dd$case_control==0] 
      
      
      
      
      ps_contr <- hist(dd$ps[dd$case_control==0], 100,plot=F)
      ps_case  <- hist(dd$ps[dd$case_control==1], 100,plot=F)
      ylimm <- c(0,max(ps_contr$counts,ps_case$counts))
      plot(ps_contr, col=rgb(t(col2rgb("green"))/255,alpha=0.5), xlim=range(dd$ps), ylim=ylimm); par(new=T)
      plot(ps_case, col=rgb(t(col2rgb("red"))/255,alpha=0.5),xlim=range(dd$ps),ylim=ylimm)
      
      summary(dd$ps[dd$case_control==0])
      summary(dd$ps[dd$case_control==1])
      
      
      swt_contr <- hist(dd$swt[dd$case_control==0], 100,plot=F)
      swt_case  <- hist(dd$swt[dd$case_control==1], 100,plot=F)
      ylimm <- c(0,max(swt_contr$counts,swt_case$counts))
      plot(swt_contr, col="skyblue", xlim=range(dd$swt), ylim=ylimm); par(new=T)
      plot(swt_case, col=rgb(t(col2rgb("red"))/255,alpha=0.4),xlim=range(dd$swt),ylim=ylimm)
      
      summary(dd$swt[dd$case_control==0])
      summary(dd$swt[dd$case_control==1])
      
      
      
      ############# balance per covariate:
        
      par(mfrow=c(3,4))
      
      
      cov_names <- c("sex","age_at_study_entry")
      tab_3a <- matrix(NA,nrow=length(cov_names),ncol=9*6+4*3, 
                       dimnames=list(cov_names,
                                     c( paste0(rep(c("vax_","sampled_","vax_swt_","sampled_swt_","vax_wt_","sampled_wt_"),each=9), 
                                            rep( c("n","n_nonmissing","mean","sd","min","Q1","median","Q3","max"),  6)),
                                        paste0(rep(c("diff_","diff_swt_","diff_wt_"),each=4), rep(c("value","CI_left","CI_right","p_value"),3) )
                                     )  ) )
      for(icov in cov_names)  tab_3a <- cov_balance(icov, dd, categorical=F, tab=tab_3a)
      
       
      
      #
      #############################
      
      
      
      
      
      
      
      
      
      
      
      ##############################
      #  
      #   Cox regression:
      dd[,paste0(iae,"_days_first_after")] <- as.numeric(difftime(dd[,paste0(iae,"_date_first_after")],as.Date("2020-08-31"),units="days"))
      dd$index_days <- as.numeric(difftime(dd$index_date,as.Date("2020-08-31"),units="days"))
      dd$survival_time <- pmin( dd[,paste0(iae,"_days_first_after")], dd$end_same_vax_days, dd$death_days, dd$study_exit_days, na.rm=T )
      #dd$survival_time <- pmin( dd[,paste0(iae,"_days")], dd$end_same_vax_days, dd$death_days, dd$study_exit_days, na.rm=T )
      dd$event_after_vax_days <- dd$survival_time - dd$index_days
      
      summary(dd$survival_time)
      summary(pmin( dd[,paste0(iae,"_days_first_after")],          dd$death_days, dd$study_exit_days, na.rm=T ))
      dd[,iae] <- as.numeric( dd$survival_time == dd[,paste0(iae,"_days_first_after")]  )
      dd[is.na(dd[,iae]),iae] <- 0
      
      summary(cox_res     <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control, data=dd))
      summary(cox_res_wt  <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control, data=dd, weight=swt))
      summary(cox_res_cov <- coxph(Surv( event_after_vax_days, get(iae)) ~ case_control + sex + age_at_study_entry + cancer + ckd + cvd + diab + hiv + 
                                     immunosup + obes + pulmon + sickle + as.numeric(number_conditions) , data=dd, weight=swt))
      
    
      
      
      
      
      #######################################
      #
      #         Survival
      # 
      #    only baseline information used!!! 
      #
      
       
      ylimm <- c(0.99975,1)
      ylimm <- c(0.998,1)
      
      
      survfit_all <- survfit( Surv( event_after_vax_days, get(iae) ) ~ case_control, data=dd)
       
      plot(survfit_all, ylim=ylimm,col=c("blue","red"),lwd=2, main= paste0(iae," unweighted"));     abline(v=0,col="gray",lty=33);grid()
      legend("bottomleft",legend=names(survfit_all$strata),col=c("blue","red"),lty=1,lwd=3,bty="n")
      
      # weighted:
      survfit_all_wt <- survfit( Surv( event_after_vax_days, get(iae) ) ~ case_control, data=dd, weight=swt)
      
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
      
      
      
      
      
      
      
      survfit_all <- survfit( Surv( event_after_vax_days_matched_min, get(paste0(iae,"_matched_min")) ) ~ case_control, data=dd)
      
      plot(survfit_all, ylim=ylimm,col=c("blue","red"),lwd=2, main=paste0(iae,"_matched_min"));      abline(v=0,col="gray",lty=33);grid()
      legend("bottomleft",legend=names(survfit_all$strata),col=c("blue","red"),lty=1,lwd=3)
      
      
      par(mfrow=c(2,2))
      for(iv in 1:max(dd$vax_n_case)){
        survfit_all <- survfit( Surv( event_after_vax_days_matched_min, get(paste0(iae,"_matched_min")) ) ~ vax_n_case + vax_brand_case + case_control, data=dd[dd$vax_n_case==iv,])
        coll <- rainbow(length(survfit_all$strata))
        #coll <- colors()[sample.int(length(colors()),length(survfit_all$strata))]
        plot(survfit_all, ylim=ylimm, col=coll,lwd=2, lty=c(2,1), main=paste0(iae,"_matched_min","; vax n=",iv));      abline(v=0,col="gray",lty=33);grid()
        legend("bottomleft",legend=names(survfit_all$strata),col=coll,lty=c(3,1),lwd=4, box.lty=0, cex=0.8)
      }
      
      par(mfrow=c(2,2))
      for(iv in 1:max(dd$vax_n_case)){
        survfit_all <- survfit( Surv( event_or_next_after_vax_days_matched_min, get(paste0(iae,"_before_next_vax_matched_min")) ) ~ vax_n_case + vax_brand_case + case_control, data=dd[dd$vax_n_case==iv,])
        coll <- rainbow(length(survfit_all$strata))
        #coll <- colors()[sample.int(length(colors()),length(survfit_all$strata))]
        plot(survfit_all, ylim=ylimm, col=coll,lwd=2, lty=c(2,1), main=paste0(iae,"_matched_min; vax n=",iv, " before ", iv+1));      abline(v=0,col="gray",lty=33);grid()
        legend("bottomleft",legend=names(survfit_all$strata),col=coll,lty=c(3,1),lwd=4, box.lty=0, cex=0.8)
      }
      
    
 
  } # end iae      
  ##########
  # restore options:
  options(old_width)
  
  
}  # end of 'subpop' loop

