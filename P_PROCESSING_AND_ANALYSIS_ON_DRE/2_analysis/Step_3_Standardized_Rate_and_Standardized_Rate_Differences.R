# Standardized Rate and Rate differences


# ---------------------------------------------------Information-------------------------------------------------------

#       Title         : Calculating Standardized Rate and Standardized Rate differences 
#       Author        : Ivonne Martin
#       Date created  : 28 Jan 2022
#       Input         : D4_persontime_monthly_poisson_RF.RData
#       Output        : standardized Rates and Standardized Rate differences tables for each DAP. 
#     
#       external Function used : 
#                       - dsr.R      # to calculate standardized rates
#                       - dsrr.R     # to calculate standardized rates difference,
#                                  # both dsr.R and dsrr.R were developed in github.com/mattkumar/dsr 
#
# ---------------------------------------------------------------------------------------------------------------------

rm(list=ls())
library(dplyr)
library(tidyr)
library(data.table)

load("Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/D4_persontime_monthly_poisson_RF.RData")
FDat1 <- as.data.frame(D4_persontime_monthly_poisson_RF)


## Removing vaccine == "UKN"

FDat <- FDat1[-which(FDat1$Vaccine1 == "UKN" | FDat1$Vaccine2 == "UKN"),]

## First remove colnames with "broad"

FD <- FDat[,-c(grep("broad",colnames(FDat)))]

idxPT <- grep("Persontime", colnames(FDat))
CN <- colnames(FD)
idxPT <- grep("Persontime",CN)
idxAESI <- c(74:113, 116:117, 122:124) ## the same as previous. So, fine


## AgeCategory, optional 

AgeCatfun <- function(par){
  if(par == "0-4" | par == "5-11" | par == "12-17"){
    agecat <- "<18"
  } else if (par == "18-24" | par == "25-29" | par == "30-39" | par == "40-49" | par == "50-59") {
    agecat <- "18-59"
  } else {agecat <- "60+"}
  return(agecat)  
}

## Vaccine Dose Category, used to categorized rowdata into background or vaccinated (either dose). 

VacDoseCat <- function(data){
  vax.cat <- vax.cat1 <- NULL
  
  if(data$year == "2020" & data$Dose1 == 0 & data$Dose2 == 0){
    vax.cat <- 0
  }
}


# Prepare the pooled data across DAP. Here the function is copied from Step1_Preprocess for Poisson Analysis

PrepVac <- function(Dat){
  p.Dat <- Dat
  
  p.Dat$AgeCat <-  relevel(as.factor(as.character(sapply(p.Dat$ageband_at_study_entry,AgeCatfun))),ref = "18-59")
  
  
  p.Dat$Gender <- as.factor(p.Dat$Gender)
  p.Dat$month <- as.factor(p.Dat$month)
  p.Dat$year <- as.factor(p.Dat$year)
  p.Dat$Dose1 <- as.factor(p.Dat$Dose1)
  p.Dat$Dose2 <- as.factor(p.Dat$Dose2)
  p.Dat$anyrisk.cat <- as.factor(p.Dat$any_risk_factors_at_date_vax)
  
  p.Dat$risk.fac.new <- as.factor(p.Dat$any_risk_factors_at_date_vax)
  p.Dat$risk.fac.new[which(p.Dat$year == "2020")] <- p.Dat$any_risk_factors_at_study_entry[which(p.Dat$year == "2020")]
  
  p.Dat$anyrisk.cat.pre <- as.factor(p.Dat$any_risk_factors_at_study_entry)
  
  # COVID19.1 
  p.Dat$COVID19.1 <- ifelse(p.Dat$COVID19 == 1,0,1)
  p.Dat$COVID19.1 <- as.factor(p.Dat$COVID19.1)
  
  p.Dat$Vaccine1 <- relevel(as.factor(as.character(p.Dat$Vaccine1)),ref = "none")
  p.Dat$Vaccine2 <- relevel(as.factor(as.character(p.Dat$Vaccine2)),ref = "none")
  
  idxrem1 <- which(p.Dat$year == "2021" & p.Dat$Dose1 == 0 & p.Dat$Dose2 == 0)
  if (length(idxrem1 != 0)){
    p.Dat.vac <- p.Dat[-which(p.Dat$year == "2021" & p.Dat$Dose1 == 0 & p.Dat$Dose2 == 0),]
  } else {p.Dat.vac = p.Dat}
  
  idxvax12 <- which(p.Dat.vac$year == "2021" | (p.Dat.vac$year == "2020" & p.Dat.vac$Dose1 == 1 & p.Dat.vac$Dose2 == 0) | (p.Dat.vac$year == "2020" & p.Dat.vac$Dose2 == 1))
  idxvax1 <- which(p.Dat.vac$Dose1 == 1 & p.Dat.vac$Dose2 == 0)
  idxvax2 <- which(p.Dat.vac$Dose2 == 1)
  
  p.Dat.vac$vax.cat1 <- p.Dat.vac$vax.cat <- rep(0, nrow(p.Dat.vac))
  p.Dat.vac$vax.cat[idxvax12] <- 1
  p.Dat.vac$vax.cat1[idxvax1] <- 1
  p.Dat.vac$vax.cat1[idxvax2] <- 2
  
  return(p.Dat.vac)
}

# FD1 now contain AESI for unvaccinated and vaccinated in 2021. All rowdata in 2021 who did not get vaccination
# has been removed. 

FD1 <- PrepVac(FD)

## Reference population. 
pop.eustat<- read.csv("Z:/folder_analysis_for_rebuttal/3_editing/g_csv/ESP_ageband.csv",sep = "")


Dat.ARS <- subset(FD1, DAP == "ARS_ALL")
Dat.BPC <- subset(FD1, DAP == "BIFAP_PC")
Dat.BPCH <- subset(FD1, DAP == "BIFAP_PC_HOSP")
Dat.CPRD <- subset(FD1, DAP == "CPRD_ALL")
Dat.PH <- subset(FD1, DAP == "PHARMO_ALL")



#------ Procedure for calculating rate. For specific AESI --------
# 1. For each DAP, the loop runs for each AESI to calculate the direct standardized rate for each vaccine brand, and dose
#    then calculate their standardized rate difference with the background. 
# 2. The loop will first check whether DAP could identify the AESI. If this is not the case, the standardized rate and 
#    the standardized rate difference will be presented as missing. 
#-----------------------------------------------------------------


dsr <- function(data, event, fu, subgroup, ..., refdata, mp, method="normal", sig=0.95, decimals ) {
  
  subgroup <- enquo(subgroup)
  event <- enquo(event)
  fu <- enquo(fu)
  stdvars <- quos(...)
  
  
  
  #Compute crude and standardized rates and variances
  all_data_st = data %>%
    left_join(refdata) %>%
    group_by(!!subgroup) %>%
    mutate(n=sum(!!event),
           d=sum(!!fu),
           cr_rate=n/d,
           cr_var=n/d^2,
           wts=pop/sum(pop),
           st_rate=sum(wts*(!!event/!!fu)),
           st_var=sum(as.numeric((wts^2)*(!!event/(!!fu )^2)))) %>%
    distinct(!!subgroup, .keep_all=TRUE) %>%
    select(!!subgroup, n, d, cr_rate, cr_var, st_rate, st_var)
  
  
  
  #Compute Confidence Intervals (CI) according to method. The default is 'normal'
  if(method=="gamma") {
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*qgamma((1-sig)/2, shape=cr_rate^2/(cr_var))/(cr_rate/cr_var),
        c_upper=mp*qgamma(1-((1-sig)/2), shape=1+cr_rate^2/(cr_var))/(cr_rate/cr_var),
        s_rate=mp*st_rate,
        s_lower=mp*qgamma((1-sig)/2, shape=st_rate^2/st_var)/(st_rate/st_var),
        s_upper=mp*qgamma(1-((1-sig)/2), shape=1+(st_rate^2/st_var))/(st_rate/st_var)) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
    
  } else if (method=="normal") {
    
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*(cr_rate+qnorm((1-sig)/2)*sqrt(cr_var)),
        c_upper=mp*(cr_rate-qnorm((1-sig)/2)*sqrt(cr_var)),
        s_rate=mp*st_rate,
        s_lower=mp*(st_rate+qnorm((1-sig)/2)*sqrt(st_var)),
        s_upper=mp*(st_rate-qnorm((1-sig)/2)*sqrt(st_var))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
  } else if (method=="lognormal") {
    
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*exp((log(cr_rate)+qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        c_upper=mp*exp((log(cr_rate)-qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        s_rate=mp*st_rate,
        s_lower=mp*exp((log(st_rate)+qnorm((1-sig)/2)*sqrt(st_var)/(st_rate))),
        s_upper=mp*exp((log(st_rate)-qnorm((1-sig)/2)*sqrt(st_var)/(st_rate)))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
  }
  
  
  #Clean up and output
  
  tmp1$c_rate  <- round(tmp1$c_rate,  digits=decimals)
  tmp1$c_lower <- round(tmp1$c_lower, digits=decimals)
  tmp1$c_upper <- round(tmp1$c_upper, digits=decimals)
  tmp1$s_rate  <- round(tmp1$s_rate, digits=decimals)
  tmp1$s_lower <- round(tmp1$s_lower, digits=decimals)
  tmp1$s_upper <- round(tmp1$s_upper, digits=decimals)
  
  colnames(tmp1) <-  c('Subgroup', 'Numerator','Denominator',
                       paste('Crude Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Crude)',sep=''),
                       paste(sig*100,'% UCL (Crude)',sep=''),
                       paste('Std Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Std)',sep=''),
                       paste(sig*100,'% UCL (Std)',sep=''))
  
  
  tmp1 <- as.data.frame(tmp1)
  
}


dsrr <- function(data, event, fu, subgroup, ..., refdata, estimate, refgroup, mp, sig = 0.95, decimals){
  
  refgroup <- enquo(refgroup)
  subgroup <- enquo(subgroup)
  event <- enquo(event)
  fu <- enquo(fu)
  stdvars <- quos(...)
  
  
  ## Compute crude and standardized rates and variances
  
  all_data_st = data %>%
    left_join(refdata) %>%
    group_by(!!subgroup) %>%
    mutate(n=sum(!!event),
           d=sum(!!fu),
           cr_rate=n/d,
           cr_var=n/d^2,
           wts=pop/sum(pop),
           st_rate=sum(wts*(!!event/!!fu)),
           st_var=sum(as.numeric((wts^2)*(!!event/(!!fu)^2)))) %>%
    distinct(!!subgroup, .keep_all=TRUE) %>%
    select(!!subgroup, n, d, st_rate, st_var)
  
  ## Select a referent
  
  ref = all_data_st %>% 
    filter(!! subgroup == !! quo_name(refgroup)) %>% 
    select(reference=!!subgroup, st_rate, st_var)
  
  ## Compute comparisons based on desired estimate
  
  if(estimate == "ratio"){
    out = all_data_st %>%
      mutate(ratio=(st_rate/ref$st_rate),
             se = sqrt(((1/(st_rate^2)) * st_var) + 
                         ((1/(ref$st_rate^2))*ref$st_var)),
             lower = exp(log(ratio + qnorm((1-sig)/2)*se)),
             upper = exp(log(ratio - qnorm((1-sig)/2)*se)),
             reference = ref$reference,
             s_rate = mp*st_rate) %>%
      select(!!subgroup, reference, s_rate, ratio, lower, upper)
    
    ## Clean up and output 
    
    out$s_rate <- round(out$s_rate, digits = decimals)
    out$lower <- round(out$lower, digits = decimals)
    out$upper <- round(out$upper, digits = decimals)
    out$ratio <- round(out$ratio, digits = decimals)
    
    colnames(out) <- c('Comparator', 'Reference',
                       paste('Std Rate (per ',mp,')', sep = ' '),
                       paste('Rate Ratio (RR)', sep=' '),
                       paste(sig*100,'% LCL(RR)', sep=' '),
                       paste(sig*100,'% UCL(RR)',sep= ' ')
    )
    
    out <- as.data.frame(out)
    
  } else if(estimate == "difference"){
    
    out = all_data_st %>% 
      mutate(ratio=(st_rate-ref$st_rate),
             se=sqrt((st_var) + (ref$st_var)),
             lower = mp*(ratio + qnorm((1-sig)/2)*se),
             upper = mp*(ratio - qnorm((1-sig)/2)*se),
             reference = ref$reference,
             s_rate=mp*st_rate,
             difference=mp*ratio) %>%
      select(!!subgroup, reference, s_rate, difference, lower, upper)
    
    ## Clean up and output
    
    out$s_rate <- round(out$s_rate, digits = decimals)
    out$lower <- round(out$lower, digits = decimals)
    out$upper <- round(out$upper, digits = decimals)
    out$difference <- round(out$difference, digits = decimals)
    
    colnames(out) <- c('Comparator','Reference',
                       paste('Std Rate (per ',mp,')', sep= ' '),
                       paste('Rate Difference (RD)', sep= ' '),
                       paste(sig*100,'% LCL (RD)',sep=' '),
                       paste(sig*100,'% UCL (RD)', sep = ' '))
    
    out <- as.data.frame(out)
  }
  
}



# Checking whether DAP has identified the AESI
check.dat <- function(data){
  if(sum(is.na(data)) == length(data)){
    val <- 0
  } else { val <- 1}
  return(val)}


format_df <- function(fml1, data){
  tapdat <- aggregate(fml1, data = data,function(x) sum(na.omit(x)))
  tapdat[,4] <- tapdat[,4]/365
  
  df_base <- data.frame(age = tapdat$ageband_at_study_entry[which(tapdat$Vaccine1 == "none")],
                        aesi = tapdat[which(tapdat$Vaccine1 == "none"),3],
                        fu=c(tapdat[which(tapdat$Vaccine1 == "none"),4]),
                        state = 'Baseline')
  df_AZ <- data.frame(age = tapdat$ageband_at_study_entry[which(tapdat$Vaccine1 == "AstraZeneca")],
                      aesi = tapdat[which(tapdat$Vaccine1 == "AstraZeneca"),3],
                      fu=tapdat[which(tapdat$Vaccine1 == "AstraZeneca"),4],
                      state = 'AZ')
  
  if(length(unique(tapdat$Vaccine1)) == 5){   # This is to cater for CPRD which does not have data on J&J
    df_JJ <- data.frame(age = tapdat$ageband_at_study_entry[which(tapdat$Vaccine1 == "J&J")],
                        aesi = tapdat[which(tapdat$Vaccine1 == "J&J"),3],
                        fu=tapdat[which(tapdat$Vaccine1 == "J&J"),4],
                        state = 'JJ')} else{
                          data.frame(age = rep(c('0-4','12-17','18-24','25-29','30-39','40-49',
                                                 '5-11','50-59','60-69','70-79','80+'),1),
                                     aesi = rep(NA,11),
                                     fu = rep(NA,11),
                                     state = "JJ")
                        }
  df_Moderna <- data.frame(age = tapdat$ageband_at_study_entry[which(tapdat$Vaccine1 == "Moderna")],
                           aesi =  tapdat[which(tapdat$Vaccine1 == "Moderna"),3],
                           fu= tapdat[which(tapdat$Vaccine1 == "Moderna"),4],
                           state = 'Moderna')
  
  ## Pfizer has missing data for certain age category at the second dose. This is to cater that situation
  
  df_PF1 <- data.frame(age =  tapdat$ageband_at_study_entry[which(tapdat$Vaccine1 == "Pfizer")],
                       aesi =  tapdat[which(tapdat$Vaccine1 == "Pfizer"),3],
                       fu = tapdat[which(tapdat$Vaccine1 == "Pfizer"),4],
                       state = 'Pfizer')
  idxrow <- which(is.na(df_PF1[,2]))
  
  if(length(idxrow) != 0 & length(idxrow) != nrow(df_PF1)){
    df_PF <- df_PF1[-idxrow,]
  } else {df_PF <- df_PF1}
  
  df_all <- rbind(df_base,df_AZ,df_JJ,df_Moderna,df_PF) 
  
  return(df_all)
}


# Calculating the standardized rates and standardized rate difference 

RateDiff <- function(temp,temp3,dat){
  
  cova1 <- c("ageband_at_study_entry","Vaccine1")
  fml1 <- as.formula(paste0("cbind(",temp,",",temp3,")","~",paste0(cova1, collapse = "+")))
  
  test.val <- check.dat(dat[,temp])
  if( test.val != 0){
    df_all12 <- format_df(fml1,dat)
    df_dose1 <- format_df(fml1, dat[which(dat$vax.cat1 != 2),])
    df_dose2 <- format_df(fml1,dat[which(dat$vax.cat1 != 1),])
    
    #population from eurostat used as reference population
    df_pop <- data.frame(age = rep(c('0-4','12-17','18-24','25-29','30-39','40-49',
                                     '5-11','50-59','60-69','70-79','80+'),1),pop = pop.eustat[,2])
    
    
    
    my_results12 <- dsrr(data=df_all12,
                         event=aesi,
                         fu=fu,
                         subgroup=state,
                         age, 
                         refdata=df_pop,
                         refgroup="Baseline",
                         estimate="difference",
                         sig=0.95,
                         mp=100000,decimals=4)
    
    stdrate12 <- dsr(data=df_all12,
                     event=aesi,
                     fu=fu,
                     subgroup=state,
                     age, 
                     refdata=df_pop,
                     method = "gamma",
                     refgroup="Baseline",
                     estimate="difference",
                     sig=0.95,
                     mp=100000,decimals=4)
    
    
    my_results1 <- dsrr(data=df_dose1,
                        event=aesi,
                        fu=fu,
                        subgroup=state,
                        age, 
                        refdata=df_pop,
                        method = "gamma",
                        refgroup="Baseline",
                        estimate="difference",
                        sig=0.95,
                        mp=100000,decimals=4)
    
    stdrate1 <- dsr(data=df_dose1,
                    event=aesi,
                    fu=fu,
                    subgroup=state,
                    age, 
                    refdata=df_pop,
                    method = "gamma",
                    refgroup="Baseline",
                    estimate="difference",
                    sig=0.95,
                    mp=100000,decimals=4)
    
    my_results2 <- dsrr(data=df_dose2,
                        event=aesi,
                        fu=fu,
                        subgroup=state,
                        age, 
                        refdata=df_pop,
                        refgroup="Baseline",
                        estimate="difference",
                        sig=0.95,
                        mp=100000,decimals=4)
    stdrate2 <- dsr(data=df_dose2,
                    event=aesi,
                    fu=fu,
                    subgroup=state,
                    age, 
                    refdata=df_pop,
                    method = "gamma",
                    refgroup="Baseline",
                    estimate="difference",
                    sig=0.95,
                    mp=100000,decimals=4)
    
    
    RD12 <- my_results12[2:5,4:6]        
    RD1 <- my_results1[2:5,4:6]
    RD2 <- my_results2[2:5,4:6]
    
    Std12 <- stdrate12[,c(1:3,7:9)]
    Std1 <- stdrate1[,c(1:3,7:9)]
    Std2 <- stdrate2[,c(1:3,7:9)]
    
    l.RD <- rbind(RD1,RD2,RD12)
    l.std <- rbind(Std1,Std2,Std12) }  else {
      l.RD <- matrix(NA,nrow = 12,ncol=3)
      colnames(l.RD) <- c("Rate Difference (RD)", "95 % LCL (RD)", "95 % UCL (RD)")
      l.std <- matrix(NA,nrow = 15,ncol = 6)
      colnames(l.std) <- c("Subgroup","Numerator","Denominator","Std Rate (per 1e+05)", "95% LCL (Std)",
                           "95% UCL (Std)")}
  return(list(l.RD,l.std))}




Calc.RD <- function(dataDap){
  list.RD <- NULL
  list.Std <- NULL
  
  for (i in 1:length(idxAESI)){
    temp <- CN[idxAESI[i]] 
    temp3 <- CN[idxPT[i+1]]
    
    output <- RateDiff(temp, temp3,dataDap)
    list.RD <- rbind(list.RD, output[[1]])
    list.Std <- rbind(list.Std, output[[2]])    
    print(i)  
  }
  
  
  
  
  list.RD.print <- cbind(rep(CN[idxAESI], each = 12),rep(c("AZ","JJ","Moderna","Pfizer"),135),
                         rep(c(rep("Dose1 - Back",4),rep("Dose2 - Back",4),rep("Dose12 - Back",4)),45),
                         list.RD)
  colnames(list.RD.print) <- c("AESI","Vaccine","Dose-Background", "RD","lower95%CI","upper95%CI")
  
  
  
  list.Std.print <- cbind(rep(CN[idxAESI],each = 15),rep(c(rep("Dose1",5), rep("Dose2",5), rep("Dose12",5)),45), list.Std)
  colnames(list.Std.print) <- c("AESI","Dose",colnames(list.Std.print)[3:8])
  list.Std.print[seq(1,675, by = 5),2] <- "" 
  
  res <- list(list.Std.print,list.RD.print)
  return(res)}



write.csv(list.RD.print, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/DAP_PHARMO_RateDiff.csv")
write.csv(list.Std.print, file ="Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/DAP_PHARMO_StandardizedRates.csv")
