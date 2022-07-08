# Preprocess the data for Poisson analysis 

# -----Information------------------------------------------
#
#   Author  : Ivonne Martin
#   Date    : 25-01-2022
#   Input   : D4_persontime_monthly_poisson_RF
#   Output  : 5 data tables for each DAP
#
#  
# ----------------------------------------------------------
#
#
# Details : 
# This preprocessing data is to align with the previous 
# analysis done in November - December 2021. Several steps in 
# the preprocessing :
# 1. Removing column names of each AESI that contains "broad", if any. 
# 2. Separate the data according to each DAP
# 3. Removing records with "UKN" vaccine. 
# 4. Define the function of merging age into 3 and 2 categories 
# 5. Define the PrepVac function to restructure the variables inside the datatable for analysis purpose 




rm(list=ls())

library(dplyr)
library(tidyr)
library(data.table)
library(MASS)

# loading the data 

load("Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/D4_persontime_monthly_poisson_RF.RData")
FDat <- as.data.frame(D4_persontime_monthly_poisson_RF)

## First remove colnames with "broad"

FD <- FDat[,-c(grep("broad",colnames(FDat)))]
idxPT <- grep("Persontime", colnames(FDat))
CN <- colnames(FD)
idxPT <- grep("Persontime",CN)
idxAESI <- c(74:113, 116:117, 122:124) ## the same as previous. So, fine

## Separate the data according to the DAP

Dat.ARS <- FD[which(FD$DAP == "ARS_ALL"),]
Dat.BIFAPPC <- FD[which(FD$DAP == "BIFAP_PC"),]
Dat.BIFAPPC_HOSP <- FD[which(FD$DAP == "BIFAP_PC_HOSP"),]
Dat.CPRD <- FD[which(FD$DAP == "CPRD_ALL"),]
Dat.PHARMO <- FD[which(FD$DAP == "PHARMO_ALL"), ]

## Removing UKN vaccine brand from Pharmo and BIFAPPC

Dat.PHARMO <- Dat.PHARMO[-which(Dat.PHARMO$Vaccine1 == "UKN" | Dat.PHARMO$Vaccine2 == "UKN"),]
Dat.BIFAPPC <- Dat.BIFAPPC[-which(Dat.BIFAPPC$Vaccine1 == "UKN" | Dat.BIFAPPC$Vaccine2 == "UKN"),]



## AgeCategory

AgeCatfun <- function(par){
  if(par == "0-4" | par == "5-11" | par == "12-17"){
    agecat <- "<18"
  } else if (par == "18-24" | par == "25-29" | par == "30-39" | par == "40-49" | par == "50-59") {
    agecat <- "18-59"
  } else {agecat <- "60+"}
  return(agecat)  
}

AgeCatfun2 <- function(par){
  if(par == "0-4" | par == "5-11" | par == "12-17" | par == "18-24" | par == "25-29" | par == "30-39" | par == "40-49" | par == "50-59"){
    agecat <- "<60"
  } else {agecat <- "60+"}
  return(agecat)  
}


# The following function PrepVac restructure several variables in each DAPs' datatable. 
# This requires releveling the factor variables, removal of unvaccinated record in 2021 as well
# as reformat the covid19 variable. 


PrepVac <- function(Dat){
  p.Dat <- Dat
  
  p.Dat$AgeCat <-  relevel(as.factor(as.character(sapply(p.Dat$ageband_at_study_entry,AgeCatfun))),ref = "18-59")
  p.Dat$AgeCat2 <- relevel(as.factor(as.character(sapply(p.Dat$ageband_at_study_entry,AgeCatfun2))), ref ="<60")
  
  p.Dat$Gender <- as.factor(p.Dat$Gender)
  p.Dat$month <- as.factor(p.Dat$month)
  p.Dat$year <- as.factor(p.Dat$year)
  p.Dat$Dose1 <- as.factor(p.Dat$Dose1)
  p.Dat$Dose2 <- as.factor(p.Dat$Dose2)
  p.Dat$anyrisk.cat <- as.factor(p.Dat$any_risk_factors_at_date_vax)
  
  p.Dat$risk.fac.new <- as.factor(p.Dat$any_risk_factors_at_date_vax)
  p.Dat$risk.fac.new[which(p.Dat$year == "2020")] <- p.Dat$any_risk_factors_at_study_entry[which(p.Dat$year == "2020")]
  
  p.Dat$anyrisk.cat.pre <- as.factor(p.Dat$any_risk_factors_at_study_entry)
  
  p.Dat$COVID19.1 <- ifelse(p.Dat$COVID19 == 1,0,1)
  p.Dat$COVID19.1 <- as.factor(p.Dat$COVID19.1)
  
  p.Dat$Vaccine1 <- relevel(as.factor(as.character(p.Dat$Vaccine1)),ref = "none")
  p.Dat$Vaccine2 <- relevel(as.factor(as.character(p.Dat$Vaccine2)),ref = "none")
  
  idxrem1 <- which(p.Dat$year == "2021" & p.Dat$Dose1 == 0 & p.Dat$Dose2 == 0)
  if (length(idxrem1 != 0)){
    p.Dat.vac <- p.Dat[-which(p.Dat$year == "2021" & p.Dat$Dose1 == 0 & p.Dat$Dose2 == 0),]}
  else {p.Dat.vac = p.Dat}
  idxvax12 <- which(p.Dat.vac$year == "2021" | (p.Dat.vac$year == "2020" & p.Dat.vac$Dose1 == 1 & p.Dat.vac$Dose2 == 0) | (p.Dat.vac$year == "2020" & p.Dat.vac$Dose2 == 1))
  idxvax1 <- which(p.Dat.vac$Dose1 == 1 & p.Dat.vac$Dose2 == 0)
  idxvax2 <- which(p.Dat.vac$Dose2 == 1)
  
  p.Dat.vac$vax.cat1 <- p.Dat.vac$vax.cat <- rep(0, nrow(p.Dat.vac))
  p.Dat.vac$vax.cat[idxvax12] <- 1
  p.Dat.vac$vax.cat1[idxvax1] <- 1
  p.Dat.vac$vax.cat1[idxvax2] <- 2
  
  return(p.Dat.vac)
}

Dat.ARSvac <- PrepVac(Dat.ARS)   
Dat.BPCvac <- PrepVac(Dat.BIFAPPC) 
Dat.BPCHvac <- PrepVac(Dat.BIFAPPC_HOSP)
Dat.CPRDvac <- PrepVac(Dat.CPRD) 
Dat.PHvac <- PrepVac(Dat.PHARMO)

## Make the AESI inline with the previous data 

obj <- readRDS("Z:/testPoisson/TestRCode/listofEstNEW.rds")

objdat <- data.frame(obj)
temp.dat <- cbind(CN[idxAESI],rownames(objdat)[seq(1,180,by=4)])

idxsel <- NULL
for (i in 1:45){
  idxsel[i] <- match(temp.dat[i,2],temp.dat[,1])
}


con.AESI <- function(data,idxsel,idxPT,idxAESI){
  idxPTused <- idxPT[2:length(idxPT)]
  
  Mat.Count.AESI <- matrix(NA,nrow(data),45)
  Mat.PT.AESI <- matrix(NA,nrow(data),45)
  
  Mat.Count.AESI <- data[,idxAESI[idxsel]]
  Mat.PT.AESI <- data[,idxPTused[idxsel]]
  new.dat <- cbind(data[,-c(idxAESI,idxPTused)],Mat.Count.AESI,Mat.PT.AESI)
  return(new.dat)}


Dat.ARSvac1 <- con.AESI(Dat.ARSvac,idxsel, idxPT,idxAESI)
Dat.BPCvac1 <- con.AESI(Dat.BPCvac,idxsel, idxPT,idxAESI)
Dat.BPCHvac1 <- con.AESI(Dat.BPCHvac,idxsel, idxPT,idxAESI)
Dat.CPRDvac1 <- con.AESI(Dat.CPRDvac,idxsel, idxPT,idxAESI)
Dat.PHvac1 <- con.AESI(Dat.PHvac,idxsel, idxPT,idxAESI)

CN2 <- colnames(Dat.ARSvac1)
idxPT2 <- grep("Persontime",colnames(Dat.ARSvac1))
idxAESI2 <- c(43:87)


### Flagging the AESI indices where DAP cannot identify the AESI.

flagging <- function(data){
  idxflag <- NULL 
  for (i in 1:length(idxAESI2)){
    temp <- CN2[idxAESI2[i]]
    data.used <- data
    sum.vec <- summary(data.used[,temp]) 
    l <- length(as.vector(sum.vec))
    if (l == 6) {
      if (prod(as.vector(sum.vec[1:6]) == rep(0,6)) == 1){
        idxflag <- c(idxflag,i)        
      }
    } else if (l == 7){
      if (sum.vec[7] == nrow(data.used)){
        idxflag <- c(idxflag,i)
      }
    }}
  return(idxflag)
}

idARSflag <- flagging(Dat.ARSvac1)
idBPCflag <- flagging(Dat.BPCvac1)
idBPCHflag <- flagging(Dat.BPCHvac1)
idCPRDflag <- flagging(Dat.CPRDvac1)
idPHflag <- flagging(Dat.PHvac1)


listflagnew <- list("ARS"=idARSflag, "BPC" = idBPCflag, "BPCH" = idBPCHflag, "CPRD" = idCPRDflag, "PH" = idPHflag)
saveRDS(listflag,"Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/listflagDAP.rds")





## The following is optional, in case Moderna and Pfizer will be combined to mRNA,
## and AZ with J&J were combined into adeno. 

comb.vax <- function(data){
  
  idxAdeno <- which(data$Vaccine1 == "AstraZeneca" | data$Vaccine2 == "AstraZeneca" | data$Vaccine1 == "J&J" | data$Vaccine2 == "J&J")
  idxmRNA <- which(data$Vaccine1 == "Moderna" | data$Vaccine2 == "Moderna" | data$Vaccine1 == "Pfizer" | data$Vaccine2 == "Pfizer")
  
  
  data$vaxtype <- rep("none",nrow(data))
  data$vaxtype[idxmRNA] <- "mRNA"
  data$vaxtype[idxAdeno] <- "adeno"
  data$vaxtype <- relevel(as.factor(as.character(data$vaxtype)),ref = "none")
  return(data)
  
}

Dat.ARSvac2 <- comb.vax(Dat.ARSvac1)
Dat.BPCvac2 <- comb.vax(Dat.BPCvac1)
Dat.BPCHvac2 <- comb.vax(Dat.BPCHvac1)
Dat.CPRDvac2 <- comb.vax(Dat.CPRDvac1)
Dat.PHvac2 <- comb.vax((Dat.PHvac1))








