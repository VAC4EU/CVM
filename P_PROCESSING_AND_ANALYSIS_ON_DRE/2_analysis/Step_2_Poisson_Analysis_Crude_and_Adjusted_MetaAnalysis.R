# Poisson Analysis : crude and adjusted, 

# --Program Information -----------------------------------------------
#       
#       Title                    : Fitting Poisson regression
#       Author (date created)    : Ivonne Martin
#       Input                    : 5 datatable from each DAP after preprocess
#       Output                   : 1 table for crude analysis, 1 table for adjusted analysis 
# 
# --Program Information -----------------------------------------------


# Function create.CI receives matrix of two columns in which columns are estimates 
# and their standard errors and rows are vaccine brands. The function returns matrix 
# of three columns (estimate, lower.95%.CI and upper.95%.CI)
# and rows of vaccine brands. 

create.CI <- function(pars){
  CI.bands <- matrix(NA,nrow = nrow(pars),ncol = 2)
  for (i in 1:nrow(pars)){
    CI.bands[i,] <- pars[i,1] + c(-1,1)*1.96*pars[i,2]
  }
  res <- cbind(pars[,1],CI.bands)
  colnames(res) <- c("Est","lower","upper")
  return(res)
}



# Check.fun  : A function to check for which DAP, the corresponding AESI is detected. 
# If it is not detected then the Poisson analysis for the corresponding DAP is ignored.

check.fun <- function(idx){
  check.vec <- c(is.na(match(idx,listflagnew$PH)) == 1 , is.na(match(idx,listflagnew$BPC)) == 1 , 
                 is.na(match(idx,listflagnew$BPCH)) == 1 , is.na(match(idx,listflagnew$CPRD)) == 1 ,
                 is.na(match(idx,listflagnew$ARS)) == 1 )
  sel <- which(check.vec == "FALSE")
  if (length(sel) != 0){
    val <- sel
  } else {
    val <- 0
  } 
  return(val)
}


# print.est : A function to print the estimate and 95%CI for the table.   

print.est <- function(pars){
  idx <- which(pars[,3] >=50)
  pars[idx,3] <- Inf 
  return(paste0(sprintf("%.2f",pars[,1]),sprintf(" (%.2f",pars[,2]),sprintf(",%.2f",pars[,3]),")"))
}


# Poisson analysis for each DAP 
ARSDap <- function(fml){
  
  Dat <- Dat.ARSvac1
  mod <- try(glm(fml, family = poisson, data = Dat,maxit = 100))              
  return(summary(mod)$coefficients)}

BPCDap <- function(fml){
  
  Dat <- Dat.BPCvac1
  mod <- try(glm(fml, family = poisson, data = Dat,maxit = 100))              
  return(summary(mod)$coefficients)}

BPCHDap <- function(fml){
  Dat <- Dat.BPCHvac1
  mod <- try(glm(fml, family = poisson, data = Dat,maxit = 100))              
  return(summary(mod)$coefficients)}

CPDap <- function(fml){
  Dat <- Dat.CPRDvac1
  mod <- try(glm(fml, family = poisson, data = Dat,maxit = 100))              
  return(summary(mod)$coefficients)}

PHDap <- function(fml){
  Dat <- Dat.PHvac1
  mod <- try(glm(fml, family = poisson, data = Dat,maxit = 100))              
  return(summary(mod)$coefficients)}

## Function to select model based on the availability of the AESI in each DAP, this is for the adjusted analysis

Sel.Mod <- function(tvalue, fml){ 
  test.val <- tvalue
  spec.row <- c(2:5)
  spec.col <- c(1:2)
  lrow <- length(spec.row)
  empty.DAP <- matrix(log(0),nrow = lrow, ncol = 3)
  empty.DAPB <- matrix(log(0),nrow = lrow, ncol = 2)
  
  
  if(length(test.val) == 1){
    if(prod(test.val == 0) == 1){
      
      Afit <- ARSDap(fml)
      B1fit <- BPCDap(fml)
      B2fit <- BPCHDap(fml)
      Cfit <- CPDap(fml)
      Pfit <- PHDap(fml)
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      B1fit[spec.row,spec.col],
                      B2fit[spec.row,spec.col],
                      rbind(Cfit[2,1:2],rep(log(0),2),Cfit[3:4,1:2]),
                      Pfit[spec.row,spec.col])
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
    }
    
    if(prod(test.val == 1) == 1){
      Afit <- ARSDap(fml)
      B1fit <- BPCDap(fml)
      B2fit <- BPCHDap(fml)
      Cfit <- CPDap(fml)
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      B1fit[spec.row,spec.col],
                      B2fit[spec.row,spec.col],
                      rbind(Cfit[2,1:2],rep(log(0),2),Cfit[3:4,1:2]),
                      empty.DAPB)
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    } 
    
    if (prod(test.val == 4) == 1){
      Afit <- ARSDap(fml)
      B1fit <- BPCDap(fml)
      B2fit <- BPCHDap(fml)
      Pfit <- PHDap(fml)
      
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      B1fit[spec.row,spec.col],
                      B2fit[spec.row,spec.col],
                      empty.DAPB,
                      Pfit[spec.row,spec.col])
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    }
  }
  ## end of 1 
  
  if(length(test.val) == 2){
    if(prod(test.val == c(1,3)) == 1){
      Afit <- ARSDap(fml)
      B1fit <- BPCDap(fml)
      Cfit <- CPDap(fml)
      
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      B1fit[spec.row,spec.col],
                      empty.DAPB,
                      rbind(Cfit[2,1:2],rep(log(0),2),Cfit[3:4,1:2]),
                      empty.DAPB)
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    }
    
    if(prod(test.val == c(2,3)) == 1){
      Afit <- ARSDap(fml)
      Cfit <- CPDap(fml)
      Pfit <- PHDap(fml)
      
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      empty.DAPB,
                      empty.DAPB,
                      rbind(Cfit[2,1:2],rep(log(0),2),Cfit[3:4,1:2]),
                      Pfit[spec.row,spec.col])
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    }
  }
  
  
  ## end of 2
  
  if (length(test.val) == 3){
    if(prod(test.val == c(1,2,3))){
      
      Afit <- ARSDap(fml)
      Cfit <- CPDap(fml)
      
      
      B.list <- cbind(Afit[spec.row,spec.col],
                      empty.DAPB,
                      empty.DAPB,
                      rbind(Cfit[2,1:2],rep(log(0),3),Cfit[3:4,1:2]),
                      empty.DAPB)
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    } 
    
  }
  
  ## end of 3 
  if (length(test.val) == 4){
    if(prod(test.val == c(1,2,4,5))){
      
      B2fit <- BPCHDap(fml)
      
      
      B.list <- cbind(empty.DAPB,
                      empty.DAPB,
                      B2fit[spec.row,spec.col],
                      empty.DAPB,
                      empty.DAPB)
      A.list <- cbind(create.CI(B.list[,1:2]),create.CI(B.list[,3:4]),create.CI(B.list[,5:6]),
                      create.CI(B.list[,7:8]),create.CI(B.list[,9:10]))
      
    } 
    
  }
  
  return(list(A.list,B.list))}


## In the crude analysis, the Sel.Mod function is slightly adjusted 
#--------------- Sel.Mod.Crude -------------------------------------


Sel.Mod.crude <- function(tvalue, fml){ ## Function to select model based on the availability of theAESI in each DAP
  test.val <- tvalue
  spec.row <- c(2:5)
  spec.col <- c(1:2)
  lrow <- length(spec.row)
  empty.DAP <- matrix(log(0),nrow = lrow, ncol = 3)
  
  if(length(test.val) == 1){
    if(prod(test.val == 0) == 1){
      
      CPfit <- create.CI(CPDap(fml)[c(2:4),spec.col])
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      create.CI(BPCDap(fml)[spec.row,spec.col]),
                      create.CI(BPCHDap(fml)[spec.row,spec.col]),
                      rbind(CPfit[1,],rep(log(0),3),CPfit[2:3,]),
                      create.CI(PHDap(fml)[spec.row,spec.col]))
    }
    
    if(prod(test.val == 1) == 1){
      CPfit <- create.CI(CPDap(fml)[c(2:4),spec.col])
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      create.CI(BPCDap(fml)[spec.row,spec.col]),
                      create.CI(BPCHDap(fml)[spec.row,spec.col]),
                      rbind(CPfit[1,],rep(log(0),3),CPfit[2:3,]),
                      empty.DAP)
    } 
    
    if (prod(test.val == 4) == 1){
      
      
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      create.CI(BPCDap(fml)[spec.row,spec.col]),
                      create.CI(BPCHDap(fml)[spec.row,spec.col]),empty.DAP,
                      create.CI(PHDap(fml)[spec.row,spec.col]))
    }
  }
  ## end of 1 
  
  if(length(test.val) == 2){
    if(prod(test.val == c(1,3)) == 1){
      CPfit <- create.CI(CPDap(fml)[c(2:4),spec.col])
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      create.CI(BPCDap(fml)[spec.row,spec.col]),empty.DAP,
                      rbind(CPfit[1,],rep(log(0),3),CPfit[2:3,]),empty.DAP)
    }
    
    if(prod(test.val == c(2,3)) == 1){
      CPfit <- create.CI(CPDap(fml)[c(2:4),spec.col])
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      empty.DAP,empty.DAP,
                      rbind(CPfit[1,],rep(log(0),3),CPfit[2:3,]),
                      create.CI(PHDap(fml)[spec.row,spec.col]))
    }
    
  }
  
  
  ## end of 2
  
  if (length(test.val) == 3){
    if(prod(test.val == c(1,2,3))){
      CPfit <- create.CI(CPDap(fml)[c(2:4),spec.col])
      A.list <- cbind(create.CI(ARSDap(fml)[spec.row,spec.col]),
                      empty.DAP, empty.DAP,
                      rbind(CPfit[1,],rep(log(0),3),CPfit[2:3,]),empty.DAP)
    } 
    
  }
  
  ## end of 3 
  if (length(test.val) == 4){
    if(prod(test.val == c(1,2,4,5))){
      A.list <- cbind(empty.DAP,
                      empty.DAP, 
                      create.CI(BPCHDap(fml)[spec.row,spec.col]), empty.DAP,empty.DAP)
    } 
    
  }
  
  return(A.list)}



#-------------------------------------------------------------------

idxusedAESI <- idxAESI2
list.est <- list.est.se <- vector("list",length(idxusedAESI))

for (i in 1:length(idxusedAESI)){
  if(i != 15){
    test.val <- check.fun(i)
    temp <- CN2[idxusedAESI[i]]
    AESIname <- substr(temp,1,nchar(temp) - 2)
    PT <- paste0("Persontime_",AESIname)
    
    # For the crude analysis, use the following cova1
    # cova1 <- c("Vaccine1",paste0("offset(log(",PT,"+1))"))

    cova1 <- c("Vaccine1","AgeCat","Gender","COVID19.1","risk.fac.new",paste0("offset(log(",PT,"+1))"))
    
    fml <- as.formula(paste0(temp,"~",paste0(cova1, collapse = "+")))
    
    AA.list <- Sel.Mod(test.val,fml)
    # AA.list <- Sel.Mod.crude(test.val, fml) ## for the crude analysis. 
  }
  list.est[[i]] <- AA.list[[1]]
  list.est.se[[i]] <- AA.list[[2]]
  print(i)
}

saveRDS(list.est, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220124_listestADJPoisson.rds")
saveRDS(list.est.se, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220124_listestseADJPoisson.rds")

## exponentiate the estimate 

list.est.exp <- lapply(list.est,exp)
list.est.se.exp <- lapply(list.est.se, exp)

print.est <- function(pars){
  idx <- which(pars[,3] >=50)
  pars[idx,3] <- pars[idx,3] 
  return(paste0(sprintf("%.2f",pars[,1]),sprintf(" (%.2f",pars[,2]),sprintf(",%.2f",pars[,3]),")"))
}


## Counting the number of nonzero events postvaccination
CountAESI <- function(Dat){
  nonzero <- vector("list",length(idxAESI2))
  for (i in 1:length(idxAESI2)){
    temp <- CN2[idxAESI2[i]]
    row.vec <- as.numeric(names(table(Dat[,temp])))
    Atab <- table(Dat[,temp],Dat$Vaccine1)
    if (nrow(Atab) == 0) {
      #nonzero[[i]] <- c("AstraZeneca" = -1,"J&J" = -1, "Moderna" = -1, "Pfizer" = -1)
      nonzero[[i]] <- c("AstraZeneca" = 0,"J&J" = 0, "Moderna" = 0, "Pfizer" = 0)
    } else if (nrow(Atab) == 1) {
      nonzero[[i]] <- c("AstraZeneca" = 0,"J&J" = 0, "Moderna" = 0, "Pfizer" = 0)
    } else if (nrow(Atab) == 2){
      nonzero[[i]] <- Atab[2:nrow(Atab),2:ncol(Atab)]*row.vec[2:length(row.vec)]
    } else {
      nonzero[[i]] <- colSums(Atab[2:nrow(Atab),2:ncol(Atab)]*row.vec[2:length(row.vec)]) 
    }
    print(i)
  }
  return(nonzero)
}

nz.ARS <- CountAESI(Dat.ARSvac1)
nz.BPC <- CountAESI(Dat.BPCvac1)
nz.BPCH <- CountAESI(Dat.BPCHvac1)
nz.CP <- CountAESI(Dat.CPRDvac1)
nz.PH <- CountAESI(Dat.PHvac1)


list.est.exp.print <- NULL

for (i in 1:length(idxAESI2)){
  l.exp <- list.est.exp[[i]]
  nz.CP1 <- nz.CP[[i]]
  list.est.exp.print <- rbind(list.est.exp.print,
                              cbind(nz.ARS[[i]], print.est(l.exp[,1:3]), nz.BPC[[i]], print.est(l.exp[,4:6]), nz.BPCH[[i]], 
                                    print.est(l.exp[,7:9]), c(nz.CP1[1:2],0,nz.CP1[3]),print.est(l.exp[,10:12]), nz.PH[[i]],print.est(l.exp[,13:15])))
  print(i)
}

list.print.new <- cbind(rep(c("AZ","J&J","Moderna","Pfizer"),45),list.est.exp.print)
rownames(list.print.new) <- rep(CN2[idxAESI2],each = 4)
colnames(list.print.new) <- c("Vaccine",rep(c("Events","IR (95%CI)"),5))
saveRDS(list.print.new, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220124_PRINTlistestADJPoisson.rds")



### Meta analysis
library(meta)

## Selection based on nonzero counts
  stud.name <- function(idx){
    if (length(idx) !=0){
      ref.name <- c("ARS","BIFAP_PC","CPRD","PHARMO")
      val <- ref.name[-idx]} else {
        val <- ref.name
      }
    return(val)
  }
  
  vaccine.brand <- function(idx){
    ref.name <- c("AZ","J&J","Moderna","Pfizer")
    return(ref.name[idx])
  }
  
  per.row.check<-function(mat,M.est,M.se,aesiname){
    list.ind <- vector("list",length = 4)
    est.fixed <- est.random <- detail.vec <- NULL
    for (i in 1:4){
      TE1 <- M.est[i,]
      seTE1 <- M.se[i,]
      idx.sel <- which(mat[i,] == 0)
      if(length(idx.sel) != 0){
        if(length(idx.sel) != 4)
        {TE <- TE1[-idx.sel]
        seTE <- seTE1[-idx.sel]
        } else {
          TE <- rep(NA,4)
          seTE <- rep(NA,4)
        }
      } else {
        TE <- TE1
        seTE <- seTE1
      }        
      
      stud.lab <- stud.name(idx.sel)
      vax.br <- vaccine.brand(i)
      
      if (length(stud.lab) != 0){
        png(file = paste0("Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/Plots/AESI_",aesiname,"_Vax_",vax.br,".png"), width = 720, height = 500, units = "px")
        m1 <- metagen(TE,seTE, sm = "IRR",backtransf = TRUE,studlab = stud.lab)
        forest(m1)
        grid.text(paste0("AESI_",aesiname,"Vax_",vax.br),.5,.9,gp = gpar(cex=1.5))
        dev.off()
        est.fixed <- rbind(est.fixed,exp(c(m1$TE.fixed,m1$lower.fixed,m1$upper.fixed)))
        est.random <- rbind(est.random,exp(c(m1$TE.random,m1$lower.random,m1$upper.random)))
        detail.vec <- rbind(detail.vec, c("fixed.pval" = m1$pval.fixed,"ran.pval" = m1$pval.random,"I2" = m1$I2,"Q.pval" = m1$pval.Q))
      }else {
        m1 <- metagen(TE,seTE, sm = "IRR",backtransf = TRUE)
        est.fixed <- rbind(est.fixed,exp(c(m1$TE.fixed,m1$lower.fixed,m1$upper.fixed)))
        est.random <- rbind(est.random,exp(c(m1$TE.random,m1$lower.random,m1$upper.random)))
        detail.vec <- rbind(detail.vec, c("fixed.pval" = m1$pval.fixed,"ran.pval" = m1$pval.random,"I2" = m1$I2,"Q.pval" = m1$pval.Q))
      }
    }
    return(list(est.fixed,est.random,detail.vec))  }
  
  
  
# Result of meta analysis
  
  Result <- vector("list",length = length(idxAESI2))
  for (i in 1:45){
    if(i!=15){
      tmp.name <- CN2[idxusedAESI[i]]
      aesiname <- substr(tmp.name,1,nchar(tmp.name) - 2)
      tmp1 <- list.est.se[[i]]
      tmp1.red <- tmp1[,-c(5:6)]
      nz.CP1 <- nz.CP[[i]]
      nz.All <- cbind(nz.ARS[[i]],nz.BPC[[i]],c(nz.CP1[1],0,nz.CP1[2:3]),nz.PH[[i]])
      
      m.est <- tmp1.red[,seq(1,ncol(tmp1.red),by=2)]
      m.se <- tmp1.red[,seq(2,ncol(tmp1.red),by=2)]
      
      Result[[i]] <- per.row.check(nz.All,m.est,m.se,aesiname)
    } else {
      Result[[i]] <- list(matrix(NA,nrow = 4,ncol = 3),matrix(NA,nrow=4,ncol=3), matrix(NA,nrow=4,ncol=4))
    }
    
    print(i)}
  
  
  ## Presenting the result with additional pooled effect size
  
  list.est.print <- readRDS("Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220124_PRINTlistestADJPoisson.rds")
  
  l.fix <- l.rnd <- NULL
  I2.stat <- c()
  for (i in 1:45){
    
    l.fix <- c(l.fix,print.est(Result[[i]][[1]]))
    l.rnd <- c(l.rnd,print.est(Result[[i]][[2]]))
    I2.stat <- rbind(I2.stat,Result[[i]][[3]][,3:4]) 
    print(i)
  }
  
  list.est.print.add <- cbind(list.est.print, "Pool.Fixed" = l.fix,"Pool.Rnd" = l.rnd)
  saveRDS(list.est.print.add, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220125_PRINTlistestADJPoissonMeta.rds")
  
  list.est.print.add.csv <- cbind(list.est.print.add,I2.stat)
  write.csv(list.est.print.add.csv, file = "Z:/folder_analysis_for_rebuttal/2_analysis/pooled_input_from_DAPs/RCodes/ObjectUsed/220125_PRINTlistestADJPoissonMeta.csv")
  