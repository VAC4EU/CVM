#-----------------------------------------------
# Calculate IR for MIS/Myocard

# input: D4_persontime_b, D4_persontime_monthly_b_BC, D4_persontime_c, D4_persontime_monthly_c_BC, D4_persontime_d, D4_persontime_monthly_d_BC
# output: RES_IR_MIS_b, RES_IR_monthly_MIS_b, RES_IR_MIS_c, RES_IR_monthly_MIS_c, RES_IR_MIS_d, RES_IR_monthly_MIS_d

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  namedataset<-paste0("D4_persontime_b",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_b",suffix[[subpop]],".RData"))
  D4_persontime_b<-get(namedataset)

for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_b[, (name_cols) := exactPoiCI(D4_persontime_b, name_count, name_pt)]
}

  nameoutput<-paste0("RES_IR_MIS_b")
  assign(nameoutput,D4_persontime_b[, !grep("^Person", D4_persontime_b) , with = FALSE])
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)
  
  fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
  rm(list=nameoutput)

  

  namedataset<-paste0("D4_persontime_monthly_b_BC",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_monthly_b_BC",suffix[[subpop]],".RData"))
  D4_persontime_monthly_b_BC<-get(namedataset)

for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_monthly_b_BC[, (name_cols) := exactPoiCI(D4_persontime_monthly_b_BC, name_count, name_pt)]
}
  nameoutput<-paste0("RES_IR_monthly_MIS_b")

  
  assign(nameoutput,D4_persontime_monthly_b_BC)
  first_cols <- c("sex", "month", "year", "Ageband")
  all_cols <- colnames(get(nameoutput))
  setcolorder(get(nameoutput), c(first_cols, all_cols[all_cols %not in% first_cols]))
  
  save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)
  
  fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
  rm(list=nameoutput)
  


#D4_persontime c----------------------------------------------

  namedataset<-paste0("D4_persontime_c",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_c",suffix[[subpop]],".RData"))
  D4_persontime_c<-get(namedataset)

for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_c[, (name_cols) := exactPoiCI(D4_persontime_c, name_count, name_pt)]
}

  nameoutput<-paste0("RES_IR_MIS_c")
  assign(nameoutput,D4_persontime_c[, !grep("^Person", names(D4_persontime_c)) , with = FALSE])

  save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)
  
  fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
  rm(list=nameoutput)
  

  namedataset<-paste0("D4_persontime_monthly_c_BC",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_monthly_c_BC",suffix[[subpop]],".RData"))
  D4_persontime_monthly_c_BC<-get(namedataset)

for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_monthly_c_BC[, (name_cols) := exactPoiCI(D4_persontime_monthly_c_BC, name_count, name_pt)]
}

  nameoutput<-paste0("RES_IR_monthly_MIS_c")
  
  
  assign(nameoutput,D4_persontime_monthly_c_BC)
  first_cols <- c("sex", "month", "year", "Ageband")
  all_cols <- colnames(get(nameoutput))
  setcolorder(get(nameoutput), c(first_cols, all_cols[all_cols %not in% first_cols]))
  
  save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)
  
  fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
  rm(list=nameoutput)



#D4_persontime d----------------------------------------------
namedataset<-paste0("D4_persontime_d",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_d",suffix[[subpop]],".RData"))
D4_persontime_d<-get(namedataset)
  
for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_d[, (name_cols) := exactPoiCI(D4_persontime_d, name_count, name_pt)]
}

nameoutput<-paste0("RES_IR_MIS_d")


assign(nameoutput,D4_persontime_d[, !grep("^Person", names(D4_persontime_d)) , with = FALSE])

save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)

fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
rm(list=nameoutput)



namedataset<-paste0("D4_persontime_monthly_d_BC",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_monthly_d_BC",suffix[[subpop]],".RData"))
D4_persontime_monthly_d_BC<-get(namedataset)


for (ev in list_outcomes_MIS) {
  name_cols <- paste0(c("IR_", "lb_", "ub_"), ev)
  name_count <- paste0(ev,"_b")
  name_pt <- paste0("Persontime_",ev)
  D4_persontime_monthly_d_BC[, (name_cols) := exactPoiCI(D4_persontime_monthly_d_BC, name_count, name_pt)]
}

nameoutput<-paste0("RES_IR_monthly_MIS_d")


assign(nameoutput,D4_persontime_monthly_d_BC)
first_cols <- c("sex", "month", "year", "Ageband")
all_cols <- colnames(get(nameoutput))
setcolorder(get(nameoutput), c(first_cols, all_cols[all_cols %not in% first_cols]))

save(nameoutput,file=paste0(thisdirexp,nameoutput,".RData"),list=nameoutput)

fwrite(get(nameoutput),file=paste0(thisdirexp,nameoutput,".csv"))
rm(list=nameoutput)
rm(namedataset)


}

