#-----------------------------------------------
# Aggregate persontime datasets by risk factors

# input: D4_persontime_risk_week, D4_persontime_benefit_week, D4_persontime_risk_year, D4_persontime_benefit_year
# output: D4_persontime_risk_week_RF, D4_persontime_benefit_week_RF, D4_persontime_risk_year_RF, D4_persontime_benefit_year_RF


for (subpop in subpopulations_non_empty) {  
  print(subpop)
  namedataset1<-paste0("D4_persontime_risk_week",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_risk_week",suffix[[subpop]],".RData"))
  
D4_persontime_risk_week <- get(namedataset1)[any_risk_factors != 0]

cols_to_sums = names(get(namedataset1))[16:length(get(namedataset1))]

assign(namedataset1, get(namedataset1)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week", "CV", "COVCANCER",
                                                          "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                                                          "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                                   .SDcols = cols_to_sums])

assign(namedataset1, melt(get(namedataset1),
                                    measure.vars = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                     "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                    variable.name = "riskfactor", value.name = "to_drop"))

assign(namedataset1,get(namedataset1)[to_drop == 1, ][, to_drop := NULL])

       assign(namedataset1,get(namedataset1)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week", "riskfactor"),
                                                   .SDcols = cols_to_sums])

all_sex <- copy(get(namedataset1))[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("Dose", "type_vax", "week", "riskfactor"),
                                         .SDcols = cols_to_sums]

all_sex <- all_sex[, sex := "both_sexes"]

nameoutput1<-paste0("D4_persontime_risk_week_RF",suffix[[subpop]])
assign(nameoutput1, rbind(get(namedataset1), all_sex))

save(nameoutput1,file=paste0(diroutput,nameoutput1,".RData"),list=nameoutput1)
rm(list=nameoutput1)
rm(list=namedataset1)
rm(namedataset1,nameoutput1)


namedataset2<-paste0("D4_persontime_benefit_week",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_benefit_week",suffix[[subpop]],".RData"))

assign(namedataset2,get(namedataset2)[any_risk_factors != 0])

cols_to_sums = names(get(namedataset2))[16:length(get(namedataset2))]

assign(namedataset2,get(namedataset2)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week", "CV", "COVCANCER",
                                                          "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                                                          "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                                   .SDcols = cols_to_sums])

assign(namedataset2,melt(get(namedataset2),
                                measure.vars = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                 "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                variable.name = "riskfactor", value.name = "to_drop"))

assign(namedataset2,get(namedataset2)[to_drop == 1, ][, to_drop := NULL])

assign(namedataset2,get(namedataset2)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week", "riskfactor"),
                                                   .SDcols = cols_to_sums])

all_sex <- copy(get(namedataset2))[, lapply(.SD, sum, na.rm=TRUE),
                                            by = c("Dose", "type_vax", "week", "riskfactor"),
                                            .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]

nameoutput2<-paste0("D4_persontime_benefit_week_RF",suffix[[subpop]])
assign(nameoutput2, rbind(get(namedataset2), all_sex))

save(nameoutput2,file=paste0(diroutput,nameoutput2,".RData"),list=nameoutput2)
rm(list=nameoutput2)
rm(list=namedataset2)
rm(namedataset2,nameoutput2)


namedataset3<-paste0("D4_persontime_risk_year",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_risk_year",suffix[[subpop]],".RData"))


assign(namedataset3, get(namedataset3)[any_risk_factors != 0, ][Dose != 0 | year == 2020, ][, year := NULL])
cols_to_sums = names(get(namedataset3))[16:length(get(namedataset3))]

assign(namedataset3, get(namedataset3)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week_fup", "CV", "COVCANCER",
                                                          "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                                                          "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                                   .SDcols = cols_to_sums])

assign(namedataset3,melt(get(namedataset3),
                                measure.vars = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                 "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                variable.name = "riskfactor", value.name = "to_drop"))

assign(namedataset3,get(namedataset3)[to_drop == 1, ][, to_drop := NULL])

assign(namedataset3,get(namedataset3)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week_fup", "riskfactor"),
                                                   .SDcols = cols_to_sums])

all_sex <- copy(get(namedataset3))[, lapply(.SD, sum, na.rm=TRUE),
                                         by = c("Dose", "type_vax", "week_fup", "riskfactor"),
                                         .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]

nameoutput3<-paste0("D4_persontime_risk_year_RF",suffix[[subpop]])
assign(nameoutput3,rbind(get(namedataset3), all_sex))

setorder(get(nameoutput3), "week_fup")

vax_dose <- unique(copy(get(nameoutput3))[, c("Dose", "type_vax", "week_fup")])
vax_dose <- vax_dose[, .SD[which.max(week_fup)], by = c("Dose", "type_vax")]
vax_dose_0 <- vax_dose[week_fup == 0,]
week_vax_dose <- vax_dose[, lapply(.SD, rep, vax_dose[, week_fup])][, week_fup := unlist(
  lapply(vax_dose[, week_fup], seq_len))]
week_vax_dose <- rbind(week_vax_dose, vax_dose_0)

sex_vect <- rep(c("0", "1", "both_sexes"), each = nrow(week_vax_dose))
week_vax_dose <- week_vax_dose[, lapply(.SD, rep, 3)][, sex := sex_vect]
riskfactors <- c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES", "COVSICKLE", "IMMUNOSUPPR",
                 "any_risk_factors")
riskfactor_vect <- rep(riskfactors, each = nrow(week_vax_dose))
empty_risk_year <- week_vax_dose[, lapply(.SD, rep, length(riskfactors))][, riskfactor := riskfactor_vect]

assign(nameoutput3, merge(empty_risk_year, get(nameoutput3), all.x = T,
                                    by = c("Dose", "type_vax", "week_fup", "sex", "riskfactor")))

for (i in names(get(nameoutput3))){
  get(nameoutput3)[is.na(get(i)), (i) := 0]
}

assign(nameoutput3, get(nameoutput3)[, (cols_to_sums) := lapply(.SD, cumsum),
                                     by = c("Dose", "type_vax", "riskfactor", "sex"),
                                     .SDcols = cols_to_sums])

save(nameoutput3,file=paste0(diroutput,nameoutput3,".RData"),list=nameoutput3)
rm(list=nameoutput3)
rm(list=namedataset3)
rm(namedataset3,nameoutput3)


namedataset4<-paste0("D4_persontime_benefit_year",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_benefit_year",suffix[[subpop]],".RData"))


assign(namedataset4,get(namedataset4)[any_risk_factors != 0][Dose != 0 | year == 2020, ][, year := NULL])
cols_to_sums = names(get(namedataset4))[16:length(get(namedataset4))]

assign(namedataset4,get(namedataset4)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week_fup", "CV", "COVCANCER",
                                                          "COVCOPD", "COVHIV", "COVCKD", "COVDIAB", "COVOBES",
                                                          "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                                   .SDcols = cols_to_sums])

assign(namedataset4,melt(get(namedataset4),
                                measure.vars = c("CV", "COVCANCER", "COVCOPD", "COVHIV", "COVCKD", "COVDIAB",
                                                 "COVOBES", "COVSICKLE", "IMMUNOSUPPR", "any_risk_factors"),
                                variable.name = "riskfactor", value.name = "to_drop"))

assign(namedataset4,get(namedataset4)[to_drop == 1, ][, to_drop := NULL])

assign(namedataset4,get(namedataset4)[, lapply(.SD, sum, na.rm=TRUE),
                                                   by = c("sex", "Dose", "type_vax", "week_fup", "riskfactor"),
                                                   .SDcols = cols_to_sums])

all_sex <- copy(get(namedataset4))[, lapply(.SD, sum, na.rm=TRUE),
                                            by = c("Dose", "type_vax", "week_fup", "riskfactor"),
                                            .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]

nameoutput4<-paste0("D4_persontime_benefit_year_RF",suffix[[subpop]])
assign(nameoutput4,rbind(get(namedataset4), all_sex))

save(nameoutput4,file=paste0(diroutput,nameoutput4,".RData"),list=nameoutput4)
rm(list=nameoutput4)
rm(list=namedataset4)
rm(namedataset4,nameoutput4)


rm(all_sex, cols_to_sums,empty_risk_year,vax_dose,vax_dose_0,week_vax_dose)


}
