#-----------------------------------------------
# Aggregate persontime datasets for MIS/Myocard

# input: D4_persontime_monthly_b, D4_persontime_monthly_c, D4_persontime_monthly_d
# output: D4_persontime_monthly_b_BC, D4_persontime_monthly_c_BC, D4_persontime_monthly_d_BC

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  namedataset1<-paste0("D4_persontime_monthly_b",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_monthly_b",suffix[[subpop]],".RData"))
  
cols_to_sums = names(get(namedataset1))[4:length(get(namedataset1))]
setnames(get(namedataset1), "ageband_at_1_jan_2021", "Ageband")

assign(namedataset1,get(namedataset1)[, c("year", "month") := tstrsplit(month, "-")])
get(namedataset1)[, month := month.name[as.integer(month)]]

all_sex <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("Ageband", "month", "year"),
                                          .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]
assign(namedataset1, rbind(get(namedataset1), all_sex))

all_month <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "Ageband", "year"),
                                           .SDcols = cols_to_sums]
all_month <- all_month[, month := "all_months"]
assign(namedataset1, rbind(get(namedataset1), all_month))

all_year <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "Ageband", "month"),
                                          .SDcols = cols_to_sums]
all_year <- all_year[, year := "all_years"]
assign(namedataset1, rbind(get(namedataset1), all_year))

all_ages <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "month", "year"),
                                           .SDcols = cols_to_sums]
all_ages <- unique(all_ages[, Ageband := "all_birth_cohorts"])

nameoutput1<-paste0("D4_persontime_monthly_b_BC",suffix[[subpop]])
assign(nameoutput1,rbind(get(namedataset1), all_ages))

save(nameoutput1,file=paste0(diroutput,nameoutput1,".RData"),list=nameoutput1)
rm(list=nameoutput1)
rm(list=namedataset1)
rm(namedataset1,nameoutput1)



namedataset2<-paste0("D4_persontime_monthly_c",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_monthly_c",suffix[[subpop]],".RData"))

cols_to_sums = names(get(namedataset2))[4:length(get(namedataset2))]
setnames(get(namedataset2), "ageband_at_1_jan_2021", "Ageband")

assign(namedataset2,get(namedataset2)[, c("year", "month") := tstrsplit(month, "-")])
get(namedataset2)[, month := month.name[as.integer(month)]]

all_sex <- copy(get(namedataset2))[, lapply(.SD, sum), by = c("Ageband", "month", "year"),
                                         .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]
assign(namedataset2,rbind(get(namedataset2), all_sex))

all_month <- copy(get(namedataset2))[, lapply(.SD, sum), by = c("sex", "Ageband", "year"),
                                           .SDcols = cols_to_sums]
all_month <- all_month[, month := "all_months"]
assign(namedataset2,rbind(get(namedataset2), all_month))

all_year <- copy(get(namedataset2))[, lapply(.SD, sum), by = c("sex", "Ageband", "month"),
                                          .SDcols = cols_to_sums]
all_year <- all_year[, year := "all_years"]
assign(namedataset2,rbind(get(namedataset2), all_year))

all_ages <- copy(get(namedataset2))[, lapply(.SD, sum), by = c("sex", "month", "year"),
                                          .SDcols = cols_to_sums]
all_ages <- unique(all_ages[, Ageband := "all_birth_cohorts"])

nameoutput2<-paste0("D4_persontime_monthly_c_BC",suffix[[subpop]])
assign(nameoutput2,rbind(get(namedataset2), all_ages))

save(nameoutput2,file=paste0(diroutput,nameoutput2,".RData"),list=nameoutput2)
rm(list=nameoutput2)
rm(list=namedataset2)
rm(namedataset2,nameoutput2)




namedataset3<-paste0("D4_persontime_monthly_d",suffix[[subpop]])
load(paste0(diroutput,"D4_persontime_monthly_d",suffix[[subpop]],".RData"))

cols_to_sums = names(get(namedataset3))[6:length(get(namedataset3))]
setnames(get(namedataset3), "ageband_at_1_jan_2021", "Ageband")

assign(namedataset3,get(namedataset3)[, c("year", "month") := tstrsplit(month, "-")])
get(namedataset3)[, month := month.name[as.integer(month)]]

all_sex <- copy(get(namedataset3))[, lapply(.SD, sum), by = c("Ageband", "month", "year", "type_vax_1", "history_covid"),
                                         .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]
assign(namedataset3,rbind(get(namedataset3), all_sex))

all_month <- copy(get(namedataset3))[, lapply(.SD, sum), by = c("sex", "Ageband", "year", "type_vax_1", "history_covid"),
                                          .SDcols = cols_to_sums]
all_month <- all_month[, month := "all_months"]
assign(namedataset3,rbind(get(namedataset3), all_month))

all_year <- copy(get(namedataset3))[, lapply(.SD, sum), by = c("sex", "Ageband", "month", "type_vax_1", "history_covid"),
                                          .SDcols = cols_to_sums]
all_year <- all_year[, year := "all_years"]
assign(namedataset3,rbind(get(namedataset3), all_year))

all_ages <- copy(get(namedataset3))[, lapply(.SD, sum), by = c("sex", "month", "year", "type_vax_1", "history_covid"),
                                          .SDcols = cols_to_sums]
all_ages <- unique(all_ages[, Ageband := "all_birth_cohorts"])

nameoutput3<-paste0("D4_persontime_monthly_d_BC",suffix[[subpop]])
assign(nameoutput3,rbind(get(namedataset3), all_ages))

save(nameoutput3,file=paste0(diroutput,nameoutput3,".RData"),list=nameoutput3)
rm(list=nameoutput3)
rm(list=namedataset3)
rm(namedataset3,nameoutput3)
rm(all_sex,all_ages,all_month,all_year)

}
