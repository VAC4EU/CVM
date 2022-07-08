#-----------------------------------------------
# Aggregate persontime datasets by month

# input: D4_persontime_risk_month
# output: D4_persontime_risk_month_RFBC

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  namedataset1<-paste0("D4_persontime_risk_month",suffix[[subpop]])
  load(paste0(diroutput,"D4_persontime_risk_month",suffix[[subpop]],".RData"))
  
cols_to_sums = names(get(namedataset1))[5:length(get(namedataset1))]
setnames(get(namedataset1), "ageband_at_study_entry", "Ageband")

assign(namedataset1,get(namedataset1)[, c("year", "month") := tstrsplit(month, "-")])
get(namedataset1)[, month := month.name[as.integer(month)]]

all_sex <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("Ageband", "month", "year", "at_risk_at_study_entry"),
                                         .SDcols = cols_to_sums]
all_sex <- all_sex[, sex := "both_sexes"]
assign(namedataset1, rbind(get(namedataset1), all_sex))

all_risk <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "Ageband", "month", "year"),
                                          .SDcols = cols_to_sums]
all_risk <- all_risk[, at_risk_at_study_entry := "total"]
assign(namedataset1, rbind(get(namedataset1), all_risk))

all_year <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "Ageband", "year", "at_risk_at_study_entry"),
                                          .SDcols = cols_to_sums]
all_year <- all_year[, month := "all_months"]
assign(namedataset1, rbind(get(namedataset1), all_year))

all_ages <- copy(get(namedataset1))[, lapply(.SD, sum), by = c("sex", "month", "year", "at_risk_at_study_entry"),
                                             .SDcols = cols_to_sums]
all_ages <- unique(all_ages[, Ageband := "all_agebands"])
assign(namedataset1, rbind(get(namedataset1), all_ages))

assign(namedataset1, bc_divide_60(get(namedataset1), c("sex", "month", "year", "at_risk_at_study_entry"),
                                  cols_to_sums, col_used = "Ageband"))

nameoutput1<-paste0("D4_persontime_risk_month_RFBC",suffix[[subpop]])
assign(nameoutput1, get(namedataset1))
save(nameoutput1,file=paste0(diroutput,nameoutput1,".RData"),list=nameoutput1)
rm(list=nameoutput1)
rm(list=namedataset1)
rm(namedataset1,nameoutput1)

}



