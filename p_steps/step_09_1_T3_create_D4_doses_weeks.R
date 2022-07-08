#-----------------------------------------------
# Create the dataset D4_doses_weeks

# input: D3_studyweeks
# output: D4_doses_weeks

for (subpop in subpopulations_non_empty) {  
  print(subpop)
  
  load(paste0(dirtemp,"D3_studyweeks",suffix[[subpop]],".RData"))

  
all_mondays <- seq.Date(as.Date("19000101","%Y%m%d"), Sys.Date(), by = "week")

monday_week <- seq.Date(from = find_last_monday(study_start, all_mondays), to = find_last_monday(study_end, all_mondays),
                        by = "week")
double_weeks <- data.table(weeks_to_join = monday_week, monday_week = monday_week)
all_days_df <- data.table(all_days = seq.Date(from = find_last_monday(study_start, monday_week), to = study_end, by = "days"))
all_days_df <- merge(all_days_df, double_weeks, by.x = "all_days", by.y = "weeks_to_join", all.x = T)
all_days_df <- all_days_df[, monday_week := nafill(monday_week, type="locf")]
all_days_df <- all_days_df[all_days >= study_start,]

cohort_to_vaxweeks <- get(paste0("D3_studyweeks",suffix[[subpop]]))
rm(list = paste0("D3_studyweeks",suffix[[subpop]]))

cohort_to_vaxweeks <- join_and_replace(cohort_to_vaxweeks, all_days_df, c("study_entry_date", "all_days"), "monday_week")
cohort_to_vaxweeks <- join_and_replace(cohort_to_vaxweeks, all_days_df, c("study_exit_date", "all_days"), "monday_week")
cohort_to_vaxweeks <- cohort_to_vaxweeks[, fup := as.numeric(difftime(study_exit_date, study_entry_date, units = "weeks") +  1)]
cohort_to_vaxweeks <- cohort_to_vaxweeks[, study_exit_date := NULL]


for (ageband in Agebands_labels) {
  nameoutput <- paste0("pop_age_", gsub("-", "_", ageband), suffix[[subpop]])
  assign(nameoutput, cohort_to_vaxweeks[ageband_at_study_entry == ageband, ])
  save(nameoutput, file = paste0(dirtemp, nameoutput,".RData"),list=nameoutput)
  rm(list=nameoutput)
}


df_events_ages <- paste0("pop_age_", gsub("-", "_", Agebands_labels))


rm(cohort_to_vaxweeks)

vect_to_rbind <- c()
for (name_temp_df in df_events_ages) {
  
  print(paste("computing for", name_temp_df))
  
  temp_df <- get(load(paste0(dirtemp, name_temp_df,suffix[[subpop]], ".RData"))[1])
  
  temp_df <- as.data.table(lapply(temp_df, rep, temp_df$fup))
  
  temp_df <- temp_df[, week := rowid(person_id, Dose, study_entry_date)]
  temp_df <- temp_df[, week := week - 1][, fup := NULL]
  temp_df <- temp_df[, start_date_of_period := study_entry_date + 7 * week]
  temp_df <- temp_df[, c("study_entry_date", "week") := NULL]
  
  tot_cohort <- copy(temp_df)
  tot_cohort <- tot_cohort[, c("sex", "ageband_at_study_entry", "Dose",
                               "type_vax", "at_risk_at_study_entry") := NULL]
  tot_cohort <- unique(tot_cohort)
  tot_cohort <- tot_cohort[, Persons_in_week := .N, by = c("start_date_of_period")]
  tot_cohort <- unique(tot_cohort[, person_id := NULL])
  
  temp_df <- temp_df[, Doses_in_week := .N,
                     by = c("sex", "ageband_at_study_entry", "Dose",
                            "type_vax", "start_date_of_period", "at_risk_at_study_entry")]
  
  
  temp_df <- unique(temp_df[, person_id := NULL])
  temp_df <- merge(temp_df, tot_cohort, by = "start_date_of_period")
  
  save(temp_df, file = paste0(dirtemp, name_temp_df,suffix[[subpop]], ".RData"))

  rm(temp_df, tot_cohort)
  rm(list=paste0( name_temp_df,suffix[[subpop]]))
}


for (name_temp_df in df_events_ages) {
  load(paste0(dirtemp, name_temp_df,suffix[[subpop]], ".RData"))
  vect_to_rbind <- append(vect_to_rbind, list(temp_df))
}

cohort_to_vaxweeks <- rbindlist(vect_to_rbind)
cohort_to_vaxweeks <- cohort_to_vaxweeks[Dose == 0, Doses_in_week := 0]

cohort_to_vaxweeks <- bc_divide_60(cohort_to_vaxweeks,
                                   c("start_date_of_period", "sex", "type_vax", "at_risk_at_study_entry", "Dose"),
                                   c("Doses_in_week", "Persons_in_week"))

older60 <- unique(cohort_to_vaxweeks[ageband_at_study_entry %in% Agebands60, c("start_date_of_period", "Persons_in_week")])
older60 <- older60[, temp_var := lapply(.SD, sum, na.rm=TRUE), by = "start_date_of_period", .SDcols = "Persons_in_week"]
older60 <- unique(older60[, ageband_at_study_entry := "60+"][, Persons_in_week := NULL])

younger60 <- unique(cohort_to_vaxweeks[ageband_at_study_entry %in% Agebands059, c("start_date_of_period", "Persons_in_week")])
younger60 <- younger60[, temp_var := lapply(.SD, sum, na.rm=TRUE), by = "start_date_of_period", .SDcols = "Persons_in_week"]
younger60 <- unique(younger60[, ageband_at_study_entry := "0-59"][, Persons_in_week := NULL])

cohort_to_vaxweeks <- merge(cohort_to_vaxweeks, older60, all.x = T, by = c("start_date_of_period", "ageband_at_study_entry"))
cohort_to_vaxweeks <- cohort_to_vaxweeks[ageband_at_study_entry == "60+", Persons_in_week := temp_var][, temp_var := NULL]
cohort_to_vaxweeks <- merge(cohort_to_vaxweeks, younger60, all.x = T, by = c("start_date_of_period", "ageband_at_study_entry"))
cohort_to_vaxweeks <- cohort_to_vaxweeks[ageband_at_study_entry == "0-59", Persons_in_week := temp_var][, temp_var := NULL]

setorder(cohort_to_vaxweeks, sex, ageband_at_study_entry, Dose,
         type_vax, start_date_of_period, at_risk_at_study_entry)

cohort_to_vaxweeks <- cohort_to_vaxweeks[, Year := year(start_date_of_period)]
cohort_to_vaxweeks <- cohort_to_vaxweeks[, Datasource := thisdatasource]

setnames(cohort_to_vaxweeks, c("start_date_of_period", "type_vax", "sex", "at_risk_at_study_entry"),
         c("Week_number", "Type_vax", "Sex", "At_Risk"))
cohort_to_vaxweeks <- cohort_to_vaxweeks[, sex := fifelse(Sex == 1, "Male", "Female")]

D4_doses_weeks <- cohort_to_vaxweeks[, .(Datasource, Year, Week_number, ageband_at_study_entry, Sex, At_Risk, Dose,
                                         Type_vax, Persons_in_week, Doses_in_week)]

thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
nameoutput <- paste0("D4_doses_weeks")
assign(nameoutput, D4_doses_weeks)
save(nameoutput, file = paste0(thisdirexp, nameoutput,".RData"),list=nameoutput)
rm(list=nameoutput)

}

rm( all_mondays, monday_week, double_weeks, all_days_df, cohort_to_vaxweeks,temp_df)