#---------------------------------------------------------------
# Apply all exclusion criteria

# input: persons_doses
# output: Flowchart_doses

#merge vaccine information to persons
for (subpop in subpopulations_non_empty){
  print(subpop)
  load(paste0(dirtemp,"persons_doses",suffix[[subpop]],".RData"))
  persons_doses_temp<-get(paste0("persons_doses", suffix[[subpop]]))

all_mondays <- seq.Date(as.Date("19000101","%Y%m%d"), Sys.Date(), by = "week")

monday_week <- seq.Date(from = find_last_monday(study_start, all_mondays), to = find_last_monday(study_end, all_mondays),
                        by = "week")
double_weeks <- data.table(weeks_to_join = monday_week, monday_week = monday_week)
all_days_df <- data.table(all_days = seq.Date(from = find_last_monday(study_start, monday_week), to = study_end, by = "days"))
all_days_df <- merge(all_days_df, double_weeks, by.x = "all_days", by.y = "weeks_to_join", all.x = T)
all_days_df <- all_days_df[, monday_week := nafill(monday_week, type="locf")]
all_days_df <- all_days_df[all_days >= study_start,]

##add the corresponding moday to each date
temp<-merge(persons_doses_temp,all_days_df, by.x=c("date"), by.y = "all_days")
monday_start_covid<-find_last_monday(start_COVID_vaccination_date,monday_week)
temp<-temp[monday_week>=monday_start_covid,]

temp2<-unique(temp[,doses_week:=.N ,by="monday_week"][,.(monday_week,doses_week)])
setnames(temp, "sex_or_birth_date_missing", "dose_not_in_persons")

selected_population <- CreateFlowChart(
  dataset = temp,
  listcriteria = c("duplicated_records", "missing_date", "date_before_start_vax", "distance_btw_1_2_doses",
                   "distance_btw_2_3_doses", "dose_after_3", "dose_after_2", "dose_not_in_persons", "birth_date_absurd",
                   "no_observation_period", "death_before_study_entry", "no_observation_period_including_study_start",
                   "insufficient_run_in", "death_before_vax", "exit_spell_before_vax", "study_end_before_vax"),
  strata = "monday_week",
  flowchartname = paste0("Flowchart_doses",suffix[[subpop]]))

thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])

  fwrite(get(paste0("Flowchart_doses",suffix[[subpop]])), paste0(thisdirexp,"Flowchart_doses.csv"))
  
  rm(list=paste0("Flowchart_doses",suffix[[subpop]]))
  rm(list=paste0("persons_doses", suffix[[subpop]]))
}




rm(persons_doses_temp, all_mondays, monday_week, double_weeks, all_days_df, temp, temp2, selected_population )