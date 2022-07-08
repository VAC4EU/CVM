# -----------------------------------------------------
# CREATE EXCLUSION CRITERIA for persons/spells

# input: D3_PERSONS, OBSERVATION_PERIODS, output_spells_category
# output: D3_selection_criteria

print('CREATE EXCLUSION CRITERIA')

load(paste0(dirtemp,"D3_PERSONS.RData"))

OBSERVATION_PERIODS <- data.table()
files<-sub('\\.csv$', '', list.files(dirinput))
for (i in 1:length(files)) {
  if (stringr::str_detect(files[i],"^OBSERVATION_PERIODS")) {  
    temp <- fread(paste0(dirinput,files[i],".csv"), colClasses = list(character="person_id"))
    OBSERVATION_PERIODS <- rbind(OBSERVATION_PERIODS, temp,fill=T)
    rm(temp)
  }
}

#CHANGE COLUMN NAMES
setnames(D3_PERSONS, c("date_birth", "date_death"), c("date_of_birth", "date_of_death"))

#CONVERT SEX to BINARY 0/1
D3_PERSONS<-D3_PERSONS[, sex := fifelse(sex_at_instance_creation == "M", 1, 0)] #1:M 0:F
#[,age_at_index_date:=age_fast(date_of_birth,index_date)][age_at_index_date<12 | age_at_index_date>55,age:=1]

D3_PERSONS<-D3_PERSONS[, sex_or_birth_date_missing:= fifelse(is.na(sex) | is.na(date_of_birth), 1, 0)]
D3_PERSONS<-D3_PERSONS[, birth_date_absurd := fifelse(year(date_of_birth) < 1899 | year(date_of_birth) > 2021, 1, 0)]

# no observation period (NA + )
PERSONS_in_OP <- merge(D3_PERSONS, OBSERVATION_PERIODS, all.x = T, by="person_id")
PERSONS_in_OP <- PERSONS_in_OP[, no_observation_period := fifelse(is.na(op_start_date), 1, 0)][, Minop_start_date := min(no_observation_period), by="person_id"][Minop_start_date==no_observation_period,]

D3_exclusion_no_observation_period <- unique(PERSONS_in_OP[,.(person_id, sex_or_birth_date_missing, birth_date_absurd, no_observation_period)])

## KEEP ONLY NEEDED VARs
D3_inclusion_from_PERSONS <- D3_PERSONS[,.(person_id,sex,date_of_birth,date_of_death)]

# if (this_datasource_has_subpopulations == TRUE)D3_selection_criteria <- vector(mode="list")

# if (this_datasource_has_subpopulations == TRUE)D3_selection_criteria <- vector(mode="list")

# OBSERVATION PERIODS -----------------------------------------------------
#new
 for (subpop in subpopulations_non_empty){
  print(subpop)
   if (this_datasource_has_subpopulations == T){
load(paste0(dirtemp,"output_spells_category.RData"))
output_spells_category<-output_spells_category[[subpop]]
   }else{
     load(paste0(dirtemp,"output_spells_category.RData"))
}

  na_date = lubridate::ymd(99991231)
  
  output_spells_category_enriched <- merge(D3_inclusion_from_PERSONS, output_spells_category, all.x = T, by="person_id")
  output_spells_category_enriched <- output_spells_category_enriched[entry_spell_category < date_of_birth + 60, entry_spell_category := date_of_birth]
  output_spells_category_enriched <- output_spells_category_enriched[, death_before_study_entry := fifelse(!is.na(date_of_death) & date_of_death < study_start, 1, 0)]
  output_spells_category_enriched <- output_spells_category_enriched[!is.na(entry_spell_category) & !is.na(exit_spell_category), no_observation_period_including_study_start := fifelse(study_start %between% list(entry_spell_category,exit_spell_category) & entry_spell_category < exit_spell_category, 0, 1)]
  output_spells_category_enriched <- output_spells_category_enriched[is.na(entry_spell_category) | is.na(exit_spell_category), no_observation_period_including_study_start := 1]
  output_spells_category_enriched <- output_spells_category_enriched[, insufficient_run_in := fifelse(entry_spell_category >= date_of_birth & entry_spell_category <= study_start - 365, 0, 1)]
  output_spells_category_enriched <- output_spells_category_enriched[, insufficient_run_in := min(fifelse(entry_spell_category == date_of_birth & year(date_of_birth) == 2019, 0 , insufficient_run_in)), by="person_id"]
  output_spells_category_enriched <- output_spells_category_enriched[, study_entry_date := study_start][, start_follow_up := start_lookback][, study_exit_date:= fifelse(no_observation_period_including_study_start == 0, pmin(date_of_death, exit_spell_category, study_end, na.rm = T), na_date)]
  output_spells_category_enriched <- output_spells_category_enriched[, no_observation_period_including_study_start := min(no_observation_period_including_study_start, na.rm = T), by="person_id"][, study_exit_date := min(study_exit_date, na.rm = T), by="person_id"]
  D3_exclusion_observation_periods_not_overlapping <- unique(output_spells_category_enriched[,.(person_id, study_entry_date, start_follow_up, study_exit_date, death_before_study_entry, insufficient_run_in, no_observation_period_including_study_start)])
  
  PERSONS_OP <- merge(D3_inclusion_from_PERSONS,
                      D3_exclusion_no_observation_period,
                      by="person_id",
                      all.x = T)
  PERSONS_OP2 <- merge(PERSONS_OP,
                       D3_exclusion_observation_periods_not_overlapping,
                       by="person_id",
                       all.x = T)
  
  coords<-c("sex_or_birth_date_missing", "birth_date_absurd", "no_observation_period", "insufficient_run_in","no_observation_period_including_study_start", "death_before_study_entry")
  selection_criteria <- PERSONS_OP2[, (coords) := replace(.SD, is.na(.SD), 0), .SDcols = coords]
  selection_criteria <- selection_criteria[!is.na(study_entry_date), ]
  
  tempname<-paste0("D3_selection_criteria",suffix[[subpop]])
  assign(tempname,selection_criteria)
  save(list=paste0("D3_selection_criteria",suffix[[subpop]]),file=paste0(dirtemp,"D3_selection_criteria",suffix[[subpop]],".RData"))
  rm(list=paste0("D3_selection_criteria",suffix[[subpop]]))
  
}
  
 rm(output_spells_category_enriched,D3_inclusion_from_PERSONS,D3_exclusion_observation_periods_not_overlapping)
rm(PERSONS_OP, PERSONS_OP2, na_date, coords)
 rm(D3_PERSONS, PERSONS_in_OP, output_spells_category,OBSERVATION_PERIODS, D3_exclusion_no_observation_period,selection_criteria)

