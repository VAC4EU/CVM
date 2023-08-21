# -----------------------------------------------------
# CREATE EXCLUSION CRITERIA and CHECK CORRECT DATE OF BIRTH

# input: PERSONS, OBSERVATION_PERIODS
# output: D3_PERSONS, D3_events_DEATH

print('PRE-PROCESSING OF PERSONS')

# import input datasets
PERSONS <- data.table()
for (file in files_ConcePTION_CDM_tables[["PERSONS"]]) {
  temp <- fread(paste0(dirinput, file, ".csv"),
                colClasses = list(character = "person_id",
                                  integer = c("day_of_birth", "month_of_birth", "year_of_birth",
                                              "day_of_death", "month_of_death", "year_of_death")))
  PERSONS <- rbind(PERSONS, temp, fill = T)
  rm(temp)
}

OBSERVATION_PERIODS <- data.table()
for (file in files_ConcePTION_CDM_tables[["OBSERVATION_PERIODS"]]) {
  temp <- fread(paste0(dirinput, file, ".csv"), colClasses = list(character = "person_id"))
  OBSERVATION_PERIODS <- rbind(OBSERVATION_PERIODS, temp, fill = T)
  rm(temp)
}

OBSERVATION_PERIODS <- OBSERVATION_PERIODS[,`:=`(op_start_date = lubridate::ymd(op_start_date),
                                                 op_end_date = lubridate::ymd(op_end_date))]
D3_PERSONS <- PERSONS[, .(person_id, day_of_birth, month_of_birth, year_of_birth, sex_at_instance_creation,
                          day_of_death, month_of_death, year_of_death)]
rm(PERSONS)

# decide if pre-processing is needed for birth date (Present year, missing day or month)
PERSONS_date_missing <- D3_PERSONS[!is.na(year_of_birth) & (is.na(day_of_birth) | is.na(month_of_birth)),]
if(nrow(PERSONS_date_missing) != 0){
  
  print('SOME PERSONS HAVE DAY OR MONTH OF BIRTH MISSING')
  
  # Get the first start of observation periods for each person
  min_OBSERVATION_PERIODS <- OBSERVATION_PERIODS[!is.na(op_end_date), ]
  min_OBSERVATION_PERIODS <- min_OBSERVATION_PERIODS[, .(op_start_date = min(op_end_date)), by = person_id]
  
  # Merge persons with a missing date(as explained above) with observation periods
  CreateDateBirth <- merge(PERSONS_date_missing, min_OBSERVATION_PERIODS, all.x = T, by = "person_id")
  
  # start of first observation periods assumed date of birth
  CreateDateBirth[, assumed_year_birth := year(op_start_date)]
  CreateDateBirth[, assumed_month_birth := month(op_start_date)]
  CreateDateBirth[, assumed_day_birth := day(op_start_date)]

  # The order is IMPORTANT: first recode just day then in case month + day
  CreateDateBirth <- CreateDateBirth[year_of_birth == assumed_year_birth & !is.na(month_of_birth) & is.na(day_of_birth),
                                     birth_day_imputed := 1]
  CreateDateBirth <- CreateDateBirth[birth_day_imputed == 1 & month_of_birth == assumed_month_birth,
                                     day_of_birth := assumed_day_birth]
  CreateDateBirth <- CreateDateBirth[birth_day_imputed == 1 & month_of_birth != assumed_month_birth, day_of_birth := 15]
  
  CreateDateBirth <- CreateDateBirth[year_of_birth != assumed_year_birth | (year_of_birth == assumed_year_birth & is.na(month_of_birth)),
                                     c("birth_month_imputed", "birth_day_imputed") := list(1, 1)]
  CreateDateBirth <- CreateDateBirth[year_of_birth == assumed_year_birth & is.na(month_of_birth),
                                     c("month_of_birth", "day_of_birth") := list(assumed_month_birth, assumed_day_birth)]
  CreateDateBirth <- CreateDateBirth[year_of_birth != assumed_year_birth,
                                     c("month_of_birth", "day_of_birth") := list(6, 30)]
  
  # Clean the dataset
  CreateDateBirth <- CreateDateBirth[, c("op_start_date", "assumed_year_birth",
                                         "assumed_month_birth", "assumed_day_birth") := NULL]
  
  # Combine the persons which didn't need the correction for birth date to the one with the imputation
  D3_PERSONS <- rbind(D3_PERSONS[is.na(year_of_birth) | (!is.na(day_of_birth) & !is.na(month_of_birth)),],
                      CreateDateBirth, fill = T)
  
  rm(CreateDateBirth, min_OBSERVATION_PERIODS)
  
  print('DATE OF BIRTH IN PERSONS ADJUSTED')
} else {
  D3_PERSONS[, c("birth_month_imputed", "birth_day_imputed") := list(0, 0)]
}

# TODO check Codebook for specification of death
PERSONS_date_missing <- D3_PERSONS[!is.na(year_of_death) & (is.na(day_of_death) | is.na(month_of_death)),]

# decide if pre-processing is needed for death date
if(nrow(PERSONS_date_missing) != 0){
  
  print('SOME PERSONS HAVE DAY OR MONTH OF DEATH MISSING (WHILE YEAR EXISTS)')
  
  max_OBSERVATION_PERIODS <- OBSERVATION_PERIODS[!is.na(op_end_date), ]
  max_OBSERVATION_PERIODS <- max_OBSERVATION_PERIODS[, .(op_end_date = max(op_end_date)), by = person_id]
  
  D3_death_date <- merge(PERSONS_date_missing, max_OBSERVATION_PERIODS, by = "person_id")
  
  # end of last observation period is assumed date of death
  D3_death_date[, assumed_day_death := day(op_end_date)]
  D3_death_date[, assumed_month_death := month(op_end_date)]
  D3_death_date[, assumed_year_death := year(op_end_date)]

  D3_death_date <- D3_death_date[year_of_death == assumed_year_death & !is.na(month_of_death) & is.na(day_of_death),
                                 death_day_imputed := 1]
  D3_death_date <- D3_death_date[death_day_imputed == 1 & month_of_death == assumed_month_death,
                                 day_of_death := assumed_day_death]
  D3_death_date <- D3_death_date[death_day_imputed == 1 & month_of_death != assumed_month_death, day_of_death := 15]
  
  D3_death_date <- D3_death_date[year_of_death == assumed_year_death & is.na(month_of_death),
                                 c("death_month_imputed", "death_day_imputed") := list(1, 1)]
  D3_death_date <- D3_death_date[death_month_imputed == 1,
                                 c("month_of_death", "day_of_death") := list(assumed_month_death, assumed_day_death)]
  
  D3_death_date <- D3_death_date[, c("op_end_date", "assumed_year_death", "assumed_month_death", "assumed_day_death") := NULL]
  D3_PERSONS <- rbind(D3_PERSONS[is.na(year_of_death) | (!is.na(day_of_death) & !is.na(month_of_death)),],
                      D3_death_date, fill = T)
  
  rm(D3_death_date)
  
  print('DATE OF DEATH IN PERSONS ADJUSTED')
} else {
  D3_PERSONS[, c("death_month_imputed", "death_day_imputed") := list(0, 0)]
}

rm(PERSONS_date_missing)


# RETRIEVE FROM PERSONS ALL DEATHS AND SAVE
print('TRANSFORM in COMPLETED DATE FOR BIRTH and DEATH')

missing_date <- "9999-12-31"

# Convert birth date and clean the dataset
D3_PERSONS[, birth_date := fifelse(is.na(year_of_birth) | is.na(month_of_birth) | is.na(day_of_birth),
                                   missing_date,
                                   paste(year_of_birth, month_of_birth, day_of_birth, sep = "-"))]
D3_PERSONS[, birth_date := lubridate::ymd(birth_date)]
D3_PERSONS[, c("year_of_birth", "month_of_birth", "day_of_birth") := NULL]

# Convert death date and clean the dataset
D3_PERSONS[, flag := rowSums(is.na(.SD)), .SDcols = c("year_of_death", "month_of_death", "day_of_death")]
D3_PERSONS[flag %in% c(1, 2), death_date := missing_date]
D3_PERSONS[flag == 0, death_date := paste(year_of_death, month_of_death, day_of_death, sep = "-")]
D3_PERSONS[, death_date := lubridate::ymd(death_date)]
D3_PERSONS <- D3_PERSONS[, c("year_of_death", "month_of_death", "day_of_death", "flag") := NULL]

# Imputation of missing values
for (i in names(D3_PERSONS)[names(D3_PERSONS) != "death_date"]){
  D3_PERSONS[is.na(get(i)), (i) := 0]
}

# Retain only unique values
D3_PERSONS <- unique(D3_PERSONS)

# Create and save D3_events_DEATH
D3_events_DEATH <- D3_PERSONS[!is.na(death_date),.(person_id, death_date)][, date := death_date][, -"death_date"]
save(D3_events_DEATH,file = paste0(direvents,"D3_events_DEATH.RData"))
rm(D3_events_DEATH)

# Save D3_PERSONS
save(D3_PERSONS,file = paste0(dirtemp,"D3_PERSONS.RData"))

rm(D3_PERSONS, OBSERVATION_PERIODS)
