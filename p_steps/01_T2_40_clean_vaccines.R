# -----------------------------------------------------------------------
# Create Quality Check criteria for vaccination records

# input: CONCEPT
# output: D3_concepts_QC_criteria, QC_dose_derived, table_QC_dose_derived

load(paste0(dirconceptsets, "COVID_VACCINES.RData"))
Covid_vaccine <- COVID_VACCINES
rm(COVID_VACCINES)

# Recode, create variables and keep only useful ones

Covid_vaccine[, vx_record_date := ymd(vx_record_date)]
Covid_vaccine[, vx_manufacturer := tolower(as.character(vx_manufacturer))]
Covid_vaccine[!str_detect(vx_dose, "^[1-9]$"), vx_dose := 0]
Covid_vaccine[, vx_dose := as.integer(vx_dose)]
Covid_vaccine[, date_curated := fifelse(!is.na(date), date, vx_record_date)]
Covid_vaccine[, manufacturer_curated := vx_manufacturer]
Covid_vaccine[manufacturer_curated %not in% covid_vaccines_ConcePTION_CDM_vocabulary, manufacturer_curated := "unk"]

key_variables <- c("person_id", "date", "vx_record_date", "vx_dose",
                   "vx_manufacturer", "date_curated", "manufacturer_curated")
Covid_vaccine <- Covid_vaccine[, ..key_variables]


### Start with exclusion criteria
# Duplicated record
setkeyv(Covid_vaccine, key_variables)
Covid_vaccine[, duplicated_records := fifelse(rowidv(Covid_vaccine, key_variables) != 1, 1, 0)]
Covid_vaccine[, removed_row := duplicated_records]

# Wrong manufacturer
Covid_vaccine[removed_row == 0, manufacturer_not_in_study := fifelse(manufacturer_curated %not in% manufacturer_in_study, 1, 0)]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "manufacturer_not_in_study")]

# Missing date
Covid_vaccine[removed_row == 0, missing_date := fifelse(is.na(date_curated), 1, 0)]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "missing_date")]

# Date before start of vaccination in DAP region
Covid_vaccine[removed_row == 0, date_before_start_vax := fifelse(date_curated < start_COVID_vaccination_date, 1, 0)]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "date_before_start_vax")]

### Distance between doses and creation of imputed doses
# Order by person_id and date to get the first vaccination date
key_variables <- c("person_id", "date_curated")
setorderv(Covid_vaccine, c(key_variables))
Covid_vaccine[removed_row == 0, min_date_curated := date_curated[1], by = person_id]

# Calculate distance from first dose
Covid_vaccine[removed_row == 0, distance_doses := as.numeric(date_curated - min_date_curated)]

# Removed them if they are less than 19 days, exclude the first row of distance == 0 since it is the "correct" dose
Covid_vaccine[removed_row == 0, distance_btw_1_2_doses := fifelse(distance_doses >= 0 & distance_doses < 19, 1, 0)]
Covid_vaccine[removed_row == 0 & distance_doses == 0, flag := seq_len(.N), by = "person_id"]
Covid_vaccine[flag == 1, distance_btw_1_2_doses := 0]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "distance_btw_1_2_doses")]
Covid_vaccine[, c("flag", "min_date_curated", "distance_doses") := NULL]

# Order by person_id and date to get the second vaccination date
setorderv(Covid_vaccine, c(key_variables))
Covid_vaccine[removed_row == 0, min_date_curated := date_curated[2], by = person_id]

# Calculate distance from second dose
Covid_vaccine[removed_row == 0, distance_doses := as.numeric(date_curated - min_date_curated)]

# Removed them if they are less than 28 days but non negative, exclude the first row of distance == 0 since it is the "correct" dose
Covid_vaccine[removed_row == 0, distance_btw_2_3_doses := fifelse(distance_doses >= 0 & distance_doses < 28, 1, 0)]
Covid_vaccine[removed_row == 0 & distance_doses == 0, flag := seq_len(.N), by = "person_id"]
Covid_vaccine[flag == 1, distance_btw_2_3_doses := 0]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "distance_btw_2_3_doses")]
Covid_vaccine[, c("flag", "min_date_curated", "distance_doses") := NULL]

# Order by person_id and date to get the third vaccination date
setorderv(Covid_vaccine, c(key_variables))
Covid_vaccine[removed_row == 0, min_date_curated := date_curated[3], by = person_id]

# Calculate distance from third dose
Covid_vaccine[removed_row == 0, distance_doses := as.numeric(date_curated - min_date_curated)]

# Removed them if they are less than 28 days but non negative, exclude the first row of distance == 0 since it is the "correct" dose
Covid_vaccine[removed_row == 0, distance_btw_3_4_doses := fifelse(distance_doses >= 0 & distance_doses < 28, 1, 0)]
Covid_vaccine[removed_row == 0 & distance_doses == 0, flag := seq_len(.N), by = "person_id"]
Covid_vaccine[flag == 1, distance_btw_3_4_doses := 0]
Covid_vaccine[, removed_row := rowSums(.SD, na.rm = T), .SDcols=c("removed_row", "distance_btw_3_4_doses")]
Covid_vaccine[, c("flag", "min_date_curated", "distance_doses") := NULL]

# Imputation of doses in chronological order
Covid_vaccine <- Covid_vaccine[removed_row == 0, dose_curated := rowid(person_id)]
Covid_vaccine <- Covid_vaccine[removed_row == 0, imputed_dose := dose_curated != vx_dose]

# Dose after third
Covid_vaccine <- Covid_vaccine[removed_row == 0, dose_after_4 := fifelse(dose_curated > max_number_doses, 1, 0)]
Covid_vaccine[, removed_row := NULL]

# Clean dataset
D3_clean_vaccines <- Covid_vaccine[, .(person_id, date, vx_record_date, vx_dose, vx_manufacturer, date_curated,
                                       dose_curated, manufacturer_curated, imputed_dose, duplicated_records,
                                       manufacturer_not_in_study, missing_date, date_before_start_vax,
                                       distance_btw_1_2_doses, distance_btw_2_3_doses, distance_btw_3_4_doses,
                                       dose_after_4)]

# Imputation of missing values
for (i in names(D3_clean_vaccines)){
  D3_clean_vaccines[is.na(get(i)), (i) := 0]
}

# Saving
save(D3_clean_vaccines, file = paste0(dirtemp, "D3_clean_vaccines.RData"))
rm(Covid_vaccine, D3_clean_vaccines)
