# -----------------------------------------------------------------------
# Create Quality Check criteria for vaccination records

# input: CONCEPT
# output: D3_concepts_QC_criteria, QC_dose_derived, table_QC_dose_derived

concepts <- import_concepts(dirtemp, vaccine__conceptssets)

# Recode, create variables and keep only useful ones

concepts <- concepts[, vx_record_date := ymd(vx_record_date)]
concepts <- concepts[, vx_manufacturer := as.character(vx_manufacturer)]

if (thisdatasource %in% c("TEST")) {
  concepts <- concepts[.(vx_manufacturer = c("MODERNA BIOTECH SPAIN S.L.", "PFIZER Srl", "ASTRAZENECA SpA", "JANSSEN CILAG SpA"),
                         to = c("Moderna", "Pfizer", "AstraZeneca", "J&J")), on = "vx_manufacturer",
                       vx_manufacturer := i.to]
}

if (thisdatasource %in% c("ARS")) {
  concepts <- concepts[.(vx_manufacturer = c("Moderna", "Pfizer", "AstraZeneca", "Janssen"),
                         to = c("Moderna", "Pfizer", "AstraZeneca", "J&J")), on = "vx_manufacturer",
                       vx_manufacturer := i.to]
}

concepts <- concepts[tolower(vx_manufacturer) == "moderna", vx_manufacturer := "Moderna"]
concepts <- concepts[tolower(vx_manufacturer) == "pfizer", vx_manufacturer := "Pfizer"]
concepts <- concepts[tolower(vx_manufacturer) == "astrazeneca", vx_manufacturer := "AstraZeneca"]
concepts <- concepts[tolower(vx_manufacturer) == "j&j" | tolower(vx_manufacturer) == "janssen", vx_manufacturer := "J&J"]

concepts <- concepts[vx_manufacturer %not in% c("Moderna", "Pfizer", "AstraZeneca", "J&J"), vx_manufacturer := "UKN"]
concepts <- concepts[, derived_date := fifelse(!is.na(date), date, vx_record_date)]

concepts <- concepts[, .(person_id, date, derived_date, vx_dose, vx_manufacturer)]

### Start with exclusion criteria

# Duplicated record
key_variables <- c("person_id", "vx_dose", "derived_date", "vx_manufacturer")
setkeyv(concepts, key_variables)
concepts <- concepts[, duplicated_records := fifelse(rowidv(concepts, key_variables) == 1, 0, 1)]
concepts <- concepts[, removed_row := duplicated_records]

# Missing date
concepts <- concepts[removed_row == 0, missing_date := as.numeric(is.na(derived_date))]
concepts[, removed_row := rowSums(.SD), .SDcols=c("removed_row", "missing_date")]

# Date before start of vaccination in DAP region
concepts <- concepts[removed_row == 0, date_before_start_vax := fifelse(derived_date >= start_COVID_vaccination_date, 0, 1)]
concepts[, removed_row := rowSums(.SD), .SDcols=c("removed_row", "date_before_start_vax")]

# Distance between doses and creation of imputed doses
key_variables <- c("person_id", "derived_date")
setorderv(concepts, c(key_variables))
concepts <- concepts[removed_row == 0, min_derived_date := derived_date[1], by = person_id]
concepts <- concepts[removed_row == 0, second_min_derived_date := derived_date[2], by = person_id]
concepts <- concepts[removed_row == 0, distance_doses := as.numeric(derived_date - min_derived_date)]
concepts <- concepts[removed_row == 0, second_distance_doses := as.numeric(derived_date - second_min_derived_date)]
concepts <- concepts[removed_row == 0, distance_btw_1_2_doses := fifelse(distance_doses >= 0 & distance_doses < 14, 1, 0)]
concepts <- concepts[removed_row == 0, distance_btw_2_3_doses := fifelse(second_distance_doses >= 0 & second_distance_doses < 90, 1, 0)]
concepts <- concepts[is.na(distance_btw_2_3_doses), distance_btw_2_3_doses := 0]
concepts <- concepts[, c("min_derived_date", "distance_doses", "second_distance_doses") := NULL]
key_variables <- "person_id"
concepts <- concepts[removed_row == 0, flag := rowidv(concepts[removed_row == 0, ], key_variables)]
concepts <- concepts[removed_row == 0 & flag == 1, distance_btw_1_2_doses := 0]
concepts[, removed_row := rowSums(.SD), .SDcols=c("removed_row", "distance_btw_1_2_doses")]
concepts <- concepts[removed_row == 0, flag := rowidv(concepts[removed_row == 0, ], key_variables)]
concepts <- concepts[removed_row == 0 & flag == 2, distance_btw_2_3_doses := 0]
concepts[, removed_row := rowSums(.SD), .SDcols=c("removed_row", "distance_btw_2_3_doses")]
concepts <- concepts[, flag := NULL]

concepts <- concepts[removed_row == 0, imputed_dose := rowid(person_id)]

# Create quality check table showing concordant dose number between original and imputed number
QC_dose_derived <- concepts[removed_row == 0,  wrong_dose := fifelse(vx_dose != imputed_dose, 1, 0)]

total_doses <- nrow(QC_dose_derived)
missing_dose_1 <- nrow(QC_dose_derived[is.na(vx_dose) & imputed_dose == 1, ])
missing_dose_2 <- nrow(QC_dose_derived[is.na(vx_dose) & imputed_dose == 2, ])
missing_dose_3 <- nrow(QC_dose_derived[is.na(vx_dose) & imputed_dose == 3, ])
dose_1_to_n <- nrow(QC_dose_derived[vx_dose == 1 & wrong_dose == 1, ])
dose_2_to_n <- nrow(QC_dose_derived[vx_dose == 2 & wrong_dose == 1, ])
dose_3_to_n <- nrow(QC_dose_derived[vx_dose == 3 & wrong_dose == 1, ])
dose_1_to_1 <- nrow(QC_dose_derived[vx_dose == 1 & wrong_dose == 0, ])
dose_2_to_2 <- nrow(QC_dose_derived[vx_dose == 2 & wrong_dose == 0, ])
dose_3_to_3 <- nrow(QC_dose_derived[vx_dose == 3 & wrong_dose == 0, ])
save(QC_dose_derived, file = paste0(dirtemp, "QC_dose_derived.RData"))

column_names <- c("", "Number of doses", "Missing first doses", "Missing second doses", "Missing third doses",
                  "Discordant first to n", "Discordant second to n", "Discordant third to n", "Concordant first doses",
                  "Concordant second doses", "Concordant third doses")
count_values <- c(total_doses, missing_dose_1, missing_dose_2, missing_dose_3, dose_1_to_n, dose_2_to_n, dose_3_to_n, 
                  dose_1_to_1, dose_2_to_2, dose_3_to_3)

table_QC_dose_derived <- data.table(a = column_names, b = c("N", count_values), c = c("%", sapply(count_values, prop_to_total)))
setnames(table_QC_dose_derived, c("a", "b", "c"), c("", thisdatasource, thisdatasource))
fwrite(table_QC_dose_derived, file = paste0(direxp, "table_QC_dose_derived.csv"))

# Dose after third
key_variables <- c("person_id", "derived_date")
setorderv(concepts, c(key_variables))
concepts <- concepts[removed_row == 0, dose_after_3 := fifelse(imputed_dose <= 3, 0, 1)]
concepts[, removed_row := rowSums(.SD), .SDcols=c("removed_row", "distance_btw_2_3_doses")]

# Dose after second
key_variables <- c("person_id", "derived_date")
setorderv(concepts, c(key_variables))
concepts <- concepts[removed_row == 0, dose_after_2 := fifelse(imputed_dose <= 2, 0, 1)]


setnames(concepts, "vx_dose", "old_dose")
setnames(concepts, "imputed_dose", "vx_dose")
setnames(concepts, "date", "old_date")
setnames(concepts, "derived_date", "date")
D3_concepts_QC_criteria <- concepts[, .(person_id, date, vx_dose, vx_manufacturer, duplicated_records, missing_date,
                                        date_before_start_vax, distance_btw_1_2_doses, distance_btw_2_3_doses,
                                        dose_after_3, dose_after_2)]

for (i in names(D3_concepts_QC_criteria)){
  D3_concepts_QC_criteria[is.na(get(i)), (i):=0]
}

save(D3_concepts_QC_criteria, file = paste0(dirtemp, "D3_concepts_QC_criteria.RData"))
rm(concepts, D3_concepts_QC_criteria, QC_dose_derived, table_QC_dose_derived)
