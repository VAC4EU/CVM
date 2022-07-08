#-----------------------------------------------
# Create D4 descriptive tables for MIS/Myocard

# input: Flowchart_basic_exclusion_criteria, Flowchart_exclusion_criteria, D4_descriptive_dataset_ageband_studystart_MIS, D4_descriptive_dataset_age_studystart_MIS, D4_followup_fromstudystart_MIS, D4_descriptive_dataset_covariate_studystart_MIS, D4_population_d, D4_descriptive_dataset_ageband_studystart_c_MIS, D4_descriptive_dataset_age_studystart_c_MIS, D4_followup_fromstudystart_MIS_c, D4_descriptive_dataset_covid_studystart_c_MIS, D4_descriptive_dataset_covariate_covid_c_MIS, D3_Vaccin_cohort, QC_code_counts_in_study_population_OUTCOME_YEAR, RES_IR_monthly_MIS_b, RES_IR_monthly_MIS_c, RES_IR_monthly_MIS_d
# output: Attrition diagram 1, Attrition diagram 2, Cohort characteristics at start of study (1-1-2020), Cohort characteristics at first COVID-19 vaccination, Cohort characteristics at first occurrence of COVID-19 prior to vaccination (cohort c), COVID-19 vaccination by dose and time period between first and second dose (days), Code counts for narrow definitions (for each event) separately, Incidence of AESI (narrow) per 100,000 PY by calendar month in 2020, Incidence of each concept (narrow) per 100,000 PY prior to vaccination and COVID-19, Incidence of each concept (narrow) per 100,000 PY after COVID-19 and prior to vaccination, Incidence of each concept (narrow) per 100,000 PY after vaccination (BRAND)

##FUNCTIONS---------------------------------
`%not in%` <- negate(`%in%`)

create_empty_table_1a <- function() {
  n6 <- numeric(6)
  row_names_1 <- c("Start population", "A_sex_or_birth_date_missing", "C_no_observation_period",
                   "D_death_before_study_entry", "E_no_observation_period_including_study_start", "end population")
  if (length(intersect(names(flow_source_1a), c("Italy_ARS", "NL_PHARMO", "UK_CPRD", "ES_BIFAP"))) == 1) {
    ext = data.table(a = row_names_1, datasource = n6)
    setnames(ext, "datasource", vect_recode_manufacturer[thisdatasource])
    return(ext)
  }
  return(data.table(a = row_names_1, Italy_ARS = n6, NL_PHARMO = n6, UK_CPRD = n6, ES_BIFAP = n6))
}

create_empty_table_1b <- function() {
  n4 <- numeric(4)
  row_names_2 <- c("Source population present at January 1st 2020", "Subjects without valid birth date",
                   "Subjects without one year of look back at 1/1/2020", "Study population")
  if (length(intersect(names(flow_source_1a), c("Italy_ARS", "NL_PHARMO", "UK_CPRD", "ES_BIFAP"))) == 1) {
    ext = data.table(a = row_names_2, datasource = n4)
    setnames(ext, "datasource", vect_recode_manufacturer[thisdatasource])
    return(ext)
  }
  return(data.table(a = row_names_2, Italy_ARS = n4, NL_PHARMO = n4, UK_CPRD = n4, ES_BIFAP = n4))
}
##-------------------------


# Table1 --------------------------------------------------------------------
for (subpop in subpopulations_non_empty) {
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  if(this_datasource_has_subpopulations == T) dirD4tables <-paste0(thisdirexp,"D4 tables/")
  suppressWarnings(if (!file.exists(dirD4tables)) dir.create(file.path(dirD4tables)))
  
  if(this_datasource_has_subpopulations == T)   dummytables_MIS <- paste0(thisdirexp,"Dummy tables for report MIS-KD/")
  suppressWarnings(if (!file.exists(dummytables_MIS)) dir.create(file.path(dummytables_MIS)))
  
flow_source <- fread(paste0(thisdirexp, "Flowchart_basic_exclusion_criteria.csv"))
flow_study <- fread(paste0(thisdirexp, "Flowchart_exclusion_criteria.csv"))

vect_recode_manufacturer <- c(TEST = "Italy_ARS", ARS = "Italy_ARS", PHARMO = "NL_PHARMO",
                              CPRD = "UK_CPRD", BIFAP = "ES_BIFAP")

if ("datasource" %not in% names(flow_source)) {
  flow_source <- flow_source[ , Datasource := thisdatasource]
  flow_study <- flow_study[ , Datasource := thisdatasource]
} else {
  setnames(flow_source, "datasource", "Datasource")
  setnames(flow_source, "datasource", "Datasource")
}

flow_source <- flow_source[ , Datasource := vect_recode_manufacturer[Datasource]]
flow_study <- flow_study[ , Datasource := vect_recode_manufacturer[Datasource]]

flow_source_totals <- unique(data.table::copy(flow_source)[, TOTAL := sum(N), by = "Datasource"][, .(Datasource, TOTAL)])
flow_source_totals <- flow_source_totals[, a := "Start population"]
flow_source_totals <- data.table::dcast(flow_source_totals, a ~ Datasource, value.var = c("TOTAL"))

flow_source <- flow_source[, row_id := rowid(Datasource)]
flow_source <- data.table::melt(flow_source, id.vars = c("row_id", "Datasource", "N"),
                                measure.vars = c("A_sex_or_birth_date_missing", "B_birth_date_absurd",
                                                 "C_no_observation_period", "D_death_before_study_entry",
                                                 "E_no_observation_period_including_study_start"), variable.name = "a")

flow_source <- flow_source[, .SD[which.min(value)], by = c("row_id", "Datasource", "N")][value == 0, ][, .(Datasource, N, a)]
flow_source <- data.table::dcast(flow_source, a ~ Datasource, value.var = c("N"))
flow_source_1a <- flow_source[a != "B_birth_date_absurd", ]
flow_source_1a <- rbind(flow_source_1a, flow_source_totals)
empty_table_1 <- create_empty_table_1a()
flow_source_1a <- rbind(empty_table_1, flow_source_1a, fill = T)
flow_source_1a <- flow_source_1a[, lapply(.SD, max, na.rm = T), by = a]

max_col <- length(names(empty_table_1))
correct_col_names <- names(empty_table_1)[2:max_col]
correct_row_names <- empty_table_1$a[2:5]
table_1a <- flow_source_1a[a == "end population",
                           (correct_col_names) := flow_source_1a[a == "Start population", 2:max_col] -
                             colSums(flow_source_1a[a %in% correct_row_names, 2:max_col])]


flow_source_1b <- flow_source[a == "B_birth_date_absurd", ]
flow_source_1b <- flow_source_1b[a == "B_birth_date_absurd", a := "Subjects without valid birth date"]
flow_source_1b <- rbind(flow_source_1b, table_1a[a == "end population", ][, a := "Source population present at January 1st 2020"])
flow_study <- flow_study[, row_id := rowid(Datasource)]
flow_study <- data.table::melt(flow_study, id.vars = c("row_id", "Datasource", "N"),
                               measure.vars = c("A_insufficient_run_in"), variable.name = "a")
flow_study <- flow_study[value == 0, a := "Subjects without one year of look back at 1/1/2020"]
flow_study <- flow_study[value == 1, a := "Study population"][, .(Datasource, N, a)]
flow_study <- data.table::dcast(flow_study, a ~ Datasource, value.var = c("N"))
flow_source_1b <- rbind(flow_source_1b, flow_study)
empty_table_2 <- create_empty_table_1b()
flow_source_1b <- rbind(empty_table_2, flow_source_1b, fill = T)
table_1b <- flow_source_1b[, lapply(.SD, max, na.rm = T), by = a]

setnames(table_1a, "a", " ")
setnames(table_1b, "a", " ")

fwrite(table_1a, file = paste0(dummytables_MIS, "Attrition diagram 1.csv"))
fwrite(table_1b, file = paste0(dummytables_MIS, "Attrition diagram 2.csv"))



# Table2 ----------------------------------------------------------------------------------------------------------


ageband_studystart <- fread(paste0(dirD4tables, "D4_descriptive_dataset_ageband_studystart_MIS.csv"))

ageband_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                     BIFAP = "ES_BIFAP")[Datasource]]

ageband_studystart[, TOTAL := sum(.SD), by = Datasource, .SDcols = paste0("AgeCat_", Agebands_labels)]

total_pop <- ageband_studystart[, a := "Study population"][, Parameters := "N"][, .(a, Parameters, Datasource, TOTAL)]
total_pop <- dcast(total_pop, a + Parameters ~ Datasource, value.var = 'TOTAL')
col_to_keep <- intersect(c("a", "Parameters", "Italy_ARS", "NL_PHARMO",
                           "UK_CPRD", "ES_BIFAP", "Test"), names(total_pop))
total_pop <- total_pop[, ..col_to_keep]


age_studystart <- fread(paste0(dirD4tables, "D4_descriptive_dataset_age_studystart_MIS.csv"))

age_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                 BIFAP = "ES_BIFAP")[Datasource]]

followup_ss <- fread(paste0(dirD4tables, "D4_followup_fromstudystart_MIS.csv"))
followup_ss[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                              BIFAP = "ES_BIFAP")[Datasource]]
pt_total <- followup_ss[, a := "Person years of follow-up"][, Parameters := "PY"][, .(a, Parameters, Datasource,
                                                                                      Followup = Followup_males + Followup_females)]
pt_total <- dcast(pt_total, a + Parameters ~ Datasource, value.var = 'Followup')
pt_total <- pt_total[, ..col_to_keep]

age_start <- copy(age_studystart)[, a := "Age in years"][, Followup := NULL]
age_start <- melt(age_start, id.vars = c("a", "Datasource"),
                  measure.vars = c("Age_min", "Age_P25", "Age_P50", "Age_mean", "Age_p75", "Age_max"),
                  variable.name = "Parameters")
age_start <- dcast(age_start, a + Parameters  ~ Datasource, value.var = 'value')
age_start[, Parameters := c(Age_min = "Min", Age_P25 = "P25", Age_P50 = "P50", Age_mean = "Mean", Age_p75 = "P75",
                            Age_max = "Max")[Parameters]]
age_start <- age_start[, ..col_to_keep]


ageband_start <- ageband_studystart[, a := "Age in categories"][, TOTAL := NULL]
ageband_start <- melt(ageband_start, id.vars = c("a", "Datasource"),
                      measure.vars = paste0("AgeCat_", Agebands_labels),
                      variable.name = "Parameters")
ageband_start <- dcast(ageband_start, a + Parameters  ~ Datasource, value.var = 'value')
names(Agebands_labels) = paste0("AgeCat_", Agebands_labels)
ageband_start[, Parameters := Agebands_labels[Parameters]]


followup_studystart <- fread(paste0(dirD4tables, "D4_followup_fromstudystart_MIS.csv"))
followup_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                      BIFAP = "ES_BIFAP")[Datasource]]
followup_studystart <- followup_studystart[, a := "Person years across age categories"]
col_to_keep <- c("a", "Datasource", paste0("Followup_", Agebands_labels))
followup_start <- followup_studystart[, ..col_to_keep]

followup_start <- melt(followup_start, id.vars = c("a", "Datasource"),
                       measure.vars = paste0("Followup_", Agebands_labels),
                       variable.name = "Parameters")
followup_start <- dcast(followup_start, a + Parameters  ~ Datasource, value.var = 'value')
names(Agebands_labels) = paste0("Followup_", Agebands_labels)
followup_start[, Parameters := Agebands_labels[Parameters]]

sex_studystart <- fread(paste0(dirD4tables, "D4_followup_fromstudystart_MIS.csv"))
sex_studystart <- sex_studystart[, .(Datasource, Followup_males, Followup_females)]
sex_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS",
                                 PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                 BIFAP = "ES_BIFAP")[Datasource]]
sex_start <- sex_studystart[, a := "Person years across sex"]
sex_start <- melt(sex_start, id.vars = c("a", "Datasource"),
                  measure.vars = c("Followup_females", "Followup_males"),
                  variable.name = "Parameters")
sex_start <- dcast(sex_start, a + Parameters  ~ Datasource, value.var = 'value')
sex_start[, Parameters := c(Followup_males = "Male", Followup_males = "Female")[Parameters]]

risk_factors_studystart <- fread(paste0(dirD4tables, "D4_descriptive_dataset_covariate_studystart_MIS.csv"))
risk_factors_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                          BIFAP = "ES_BIFAP")[Datasource]]
risk_factors_start <- risk_factors_studystart[, a := "At risk population at January 1-2020"]
risk_factors_start <- melt(risk_factors_start, id.vars = c("a", "Datasource"),
                           measure.vars = c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes",
                                            "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"),
                           variable.name = "Parameters")
risk_factors_start <- dcast(risk_factors_start, a + Parameters  ~ Datasource, value.var = 'value')
risk_factors_start[, Parameters := c(CV = "Cardiovascular disease", Cancer = "Cancer", CLD = "Chronic lung disease",
                                     HIV = "HIV", CKD = "Chronic kidney disease", Diabetes = "Diabetes",
                                     Obesity = "Severe obesity", Sicklecell = "Sickle cell disease",
                                     immunosuppressants = "Use of immunosuppressants",
                                     any_risk_factors = "Any risk factors")[Parameters]]

table2 <- rbind(total_pop, pt_total, age_start, ageband_start, followup_start, sex_start, risk_factors_start)
daps <- intersect(c("Italy_ARS", "NL_PHARMO", "UK_CPRD", "ES_BIFAP", "Test"), names(table2))
daps_perc <- paste("perc", daps, sep="_")
col_order <- c(rbind(daps, daps_perc))
table2 <- table2[, (daps_perc) := character(nrow(table2))]
total_pop <- total_pop[, ..daps]
pt_total <- pt_total[, ..daps]
table2 <- table2[a %in% c("Age in categories", "At risk population at January 1-2020"),
                 (daps_perc) := round(.SD / as.numeric(total_pop) * 100, 1), .SDcols = daps]

table2 <- table2[a %in% c("Person years across age categories", "Person years across sex"),
                 (daps_perc) := round(.SD / as.numeric(pt_total) * 100, 1), .SDcols = daps]

table2 <- table2[a %in% c("Age in categories", "Person years across sex", "At risk population at January 1-2020", 
                          "Person years across age categories"), (daps_perc) := lapply(.SD, paste0, "%"), .SDcols = daps_perc]

empty_df <- table2[0,]
empty_df <- rbindlist(list(empty_df, list("", "", "N", "%")))

table2 <- rbindlist(list(empty_df, table2))

setcolorder(table2, c("a", "Parameters", col_order))

setnames(table2, "a", " ")

fwrite(table2, file = paste0(dummytables_MIS, "Cohort characteristics at start of study (1-1-2020).csv"))




# table2b ------------------------------------------------------------------

load(file = paste0(diroutput, "D4_population_d",suffix[[subpop]],".RData"))
D4_population_d<-get(paste0("D4_population_d",suffix[[subpop]]))
N_fup_pop <- D4_population_d[, .(person_id, sex, date_vax1, type_vax_1, fup_days, age_at_1_jan_2021, CV_at_date_vax_1,
                                  COVCANCER_at_date_vax_1, COVCOPD_at_date_vax_1, COVHIV_at_date_vax_1,
                                  COVCKD_at_date_vax_1, COVDIAB_at_date_vax_1, COVOBES_at_date_vax_1,
                                  COVSICKLE_at_date_vax_1, immunosuppressants_at_date_vax_1, at_risk_at_date_vax_1)]
setnames(N_fup_pop, c("date_vax1", "type_vax_1", "fup_days","age_at_1_jan_2021", "CV_at_date_vax_1",
                      "COVCANCER_at_date_vax_1", "COVCOPD_at_date_vax_1", "COVHIV_at_date_vax_1",
                      "COVCKD_at_date_vax_1", "COVDIAB_at_date_vax_1", "COVOBES_at_date_vax_1",
                      "COVSICKLE_at_date_vax_1", "immunosuppressants_at_date_vax_1", "at_risk_at_date_vax_1"),
         c("date_vax", "type_vax", "fup_vax","age_at_date_vax", "CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes",
           "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"))

# N_fup_pop <- melt(N_fup_pop, measure = list(c("date_vax1", "date_vax2"),
#                                             c("type_vax_1", "type_vax_2"),
#                                             c("fup_vax1", "fup_vax2"),
#                                             c("age_at_date_vax_1", "age_at_date_vax_2")),
#                   value.name = c("date_vax", "type_vax", "fup_vax", "age_at_date_vax"), na.rm = T)[, variable := NULL]
N_fup_pop[type_vax == "J&J", type_vax := "Janssen"]
N_fup_pop[, fup_vax := as.numeric(fup_vax) / 365.25]

vax_man <- unique(N_fup_pop[, type_vax])

vax_man_tot <- c("Pfizer", "Moderna", "AstraZeneca", "Janssen", "UKN")

vax_man <- intersect(vax_man_tot, vax_man)
vax_man_perc <- paste("perc", vax_man, sep = "_")
col_order <- c(rbind(vax_man, vax_man_perc))
cols_to_keep = c("a", "Parameters", col_order)

N_pop <- N_fup_pop[, .N, by = "type_vax"]
N_pop_by_vax <- setNames(copy(N_pop)$N, as.character(copy(N_pop)$type_vax))
N_pop_by_vax <- N_pop_by_vax[names(N_pop_by_vax) %in% vax_man]
N_pop_by_vax <- N_pop_by_vax[order(match(names(N_pop_by_vax), vax_man))]
total_pop <- N_pop[, sum(N)]
N_pop <- dcast(N_pop, . ~ type_vax, value.var = "N")[, . := NULL]
totals_man_d <- copy(N_pop)
N_pop <- N_pop[, Parameters := "N"][, a := "Study population"]
setnafill(N_pop, cols = c(vax_man), fill = 0)
N_pop <- N_pop[, (vax_man_perc) := round(.SD / as.numeric(total_pop) * 100, 1), .SDcols = vax_man]
N_pop <- N_pop[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
N_pop <- N_pop[, ..cols_to_keep]

fup_pop <- N_fup_pop[, sum(fup_vax), by = "type_vax"][, V1 := round(V1, 0)]
pt_total <- fup_pop[, sum(V1)]
fup_pop <- dcast(fup_pop, . ~ type_vax, value.var = "V1")[, . := NULL]
fup_pop <- fup_pop[, Parameters := "PY"][, a := "Person-years of follow-up"]
setnafill(fup_pop, cols = c(vax_man), fill = 0)
fup_pop <- fup_pop[, (vax_man_perc) := round(.SD / as.numeric(pt_total) * 100, 1), .SDcols = vax_man]
fup_pop <- fup_pop[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
fup_pop <- fup_pop[, ..cols_to_keep]


min_month <- N_fup_pop[, min(date_vax), by = "type_vax"][, V1 := month(V1)]
min_month <- dcast(min_month, . ~ type_vax, value.var = "V1")[, . := NULL]
min_month <- min_month[, Parameters := ""][, a := "Month of first vaccination"]
min_month <- min_month[, (vax_man_perc) := character(nrow(min_month))]
min_month <- min_month[, ..cols_to_keep]

year_month_pop <- N_fup_pop[, c("year", "month") := list(lubridate::year(date_vax), lubridate::month(date_vax))]
year_month_pop <- year_month_pop[, .N, by = c("type_vax", "year", "month")]
year_month_pop <- dcast(year_month_pop, year + month ~ type_vax, value.var = "N")
setorder(year_month_pop, month, year)
year_month_pop <- year_month_pop[, Parameters := "N"][, a := paste(year, month.name[month])]
setnafill(year_month_pop, cols = c(vax_man), fill = 0)
round_sum <- function(x) {round(x / sum(x) * 100, 1)}
year_month_pop <- year_month_pop[, (vax_man_perc) := lapply(.SD, round_sum), .SDcols = vax_man]
year_month_pop <- year_month_pop[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
year_month_pop <- year_month_pop[, ..cols_to_keep]

age_pop <- copy(N_fup_pop)[, .(Age_P25 = round(quantile(age_at_date_vax, probs = 0.25)),
                               Age_P50 = round(quantile(age_at_date_vax, probs = 0.50)),
                               Age_p75 = round(quantile(age_at_date_vax, probs = 0.75)),
                               Age_mean = round(mean(age_at_date_vax)),
                               Age_min = min(age_at_date_vax),
                               Age_max = max(age_at_date_vax)),
                           by = "type_vax"]
age_pop <- melt(age_pop, id.vars = "type_vax",
                measure.vars = c("Age_min", "Age_P25", "Age_P50", "Age_mean", "Age_p75", "Age_max"),
                variable.name = "child", value.name = "dob")
age_pop <- dcast(age_pop, child ~ type_vax, value.var = "dob")
age_pop <- age_pop[, Parameters := c("Min", "P25", "P50", "Mean", "P75", "Max")][, a := "Age in years"]
age_pop <- age_pop[, (vax_man_perc) := character(nrow(age_pop))]
age_pop <- age_pop[, ..cols_to_keep]

age_cat <- N_fup_pop[, age_at_date_vax := as.character(cut(age_at_date_vax, breaks = Agebands, labels = Agebands_labels))]


N_age_cat <- age_cat[, .N, by = c("type_vax", "age_at_date_vax")]
N_age_cat <- dcast(N_age_cat, age_at_date_vax ~ type_vax, value.var = "N")
setnames(N_age_cat, "age_at_date_vax", "Parameters")
N_age_cat <- N_age_cat[, a := "Age in categories"]
setnafill(N_age_cat, cols = c(vax_man), fill = 0)
N_age_cat <- N_age_cat[, (vax_man_perc) := lapply(.SD, round_sum), .SDcols = vax_man]
N_age_cat <- N_age_cat[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
N_age_cat <- N_age_cat[, ..cols_to_keep]

fup_age_cat <- age_cat[, sum(fup_vax), by = c("type_vax", "age_at_date_vax")][, V1 := round(V1, 0)]
fup_age_cat <- dcast(fup_age_cat, age_at_date_vax ~ type_vax, value.var = "V1")
setnames(fup_age_cat, "age_at_date_vax", "Parameters")
fup_age_cat <- fup_age_cat[, a := "Person years across age categories"]
setnafill(fup_age_cat, cols = c(vax_man), fill = 0)
fup_age_cat <- fup_age_cat[, (vax_man_perc) := lapply(.SD, round_sum), .SDcols = vax_man]
fup_age_cat <- fup_age_cat[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
fup_age_cat <- fup_age_cat[, ..cols_to_keep]

sex_pop <- copy(N_fup_pop)[, sum(fup_vax), by = c("type_vax", "sex")][, V1 := round(V1, 0)]
sex_pop[type_vax == "J&J", type_vax := "Janssen"]
sex_pop[, sex := as.character(sex)][, sex := c("0" = "Female", "1" = "Male")[sex]]
sex_pop <- dcast(sex_pop, sex ~ type_vax, value.var = "V1")
setnames(sex_pop, "sex", "Parameters")
sex_pop <- sex_pop[, a := "Person years across sex"]
setnafill(sex_pop, cols = c(vax_man), fill = 0)
sex_pop <- sex_pop[, (vax_man_perc) := lapply(.SD, round_sum), .SDcols = vax_man]
sex_pop <- sex_pop[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
sex_pop <- sex_pop[, ..cols_to_keep]

load(file = paste0(diroutput, "D4_population_d",suffix[[subpop]],".RData"))
D4_population_d<-get(paste0("D4_population_d",suffix[[subpop]]))

positive_before_vax <- D4_population_d[, .(person_id, history_covid, type_vax_1)]
positive_before_vax[type_vax_1 == "J&J", type_vax_1 := "Janssen"]


positive_before_vax <- positive_before_vax[, .N, by = c("type_vax_1", "history_covid")]

positive_before_vax <- dcast(positive_before_vax, history_covid ~ type_vax_1, value.var = "N")
positive_before_vax <- positive_before_vax[, Parameters := "N, %"]
positive_before_vax <- positive_before_vax[, a := "COVID-19 diagnosis/test in past year"]


setnafill(positive_before_vax, cols = c(vax_man), fill = 0)
positive_before_vax <- positive_before_vax[, (vax_man_perc) := lapply(.SD, round_sum), .SDcols = vax_man]
positive_before_vax <- positive_before_vax[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
positive_before_vax <- positive_before_vax[history_covid == 1, ]
positive_before_vax <- positive_before_vax[, history_covid := NULL]
positive_before_vax <- positive_before_vax[, ..cols_to_keep]

risk_factors <- copy(N_fup_pop)[, c("person_id", "type_vax", "CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes",
                                    "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")]
cols_chosen <- c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes", "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors")
risk_factors <- risk_factors[, lapply(.SD, sum, na.rm = T), by = "type_vax", .SDcols = cols_chosen]
risk_factors <- melt(risk_factors, id.vars = "type_vax",
                     measure.vars = cols_chosen,
                     variable.name = "Parameters", value.name = "dob")
risk_factors <- dcast(risk_factors, Parameters ~ type_vax, value.var = "dob")
risk_factors <- risk_factors[, a := "At risk population at date of vaccination"]
round_coverage <- function(x, y){
  round(x / as.numeric(y) * 100, 1)
}
risk_factors[, (vax_man_perc) := Map(round_coverage, .SD, N_pop_by_vax), .SDcols = vax_man]
risk_factors <- risk_factors[, (vax_man_perc) := lapply(.SD, paste0, "%"), .SDcols = vax_man_perc]
risk_factors <- risk_factors[, ..cols_to_keep]
risk_factors[, Parameters := c(CV = "Cardiovascular disease", Cancer = "Cancer", CLD = "Chronic lung disease",
                               HIV = "HIV", CKD = "Chronic kidney disease", Diabetes = "Diabetes",
                               Obesity = "Severe obesity", Sicklecell = "Sickle cell disease",
                               immunosuppressants = "Use of immunosuppressants",
                               any_risk_factors = "Any risk factors")[Parameters]]

table3_4_5_6 <- rbind(N_pop, fup_pop, min_month, year_month_pop, age_pop, N_age_cat, fup_age_cat, sex_pop,
                      positive_before_vax, risk_factors)
setnames(table3_4_5_6, "a", " ")

final_name_table3_4_5_6 <- c(TEST = "table 2", ARS = "table 2", PHARMO = "table 3",
                             CPRD = "table 3", BIFAP = "table 4")[[thisdatasource]]

vect_recode_manufacturer <- c(TEST = "Italy_ARS", ARS = "Italy_ARS", PHARMO = "Netherlands-PHARMO",
                              CPRD = "UK_CPRD", BIFAP = "ES_BIFAP")

empty_df <- table3_4_5_6[0,]
empty_df <- rbindlist(list(empty_df, as.list(c("", "", unlist(rep(c("N", "%"), length(vax_man)))))))

table3_4_5_6 <- rbindlist(list(empty_df, table3_4_5_6))

fwrite(table3_4_5_6, file = paste0(dummytables_MIS, final_name_table3_4_5_6,
                                   " Cohort characteristics at first COVID-19 vaccination ", 
                                   vect_recode_manufacturer[[thisdatasource]],".csv"))



# table5 ------------------------------------------------------------------


ageband_studystart_c <- fread(paste0(dirD4tables, "D4_descriptive_dataset_ageband_studystart_c_MIS.csv"))

ageband_studystart_c[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                       BIFAP = "ES_BIFAP")[Datasource]]

ageband_studystart_c[, TOTAL := sum(.SD), by = Datasource, .SDcols = paste0("AgeCat_", Agebands_labels)]

total_pop_c <- ageband_studystart_c[, a := "Study population"][, Parameters := "N"][, .(a, Parameters, Datasource, TOTAL)]
total_pop_c <- dcast(total_pop_c, a + Parameters ~ Datasource, value.var = 'TOTAL')
col_to_keep <- intersect(c("a", "Parameters", "Italy_ARS", "NL_PHARMO",
                           "UK_CPRD", "ES_BIFAP", "Test"), names(total_pop_c))
total_pop_c <- total_pop_c[, ..col_to_keep]


age_studystart_c <- fread(paste0(dirD4tables, "D4_descriptive_dataset_age_studystart_c_MIS.csv"))

age_studystart_c[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                   BIFAP = "ES_BIFAP")[Datasource]]

pt_total_c <- age_studystart_c[, a := "Person years of follow-up"][, Parameters := "PY"][, .(a, Parameters, Datasource, Followup)]
pt_total_c <- dcast(pt_total_c, a + Parameters ~ Datasource, value.var = 'Followup')
pt_total_c <- pt_total_c[, ..col_to_keep]

age_start_c <- copy(age_studystart_c)[, a := "Age in years"][, Followup := NULL]
age_start_c <- melt(age_start_c, id.vars = c("a", "Datasource"),
                    measure.vars = c("Age_min", "Age_P25", "Age_P50", "Age_mean", "Age_p75", "Age_max"),
                    variable.name = "Parameters")
age_start_c <- dcast(age_start_c, a + Parameters  ~ Datasource, value.var = 'value')
age_start_c[, Parameters := c(Age_min = "Min", Age_P25 = "P25", Age_P50 = "P50", Age_mean = "Mean", Age_p75 = "P75",
                              Age_max = "Max")[Parameters]]
age_start_c <- age_start_c[, ..col_to_keep]


ageband_start_c <- ageband_studystart_c[, a := "Age in categories"][, TOTAL := NULL]
ageband_start_c <- melt(ageband_start_c, id.vars = c("a", "Datasource"),
                        measure.vars = paste0("AgeCat_", Agebands_labels),
                        variable.name = "Parameters")
ageband_start_c <- dcast(ageband_start_c, a + Parameters  ~ Datasource, value.var = 'value')
names(Agebands_labels) = paste0("AgeCat_", Agebands_labels)
ageband_start_c[, Parameters := Agebands_labels[Parameters]]

followup_studystart <- fread(paste0(dirD4tables, "D4_followup_fromstudystart_MIS_c.csv"))
followup_studystart[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                      BIFAP = "ES_BIFAP")[Datasource]]
followup_studystart <- followup_studystart[, a := "Person years across age categories"]
setcolorder(followup_studystart, c("a", "Datasource"))

followup_start <- melt(followup_studystart, id.vars = c("a", "Datasource"),
                       measure.vars = colnames(followup_studystart)[colnames(followup_studystart) %not in% c("a", "Datasource")],
                       variable.name = "Parameters")
followup_start <- dcast(followup_start, a + Parameters  ~ Datasource, value.var = 'value')
names(Agebands_labels) = paste0("Followup_", Agebands_labels)
followup_start[, Parameters := Agebands_labels[Parameters]]

D4_descriptive_dataset_covid_studystart_c_MIS <- fread(paste0(dirD4tables, "D4_descriptive_dataset_covid_studystart_c_MIS.csv"))
D4_descriptive_dataset_covid_studystart_c_MIS[, Datasource := c(TEST = "Test", ARS = "Italy_ARS",
                                                                PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                                                BIFAP = "ES_BIFAP")[Datasource]]


covid_month <- D4_descriptive_dataset_covid_studystart_c_MIS[, a := "Month of first diagnosis"]
x<-colnames(covid_month)
cols_covid<-x[grepl("-", x)]

covid_month <- melt(covid_month, id.vars = c("a", "Datasource"),
                    measure.vars = cols_covid,
                    variable.name = "Parameters")
covid_month <- dcast(covid_month, a + Parameters  ~ Datasource, value.var = 'value')
#covid_month[, Parameters := c(Sex_male = "Male", Sex_female = "Female")[Parameters]]


risk_factors_studystart_c <- fread(paste0(dirD4tables, "D4_descriptive_dataset_covariate_covid_c_MIS.csv"))
risk_factors_studystart_c[, Datasource := c(TEST = "Test", ARS = "Italy_ARS", PHARMO = "NL_PHARMO", CPRD = "UK_CPRD",
                                            BIFAP = "ES_BIFAP")[Datasource]]
risk_factors_start_c <- risk_factors_studystart_c[, a := "At risk population at first covid diagnosis"]
risk_factors_start_c <- melt(risk_factors_start_c, id.vars = c("a", "Datasource"),
                             measure.vars = c("CV", "Cancer", "CLD", "HIV", "CKD", "Diabetes",
                                              "Obesity", "Sicklecell", "immunosuppressants", "any_risk_factors"),
                             variable.name = "Parameters")
risk_factors_start_c <- dcast(risk_factors_start_c, a + Parameters  ~ Datasource, value.var = 'value')
risk_factors_start_c[, Parameters := c(CV = "Cardiovascular disease", Cancer = "Cancer", CLD = "Chronic lung disease",
                                       HIV = "HIV", CKD = "Chronic kidney disease", Diabetes = "Diabetes",
                                       Obesity = "Severe obesity", Sicklecell = "Sickle cell disease",
                                       immunosuppressants = "Use of immunosuppressants",
                                       any_risk_factors = "Any risk factors")[Parameters]]

table5 <- rbind(total_pop_c, pt_total_c,covid_month, age_start_c, ageband_start_c, followup_start, risk_factors_start_c)
daps <- intersect(c("Italy_ARS", "NL_PHARMO", "UK_CPRD", "ES_BIFAP", "Test"), names(table5))
daps_perc <- paste("perc", daps, sep="_")
col_order <- c(rbind(daps, daps_perc))
table5 <- table5[, (daps_perc) := character(nrow(table5))]
total_pop_c <- total_pop_c[, ..daps]
pt_total_c <- pt_total_c[, ..daps]
table5 <- table5[a %in% c("Age in categories", "Month of first diagnosis", "At risk population at January 1-2020"),
                 (daps_perc) := round(.SD / as.numeric(total_pop_c) * 100, 1), .SDcols = daps]

table5 <- table5[a == "Person years across age categories",
                 (daps_perc) := round(.SD / as.numeric(pt_total_c) * 100, 1), .SDcols = daps]

table5 <- table5[a %in% c("Age in categories", "Person years across sex", "At risk population at January 1-2020", 
                          "Person years across age categories"), (daps_perc) := lapply(.SD, paste0, "%"), .SDcols = daps_perc]

empty_df <- table5[0,]
empty_df <- rbindlist(list(empty_df, list("", "", "N", "%")))

table5 <- rbindlist(list(empty_df, table5))

setcolorder(table5, c("a", "Parameters", col_order))

setnames(table5, "a", " ")

fwrite(table5, file = paste0(dummytables_MIS, "Cohort characteristics at first occurrence of COVID-19 prior to vaccination (cohort c).csv"))





# Table6 ----------------------------------------------------------------------------------------------------------

load(file = paste0(dirtemp, "D3_Vaccin_cohort",suffix[[subpop]],".RData"))
D3_Vaccin_coh<-get(paste0("D3_Vaccin_cohort",suffix[[subpop]]))
empty_table_7 <- data.table(a = character(0), Parameters = character(0), N = numeric(0))

names_vect <- c("Pfizer", "Moderna", "AstraZeneca", "Janssen", "UKN")
recode_rows <- paste(names_vect, "dose 1")
names(recode_rows) <- names_vect

vaccinated_persons <- D3_Vaccin_coh[, .(person_id, date_vax1, date_vax2, type_vax_1, type_vax_2)]
vaccinated_persons <- vaccinated_persons[date_vax1 <= study_end, ]
vaccinated_persons <- vaccinated_persons[date_vax2 > study_end, c("date_vax2", "type_vax_2") := NA]
vaccinated_persons <-vaccinated_persons[type_vax_1 == "J&J", type_vax_1 := "Janssen"]

names_vect_used <- names_vect[names_vect %in% vaccinated_persons[, type_vax_1] |
                                names_vect %in% vaccinated_persons[, type_vax_2]]

Totals_dose_1 <- vaccinated_persons[, .N, by = "type_vax_1"]
Totals_dose_1 <- Totals_dose_1[, index := 1]
Totals_dose_1 <- Totals_dose_1[, a := recode_rows[type_vax_1]][, Parameters := "Persons"]

Totals <- Totals_dose_1[, sum(N)]
Totals_df <- data.table::data.table(N = Totals, a = "Total population", Parameters = "Persons")
base_table_7 <- rbindlist(list(empty_table_7, Totals_df), use.names=TRUE)

df_with_second_doses  <- vaccinated_persons[!is.na(date_vax2), ]

dissimilar_doses <- df_with_second_doses[type_vax_1 != type_vax_2, ]
dissimilar_doses <- dissimilar_doses[, .N, by = "type_vax_1"]
dissimilar_doses <- merge(data.table(type_vax_1 = names_vect_used), dissimilar_doses, all.x = T, by = "type_vax_1")
dissimilar_doses <- dissimilar_doses[is.na(N), N := 0]
dissimilar_doses <- dissimilar_doses[, a := "Other vaccine dose 2"][, Parameters := "Persons"][, index := 3]

distance_doses <- df_with_second_doses[type_vax_1 == type_vax_2, ]
Totals_dose_2 <- distance_doses[, .N, by = "type_vax_1"]
recode_rows <- paste(names_vect, "dose 2")
names(recode_rows) <- names_vect
if ("Janssen" %not in% Totals_dose_2$type_vax_1) {
  Totals_dose_2 <- rbind(Totals_dose_2, list("Janssen", 0))
}
Totals_dose_2 <- Totals_dose_2[, a := recode_rows[type_vax_1]][, Parameters := "Persons"][, index := 2]

distance_doses <- distance_doses[, distance := correct_difftime(date_vax2 - 1, date_vax1)]
distance_doses <- distance_doses[, .(Min = min(distance), P25 = quantile(distance, 0.25),
                                     P50 = quantile(distance, 0.50), P75 = quantile(distance, 0.75),
                                     Max = max(distance)), by = "type_vax_1"]
cols_difftime <- sapply(distance_doses, is.difftime)
cols_difftime <- names(cols_difftime)[cols_difftime]
distance_doses <- distance_doses[, (cols_difftime) := lapply(.SD, function(x) round(as.numeric(x))), .SDcols = cols_difftime]
vect_measures <- c("Min", "P25", "P50", "P75", "Max")
distance_doses <- data.table::melt(distance_doses, id.vars = "type_vax_1",
                                   measure.vars = vect_measures, variable.name = "Parameters",
                                   value.name = "N")
recode_rows <- paste("Amongst persons with", names_vect, "dose 2 distance")
names(recode_rows) <- names_vect
rows_to_index <- seq(4, 8)
names(rows_to_index) <- vect_measures
distance_doses <- distance_doses[, a := recode_rows[type_vax_1]][, index := rows_to_index[Parameters]]

part2_table_7 <- rbindlist(list(Totals_dose_1, Totals_dose_2, dissimilar_doses, distance_doses), use.names=TRUE)
data.table::setorder(part2_table_7, type_vax_1, index)

dose_1_to_join <- Totals_dose_1[, .(type_vax_1, tot_type_1 = N)]

table_7 <- rbindlist(list(base_table_7, part2_table_7), fill = TRUE)
table_7 <- table_7[index == 1, Perc := N / Totals]
table_7 <- merge(table_7, dose_1_to_join, by = "type_vax_1")
table_7 <- table_7[data.table::between(index, 2, 3), Perc := N / tot_type_1]
table_7 <- table_7[, Perc := paste0(round(Perc * 100, 1), "%")]
table_7 <- table_7[Perc == "NA%", Perc := ""]

vect_recode_manufacturer <- c(TEST = "IT-ARS", ARS = "IT-ARS", PHARMO = "NL-PHARMO",
                              CPRD = "UK_CPRD", BIFAP = "ES_BIFAP")
correct_datasource <- vect_recode_manufacturer[thisdatasource]
table_7 <- table_7[, .(a, Parameters, N, Perc)]

empty_df <- table_7[0,]
empty_df <- rbindlist(list(empty_df, list("", "", "N", "%")))
table_7 <- rbindlist(list(empty_df, table_7))

setnames(table_7, c("a", "N", "Perc"), c("", correct_datasource, correct_datasource))

fwrite(table_7, file = paste0(dummytables_MIS, "COVID-19 vaccination by dose and time period between first and second dose (days).csv"))






# Table7 ----------------------------------------------------------------------------------------------------------

load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))

list_outcomes_observed<-get(paste0("list_outcomes_observed",suffix[[subpop]]))
list_outcomes_observed <- intersect(list_outcomes_observed, list_outcomes_MIS)
list_outcomes_observed <- list_outcomes_observed[str_detect(list_outcomes_observed, "narrow")]

table_7 <- data.table::data.table(meaning_of_first_event = character(), coding_system_of_code_first_event = character(),
                                  code_first_event = character(), count_n = character(), Event = character())

for (outcome in list_outcomes_observed) {
  for (year in c(2020, 2021)) {
    if (file.exists(paste0(thisdirexp, "QC_code_counts_in_study_population_", outcome, "_", year, ".csv"))) {
      temp_df <- data.table::fread(paste0(thisdirexp, "QC_code_counts_in_study_population_", outcome, "_", year, ".csv"),
                                   colClasses = list(character = "code_first_event"))
      event <- strsplit(outcome, "_")[[1]][1]
      temp_df <- temp_df[, Event := event]
      table_7 <- data.table::rbindlist(list(table_7, temp_df))
    }
  }
}

setnames(table_7, c("meaning_of_first_event", "coding_system_of_code_first_event", "code_first_event"),
         c("meaning_of_event", "Coding_system", "code"))

table_7 <- table_7[, count_n := as.integer(count_n)]
table_7 <- table_7[, .(sum = sum(count_n)), by = c("Event", "Coding_system", "code", "meaning_of_event")]

table_7 <- table_7[, Code := code][, DAP := thisdatasource]
table_7 <- table_7[, .(DAP, Event, Coding_system, Code, meaning_of_event, sum)]

fwrite(table_7, file = paste0(dummytables_MIS, "Code counts for narrow definitions (for each event) separately.csv"))



# Table8 ----------------------------------------------------------------------------------------------------------

load(paste0(thisdirexp,"RES_IR_monthly_MIS_b.RData"))
load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))

list_outcomes_observed <-get(paste0("list_outcomes_observed",suffix[[subpop]]))
RES_IR_monthly_MIS_b <-get(paste0("RES_IR_monthly_MIS_b"))
  
list_outcomes_observed <- intersect(list_outcomes_observed, list_outcomes_MIS)
list_outcomes_observed <- list_outcomes_observed[str_detect(list_outcomes_observed, "narrow")]

list_risk <- list_outcomes_observed
vect_recode_AESI <- list_outcomes_observed
names(vect_recode_AESI) <- c(as.character(seq_len(length(list_outcomes_observed))))

colA = paste0(list_risk, "_b")
colB = paste0("IR_", list_risk)
colC = paste0("lb_", list_risk)
colD = paste0("ub_", list_risk)

PT_monthly <- data.table::melt(RES_IR_monthly_MIS_b, measure = list(colA, colB, colC, colD),
                               variable.name = "AESI", value.name = c("Cases", "IR", "lb", "ub"), na.rm = F)

PT_monthly <- PT_monthly[, DAP := thisdatasource][ , AESI := vect_recode_AESI[AESI]]
PT_monthly <- PT_monthly[, .(DAP, sex, month, year, Ageband, AESI, Cases, IR, lb, ub)]


table_10 <- PT_monthly[year == 2020 & sex == "both_sexes" & Ageband == "all_birth_cohorts" & month != "all_months"
                       & !stringr::str_detect(AESI, "broad"), ]
table_10 <- table_10[, c("year", "sex", "Ageband") := NULL]

setcolorder(table_10, c("DAP", "AESI", "month", "Cases", "IR", "lb", "ub"))

setnames(table_10, c("month", "Cases", "IR", "lb", "ub"),
         c("Month in 2020", "Person years", "IR narrow", "LL narrow", "UL narrow"))

fwrite(table_10, file = paste0(dummytables_MIS, "Incidence of AESI (narrow) per 100,000 PY by calendar month in 2020.csv"))




# table_9 ----------------------------------------------------------------------------------------------------------

table_12 <- PT_monthly[year == "all_years" & sex != "both_sexes" & Ageband != "all_birth_cohorts" & month == "all_months"
                       & !stringr::str_detect(AESI, "broad"), ]
table_12 <- table_12[, c("year", "month") := NULL]

vect_recode_gender <- c("Male", "Female")
names(vect_recode_gender) <- c(1, 0)
table_12 <- table_12[ , sex := vect_recode_gender[sex]]

setcolorder(table_12, c("DAP", "AESI", "sex", "Ageband", "Cases", "IR", "lb", "ub"))

setorder(table_12, DAP, AESI, Ageband, sex)

setnames(table_12, c("Ageband", "sex", "Cases", "IR", "lb", "ub"),
         c("Age in 2020", "Sex", "Person years", "IR narrow", "LL narrow", "UL narrow"))

fwrite(table_12, file = paste0(dummytables_MIS, "Incidence of each concept (narrow) per 100,000 PY prior to vaccination and COVID-19.csv"))





# table_10 ----------------------------------------------------------------------------------------------------------

load(paste0(thisdirexp,"RES_IR_monthly_MIS_c.RData"))
load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))

RES_IR_monthly_MIS_c<-get(paste0("RES_IR_monthly_MIS_c"))
list_outcomes_observed<-get(paste0("list_outcomes_observed",suffix[[subpop]]))

list_outcomes_observed <- intersect(list_outcomes_observed, list_outcomes_MIS)
list_outcomes_observed <- list_outcomes_observed[str_detect(list_outcomes_observed, "narrow")]

list_risk <- list_outcomes_observed
vect_recode_AESI <- list_outcomes_observed
names(vect_recode_AESI) <- c(as.character(seq_len(length(list_outcomes_observed))))

colA = paste0(list_risk, "_b")
colB = paste0("IR_", list_risk)
colC = paste0("lb_", list_risk)
colD = paste0("ub_", list_risk)

PT_monthly <- data.table::melt(RES_IR_monthly_MIS_c, measure = list(colA, colB, colC, colD),
                               variable.name = "AESI", value.name = c("Cases", "IR", "lb", "ub"), na.rm = F)

PT_monthly <- PT_monthly[, DAP := thisdatasource][ , AESI := vect_recode_AESI[AESI]]
PT_monthly <- PT_monthly[, .(DAP, sex, month, year, Ageband, AESI, Cases, IR, lb, ub)]

table_12 <- PT_monthly[year == "all_years" & sex != "both_sexes" & Ageband != "all_birth_cohorts" & month == "all_months"
                       & !stringr::str_detect(AESI, "broad"), ]
table_12 <- table_12[, c("year", "month") := NULL]

vect_recode_gender <- c("Male", "Female")
names(vect_recode_gender) <- c(1, 0)
table_12 <- table_12[ , sex := vect_recode_gender[sex]]

setcolorder(table_12, c("DAP", "AESI", "sex", "Ageband", "Cases", "IR", "lb", "ub"))

setorder(table_12, DAP, AESI, Ageband, sex)

setnames(table_12, c("Ageband", "sex", "Cases", "IR", "lb", "ub"),
         c("Age in 2020", "Sex", "Person years", "IR narrow", "LL narrow", "UL narrow"))

fwrite(table_12, file = paste0(dummytables_MIS, "Incidence of each concept (narrow) per 100,000 PY after COVID-19 and prior to vaccination.csv"))





# table_11 ---------------------------------------------------------------------------

load(paste0(thisdirexp,"RES_IR_monthly_MIS_d.RData"))
load(paste0(dirtemp,"list_outcomes_observed",suffix[[subpop]],".RData"))

RES_IR_monthly_MIS_d<-get(paste0("RES_IR_monthly_MIS_d"))
list_outcomes_observed<-get(paste0("list_outcomes_observed",suffix[[subpop]]))

list_outcomes_observed <- intersect(list_outcomes_observed, list_outcomes_MIS)
list_outcomes_observed <- list_outcomes_observed[str_detect(list_outcomes_observed, "narrow")]

list_risk <- list_outcomes_observed
vect_recode_AESI <- list_outcomes_observed
names(vect_recode_AESI) <- c(as.character(seq_len(length(list_outcomes_observed))))

colA = paste0(list_risk, "_b")
colB = paste0("IR_", list_risk)
colC = paste0("lb_", list_risk)
colD = paste0("ub_", list_risk)

PT_monthly <- data.table::melt(RES_IR_monthly_MIS_d, measure = list(colA, colB, colC, colD),
                               variable.name = "AESI", value.name = c("Cases", "IR", "lb", "ub"), na.rm = F)

PT_monthly <- PT_monthly[, DAP := thisdatasource][ , AESI := vect_recode_AESI[AESI]]
PT_monthly <- PT_monthly[, .(DAP, sex, month, year, Ageband, type_vax_1, history_covid, AESI, Cases, IR, lb, ub)]

for (vax_m in recode(vax_man, Janssen = "J&J")) {
  table_12 <- PT_monthly[year == "all_years" & sex != "both_sexes" & Ageband != "all_birth_cohorts" & month == "all_months"
                         & !stringr::str_detect(AESI, "broad") & type_vax_1 == vax_m, ]
  table_12 <- table_12[, c("year", "month", "type_vax_1") := NULL]
  
  vect_recode_gender <- c("Male", "Female")
  names(vect_recode_gender) <- c(1, 0)
  table_12 <- table_12[ , sex := vect_recode_gender[sex]]
  
  setcolorder(table_12, c("DAP", "AESI", "sex", "Ageband", "history_covid", "Cases", "IR", "lb", "ub"))
  
  setorder(table_12, DAP, AESI, Ageband, sex)
  
  setnames(table_12, c("Ageband", "sex", "Cases", "IR", "lb", "ub"),
           c("Age in 2020", "Sex", "Person years", "IR narrow", "LL narrow", "UL narrow"))
  
  fwrite(table_12, file = paste0(dummytables_MIS, "Incidence of each concept (narrow) per 100,000 PY after vaccination (",
                                 vax_m,").csv"))
}

}
