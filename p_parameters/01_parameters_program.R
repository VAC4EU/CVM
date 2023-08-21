# #set directory with input data

setwd("..")
setwd("..")
# dirbase <- getwd()
# dirinput <- paste0(dirbase,"/CDMInstances/PASS_COVIDVACCINES2205/")
# dirpregnancy <- paste0(dirbase,"/StudyScripts/pregnancy_20221017/g_output/")

dirinput <- paste0(thisdir,"/i_input/")
dirinput <- paste0(thisdir,"/i_input_subpop/")
dirpregnancy <- ""

# SCRI_list_variables <-  c("E_GOUT_COV", "C_MYOCARD_AESI", "C_PERICARD_AESI", "SO_OTITISEXT_COV", "C_VALVULAR_COV")

set_and_create_dir <- function(x) {
  x <- paste0(thisdir, x)
  dir.create(file.path(x), showWarnings = F)
  return(x)
}

# Sanitize dirpregnancy
if (dirpregnancy == "" | is.null(dirpregnancy) | is.na(dirpregnancy)) {
  skip_pregnancy = T
} else {
  dirpregnancy <- paste0(gsub("/$", "", dirpregnancy), "/")
  if (file.exists(paste0(dirpregnancy, "D3_pregnancy_final.RData"))) {
    skip_pregnancy = F
  } else {
    stop("there is no D3_pregnancy_final inside the folder specified in dirpregnancy")
  }
}

# set other directories
diroutput <- set_and_create_dir("/g_output/")
dirtemp <- set_and_create_dir("/g_intermediate/")
dirconceptsets <- set_and_create_dir("/g_intermediate/conceptset_datasets/")
diritemsets <- set_and_create_dir("/g_intermediate/itemset_datasets/")
dirpromptsets <- set_and_create_dir("/g_intermediate/promptset_datasets/")
direxp <- set_and_create_dir("/g_export/")
dirmacro <- set_and_create_dir("/p_macro/")
dirpargen <- set_and_create_dir("/g_parameters/")
direvents <- set_and_create_dir("/g_intermediate/events/")
dircomponents <- set_and_create_dir("/g_intermediate/components/")
dirTD <- set_and_create_dir("/g_intermediate/TD/")

rm(set_and_create_dir)

# load packages
read_library <- function(...) {
  x <- c(...)
  invisible(lapply(x, library, character.only = TRUE))
}

list.of.packages <- c("MASS", "haven", "tidyverse", "lubridate", "AdhereR", "stringr", "purrr", "readr", "dplyr",
                      "survival", "rmarkdown", "ggplot2", "data.table", "qpdf", "parallel", "readxl", "gtsummary",
                      "labelled", "huxtable", "metafor", "markdown", "R.utils", "SCCS", "RcppAlgos")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, require, character.only = T))

rm(read_library, new.packages, list.of.packages)

# load macros
source(paste0(dirmacro,"CreateConceptSetDatasets_v20.R"))
source(paste0(dirmacro,"CreateItemsetDatasets_v03.R"))
source(paste0(dirmacro,"MergeFilterAndCollapse_v5.R"))
source(paste0(dirmacro,"CreateSpells_v15.R"))
source(paste0(dirmacro,"CreateFlowChart.R"))
source(paste0(dirmacro,"CleanOutcomes.R"))
source(paste0(dirmacro,"CreateAgebandIntervals.R"))
source(paste0(dirmacro,"CreateTimeIntervals.R"))
source(paste0(dirmacro,"CheckAndPrepareDates.R"))
source(paste0(dirmacro,"CalculateSubtractionDenominator.R"))
source(paste0(dirmacro,"CalculateNumeratorAggregated.R"))
source(paste0(dirmacro,"SplitSpellsAgeBands.R"))
source(paste0(dirmacro,"CalculateNumeratorNotRecurrent.R"))
source(paste0(dirmacro,"SetToInteger.R"))
source(paste0(dirmacro,"CountPersonTimeV13.8.R"))
source(paste0(dirmacro,"df_to_list_of_list.R"))
source(paste0(dirmacro,"list_of_list_to_df.R"))
source(paste0(dirmacro,"dsr.R"))
source(paste0(dirmacro,"launch_step.R"))
source(paste0(dirmacro,"standardsccs2.R"))
source(paste0(dirmacro,"sensitivity_functions.R"))
source(paste0(dirmacro,"CompareListsOfCodes.R"))

#other parameters
# TODO remove this
firstjan2021 <- ymd(20210101)
#---------------------------------------
# understand which datasource the script is querying

study_start <- ymd(20190101)
start_lookback <- ymd(20180101)

CDM_SOURCE<- fread(paste0(dirinput,"CDM_SOURCE.csv"))

thisdatasource <- as.character(CDM_SOURCE[1,3])
instance_creation <- ymd(CDM_SOURCE[1,"date_creation"])
recommended_end_date <- ymd(CDM_SOURCE[1,"recommended_end_date"])
study_end <- min(instance_creation, recommended_end_date, na.rm = T)

rm(recommended_end_date, CDM_SOURCE)








start_COVID_vaccination_date <- fifelse(thisdatasource == 'CPRD', ymd(20201206), ymd(20201227))

start_COVID_diagnosis_date <- case_when((thisdatasource == 'TEST') ~ ymd(20200131),
                                        (thisdatasource == 'ARS') ~ ymd(20200131),
                                        (thisdatasource == 'PHARMO') ~ ymd(20200227),
                                        (thisdatasource == 'CPRD') ~ ymd(20200123),
                                        (thisdatasource == 'BIFAP') ~ ymd(20200131),
                                        (thisdatasource == 'SIDIAP') ~ ymd(20200131),
                                        (thisdatasource == 'UOSL') ~ ymd(20200226),
                                        TRUE ~ ymd(20200123))

###################################################################
# CREATE EMPTY FILES
###################################################################

files <- sub('\\.csv$', '', list.files(dirinput))

if (!any(str_detect(files,"^SURVEY_ID"))) {
  print("Creating empty SURVEY_ID since none were found")
  fwrite(data.table(person_id = character(0), survey_id = character(0), survey_date = character(0),
                    survey_meaning = character(0)),
         paste0(dirinput, "SURVEY_ID", ".csv"))
}

if (!any(str_detect(files,"^SURVEY_OBSERVATIONS"))) {
  print("Creating empty SURVEY_OBSERVATIONS since none were found")
  fwrite(data.table(person_id = character(0), so_date = character(0), so_source_table = character(0),
                    so_source_column = character(0), so_source_value = character(0), so_unit = character(0),
                    survey_id = character(0)),
         paste0(dirinput, "SURVEY_OBSERVATIONS", ".csv"))
}

if (!any(str_detect(files,"^MEDICINES"))) {
  print("Creating empty MEDICINES since none were found")
  fwrite(data.table(person_id = character(0), medicinal_product_id = integer(0),
                    medicinal_product_atc_code = character(0), date_dispensing = integer(0),
                    date_prescription = logical(0), disp_number_medicinal_product = numeric(0),
                    presc_quantity_per_day = logical(0), presc_quantity_unit = logical(0),
                    presc_duration_days = logical(0), product_lot_number = logical(0),
                    indication_code = logical(0), indication_code_vocabulary = logical(0),
                    meaning_of_drug_record = character(0), origin_of_drug_record = character(0),
                    prescriber_speciality = logical(0), prescriber_speciality_vocabulary = logical(0),
                    visit_occurrence_id = character(0)),
         paste0(dirinput, "MEDICINES_FED", ".csv"))
}

rm(files)

#############################################
#SAVE METADATA TO direxp
#############################################

file.copy(paste0(dirinput,'/METADATA.csv'), direxp, overwrite = T)
file.copy(paste0(dirinput,'/CDM_SOURCE.csv'), direxp, overwrite = T)
file.copy(paste0(dirinput,'/INSTANCE.csv'), direxp, overwrite = T)

#############################################
#SAVE to_run.R TO direxp
#############################################

file.copy(paste0(thisdir,'/to_run.R'), direxp, overwrite = T)

#study_years_datasource

# study_years <- c("2019", "2020", "2021")

# TODO should add 2018?
ComponentAnalysisYears <- c("2018", "2019")

days <- ifelse(thisdatasource %in% c("ARS","TEST"), 180, 1)

#############################################
#RECODING FOR OUTPUT TABLES
#############################################

vect_recode_dap <- c(TEST = "Italy_ARS",
                     ARS = "Italy_ARS",
                     PHARMO = "NL_PHARMO",
                     CPRD = "UK_CPRD",
                     BIFAP = "ES_BIFAP")

vect_recode_dap <- data.table(ori = names(vect_recode_dap), new = vect_recode_dap)

export_dap_name <- as.character(as.data.table(thisdatasource)[vect_recode_dap,
                                                              on = .(thisdatasource = ori),
                                                              "thisdatasource" := .(i.new)])

# TODO remove?
vect_new_severity <- c("covid_severity_1_plus", "covid_severity_2_plus", "covid_severity_3_plus",
                       "covid_severity_4_plus", "covid_severity_1", "covid_severity_2", "covid_severity_3",
                       "covid_severity_4")

### Set gtsummary theme
suppressMessages(gtsummary::theme_gtsummary_language("en", big.mark = ""))

split_by_default <- c("sex", "dose", "type_vax")

#############################################
#FUNCTION TO COMPUTE AGE
#############################################

# TODO check agebands
Agebands = c(-1, 4, 11, 17, 24, 29, 39, 49, 59, 69, 79, Inf)
Agebands_countpersontime = c(0, 4, 11, 17, 24, 29, 39, 49, 59, 69, 79)
Agebands_labels = c("0-4","5-11","12-17","18-24","25-29", "30-39", "40-49","50-59","60-69", "70-79","80+")
names(Agebands_countpersontime) <- Agebands_labels

age_fast = function(from, to) {
  from_lt = as.POSIXlt(from)
  to_lt = as.POSIXlt(to)
  
  age = to_lt$year - from_lt$year
  
  ifelse(to_lt$mon < from_lt$mon |
           (to_lt$mon == from_lt$mon & to_lt$mday < from_lt$mday),
         age - 1, age)
}

`%not in%` = Negate(`%in%`)

substrRight <- function(x, n){
  char_x <- nchar(x)
  substr(x, char_x - n + 1, char_x)
}

find_last_monday <- function(tmp_date, monday_week) {
  
  tmp_date <- as.Date(lubridate::ymd(tmp_date))
  
  while (tmp_date %not in% monday_week) {
    tmp_date <- tmp_date - 1
  }
  return(tmp_date)
}

find_first_monday_year <- function(tmp_date, monday_week) {
  
  tmp_date <- as.Date(lubridate::ymd(tmp_date))
  
  while (tmp_date %not in% monday_week) {
    tmp_date <- tmp_date + 1
  }
  return(tmp_date)
}

correct_difftime <- function(t1, t2, t_period = "days") {
  return(difftime(t1, t2, units = t_period) + 1)
}

calc_precise_week <- function(time_diff) {
  # correction in case a person exit the same date it enter
  time_diff <- fifelse(time_diff == 1, time_diff + 1, time_diff)
  weeks_frac <- time_length(time_diff - 1, "week")
  fifelse(weeks_frac%%1==0, weeks_frac, floor(weeks_frac) + 1)
}

exactPoiCI <- function (df, X, PT, conf.level = 0.95) {
  alpha <- 1 - conf.level
  IR <- df[, get(X)]
  upper <- df[, 0.5 * qchisq((1-(alpha/2)), 2*(get(X)+1))]
  lower <- df[, 0.5 * qchisq(alpha/2, 2*get(X))]
  temp_list <- lapply(list(IR, lower, upper), `/`, df[, get(PT)/365.25])
  temp_list <- lapply(temp_list, `*`, 100000)
  temp_list <- lapply(temp_list, function(x) {fifelse(x == Inf, 0, x)})
  return(lapply(temp_list, round, 2))
}

divide_period_per_event <- function(.data, .x, start_period, end_period) {
  
  # Set COVID19 to 0 for all period ending before the event and person without the event
  temp_0 <- .data[get(end_period) < get(.x) | is.na(get(.x)), ][, COVID19 := 0]
  
  # Set COVID19 to 1 for all period starting after the event
  temp_3 <- .data[get(start_period) > get(.x), ][, COVID19 := 1]
  
  # Get all periods containing events
  temp_1 <- .data[data.table::between(get(.x), get(start_period), get(end_period)), ]
  
  # Create a copy of the periods with events and set events to 0.
  # Then input end of the period as the day before the event and remove periods with eent in the first day
  temp_2 <- copy(temp_1)[, COVID19 := 0][, (end_period) := get(.x) - 1][get(end_period) >= get(start_period), ]
  
  # Set events to 1 for periods with one of them
  # Then input start of the period as the day of the event
  temp_1 <- temp_1[, COVID19 := 1][, (start_period) := get(.x)]
  
  # Combine all dataset and return
  return(rbind(temp_0, temp_1, temp_2, temp_3))
}

generate_formulas <- function(left, right_string) {
  return(formula(paste(left, paste0("'", right_string, "'"), sep = " ~ ")))
}

gt_dichotomous <- function(.data, col_to_print, values_to_print, str_to_print, other = F, remove_na = F) {
  
  values_to_print <- as.list(values_to_print)
  tbl_list <- list()
  
  for (level in values_to_print) {
    
    if (other) {
      other_levels <- setdiff(levels(.data[, get(col_to_print)]), basic_vx_manufacturer)
      tbl_out <- copy(.data)[get(col_to_print) %in% other_levels, (col_to_print) := "other"]
      level <- "other"
    } else {
      tbl_out <- copy(.data)
    }
    
    tbl_out <- tbl_out[is.na(get(col_to_print)) | get(col_to_print) != level, (col_to_print) := paste("not", level)]
    
    tbl_out <- tbl_out %>%
      tbl_summary(label = generate_formulas(col_to_print, sprintf(str_to_print, level)),
                  value = generate_formulas(col_to_print, eval(level)),
                  by = "DAP",
                  type = generate_formulas(col_to_print, "dichotomous"),
                  include = eval(col_to_print), digits = everything() ~ c(0, 2)) %>%
      modify_header(all_stat_cols(T) ~ header_string) %>%
      modify_footnote(all_stat_cols(FALSE) ~ NA)
    
    tbl_list <- append(tbl_list, list(tbl_out))
  }
  
  return(tbl_stack(tbl_list))
}

tbl_PT_IR_dichotomous <- function(.data, .string) {
  
  # Table with count and persontime
  first_tbl <- copy(.data)[statistic %in% c("counts", "persontime"), ]
  first_tbl[, statistic := factor(statistic, levels = c("counts", "persontime"))]

  first_tbl <- first_tbl %>%
    tbl_custom_summary(label = var ~ .string,
                       by = "statistic",
                       value = list(var = current_var),
                       stat_fns = everything() ~ function(data, ...) dplyr::tibble(value = data$value),
                       statistic = everything() ~ "{value}",
                       type = var ~ "dichotomous",
                       include = var) %>%
    modify_header(label = "**{current_var}**", stat_1 = "**counts**", stat_2 = "**persontime**") %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  # Table with IR
  second_tbl <- copy(.data)[statistic %in% c("IR", "lb", "ub"), ]
  second_tbl[, value := as.numeric(value)]
  second_tbl[, statistic := factor(statistic, levels = c("IR", "lb", "ub"))]
  
  second_tbl <- second_tbl %>%
    tbl_custom_summary(label = var ~ .string,
                       by = "statistic",
                       value = list(var = current_var),
                       stat_fns = everything() ~ function(data, ...) dplyr::tibble(value = data$value),
                       statistic = everything() ~ "{value}",
                       type = var ~ "dichotomous",
                       digits = var ~ 2,
                       include = var) %>%
    modify_header(label = "**{current_var}**", stat_1 = "**IR**", stat_2 = "**lb**", stat_3 = "**ub**") %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  return(tbl_merge(list(first_tbl, second_tbl), tab_spanner = F))
}

tbl_PT_IR_categorical <- function(.data, str_to_print, col_to_print) {
  
  # Table with count and persontime
  first_tbl <- copy(.data)[statistic %in% c("counts", "persontime"), ]
  first_tbl[, statistic := factor(statistic, levels = c("counts", "persontime"))]
  
  first_tbl <- first_tbl %>%
    tbl_custom_summary(label = generate_formulas(col_to_print, str_to_print),
                       by = "statistic",
                       stat_fns = everything() ~ function(data, ...) dplyr::tibble(value = data$value),
                       statistic = everything() ~ "{value}",
                       type = generate_formulas(col_to_print, "categorical"),
                       digits = formula(paste(col_to_print, 0, sep = " ~ ")),
                       include = eval(col_to_print)) %>%
    modify_header(label = "**{current_var}**", stat_1 = "**counts**", stat_2 = "**persontime**") %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  # Table with IR
  second_tbl <- copy(.data)[statistic %in% c("IR", "lb", "ub"), ]
  second_tbl[, statistic := factor(statistic, levels = c("IR", "lb", "ub"))]
  
  second_tbl <- second_tbl %>%
    tbl_custom_summary(label = generate_formulas(col_to_print, str_to_print),
                       by = "statistic",
                       stat_fns = everything() ~ function(data, ...) dplyr::tibble(value = data$value),
                       statistic = everything() ~ "{value}",
                       type = generate_formulas(col_to_print, "categorical"),
                       digits = formula(paste(col_to_print, 2, sep = " ~ ")),
                       include = eval(col_to_print)) %>%
    modify_header(label = "**{current_var}**", stat_1 = "**IR**", stat_2 = "**lb**", stat_3 = "**ub**") %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  return(tbl_merge(list(first_tbl, second_tbl), tab_spanner = F))
}

create_table_characteristic_population <- function(study_pop, persontime = NULL, covariates, at_time,
                                                   agebands = Agebands_countpersontime) {
  
  if (at_time == "baseline") {
    col_by = "DAP"
    header_string <- "**{paste0(thisdatasource, suffix[[subpop]])}**"
    
    # Create a column for the DAP name
    study_pop[, DAP := paste0(thisdatasource, suffix[[subpop]])]
    # Create a fixed variable to 1 and add the column for the DAP
    tot_pop <- study_pop[, .(person_id, total = 1, DAP)]
    # Select only the persontime, divide it by 365.25 to get PT in years and add the column for the DAP
    tot_PT <- study_pop[, .(person_id, Persontime = correct_difftime(study_exit_date, spell_start_date) / 365.25, DAP)]
    tot_PT <- tot_PT[, .(Persontime = round(sum(Persontime), 0)), by = col_by]
    # Calculate age at start follow-up
    study_pop[, age := age_fast(date_of_birth, start_followup_study)]
    # Select only variables of interest and add the DAP name
    pop_age_sex <- study_pop[, .(person_id, age, sex, DAP)]
    # Keep only time which we need and then remove the type of date
    covariates <- covariates[type_of_date == "baseline", ][, type_of_date := NULL]
    # Create study_pop for covariates
    study_pop_for_covariates <- study_pop[, .(person_id, sex, DAP)]
  } else if (at_time == "vax1") {
    col_by = "type_vax_1"
    header_string <- "**at first {level}**"
    
    # Select only persons with at least 1 vaccination
    study_pop <- study_pop[!is.na(date_vax_1), ]
    # Create a fixed variable to 1 and retain the column for the manufacturer
    tot_pop <- study_pop[, .(person_id, type_vax_1, total = 1)]
    tot_pop[, type_vax_1 := factor(type_vax_1, levels = manufacturer_in_study)]
    # Select only the persontime, divide it by 365.25 to get PT in years and retain the column for the manufacturer
    tot_PT <- study_pop[, .(Persontime = correct_difftime(study_exit_date, date_vax_1) / 365.25, type_vax_1)]
    tot_PT <- tot_PT[, .(Persontime = round(sum(Persontime), 0)), by = col_by]
    tot_PT[, type_vax_1 := factor(type_vax_1, levels = manufacturer_in_study)]
    # Calculate age at start follow-up
    study_pop[, age := age_fast(date_of_birth, date_vax_1)]
    # Select only variables of interest and add the manufacturer
    pop_age_sex <- study_pop[, .(person_id, age, sex, type_vax_1)]
    pop_age_sex[, type_vax_1 := factor(type_vax_1, levels = manufacturer_in_study)]
    # Keep only time which we need and then remove the type of date
    covariates <- covariates[type_of_date == "vax1", ][, type_of_date := NULL]
    # Create study_pop for covariates
    study_pop_for_covariates <- study_pop[, .(person_id, sex, type_vax_1)]
    study_pop_for_covariates[, type_vax_1 := factor(type_vax_1, levels = manufacturer_in_study)]
  } 
  
  ### Total population
  # Create the table which contains the total population
  tot_pop <- tot_pop %>%
    tbl_summary(label = list(total ~ "Study population"), by = all_of(col_by), include = total, percent = "row") %>%
    modify_header(all_stat_cols(T) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  defaultW <- getOption("warn") 
  options(warn = -1)
  ### Total PT
  # Create the table which contains the total population
  tot_PT <- tot_PT %>%
    tbl_summary(label = list(Persontime ~ "follow-up (years)"), by = all_of(col_by),
                type = Persontime ~ "continuous", digits = Persontime ~ 0,
                statistic = Persontime ~ "{min} (PY)") %>%
    modify_header(all_stat_cols(FALSE) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  ### Age and sex
  # Calculate age at start follow-up
  pop_age_sex[, ageband := Agebands_labels[findInterval(pop_age_sex[, age], agebands,
                                                        rightmost.closed = T, left.open = T)]]
  
  # Transform ageband to a factor to set the ordering the agebands
  pop_age_sex[, ageband := factor(ageband, levels = names(agebands))]
  
  # Recode sex
  pop_age_sex[, sex := fcase(sex == "F", "Female", sex == "M", "Male", sex == "O", "Other")]
  
  # Table which age basic statistics
  age_sex_characteristics <- pop_age_sex %>%
    tbl_summary(label = list(age ~ "Age in years",
                             ageband ~ "Age in categories",
                             sex ~ "Persons"),
                by = all_of(col_by),
                type = age ~ "continuous2",
                statistic = age ~ c("{min}", "{p25}", "{median}", "{mean}", "{p75}", "{max}"),
                include = c(age, ageband, sex)) %>%
    modify_header(all_stat_cols(FALSE) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  options(warn = defaultW)
  
  ### Covariates
  # covariates[ , PREGNANCY := round(runif(nrow(covariates)))]
  
  # Create vector with risk factors and medicines defined as variables (PREGNANCY should be counted as special case)
  risk_factors <- check_columns_exist(covariates, setdiff(c(COV_variables, "COVID"), "PREGNANCY"))
  medicines <- check_columns_exist(covariates, DP_variables)
  pregnancy <- check_columns_exist(covariates, "PREGNANCY")
  
  # Merge study population and covariates to add sex which we will need for PREGNANCY
  pop_covariates <- merge(study_pop_for_covariates, covariates, by = "person_id")
  
  # Change column names to final names in tables
  
  # Get the column to recode and the vector with the names of the covariate
  cols <- c(risk_factors, medicines, pregnancy)
  covariate_names <- c(str_match(setdiff(cols, pregnancy), "_(.*?)_")[, 2], pregnancy)
  covariate_names[is.na(covariate_names)] <- "COVID"
  
  # Recode 1 to name of covariate
  pop_covariates <- pop_covariates[ , (cols) := Map(function(x, single_cov) ifelse(x == 1, single_cov, F),
                                                    .SD, covariate_names),
                                    .SDcols = cols]
  
  param <- setdiff(covariate_names, pregnancy)
  names(param) <- setdiff(cols, pregnancy)
  
  # Table which covariates counts
  covariate_characteristics <- pop_covariates %>%
    tbl_summary(label = lapply(Map(paste, names(param), paste0('"', param, '"'), sep = " ~ ", USE.NAMES = F), formula),
                by = all_of(col_by),
                type = everything() ~ "dichotomous",
                value = as.list(param),
                # statistic = everything() ~ "{n}",
                include = c(all_of(risk_factors), all_of(medicines))) %>%
    modify_header(all_stat_cols(FALSE) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  # Table for pregnancy. That's special because we need the percentage in relation to only females
  pregnancy_characteristics <- pop_covariates[sex == "F", ] %>%
    tbl_summary(label = all_of(pregnancy) ~ "PREGNANCY",
                by = all_of(col_by),
                type = all_of(pregnancy) ~ "dichotomous",
                value = list(all_of(pregnancy) ~ "PREGNANCY"),
                include = all_of(pregnancy)) %>%
    modify_header(all_stat_cols(FALSE) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  
  return(tbl_stack(list(tot_pop, tot_PT, age_sex_characteristics, pregnancy_characteristics,
                        covariate_characteristics)))
}

save_tbl_summary <- function(dir_export, tbl_obj, tbl_name, additional_folder = NULL) {
  
  suppressWarnings(if(!file.exists(paste0(dir_export, "csv 3-4-5-6/", additional_folder))) dir.create(file.path(paste0(dir_export, "csv 3-4-5-6/", additional_folder))))
  suppressWarnings(if(!file.exists(paste0(dir_export, "rtf 3-4-5-6/", additional_folder))) dir.create(file.path(paste0(dir_export, "rtf 3-4-5-6/", additional_folder))))
  suppressWarnings(if(!file.exists(paste0(dir_export, "html 3-4-5-6/", additional_folder))) dir.create(file.path(paste0(dir_export, "html 3-4-5-6/", additional_folder))))
  suppressWarnings(if(!file.exists(paste0(dir_export, "rdata 3-4-5-6/", additional_folder))) dir.create(file.path(paste0(dir_export, "rdata 3-4-5-6/", additional_folder))))
  
  tbl_obj$inputs <- NULL
  
  final_csv <- tbl_obj %>%
    gtsummary::as_tibble()
  write.csv(final_csv, paste0(paste0(dir_export, "csv 3-4-5-6/", additional_folder), tbl_name, ".csv"))
  
  tbl_obj <- tbl_obj %>%
    as_hux_table()
  
  save(tbl_obj, file = paste0(paste0(dir_export, "rdata 3-4-5-6/", additional_folder), tbl_name, ".RData"))
  
  huxtable::quick_html(tbl_obj, file = paste0(paste0(dir_export, "html 3-4-5-6/", additional_folder), tbl_name, ".html"),
                       open = FALSE)
  huxtable::quick_rtf(tbl_obj, file = paste0(paste0(dir_export, "rtf 3-4-5-6/", additional_folder), tbl_name, ".rtf"),
                      open = FALSE)
}

smart_save <- function(df, folder, subpop = "") {
  qsave(df, paste0(folder, deparse(substitute(df)), suffix[[subpop]], ".qs"), nthreads = parallel::detectCores())
}

smart_load <- function(df, folder, subpop = "") {
  qread(paste0(folder, deparse(substitute(df)), suffix[[subpop]], ".qs"), nthreads = parallel::detectCores())
}

check_columns_exist <- function(start_df, columns) {
  colnames(start_df)[grepl(paste(columns, collapse = "|"), colnames(start_df))]
}

split_and_save <- function(.data, col_to_split) {
  col_to_split_levels <- unique(.data[, get(col_to_split)])
  tmp <- .data[get(col_to_split) == lv, ]
  for (lv in col_to_split_levels) {
    save(tmp, file = paste0(dirtemp, "TEMP_study_population_", col_to_split, "_", lv, ".RData"),
         list = "study_population")
  }
  rm(tmp)
  return(list(cols = col_to_split, levels = col_to_split_levels))
}

load_and_combine <- function(cols_splitted, levels_splitted) {
  
  persontime_list <- lapply(levels_splitted, function(x) {
    df <- get(load(paste0(dirtemp, "TEMP_persontime_", cols_splitted, "_", x, ".RData"))[[1]])
    return(df)
  })
  
  return(rbindlist(persontime_list))
}
