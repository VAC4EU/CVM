dirinput_pp <- paste0(thisdir, paste("/../2_analysis/pooled_input_from_DAPs/"))


base_dir <- "Z:/inbox"

list_of_submitted_folders <- vector(mode="list")

# step 0: included data from submissions of DAPs (comment line if DAP is not present)
list_of_submitted_folders[["BIFAP_PC"]]<-c("Z:/inbox/transfer-2021-11-21-02-45-pm/g_export_PC/")
list_of_submitted_folders[["BIFAP_PC_HOSP"]]<-c("Z:/inbox/transfer-2021-11-21-02-47-pm/g_export_PC_HOSP/")
list_of_submitted_folders[["ARS"]]<-c("Z:/inbox/transfer-2021-11-25-11-34-pm/g_export/")
list_of_submitted_folders[["PHARMO"]]<-c("Z:/inbox/transfer-2021-11-24-07-40-pm/Results of 6.5_24112021/g_export/")
list_of_submitted_folders[["CPRD"]]<-c("Z:/inbox/transfer-2021-11-22-09-13-pm/")



list_of_DAP <- names(list_of_submitted_folders)


# AESI to be excluded per each DAP

excluded_AESI <- vector(mode="list")


excluded_AESI[['BIFAP_PC']] <- c("DM1_narrow","DM1_broad","TRANSMYELITIS_narrow", "TRANSMYELITIS_broad", "TRANSMYELITIS","ANOSMIA_narrow", "ANOSMIA_broad", "ANOSMIA")
excluded_AESI[['BIFAP_PC_HOSP']] <- c("DM1_narrow","DM1_broad","TRANSMYELITIS_narrow", "TRANSMYELITIS_broad", "TRANSMYELITIS","ANOSMIA_narrow", "ANOSMIA_broad", "ANOSMIA")
excluded_AESI[['PHARMO']] <- c("DM1_narrow","DM1_broad","DEATH","ACUASEARTHRITIS", "ADEM", "AKI", "ALI",  "ARD",  "DIC" , "ERYTH", "MICROANGIO", "MISCC", "NARCOLEPSY", "SOCV", "STRCARD",  "Sinusthrom", "TMA", "TRANSMYELITIS", "CVSTNoTP", "CVSTTP", "COVID_L4plus", "COVID_L5plus")
excluded_AESI[['CPRD']] <- c("DM1_narrow","DM1_broad")
excluded_AESI[['ARS']] <- c("DM1_narrow","DM1_broad")


# dirinput_pp_subfolder1<- paste0(dirinput_pp, "D4 tables")
# dirinput_pp_subfolder2<- paste0(dirinput_pp, "dashboard tables")
# dirinput_pp_subfolder3<- paste0(dirinput_pp, "Dummy tables for report")
# dirinput_pp_subfolder4<- paste0(dirinput_pp, "Dummy tables for report MIS-KD")
# dirinput_pp_subfolder5<- paste0(dirinput_pp, "Dummy tables October")

suppressWarnings(if (!file.exists(dirinput_pp)) dir.create(file.path( dirinput_pp)))
# suppressWarnings(if (!file.exists(dirinput_pp_subfolder1)) dir.create(file.path( dirinput_pp_subfolder1)))
# suppressWarnings(if (!file.exists(dirinput_pp_subfolder2)) dir.create(file.path( dirinput_pp_subfolder2)))
# suppressWarnings(if (!file.exists(dirinput_pp_subfolder3)) dir.create(file.path( dirinput_pp_subfolder3)))
# suppressWarnings(if (!file.exists(dirinput_pp_subfolder4)) dir.create(file.path( dirinput_pp_subfolder4)))
# suppressWarnings(if (!file.exists(dirinput_pp_subfolder5)) dir.create(file.path( dirinput_pp_subfolder5)))

# files_to_join <- c("Dummy tables for report/Attrition diagram 1.csv",
#                    "Dummy tables for report/Attrition diagram 2.csv",
#                    "Dummy tables for report/Cohort characteristics at start of study (1-1-2020).csv",
#                    "Dummy tables for report/Doses of COVID-19 vaccines and distance between first and second dose.csv",
#                    "Dummy tables for report/Number of incident cases entire study period.csv",
#                    "Dummy tables October/Table 5, Attrition diagram 1.csv",
#                    "Dummy tables October/Table 6, Cohort characteristics at start of study (1-1-2020).csv",
#                    "Dummy tables October/Table 11, COVID-19 vaccination by dose and time period between first and second dose (days).csv",
#                    "Dummy tables October/Table 12, Number of incident cases entire study period.csv",
#                    "Dummy tables October/Table XX, Incidence rates of AESI by vaccine and datasource.csv")
# files_slated <- list(c("Dummy tables for report/table 3 Cohort characteristics at first COVID-19 vaccination Italy_ARS.csv",
#                        "Dummy tables for report/table 4 Cohort characteristics at first COVID-19 vaccination Netherlands-PHARMO.csv",
#                        "Dummy tables for report/table 5 Cohort characteristics at first COVID-19 vaccination UK_CPRD.csv",
#                        "Dummy tables for report/table 6 Cohort characteristics at first COVID-19 vaccination ES_BIFAP.csv"),
#                      c("Dummy tables October/Table 7, Cohort characteristics at first dose of COVID-19 vaccine Italy_ARS.csv",
#                        "Dummy tables October/Table 8, Cohort characteristics at first dose of COVID-19 vaccine Netherlands-PHARMO.csv",
#                        "Dummy tables October/Table 9, Cohort characteristics at first dose of COVID-19 vaccine UK_CPRD.csv",
#                        "Dummy tables October/Table 10, Cohort characteristics at first dose of COVID-19 vaccine ES_BIFAP.csv"),
#                      c("Dummy tables October/table 7, Cohort characteristics at second dose of COVID-19 vaccine Italy_ARS.csv",
#                        "Dummy tables October/table 8, Cohort characteristics at second dose of COVID-19 vaccine Netherlands-PHARMO.csv",
#                        "Dummy tables October/table 9, Cohort characteristics at second dose of COVID-19 vaccine UK_CPRD.csv",
#                        "Dummy tables October/table 10, Cohort characteristics at second dose of COVID-19 vaccine ES_BIFAP.csv"))
# names_files_slated <- c("/Table 3_4_5_6 Cohort characteristics at first COVID-19 vaccination (dashboard)",
#                         "/Table 7_8_9_10 Cohort characteristics at first dose of COVID-19 vaccine (dashboard)",
#                         "/Table 7_8_9_10 Cohort characteristics at second dose of COVID-19 vaccine (dashboard)")

# load packages
if (!require("haven")) install.packages("haven")
library(haven)
if (!require("tidyverse")) install.packages("tidyverse")
library(dplyr)
if (!require("lubridate")) install.packages("lubridate")
library(lubridate)
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("AdhereR")) install.packages("AdhereR")
library(AdhereR)
if (!require("stringr")) install.packages("stringr")
library(stringr)
if (!require("purrr")) install.packages("purrr")
library(purrr)
if (!require("readr")) install.packages("readr")
library(readr)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("survival")) install.packages("survival")
library(survival)

#############
`%not in%` <- negate(`%in%`)
