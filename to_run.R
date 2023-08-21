#-------------------------------
# CVM script - READINESS + SCCS + SCRI + COHORT

# authors Readiness: Rosa Gini, Davide Messina, Anna Schultze

# v 3.3.1 - 19 July 2023
# Patch for V_MICROANGIO_AESI

# v 3.3.0 - 18 July 2023
# Updated codelist
# Fixed calculation of variables which are concepts too (V_MICROANGIO_AESI)
# UOSL optimization regarding persontime calculation
# Concepts in TD steps use only events with defined date before study_end

# v 3.2.1 - 03 July 2023
# Correct subpopulations for BIFAP

# v 3.2.1 - 03 July 2023
# Correct subpopulations for BIFAP

# v 3.2.0 - 27 June 2023
# E_DM1_AESI change in the parameter

# v 3.1.1 - 26 June 2023
# Fixed covariates
# Exclude codes
# New and updated codelists and variable metadata
# New algo E_DM1_AESI
# ICD10DA as different vocabulary than ICD10/ICD10CM

# v 3.1.0 - 29 May 2023
# sampling and matching testing

# v 3.0.9 - 26 May 2023
# Added V_THROMBOSISARTERIAL_AESI concept
# Correction for TTS algorithm

# v 3.0.8 - 12 May 2023
# Fixed component B of B_TTS_AESI

# v 3.0.7 - 05 May 2023
# Fixed DEATH missing IRs (this time for real)

# v 3.0.6 - 02 May 2023
# SCRI bugfixes
# Added vocabulary "free_text" as copy of "Free_text"
# PERSONS will be imported with day, month and year columns forced to integer
# Fixed DEATH missing IRs
# Fix for UOSL itemset
# Components for NCO

# v 3.0.5 - 28 April 2023
# Fix for COVID and pregnacy covariate

# v 3.0.4 - 27 April 2023
# Fixed additional error in SCRI script

# v 3.0.2 - 26 April 2023
# Added in to_run_short creation of event DEATH

# v 3.0.1 - 26 April 2023
# Fixed parameter importation
# Recreated datasources_SCRI_SCCS_COHORT
# NOT recommended_end_date but we should use ONLY the study_end in SCRI

# v 3.0.0 - 25 April 2023
# Updated SCRI
# Updated SCCS
# Added COHORT

# v 2.2.0 - 21 April 2023
# OUTCOMES o_deathsudden_aesi and death have been added
# Calculating and reporting IRs before in 2019 and 2020 separately
# Fixed for covid in CPD and PHARMO

# v 2.2.alpha - 19 April 2023
# Fixed VACCINES covariate
# PEDIANET with prescription
# Bugfixes

# v 2.1.1 - 07 April 2023
# Fixed TD computation in case all conceptset of a covariate are empty

# v 2.1.0 - 06 April 2023
# New codelist
# Exact matching of codes
# Itemset for UOSL
# Improved memory usage in CountPersonTime
# Creation of TD datasets
# Included levels of COVID
# Included components
# Bugfixes and general improvement

# v 2.0.3 - 23 January 2022
# New parameter specification for UOSL
# Minor bugfixes

# v 2.0.2 - 21 November 2022
# Fixed selection criteria in case there is a vax before the start of the spell
# Fixed criteria no_spells
# Added COVID to table 3-4
# Added possibility to divide too large conceptsets

# v 2.0.2 - 04 November 2022
# Fixed selection criteria in case there is a vax before the start of the spell
# Fixed criteria no_spells

# v 2.0.1 - 28 October 2022
# Fixed stability issue regarding the importation of D3_pregnancy_final
# Countpersontime for background IR from month to year

# v 2.0 - 27 October 2022
# Readiness
# updated codelist and variable names to adapt to the VAC4EU standards
# Major changes in most of the steps

# v 1.4 - 09 June 2022
# fixed bug about covid severity
# additional covid severity levels

# v 1.3 - 06 June 2022
# fixed drug proxies (except PEDIANET)
# fixed covid itemset for PEDIANET
# bugfix for end of cohort d
# fixed folder of final table in case of subpopulations

# v 1.2 - 01 June 2022
# mapped codelists of diagnosis and drug proxies to VAC4EU codelists
# added time dependent age in IR
# add MEDICINES if does not exist
# Itemset for PEDIANET

# v 1.1 - 27 May 2022
# Completed covid severity and relative IR

# v1.0.0 - 9 March 2022
# First release of SCCS


rm(list=ls(all.names=TRUE))

#set the directory where the file is saved as the working directory
if (!require("rstudioapi")) install.packages("rstudioapi")
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##%######################################################%##
#                                                          #
####                     PARAMETERS                     ####
#                                                          #
##%######################################################%##

source(paste0(thisdir,"/p_parameters/01_parameters_program.R"))
source(paste0(thisdir,"/p_parameters/02_parameters_CDM.R"))
source(paste0(thisdir,"/p_parameters/03_concept_sets.R"))
source(paste0(thisdir,"/p_parameters/04_itemsets.R"))
source(paste0(thisdir,"/p_parameters/05_subpopulations_restricting_meanings.R"))
source(paste0(thisdir,"/p_parameters/06_variable_lists.R"))
source(paste0(thisdir,"/p_parameters/07_algorithms.R"))
source(paste0(thisdir,"/p_parameters/08_SCRI_parameters.R"))
source(paste0(thisdir,"/p_parameters/09_SCCS_parameters.R"))
source(paste0(thisdir,"/p_parameters/10_cohort_parameters.R"))
source(paste0(thisdir,"/p_parameters/11_design_parameters.R"))
source(paste0(thisdir,"/p_parameters/99_saving_all_parameters.R"))


##%######################################################%##
#                                                          #
####                    MAIN SCRIPT                     ####
#                                                          #
##%######################################################%##

launch_step("p_steps/01_T2_10_create_persons.R")
launch_step("p_steps/01_T2_20_apply_CreateSpells.R")
launch_step("p_steps/01_T2_31_CreateConceptSetDatasets.R")
launch_step("p_steps/01_T2_32_CreateItemSetDatasets.R")
launch_step("p_steps/01_T2_33_CreatePromptSetDatasets.R")
launch_step("p_steps/01_T2_40_clean_vaccines.R")
launch_step("p_steps/01_T2_41_apply_criteria_for_doses.R")
launch_step("p_steps/01_T2_50_clean_spells.R")
launch_step("p_steps/01_T2_60_selection_criteria_from_PERSON_to_study_population.R")

launch_step("p_steps/02_T3_10_create_study_population.R")

launch_step("p_steps/03_T2_10_create_D3_outcomes_simple_algorithm.R")
launch_step("p_steps/03_T2_11_create_D3_outcomes_complex_algorithm.R")
launch_step("p_steps/03_T2_12_create_D3_event_outcomes_ALL.R")
launch_step("p_steps/03_T2_20_create_D3_covid_episodes.R")
launch_step("p_steps/03_T2_21_COVID_severity_hospitalised.R")
launch_step("p_steps/03_T2_22_COVID_severity_ICU.R")
launch_step("p_steps/03_T2_23_COVID_severity_DEATH.R")
launch_step("p_steps/03_T2_24_TD_COVID_severity_levels.R")

launch_step("p_steps/03_T2_30_create_covariates.R")
launch_step("p_steps/03_T2_40_create_components.R")

launch_step("p_steps/04_T2_10_create_total_study_population.R")
launch_step("p_steps/04_T2_20_SCRI.R")

### Calculation of Time Dependent variables
launch_step("p_steps/04_T2_30_create_study_population_cohort.R")
launch_step("p_steps/04_T2_40_create_TD_datasets.R")
launch_step("p_steps/04_T2_41_create_TD_NUMBER_CONDITIONS.R")

launch_step("p_steps/05_T3_10_count_events_windows.R")
launch_step("p_steps/05_T3_11_aggregate_events_windows.R")
launch_step("p_steps/05_T3_20_create_person_time_monthly.R")
launch_step("p_steps/05_T3_21_aggregate_person_time_monthly.R")
launch_step("p_steps/05_T3_30_create_person_time_background.R")
launch_step("p_steps/05_T3_31_aggregate_person_time_background.R")

launch_step("p_steps/06_T4_10_create_D5_IR_background.R")
launch_step("p_steps/06_T4_20_create_D5_IR_background_std.R")

launch_step("p_steps/07_T5_10_final_tables.R")
