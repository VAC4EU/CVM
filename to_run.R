#-------------------------------
# CVM script - READINESS + SCCS + SCRI + COHORT

# authors Readiness: Rosa Gini, Davide Messina, Anna Schultze

# WP4 SCCS sensitivity code for publication - 14 March 2024
# SCCS additional analyses for publication
# Based on Readiness v3.0.1

# v 3.0.1 - 26 April 2023
# Fixed parameter importation
# Recreated datasources_SCRI_SCCS_COHORT

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

# v 2.0.2 - 23 January 2022
# Bugfixes
# New covid itemset specification for UOSL
# Based on Readiness v 2.0.3

# v 2.0.1 - 19 December 2022
# Bugfix: some SCRI parameters where not calculated
# Based on Readiness v 2.0.3.Beta

# v 2.0 - 15 December 2022
# Myocarditis
# Based on Readiness v 2.0.3.Beta

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

launch_step("p_steps/05_T3_10_count_events_windows.R")
launch_step("p_steps/05_T3_11_aggregate_events_windows.R")
launch_step("p_steps/05_T3_20_create_person_time_monthly.R")
launch_step("p_steps/05_T3_21_aggregate_person_time_monthly.R")
launch_step("p_steps/05_T3_30_create_person_time_background.R")
launch_step("p_steps/05_T3_31_aggregate_person_time_background.R")

launch_step("p_steps/06_T4_10_create_D5_IR_background.R")
launch_step("p_steps/06_T4_20_create_D5_IR_background_std.R")

launch_step("p_steps/07_T5_10_final_tables.R")

### Calculation of Time Dependent variables
launch_step("p_steps/04_T2_30_create_study_population_cohort.R")
launch_step("p_steps/04_T2_40_create_TD_datasets.R")
launch_step("p_steps/04_T2_41_create_TD_NUMBER_CONDITIONS.R")

if (thisdatasource %in% datasources_SCRI_SCCS_COHORT) {
  
  # SCCS
  launch_step("p_steps/08_T3_1_clean_data.R", print.eval = TRUE)
  launch_step("p_steps/08_T3_2_select_population.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_3_describe_event_dependency.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_4a_scri_pre.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_4b_scri_post.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_4c_standard_sccs.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_4d_extended_sccs.R", print.eval = TRUE)
  
  # sensitivity analyses for peer review 
  launch_step("p_steps/08_T3_5_clean_data_SENS.R", print.eval = TRUE)
  launch_step("p_steps/08_T3_6_select_population_SENS.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_71_standard_sccs_SENS.R", print.eval = TRUE)
  launch_step("p_steps/08_T4_7b_extended_sccs_SENS.R", print.eval = TRUE)

}