#-------------------------------
# CVM script

# authors: Rosa Gini, Olga Paoletti, Davide Messina, Giorgio Limoncella
# authors: Anna Schultze, Svetlana Belitser; Ema Alsina, Sophie Bots, Ivonne Martens 

# v 1.0.6 - 25 March 2022
# Fixed a bug in step 12_2 and scri tools

# v 1.0.5 - 21 March 2022
# Fixed a bug in step 07_5
# Implemented myocarditis cohort calculation
# Fixed a bug in step 12_1

# v 1.0.4 - 15 March 2022
# Fixed a bug in Poisson dataset creation
# Fixed a bug in step 12_2.
# Parallel+computing for step 12_2
# SCRI models will now be saved in g_intermediate

# v 1.0.3 - 08 March 2022
# Fixed scri_tools
# Added region 13 for BIFAP
# Memory optimization for step 06_6 and 07_9

# v 1.0.2 - 05 March 2022
# Added CI's in plots and fixed some issues
# Bugfixes in scri_tools, 12_1 and 12_2

# v 1.0.1 - 19 February 2022
# Added free text search option for concepts
# Fixed a bug in step 12_1 and 12_2
# Version 1.0 is comparable as version 6.9 in ECVM

# v 1.0 - 14 February 2022
# Modified step 12_2, created macro scri_tools
# This is a fork of the ECVM project (https://github.com/ARS-toscana/ECVM)
# Version 1.0 is comparable as version 6.9 in ECVM

rm(list=ls(all.names=TRUE))

#set the directory where the file is saved as the working directory
if (!require("rstudioapi")) install.packages("rstudioapi")
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#load parameters
source(paste0(thisdir,"/p_parameters/01_parameters_program.R"))
source(paste0(thisdir,"/p_parameters/02_parameters_CDM.R"))
source(paste0(thisdir,"/p_parameters/03_concept_sets.R"))
source(paste0(thisdir,"/p_parameters/04_itemsets.R"))
source(paste0(thisdir,"/p_parameters/05_subpopulations_restricting_meanings.R"))
source(paste0(thisdir,"/p_parameters/06_algorithms.R"))
source(paste0(thisdir,"/p_parameters/07_scri_inputs.R"))


#run scripts

# 01 RETRIEVE RECORDS FRM CDM

system.time(source(paste0(thisdir,"/p_steps/step_01_1_T2.1_create_conceptset_datasets.R")))
system.time(source(paste0(thisdir,"/p_steps/step_01_2_T2.1_create_spells.R")))

system.time(source(paste0(thisdir,"/p_steps/step_01_3_T2.1_create_dates_in_PERSONS.R")))

system.time(source(paste0(thisdir,"/p_steps/step_01_4_T2.1_create_prompt_and_itemset_datasets.R")))

#02 quality checks

system.time(source(paste0(thisdir,"/p_steps/step_02_1_T2_create_QC_criteria.R")))
system.time(source(paste0(thisdir,"/p_steps/step_02_2_T3_apply_QC_exclusion_criteria.R")))


#03 create exclusion criteria
system.time(source(paste0(thisdir,"/p_steps/step_03_1_T2_create_exclusion_criteria.R")))
system.time(source(paste0(thisdir,"/p_steps/step_03_2_T2_merge_persons_concept.R")))

#04 apply exclusion criteria
system.time(source(paste0(thisdir,"/p_steps/step_04_1_T3_apply_exclusion_criteria.R")))
system.time(source(paste0(thisdir,"/p_steps/step_04_2_T3_apply_quality_check_exclusion_criteria_doses.R")))
##use flowchart (apply also quality checks)

#10 describing the datasets
system.time(source(paste0(thisdir,"/p_steps/step_10_1_FlowChart_description.R")))

#05 create D3 for doses and coverage
system.time(source(paste0(thisdir,"/p_steps/step_05_1_T2.2_components.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_2_T2.2_secondary_components.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_3_T2_create_events_ALL_OUTCOMES.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_4_QC_count_codes_ALL_OUTCOMES.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_5_QC_apply_component_strategy.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_6_T2.2_covariates_at_baseline.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_7_T2.2_DP_at_baseline.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_8_T2.3_baseline_characteristics.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_9_T2.3_ALL_covariates_at_baseline_V2.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_10_T2.2_components_COVID_severity.R")))
system.time(source(paste0(thisdir,"/p_steps/step_05_11_T2.3_algorithms_COVID_severity.R")))

#06 create D3 for doses and coverage
system.time(source(paste0(thisdir,"/p_steps/step_06_1_T2_create_D3_datasets.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_2_T2.2_covariates_at_vaccination.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_3_T2.2_DP_at_vaccination.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_4_T2.3_vaccination_characteristics.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_5_T2.3_ALL_covariates_at_vaccination_V2.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_6_T2_create_D3_datasets.R")))

#MIS and Myocarditis section
# create D3 MIS population

system.time(source(paste0(thisdir,"/p_steps/step_06_7_MIS_population.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_8_T2.2_covariates_at_covid.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_9_T2.2_DP_at_covid.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_10_T2.3_covid_characteristics.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_11_T2.3_ALL_covariates_at_covid_V2.R")))
system.time(source(paste0(thisdir,"/p_steps/step_06_12_MIS_population_c.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_9_T3_create_person_time_MIS_year.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_10_T3_aggregate_monthly_MIS.R")))
system.time(source(paste0(thisdir,"/p_steps/step_08_2_T4_IR_MIS.R")))

#descriptive
system.time(source(paste0(thisdir,"/p_steps/step_09_4_T3_create_D4_descriptive_tables_MIS.R")))
system.time(source(paste0(thisdir,"/p_steps/step_11_2_T4_create_dummy_tables_MIS_KD.R")))

##end MIS and Myocarditis section------------

#07 create persontime
system.time(source(paste0(thisdir,"/p_steps/step_07_1_T3_create_person_time_risks.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_2_T3_create_person_time_risks_year.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_3_T3_create_person_time_benefits.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_4_T3_create_person_time_benefits_year.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_5_T3_aggregate_sex_birth_cohort.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_6_T3_aggregate_sex_risk_factor.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_7_T3_create_person_time_vax_cohort.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_8_T3_aggregate_monthly.R")))



#08 Calculate Incidence Rates
system.time(source(paste0(thisdir,"/p_steps/step_08_1_T4_IR.R")))


#09 create D4 for doses and coverage
system.time(source(paste0(thisdir,"/p_steps/step_09_1_T3_create_D4_doses_weeks.R")))
system.time(source(paste0(thisdir,"/p_steps/step_09_2_T3_create_D4_descriptive_tables.R")))
system.time(source(paste0(thisdir,"/p_steps/step_09_3_T3_create_dashboard_tables.R")))
system.time(source(paste0(thisdir,"/p_steps/step_09_4_T3_create_D4_descriptive_tables_MIS.R")))


# system.time(source(paste0(thisdir,"/p_steps/step_10_2_Coverage_description.R")))
# system.time(source(paste0(thisdir,"/p_steps/step_10_3_Doses_description.R")))
# system.time(source(paste0(thisdir,"/p_steps/step_10_4_benefit_description.R")))


#11 Create descriptive tables
system.time(source(paste0(thisdir,"/p_steps/step_11_1_T4_create_dummy_tables.R")))
system.time(source(paste0(thisdir,"/p_steps/step_11_2_T4_create_dummy_tables_MIS_KD.R")))
system.time(source(paste0(thisdir,"/p_steps/step_11_3_T4_create_dummy_tables_October.R")))


#POISSON section
system.time(source(paste0(thisdir,"/p_steps/step_06_13_Poisson_population.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_11_T3_create_person_time_poisson.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_12_T3_aggregate_monthly_Poisson.R")))


#SCRI section
#create D3 MIS population
system.time(source(paste0(thisdir,"/p_steps/step_06_14_SCRI_population.R")))
system.time(source(paste0(thisdir,"/p_steps/step_12_1_create_scri_dataset.R")))
system.time(source(paste0(thisdir,"/p_steps/step_12_2_run_scri.R")))
system.time(source(paste0(thisdir,"/p_steps/step_12_3_count_subgroup_numbers.R")))
# system.time(source(paste0(thisdir,"/p_steps/step_12_4_prepare_meta_dataset.R")))

#MYOCARD cohort
system.time(source(paste0(thisdir,"/p_steps/step_06_15_Cohort_population.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_13_T3_create_person_time_windows.R")))
system.time(source(paste0(thisdir,"/p_steps/step_07_14_T3_aggregate_windows.R")))
system.time(source(paste0(thisdir,"/p_steps/step_08_3_T4_IR_windows.R")))
