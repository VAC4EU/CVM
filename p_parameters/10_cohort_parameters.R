# Program Information  ----------------------------------------------------

# Program:      cohort_inputs.R 
# Description:  sets a series of "parameters" (i.e, input values) 
# Requirements: none
#
# INPUT LISTS  -----------------------------------------------------------

# functions list:
func_names <- c("cohort_tools.R")
func_dir   <- paste0(thisdir,"/p_steps/")

# Load functions
source(paste0(dirmacro, "cohort_tools.R"))

# names of input data 
#raw_data <- "D3_study_population_SCRI.RData"
#raw_data_name <- "D3_study_population_SCRI"

# SCRI variables
SCRI_variables_vocabulary <- data.table(vac4eu = c("E_GOUT_AESI", "C_MYOCARD_AESI", "C_PERICARD_AESI",
                                                   "SO_OTITISEXT_AESI", "C_VALVULAR_AESI"),
                                        scri = c("gout", "myocarditis", "pericarditis",
                                                 "otitis_externa", "valvular_heart_disease"))

#############################
#
# events and outcome variables:
#
ae_events <-  c("myocarditis") #, "pericarditis", "myopericarditis", "otitis_externa", "valvular_heart_disease")


#############################
#
#   dataset variable
#
id_original <- "person_id"
id          <- "pat_n"


###########################
#
# define whether to sample and/or match again:
#
lnew_sampling <- F
lnew_matching <- F





