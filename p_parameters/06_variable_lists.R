#---------------------------------------------------------------------------
# create few lists, each containing names of variables of interest for different D3s

# input: the list of variable names associated to the algorithms, created by the PI: the name of the file was stored in File_variables_ALG_DP_ROC20 in step 03_concept_sets

#output:
# VAR_codelist - dataset with instructions (also used in step 07)
# OUTCOME_variables - list of names of variables that are outcomes
# recurrent_OUTCOME_variables - subset of OUTCOME_variables inclusing only outcomes that are recurrent
# CONTROL_variables - list of names of variables that are negative controls
# COV_variables_raw - list of names of variables that are covariates
# DP_variables - subset of COV_variables_raw including only drug proxies
# COV_variables - subset of COV_variables_raw including only diagnosis
# VACCINES_variable - name of the variable containing the vaccines

VAR_codelist <- readxl::read_excel(File_variables_ALG_DP_ROC20, sheet = "Variables")
VAR_codelist <- unique(as.data.table(VAR_codelist))

# TODO remove before release?
VAR_codelist <- VAR_codelist[Varname == "D_Diverticulitis_AESI", Varname := "D_DIVERTICULITIS_AESI"]

# Divide drug and diagnosis
DRUG_codelist <- VAR_codelist[(DP), ]
VAR_codelist <- VAR_codelist[!(DP), ]

# list of variable names of OUTCOMES
OUTCOME_variables <- c(VAR_codelist[(AESI), Varname], DRUG_codelist[(AESI), Varname])
# OUTCOME_variables <- OUTCOME_variables[OUTCOME_variables %in% SCRI_list_variables]

# list of variable names of CONTROLS
CONTROL_variables <- VAR_codelist[(NEG), Varname]
# CONTROL_variables <- CONTROL_variables[CONTROL_variables %in% SCRI_list_variables]
# list of variable names of COVARIATES
COV_variables_raw <- c(VAR_codelist[(COV), Varname], DRUG_codelist[(COV), Varname])

# Variable for COVID_VACCINES
VACCINES_variable <- "COVID_VACCINES"

# TODO check recurrent events
recurrent_OUTCOME_variables <- c("Im_ANAPHYLAXIS_AESI")

# Creating DP_variables from COV manually
# list of variable names of DRUG_PROXIES
DP_variables <- COV_variables_raw[grepl("^DP_", COV_variables_raw)]

# list of variable names of COVARIATES without DRUG_PROXIES
COV_variables <- setdiff(COV_variables_raw, DP_variables) 

# TODO test if needed
# variables_of_our_study <- c(VAR_codelist[, Varname], DRUG_codelist[, Drug_proxie])
# auxiliary_variables <- c(VAR_codelist[(Algorithm_input) & !((COV) | (NEG) | (AESI)), Varname],
#                          DRUG_codelist[(Algorithm_input) & !((COV) | (AESI)), Drug_proxie])
