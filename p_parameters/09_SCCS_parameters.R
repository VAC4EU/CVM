# Program Information  ----------------------------------------------------
#
# Program:      00_sensitivity_functions
# Description:  load packages required to support pipeline 
#               run custom functions required 
#               *utility functions (to plot basic summaries)
#               *event_dependency_table
#               *event_dependency_plot
#               *event_information_table 
#               *sccs_analysis 
#
# Dependencies: Function 5 sccs_analysis requires and loads packages tidyverse, SCCS
#
# 0. HOUSEKEEPING ------------------------------------------------------------

# specify the input data name 
raw_data_name <- "D3_study_population_SCRI"

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp, "sccs_sensitivity")),           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir, "/log_files/sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)
for (subpop in subpopulations_non_empty) {
  dir.create(file.path(paste0(direxpsubpop[[subpop]], "sccs_sensitivity")),         showWarnings = FALSE, recursive = TRUE)
}

# study_end <- ymd(20261231)