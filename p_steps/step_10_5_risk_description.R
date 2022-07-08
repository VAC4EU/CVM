##############################################################################
##########################     risk description     ##########################
##############################################################################

# input: RISK_BIRTHCOHORTS_CALENDARTIME, RISK_BIRTHCOHORTS_TIMESINCEVACCINATION
# output: HTML_risk_description

# library
if (!require("rmarkdown")) install.packages("rmarkdown")
library(rmarkdown )

# output file path
PathOutputFolder=paste0(thisdir,"/g_describeHTML")

# creating output file folder 
suppressWarnings(if (!file.exists(PathOutputFolder)) dir.create(file.path( PathOutputFolder)))

for (subpop in subpopulations_non_empty ) {
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  if(this_datasource_has_subpopulations == T) dirdashboard <-paste0(thisdirexp,"dashboard tables/")
  
#loading the datasets
RISK_BIRTHCOHORTS_CALENDARTIME <- fread(paste0(dirdashboard, "RISK_BIRTHCOHORTS_CALENDARTIME.csv"))
RISK_BIRTHCOHORTS_TIMESINCEVACCINATION <- fread( paste0(dirdashboard, "RISK_BIRTHCOHORTS_TIMESINCEVACCINATION.csv"))

#rendering the file
render(paste0(dirmacro,"risk_Description.Rmd"),           
       output_dir=PathOutputFolder,
       output_file=paste0("HTML_risk_description",suffix[[subpop]]), 
       params=list(RISK_BIRTHCOHORTS_CALENDARTIME = RISK_BIRTHCOHORTS_CALENDARTIME,
                   RISK_BIRTHCOHORTS_TIMESINCEVACCINATION = RISK_BIRTHCOHORTS_TIMESINCEVACCINATION)) 
}
