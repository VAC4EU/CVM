##############################################################################
##########################    benefit description   ##########################
##############################################################################

# input: BENEFIT_BIRTHCOHORTS_CALENDARTIME, BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION
# output: HTML_benefit_description

# library
if (!require("rmarkdown")) install.packages("rmarkdown")
library(rmarkdown )

# output file path
PathOutputFolder=paste0(thisdir,"/g_describeHTML")

# creating output file folder 
suppressWarnings(if (!file.exists(PathOutputFolder)) dir.create(file.path( PathOutputFolder)))

#loading the datasets

for (subpop in subpopulations_non_empty ) {
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  if(this_datasource_has_subpopulations == T) dirdashboard <-paste0(thisdirexp,"dashboard tables/")
  
BENEFIT_BIRTHCOHORTS_CALENDARTIME <- fread(paste0(dirdashboard, "BENEFIT_BIRTHCOHORTS_CALENDARTIME.csv"))
BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION <- fread(paste0(dirdashboard, "BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION.csv"))
BENEFIT_BIRTHCOHORTS_CALENDARTIME[birth_cohort == "all_birth_cohorts" & vx_manufacturer == "Pfizer"  & COVID == "COVID"]

#rendering the file
render(paste0(dirmacro,"benefit_Description.Rmd"),           
       output_dir=PathOutputFolder,
       output_file=paste0("HTML_benefit_description",suffix[[subpop]]), 
       params=list(BENEFIT_BIRTHCOHORTS_CALENDARTIME = BENEFIT_BIRTHCOHORTS_CALENDARTIME,
                   BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION = BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION)) 
}


