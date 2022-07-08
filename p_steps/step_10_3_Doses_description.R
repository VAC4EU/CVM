##############################################################################
##########################    Doses description     ##########################
##############################################################################

# input: DOSES_BIRTHCOHORTS
# output: HTML_DOSES_BIRTHCOHORTS_description

# library
if (!require("rmarkdown")) install.packages("rmarkdown")
library(rmarkdown )

# output file path
PathOutputFolder=paste0(thisdir,"/g_describeHTML")

# creating output file folder 
suppressWarnings(if (!file.exists(PathOutputFolder)) dir.create(file.path( PathOutputFolder)))

#loading the flowchart


for (subpop in subpopulations_non_empty ) {
  
  thisdirexp <- ifelse(this_datasource_has_subpopulations == FALSE,direxp,direxpsubpop[[subpop]])
  
  if(this_datasource_has_subpopulations == T) dirdashboard <-paste0(thisdirexp,"dashboard tables/")
  
  DOSES_BIRTHCOHORTS<- fread(paste0(dirdashboard,"DOSES_BIRTHCOHORTS.csv"))
  
  #rendering the file
  render(paste0(dirmacro,"DOSES_BIRTHCOHORTS_description.Rmd"),           
         output_dir=PathOutputFolder,
         output_file=paste0("HTML_DOSES_BIRTHCOHORTS_description",suffix[[subpop]]), 
         params=list(Dataset=DOSES_BIRTHCOHORTS))
}