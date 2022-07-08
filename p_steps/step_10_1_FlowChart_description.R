##############################################################################
##########################  FlowChart description   ##########################
##############################################################################

# input: Flowchart_doses
# output: HTML_Flowchart_doses_description

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
  
  Flowchart_doses<- fread(paste0(thisdirexp,"Flowchart_doses.csv"))
  
  #rendering the file
  render(paste0(dirmacro,"FlowChart_Description.Rmd"),           
         output_dir=PathOutputFolder,
         output_file="HTML_Flowchart_doses_description", 
         params=list(FlowChart = Flowchart_doses))  
}

