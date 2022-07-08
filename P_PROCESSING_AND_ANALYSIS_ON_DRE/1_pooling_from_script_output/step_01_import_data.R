# this step retrieves data submitted by DAPs, as described
# in the parameters.R script
# append the data for further statistical analysis

# input: inbox folders

# output: 2_analysis\pooled_input_from_DAPs

# authors: Rosa Gini and Davide Messina

# version 1.0

# 17 Jan 2022

rm(list=ls(all.names=TRUE))

#set the directory where the file is saved as the working directory
if (!require("rstudioapi")) install.packages("rstudioapi")
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))
thisdir<-setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source(paste0(thisdir,"/parameters.R"))

## List files to be appended

namefile <- "D4_persontime_monthly_poisson_RF"

files_list <- unique(unlist(lapply(list_of_DAP, function(x) {list.files(list_of_submitted_folders[[x]], recursive = T, include.dirs = T)})))
files_list <- files_list[stringr::str_detect(files_list,namefile)]
files_list <- files_list[stringr::str_detect(files_list,"\\.RData")]

data <- data.table()
for (dap in list_of_DAP){
  for (f in files_list){
    if( file.exists( paste0(list_of_submitted_folders[[dap]],f)) ){
      print(paste(dap,f))
      input <- get(load(paste0(list_of_submitted_folders[[dap]],f)))
      # exclude columns referring to AESI that the DAP cannot analyse
      columnstobeincluded <- colnames(input)
      for (AESI in excluded_AESI[[dap]]){
        columnstobeincluded <- columnstobeincluded[!str_detect(columnstobeincluded,AESI)]
      }
      input <- input[,..columnstobeincluded]
      input <- input[,datasource := dap]
      data <- rbind(data, input, fill=T)
    }
  }
}

assign(namefile,data)
fwrite(data, file= paste0(dirinput_pp, namefile, ".csv"))
save(list = namefile,file= paste0(dirinput_pp, namefile,".RData"))

# save(timestamp_df, file= paste0(dirinput_pp,"timestamp_df.RData"))

# rm(data, input)


