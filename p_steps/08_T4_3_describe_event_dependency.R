# Program Information  ----------------------------------------------------
#
# Program:      03_describe_event_dependency.R 
# Description:  create descriptions of event dependency of events. Run once irrespective of design. 
#               the plot is made for all exposed individuals with the outcome, so essentially in the SCCS population
# Dependencies: 01_clean_data.R  
#               02_select_population.R 
#               00_sensitivity_functions 
#               tidyverse, ggplot2
# 
# 0. HOUSEKEEPING ---------------------------------------------------------

# Start the subpopulation
for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  # 1. READ IN DATA ------------------------------------------------------------
  load(paste0(dirtemp, "sccs_sensitivity/", "sccs_population", suffix[[subpop]], ".RData"))
  
  # TABULATE EXPOSURES AND EVENTS -------------------------------------------
  # use the function event_dependency_table to plot basic information on the relationship between dose and event 
  
  # for the overall population 
  table2b.interim <- lapply(split(sccs_population, sccs_population$type_vax1), event_dependency_table)
  table2b <- do.call(rbind, table2b.interim )
  
  table2b <- cbind(Names = rownames(table2b), table2b)
  table2b$Names <- substring(table2b$Names,1,nchar(table2b$Names)-2)
  
  write.csv(table2b, paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", "table_2b_event_dependency_by_type.csv"),
            row.names = FALSE)
  
  # by vaccine type 
  table2b.interim <- lapply(split(sccs_population, sccs_population$type_vax1), event_dependency_table)
  table2b <- do.call(rbind, table2b.interim )
  
  table2b <- cbind(Names = rownames(table2b), table2b)
  table2b$Names <- substring(table2b$Names,1,nchar(table2b$Names)-2)
  
  write.csv(table2b, paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", "table_2b_event_dependency_by_type.csv"),
            row.names = FALSE)
  
  # EXPOSURE-CENTRED INTERVAL PLOTS -----------------------------------------
  # use the function event_dependency_plot to make exposure-centered interval plots by vaccine type 
  
  dose1_plots <- lapply(split(sccs_population, sccs_population$type_vax1), event_dependency_plot, 1)
  dose2_plots <- lapply(split(sccs_population, sccs_population$type_vax1), event_dependency_plot, 2)
  dose3_plots <- lapply(split(sccs_population, sccs_population$type_vax1), event_dependency_plot, 3)
  
  # save dose 1 plots 
  dose1_names <- names(dose1_plots)
  for (i in dose1_names) {
    
    plot <- dose1_plots[[i]]
    file_name = paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", "plot_d1_event_dependency_plot_", i, ".pdf")
    ggsave(file_name, plot = plot)
    while (!is.null(dev.list()))  dev.off()
    
  }
  
  # save dose 2 plots 
  dose2_names <- names(dose2_plots)
  for (i in dose2_names) {
    
    plot <- dose2_plots[[i]]
    file_name = paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", "plot_d2_event_dependency_plot_", i, ".pdf")
    ggsave(file_name, plot = plot)
    while (!is.null(dev.list()))  dev.off()
    
  }
  
  # save dose 3 plots 
  dose3_names <- names(dose3_plots)
  for (i in dose3_names) {
    
    plot <- dose3_plots[[i]]
    file_name = paste0(direxpsubpop[[subpop]], "sccs_sensitivity/", "plot_d3_event_dependency_plot_", i, ".pdf")
    ggsave(file_name, plot = plot)
    while (!is.null(dev.list()))  dev.off()
    
  }
  
}

