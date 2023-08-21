for (subpop in subpopulations_non_empty) {
  print(subpop)
  
  # Table 1
  print("Now creating: Table 1")
  flow_source <- fread(paste0(direxpsubpop[[subpop]], "Flowchart_exclusion_criteria.csv"))
  
  excl_criteria <- c("sex_or_birth_date_is_not_defined", "birth_date_absurd", "partial_date_of_death", "no_spells",
                     "all_spells_start_after_ending", "no_spell_overlapping_the_study_period",
                     "no_spell_longer_than_365_days", "all_spells_include_vax1_but_less_than_365_days_from_it",
                     "higher_doses_included_but_lower_doses_missing")
  
  sex_birth_death_entry_exit_missing <- check_columns_exist(flow_source, c("sex_or_birth_date_is_not_defined", "birth_date_absurd",
                                                                           "partial_date_of_death", "no_spells", "all_spells_start_after_ending"))
  death_exit_before_2019 <- check_columns_exist(flow_source, c("no_spell_overlapping_the_study_period"))
  less_365 <- check_columns_exist(flow_source, c("no_spell_longer_than_365_days"))
  no_lookback_vax1 <- check_columns_exist(flow_source, c("all_spells_include_vax1_but_less_than_365_days_from_it"))
  lower_doses <- check_columns_exist(flow_source, c("higher_doses_included_but_lower_doses_missing"))
  
  flow_source[, excl_criteria := names(flow_source)[apply(flow_source, 1, which.min)]]
  which_all_1 <- apply(flow_source[, .SD, .SDcols = c(sex_birth_death_entry_exit_missing, death_exit_before_2019,
                                                      less_365, lower_doses)], 1, function(x) all(as.logical(x)))
  flow_source[which_all_1, excl_criteria := "Final study population"]
  
  flow_source[excl_criteria %in% sex_birth_death_entry_exit_missing, order_col := 1]
  flow_source[excl_criteria %in% sex_birth_death_entry_exit_missing, excl_criteria := "Sex or birth date missing or absurd or no dates of entry or exit"]
  
  flow_source[excl_criteria %in% death_exit_before_2019, order_col := 2]
  flow_source[excl_criteria %in% death_exit_before_2019, excl_criteria := "Exit from the data source before 1/1/2019"]
  
  flow_source <- flow_source[, .(excl_criteria, order_col, N)]
  
  flow_source <- rbind(flow_source, data.frame(excl_criteria = "Persons in the instance of the data source",
                                               N = t(colSums(flow_source[, -c(1, 2)])), order_col = 0))
  
  
  flow_source <- rbind(flow_source,
                       data.frame(excl_criteria = "Persons in the data source at or after 1/1/2019",
                                  N = t(colSums(flow_source[excl_criteria %in% c(less_365, no_lookback_vax1, lower_doses,
                                                                                 "Final study population"),
                                                            -c(1, 2)])), order_col = 3))
  
  flow_source[excl_criteria %in% less_365, order_col := 4]
  flow_source[excl_criteria %in% less_365, excl_criteria := "Less than 365 days history at any point in time after 1.1.2019 (and not born after 1.1.2019)"]
  
  
  flow_source[excl_criteria %in% no_lookback_vax1, order_col := 5]
  
  flow_source[excl_criteria %in% lower_doses, order_col := 6]
  flow_source[excl_criteria == "Final study population", order_col := nrow(flow_source)]
  
  flow_source <- flow_source[, .(excl_criteria, order_col, N)]
  
  flow_source <- flow_source[, lapply(.SD, sum), by = c("excl_criteria", "order_col")]
  
  setorder(flow_source, order_col)
  flow_source <- flow_source[, .(excl_criteria, N)]
  
  setnames(flow_source, c("excl_criteria", "N"), c("Exclusion criteria", paste0(thisdatasource, suffix[[subpop]])))
  
  fwrite(flow_source, file = paste0(dirtablesubpop[[subpop]], "Table 1 - Attrition diagram.csv"))
  save(flow_source, file = paste0(dirtablesubpop[[subpop]], "Table 1 - Attrition diagram.RData"))
  
  
  
  ### Table 2
  print("Now creating: Table 2")
  load(paste0(dirtemp, "D3_events_ALL_OUTCOMES", suffix[[subpop]], ".RData"))
  events_ALL_OUTCOMES <- get(paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  rm(list = paste0("D3_events_ALL_OUTCOMES", suffix[[subpop]]))
  
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, .N, by = c("type_outcome", "meaning_renamed",
                                                          "codvar", "event_record_vocabulary")]
  
  events_ALL_OUTCOMES <- events_ALL_OUTCOMES[, .(DAP = paste0(thisdatasource, suffix[[subpop]]), "AESI/NCO" = type_outcome,
                                                 "Coding system" = event_record_vocabulary, "Code" = codvar,
                                                 "Meaning" = meaning_renamed, "Count" = N)]
  
  fwrite(events_ALL_OUTCOMES, file = paste0(dirtablesubpop[[subpop]],
                                            "Table 2.1-x - Code count for events during follow-up (per AESI-NCO) and by DAP (as excel).csv"))
  save(events_ALL_OUTCOMES, file = paste0(dirtablesubpop[[subpop]], "Table 2.1-x - Code count for events during follow-up (per AESI-NCO) and by DAP (as excel).RData"))
  
  
  
  
  ### Table 3
  print("Now creating: Table 3")
  
  # Load D3_Total_study_population
  load(paste0(dirtemp, "D3_Total_study_population", suffix[[subpop]], ".RData"))
  study_population <- get(paste0("D3_Total_study_population", suffix[[subpop]]))
  rm(list = paste0("D3_Total_study_population", suffix[[subpop]]))
  
  # # Load D3_Total_study_population
  # load(paste0(diroutput, "D4_persontime_risk_month_aggregated", suffix[[subpop]], ".RData"))
  # PT_aggr <- get(paste0("D4_persontime_risk_month_aggregated", suffix[[subpop]]))
  # rm(list = paste0("D4_persontime_risk_month_aggregated", suffix[[subpop]]))
  
  # Load covariates
  load(paste0(dirtemp,"D3_covariates_ALL",suffix[[subpop]],".RData"))
  covariates <- get(paste0("D3_covariates_ALL", suffix[[subpop]]))
  rm(list = paste0("D3_covariates_ALL", suffix[[subpop]]))
  
  table_3 <- create_table_characteristic_population(study_pop = study_population,
                                                    covariates = covariates, at_time = "baseline",
                                                    agebands = Agebands_countpersontime)
  
  table3_name <- "Table 3 - Baseline characteristics at start of follow-up for entire study population"
  
  save_tbl_summary(dirtablesubpop[[subpop]], table_3, table3_name)
  
  ### Table 4
  print("Now creating: Table 4")
  table_4 <- create_table_characteristic_population(study_pop = study_population,
                                                    covariates = covariates, at_time = "vax1",
                                                    agebands = Agebands_countpersontime)
  table_4 <- tbl_merge(list(table_3, table_4), tab_spanner = F)
  
  table4_name <- "Table 4.1-10 - Characteristics at first COVID-19 vaccination per DAP (each DAP separate table)"
  
  save_tbl_summary(dirtablesubpop[[subpop]], table_4, table4_name)
  
  
  
  ### Table 5
  print("Now creating: Table 5")
  
  # Generate empty list which will contains all the tables to combined
  tbl_list <- list()
  
  # Use as header the name of the DAP with attacched the subpopulation
  header_string <- "**{paste0(thisdatasource, suffix[[subpop]])}**"
  
  # Set the DAP name
  study_population[, DAP := paste0(thisdatasource, suffix[[subpop]])]
  
  ### Total population
  # Select columns of interest and creates a temporary variable to count the total
  tot_pop <- study_population[, .(person_id, total = 1, DAP)]
  
  # Create the table which contains the total population
  tot_pop <- tot_pop %>%
    tbl_summary(label = list(total ~ "Study population"), by = "DAP", include = total, percent = "row") %>%
    modify_header(all_stat_cols(T) ~ header_string) %>%
    modify_footnote(all_stat_cols(FALSE) ~ NA)
  tbl_list <- append(tbl_list, list(tot_pop))
  
  # Cycle for each brand
  for (basic_vx_manufacturer in manufacturer_in_study) {
    ### Percentage of 1^ doses wrt total population
    # Select id, vax1 brand and DAP name
    study_pop_vax1 <- study_population[, .(person_id, type_vax_1, DAP)]
    study_pop_vax1 <- gt_dichotomous(study_pop_vax1, "type_vax_1", basic_vx_manufacturer, "%1$s dose 1")
    tbl_list <- append(tbl_list, list(study_pop_vax1))
    
    # Filter to get only vax1 as basic_vx_manufacturer
    study_pop_vax1 <- copy(study_population)[type_vax_1 == basic_vx_manufacturer, ]
    
    # Cycle for each doses we want the distances and percentage of.
    for (i in seq(2, 3)) {
      
      # Create names of columns which we will need in this cycle
      current_type <- paste("type_vax", i, sep = "_")
      previous_type <- paste("type_vax", i - 1, sep = "_")
      current_date <- paste("date_vax", i, sep = "_")
      previous_date <- paste("date_vax", i - 1, sep = "_")
      
      ### Percentage of 2^ doses wrt basic_vx_manufacturer 1^ doses
      # Type_vax_2 as factors
      study_pop_vax1_n <- study_pop_vax1[, (current_type) := factor(get(current_type), levels = manufacturer_in_study)]
      
      # Select variables of interest
      study_pop_vax1_n <- study_pop_vax1_n[, c("person_id", current_type, "DAP"), with = F]
      
      # Compute tables
      if (nrow(study_pop_vax1_n) > 1) {
        basic_tbl_vax_n <- gt_dichotomous(study_pop_vax1_n, current_type, basic_vx_manufacturer, paste("%1$s dose", i))
        other_tbl_vax_n <- gt_dichotomous(study_pop_vax1_n, current_type, basic_vx_manufacturer, paste("%1$s dose", i),
                                          other = T)
        rest_tbl_vax_n <- gt_dichotomous(study_pop_vax1_n, current_type,
                                         setdiff(manufacturer_in_study, basic_vx_manufacturer), paste("%1$s dose", i))
        
        # Add the tables to the final list
        tbl_list <- append(tbl_list, list(basic_tbl_vax_n, other_tbl_vax_n, rest_tbl_vax_n))
        rm(basic_tbl_vax_n, other_tbl_vax_n, rest_tbl_vax_n)
      }
      
      ### Distances 1^ and 2^ doses homologous wrt basic_vx_manufacturer
      # Compute distance of 1^ and 2^ doses for everyone
      study_pop_vax1[, distance_1_n := difftime(get(current_date), date_vax_1, units = "days")]
      
      # Select homologous doses for basic_vx_manufacturer and variables of interest
      vax_n_distance <- copy(study_pop_vax1)[get(current_type) == basic_vx_manufacturer, ]
      vax_n_distance <- vax_n_distance[, .(person_id, distance_1_n, DAP)]
      
      if (nrow(vax_n_distance) > 1) {
        # Remove warning in case we don't have observed this combination of manufacturer/doses
        defaultW <- getOption("warn") 
        options(warn = -1)
        
        # Table with distances basic statistics
        vax_n_distance <- vax_n_distance %>%
          tbl_summary(label = list(distance_1_n ~ paste("Amongst persons with", basic_vx_manufacturer,
                                                        "dose", i, "distance")),
                      by = "DAP",
                      type = distance_1_n ~ "continuous2",
                      statistic = distance_1_n ~ c("{min}", "{p25}", "{p50}", "{p75}", "{max}"),
                      include = distance_1_n) %>%
          modify_header(all_stat_cols(FALSE) ~ header_string) %>%
          modify_footnote(all_stat_cols(FALSE) ~ NA)
        # Reinstate warnings
        options(warn = defaultW)
        
        # Add the table to the final list
        tbl_list <- append(tbl_list, list(vax_n_distance))
        rm(vax_n_distance)
      }
      
      ### Distances 1^ and 2^ doses heterologous wrt basic_vx_manufacturer
      # Find the other levels wrt basic_vx_manufacturer
      other_levels <- setdiff(manufacturer_in_study, basic_vx_manufacturer)
      
      # Select only observation with vax1 as basic_vx_manufacturer and va2 as other manufacturers
      vax_n_distance_other <- copy(study_pop_vax1)[get(current_type) %in% other_levels, ]
      
      # Retain variables of interest
      vax_n_distance_other <- vax_n_distance_other[, .(person_id, distance_1_n, DAP)]
      
      if (nrow(vax_n_distance_other) > 1) {
        
        # Remove warning in case we don't have observed this combination of manufacturer/doses
        defaultW <- getOption("warn") 
        options(warn = -1)
        
        # Table with distances basic statistics
        vax_n_distance_other <- vax_n_distance_other %>%
          tbl_summary(label = list(distance_1_n ~ paste("Amongst persons with other dose", i, "distance")),
                      by = "DAP",
                      type = distance_1_n ~ "continuous2",
                      statistic = distance_1_n ~ c("{min}", "{p25}", "{p50}", "{p75}", "{max}"),
                      include = distance_1_n) %>%
          modify_header(all_stat_cols(FALSE) ~ header_string) %>%
          modify_footnote(all_stat_cols(FALSE) ~ NA)
        # Reinstate warnings
        options(warn = defaultW)
        
        # Add the table to the final list
        tbl_list <- append(tbl_list, list(vax_n_distance_other))
        rm(vax_n_distance_other)
      }
    }
  }
  
  table_5 <- tbl_stack(tbl_list)
  
  table5_name <- "Table 5 - Covid-19 Vaccination doses and distances"
  
  save_tbl_summary(dirtablesubpop[[subpop]], table_5, table5_name)
  
  
  
  
  
  
  # Load the IR
  load(paste0(dirD4D5subpop[[subpop]], "D5_IR_background", ".RData"))
  RES_IR <- get(paste0("D5_IR_background"))
  rm(list = paste0("D5_IR_background"))
  
  # Load the standardized IR
  load(paste0(dirD4D5subpop[[subpop]], "D5_IR_background_std", ".RData"))
  RES_IR_std <- get(paste0("D5_IR_background_std"))
  rm(list = paste0("D5_IR_background_std"))
  
  ### Table 6
  print("Now creating: Table 6")
  
  # Remove all columns not useful anymore
  RES_IR <- RES_IR[, Persontime := NULL]
  
  # Retain only year 2019/2020
  RES_IR <- RES_IR[year %in% c("2019", "2020"), ]
  
  # Prepare for the melt by creating the list of columns needed
  list_AESI_NCO <- c(OUTCOME_variables, "DEATH", CONTROL_variables)
  
  colA = paste0(list_AESI_NCO, "_b")
  colB = paste0("Persontime_", list_AESI_NCO)
  colC = paste0("IR_", list_AESI_NCO)
  colD = paste0("lb_", list_AESI_NCO)
  colE = paste0("ub_", list_AESI_NCO)
  
  RES_IR <- data.table::melt(RES_IR, measure = list(colA, colB, colC, colD, colE),
                             variable.name = "var", variable.factor = F,
                             value.name = c("counts", "persontime", "IR", "lb", "ub"))
  
  colF = paste0("IR_std_", list_AESI_NCO)
  colG = paste0("lb_std_", list_AESI_NCO)
  colH = paste0("ub_std_", list_AESI_NCO)
  
  RES_IR_std <- data.table::melt(RES_IR_std, measure = list(colF, colG, colH),
                                 variable.name = "var", variable.factor = F,
                                 value.name = c("IR", "lb", "ub"))
  
  # Convert numeric variable to correct variable names
  vect_recode_AESI_NCO <- list_AESI_NCO
  names(vect_recode_AESI_NCO) <- as.character(seq_along(vect_recode_AESI_NCO))
  RES_IR <- RES_IR[ , var := vect_recode_AESI_NCO[var]]
  RES_IR_std <- RES_IR_std[ , var := vect_recode_AESI_NCO[var]]
  
  # Cycle for each AESI/NCO
  for (current_var in list_AESI_NCO) {
    
    # Select current variable
    RES_IR_var <- copy(RES_IR)[var == current_var, ]
    RES_IR_std_var <- copy(RES_IR_std)[var == current_var, ]
    
    # Create list wich will contains the tables
    tbl_list <- list()
    
    # Create vector to cycle
    year_to_select <- c(2019, 2020, 2020)
    covid_to_select <- c(0, 0, 1)
    
    ### Divide for before and after COVID
    for (i in seq_len(3)) {
      
      # Select only PT without COVID
      RES_IR_no_covid <- copy(RES_IR_var)[COVID19 == covid_to_select[i], ][, COVID19 := NULL]
      RES_IR_no_covid <- RES_IR_no_covid[year == year_to_select[i], ][, year := NULL]
      RES_IR_std_no_covid <- copy(RES_IR_std_var)[COVID19 == covid_to_select[i], ][, COVID19 := NULL]
      RES_IR_std_no_covid <- RES_IR_std_no_covid[year == year_to_select[i], ][, year := NULL]
      
      # Select the total for sex and age
      total_RES_IR_no_covid <- RES_IR_no_covid[sex == "total" & Ageband == "total", ][, c("sex", "Ageband") := NULL]
      
      total_RES_IR_no_covid <- data.table::melt(total_RES_IR_no_covid,
                                                measure = list(c("counts", "persontime", "IR", "lb", "ub")),
                                                variable.name = "statistic", variable.factor = F,
                                                value.name = c("value"))
      
      total_RES_IR_no_covid <- tbl_PT_IR_dichotomous(total_RES_IR_no_covid, "Total (all ages/gender)")
      
      RES_IR_std_no_covid[, names(RES_IR_std_no_covid) := lapply(.SD, as.character)]
      RES_IR_std_no_covid <- RES_IR_std_no_covid[, c("counts", "persontime") := list("", "")]
      total_RES_IR_std_no_covid <- data.table::melt(RES_IR_std_no_covid,
                                                    measure = list(c("counts", "persontime", "IR", "lb", "ub")),
                                                    variable.name = "statistic", variable.factor = F,
                                                    value.name = c("value"))
      
      total_RES_IR_std_no_covid <- tbl_PT_IR_dichotomous(total_RES_IR_std_no_covid, "Total age Standardized")
      
      # Select the total for age and remove the total for sex
      sex_RES_IR_no_covid <- RES_IR_no_covid[sex != "total" & Ageband == "total", ][, Ageband := NULL]
      
      sex_RES_IR_no_covid <- data.table::melt(sex_RES_IR_no_covid,
                                              measure = list(c("counts", "persontime", "IR", "lb", "ub")),
                                              variable.name = "statistic", variable.factor = F,
                                              value.name = c("value"))
      
      sex_RES_IR_no_covid <- tbl_PT_IR_categorical(sex_RES_IR_no_covid, "Gender specific", "sex")
      
      # Select the total for sex and remove the total for ageband
      age_RES_IR_no_covid <- RES_IR_no_covid[sex == "total" & Ageband != "total", ][, sex := NULL]
      
      age_RES_IR_no_covid <- data.table::melt(age_RES_IR_no_covid,
                                              measure = list(c("counts", "persontime", "IR", "lb", "ub")),
                                              variable.name = "statistic", variable.factor = F,
                                              value.name = c("value"))
      age_RES_IR_no_covid[, Ageband := factor(Ageband, levels = Agebands_labels)]
      
      age_RES_IR_no_covid <- tbl_PT_IR_categorical(age_RES_IR_no_covid, "Age specific (years)", "Ageband")
      
      tbl_list <- append(tbl_list, list(tbl_stack(list(total_RES_IR_no_covid, total_RES_IR_std_no_covid,
                                                       sex_RES_IR_no_covid, age_RES_IR_no_covid))))
    }
    
    table_6 <- tbl_stack(tbl_list, group_header = c("2019 before COVID-19 infection",
                                                    "2020 before COVID-19 infection",
                                                    "2020 after COVID-19 infection before vaccination"))
    
    table6_name <- paste0("Table 6.1-", which(list_AESI_NCO == current_var),
                          " - Background Incidence rates of ", current_var)
    
    save_tbl_summary(dirtablesubpop[[subpop]], table_6, table6_name, "table 6/")
  }
}
# # Covariates to a single column
# pop_monthly <- data.table::melt(pop_covariates, measure = c(risk_factors, medicines, pregnancy),
#                                 variable.name = "has_covariate", variable.factor = F,
#                                 value.name = c("covariate"))
# pop_monthly <- pop_monthly[, fifelse(covariate == 1, covariate:)]
