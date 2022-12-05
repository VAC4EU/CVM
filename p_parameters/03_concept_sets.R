# the output of this step are several lists and lists of lists to be used by CreateConceptsetDatasets

# concept_set_domains
# level 1 name of the conceptset (including in the name whether it is narrow or possible)
# level 2 domain

# concept_set_codes_our_study,
# level 1 name of the conceptset (including in the name whether it is narrow or possible)
# level 2 coding system
# level 3 list of codes incuded in that conceptset for that coding system

# concept_set_codes_our_study_excl
# level 1 name of the conceptset (including in the name whether it is narrow or possible)
# level 2 coding system
# level 3 list of codes to be excluded from that conceptset for that coding system

# the lists of all conceptsets 
# concept_sets_of_our_study

# input: the VAC4EU spreadsheets, restricted to the conceptsets associated with this study

### IN CASE A CONCEPT IS TOO BIG IN A DAP
datasource_needing_split_conceptsets <- c("CPRD", "TEST", "ARS")
CONCEPTSETS_to_be_split <- if(thisdatasource %in% datasource_needing_split_conceptsets) c("DP_COVCARDIOCEREBROVAS") else c()
numbers_split <- c(15)

OUT_codelist <- fread(paste0(thisdir,"/p_parameters/archive_parameters/20221004_V2_ALL_full_codelist.csv"))
OUT_codelist <- OUT_codelist[, .(coding_system, code, type, tags, event_abbreviation, system)]
OUT_codelist <- OUT_codelist[, Varname := paste(system, event_abbreviation, type, sep = "_")]
File_variables_ALG_DP_ROC20 <- paste0(thisdir,"/p_parameters/archive_parameters/Variables_ALG_DP_ROC20_July22.xlsx")

VAR_list <- as.data.table(readxl::read_excel(File_variables_ALG_DP_ROC20, sheet = "Variables"))[, .(Varname)]

# Adding I_COVID19DX_AESI manually
VAR_list <- rbindlist(list(VAR_list, data.table(Varname = "I_COVID19DX_AESI")))
OUT_codelist <- merge(VAR_list, OUT_codelist, all.x = T, by = "Varname")
rm(VAR_list)

# TODO ok for release?
OUT_codelist <- OUT_codelist[code != "" & !is.na(code), ][, event_abbreviation := toupper(event_abbreviation)]
OUT_codelist <- OUT_codelist[tags != ""][tags == "possbie", tags := "possible"]
OUT_codelist <- OUT_codelist[coding_system %not in% c("MEDCODEID", "MedCodeId")]

concept_set_codes_our_study <- df_to_list_of_list(OUT_codelist, codying_system_recode = "auto", type_col = "type")
rm(OUT_codelist)

concept_set_domains<- vector(mode="list")
for (concept in names(concept_set_codes_our_study)) {
  concept_set_domains[[concept]] = "Diagnosis"
}

DRUG_codelist <- as.data.table(readxl::read_excel(File_variables_ALG_DP_ROC20, sheet = "DrugProxies",
                                                  .name_repair = ~ vctrs::vec_as_names(..., repair = "universal", quiet = TRUE)))

DRUG_codelist <- DRUG_codelist[, ATC.codes := strsplit(ATC.codes, ",")]

DRUG_codelist_list <- df_to_list_of_list(DRUG_codelist, code_col = "ATC.codes", concepts_col = "Drug_proxie",
                                    codying_system_col = F, codying_system_recode = "auto")

for (concept in names(DRUG_codelist_list)) {
  concept_set_domains[[concept]] = "Medicines"
}

concept_set_codes_our_study <- c(concept_set_codes_our_study, DRUG_codelist_list)

if (!is.null(CONCEPTSETS_to_be_split)) {
  
  # Select only concepts to be split
  list_CONCEPTSETS_to_be_split <- concept_set_codes_our_study[names(concept_set_codes_our_study) %in% CONCEPTSETS_to_be_split]
  
  # Remove the above mentioned concept from the general list
  concept_set_codes_our_study <- concept_set_codes_our_study[names(concept_set_codes_our_study) %not in% CONCEPTSETS_to_be_split]
  
  df_CONCEPTSETS_to_be_split <- list_of_list_to_df(list_CONCEPTSETS_to_be_split)
  df_CONCEPTSETS_to_be_split[, group := cut(seq_along(code), numbers_split, labels = FALSE),
                             by = c("coding_system", "event_abbreviation", "tags")]
  
  df_CONCEPTSETS_to_be_split_NA <- df_CONCEPTSETS_to_be_split[is.na(tags), ]
  df_CONCEPTSETS_to_be_split_not_NA <- df_CONCEPTSETS_to_be_split[!is.na(tags), ]
  df_CONCEPTSETS_to_be_split_NA <- df_CONCEPTSETS_to_be_split_NA[CJ(coding_system = coding_system,
                                                                    event_abbreviation = event_abbreviation, tags = tags,
                                                                    group = seq_len(numbers_split), unique = TRUE),
                                                                 on = .(coding_system, event_abbreviation, tags, group)]
  df_CONCEPTSETS_to_be_split_not_NA <- df_CONCEPTSETS_to_be_split_not_NA[CJ(coding_system = coding_system,
                                                                            event_abbreviation = event_abbreviation, tags = tags,
                                                                            group = seq_len(numbers_split), unique = TRUE),
                                                                         on = .(coding_system, event_abbreviation, tags, group)]
  df_CONCEPTSETS_to_be_split <- rbind(df_CONCEPTSETS_to_be_split_not_NA, df_CONCEPTSETS_to_be_split_NA)
  
  df_CONCEPTSETS_to_be_split[!is.na(tags), event_abbreviation := paste(event_abbreviation, tags, group, sep = "_")]
  df_CONCEPTSETS_to_be_split[is.na(tags), event_abbreviation := paste(event_abbreviation, group, sep = "_")]
  df_CONCEPTSETS_to_be_split[, c("group", "tags") := NULL]
  
  list_CONCEPTSETS_splitted <- lapply(split(df_CONCEPTSETS_to_be_split, by = "event_abbreviation", keep.by = F),
                                      split, by = "coding_system", keep.by = F)
  
  list_CONCEPTSETS_splitted <- lapply(list_CONCEPTSETS_splitted, sapply, unlist, use.names = F, simplify = F)
  
  concept_set_codes_our_study <- c(concept_set_codes_our_study, list_CONCEPTSETS_splitted)
  
  
  # Select only concepts to be split
  list_concept_set_domains_to_be_split <- concept_set_domains[names(concept_set_domains) %in% CONCEPTSETS_to_be_split]
  
  # Remove the above mentioned concept from the general list
  concept_set_domains <- concept_set_domains[names(concept_set_domains) %not in% CONCEPTSETS_to_be_split]
  
  df_concept_set_domains_to_be_split <- list()
  for (i in names(list_concept_set_domains_to_be_split)) {
    df_concept_set_domains_to_be_split <- append(df_concept_set_domains_to_be_split,
                                                 list(data.table(domain = list_concept_set_domains_to_be_split[[i]], concept = i)))
  }
  df_concept_set_domains_to_be_split <- data.table::rbindlist(df_concept_set_domains_to_be_split)
  
  df_concept_set_domains_to_be_split <- df_concept_set_domains_to_be_split[CJ(concept = concept,
                                                                              group = seq_len(numbers_split),
                                                                              unique = TRUE),
                                                                           on = .(concept)]
  df_concept_set_domains_to_be_split[, concept := paste(concept, group, sep = "_")]
  df_concept_set_domains_to_be_split <- df_concept_set_domains_to_be_split[, group := NULL]
  list_concept_set_domains_splitted <- lapply(split(df_concept_set_domains_to_be_split, by = "concept", keep.by = F),
                                              unlist, use.names = F)
    
  concept_set_domains <- c(concept_set_domains, list_concept_set_domains_splitted)
  
}

concept_sets_of_our_study <- names(concept_set_codes_our_study)
rm(DRUG_codelist_list)

concept_set_codes_our_study[["COVID_VACCINES"]][["ATC"]] <- c("J07BX03")
concept_set_domains[["COVID_VACCINES"]] = "VaccineATC"
vaccine_conceptssets <- c("COVID_VACCINES")

concept_set_codes_our_study_excl <- vector(mode="list")
concept_set_codes_our_study_excl[["DP_VACCINES"]] <- concept_set_codes_our_study[["COVID_VACCINES"]]

rm(concept)
