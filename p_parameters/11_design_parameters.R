manufacturer_in_study <- covid_vaccines_ConcePTION_CDM_vocabulary[covid_vaccines_ConcePTION_CDM_vocabulary %in%
                                                                    c("pfizer", "moderna", "astrazeneca", "janssen", "novavax")]
manufacturer_in_study <- c(manufacturer_in_study, "unk")

max_number_doses <- 4

# for debugging purpose set a predefined subpop (it does not affect the script results)
subpop <- subpopulations_non_empty[1]

datasources_SCRI_SCCS_COHORT <- c("TEST", "ARS", "BIFAP", "FISABIO", "SIDIAP", "CPRD")
