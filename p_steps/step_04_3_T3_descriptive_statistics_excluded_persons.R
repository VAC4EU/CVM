#---------------------------------------------------------------
# Descriptive statistics on exclusion criteria

# input: D3_selection_criteria, output_spells_category
# output: pop_excluded, number_criteria_excluded, plots

load(paste0(dirtemp,"D3_selection_criteria.RData"))

base_pop <- D3_selection_criteria[sex_or_birth_date_missing == 0 & birth_date_absurd == 0 & 
                                    no_observation_period == 0 & death_before_study_entry == 0, ]

base_pop <- base_pop[, date_of_birth := findInterval(year(date_of_birth), c(1940, 1950, 1960, 1970, 1980, 1990))]
base_pop$date_of_birth <- as.character(base_pop$date_of_birth)
base_pop <- base_pop[.(date_of_birth = c("0", "1", "2", "3", "4", "5", "6"),
                       to = c("<1940", "1940-1949", "1950-1959", "1960-1969", "1970-1979", "1980-1989", "1990+")),
                     on = "date_of_birth", date_of_birth := i.to]

population_ins_run_in <- copy(base_pop)[no_observation_period_including_study_start == 0 & insufficient_run_in == 1, ]
population_ins_run_in <- population_ins_run_in[, .(person_id, sex, date_of_birth)][ , variable := "ins_run_in"]

population_no_obs_start <- copy(base_pop)[no_observation_period_including_study_start == 1, ]
population_no_obs_start <- population_no_obs_start[, .(person_id, sex, date_of_birth)][ , variable := "no_obs_in_start"]

population_discarded <- unique(rbind(population_ins_run_in, population_no_obs_start)[ , variable := NULL])

population_filtered <- copy(base_pop)[no_observation_period_including_study_start == 0 & insufficient_run_in == 0, ]
population_filtered <- population_filtered[, .(person_id, sex, date_of_birth)][ , variable := "filtered_pop"]

rm(base_pop)

population_all <- rbind(population_ins_run_in, population_no_obs_start, population_filtered)

population_sex <- copy(population_all)[, .N, by = c("sex", "variable")]
population_sex <- population_sex[, value := paste0(N, " (", round(N / sum(N) * 100, 1), "%)"), by = c("variable")]
population_date <- population_all[, .N, by = c("date_of_birth", "variable")]
population_date <- population_date[, value := paste0(N, " (", round(N / sum(N) * 100, 1), "%)"), by = c("variable")]

population_sex <- dcast(population_sex, sex ~ variable, value.var = "value")
population_date <- dcast(population_date, date_of_birth ~ variable, value.var = "value")

for (i in names(population_sex)){population_sex[is.na(get(i)), (i) := "0 (0%)"]}
for (i in names(population_date)){population_date[is.na(get(i)), (i) := "0 (0%)"]}

setnames(population_sex, "sex", "")
setnames(population_date, "date_of_birth", "")

pop_excluded <- rbind(population_sex, population_date)

to_ord <- c(1, 2, 4, 5, 6, 7, 8, 9, 3)
setorder(pop_excluded[, .r := to_ord], .r)[, .r := NULL]

fwrite(pop_excluded, file = paste0(direxp, "pop_excluded.csv"))

rm(D3_selection_criteria, population_ins_run_in, population_no_obs_start, population_sex,
   population_date, pop_excluded, population_all, population_filtered)

load(paste0(dirtemp,"output_spells_category.RData"))

out_spells <- merge(population_discarded, output_spells_category[, op_meaning := NULL], by = "person_id")
number_criteria_excluded <- data.table::data.table(incorrect_spell = copy(out_spells)[entry_spell_category > exit_spell_category, .N])
out_spells <- out_spells[entry_spell_category < exit_spell_category, ]
out_spells <- out_spells[, after_20200101 := fifelse(study_start <= entry_spell_category |
                                                       study_start <= exit_spell_category, 1, 0)]
first2020 <- out_spells[out_spells[after_20200101 == 1, .I[which.min(num_spell)], by = "person_id"]$V1]
last2019 <- out_spells[out_spells[after_20200101 == 0, .I[which.max(num_spell)], by = "person_id"]$V1]

out_spells <- data.table::rbindlist(list(first2020, last2019))
out_spells <- out_spells[, N_spell := .N, by = "person_id"]
number_criteria_excluded <- number_criteria_excluded[, single_spell := copy(out_spells)[N_spell == 1, .N]]
number_criteria_excluded <- number_criteria_excluded[, double_spell := copy(out_spells)[N_spell == 2, .N]]

fwrite(number_criteria_excluded, file = paste0(direxp, "number_criteria_excluded.csv"))

out_spells_start_2019 <- copy(out_spells)[exit_spell_category < study_start, ]
out_spells_start_2019 <- out_spells_start_2019[, exit_spell_category := as.Date(exit_spell_category)]
fifth_percentile <- as.Date(out_spells_start_2019[, quantile(exit_spell_category, probs = 0.05, type = 1)[[1]]] + 5)
out_spells_start_2019 <- out_spells_start_2019[exit_spell_category >= fifth_percentile, ]

g <- ggplot(out_spells_start_2019, aes(exit_spell_category)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = 1) +
  geom_density(alpha=0.2) +
  labs(title = "Density plot of spell end date", 
       subtitle = "End date of last spell <2020",
       caption = "Source: out_spells_start_2019",
       x = "Spell end date")
suppressMessages(ggsave(paste0(direxp, "Last spells end date.png"), plot = g,
                        units = c("cm"), dpi = 600))

out_spells <- out_spells[N_spell == 2, ][, .(entry_spell_category = max(entry_spell_category),
                                             exit_spell_category = min(exit_spell_category)), by = "person_id"]
out_spells <- out_spells[, distance := as.numeric(entry_spell_category - exit_spell_category, units="days")]

g <- ggplot(out_spells, aes(distance)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = 1) +
  geom_density(alpha=0.2) +
  scale_x_continuous(limits = c(out_spells[, min(distance)] - 5,
                                out_spells[, quantile(distance, probs = 0.95)] + 5),
                     oob = scales::oob_keep) +
  labs(title = "Density plot of distances", 
       subtitle = "Distances between last spell <2020 and first spell >=2020",
       caption = "Source: out_spells",
       x = "Distance")
suppressMessages(ggsave(paste0(direxp, "Density_plot_distance_spells.png"), plot = g,
                        units = c("cm"), dpi = 600))