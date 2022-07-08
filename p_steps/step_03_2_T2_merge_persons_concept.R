# -----------------------------------------------------
# Merge persons and concepts

# input: output_spells_category, D3_selection_criteria, D3_concepts_QC_criteria
# output: persons_doses

#merge vaccine information to persons
 for (subpop in subpopulations_non_empty){
   print(subpop)
   
   if (this_datasource_has_subpopulations == T){
     load(paste0(dirtemp,"output_spells_category.RData"))
     output_spells_category<-output_spells_category[[subpop]]
   }else{
     load(paste0(dirtemp,"output_spells_category.RData"))
   }
  
  load(paste0(dirtemp,"D3_selection_criteria", suffix[[subpop]] ,".RData"))
  selection_criteria<-get(paste0("D3_selection_criteria", suffix[[subpop]]))
  
  load(paste0(dirtemp,"D3_concepts_QC_criteria.RData"))

persons_doses<-merge(selection_criteria,D3_concepts_QC_criteria, by=c("person_id"),all=T)

persons_doses<-persons_doses[is.na(sex_or_birth_date_missing),sex_or_birth_date_missing:=1]

temp <- selection_criteria[, .(person_id, date_of_death)]
temp1 <- D3_concepts_QC_criteria[, .(person_id, date, vx_dose)]
names(output_spells_category)
temp2 <- output_spells_category[, .(person_id, entry_spell_category, exit_spell_category)]
temp2 <- temp2[study_start %between% list(entry_spell_category,exit_spell_category) & entry_spell_category < exit_spell_category, ]

temp_tot <- merge(temp, temp1)
temp_tot <- merge(temp_tot, temp2)
temp_tot <- temp_tot[date_of_death < date, death_before_vax := 1]
temp_tot <- temp_tot[exit_spell_category < date, exit_spell_before_vax := 1]
temp_tot <- temp_tot[study_end < date, study_end_before_vax := 1]
temp_tot <- temp_tot[, .(person_id, vx_dose, death_before_vax, exit_spell_before_vax, study_end_before_vax)]

for (i in names(temp_tot)){
  temp_tot[is.na(get(i)), (i):=0]
}

persons_doses <- merge(persons_doses, unique(temp_tot), by = c("person_id", "vx_dose"),all.x=T)


assign(paste0("persons_doses",suffix[[subpop]]), persons_doses)
save(list=paste0("persons_doses",suffix[[subpop]]),file=paste0(dirtemp,"persons_doses",suffix[[subpop]],".RData"))
rm(list=paste0("persons_doses",suffix[[subpop]]))
rm(list=paste0("D3_selection_criteria", suffix[[subpop]]))
}

if(this_datasource_has_subpopulations==T) rm(persons_doses)
rm(selection_criteria, D3_concepts_QC_criteria , temp, temp1, output_spells_category, temp2, temp_tot)
