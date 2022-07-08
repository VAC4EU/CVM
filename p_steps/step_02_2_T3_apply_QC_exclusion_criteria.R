# -----------------------------------------------------------------------
# Apply Quality Check criteria for vaccination records

# input: D3_concepts_QC_criteria
# output: Flowchart_QC_criteria, selected_doses

load(paste0(dirtemp, "D3_concepts_QC_criteria.RData"))

selected_doses<- CreateFlowChart(
  dataset = D3_concepts_QC_criteria,
  listcriteria = c("duplicated_records", "missing_date", "date_before_start_vax", "distance_btw_1_2_doses",
                   "distance_btw_2_3_doses", "dose_after_3", "dose_after_2"),
  flowchartname = "Flowchart_QC_criteria")

suppressWarnings(
  if(!(file.exists(direxp))) {dir.create(file.path(direxp))}
)
suppressWarnings(
  if(!(file.exists(diroutput))) {dir.create(file.path(diroutput))}
)

fwrite(Flowchart_QC_criteria, paste0(direxp,"Flowchart_QC_criteria.csv"))

save(selected_doses, file = paste0(dirtemp, "selected_doses.RData"))

# rm(PERSONS, OBSERVATION_PERIODS)
rm(selected_doses, D3_concepts_QC_criteria, Flowchart_QC_criteria)
