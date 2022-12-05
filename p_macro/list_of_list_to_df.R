list_of_list_to_df <- function(l) {
  
  out <- list()
  for (j in names(l)) {
    for (i in names(l[[j]])) {
      out <- append(out, list(data.table(code = l[[j]][[i]], coding_system = i, event_abbreviation = j)))
    }
  }
  
  out <- data.table::rbindlist(out)
  
  if (nrow(out[grepl("_narrow|_possible", event_abbreviation), ]) > 1) {
    out[grepl("_narrow|_possible", event_abbreviation),
        c("event_abbreviation", "tags") := tstrsplit(event_abbreviation, "_(?!.*_)", perl = T, type.convert = T)]
  } else {
    out[grepl("_narrow|_possible", event_abbreviation),
        tags := NA_character_]
  }
  
  return(out)
}