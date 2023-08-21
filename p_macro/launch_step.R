launch_step <- function(x, print.eval = F) {
  
  if ("dirpargen" %in% ls(envir = .GlobalEnv)) {
    dirpargen <- get("dirpargen", envir = .GlobalEnv)
  } else {
    stop("Folder of the generated parameters is not called 'dirpargen'")
  }
  
  exec_time <- system.time(source(paste0(thisdir,"/", x), print.eval = print.eval))
  print(exec_time)
  
  rm(list=ls(envir = .GlobalEnv), envir = .GlobalEnv)
  load(paste0(dirpargen, "parameters.RData"), envir = .GlobalEnv)
  rm(list=ls(), inherits = T)
  
}