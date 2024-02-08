# Necessary functions for data aggregation and exploration

library(googledrive)

drive_file_by_id <- function(id=character(), dir="temp", vsizip=FALSE) {
  d <- googledrive::as_dribble(googledrive::as_id(id))
  p <- file.path(dir, d$name)
  if(!file.exists(p)){
    message(paste0('downloading ', p, '...'))
    googledrive::drive_download(file = d, path = p)
  } else {
    message(paste0(p, ' already exists and will be used...'))
  }

  if(vsizip){
    return(file.path("/vsizip",p))
  } else {
    return(p)
  }
}
