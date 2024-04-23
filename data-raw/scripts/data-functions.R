# Define coordinate reference system used for the project

#project_crs <- "EPSG:3310" # NAD83 California Teale Albers
#project_crs <- "EPSG:5070" # NAD83 CONUS Albers
project_crs <- "ESRI:102039" # NAD83 CONUS Albers USGS Version

# Necessary functions for data aggregation and exploration

library(googledrive)

drive_file_by_id <- function(id=character(), dir=here::here("data-raw", "temp"), vsizip=FALSE) {
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

# Linear interpolation function used in multiple places

linterp <- function(x, x1, x2, y1, y2){
  y1 + ((x-x1)/(x2-x1)) * (y2-y1)
}


# Spatial processing functions
# perpendicular transects used in both valley width and levees scripts

library(sf)

perpendicular_transect <- function(pair=st_multipoint(), length=numeric()) {
  x1 <- pair[1]
  x2 <- pair[2]
  y1 <- pair[3]
  y2 <- pair[4]
  midpoint_x <- (x1+x2)/2
  midpoint_y <- (y1+y2)/2
  orig_bearing <- atan2((y2-y1), (x2-x1))
  perp_bearing <- orig_bearing + pi/2
  radius <- length/2
  perp_x1 <- midpoint_x + radius*cos(perp_bearing)
  perp_y1 <- midpoint_y + radius*sin(perp_bearing)
  perp_x2 <- midpoint_x - radius*cos(perp_bearing)
  perp_y2 <- midpoint_y - radius*sin(perp_bearing)
  return(st_linestring(c(st_point(c(perp_x1,perp_y1)),
                         st_point(c(perp_x2, perp_y2)))))
}

perpendicular_transects <- function(ls=st_linestring(), length=numeric()) {
  if ("sf" %in% class(ls)){
    ls <- st_as_sfc(ls, crs=st_crs(ls))
  }
  if ("sfc" %in% class(ls)){
    m <- t(st_zm(ls[[1]]))
  }
  if ("sfg" %in% class(ls)){
    m <- t(ls)
  }
  point_pairs <- lapply(seq(1, length(m)-3, 2), function(i) c(st_point(c(m[i], m[i+1])), st_point(c(m[i+2], m[i+3]))))
  perp_lines <- lapply(point_pairs, function(pair) perpendicular_transect(pair, length))
  return(st_sfc(perp_lines, crs=st_crs(ls)))
}
