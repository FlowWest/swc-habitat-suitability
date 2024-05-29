#' @name const_proj_crs
#' @title Project CRS
#' @examples
#' sf::st_sfc(crs="+proj=latlong +datum=WGS84") |>
#'   sf::st_transform(const_proj_crs())
#' sf::st_crs(const_proj_crs())
#' @export
const_proj_crs <- function() "ESRI:102039" # NAD83 CONUS Albers USGS Version

#' Download Data from Google Drive
#'
#' @param id The character ID from a Google Drive file. Found in the middle of the share link: `https://drive.google.com/file/d/XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/view`
#' @param dir Path to the local directory in which to download the file
#' @param vsizip Whether to append the "vsizip" parameter which allows GEOS access within ZIPs without extracting first
#'
#' @returns The path to the downloaded file
#' @md
#'
#' @export
drive_file_by_id <- function(id=character(), dir=here::here("data-raw", "temp"), vsizip=F) {
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

#' Simple Linear Interpolation
#'
#' @param x The index at which to predict `y`
#' @param x1 An index `x1 <= x` with known value `y1`
#' @param x2 An index `x1 >= x` with known value `y2`
#' @param y1 The value at index `x1`
#' @param y2 The value at index `x2`
#'
#' @returns The predicted value `y` at index `x`
#' @md
#'
#' @export
linterp <- function(x, x1, x2, y1, y2){
  y1 + ((x-x1)/(x2-x1)) * (y2-y1)
}

#' Perpendicular Transect (helper function for perpendicular_transects)
#'
#' @param pair An `sf` multipoint containing two points along a stream or profile line
#' @param length The length of transect to create
#'
#' @returns An `sf` linestring transect, halfway between the two input points in `pair`, angled perpendicular to their alignment, measuring length `length`.
#' @md
#'
#' @export
perpendicular_transect <- function(pair=sf::st_multipoint(), length=numeric()) {
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
  return(sf::st_linestring(c(sf::st_point(c(perp_x1,perp_y1)),
                             sf::st_point(c(perp_x2, perp_y2)))))
}

#' Create Perpendicular Transects from Line
#'
#' @param ls An `sfc` linestring geometry collection defining a stream or profile line along which to generate transects.
#' @param length The length of each transect to create, which typically should be wide enough to span the streamway or valley in question.
#'
#' @returns An `sfc` linestring geometry collection of perpendicular transect lines, one for each pair of subsequent points in `ls`
#' @md
#'
#' @export
perpendicular_transects <- function(ls=sf::st_linestring(), length=numeric()) {
  if ("sf" %in% class(ls)){
    ls <- sf::st_as_sfc(ls, crs=sf::st_crs(ls))
  }
  if ("sfc" %in% class(ls)){
    m <- t(sf::st_zm(ls[[1]]))
  }
  if ("sfg" %in% class(ls)){
    m <- t(ls)
  }
  point_pairs <- lapply(seq(1, length(m)-3, 2),
                        function(i) c(sf::st_point(c(m[i], m[i+1])),
                                      sf::st_point(c(m[i+2], m[i+3]))))
  perp_lines <- lapply(point_pairs,
                       function(pair) perpendicular_transect(pair, length))
  return(sf::st_sfc(perp_lines, crs=sf::st_crs(ls)))
}

#' Cumulative Mean ignoring NA values
#'
#' @param x An integer or numeric vector
#'
#' @return A vector the same length as `x`
#' @md
#' @examples
#' v <- c(1, 2, 3, NA, NA, 1, 2, 6, NA)
#' cummean.na.rm(v)
#'
#' @export
cummean.na.rm <- function(x) {
  idx <- cumsum(!is.na(x))
  x_filtered <- x[!is.na(x)]
  dplyr::cummean(x_filtered)[idx]
}

#' Semi Inverse Hyperbolic Sine Transform
#'
#' @export
semiIHS <- function(x) asinh(x / 2)

#' Semi Inverse Hyperbolic Sine Transform (Inverse)
#'
#' @export
semiIHS_inv <- function(y) 2 * sinh(y)

#' Semi Inverse Hyperbolic Sine Transform * 100
#'
#' @export
semiIHS00 <- function(x) asinh(x * 100 / 2)

#' Semi Inverse Hyperbolic Sine Transform * 100 (Inverse)
#'
#' @export
semiIHS00_inv <- function(y) 2 * sinh(y) / 100


#' Semi Inverse Hyperbolic Sine Transform (scales transformation for ggplot)
#'
#' @export
trans_semiIHS <- scales::trans_new("semiIHS",
                                   transform = function(x) asinh(x / 2),
                                   inverse = function(y) 2 * sinh(y),
                                   breaks = function(i) scales::breaks_log(n=5, base=10)(pmax(i,0.01)),
                                   #minor_breaks = scales::minor_breaks_n(n = 0),
                                   domain=c(0, Inf),
                                   format = scales::label_comma())

#' Order of Magnitude Range
#'
#' @param from An integer defining the start of the sequence
#' @param to An integer defining the end of the sequence
#'
#' @returns An integer vector sequence
#' @md
#'
#' @examples
#' oom_range(9, 47000)
#'
#' @export
oom_range <- function(from, to) {
  from_mag <- floor(log10(if (from > 0) from else 1))
  to_mag <- ceiling(log10(to))
  magnitudes <- seq(from_mag, to_mag + 1, 1)
  expanded <- as.vector(t(10^seq(magnitudes) %o% seq(1,9,1)))
  result <- expanded[which(expanded>=signif(from,1) & expanded<=signif(to,1))]
  return(if (from == 0) c(0, result) else result)
}

#' Hex Color Blend
#'
#' @param x A character string describing a color in the form `"#000000"`
#' @param y A character string describing a color in the form `"#000000"`
#'
#' @returns A character string describing a color in the form `"#000000"`
#' @md
#'
#' @examples
#' hex_color_blend("#55ad89", "#c3bc3f")
#'
#' @export
hex_color_blend <- function(x, y){
  hex_r <- as.character.hexmode((strtoi(substr(x,2,3), base=16L) + strtoi(substr(y,2,3), base=16L)) / 2)
  hex_g <- as.character.hexmode((strtoi(substr(x,4,5), base=16L) + strtoi(substr(y,4,5), base=16L)) / 2)
  hex_b <- as.character.hexmode((strtoi(substr(x,6,7), base=16L) + strtoi(substr(y,6,7), base=16L)) / 2)
  return(paste0("#", hex_r, hex_g, hex_b))
}
