#' NHD Flowline (ComID) Geometries
#'
#' @examples
#' flowline_geom |> pillar::glimpse()
#' @export
"flowline_geom"

#' NHD Flowline (ComID) Geometries - Projected
#'
#' @examples
#' flowline_geom_proj |> pillar::glimpse()
#' @export
"flowline_geom_proj"

#' NHD Flowline (ComID) Attributes
#'
#' @examples
#' flowline_attr |> pillar::glimpse()
#' @export
"flowline_attr"

#' WUA/LF by ComID by Flow (training dataset)
#'
#' @examples
#' wua_hydraulic |> pillar::glimpse()
#' @export
"wua_hydraulic"

#' WUA/LF by ComID by Flow (predicted by statistical model)
#'
#' @examples
#' wua_predicted |> pillar::glimpse()
#' @export
"wua_predicted"
