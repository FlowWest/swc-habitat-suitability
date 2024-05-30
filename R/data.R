#' @name flowline_geom
#' @title NHD Flowline (ComID) Geometries
#' @examples
#' flowline_geom |> pillar::glimpse()
"flowline_geom"

#' @name flowline_geom_proj
#' @title NHD Flowline (ComID) Geometries - Projected
#' @examples
#' flowline_geom_proj |> pillar::glimpse()
"flowline_geom_proj"

#' @name flowline_attr
#' @title NHD Flowline (ComID) Attributes
#' @examples
#' flowline_attr |> pillar::glimpse()
"flowline_attr"

#' @name wua_hydraulic
#' @title WUA/LF by ComID by Flow (training dataset)
#' @examples
#' wua_hydraulic |> pillar::glimpse()
"wua_hydraulic"

#' @name wua_predicted
#' @title WUA/LF by ComID by Flow (predicted by statistical model)
#' @examples
#' wua_predicted |> pillar::glimpse()
"wua_predicted"
