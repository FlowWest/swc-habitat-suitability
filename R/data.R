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

#' @name wua_hydraulic_interp
#' @title WUA/LF by ComID by Flow (training dataset, interpolated)
#' @examples
#' wua_hydraulic_interp |> pillar::glimpse()
"wua_hydraulic"

#' @name wua_predicted
#' @title WUA/LF by ComID by Flow (predicted by statistical model)
#' @examples
#' wua_predicted |> pillar::glimpse()
"wua_predicted"

#' @name cv_mainstems
#' @title CVPIA Mainstems (geometry)
#' @examples
#' cv_mainstems |> ggplot2::ggplot() + ggplot2::geom_sf()
"cv_mainstems"

#' @name cv_watersheds
#' @title CVPIA Watersheds (geometry)
#' @examples
#' cv_watersheds |> ggplot2::ggplot() + ggplot2::geom_sf()
"cv_watersheds"

#' @name streamgage_attr
#' @title CDEC Mainstem Streamgage Attributes
"streamgage_attr"

#' @name streamgage_duration_rating_curves
#' @title CDEC Mainstem Streamgage Duration Rating Curves
"streamgage_duration_rating_curves"
