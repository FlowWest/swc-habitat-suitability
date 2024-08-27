#' Prep depth and velocity raster stacks
#'
#' @param filenames A data frame or tibble containing the file path to each depth and velocity raster at each flow. There should be a numeric column called `flow_cfs`, and two character columns called `depth` and `velocity` each containing valid file paths to rasters (e.g. `.tif` or `.vrt` format).
#'
#' @return A named list of two `terra` `SpatRaster` raster stacks (`depth` and `velocity`), each containing depth or velocity grids at each modeled flow.
#' @md
#' @export
#'
#' @examples
raster_prep_grid <- function(filenames){

  dep_r <- terra::rast(filenames$depth)
  terra::set.names(dep_r, filenames$flow_cfs)

  vel_r <- terra::rast(filenames$velocity)
  terra::set.names(vel_r, filenames$flow_cfs)

  return(list("depth" = dep_r,
              "velocity" = vel_r))

}

#' Raster Habitat Suitability Index (HSI) function: HQT
#'
#' @description
#'
#' @param d a `terra` raster layer for depth in feet
#' @param v a `terra` raster layer for velocity in feet per second
#'
#' @return A `terra` raster layer of habitat suitability index values between 0 and 1.
#' @md
#' @export
#'
#' @examples

raster_dvhsi_hqt <- function(d, v) {
  #(d>1.0 & d<=3.28) & (v>0 & v<=1.5)
  terra::lapp(c(d, v), raster_dvhsi_hqt)
}

#' Raster Habitat Suitability Index (HSI) function: LYR
#'
#' @description
#'
#' @param d a `terra` raster layer for depth in feet
#' @param v a `terra` raster layer for velocity in feet per second
#'
#' @return A `terra` raster layer of habitat suitability index values between 0 and 1.
#' @md
#' @export
#'
#' @examples
raster_dvhsi_lyr <- function(d, v) {
  #(d>0.5 & d<=5.2) & (v>0 & v<=4.0)
  terra::lapp(c(d, v), raster_dvhsi_lyr)
}

#' Raster Habitat Suitability Index (HSI) function: Fall-Run Chinook Spawning
#'
#' @description
#'
#' @param d a `terra` raster layer for depth in feet
#' @param v a `terra` raster layer for velocity in feet per second
#'
#' @return A `terra` raster layer of habitat suitability index values between 0 and 1.
#' @md
#' @export
#'
#' @examples
raster_dvhsi_spawning <- function(d, v, run="max") {
  terra::lapp(c(d, v), function(d, v) vector_dvhsi_spawning(d, v, run=run))
}

#' Calculate HSI based on HEC-RAS 2D rasters
#'
#' @param rasters A named list or `SpatRasterDataset` of `depth` and `velocity` raster stacks (in `terra` format) containing model results for each flow, as output by `raster_prep_grid()`
#' @param group_polygons An `sf` polygon feature collection with one polygon for each analysis reach. These may have been delineated previously by splitting the hydraulic model domain into subsections.
#' @param mask_polygons Optional: An `sf` polygon feature or feature collection containing a mask to be eliminated from the study area; for example, a baseflow channel. Only works if group polygons are also supplied.
#' @param mask_raster Optional: A `terra` raster of 1 or 0 to be multipled for the purpose of masking out a baseflow channel.
#' @param hsi_func A `terra::lapp`-compatible function that takes depth (ft) and velocity (ft/s) raster layers and returns habitat suitability index values between 0 and 1. Function arguments should be `d` and `v` in US units. Defaults to the `raster_dvhsi_hqt` function in this package.
#' @param .group_var The unquoted name of the attribute used to identify groups in the `group_polygons` layer. Defaults to `comid`
#' @param parallel Toggles parallel processing. Currently under development, keep `FALSE` for now.
#'
#' @return A `tbl_df` data frame with one row per group (`comid`) per flow (`flow_cfs`), containing the total inundated area (`area_tot`) in square feet, the suitable habitat area (`area_wua`) in square feet, and the ratio of suitable to total area (`area_pct`).
#' @md
#' @export
#'
#' @examples
raster_summarize_hsi <- function(rasters,
                                 group_polygons = NULL,
                                 mask_polygons = NULL,
                                 mask_raster = NULL,
                                 hsi_func = raster_dvhsi_hqt,
                                 .group_var = comid,
                                 parallel = FALSE){

  if(!is.null(group_polygons)){
    grp <- group_polygons |>
      group_by({{.group_var}}) |>
      summarize() |>
      st_union(by_feature=T)

    if(!is.null(mask_polygons)) {
      mask <- mask_polygons |>
        st_transform(st_crs(group_polygons))
      grp <- grp |>
        st_difference(mask)
    }

    grp_vect <- terra::vect(grp)
  }

  #cell_area_ft2 <- prod(terra::res(rasters$depth)) * (terra::linearUnits(rasters$depth) / 0.3048)^2 # International Ft (Oregon, etc.)
  cell_area_ft2 <- prod(terra::res(rasters$depth)) * (terra::linearUnits(rasters$depth) * 3937 / 1200)^2 # US Survey Feet (California)

  if (parallel) {

    dep <- terra::wrap(rasters$depth)
    vel <- terra::wrap(rasters$velocity)

    if(!is.null(mask_raster)) {
      msk <- terra::wrap(mask_raster)
    } else {
      msk <- NULL
    }

    if (future::availableCores()>2) {
      future::plan(future::multisession, workers = future::availableCores() - 1)
    } else {
      future::plan(future::sequential)
    }

  } else {

    dep <- rasters$depth
    vel <- rasters$velocity

    if(!is.null(mask_raster)) {
      msk <- mask_raster
    } else {
      msk <- NULL
    }
  }

  # TODO: Fix issue where parallelized version is not able to pass through the non-standard symbol for comid name. Need to quo/unquo to pass. Can't properly use parallel processing until this is fixed.
  # .group_var_str = rlang::as_string(.group_var)
  # message(.group_var_str)

  run_for_flow <- function(flow, dep, vel, msk = NULL) {

    if (parallel) {
      d <- terra::unwrap(dep)[[as.character(flow)]]
      v <- terra::unwrap(vel)[[as.character(flow)]]
    } else {
      d <- dep[[as.character(flow)]]
      v <- vel[[as.character(flow)]]
    }

    message(paste(flow, " - calculate inundation extent"))
    ext <- (d > 0)

    message(paste(flow, " - apply habitat suitability function"))
    hsi <- hsi_func(d, v)
    #hsi <- terra::lapp(c(d, v), hsi_func)

    if(!is.null(msk)) {
      message(paste(flow, " - apply mask"))
      hsi <- hsi * msk
    }

    if(!is.null(group_polygons)){

      message(paste(flow, " - summarize inundation extent by group"))
      area_tot <- terra::zonal(ext, grp_vect, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      message(paste(flow, " - summarize suitable habitat area by group"))
      area_wua <- terra::zonal(hsi, grp_vect, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      return(tibble({{.group_var}} := grp |> pull({{.group_var}}),
                    area_tot,
                    area_wua,
                    area_pct = area_wua / area_tot))

    } else {

      message(paste(flow, " - summarize inundation extent"))
      area_tot <- terra::global(ext, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      message(paste(flow, " - summarize suitable habitat area"))
      area_wua <- terra::global(hsi, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      return(list(area_tot = area_tot,
                  area_wua = area_wua,
                  area_pct = area_wua / area_tot))

    }
  }

  if (parallel) {

    result <- names(rasters$depth) |>
      setNames(names(rasters$depth)) |>
      future.apply::future_lapply(function(x) run_for_flow(x, dep, vel, msk), future.seed = 47) |>
      bind_rows(.id = "flow_cfs") |>
      mutate(flow_cfs = as.numeric(flow_cfs))

  } else {

    result <- names(rasters$depth) |>
      setNames(names(rasters$depth)) |>
      lapply(function(x) run_for_flow(x, dep, vel, msk)) |>
      bind_rows(.id = "flow_cfs") |>
      mutate(flow_cfs = as.numeric(flow_cfs))
  }

  return(result)

}

#' Calculate First Inundating Flow
#'
#' @param rasters A named list or `SpatRasterDataset` of `depth` and `velocity` raster stacks (in `terra` format) containing model results for each flow, as output by `raster_prep_grid()`. Only depth is strictly necessary as velocity is not used in this function. Layer names should be the flow in cfs.
#'
#' @return A `terra` `SpatRaster` in which each grid cell is identified by the first modeled flow where depth exceeded 0 ft.
#' @md
#' @export
#'
#' @examples
raster_first_inundating_flow <- function(rasters) {

  rst <- rasters[[1]]

  get_inundation_bnd <- function(q) {
    message(q, " - get inundation extent")
    (terra::ifel(rst[[q]] > 0, 1, NA) * as.numeric(q))
  }

  inundation_bnds <-
    names(rst) |>
    lapply(get_inundation_bnd) |>
    terra::rast()

  message("calculate first inundating flow raster")
  first_inundating_flow <- inundation_bnds |>
    terra::app(fun = "min", na.rm=TRUE)
  terra::set.names(first_inundating_flow, "flow_cfs")

  return(first_inundating_flow)

}
