library(tidyverse)
library(sf)
library(terra)

library(future) # parallel processing backend
library(future.apply) # parallelized versions of lapply, sapply, ...
library(furrr) # parallelized versions of purrr functions: map, map2, pmap,...

################################################################################
# VECTOR MESH FUNCTIONS
# Use for native SRH-2D model outputs
################################################################################

# Prep Mesh --------------------------------------------------------------------
# Input mesh polygons (e.g. voronoi)
# Output with cell areas calculated and baseflow channel attributes

prep_mesh <- function(mesh,
                      baseflow_polygon,
                      id_var = vid,
                      group_polygons = NULL,
                      group_var = NULL) {

  bfc <- baseflow_polygon |>
    transmute(is_bfc = TRUE)

  if(st_crs(bfc) != st_crs(mesh)) {
    bfc <- bfc |> st_transform(st_crs(mesh))
  }

  if(!missing(group_polygons)) {
    groups <- group_polygons |>
      select({{group_var}})
  }

  result <- mesh |>
    st_union(bfc) |>
    st_intersection(groups) |>
    transmute({{id_var}},
              {{group_var}},
              bfc_mask = if_else(is_bfc, 0, 1),
              area_ft2 = st_area(geometry) |> units::set_units("ft2") |> units::drop_units())

  return(result)

}

# HSI Function -----------------------------------------------------------------
dvhsi <- function(d, v) {
  if_else(d>0.5 & d<=5.2 & v>0 & v<=4.0, 1, 0)
}

# Process Mesh -----------------------------------------------------------------
# Input table of depth and velocity by mesh cell ID at a range of flows
# Columns: cell, flow, depth, velocity
# Input a vectorized function that calculates HSI via depth and velocity

process_mesh <- function(mesh,
                         flow_velocity_table,
                         hsi_func = dvhsi,
                         depth_var = depth,
                         velocity_var = velocity) {
  result <- mesh |>
    inner_join(flow_velocity_table, by=join_by(cell), relationship="one-to-many") |>
    mutate(hsi = hsi_func({{depth_var}}, {{velocity_var}}),
           wua = area_ft2 * hsi * bfc_mask)
  return(result)
}

summarize_mesh <- function(mesh, group_var = NULL) {
  result <- mesh |>
    group_by({{group_var}}) |>
    summarize(area_tot = sum(area_ft2),
              area_wua = sum(wua_ft2)) |>
    mutate(area_pct = area_wua / area_tot)
  return(result)
}


################################################################################
# RASTER FUNCTIONS
# Use for native HEC-RAS 2D model outputs
################################################################################

# Prep Flow Raster Stacks ------------------------------------------------------

# Input: a data frame of raster filenames with columns flow_cfs, depth, velocity
prep_rast <- function(filenames){

  dep_r <- terra::rast(filenames$depth)
  terra::set.names(dep_r, filenames$flow_cfs)

  vel_r <- terra::rast(filenames$velocity)
  terra::set.names(vel_r, filenames$flow_cfs)

  return(list("depth" = dep_r,
              "velocity" = vel_r))

}

# Depth/velocity HSI functions to be selected from -----------------------------

dvhsi_hqt <- function(d, v) (d>1.0 & d<=3.28) & (v>0 & v<=1.5)
dvhsi_lyr <- function(d, v) (d>0.5 & d<=5.2) & (v>0 & v<=4.0)

# Main HSI Calculation Function ------------------------------------------------
# input should be the output from prep_rast()
raster_hsi <- function(rasters,
                       group_polygons = NULL,
                       #mask_raster = NULL,
                       hsi_func = dvhsi_hqt,
                       .id_var = comid,
                       parallel = FALSE){

  if(!is.null(group_polygons)){
    grp <- terra::vect(group_polygons)
  }

  cell_area_ft2 <- prod(terra::res(rasters$depth)) * (terra::linearUnits(rasters$depth)/0.3048)^2

  if (parallel) {

    dep <- terra::wrap(rasters$depth)
    vel <- terra::wrap(rasters$velocity)

    if (future::availableCores()>1) {
        future::plan(future::multisession, workers = future::availableCores())
    } else {
        future::plan(future::sequential)
    }

  }

  run_for_flow <- function(flow) {

    if (parallel) {
      d <- terra::unwrap(dep)[[as.character(flow)]]
      v <- terra::unwrap(vel)[[as.character(flow)]]
    } else {
      d <- rasters$depth[[as.character(flow)]]
      v <- rasters$velocity[[as.character(flow)]]
    }

    message(paste(flow, " - calculate inundation extent"))
    ext <- (d > 0)

    message(paste(flow, " - apply habitat suitability function"))
    hsi <- hsi_func(d, v) # * mask_raster

    if(!is.null(group_polygons)){

      message(paste(flow, " - summarize inundation extent by group"))
      area_tot <- terra::zonal(ext, grp, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      message(paste(flow, " - summarize suitable habitat area by group"))
      area_wua <- terra::zonal(hsi, grp, "sum", na.rm=TRUE)[[1]] * cell_area_ft2

      return(tibble({{.id_var}} := group_polygons |> pull({{.id_var}}),
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
      future.apply::future_lapply(run_for_flow, future.seed = 47) |>
      bind_rows(.id = "flow_cfs")

  } else {

    result <- names(rasters$depth) |>
      setNames(names(rasters$depth)) |>
      lapply(run_for_flow) |>
      bind_rows(.id = "flow_cfs")
  }

  return(result)

}

# Prep Flowlines Helper Function -----------------------------------------------
# Extract NHD flowlines that correspond to polys and calculate length

prep_flowlines <- function(group_polygons, .id_var = comid) {

  flowlines <- readRDS(here::here("data-raw", "results", "flowline_geometries_proj.Rds"))
  # TODO: Update this reference to flowline data to use the packaged flowline data

  flowlines |>
    filter(comid %in% (group_polygons |> pull({{.id_var}}))) |>
    st_transform(st_crs(group_polygons)) |>
    st_intersection(group_polygons |>
                      summarize() |>
                      st_union(by_feature = F)) |>
    mutate(length_ft = st_length(geometry) |>
             units::set_units("ft") |> units::drop_units())
}

# Postprocess raster, adding reach lengths -------------------------------------

post_rast <- function(result_tbl, group_polygons, .id_var = comid) {

  flowlines_prepped <- group_polygons |> prep_flowlines(.id_var = {{.id_var}})

  result_tbl |>
    left_join(st_drop_geometry(flowlines_prepped),
            by=join_by({{.id_var}}),
            relationship="many-to-one") |>
    mutate(ind_per_lf = area_tot / length_ft,
           wua_per_lf = area_wua / length_ft) |>
    arrange(comid, flow_cfs)

}

# First Inundating Flow Raster Function ----------------------------------------

# input should be the output from prep_rast()
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
