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

# Import SRH2D tabular results -------------------------------------------------

vector_import_srh2d <- function(csv_tbl,
                         units = "ft",
                         units_xy = units,
                         units_z = units,
                         .filename_var = filename,
                         has_vid_col = FALSE) {

  conv_fact_xy <- case_when(units_xy == "ft" ~ 1,
                            units_xy == "m" ~ 3.28084)

  conv_fact_z <- case_when(units_z == "ft" ~ 1,
                           units_z == "m" ~ 3.28084)


  cols <- c("x"="n", "y"="n", "z"="n",
           "wse"="n", "wdepth"="n",
           "vel_x"="n", "vel_y"="n", "vel_mag"="n",
           "froude"="n", "stress"="n")

  if (has_vid_col) {
    cols <- c("vid"="n", cols)
  }

  csv_tbl |>
    mutate(tbl = map({{.filename_var}}, function(x) {
      message("import SRH-2D results from ", x)
      csv_in <- read_csv(x, col_names=names(cols), col_types=cols, skip=1)
      if (has_vid_col) {
        return(csv_in)
      } else {
        return(csv_in |> mutate(vid = row_number()))
      }
    })) |>
    unnest(tbl) |>
    mutate(across(everything(), function(n) if_else(n==-999, NA, n))) |>
    mutate(across(c("x", "y", "vel_x", "vel_y", "vel_mag"), function(n) n * conv_fact_xy)) |>
    mutate(across(c("z", "wse", "wdepth"), function(n) n * conv_fact_z))
}

# HSI Function -----------------------------------------------------------------
vector_dvhsi_hqt <- function(d, v) {
  if_else(d>1.0 & d<=3.28 & v>0 & v<=1.5, 1, 0)
}
vector_dvhsi_lyr <- function(d, v) {
  if_else(d>0.5 & d<=5.2 & v>0 & v<=4.0, 1, 0)
}

# Calculate HSI for each cell --------------------------------------------------

vector_calculate_hsi <- function(data, hsi_func = vector_dvhsi_hqt) {

  data |>
    select(flow_cfs, vid, wdepth, vel_mag) |>
    mutate(inundated = wdepth > 0,
           hsi = hsi_func(wdepth, vel_mag))

}

# Prep Mesh --------------------------------------------------------------------
# Input mesh polygons (e.g. voronoi)
# Output with cell areas calculated and baseflow channel attributes

vector_prep_mesh <- function(mesh,
                      group_polygons = NULL,
                      .group_var = comid,
                      .id_var = vid) {

  groups <- group_polygons |>
    group_by({{.group_var}}) |>
    summarize() |>
    st_union(by_feature = TRUE)

  result <- mesh |>
    mutate(area_ft2 = st_area(geometry) |> units::set_units("ft2") |> units::drop_units()) |>
    st_centroid() |>
    st_intersection(groups) |> # TODO Check performance difference between st_join and st_intersection
    st_drop_geometry() |>
    transmute({{.id_var}},
              {{.group_var}},
              area_ft2)

  return(result)

}

vector_summarize_hsi <- function(mesh_tbl, hsi_tbl, .group_var = comid, .id_var=vid) {
  mesh_tbl |>
    inner_join(hsi_tbl, by=join_by({{.id_var}}), relationship="one-to-many") |>
    group_by({{.group_var}}, flow_cfs) |>
    summarize(area_tot = sum(if_else(wdepth > 0, area_ft2, 0)),
              area_wua = sum(area_ft2 * hsi)) |>
    mutate(area_pct = area_wua / area_tot)
}

################################################################################
# RASTER FUNCTIONS
# Use for native HEC-RAS 2D model outputs
################################################################################

# Prep Flow Raster Stacks ------------------------------------------------------

# Input: a data frame of raster filenames with columns flow_cfs, depth, velocity
raster_prep_grid <- function(filenames){

  dep_r <- terra::rast(filenames$depth)
  terra::set.names(dep_r, filenames$flow_cfs)

  vel_r <- terra::rast(filenames$velocity)
  terra::set.names(vel_r, filenames$flow_cfs)

  return(list("depth" = dep_r,
              "velocity" = vel_r))

}

# Depth/velocity HSI functions to be selected from -----------------------------

raster_dvhsi_hqt <- function(d, v) (d>1.0 & d<=3.28) & (v>0 & v<=1.5)
raster_dvhsi_lyr <- function(d, v) (d>0.5 & d<=5.2) & (v>0 & v<=4.0)

# Main HSI Calculation Function ------------------------------------------------
# input should be the output from prep_rast()
raster_summarize_hsi <- function(rasters,
                       group_polygons = NULL,
                       hsi_func = raster_dvhsi_hqt,
                       .group_var = comid,
                       parallel = FALSE){

  if(!is.null(group_polygons)){
    grp <- group_polygons |>
      group_by({{.group_var}}) |>
      summarize() |>
      st_union(by_feature=T)

    grp_vect <- terra::vect(grp)
  }

  gggarea_ft2 <- prod(terra::res(rasters$depth)) * (terra::linearUnits(rasters$depth)/0.3048)^2

  if (parallel) {

    dep <- terra::wrap(rasters$depth)
    vel <- terra::wrap(rasters$velocity)

    if (future::availableCores()>2) {
        future::plan(future::multisession, workers = future::availableCores() - 1)
    } else {
        future::plan(future::sequential)
    }

  } else {

    dep <- rasters$depth
    vel <- rasters$velocity

  }

  # TODO: Fix issue where parallelized version is not able to pass through the non-standard symbol for comid name. Need to quo/unquo to pass. Can't properly use parallel processing until this is fixed.
  # .group_var_str = rlang::as_string(.group_var)
  # message(.group_var_str)

  run_for_flow <- function(flow, dep, vel) {

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
    #hsi <- hsi_func(d, v)
    hsi <- terra::lapp(c(d, v), hsi_func)

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
      future.apply::future_lapply(function(x) run_for_flow(x, dep, vel), future.seed = 47) |>
      bind_rows(.id = "flow_cfs") |>
      mutate(flow_cfs = as.numeric(flow_cfs))

  } else {

    result <- names(rasters$depth) |>
      setNames(names(rasters$depth)) |>
      lapply(function(x) run_for_flow(x, dep, vel)) |>
      bind_rows(.id = "flow_cfs") |>
      mutate(flow_cfs = as.numeric(flow_cfs))
  }

  return(result)

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


################################################################################
# ADDITIONAL FUNCTIONS used for both raster and vector
################################################################################

# Prep Flowlines Helper Function -----------------------------------------------
# Extract NHD flowlines that correspond to polys and calculate length

prep_flowlines <- function(group_polygons, .group_var = comid) {

  flowlines <- readRDS(here::here("data-raw", "results", "flowline_geometries_proj.Rds"))
  # TODO: Update this reference to flowline data to use the packaged flowline data

  flowlines |>
    filter({{.group_var}} %in% (group_polygons |> pull({{.group_var}}))) |>
    st_transform(st_crs(group_polygons)) |>
    st_intersection(group_polygons |>
                      summarize() |>
                      st_union(by_feature = F)) |>
    mutate(length_ft = st_length(geometry) |>
             units::set_units("ft") |> units::drop_units())
}

# Postprocess raster, adding reach lengths -------------------------------------

postprocess <- function(result_tbl, group_polygons, .group_var = comid) {

  flowlines_prepped <- group_polygons |> prep_flowlines(.group_var = {{.group_var}})

  result_tbl |>
    left_join(st_drop_geometry(flowlines_prepped),
              by=join_by({{.group_var}}),
              relationship="many-to-one") |>
    mutate(ind_per_lf = area_tot / length_ft,
           wua_per_lf = area_wua / length_ft) |>
    arrange(comid, flow_cfs)

}
