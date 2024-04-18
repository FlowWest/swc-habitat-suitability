################################################################################
# VECTOR MESH FUNCTIONS
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
################################################################################

# Prep Flow Raster Stacks ------------------------------------------------------

# Input: a named list of raster filenames, named by flow
prep_rast <- function(filenames){

  dep_r <- terra::rast(filenames$depth)
  terra::set.names(dep_r, filenames$flow_cfs)

  vel_r <- terra::rast(filenames$velocity)
  terra::set.names(vel_r, filenames$flow_cfs)

  #return(list("depth" = terra::wrap(dep_r),
  #            "velocity" = terra::wrap(vel_r)))

  return(list("depth" = dep_r, #terra::wrap(dep_r),
              "velocity" = vel_r))#,terra::wrap(vel_r)))

}

# HSI Function -----------------------------------------------------------------
dvhsi <- function(d, v) (d>0.5 & d<=5.2) & (v>0 & v<=4.0)
# ------------------------------------------------------------------------------

raster_hsi <- function(rasters,
                       group_polygons,
                       mask_raster,
                       hsi_func){

  if(!missing(group_polygons)){
    grp <- terra::vect(group_polygons)
  }

  #dep <- terra::unwrap(rasters$depth)
  #vel <- terra::unwrap(rasters$velocity)

  cell_area_ft2 <- prod(terra::res(rasters$depth)) * (terra::linearUnits(rasters$depth)/0.3048)^2

  run_for_flow <- function(flow) {
    d <- rasters$depth[[as.character(flow)]]
    v <- rasters$velocity[[as.character(flow)]]

    message(paste(flow, " - calculate inundation extent"))
    ext <- (d > 0)
    message(paste(flow, " - apply habitat suitability function"))
    hsi <- hsi_func(d, v) * mask_raster

    if(!missing(group_polygons)){
      message(paste(flow, " - summarize inundation extent"))
      area_tot <- terra::zonal(ext, grp, "sum", na.rm=TRUE)[[1]] * cell_area_ft2
      message(paste(flow, " - summarize suitable habitat area"))
      area_wua <- terra::zonal(hsi, grp, "sum", na.rm=TRUE)[[1]] * cell_area_ft2
    } else {
      message(paste(flow, " - summarize inundation extent"))
      area_tot <- terra::global(ext, "sum", na.rm=TRUE)[[1]] * cell_area_ft2
      message(paste(flow, " - summarize suitable habitat area"))
      area_wua <- terra::global(hsi, "sum", na.rm=TRUE)[[1]] * cell_area_ft2
    }

    area_pct <- area_wua / area_tot

    return(list("area_tot"=area_tot,
                "area_wua"=area_wua,
                "area_pct"=area_pct))
  }

  result <- names(rasters$depth) |> lapply(run_for_flow) |> bind_rows()

  return(result)

}





























