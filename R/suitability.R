#' Flowline length helper function: Extract NHD flowlines that correspond to groups and calculate length
#'
#' @param group_polygons An `sf` polygon feature collection with one polygon for each analysis reach. These may have been delineated previously by splitting the hydraulic model domain into subsections.
#' @param .group_var The unquoted name of the attribute used to identify groups in the `group_polygons` layer. Defaults to `comid`
#'
#' @return A `tbl_df` data frame, one row per group (`comid`), containing the calculated length of the associated NHDPlusV2 flowline (joined by `comid`) clipped to the group polygon extent.
#' @md
#'
#' @examples
#' @export
suitability_measure_flowline_lengths <- function(group_polygons, .group_var = comid) {

  habistat::flowline_geom_proj |>
    filter({{.group_var}} %in% (group_polygons |> pull({{.group_var}}))) |>
    st_transform(st_crs(group_polygons)) |>
    st_intersection(group_polygons |>
                      summarize() |>
                      st_union(by_feature = F)) |>
    mutate(length_ft = st_length(geometry) |>
             units::set_units("ft") |> units::drop_units()) |>
    st_drop_geometry()
}

#' Postprocessor for both SRH-2D and HEC-RAS results, applying linear ft normalization by group
#'
#' @param result_tbl A `tbl_df` data frame with one row per group (`comid`) per flow (`flow_cfs`), containing the total inundated area (`area_tot`) in square feet, the suitable habitat area (`area_wua`) in square feet, and the ratio of suitable to total area (`area_pct`). This is the format returned by `raster_summarize_hsi()` and `vector_summarize_hsi()`.
#' @param group_polygons An `sf` polygon feature collection with one polygon for each analysis reach. These may have been delineated previously by splitting the hydraulic model domain into subsections.
#' @param .group_var The unquoted name of the attribute used to identify groups in the `group_polygons` layer. Defaults to `comid`
#'
#' @return A `tbl_df` data frame with one row per group (`comid`) per flow (`flow_cfs`), containing the total inundated area (`area_tot`) in square feet, the suitable habitat area (`area_wua`) in square feet, the ratio of suitable to total area (`area_pct`), the total inundated area per linear ft (`ind_per_lf`), and the suitable habitat area per linear ft (`wua_per_lf`).
#' @md
#' @export
#'
#' @examples
#' @export
suitability_postprocess <- function(result_tbl, group_polygons, .group_var = comid) {

  flowline_lengths <- group_polygons |>
    suitability_measure_flowline_lengths(.group_var = {{.group_var}})

  result_tbl |>
    left_join(flowline_lengths,
              by=join_by({{.group_var}}),
              relationship="many-to-one") |>
    mutate(ind_per_lf = area_tot / length_ft,
           wua_per_lf = area_wua / length_ft) |>
    arrange(comid, flow_cfs)

}
