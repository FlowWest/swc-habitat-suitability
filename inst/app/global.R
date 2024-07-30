library(tidyverse)
library(sf)
library(shiny)
library(leaflet)
library(leafgl)
library(habistat)

get_data <- function(...) {
  # use this for loading data from R package such that it can be easily garbage-collected
  # get_data(wua_predicted, package = "habistat")
  # gc()
  e <- new.env()
  name <- data(..., envir = e)[1]
  e[[name]]
}

# DATA IMPORT ------------------------------------------------------------------

predictions <- get_data(wua_predicted, package = "habistat") |>
  #filter(habitat == "rearing") |>
  mutate(model_id = if_else((habitat=="rearing" & model_bfc),
                            paste0(model_name, "_ph_bfc_rm"), # post-hoc baseflow channel removal
                            model_name)) |>
  select(comid, flow_cfs, habitat, model_id, wua_per_lf_pred, river_cvpia, watershed_level_3, reach_length_ft) |>
  pivot_wider(names_from = model_id,
              values_from = wua_per_lf_pred,
              names_glue = c("wua_per_lf_pred_{model_id}")) |>
  left_join(get_data(wua_hydraulic_interp, package = "habistat") |>
              filter(!bfc) |> # only showing the actuals with prior bfc removed for now
              select(comid, flow_cfs, wua_per_lf_actual = wua_per_lf),
            by = join_by(comid, flow_cfs)) |>
  left_join(get_data(flowline_attr, package = "habistat") |>
              transmute(comid,
                        chan_width_ft = chan_width_m/0.3048,
                        baseflow_cfs = nf_bfl_dry_cfs),
            by = join_by(comid))

gc() # garbage collect after loading from habistat package

all_flows <- unique(predictions$flow_cfs)

attr <- get_data(flowline_attr, package = "habistat") |>
  filter(comid %in% predictions$comid) |>
  # stream size filter -- make sure this matches model_cleaned.Rmd
  filter(((stream_order >= 4) & (da_area_sq_km > 1)) | ((stream_order >= 3) & (da_area_sq_km >= 50)))

gc() # garbage collect after loading from habistat package

geom <- get_data(flowline_geom, package = "habistat") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  inner_join(attr |> transmute(comid, gnis_name), by = join_by(comid)) |>
  mutate(object_id = paste0("comid_", comid))

gc() # garbage collect after loading from habistat package

#predictions_summary <- predictions |>
#  group_by(watershed_level_3, flow_cfs, reach_length_ft) |>
#  summarize(across(starts_with("wua_per_lf_pred_"), function(x) sum(x * reach_length_ft) / sum(reach_length_ft)))

watersheds <- get_data(cv_watersheds, package = "habistat") |>
  group_by(watershed_level_1, watershed_level_2, watershed_level_3) |> #, range_pisces) |>
  summarize() |>
  ungroup() |>
  st_union(by_feature=T) |>
  st_transform("+proj=longlat +datum=NAD83") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  mutate(watershed_id = paste0("watershed_", row_number()),
         watershed_label = paste0(watershed_level_1, " Basin",
                                  if_else(watershed_level_2 != watershed_level_1,
                                          paste0("<br />", watershed_level_2, " Watershed"),
                                          " (Direct Drainage)"),
                                  if_else(watershed_level_3 != watershed_level_2,
                                          paste0("<br />", watershed_level_3, " Subwatershed"),
                                          "")))

gc()

# PLOT STYLES AND PALETTES -----------------------------------------------------

pal <- function(x, n = 10,
                palette_function = viridisLite::cividis,
                palette_args = list(direction = -1)) {
  #palette <- palette_function(n)
  palette <- do.call(palette_function, c(list(n=n), palette_args))
  palette[cut(x, n)]
}

theme_set(theme_minimal())

ihs <- scales::trans_new("ihs",
                         transform = function(x) asinh(x),
                         inverse = function(y) sinh(y),
                         breaks = function(i) scales::breaks_log(n=5, base=10)(pmax(i,0.01)),
                         #minor_breaks = scales::minor_breaks_n(n = 0),
                         domain=c(0, Inf),
                         format = scales::label_comma())

# before deploying shiny app, install package from github:
# options(timeout=600)
# devtools::install_github("FlowWest/swc-habitat-suitability@training-data-no-bfc-removal")
