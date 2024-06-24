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

#predictions_table <- habistat::wua_predicted |>
#  inner_join(habistat::flowline_attr |> transmute(comid, gnis_name, chan_width_ft = chan_width_m/0.3048),
#            by=join_by(comid), relationship="many-to-one") #|>
#  #filter(model_type == "RF")
#  # TODO Filter for just RF

predictions <- get_data(wua_predicted, package = "habistat") |>
  mutate(model_bfc = if_else(model_bfc, "p", "a")) |> # baseflow channel removed pre (a) or post (p) prediction
  mutate(model_id = paste0(model_name, "_", model_bfc)) |>
  select(comid, flow_cfs, model_id, wua_per_lf_pred, river_cvpia, watershed_level_3) |>
  pivot_wider(names_from = model_id,
              values_from = wua_per_lf_pred,
              #values_fn = first,
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

geom <- get_data(flowline_geom, package = "habistat") |>
  st_set_crs("+proj=longlat +datum=WGS84") |> # for display purposes only
  mutate(object_id = paste0("comid_", comid)) |>
  inner_join(get_data(flowline_attr, package = "habistat") |>
               transmute(comid, gnis_name), # chan_width_ft = chan_width_m/0.3048),
             by = join_by(comid))

gc() # garbage collect after loading from habistat package

# attr <- get_data(flowline_attr, package = "habistat") |>
#   filter(comid %in% predictions$comid)
#
# gc() # garbage collect after loading from habistat package

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
