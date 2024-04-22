library(tidyverse)
library(sf)
library(shiny)
library(leaflet)

# DATA IMPORT ------------------------------------------------------------------

flowline_attributes <- readRDS("../../data/flowline_attributes.Rds")

flowlines_gcs <- readRDS("../../data/flowline_geometries_leaflet.Rds") |>
  mutate(object_id = paste0("comid_", comid))

predictions_table <- readRDS("../../data/predictions_table.Rds") |>
  left_join(flowline_attributes |> transmute(comid, gnis_name, chan_width_ft = chan_width_m/0.3048),
            by=join_by(comid), relationship="many-to-one")

all_flows <- unique(predictions_table$flow_cfs)

# PLOT STYLES AND PALETTES -----------------------------------------------------

pal <- function(x, n = 10,
                palette_function = viridisLite::cividis,
                palette_args = list(direction = -1)) {
  #palette <- palette_function(n)
  palette <- do.call(palette_function, c(list(n=n), palette_args))
  palette[cut(x, n)]
}

theme_set(theme_minimal())

ihs <- trans_new("ihs",
                 transform = function(x) asinh(x),
                 inverse = function(y) sinh(y),
                 breaks = function(i) scales::breaks_log(n=5, base=10)(pmax(i,0.01)),
                 #minor_breaks = scales::minor_breaks_n(n = 0),
                 domain=c(0, Inf),
                 format = scales::label_comma())
