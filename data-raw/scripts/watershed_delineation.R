library(tidyverse)
library(sf)
googledrive::drive_auth()

st_identity <- function(x, y, ...) {
  st_intersection(x, y, ...) |>
  bind_rows(st_difference(x, st_union(y,...),...) |> st_collection_extract("POLYGON"))
}

st_coalesce <- function(x, y, ...) {
  bind_rows(x, st_difference(y, st_union(x, ...),...))
}

#st_identity <- function(x, y, model="open", dimensions="polygon",...) {
#  st_intersection(x, y, model=model, dimensions=dimensions,...) |>
#    bind_rows(st_difference(x, st_union(y, model=model, dimensions=dimensions,...), model=model, dimensions=dimensions,...) |> st_collection_extract#("POLYGON"))
#}
#
#st_coalesce <- function(x, y, model="open", dimensions="polygon",...) {
#  bind_rows(x, st_difference(y, st_union(x, model=model, dimensions=dimensions, ...), model=model, dimensions=dimensions,...))
#}
#

catchments <- habistat::drive_file_by_id("1ncwKAUNoJUNkPLEy6NzrUKCYqVG681p-") |>
  read_sf() |>
  mutate(HUC_4 = substr(HUC_8, 1, 4), HUC_6 = substr(HUC_8, 1, 6))
catchments

huc2poly <- function(filter_expr) {
  fc <- catchments |>
    filter({{filter_expr}}) |>
    summarize() |>
    st_union(by_feature = FALSE)
  return(fc[[1]])
}

watersheds_nested <-
  tribble(~level, ~parent, ~river, ~geometry,
          1, NA,                  "Sacramento River",        huc2poly(HUC_4 == 1802 & HUC_8 != 18020001),
          #2, "Sacramento River",  "North Delta",             huc2poly(HUC_10 %in% c(1802016305, 1802016306, 1802016307)),
          2, "Sacramento River",  "Upper Sacramento River",  huc2poly(HUC_8 %in% c(18020005)),
          2, "Sacramento River",  "McCloud River",           huc2poly(HUC_8 %in% c(18020004)),
          2, "Sacramento River",  "Pit River",               huc2poly(HUC_8 %in% c(18020003, 18020002)),
          2, "Sacramento River",  "Cow Creek",               huc2poly(HUC_8 %in% c(18020151)),
          2, "Sacramento River",  "Clear Creek",             huc2poly(HUC_10 %in% c(1802015401)),
          2, "Sacramento River",  "Bear Creek",              huc2poly(HUC_10 %in% c(1802015404)),
          2, "Sacramento River",  "Cottonwood Creek",        huc2poly(HUC_8 %in% c(18020152)),
          2, "Sacramento River",  "Battle Creek",            huc2poly(HUC_8 %in% c(18020153)),
          2, "Sacramento River",  "Paynes Creek",            huc2poly(HUC_10 %in% c(1802015501)),
          2, "Sacramento River",  "Antelope Creek",          huc2poly(HUC_10 %in% c(1802015601)),
          2, "Sacramento River",  "Elder Creek",             huc2poly(HUC_10 %in% c(1802015602)),
          2, "Sacramento River",  "Mill Creek",              huc2poly(HUC_10 %in% c(1802015603)),
          2, "Sacramento River",  "Thomes Creek",            huc2poly(HUC_10 %in% c(1802015604, 1802015605)),
          2, "Sacramento River",  "Deer Creek",              huc2poly(HUC_10 %in% c(1802015702)),
          2, "Sacramento River",  "Big Chico Creek",         huc2poly(HUC_10 %in% c(1802015705)),
          2, "Sacramento River",  "Stony Creek",             huc2poly(HUC_8 %in% c(18020115)),
          2, "Sacramento River",  "Butte Creek",             huc2poly(HUC_10 %in% c(1802015801, 1802015802, 1802015803, 1802015804)),
          2, "Sacramento River",  "Feather River",           huc2poly(HUC_8 %in% c(18020159, 18020121, 18020122, 18020123, 18020124, 18020125, 18020126)),
          3, "Feather River",     "Yuba River",              huc2poly(HUC_8 %in% c(18020125)),
          3, "Feather River",     "Bear River",              huc2poly(HUC_8 %in% c(18020126)),
          2, "Sacramento River",  "American River",          huc2poly(HUC_8 %in% c(18020111, 18020128, 18020129)),
          1, NA,                  "San Joaquin River",       huc2poly(HUC_4 == 1804),
          #2, "San Joaquin River", "South Delta",             huc2poly(HUC_10 %in% c(1804000307, 1804000308, 1804000309)),
          2, "San Joaquin River", "Mokelumne River",         huc2poly(HUC_8 %in% c(18040012, 18040013)),
          3, "Mokelumne River",   "Cosumnes River",          huc2poly(HUC_8 %in% c(18040013)),
          2, "San Joaquin River", "Calaveras River",         huc2poly(HUC_8 %in% c(18040011)),
          2, "San Joaquin River", "Stanislaus River",        huc2poly(HUC_8 %in% c(18040010)),
          2, "San Joaquin River", "Tuolumne River",          huc2poly(HUC_8 %in% c(18040009)),
          2, "San Joaquin River", "Merced River",            huc2poly(HUC_8 %in% c(18040008)),
          2, "San Joaquin River", "Upper San Joaquin River", huc2poly(HUC_10 %in% c(1804000103) | HUC_8 %in% c(18040006)),
        ) |>
  st_as_sf(crs=st_crs(catchments)) |>
  st_transform(habistat::const_proj_crs()) |>
  mutate(watershed_level_1 = if_else(level <= 1, river, NA),
         watershed_level_2 = if_else(level <= 2, river, NA),
         watershed_level_3 = if_else(level <= 3, river, NA)) |>
  fill(starts_with("watershed_level_"))
#|>
  #mutate(wid = row_number())

clip_out_children <- function(g, n) {

  cwdf <- watersheds_nested |>
    filter(parent == n)

  if (nrow(cwdf) > 0) {
    child_watersheds <- cwdf |>
      summarize() |>
      st_union()

    res <- st_difference(g, child_watersheds[[1]])

    } else {

    res <- g

    }

  return(res)
}

watersheds_clip <- watersheds_nested |>
  mutate(geometry = st_sfc(map2(geometry, river, clip_out_children), crs = habistat::const_proj_crs()))

watersheds_clip |> head(1) |> plot()

# IMPORT CVPIA MAINSTEMS -------------------------------------------------------

cv_mainstems <-
  read_sf(here::here("data-raw", "source", "rearing_spatial_data", "nhdplusv2_comid_habitat_xw.shp.zip")) |>
  janitor::clean_names() |>
  st_zm() |>
  st_transform(habistat::const_proj_crs()) |>
  select(river_group = watershed, river_cvpia = river, comid, habitat) |>
  filter(!is.na(river_group) & !is.na(river_cvpia))

# ADD UC DAVIS PISCES ----------------------------------------------------------

pisces_ranges <- tribble(~sp_id, ~species, ~dataset, ~filename,
                         "WR",  "Central Valley Winter Run Chinook Salmon",    "Historical", "Oncorhynchus_tshawytscha_SOT05_historical_expert_16.shp.zip",
                         "WR",  "Central Valley Winter Run Chinook Salmon",    "Extant",     "Oncorhynchus_tshawytscha_SOT05_extant_1.shp.zip",
                         "WR",  "Central Valley Winter Run Chinook Salmon",    "Observed",   "Oncorhynchus_tshawytscha_SOT05_observed_2.shp.zip",
                         "SR",  "Central Valley Spring Run Chinook Salmon",    "Historical", "Oncorhynchus_tshawytscha_SOT06_historical_expert_16.shp.zip",
                         "SR",  "Central Valley Spring Run Chinook Salmon",    "Extant",     "Oncorhynchus_tshawytscha_SOT06_extant_1.shp.zip",
                         "SR",  "Central Valley Spring Run Chinook Salmon",    "Observed",   "Oncorhynchus_tshawytscha_SOT06_observed_2.shp.zip",
                         "LFR", "Central Valley Late Fall Run Chinook Salmon", "Historical", "Oncorhynchus_tshawytscha_SOT07_historical_expert_16.shp.zip",
                         "LFR", "Central Valley Late Fall Run Chinook Salmon", "Extant",     "Oncorhynchus_tshawytscha_SOT07_extant_1.shp.zip",
                         "FR",  "Central Valley Fall Run Chinook Salmon",      "Historical", "Oncorhynchus_tshawytscha_SOT08_historical_expert_16.shp.zip",
                         "FR",  "Central Valley Fall Run Chinook Salmon",      "Extant",     "Oncorhynchus_tshawytscha_SOT08_extant_1.shp.zip") |>
  mutate(filepath = file.path("/vsizip", here::here("data-raw/source/ucd_pisces_ranges", filename))) |>
  mutate(result = map(filepath, read_sf)) |>
  unnest(result) |>
  janitor::clean_names() |>
  st_sf() |>
  st_transform(habistat::const_proj_crs())

range_extant <-
  pisces_ranges |>
  filter(dataset=="Extant") |>
  summarize() |>
  st_union() |>
  st_as_sf() |>
  mutate(range = "extant")

range_historical <-
  pisces_ranges |>
  filter(dataset=="Historical") |>
  summarize() |>
  st_union() |>
  st_difference(st_union(range_extant)) |>
  st_as_sf() |>
  mutate(range = "historical")

ranges <- bind_rows(range_historical, range_extant)

# IMPORT BYPASS AND DELTA ------------------------------------------------------

cv_bypasses <-
  read_sf(here::here("data-raw", "source", "rearing_spatial_data", "yolo_sutter_bypass_extents.shp.zip")) |>
  st_transform(habistat::const_proj_crs()) |>
  mutate(bypass = if_else(str_detect(area_name, "Butte Sink"), "Butte Sink", paste(bypass, "Bypass")),
         bypass_parent = if_else(str_detect(area_name, "Butte Sink"), "Butte Creek", "Sacramento River"),
         bypass_level_1 = "Sacramento River",
         bypass_level_2 = if_else(str_detect(area_name, "Butte Sink"), "Butte Creek", bypass),
         bypass_level_3 = bypass) |>
  group_by(bypass, bypass_parent, bypass_level_1, bypass_level_2, bypass_level_3) |>
  summarize() |>
  st_union(by_feature = T)

cv_bypasses |> saveRDS(here::here("data-raw", "results", "cv_bypasses.Rds"))

cv_delta <-
  read_sf(here::here("data-raw", "source", "rearing_spatial_data", "Legal_Delta_Boundary.shp.zip")) |>
  st_transform(habistat::const_proj_crs()) |>
  select(geometry) |>
  st_intersection(
    tribble(~delta_parent, ~delta, ~geometry,
            "Sacramento River", "North Delta", huc2poly(HUC_4 == 1802),
            "San Joaquin River", "South Delta", huc2poly(HUC_4 == 1804)) |>
      st_as_sf(crs=st_crs(catchments)) |>
      st_transform(habistat::const_proj_crs())) |>
  mutate(delta_level_1 = delta_parent,
         delta_level_2 = delta,
         delta_level_3 = delta)

cv_delta |> saveRDS(here::here("data-raw", "results", "cv_delta.Rds"))

cv_bypasses_delta <-
  st_coalesce(cv_bypasses, cv_delta) |>
  transmute(bypass_delta_name = coalesce(bypass, delta),
            bypass_delta_parent = coalesce(bypass_parent, delta_parent),
            bypass_delta_level_1 = coalesce(bypass_level_1, delta_level_1),
            bypass_delta_level_2 = coalesce(bypass_level_2, delta_level_2),
            bypass_delta_level_3 = coalesce(bypass_level_3, delta_level_3))

cv_bypasses_delta |> plot()

# COMBINE ----------------------------------------------------------------------

polygon_filter <- function(g) {
  if ("GEOMETRYCOLLECTION" %in% class(g)) {
    (g |>
       st_collection_extract("POLYGON") |>
       st_combine())[[1]]
  } else {
    g
  }
}

cv_watersheds <-
  watersheds_clip |>
  # intersect the pisces ranges
  st_identity(ranges) |>
  mutate(range = coalesce(range, "none")) |>
  rename(range_pisces = range,
         watershed_cvpia = river,
         watershed_cvpia_parent = parent) #|>
  # merge in delta and bypasses, replacing overlapping areas
  st_identity(cv_bypasses_delta) |>
  mutate(watershed_cvpia = coalesce(bypass_delta_name, watershed_cvpia),
         watershed_cvpia_parent = coalesce(bypass_delta_parent, watershed_cvpia_parent),
         watershed_level_1 = coalesce(bypass_delta_level_1, watershed_level_1),
         watershed_level_2 = coalesce(bypass_delta_level_2, watershed_level_2),
         watershed_level_3 = coalesce(bypass_delta_level_3, watershed_level_3)) |>
  select(-starts_with("bypass_delta_")) |>
  # eliminate dangling linestrings
  mutate(geometry = st_sfc(lapply(geometry, polygon_filter), crs=habistat::const_proj_crs())) |>
  st_collection_extract(type = "POLYGON")

#cv_watersheds |> plot()

cv_watersheds |> saveRDS(here::here("data-raw", "results", "cv_watersheds.Rds"))
cv_watersheds |> usethis::use_data(overwrite = T)

#mainstems |> st_join(cv_watersheds, largest=T, left=T)

cv_mainstems |> saveRDS(here::here("data-raw", "results", "cv_mainstems.Rds"))
cv_mainstems |> usethis::use_data(overwrite = T)

ggplot() +
  geom_sf(data = cv_watersheds, aes(fill = watershed_cvpia)) +
  geom_sf(data = cv_mainstems, aes(color = river_cvpia)) +
  scale_fill_hue(aesthetics = c("fill", "color"))

ggplot() +
  geom_sf(data = cv_watersheds |> filter(range_pisces == "extant" & watershed_level_1 == "Sacramento River"), aes(fill = watershed_cvpia))

