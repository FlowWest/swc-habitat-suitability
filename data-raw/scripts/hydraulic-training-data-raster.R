library(tidyverse)
library(sf)
library(terra)

#source(here::here("data-raw", "scripts", "data-functions.R"))
#source(here::here("data-raw", "scripts", "suitability-functions.R"))
library(habistat)

flowlines <- readRDS(here::here("data-raw", "results", "flowline_geometries_proj.Rds"))

googledrive::drive_auth() # run this before proceeding with rest of script

# TUOLUMNE BASSO-LA GRANGE -----------------------------------------------------

basso_dir <- here::here("data-raw", "temp", "tuolumne-basso")
dir.create(basso_dir, recursive=T)
drive_file_by_id("1nmczhbtWnCfmFaFPTol6rmRjnRh3ovCB") |>
  archive::archive_extract(basso_dir)

basso_filenames <-
  tribble(~flow_cfs, ~depth, ~velocity,
          80, "Depth (02JAN1900 00 00 00).vrt", "Velocity (02JAN1900 00 00 00).vrt",
          110, "Depth (02JAN1900 16 00 00).vrt", "Velocity (02JAN1900 16 00 00).vrt",
          150, "Depth (03JAN1900 08 00 00).vrt", "Velocity (03JAN1900 08 00 00).vrt",
          300, "Depth (03JAN1900 18 00 00).vrt", "Velocity (03JAN1900 18 00 00).vrt",
          500, "Depth (04JAN1900 04 00 00).vrt", "Velocity (04JAN1900 04 00 00).vrt",
          600, "Depth (04JAN1900 14 00 00).vrt", "Velocity (04JAN1900 14 00 00).vrt",
          900, "Depth (05JAN1900 00 00 00).vrt", "Velocity (05JAN1900 00 00 00).vrt",
          1100, "Depth (05JAN1900 10 00 00).vrt", "Velocity (05JAN1900 10 00 00).vrt",
          1450, "Depth (05JAN1900 20 00 00).vrt", "Velocity (05JAN1900 20 00 00).vrt",
          3000, "Depth (06JAN1900 06 00 00).vrt", "Velocity (06JAN1900 06 00 00).vrt",
          3200, "Depth (06JAN1900 16 00 00).vrt", "Velocity (06JAN1900 16 00 00).vrt",
          4500, "Depth (07JAN1900 02 00 00).vrt", "Velocity (07JAN1900 02 00 00).vrt",
          6000, "Depth (07JAN1900 12 00 00).vrt", "Velocity (07JAN1900 12 00 00).vrt",
          7050, "Depth (07JAN1900 22 00 00).vrt", "Velocity (07JAN1900 22 00 00).vrt",
          7600, "Depth (08JAN1900 08 00 00).vrt", "Velocity (08JAN1900 08 00 00).vrt",
          9600, "Depth (08JAN1900 18 00 00).vrt", "Velocity (08JAN1900 18 00 00).vrt") |>
  mutate(across(c(depth, velocity), function(x) file.path(basso_dir, x)))

basso_groups <-
  drive_file_by_id("1f_9Z8O64_t6mqicBW96KRBJ5ffOHqwen", vsizip=T) |>
  st_read(as_tibble=T)

basso_rast <- basso_filenames |>
  raster_prep_grid()

outpath <- here::here("data-raw", "results", "fsa_basso.Rds")

if(!file.exists(outpath)) {

  basso_result <- basso_rast |>
    raster_summarize_hsi(basso_groups, .group_var = comid) |>
    suitability_postprocess(basso_groups, .group_var = comid)

  basso_result |> saveRDS(outpath)

} else {

  basso_result <- readRDS(outpath)

}

#readRDS(here::here("data-raw", "results", "fsa_basso.Rds")) |>
#  mutate(flow_cfs = as.numeric(flow_cfs)) |>
#  saveRDS(here::here("data-raw", "results", "fsa_basso.Rds"))

# DEER CREEK -------------------------------------------------------------------

deer_dir <- here::here("data-raw", "temp", "deer-creek")
dir.create(deer_dir, recursive=T)
drive_file_by_id("1rmMw6PXJGS0-ui52eaotABCSlJsZOzvr") |>
  archive::archive_extract(deer_dir)

deer_filenames <-
  tribble(~flow_cfs, ~timestep,
          100, "15NOV2018 06 00 00",
          250, "15NOV2018 16 00 00",
          300, "16NOV2018 02 00 00",
          400, "16NOV2018 12 00 00",
          500, "16NOV2018 22 00 00",
          600, "17NOV2018 08 00 00",
          1000, "17NOV2018 18 00 00",
          3000, "18NOV2018 04 00 00",
          5000, "18NOV2018 14 00 00",
          6000, "19NOV2018 00 00 00",
          7000, "19NOV2018 10 00 00",
          9000, "19NOV2018 20 00 00",
          10000, "20NOV2018 06 00 00",
          11000, "20NOV2018 16 00 00",
          12000, "21NOV2018 02 00 00",
          13000, "21NOV2018 12 00 00",
          14000, "21NOV2018 22 00 00",
          15000, "22NOV2018 08 00 00") |>
  mutate(depth = file.path(deer_dir, glue::glue("Depth ({timestep}).vrt")),
         velocity = file.path(deer_dir, glue::glue("Velocity ({timestep}).vrt")))

deer_groups <-
  drive_file_by_id("1abatoUypK-oigD6zS6fUI7qmvEewZsSg", vsizip=T) |>
  st_read(as_tibble=T) |>
  janitor::clean_names() |>
  filter(comid > 0)

deer_rast <- deer_filenames |>
  raster_prep_grid()

outpath <- here::here("data-raw", "results", "fsa_deer.Rds")

if(!file.exists(outpath)) {

  deer_result <- deer_rast |>
    raster_summarize_hsi(deer_groups, .group_var = comid, parallel = FALSE) |>
    suitability_postprocess(deer_groups, .group_var = comid)

  deer_result |> saveRDS(outpath)

} else {

  deer_result <- readRDS(outpath)

}
