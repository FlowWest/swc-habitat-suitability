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

# ### ORIGINAL VERSION - no longer used
#
# outpath <- here::here("data-raw", "results", "fsa_basso.Rds")
#
# if(!file.exists(outpath)) {
#
#   basso_result <- basso_rast |>
#     raster_summarize_hsi(basso_groups, .group_var = comid, hsi_func = raster_dvhsi_hqt) |>
#     suitability_postprocess(basso_groups, .group_var = comid)
#
#   basso_result |> saveRDS(outpath)
#
# } else {
#
#   basso_result <- readRDS(outpath)
#
# }

### VARIATION FOR REARING - BASEFLOW REMOVAL

# define baseflow mask as 150 cfs (COMID 2823750 observed dry season baseflow is 166 cfs)
basso_baseflow_mask <- terra::ifel(is.na(basso_rast$depth[["150"]]), 1, 0)

outpath <- here::here("data-raw", "results", "fsa_basso_nbfc.Rds")

if(!file.exists(outpath)) {

  basso_result_nbfc <- basso_rast |>
    raster_summarize_hsi(basso_groups, .group_var = comid, mask_raster = basso_baseflow_mask, hsi_func = raster_dvhsi_hqt) |>
    suitability_postprocess(basso_groups, .group_var = comid)

  basso_result_nbfc |> saveRDS(outpath)

} else {

  basso_result_nbfc <- readRDS(outpath)

}

### VARIATION FOR SPAWNING - ALTERNATIVE HSI

outpath <- here::here("data-raw", "results", "fsa_basso_spawning.Rds")

if(!file.exists(outpath)) {

  basso_result_spawning <- basso_rast |>
    raster_summarize_hsi(basso_groups, .group_var = comid, hsi_func = raster_dvhsi_spawning) |>
    suitability_postprocess(basso_groups, .group_var = comid)

  basso_result_spawning |> saveRDS(outpath)

} else {

  basso_result_spawning <- readRDS(outpath)

}

# DEER CREEK -------------------------------------------------------------------

deer_dir <- here::here("data-raw", "temp", "deer-creek")
dir.create(deer_dir, recursive=T)
drive_file_by_id("1rmMw6PXJGS0-ui52eaotABCSlJsZOzvr", dir=deer_dir) |>
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

# ### ORIGINAL VERSION - no longer used
#
# outpath <- here::here("data-raw", "results", "fsa_deer.Rds")
#
# if(!file.exists(outpath)) {
#
#   deer_result <- deer_rast |>
#     raster_summarize_hsi(deer_groups, .group_var = comid, hsi_func = raster_dvhsi_hqt) |>
#     suitability_postprocess(deer_groups, .group_var = comid)
#
#   deer_result |> saveRDS(outpath)
#
# } else {
#
#   deer_result <- readRDS(outpath)
#
# }

### VARIATION FOR REARING - BASEFLOW REMOVAL

# define baseflow mask as 150 cfs (COMID 8020924 observed dry season baseflow is 95 cfs)
deer_baseflow_mask <- terra::ifel(is.na(deer_rast$depth[["100"]]), 1, 0)

outpath <- here::here("data-raw", "results", "fsa_deer_nbfc.Rds")

if(!file.exists(outpath)) {

  deer_result_nbfc <- deer_rast |>
    raster_summarize_hsi(deer_groups, .group_var = comid, mask_raster = deer_baseflow_mask, hsi_func = raster_dvhsi_hqt) |>
    suitability_postprocess(deer_groups, .group_var = comid)

  deer_result_nbfc |> saveRDS(outpath)

} else {

  deer_result_nbfc <- readRDS(outpath)

}

### VARIATION FOR SPAWNING - ALTERNATIVE HSI

outpath <- here::here("data-raw", "results", "fsa_deer_spawning.Rds")

if(!file.exists(outpath)) {

  deer_result_spawning <- deer_rast |>
    #list(depth=deer_rast$depth$`10000`, velocity=deer_rast$velocity$`10000`) |>
    raster_summarize_hsi(deer_groups, .group_var = comid, hsi_func = raster_dvhsi_spawning) |>
    suitability_postprocess(deer_groups, .group_var = comid)

  deer_result_spawning |> saveRDS(outpath)

} else {

  deer_result_spawning <- readRDS(outpath)

}

# MOKELUMNE RIVER (EBMUD) ------------------------------------------------------

moke_dir <- here::here("data-raw", "temp", "mokelumne-river")
dir.create(moke_dir, recursive=T)

moke_ids <- c(R2 = "1BLE0U1lZkqgSb9QM0QlPo8VpmSXOnshe",
              R3 = "1m_ivN7Sr05agsp8EMg6Zn6MFMJoM8ifS", #1__pZkEuGf6J4bHUcZbY2XWx6Dl9vd35w",
              R4 = "1i2GUrr0Y2AoI9jNIYYbUqU4tohOD1Lwu", #160FV8cy3TceAo_1GxqhDs0VTgGS1kvAu",
              R5 = "1Dhmyq-9zbYn_T-jTXu5_hYdM-sxs3r0y", #19qrFSDLzpcdejPENh4qCUlO4sx6_IS_C",
              R6 = "1c9QsKaK7F6DVIiNFpw3YnLCDjUV-Eav8", #1jq7-hNmJXUNWa_TLqTe0WbI_3-P1gtYo",
              R7 = "1c_Phhz-46YpCCdIwVYGLJAp-xWKH43Y8") #1wyG0nk2ZlZivhYwnfakiIrySeAJQJi6k")
              #R2_4410 = "1dwS1myvE9eQZ2jCwr0kM5Mtau9N1ZAqB",
              #R2_8810 = "1Y3FCsy8AvNy43lOOkqTeS0-xJ9m7qrlH",
              #R7new =  "1CPyNYb0iMnLHhrJhhgGNYVUSl8Y5CQ0P" / "1oenUfh3vDva94U0Slj5naT4TTs8KPzS5"
lapply(names(moke_ids), function(x) {
  dir.create(file.path(moke_dir, x), recursive=T)
  drive_file_by_id(moke_ids[x], dir = moke_dir) |>
    archive::archive_extract(file.path(moke_dir, x))
  })

moke_groups <-
  drive_file_by_id("1lgD57oZyKP_vEUHVGuVEiNZ-L0WkFEvs", dir = moke_dir) |>
  read_sf() |>
  janitor::clean_names() |>
  filter(comid > 0)

moke_filenames <-
  read_csv(here::here("data-raw", "source", "hydraulic_model_data", "mokelumne_ebmud_2021", "mok_model_timestamps.csv")) |>
  mutate(depth = file.path(moke_dir, dir_id, glue::glue("Depth ({time_stamp_exact}).vrt")),
         velocity = file.path(moke_dir, dir_id, glue::glue("Velocity ({time_stamp_exact}).vrt"))) |>
  filter(dir_id %in% names(moke_ids))

moke_rasters <-
  moke_filenames |>
  select(dir_id, flow_cfs = q_in, depth, velocity) |>
  nest(filenames = c(flow_cfs, depth, velocity)) |>
  mutate(rasters = map(filenames, raster_prep_grid))

moke_reaches <-
  drive_file_by_id("1MO6rIySOJ1DcsvXh8YvNHF5ISyXhgNaX", dir = moke_dir) |>
  read_sf() |>
  janitor::clean_names() |>
  mutate(dir_id = str_replace(name, "Mok", "")) |>
  select(dir_id, geometry)

moke_groups_by_reach <- moke_groups |>
  st_intersection(moke_reaches) |>
  select(dir_id, comid, geometry) |>
  nest(groups = c(comid, geometry)) |>
  deframe()

moke_rasters_cropped <- moke_rasters |>
  mutate(rasters_cropped = pmap(list(dir_id, rasters), function(x, ras) {
    e <- terra::ext(terra::vect(moke_groups_by_reach[[x]]))
    list(depth = ras$depth |> terra::crop(e),
         velocity = ras$velocity |> terra::crop(e))
  })) |>
    select(dir_id, rasters_cropped) |>
    deframe()

moke_results_nbfc <- list()
moke_results_spawning <- list()

for (r in moke_rasters$dir_id) {

  message(r)

  #moke_rast <- moke_rasters$rasters[[which(moke_rasters$dir_id==r)]]
  moke_rast <- moke_rasters_cropped[[r]]

  ### VARIATION FOR REARING - BASEFLOW REMOVAL ---------------------------------

  # use the minimum flow available: 125 or 100 depending on the reach
  moke_baseflow_mask <- terra::ifel(is.na(moke_rast$depth[[1]]), 1, 0)

  outpath <- here::here("data-raw", "temp", glue::glue("fsa_moke_nbfc_{r}.Rds"))

  if(!file.exists(outpath)) {

    moke_results_nbfc[[r]] <- moke_rast |>
      raster_summarize_hsi(moke_groups_by_reach[[r]], .group_var = comid, mask_raster = moke_baseflow_mask, hsi_func = raster_dvhsi_hqt) |>
      suitability_postprocess(moke_groups_by_reach[[r]], .group_var = comid)

    moke_results_nbfc[[r]] |> saveRDS(outpath)

  } else {

    moke_results_nbfc[[r]] <- readRDS(outpath)

  }

  ### VARIATION FOR SPAWNING - ALTERNATIVE HSI ---------------------------------

  outpath <- here::here("data-raw", "temp", glue::glue("fsa_moke_spawning_{r}.Rds"))

  if(!file.exists(outpath)) {

    moke_results_spawning[[r]] <- moke_rast |>
      raster_summarize_hsi(moke_groups_by_reach[[r]], .group_var = comid, hsi_func = raster_dvhsi_spawning) |>
      suitability_postprocess(moke_groups_by_reach[[r]], .group_var = comid)

    moke_results_spawning[[r]] |> saveRDS(outpath)

  } else {

    moke_results_spawning[[r]] <- readRDS(outpath)

  }

}

combine_reaches <- function(data) {
  tibble(dir_id = names(data),
         result = data) |>
  unnest(result) |>
  filter(!is.nan(area_tot)) |>
  group_by(comid) %>%
  complete(flow_cfs = unique(.$flow_cfs)) |>
  arrange(comid, flow_cfs) |>
  mutate(across(c(area_tot, area_wua, length_ft),
                function(x) zoo::na.approx(x, x = flow_cfs, rule = 2))) |>
  group_by(flow_cfs, comid) |>
  summarize(area_tot = sum(area_tot),
            area_wua = sum(area_wua),
            area_pct = sum(area_wua) / sum(area_tot),
            length_ft = sum(length_ft),
            ind_per_lf = sum(area_tot) / sum(length_ft),
            wua_per_lf = sum(area_wua) / sum(length_ft)) |>
  ungroup() |>
  arrange(comid, flow_cfs)
}

moke_result_nbfc_combined <- combine_reaches(moke_results_nbfc)
moke_result_spawning_combined <- combine_reaches(moke_results_spawning)

moke_result_nbfc_combined |>
  saveRDS(here::here("data-raw", "results", glue::glue("fsa_moke_nbfc.Rds")))

moke_result_spawning_combined |>
  saveRDS(here::here("data-raw", "results", glue::glue("fsa_moke_spawning.Rds")))

