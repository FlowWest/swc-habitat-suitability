## Flow Summary
library(tidyverse)

options(scipen=999)

googledrive::drive_auth()

source('data-functions.R')

fcodes <-
  drive_file_by_id("1DwrieA-jrEsY8hnuEMFAqeAjsu7sB1D7") |>
  foreign::read.dbf() |>
  as_tibble() |>
  janitor::clean_names() |>
  rename(fcode = f_code, fcode_desc = descriptio) |>
  mutate(ftype_desc = map(fcode_desc, function(x) str_split(x, ":", simplify=TRUE)[1])) |>
  unnest(ftype_desc) |>
  select(fcode, ftype_desc)

# flowline shapefile attribute table
drive_file_by_id("1UiG8AeMr6mFOw7Jx--LyNRzez7GsDhzK") |>
  archive::archive_extract(dir = "temp", file = "NHDFlowline.dbf")

flowline_table <-
  foreign::read.dbf("temp/NHDFlowline.dbf") |>
  janitor::clean_names() |>
  select(comid, reachcode, gnis_id, gnis_name, lengthkm, ftype, fcode) |>
  mutate(huc_8 = substr(reachcode, 1, 8),
         huc_10 = substr(reachcode, 1, 10),
         huc_12 = substr(reachcode, 1, 12)) |>
  inner_join(fcodes |> select(fcode, ftype_desc))

nf <-
  drive_file_by_id("1P9vH9VXbPV9jYuFq35yF3elBRatYKQmU", vsizip=F) |>
  archive::archive_read() |>
  read_csv() |>
  filter(comid %in% flowline_table$comid)

flowline_attributes <- readRDS('../data/flowline_attributes.Rds')

# this summary table will be mapped to rivers in order to understand flow ranges of each river  that is modeled
flow_summary <- nf |>
  filter(!is.na(gage_id)) |>
  filter(ffm %in% c('peak_10', 'ds_mag_50')) |>
  select(-observed_years, -source, -alteration, -observed_year_start, -observed_year_end) |>
  inner_join(flowline_attributes) |>
  select(gage_id, ffm, p10, p25, p50, p75, p90, gnis_name) |>
  mutate(gnis_name = as.character(gnis_name)) |>
  group_by(gage_id, ffm, gnis_name) |>
  summarise_if(is.numeric, list(min = min, max = max, mean = mean)) |>
  glimpse()


watersheds <- DSMflow::watershed_ordering |> pull(watershed)

# NOTE: there are missing watersheds and it's difficult to know which gage would be best to use
flow_summary_watersheds <- flow_summary |> filter(gnis_name %in% watersheds)

# compare to CVPIA --------------------------------------------------------
# For the time being, let's use CALSIM 2018/2019 flows as thresholds since they are
# realistic with the current functionings of the Central Valley Rivers

flow_2018_2019 <- DSMflow::flows_cfs$biop_itp_2018_2019  |>
  pivot_longer(cols = c(`Upper Sacramento River`:`San Joaquin River`), values_to = 'cfs', names_to = "watershed") |>
  filter(!(watershed %in% c('Sutter Bypass', 'Yolo Bypass'))) |>
  group_by(watershed) |>
  summarise_if(is.numeric, list(min_cfs = min, max_cfs = max, mean_cfs = mean)) |>
  glimpse()

saveRDS(flow_2018_2019, "../data/watershed_flow_summary.RDS")

