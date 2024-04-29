Flow Data Analysis for Inundation Duration HSI Component
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-03-07

``` r
library(tidyverse)
library(sf)
library(stars)
library(lubridate)

knitr::opts_chunk$set(eval=TRUE, fig.width=6.5, fig.height=4, dpi=300)

theme_set(theme_minimal())
```

``` r
flowlines <- readRDS("../data/flowline_geometries.Rds") |>
  st_transform("ESRI:102039")

flowline_attributes <- readRDS("../data/flowline_attributes.Rds")

flow_to_suitable_area <- readRDS("../data/fsa_combined.Rds")
```

``` r
# functions from interoperable flows

retrieve_cdec_csv <- function(sta=character(), sen=character(), dur=character(), 
                     start_date=ymd("1900-10-01"), end_date=ymd("2023-09-30"),
                     dir="temp") {
  name <- str_to_upper(paste0(sta, "_", sen, "_", dur))
  filename <- file.path(dir, paste0(name,".csv.gz"))
  if(!file.exists(filename)){
    message(paste0("downloading to ",filename))
    dir.create(dir, recursive = TRUE)
    data_raw <-
      httr::GET(
      url="https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet",
      query=list(Stations=sta, SensorNums=sen, dur_code=dur, 
                 Start=format(start_date,"%Y-%m-%d"), End=format(end_date,"%Y-%m-%d"))) |> 
      httr::content("raw") 
    gzf <- gzcon(file(filename, "wb"))
    data_raw |> writeBin(gzf)
    close(gzf)
  } else {
    message(paste0(filename, " exists, loading"))
  }
  return(filename)
}

parse_cdec_csv <- function(filename) {
  read_csv(filename, 
           col_select = c(date_time = "DATE TIME", value = "VALUE"),
           col_types = list(col_datetime(), col_double())) |>
    janitor::clean_names() |>
    mutate(water_year = if_else(month(date_time)>=10, year(date_time)+1, year(date_time)))
}
```

Duration suitability criteria from the CBEC LYR analysis are:

**Valley Lowland**

- 1 to 17 days: 0.66
- 18 to 24 days: 1.00
- 25+ days: 0.66

**Valley Foothill**

- 1 to 9 days: 1.00
- 10+ days: 0.66

``` r
# refer to https://cdec.water.ca.gov/webgis/?appid=cdecstation
# refer to https://docs.google.com/spreadsheets/d/1OaohqeYbp-dTyIim3mFDYP11GEqIZWqSYIfT67b5-eE/edit

# select a representative best flow gage to use for each reach, ideally available 25+ years
# refer to https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=### start 10/1/1997 - 09/30/2023
# drainage areas from StreamStats gage-pages
# mean annual flow is joined in based on the comid that the gage is located in

fsa_combined <- readRDS("../data/fsa_combined.Rds")

# get list of modelled flows for each dataset
training_flows <- fsa_combined |>
  select(dataset, comid, flow_cfs) |>
  group_by(dataset, flow_cfs) |> 
  summarize() %>%
  split(.$dataset) |>
  lapply(function(df) df |> pull(flow_cfs))

# get flow information for gaged reaches, and download cdec flow time series
training_reach_gages <- tribble(
  ~river, ~usgs_id, ~cdec_sta, ~cdec_sen, ~cdec_dur, ~start_date, ~end_date, ~da_gauge_sq_mi, ~gauge_comid, ~flows,
  "Deer Creek",       11383500, "DVD", 20, "E", mdy("10/1/1997"), mdy("9/30/2023"), 208, 8020924, training_flows[["Deer Creek"]],
  "Lower Yuba River", 11421000, "MRY", 20, "E", mdy("10/1/1997"), mdy("9/30/2023"), 1339, 7981844, training_flows[["Lower Yuba River"]],
  "Stanislaus River", 11303000, "RIP", 20, "E", mdy("10/1/1997"), mdy("9/30/2023"), 1075, 2819818, training_flows[["Stanislaus River"]]
) |>
  mutate(da_gauge_sq_km = da_gauge_sq_mi * 1.609344^2) |>
  mutate(filename = pmap_chr(list(cdec_sta, cdec_sen, cdec_dur), 
                             function(x, y, z) retrieve_cdec_csv(sta=x, sen=y, dur=z))) |>
  left_join(flowline_attributes |> 
              transmute(comid, maf_gauge = erom_q_ma_cfs), by=join_by(gauge_comid == comid))

# get historical mean daily flow for gages
training_cdec_data <- training_reach_gages |> 
  select(river, filename) |> 
  deframe() |>
  lapply(parse_cdec_csv) |>
  bind_rows(.id = "river") |>
  mutate(q_gauge = as.numeric(value),
         q_gauge = case_when(q_gauge>=0 ~ q_gauge)) |>
  drop_na() |>
  group_by(river, date = date(date_time)) |>
  summarize(q_gauge = mean(q_gauge)) |>
  # just data for November through June, as was done with LYR analysis
  filter(month(date) %in% c(11, 12, 1, 2, 3, 4, 5, 6)) |>
  glimpse()
```

    ## Rows: 15,982
    ## Columns: 3
    ## Groups: river [3]
    ## $ river   <chr> "Deer Creek", "Deer Creek", "Deer Creek", "Deer Creek", "Deer …
    ## $ date    <date> 1997-11-01, 1997-11-02, 1997-11-03, 1997-11-04, 1997-11-05, 1…
    ## $ q_gauge <dbl> 68.28125, 65.77083, 64.25000, 64.32292, 63.40625, 66.07292, 83…

``` r
#  function to calculate exceedence by flow
calculate_exceedence <- function(q_gauge, q_crit) {
  
  n_days <- length(q_gauge)
  
  exceedence_intervals <- 
    tibble(exceeds_q = q_gauge > q_crit,
           exceedence_event = with(rle(exceeds_q), rep(seq_along(lengths), lengths))) |>
    group_by(exceedence_event) |>
    mutate(cumulative_exceedence = if_else(exceeds_q, seq_along(exceedence_event), 0)) |>
    filter(exceeds_q) |>
    summarize(duration_days = max(cumulative_exceedence)) |>
    mutate(durhsi_vl = case_when(duration_days == 0  ~ 0,
                                 duration_days >= 1  ~ 0.66,
                                 duration_days >= 18 ~ 1.00,
                                 duration_days >= 25 ~ 0.66),
           durhsi_vf = case_when(duration_days < 9  ~ 0,
                                 duration_days >= 9  ~ 1.00,
                                 duration_days >= 10 ~ 0.66))
  exceedence_intervals |>
    summarize(n_inundations = n(),
              avg_days_inundated = mean(duration_days),
              avg_durhsi_vl = mean(durhsi_vl),
              avg_durhsi_vf = mean(durhsi_vf),
              tot_days_inundated = sum(duration_days),
              tot_days_weighted_vl = sum(duration_days * durhsi_vl),
              tot_days_weighted_vf = sum(duration_days * durhsi_vf)) |>
    mutate(avg_daily_durhsi_vl = tot_days_weighted_vl / n_days,
           avg_daily_durhsi_vf = tot_days_weighted_vf / n_days) |>
    as.list()
}

# example usage
calculate_exceedence(training_cdec_data$q_gauge, 500)
```

    ## $n_inundations
    ## [1] 179
    ## 
    ## $avg_days_inundated
    ## [1] 54.98324
    ## 
    ## $avg_durhsi_vl
    ## [1] 0.66
    ## 
    ## $avg_durhsi_vf
    ## [1] 0.3351955
    ## 
    ## $tot_days_inundated
    ## [1] 9842
    ## 
    ## $tot_days_weighted_vl
    ## [1] 6495.72
    ## 
    ## $tot_days_weighted_vf
    ## [1] 9499
    ## 
    ## $avg_daily_durhsi_vl
    ## [1] 0.4064397
    ## 
    ## $avg_daily_durhsi_vf
    ## [1] 0.5943562

``` r
#training_cdec_data |>
#  nest(gauge_ts = c(date, q_gauge)) |>
#  expand_grid(model_q = seq(300,15000,100)) |>
#  mutate(result = map2(gauge_ts, model_q, 
#                       function(df, q) calculate_exceedence(df$q_gauge, q))) |>
#  select(-gauge_ts) |>
#  unnest_wider(result)

durhsi_by_model_q <- 
  training_reach_gages |> 
  left_join(training_cdec_data, by=join_by(river), relationship="one-to-many") |>
  select(river, flows, date, q_gauge) |>
  nest(gauge_ts = c(date, q_gauge)) |>
  unnest(flows) |> rename(model_q = flows) |>
  mutate(result = map2(gauge_ts, model_q, 
                       function(df, q) calculate_exceedence(df$q_gauge, q))) |>
  select(-gauge_ts) |>
  unnest_wider(result) |>
  glimpse()
```

    ## Rows: 55
    ## Columns: 11
    ## $ river                <chr> "Deer Creek", "Deer Creek", "Deer Creek", "Deer C…
    ## $ model_q              <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 5000, 6…
    ## $ n_inundations        <int> 55, 72, 76, 66, 54, 54, 19, 2, 0, 0, 0, 0, 0, 0, …
    ## $ avg_days_inundated   <dbl> 40.127273, 11.152778, 8.815789, 6.727273, 5.85185…
    ## $ avg_durhsi_vl        <dbl> 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, N…
    ## $ avg_durhsi_vf        <dbl> 0.4727273, 0.3055556, 0.2105263, 0.1515152, 0.185…
    ## $ tot_days_inundated   <dbl> 2207, 803, 670, 444, 316, 216, 57, 3, 0, 0, 0, 0,…
    ## $ tot_days_weighted_vl <dbl> 1456.62, 529.98, 442.20, 293.04, 208.56, 142.56, …
    ## $ tot_days_weighted_vf <dbl> 2130, 652, 473, 298, 198, 109, 18, 0, 0, 0, 0, 0,…
    ## $ avg_daily_durhsi_vl  <dbl> 0.4444980165, 0.1617271895, 0.1349404944, 0.08942…
    ## $ avg_daily_durhsi_vf  <dbl> 0.649984742, 0.198962466, 0.144339335, 0.09093683…

*To apply this output:*

- Generate a summarized depth result raster grid or mesh table that
  indicates the first flow (model_q) at which the cell is inundated.
- Join the durhsi_by_model_q table so that we have a raster grid or mesh
  table of variable “first inundated flow” avg_durhsi values
- At each flow in the flow-to-suitable-area calculation script, multiply
  the grid/vector of “first inundated avg durhsi” value for the cell
  with a grid/vector of 1 = inundated or 0 = not inundated. The result
  is the “Duration HSI” component that is multiplied with the
  depth/velocity based HSI

For example, if a cell is first inundated at 500 cfs:

- In the 100cfs raster, DurHSI = 0
- In the 250cfs raster, DurHSI = 0
- In the 500cfs raster, DurHSI = the avg_durhsi value for 500 cfs
- In the 1000cfs raster, DurHSI = the avg_durhsi value for 500 cfs

The avg_durhsi is based on the typical duration of an individual
inundation period – this is the representative DurHSI value to use in
the HSI calculations. The avg_daily_durhsi also accounts for how
frequently inundation occurs and may be useful in other contexts.

``` r
durhsi_by_model_q |> saveRDS("../data/durhsi_by_model_q.Rds")
```

``` r
knitr::knit_exit()
```
