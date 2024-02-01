Width Datasets
================
Maddee Rubenson (FlowWest)
2024-02-01

## River Channel Width Dataset Exploration

### Global River Width Database

- <http://hydro.iis.u-tokyo.ac.jp/~yamadai/GWD-LR/>
- downloadable link here:<https://zenodo.org/records/1269595>

#### Vector files that overlap with the Central Valley

``` r
files_in_cv <- c('NK10.shp', 'NJ10.shp')

nk10 <- drive_file_by_id('1H-oAmip5Pp2H0sQJnGRT_srypQGy8dWE', vsizip=T) |>
  st_read()
```

    ## temp/NK10.zip already exists and will be used...

    ## Reading layer `NK10' from data source `/vsizip/temp/NK10.zip' using driver `ESRI Shapefile'
    ## Simple feature collection with 52767 features and 10 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: -124.4289 ymin: 40 xmax: -121.7742 ymax: 44
    ## Geodetic CRS:  WGS 84

``` r
nj10 <- drive_file_by_id('1xq-0uoJbCA5fQE38OsOK1mi9k3boXsmF', vsizip = T) |> 
  st_read()
```

    ## temp/NJ10.zip already exists and will be used...

    ## Reading layer `NJ10' from data source `/vsizip/temp/NJ10.zip' using driver `ESRI Shapefile'
    ## Simple feature collection with 29601 features and 10 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: -122.2447 ymin: 37.3486 xmax: -120.972 ymax: 40
    ## Geodetic CRS:  WGS 84

``` r
all_width_data_in_cv <- nk10 |> 
  bind_rows(nj10) |> 
  janitor::clean_names() |> 
  filter(lake_flag == 0) |> 
  glimpse()
```

    ## Rows: 62,126
    ## Columns: 11
    ## $ utm_east    <int> 594210, 594510, 594540, 594570, 594600, 594630, 594660, 59…
    ## $ utm_north   <dbl> 4683810, 4683240, 4683210, 4683180, 4683180, 4683150, 4683…
    ## $ width_m     <int> 3131, 3167, 3134, 3092, 2988, 2925, 2923, 2922, 2942, 2995…
    ## $ nchannels   <dbl> 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1…
    ## $ segment_id  <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ segment_ind <dbl> 592, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622…
    ## $ lake_flag   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ lon         <dbl> -121.8571, -121.8533, -121.8530, -121.8526, -121.8522, -12…
    ## $ lat         <dbl> 42.30069, 42.29552, 42.29525, 42.29511, 42.29497, 42.29470…
    ## $ elev_m      <int> 1262, 1268, 1268, 1268, 1268, 1268, 1268, 1268, 1268, 1268…
    ## $ geometry    <LINESTRING [°]> LINESTRING (-121.857 42.300..., LINESTRING (-12…

``` r
ggplot() +
    geom_sf(data = watersheds) +
    geom_sf(data = all_width_data_in_cv, aes(color = width_m))
```

![](width-datasets-exploration_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(all_width_data_in_cv, aes(x = width_m)) +
  geom_freqpoly(bins = 10)
```

![](width-datasets-exploration_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->
