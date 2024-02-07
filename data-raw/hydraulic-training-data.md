Hydraulic Training Data
================

``` r
library(tidyverse)
library(sf)
library(stars)
```

Define habitat suitability functions

``` r
# simple linear interpolation function
linterp <- function(x, x1, x2, y1, y2){
  y1 + ((x-x1)/(x2-x1)) * (y2-y1)
}

# vectorized function to calculate the depth*velocity suitability factor
# frac=FALSE uses simple 1/0 thresholds; frac=TRUE uses ranges defined in Yuba tables 8-9
dvhsi <- function(d, v, frac=FALSE) {
  if(frac){
    dhsi <- case_when(
      d<=0.5 ~ 0,
      d<=0.9 ~ linterp(d, 0.5, 0.9, 0, 1),
      d<=4.0 ~ 1,
      d<=5.2 ~ linterp(d, 4.0, 5.2, 1, 0),
      d>5.2 ~ 0
    )
    vhsi <- case_when(
      v<=0.1 ~ linterp(v, 0.0, 0.1, 0, 1),
      v<=0.8 ~ 1,
      v<=1.8 ~ linterp(v, 0.8, 1.8, 1, 0.35),
      v<=4.0 ~ linterp(v, 1.8, 4.0, 0.35, 0),
      v>4.0 ~ 0
    )
    return(sqrt(dhsi*vhsi)) 
  }else{
    return(if_else(d>0.5 & d<=5.2 & v>0.1 & v<=4.0, 1, 0))
  }
}
```

## Stanislaus

Import SRH2D model data for Stanislaus (unlike HEC-RAS, the SRH2D
outputs are natively vector format)

``` r
# SRH2D model domain, dissolved from SRH2D mesh faces using QGIS
stan_domain <- st_read("/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Domain.shp.zip", as_tibble=T)
```

    ## Reading layer `StanMesh072313_Domain' from data source 
    ##   `/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Domain.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1 feature and 1 field
    ## Geometry type: POLYGON
    ## Dimension:     XYZ
    ## Bounding box:  xmin: 655231.1 ymin: 4169204 xmax: 705835.6 ymax: 4188522
    ## z_range:       zmin: 0 zmax: 0
    ## Projected CRS: NAD83 / UTM zone 10N

``` r
# SRH2D model domain, manually split into polygons aligning with COMID segments
stan_comid <- st_read("/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Domain_COMID.shp.zip", as_tibble=T) |>
  st_zm() |> janitor::clean_names()
```

    ## Reading layer `StanMesh072313_Domain_COMID' from data source 
    ##   `/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Domain_COMID.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 54 features and 1 field
    ## Geometry type: POLYGON
    ## Dimension:     XYZ
    ## Bounding box:  xmin: 655231.1 ymin: 4169204 xmax: 705835.6 ymax: 4188522
    ## z_range:       zmin: 0 zmax: 0
    ## Projected CRS: NAD83 / UTM zone 10N

``` r
# SRH2D mesh vertices, converted from 2DM using QGIS 
stan_vertices <- st_read("/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Vertices.shp.zip", as_tibble=T) |>
  mutate(vid = row_number()) |>
  select(vid)
```

    ## Reading layer `StanMesh072313_Vertices' from data source 
    ##   `/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Vertices.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 241395 features and 1 field
    ## Geometry type: POINT
    ## Dimension:     XYZ
    ## Bounding box:  xmin: 655231.1 ymin: 4169204 xmax: 705835.6 ymax: 4188522
    ## z_range:       zmin: 2.461469 zmax: 151.0898
    ## Projected CRS: NAD83 / UTM zone 10N

``` r
# Thiessen (aka Voronoi) polygons for mesh vertices, generated using QGIS
stan_thiessen <- st_read("/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Thiessen.shp.zip", as_tibble=T) |>
  #mutate(vid = row_number()) # row order doesn't match the SRH2D outputs so need to spatial join
  st_join(stan_vertices, join=st_nearest_feature) |>
  select(vid) |>
  arrange(vid)
```

    ## Reading layer `StanMesh072313_Thiessen' from data source 
    ##   `/vsizip/hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_Thiessen.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 241395 features and 1 field
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 655231.1 ymin: 4169204 xmax: 705835.6 ymax: 4188522
    ## Projected CRS: NAD83 / UTM zone 10N

``` r
# confirm correct join via:
# ggplot() + geom_sf(data=stan_thiessen|>filter(vid<50)) + geom_sf(data=stan_vertices|>filter(vid<50)) 

# Bed elevations, extracted from mesh using QGIS "Export time series values from points of a mesh dataset"
stan_elev <- read_csv("hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_BedElevation.csv.gz") |>
  janitor::clean_names() |>
  mutate(vid = row_number()) |>
  select(vid, bed_elevation)

# Material IDs (Manning's roughness classes) extracted from mesh using QGIS "Export time series values from points of a mesh dataset"
stan_material <- 
  read_csv("hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_MaterialID.csv.gz") |>
  janitor::clean_names() |>
  mutate(vid = row_number()) |>
  select(vid, material_id)
# alternate approach with spatial join yields the same result
# stan_material <- 
#   read_csv("hydraulic_model_data/stanislaus_srh2d_2013/StanMesh072313_MaterialID.csv.gz") |>
#   janitor::clean_names() |>
#   st_as_sf(coords=c("x","y"), crs=st_crs(stan_thiessen)) |>
#   st_join(stan_thiessen, join=st_nearest_feature)

# SRH2D results by vertex point, incl. depth, velocity, shear stress, froude
stan_result_filenames <- c(
   "500" = "hydraulic_model_data/stanislaus_srh2d_2013/500cfs_072313.csv.gz",
   "750" = "hydraulic_model_data/stanislaus_srh2d_2013/750cfs_072413.csv.gz",
  "1000" = "hydraulic_model_data/stanislaus_srh2d_2013/1000cfs_072413.csv.gz",
  "1250" = "hydraulic_model_data/stanislaus_srh2d_2013/1250cfs_090313.csv.gz",
  "1500" = "hydraulic_model_data/stanislaus_srh2d_2013/1500cfs_071013.csv.gz",
  "1750" = "hydraulic_model_data/stanislaus_srh2d_2013/1750cfs_090313.csv.gz",
  "2250" = "hydraulic_model_data/stanislaus_srh2d_2013/2250cfs_121713.csv.gz",
  "3000" = "hydraulic_model_data/stanislaus_srh2d_2013/3000cfs_061913.csv.gz",
  "5000" = "hydraulic_model_data/stanislaus_srh2d_2013/5000cfs_071113.csv.gz"
)
stan_result_cols <- c("x_m"="n", "y_m"="n", "z_m"="n", 
                      "wse_m"="n", "wdepth_m"="n", 
                      "vel_x"="n", "vel_y"="n", "vel_mag"="n", 
                      "froude"="n", "stress"="n")

stan_result <- 
  names(stan_result_filenames) |>
  lapply(function(x) read_csv(stan_result_filenames[x], col_names=names(stan_result_cols), col_types=stan_result_cols, skip=1) |> 
           mutate(discharge_cfs = as.numeric(x), vid=row_number())) |>
  bind_rows() |>
  mutate(across(everything(), function(x) if_else(x==-999,NA,x))) |>
  mutate(depth_ft = wdepth_m*3.28084, velocity_fps = vel_mag*3.28084) |>
  select(discharge_cfs, vid, depth_ft, velocity_fps) |>
  left_join(stan_elev, by=join_by(vid)) |> 
  left_join(stan_material, by=join_by(vid)) |>
  mutate(
    cover_hs = 1, # cover habitat suitability, for now, actually need to get it from stan_material
    hsi_simp = dvhsi(depth_ft, velocity_fps, frac=F) * cover_hs,
    hsi_frac = dvhsi(depth_ft, velocity_fps, frac=T) * cover_hs,
  ) |>
  glimpse()
```

    ## Rows: 2,172,555
    ## Columns: 9
    ## $ discharge_cfs <dbl> 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 5…
    ## $ vid           <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1…
    ## $ depth_ft      <dbl> 0.0000000, 1.1168896, 3.6240659, 3.9757284, 3.9583468, 3…
    ## $ velocity_fps  <dbl> 0.0000000, 0.6719190, 1.5923429, 1.8671964, 1.2991604, 0…
    ## $ bed_elevation <dbl> 54.86779, 53.10716, 52.34171, 52.23186, 52.23380, 52.252…
    ## $ material_id   <dbl> 39, 28, 28, 28, 28, 25, 25, 28, 28, 28, 39, 39, 28, 28, …
    ## $ cover_hs      <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ hsi_simp      <dbl> 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,…
    ## $ hsi_frac      <dbl> 0.0000000, 1.0000000, 0.6964030, 0.5825029, 0.8219159, 1…

Import polygons created for each COMID reach, and calculate suitable
area for each

Proof of concept, first comid, 500 cfs:

``` r
test <- stan_thiessen |> 
  left_join(filter(stan_result, discharge_cfs==500), by=join_by(vid)) |>
  st_intersection(stan_comid[1]) |>
    mutate(area_m2 = units::drop_units(st_area(geometry)),
           wua_simp = hsi_simp * area_m2,
           wua_frac = hsi_frac * area_m2) |>
  glimpse()
```

    ## Warning: attribute variables are assumed to be spatially constant throughout
    ## all geometries

    ## Rows: 244,606
    ## Columns: 14
    ## $ vid           <int> 194330, 194471, 194472, 194473, 194475, 194476, 194477, …
    ## $ discharge_cfs <dbl> 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 5…
    ## $ depth_ft      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ velocity_fps  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ bed_elevation <dbl> 13.83647, 13.84992, 13.86708, 13.86089, 13.86002, 13.826…
    ## $ material_id   <dbl> 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, …
    ## $ cover_hs      <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ hsi_simp      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ hsi_frac      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ comid         <dbl> 2819852, 2819852, 2819852, 2819852, 2819852, 2819852, 28…
    ## $ geometry      <POLYGON [m]> POLYGON ((661138 4174738, 6..., POLYGON ((661145…
    ## $ area_m2       <dbl> 52.027195, 219.927514, 244.153508, 293.728856, 314.55599…
    ## $ wua_simp      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ wua_frac      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…

``` r
#test |> ggplot() + geom_sf(aes(fill=hsi_simp), color=NA)
```

Batch process, saving flow-to-suitable-area (fsa) by comid to file

``` r
stan_calc_hsi <- function(g, q){
  stan_thiessen |>
    left_join(filter(stan_result, discharge_cfs==q), by=join_by(vid)) |> 
    st_intersection(g) |>
    mutate(area_m2 = units::drop_units(st_area(geometry)),
           wua_simp = hsi_simp * area_m2,
           wua_frac = hsi_frac * area_m2) |>
    summarize(area_m2 = sum(area_m2),
              wua_simp = sum(wua_simp),
              wua_frac = sum(wua_frac)) |>
    mutate(pcthab_simp = wua_simp / area_m2,
           pcthab_frac = wua_frac / area_m2) |>
    st_drop_geometry()
}

if(!file.exists("../data/fsa_stanislaus.Rds")) {
  fsa_stanislaus <- 
    stan_comid |>
    expand_grid(discharge_cfs = as.numeric(names(stan_result_filenames))) |>
    mutate(result = map2(geometry, discharge_cfs, function(g, q) stan_calc_hsi(st_sfc(g, crs=st_crs(stan_comid)), q))) |>
    unnest_wider(result) |>
    st_as_sf() |>
    glimpse()
  fsa_stanislaus |> saveRDS("../data/fsa_stanislaus.Rds")
} else {
  fsa_stanislaus <- readRDS("../data/fsa_stanislaus.Rds") |> glimpse()
}
```

    ## Rows: 486
    ## Columns: 8
    ## $ comid         <dbl> 2819852, 2819852, 2819852, 2819852, 2819852, 2819852, 28…
    ## $ geometry      <POLYGON [m]> POLYGON ((661091.5 4174815,..., POLYGON ((661091…
    ## $ discharge_cfs <dbl> 500, 750, 1000, 1250, 1500, 1750, 2250, 3000, 5000, 500,…
    ## $ area_m2       <dbl> 2973283.79, 2973283.79, 2973283.79, 2973283.79, 2973283.…
    ## $ wua_simp      <dbl> 184235.315, 182968.634, 181451.715, 102791.552, 73867.99…
    ## $ wua_frac      <dbl> 124882.432, 123558.843, 121039.274, 62002.932, 46502.493…
    ## $ pcthab_simp   <dbl> 0.06196358, 0.06153756, 0.06102738, 0.03457173, 0.024843…
    ## $ pcthab_frac   <dbl> 0.04200152, 0.04155636, 0.04070895, 0.02085335, 0.015640…

``` r
#stan_hsi |> ggplot() + geom_sf(aes(fill=pcthab_simp), color=NA) + facet_wrap(~discharge_cfs)
```

``` r
fsa_stanislaus |> 
  ggplot() + geom_line(aes(x = discharge_cfs, y = pcthab_frac, color=factor(comid)))
```

![](hydraulic-training-data_files/figure-gfm/stan-plot-hsi-1.png)<!-- -->
