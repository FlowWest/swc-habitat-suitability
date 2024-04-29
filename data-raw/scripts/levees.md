Levees
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-02-16

``` r
# transect functions developed previously in terrain-attributes.Rmd

perpendicular_transect <- function(pair=st_multipoint(), length=numeric()) {
  x1 <- pair[1]
  x2 <- pair[2]
  y1 <- pair[3]
  y2 <- pair[4]
  midpoint_x <- (x1+x2)/2
  midpoint_y <- (y1+y2)/2
  orig_bearing <- atan2((y2-y1), (x2-x1))
  perp_bearing <- orig_bearing + pi/2
  radius <- length/2
  perp_x1 <- midpoint_x + radius*cos(perp_bearing)
  perp_y1 <- midpoint_y + radius*sin(perp_bearing)
  perp_x2 <- midpoint_x - radius*cos(perp_bearing)
  perp_y2 <- midpoint_y - radius*sin(perp_bearing)
  return(st_linestring(c(st_point(c(perp_x1,perp_y1)), 
                         st_point(c(perp_x2, perp_y2)))))
}

perpendicular_transects <- function(ls=st_linestring(), length=numeric()) {
  if ("sf" %in% class(ls)){
    ls <- st_as_sfc(ls, crs=st_crs(ls))
  }
  if ("sfc" %in% class(ls)){
    m <- t(st_zm(ls[[1]]))
  } 
  if ("sfg" %in% class(ls)){
    m <- t(ls)
  }
  point_pairs <- lapply(seq(1, length(m)-3, 2), function(i) c(st_point(c(m[i], m[i+1])), st_point(c(m[i+2], m[i+3]))))
  perp_lines <- lapply(point_pairs, function(pair) perpendicular_transect(pair, length))
  return(st_sfc(perp_lines, crs=st_crs(ls)))
}
```

``` r
flow <- 
  drive_file_by_id("1UiG8AeMr6mFOw7Jx--LyNRzez7GsDhzK", vsizip=T) |>
  st_read() |>
  janitor::clean_names() |>
  select(comid) |>
  st_transform(project_crs) |>
  st_zm()
```

    ## ! Using an auto-discovered, cached token.

    ##   To suppress this message, modify your code or options to clearly consent to
    ##   the use of a cached token.

    ##   See gargle's "Non-interactive auth" vignette for more details:

    ##   <https://gargle.r-lib.org/articles/non-interactive-auth.html>

    ## ℹ The googledrive package is using a cached token for 'slewis@flowwest.com'.

    ## Auto-refreshing stale OAuth token.

    ## temp/NHDFlowline.zip already exists and will be used...

    ## Reading layer `NHDFlowline' from data source `/vsizip/temp/NHDFlowline.zip' using driver `ESRI Shapefile'
    ## Simple feature collection with 178868 features and 14 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XYZM
    ## Bounding box:  xmin: -124.4096 ymin: 32.50006 xmax: -114.5885 ymax: 43.33627
    ## z_range:       zmin: 0 zmax: 0
    ## m_range:       mmin: 0 mmax: 100
    ## Geodetic CRS:  NAD83

``` r
levees <- 
  st_read("/vsizip/levees/nld_ca_levees.shp.zip") |>
  janitor::clean_names() |>
  st_transform(project_crs)
```

    ## Reading layer `nld_ca_levees' from data source 
    ##   `/vsizip/levees/nld_ca_levees.shp.zip' using driver `ESRI Shapefile'
    ## Simple feature collection with 1766 features and 58 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: -124.3296 ymin: 32.54026 xmax: -114.5009 ymax: 41.99989
    ## Geodetic CRS:  WGS 84

``` r
# version 1: simple minimum distance
if(FALSE){
  levee_attr_distance <- flow |> 
    mutate(distance_to_levee = st_distance(geometry, levees)) |>
    st_drop_geometry()
}

# version two: overlapping buffers
if(FALSE){
  radius <- 1000 #ft
  
  levees_buffer <- levees |> 
    st_buffer(dist = radius,
              endCapStyle="FLAT",
              joinStyle="ROUND")
  
  levee_attr_pct_buffer <- 
    st_buffer(dist = radius,
              endCapStyle="FLAT",
              joinStyle="ROUND") |>
    mutate(denom = st_area(geometry)) |>
    st_intersection(levees_buffer) |>
    mutate(num = st_area(geometry)) |>
    mutate(pct_leveed = num/denom) |>
    st_drop_geometry()
}

# version three: using transects
ls <- (flow |> filter(comid==15034139) |> pull(geometry))[[1]]
ls <- (flow |> head(1) |> pull(geometry))[[1]]

frac_leveed <- function(ls=st_linestring(), radius=1000){
  # crs of st_linestring must match crs of levees layer
  
  ls_densified <- ls |> smoothr::densify(n=3)
  transects <- perpendicular_transects(ls, radius) |> st_sfc(crs=st_crs(levees)) |> st_sf()
  n_transects <- nrow(transects)
  transect_lengths <- st_length(transects) |> units::set_units("ft") |> units::drop_units()
  
  intersected <- st_filter(transects, levees)
  n_intersected <- nrow(intersected)
  frac_leveed_longitudinal <- n_intersected / n_transects
  
  # transects_split <- st_collection_extract(lwgeom::st_split(transects, levees),"LINESTRING")
  # within_levees <- st_filter(transects_split, ls)
  # lengths_within_levees <- st_length(within_levees) |> units::set_units("ft") |> units::drop_units()
  # frac_leveed_lateral <- sum(lengths_within_levees) / sum(transect_lengths)
  
  if (n_intersected>0) {
    intersected_split <- st_collection_extract(lwgeom::st_split(intersected, levees),"LINESTRING")
    intersected_within_levees <- st_filter(intersected_split, ls)
    lengths_within_levees <- st_length(intersected_within_levees) |> units::set_units("ft") |> units::drop_units()
    lateral_levee_confinement_ft <- median(lengths_within_levees) 
  } else {
    lateral_levee_confinement_ft <- NA
  }
  
  return(list("frac_leveed_longitudinal" = frac_leveed_longitudinal,
              "lateral_levee_confinement_ft" = lateral_levee_confinement_ft))
}

if(!file.exists("../data/attr_frac_leveed.Rmd")) {
  attr_frac_leveed <- 
    flow |> 
    # head(2) |>
    # filter(comid==15034139) |>
    mutate(result = map(geometry, frac_leveed)) |> 
    unnest_wider(result) |>
    select(-geometry)
  
  attr_frac_leveed |> saveRDS("../data/attr_frac_leveed.Rmd")

} else {

  attr_frac_leveed <- readRDS("../data/attr_frac_leveed.Rmd")
}

attr_frac_leveed |> filter(frac_leveed_longitudinal>0) |> glimpse()
```

    ## Rows: 8,908
    ## Columns: 3
    ## $ comid                        <int> 22226680, 22226632, 22226660, 24680675, 2…
    ## $ frac_leveed_longitudinal     <dbl> 0.21621622, 1.00000000, 0.62500000, 0.163…
    ## $ lateral_levee_confinement_ft <dbl> 1751.7058, 2907.2030, 2616.2554, 1768.276…

``` r
flow |>
  inner_join(attr_frac_leveed |> filter(frac_leveed_longitudinal>0), by=join_by(comid)) |>
  plot()
```

![](levees_files/figure-gfm/levees-1.png)<!-- -->
