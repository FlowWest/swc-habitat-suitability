SWC Habitat Model Data Discovery
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-01-04

- [Data Import](#data-import)
- [Import HUC Watersheds](#import-huc-watersheds)
- [Import Supplemental Tables](#import-supplemental-tables)
  - [Import flowline geometry](#import-flowline-geometry)
  - [Import catchments](#import-catchments)
- [Stream Widths?](#stream-widths)

*Selected HUCs*

- Lower Yuba 18020107
- Upper Yuba 18020125

*Data Source*

- NHDPlusV2 dataset <https://nhdplus.com/NHDPlus/NHDPlusV2_home.php>
- As retrievedfrom
  <https://www.epa.gov/waterdata/nhdplus-california-data-vector-processing-unit-18>
- User Guide
  <https://www.epa.gov/system/files/documents/2023-04/NHDPlusV2_User_Guide.pdf>
- (This is different from NHDPlusHR which is generated at a higher
  resolution and more computationally intensive)

## Data Import

## Import HUC Watersheds

``` r
selected_huc_8 <- c("18020107", "18020125")

# HUC-12 watersheds and higher level heirarchies
watersheds <- st_read("nhdplus/WBD_Subwatershed.shp") |> 
  janitor::clean_names() |>
  filter(huc_8 %in% selected_huc_8) |>
  st_transform(project_crs)
```

    ## Reading layer `WBD_Subwatershed' from data source 
    ##   `C:\Users\skylerlewis\Github\swc-habitat-suitability\data-raw\nhdplus\WBD_Subwatershed.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 4564 features and 21 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -124.5351 ymin: 32.133 xmax: -114.6198 ymax: 43.34273
    ## Geodetic CRS:  NAD83

``` r
watersheds |> ggplot() + geom_sf()
```

![](data-discovery_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Import Supplemental Tables

``` r
# identifiers of each flowline (COMID)
flowline_comid_reachcode_crosswalk <- 
  foreign::read.dbf("nhdplus/NHDReachCode_Comid.dbf") |> 
  as_tibble() |> 
  janitor::clean_names() |>
  mutate(huc_8 = substr(reachcode, 1, 8),
         huc_10 = substr(reachcode, 1, 10),
         huc_12 = substr(reachcode, 1, 12)) |>
  filter(huc_8 %in% selected_huc_8)

# slopes and endpoint elevations
flowline_slopes <- 
  foreign::read.dbf("nhdplus/NHDPlusAttributes/elevslope.dbf") |> 
  as_tibble() |> 
  janitor::clean_names() |>
  mutate(slope = if_else(slope==-9998, NA, slope))

# cumulative upstream area
flowline_cumarea <- 
  foreign::read.dbf("nhdplus/NHDPlusAttributes/CumulativeArea.dbf") |> 
  as_tibble() |> 
  janitor::clean_names() |>
  rename(comid = com_id)

# flow routing attributes as described in https://www.usgs.gov/national-hydrography/value-added-attributes-vaas
flowline_vaattr <- 
  foreign::read.dbf("nhdplus/NHDPlusAttributes/PlusFlowlineVAA.dbf") |> 
  as_tibble() |> 
  janitor::clean_names()

# fcode table
fcodes <- 
  foreign::read.dbf("nhdplus/NHDFCode.dbf") |> 
  as_tibble() |> 
  janitor::clean_names() |>
  rename(fcode = f_code, fcode_desc = descriptio)
```

### Import flowline geometry

``` r
flowlines <- 
  st_read("nhdplus/Hydrography/NHDFlowline.shp") |>
  janitor::clean_names() |>
  mutate(huc_8 = substr(reachcode, 1, 8),
         huc_10 = substr(reachcode, 1, 10),
         huc_12 = substr(reachcode, 1, 12)) |>
  filter(huc_8 %in% selected_huc_8) |>
  inner_join(fcodes |> select(fcode, fcode_desc)) |>
  inner_join(flowline_slopes |> select(comid, slope, maxelevsmo, minelevsmo, slopelenkm), by = join_by(comid)) |>
  inner_join(flowline_cumarea, by = join_by(comid)) |>
  arrange(comid) |>
  st_transform(project_crs)
```

    ## Reading layer `NHDFlowline' from data source 
    ##   `C:\Users\skylerlewis\Github\swc-habitat-suitability\data-raw\nhdplus\Hydrography\NHDFlowline.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 178868 features and 14 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XYZM
    ## Bounding box:  xmin: -124.4096 ymin: 32.50006 xmax: -114.5885 ymax: 43.33627
    ## z_range:       zmin: 0 zmax: 0
    ## m_range:       mmin: 0 mmax: 100
    ## Geodetic CRS:  NAD83

    ## Joining with `by = join_by(fcode)`

``` r
# plot showing slopes
flowlines |> 
  st_zm() |>
  filter(gnis_name %in% c("Yuba River", "South Yuba River", "Middle Yuba River", "North Yuba River")) |>
  ggplot() + 
  geom_sf(data=st_zm(flowlines), aes(color = slope)) +
  geom_sf(aes(color = slope), linewidth=1) + 
  scale_color_viridis_c(trans = "log")
```

![](data-discovery_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Import catchments

``` r
# catchment associated with each flowline reach (COMID)
catchments <- 
  st_read("nhdplus/Catchment.shp") |> 
  janitor::clean_names() |>
  mutate(comid = as.numeric(featureid)) |>
  inner_join(flowlines |> st_drop_geometry() |> select(comid)) |>
  arrange(comid) |>
  st_transform(project_crs) 
```

    ## Reading layer `Catchment' from data source 
    ##   `C:\Users\skylerlewis\Github\swc-habitat-suitability\data-raw\nhdplus\Catchment.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 140835 features and 4 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -124.4098 ymin: 32.13295 xmax: -114.6198 ymax: 43.34269
    ## Geodetic CRS:  NAD83

    ## Joining with `by = join_by(comid)`

``` r
# examples
ggplot() + geom_sf(data = catchments, color="orange") + 
  geom_sf(data = watersheds, color="red", fill=NA) + 
  geom_sf(data = st_zm(flowlines), color="blue") 
```

![](data-discovery_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Stream Widths?

#### Points along line

using st_line_sample or st_segmentize

``` r
# approach with one point every 100 ft
sample_points_100ft <- 
  flowlines |> 
  filter(gnis_name %in% c("Yuba River", "South Yuba River", "Middle Yuba River", "North Yuba River")) |>
  st_line_sample(density = 1 / 100) 

# approach with just the midpoint of each segment
sample_points_midpt <- 
  flowlines |>
  filter(gnis_name %in% c("Yuba River", "South Yuba River", "Middle Yuba River", "North Yuba River")) |>
  #st_line_sample(sample = c(0, 0.25, 0.5, 0.75, 1))
  st_line_sample(sample = c(0.5))

sample_points_midpt |> ggplot() + geom_sf()
```

![](data-discovery_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
