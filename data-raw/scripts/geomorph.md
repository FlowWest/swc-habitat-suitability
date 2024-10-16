Predictor Data Prep: Geomorphology and Soils
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-08-16

- [Soils](#soils)
- [UCD Geomorphology classes
  (experimental)](#ucd-geomorphology-classes-experimental)
  - [version 1](#version-1)
  - [version 2](#version-2)
- [UCD Geomorph Classes](#ucd-geomorph-classes)

## Soils

``` r
soils <- 
  st_read(file.path("/vsizip", here::here("data-raw", "source", "soils", "Cal_STATSGO2.shp.zip")), as_tibble=T) |>
  st_transform(project_crs) |>
  janitor::clean_names() |> 
  transmute(soil_drainage = factor(drainage, levels=c("Very poorly drained",
                                                      "Poorly drained",
                                                      "Somewhat poorly drained",
                                                      "Well drained",
                                                      "Somewhat excessively drained",
                                                      "Excessively drained")),
            soil_hsg = factor(hydro_g, levels=c("A","B","C","D")),
            soil_salt_accum = factor(salt, levels=c("acid","nonacid","euic","calcareous")),
            soil_climate = factor(climate, levels=c("frigid","mesic","thermic","hyperthermic","isofrigid","isomesic","isothermic","isohyperthermic")),
            soil_minerology = factor(minerology, levels=c("superactive","active","semiactive","subactive")),
            soil_volcanic = if_else(str_detect(texture,"ashy") | str_detect(texture,"cindery") | str_detect(texture,"medial"), 1, 0),
            soil_texture_simple = case_when(str_detect(word(texture,1), "ashy") ~ "ashy",
                                            str_detect(word(texture,1), "medial") ~ "medial",
                                            str_detect(word(texture,1), "clayey") ~ "clayey",
                                            str_detect(word(texture,1), "clayey") ~ "silty",
                                            str_detect(word(texture,1), "clayey") ~ "sandy",
                                            str_detect(word(texture,1), "clayey") ~ "loamy",
                                            str_detect(word(texture,1), "clayey") ~ "fine"
                                            ) |> factor())
```

    ## Reading layer `Cal_STATSGO2' from data source 
    ##   `/vsizip/C:/Users/skylerlewis/Github/swc-habitat-suitability/data-raw/source/soils/Cal_STATSGO2.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 4235 features and 19 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -13849240 ymin: 3833661 xmax: -12705030 ymax: 5162405
    ## Projected CRS: WGS 84 / Pseudo-Mercator

``` r
soils_comid <-
  flowlines_gcs |>
  select(comid) |>
  st_point_on_surface() |>
  st_transform(st_crs(soils)) |>
  st_join(soils) |>
  st_drop_geometry() |>
  glimpse()
```

    ## Warning: st_point_on_surface assumes attributes are constant over geometries

    ## Warning in st_point_on_surface.sfc(st_geometry(x)): st_point_on_surface may not
    ## give correct results for longitude/latitude data

    ## Rows: 178,868
    ## Columns: 8
    ## $ comid               <int> 20245062, 24085230, 22226684, 22226720, 22226732, …
    ## $ soil_drainage       <fct> NA, NA, Well drained, Well drained, Well drained, …
    ## $ soil_hsg            <fct> NA, NA, B, B, B, B, B, D, D, D, D, D, D, D, D, D, …
    ## $ soil_salt_accum     <fct> NA, NA, nonacid, nonacid, NA, nonacid, nonacid, NA…
    ## $ soil_climate        <fct> NA, NA, mesic, mesic, isomesic, mesic, mesic, isom…
    ## $ soil_minerology     <fct> NA, NA, superactive, superactive, NA, superactive,…
    ## $ soil_volcanic       <dbl> NA, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ soil_texture_simple <fct> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…

``` r
soils_comid |> saveRDS(here::here("data-raw", "results", "attr_soils.Rds"))
```

## UCD Geomorphology classes (experimental)

This is a very rough generalization of spot classifications to the
entire basin. Use with caution!

``` r
geomorph_site_data <- 
  tibble(result=rjson::fromJSON(file=here::here("data-raw", "source", "ucd_geomorph", "geoSites.json"))) |> 
  unnest_wider(result) |> 
  mutate(geomorph_region = map_chr(geoClass, function(x) x$geoRegion$abbreviation),
         geomorph_class = map_chr(geoClass, function(x) x$name),
         geomorph_class_median_attributes = map(geoClass, function(x) x$medianAttributes),
         geometry = map(geometry, function(x) c(x$coordinates[[2]], x$coordinates[[1]]) |> st_point()) |> st_sfc(crs="EPSG:4326")) |>
  janitor::clean_names() |> 
  st_as_sf() |>
  filter(geomorph_region=="SAC") |>
  select(identity, geometry, geomorph_class) |>
  mutate(identity = str_replace(identity, "SAC_", ""),
         geomorph_class_num=str_replace(geomorph_class, "SAC-","") |> as.numeric(),
         geomorph_class = factor(geomorph_class, 
                                 levels = paste0("SAC-",seq(1,10,1)),
                                 labels = c("SAC-1"  = "Unconfined, boulder-bedrock, bed undulating",
                                            "SAC-2"  = "Confined, boulder, high gradient, step-pool/cascade",
                                            "SAC-3"  = "Confined, boulder-bedrock, uniform",
                                            "SAC-4"  = "Confined, boulder-bedrock, low-gradient step-pool",
                                            "SAC-5"  = "Confined, gravel-cobble, uniform",
                                            "SAC-6"  = "Partly-confined, low width-to-depth ratio, gravel-cobble, riffle-pool",
                                            "SAC-7"  = "Partly-confined, cobble-boulder, uniform",
                                            "SAC-8"  = "Partly-confined, high width-to-depth ratio, gravel-cobble, riffle-pool",
                                            "SAC-9"  = "Unconfined, low width-to-depth ratio, gravel",
                                            "SAC-10" = "Unconfined, gravel-cobble, riffle-pool"))) |>
  inner_join(read_csv(here::here("data-raw", "source", "ucd_geomorph", "geomorph_site_attributes.csv")), by=join_by(identity)) |>
  st_transform(project_crs) |> 
  st_join(flowlines |> select(comid), join=st_nearest_feature)
```

    ## Rows: 288 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): identity
    ## dbl (13): Ac, s, d, w, w/d, d/D50, CVd, CVw, k, D50, D84, Cv, Ls
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
geomorph_site_data |> 
  ggplot() + 
  geom_sf(aes(color=geomorph_class)) + 
  geom_sf(data=flowlines_sac_valley_major, color="darkgray")
```

![](geomorph_files/figure-gfm/geomorph-data-import-1.png)<!-- -->

``` r
geomorph_site_data |> select(comid, starts_with("geomorph_")) |>
  saveRDS(here::here("data-raw", "results", "geomorph_sites_ucd.Rds"))

flowlines_gcs |> 
  inner_join(geomorph_site_data |> st_drop_geometry()) |> 
  ggplot() + geom_sf(aes(color=geomorph_class)) + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

    ## Joining with `by = join_by(comid)`

![](geomorph_files/figure-gfm/geomorph-data-import-2.png)<!-- -->

``` r
geomorph_attr <- 
  geomorph_site_data |>
  left_join(flowline_attributes |> select(-starts_with("geomorph_")), by=join_by(comid)) |>
  mutate(geomorph_confined = str_split_i(geomorph_class, ", ", 1) |> as_factor(),
         geomorph_uniform = str_detect(geomorph_class, "uniform"),
         geomorph_steppool = str_detect(geomorph_class, "step-pool"), 
         geomorph_riffles = str_detect(geomorph_class, "riffle"),
         geomorph_gravel = str_detect(geomorph_class, "gravel"),
         geomorph_spawning = (geomorph_riffles | geomorph_steppool)) 

geomorph_attr |>
  ggplot() +
  geom_point(aes(x = slope*100, y = da_area_sq_km, color = geomorph_riffles)) + 
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  theme(panel.grid.minor = element_blank()) 
```

![](geomorph_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
geomorph_attr |>
  ggplot() +
  geom_point(aes(x = slope*100, y = da_area_sq_km, color = geomorph_class)) + 
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

![](geomorph_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

### version 1

### version 2

``` r
geomorph_training_data <- 
  geomorph_site_data |>
  left_join(flowline_attributes |> select(-starts_with("geomorph_")), by=join_by(comid), relationship="many-to-many") |>
  select(comid, geomorph_class, geomorph_class_num,
         slope, da_area_sq_km, da_elev_min, da_elev_max, da_elev_rel, da_ppt_mean_mm, mtpi30_min, hqt_gradient_class, hyd_cls) |>
  drop_na() |>
  filter(slope>1E-5) |>
  st_drop_geometry()

pca <- princomp(scale(geomorph_training_data |> select(-comid, -starts_with("geomorph_"), -hqt_gradient_class, -hyd_cls)))
factoextra::fviz_eig(pca, choice="variance", addlabels = TRUE)
```

![](geomorph_files/figure-gfm/geomorph-fill-gaps-train-1.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(1,2), col.var = "cos2", repel = TRUE)
```

![](geomorph_files/figure-gfm/geomorph-fill-gaps-train-2.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(1,3), col.var = "cos2", repel = TRUE)
```

![](geomorph_files/figure-gfm/geomorph-fill-gaps-train-3.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(1,4), col.var = "cos2", repel = TRUE)
```

![](geomorph_files/figure-gfm/geomorph-fill-gaps-train-4.png)<!-- -->

``` r
geomorph_rec <- 
  recipe(geomorph_class ~ slope + da_area_sq_km + da_elev_min + da_elev_max + da_elev_rel + da_ppt_mean_mm + mtpi30_min + hqt_gradient_class + hyd_cls,
         data=geomorph_training_data) |>
  #step_log(all_numeric_predictors()) |>#, -mtpi30_min) |>
  step_mutate_at(all_numeric_predictors(), fn = asinh) |>
  step_normalize(all_numeric_predictors()) |>
  step_pca(all_numeric_predictors(), num_comp = 4) |>
  step_dummy(hqt_gradient_class) |>
  step_dummy(hyd_cls) |>
  step_naomit()

geomorph_spec <- rand_forest(mode = "classification", trees = 256)

geomorph_fit <- 
  workflow() |>
  add_recipe(geomorph_rec) |>
  add_model(geomorph_spec) |>
  fit(data=geomorph_training_data)
```

``` r
geomorph_prediction_data <- flowlines_gcs |>
  left_join(flowline_attributes, by=join_by(comid)) |>
  filter(substr(reachcode,1, 4) %in% c("1802", "1803", "1804")) |>
  select(comid, 
         slope, da_area_sq_km, da_elev_min, da_elev_max, da_elev_rel, da_ppt_mean_mm, hqt_gradient_class, hyd_cls, mtpi30_min) |>
  drop_na() |>
 # filter(da_elev_min > 0) |>
  st_drop_geometry()
  
geomorph_pred <- 
  geomorph_prediction_data |>
  mutate(geomorph_class = predict(geomorph_fit, geomorph_prediction_data)[[".pred_class"]]) |>
  mutate(geomorph_confined = str_split_i(geomorph_class, ", ", 1) |> as_factor(),
         geomorph_uniform = str_detect(geomorph_class, "uniform"),
         geomorph_steppool = str_detect(geomorph_class, "step-pool"), 
         geomorph_riffles = str_detect(geomorph_class, "riffle"),
         geomorph_gravel = str_detect(geomorph_class, "gravel"),
         geomorph_spawning = (geomorph_riffles | geomorph_steppool))
  
flowlines_gcs |> 
  inner_join(geomorph_pred, by=join_by(comid)) |> 
  ggplot() + geom_sf(aes(color=geomorph_class)) + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

![](geomorph_files/figure-gfm/geomorph-fill-gaps-predict-1.png)<!-- -->

``` r
# validation matrix
geomorph_training_data |>
  select(comid, geomorph_class_actual=geomorph_class) |>
  inner_join(geomorph_pred |> select(comid, geomorph_class_pred=geomorph_class)) |>
  group_by(geomorph_class_actual, geomorph_class_pred) |>
  tally() |> 
  mutate(n=coalesce(n,0)) |>
  spread(geomorph_class_actual, n) |>
  knitr::kable()
```

    ## Joining with `by = join_by(comid)`

| geomorph_class_pred                                                    | Unconfined, boulder-bedrock, bed undulating | Confined, boulder, high gradient, step-pool/cascade | Confined, boulder-bedrock, uniform | Confined, boulder-bedrock, low-gradient step-pool | Confined, gravel-cobble, uniform | Partly-confined, low width-to-depth ratio, gravel-cobble, riffle-pool | Partly-confined, cobble-boulder, uniform | Partly-confined, high width-to-depth ratio, gravel-cobble, riffle-pool | Unconfined, low width-to-depth ratio, gravel | Unconfined, gravel-cobble, riffle-pool |
|:-----------------------------------------------------------------------|--------------------------------------------:|----------------------------------------------------:|-----------------------------------:|--------------------------------------------------:|---------------------------------:|----------------------------------------------------------------------:|-----------------------------------------:|-----------------------------------------------------------------------:|---------------------------------------------:|---------------------------------------:|
| Unconfined, boulder-bedrock, bed undulating                            |                                           3 |                                                  NA |                                 NA |                                                NA |                               NA |                                                                    NA |                                       NA |                                                                     NA |                                           NA |                                     NA |
| Confined, boulder, high gradient, step-pool/cascade                    |                                          NA |                                                  17 |                                  2 |                                                 3 |                                1 |                                                                    NA |                                       NA |                                                                     NA |                                            1 |                                      1 |
| Confined, boulder-bedrock, uniform                                     |                                          NA |                                                   3 |                                 23 |                                                 1 |                               NA |                                                                     1 |                                        5 |                                                                      1 |                                           NA |                                     NA |
| Confined, boulder-bedrock, low-gradient step-pool                      |                                          NA |                                                   2 |                                  3 |                                                27 |                               NA |                                                                    NA |                                        1 |                                                                      3 |                                           NA |                                     NA |
| Confined, gravel-cobble, uniform                                       |                                          NA |                                                  NA |                                  2 |                                                 1 |                               33 |                                                                     2 |                                       NA |                                                                     NA |                                            1 |                                     NA |
| Partly-confined, low width-to-depth ratio, gravel-cobble, riffle-pool  |                                           1 |                                                   2 |                                  4 |                                                 1 |                                3 |                                                                    39 |                                        3 |                                                                      6 |                                            2 |                                      2 |
| Partly-confined, cobble-boulder, uniform                               |                                          NA |                                                  NA |                                  2 |                                                NA |                                2 |                                                                    NA |                                       20 |                                                                      4 |                                            1 |                                     NA |
| Partly-confined, high width-to-depth ratio, gravel-cobble, riffle-pool |                                          NA |                                                  NA |                                 NA |                                                NA |                                2 |                                                                    NA |                                       NA |                                                                      8 |                                           NA |                                     NA |
| Unconfined, low width-to-depth ratio, gravel                           |                                          NA |                                                  NA |                                 NA |                                                NA |                               NA |                                                                     2 |                                       NA |                                                                     NA |                                           20 |                                     NA |
| Unconfined, gravel-cobble, riffle-pool                                 |                                          NA |                                                  NA |                                 NA |                                                NA |                               NA |                                                                    NA |                                       NA |                                                                     NA |                                           NA |                                     12 |

``` r
ggplot() +
  geom_point(data=geomorph_pred, aes(y = slope, x = da_area_sq_km, color = geomorph_class), size=0.25) + 
  geom_point(data=geomorph_attr, aes(y = slope, x = da_area_sq_km, fill = geomorph_class), shape=21, size=2) + 
  scale_y_log10(labels  = scales::label_percent()) + scale_x_log10(labels  = scales::label_comma()) + annotation_logticks() +
  theme(panel.grid.minor = element_blank())  + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color")) +
  theme(legend.position = "none") + ylab("Slope (%)") + xlab("Drainage Area (km2)")
```

    ## Warning in scale_x_log10(labels = scales::label_comma()): log-10 transformation
    ## introduced infinite values.

![](geomorph_files/figure-gfm/geomorph-fill-gaps-predict-2.png)<!-- -->

``` r
ggplot() +
  geom_point(data=geomorph_pred, aes(x = da_area_sq_km, y = da_elev_min/0.3048, color = geomorph_class), size=0.25) + 
  geom_point(data=geomorph_attr, aes(x = da_area_sq_km, y = da_elev_min/0.3048, fill = geomorph_class), shape=21, size=2) + 
  scale_y_sqrt(labels = scales::label_comma()) + scale_x_log10(labels  = scales::label_comma())  + annotation_logticks(sides="b") +
  theme(panel.grid.minor = element_blank())  + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color")) +
  theme(legend.position = "none") + xlab("Drainage Area (km2)") + ylab("Elevation (ft)")
```

    ## Warning in transformation$transform(x): NaNs produced

    ## Warning in scale_y_sqrt(labels = scales::label_comma()): sqrt transformation
    ## introduced infinite values.

    ## Warning in scale_x_log10(labels = scales::label_comma()): log-10 transformation
    ## introduced infinite values.

    ## Warning: Removed 423 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](geomorph_files/figure-gfm/geomorph-fill-gaps-predict-3.png)<!-- -->

``` r
ggplot() + 
  geom_sf(data = left_join(flowlines, geomorph_pred), aes(color = geomorph_class)) + 
  geom_sf(data = geomorph_site_data, aes(fill = geomorph_class), shape=21) + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

    ## Joining with `by = join_by(comid)`

    ## Warning: Removed 147795 rows containing missing values or values outside the scale range
    ## (`geom_sf()`).

![](geomorph_files/figure-gfm/geomorph-export2-1.png)<!-- -->

``` r
geomorph_pred |> select(comid, starts_with("geomorph_")) |>  #select(comid, geomorph_class) |> 
  saveRDS(here::here("data-raw", "results", "attr_geomorph_class.Rds"))
```

## UCD Geomorph Classes

Guillon, H., Byrne, C.F., Lane, B.A., Sandoval-Solis, S., and Pasternack
G.B. (2020). Machine Learning Predicts Reach-Scale Channel Types from
Coarse-Scale Geospatial Data in a Large River Basin. J. of Water
Resources Research. <https://doi.org/10.1029/2019WR026691> \[Article\]

Guillon, Hervé et al. (2020). Channel types predictions for the
Sacramento River basin \[Dataset\]. Dryad.
<https://doi.org/10.25338/B8031W>

Sandoval-Solis, S., Lane, B.A., Pasternack G.B., Byrne, C.F. &
Pasternack, G.B. Appendix F. Geomorphic Classification of California
Rivers. in California Environmental Flows Framework.
<https://ceff.ucdavis.edu/sites/g/files/dgvnsk5566/files/media/documents/Appendix_F%20Geomorphic%20Classification%20of%20CA.pdf>

``` r
geomorph_nhd_ucd <- read_sf(here::here("data-raw", "source", "ucd_geomorph", "SAC_channel-types_predictions_v1.shp.zip")) |>
  transmute(comid = COMID,
            geomorph_class_num=str_replace(group, "SAC","") |> as.numeric(),
            geomorph_class = factor(group, 
              levels = paste0("SAC",sprintf("%02d", seq(1,10,1))),
              labels = c("SAC01"  = "Unconfined, boulder-bedrock, bed undulating",
              "SAC02"  = "Confined, boulder, high gradient, step-pool/cascade",
              "SAC03"  = "Confined, boulder-bedrock, uniform",
              "SAC04"  = "Confined, boulder-bedrock, low-gradient step-pool",
              "SAC05"  = "Confined, gravel-cobble, uniform",
              "SAC06"  = "Partly-confined, low width-to-depth ratio, gravel-cobble, riffle-pool",
              "SAC07"  = "Partly-confined, cobble-boulder, uniform",
              "SAC08"  = "Partly-confined, high width-to-depth ratio, gravel-cobble, riffle-pool",
              "SAC09"  = "Unconfined, low width-to-depth ratio, gravel",
              "SAC010" = "Unconfined, gravel-cobble, riffle-pool")))

ggplot() + 
  geom_sf(data = geomorph_nhd_ucd, aes(color = geomorph_class)) + 
  geom_sf(data = geomorph_site_data, aes(fill = geomorph_class), shape=21) + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

![](geomorph_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Combine the official UCD dataset for Sac with our rougher classification
for San Joaquin

``` r
geomorph_final <- 
  geomorph_pred |>
  #geomorph_pred_v2 |>
  select(comid, geomorph_class_pred = geomorph_class) |>
  left_join(geomorph_nhd_ucd |>
              st_drop_geometry() |>
              select(comid, geomorph_class_ucd = geomorph_class), by=join_by(comid)) |>
  transmute(comid,
            geomorph_class = coalesce(geomorph_class_ucd, geomorph_class_pred)) |>
  mutate(geomorph_confined = str_split_i(geomorph_class, ", ", 1) |> as_factor(),
         geomorph_uniform = str_detect(geomorph_class, "uniform"),
         geomorph_steppool = str_detect(geomorph_class, "step-pool"), 
         geomorph_riffles = str_detect(geomorph_class, "riffle"),
         geomorph_gravel = str_detect(geomorph_class, "gravel"),
         geomorph_spawning = (geomorph_riffles | geomorph_steppool))

# # validation matrix
# geomorph_final |>
#   group_by(geomorph_class_ucd, geomorph_class_pred) |>
#   tally() |> 
#   mutate(n=coalesce(n,0)) |>
#   spread(geomorph_class_ucd, n) |>
#   knitr::kable()

ggplot() + 
  geom_sf(data = inner_join(flowlines, geomorph_final), aes(color = geomorph_class)) + 
  geom_sf(data = geomorph_site_data, aes(fill = geomorph_class), shape=21) + 
  scale_fill_brewer(type="qual", palette="Paired", aesthetics = c("fill", "color"))
```

    ## Joining with `by = join_by(comid)`

![](geomorph_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
geomorph_final |> select(comid, starts_with("geomorph_")) |>  #select(comid, geomorph_class) |> 
  saveRDS(here::here("data-raw", "results", "attr_geomorph_class.Rds"))
```
