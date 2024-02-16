Exploratory Modeling
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-02-16

- [Exploration and PCA of training
  data](#exploration-and-pca-of-training-data)
- [Set up training data](#set-up-training-data)
- [Split into training and testing
  datasets](#split-into-training-and-testing-datasets)
- [Modeling](#modeling)
  - [Linear regression](#linear-regression)
  - [Linear regression with regularized (lasso) feature
    selection](#linear-regression-with-regularized-lasso-feature-selection)
  - [Random Forest Regrssion](#random-forest-regrssion)
- [Random Forest model prediction
  output](#random-forest-model-prediction-output)
- [Predictions summarized by DSMHabitat
  reach](#predictions-summarized-by-dsmhabitat-reach)

``` r
library(tidyverse)
library(sf)
library(stars)
library(tidymodels)
library(broom.mixed) # tidy output of mixed model results
library(dotwhisker) # visualize regression results

knitr::opts_chunk$set(eval=TRUE)
```

``` r
flowlines <- readRDS("../data/flowline_geometries.Rds") |>
  st_transform("ESRI:102039")

flowline_attributes <- readRDS("../data/flowline_attributes.Rds")

train_data <- flowlines |> st_drop_geometry() |>
  left_join(readRDS("../data/flowline_attributes.Rds"), by=join_by("comid"), relationship="one-to-one") |>
  inner_join(readRDS("../data/fsa_combined.Rds"), by=join_by("comid"), relationship="one-to-many") |> 
  glimpse()
```

    ## Rows: 738
    ## Columns: 120

    ## Warning in grepl(",", levels(x), fixed = TRUE): input string 1 is invalid in
    ## this locale

    ## Warning in grepl(",", levels(x), fixed = TRUE): input string 2 is invalid in
    ## this locale

    ## Warning in grepl(",", levels(x), fixed = TRUE): input string 3 is invalid in
    ## this locale

    ## Warning in grepl(",", levels(x), fixed = TRUE): input string 4 is invalid in
    ## this locale

    ## Warning in grepl(",", levels(x), fixed = TRUE): input string 6 is invalid in
    ## this locale

    ## $ comid                        <dbl> 12071480, 12071480, 12071480, 12071480, 1…
    ## $ reachcode                    <fct> 18020157005039, 18020157005039, 180201570…
    ## $ gnis_id                      <fct> 1655075, 1655075, 1655075, 1655075, 16550…
    ## $ gnis_name                    <fct> "Deer Creek", "Deer Creek", "Deer Creek",…
    ## $ lengthkm                     <dbl> 3.414, 3.414, 3.414, 3.414, 3.414, 3.414,…
    ## $ ftype                        <fct> StreamRiver, StreamRiver, StreamRiver, St…
    ## $ fcode                        <int> 46006, 46006, 46006, 46006, 46006, 46006,…
    ## $ huc_8                        <chr> "18020157", "18020157", "18020157", "1802…
    ## $ huc_10                       <chr> "1802015700", "1802015700", "1802015700",…
    ## $ huc_12                       <chr> "180201570050", "180201570050", "18020157…
    ## $ ftype_desc                   <chr> "Stream/River", "Stream/River", "Stream/R…
    ## $ hydro_seq                    <dbl> 10012820, 10012820, 10012820, 10012820, 1…
    ## $ reach_code                   <fct> 18020157005039, 18020157005039, 180201570…
    ## $ stream_level                 <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ stream_order                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
    ## $ us_length_km                 <dbl> 389.906, 389.906, 389.906, 389.906, 389.9…
    ## $ ds_length_km                 <dbl> 346.603, 346.603, 346.603, 346.603, 346.6…
    ## $ da_area_sq_km                <dbl> 573.3900, 573.3900, 573.3900, 573.3900, 5…
    ## $ reach_length_km              <dbl> 3.414, 3.414, 3.414, 3.414, 3.414, 3.414,…
    ## $ slope                        <dbl> 0.00272993, 0.00272993, 0.00272993, 0.002…
    ## $ elev_min                     <dbl> 51.80, 51.80, 51.80, 51.80, 51.80, 51.80,…
    ## $ elev_max                     <dbl> 61.12, 61.12, 61.12, 61.12, 61.12, 61.12,…
    ## $ stream_power                 <dbl> 1.565314563, 1.565314563, 1.565314563, 1.…
    ## $ da_ppt_mean_mm               <dbl> 1437.351, 1437.351, 1437.351, 1437.351, 1…
    ## $ loc_ppt_mean_mm              <dbl> 584.914, 584.914, 584.914, 584.914, 584.9…
    ## $ vogel_q_ma_cfs               <dbl> 593.8086, 593.8086, 593.8086, 593.8086, 5…
    ## $ vogel_v_ma_fps               <dbl> 1.63041, 1.63041, 1.63041, 1.63041, 1.630…
    ## $ erom_q_ma_cfs                <dbl> 347.468, 347.468, 347.468, 347.468, 347.4…
    ## $ erom_v_ma_fps                <dbl> 1.57425, 1.57425, 1.57425, 1.57425, 1.574…
    ## $ sinuosity                    <dbl> 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.27,…
    ## $ da_avg_slope                 <dbl> 25.97, 25.97, 25.97, 25.97, 25.97, 25.97,…
    ## $ da_elev_mean                 <dbl> 1206.24, 1206.24, 1206.24, 1206.24, 1206.…
    ## $ da_elev_min                  <dbl> 51.80, 51.80, 51.80, 51.80, 51.80, 51.80,…
    ## $ da_elev_max                  <dbl> 2393.52, 2393.52, 2393.52, 2393.52, 2393.…
    ## $ da_elev_rel                  <dbl> 2341.72, 2341.72, 2341.72, 2341.72, 2341.…
    ## $ loc_bfi                      <dbl> 42.5746, 42.5746, 42.5746, 42.5746, 42.57…
    ## $ da_bfi                       <dbl> 63.3563, 63.3563, 63.3563, 63.3563, 63.35…
    ## $ loc_twi                      <dbl> 977.304, 977.304, 977.304, 977.304, 977.3…
    ## $ da_twi                       <dbl> 676.9103, 676.9103, 676.9103, 676.9103, 6…
    ## $ loc_pct_clay                 <dbl> 13.3988, 13.3988, 13.3988, 13.3988, 13.39…
    ## $ da_pct_clay                  <dbl> 21.6023, 21.6023, 21.6023, 21.6023, 21.60…
    ## $ loc_pct_sand                 <dbl> 48.7830, 48.7830, 48.7830, 48.7830, 48.78…
    ## $ da_pct_sand                  <dbl> 38.4830, 38.4830, 38.4830, 38.4830, 38.48…
    ## $ loc_permeability             <dbl> 7.306, 7.306, 7.306, 7.306, 7.306, 7.306,…
    ## $ loc_bedrock_depth            <dbl> 151.2547, 151.2547, 151.2547, 151.2547, 1…
    ## $ da_permeability              <dbl> 6.1970, 6.1970, 6.1970, 6.1970, 6.1970, 6…
    ## $ da_bedrock_depth             <dbl> 105.5002, 105.5002, 105.5002, 105.5002, 1…
    ## $ loc_k_erodibility            <dbl> 0.2617, 0.2617, 0.2617, 0.2617, 0.2617, 0…
    ## $ da_k_erodibility             <dbl> 0.2964, 0.2964, 0.2964, 0.2964, 0.2964, 0…
    ## $ loc_precip                   <dbl> 582.4351, 582.4351, 582.4351, 582.4351, 5…
    ## $ da_precip                    <dbl> 1417.543, 1417.543, 1417.543, 1417.543, 1…
    ## $ loc_runoff                   <dbl> 521, 521, 521, 521, 521, 521, 521, 521, 5…
    ## $ da_runoff                    <dbl> 822.3441, 822.3441, 822.3441, 822.3441, 8…
    ## $ bf_width_m                   <dbl> 31.47, 31.47, 31.47, 31.47, 31.47, 31.47,…
    ## $ bf_depth_m                   <dbl> 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72,…
    ## $ bf_w_d_ratio                 <dbl> 18.29651, 18.29651, 18.29651, 18.29651, 1…
    ## $ mean_ndvi                    <dbl> 0.4869165, 0.4869165, 0.4869165, 0.486916…
    ## $ bio_aq_rank_sw               <int> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,…
    ## $ species                      <dbl> 199, 199, 199, 199, 199, 199, 199, 199, 1…
    ## $ species_fish                 <dbl> 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 1…
    ## $ species_crust                <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
    ## $ species_herps                <dbl> 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,…
    ## $ species_inverts              <dbl> 116, 116, 116, 116, 116, 116, 116, 116, 1…
    ## $ species_mollusks             <dbl> 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,…
    ## $ species_plants               <dbl> 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 2…
    ## $ species_birds                <dbl> 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 1…
    ## $ species_mammals              <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
    ## $ species_mollusks_crust       <dbl> 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1…
    ## $ species_endemic              <dbl> 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 4…
    ## $ species_endemic_fish         <dbl> 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 1…
    ## $ species_endemic_crust        <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ species_endemic_herps        <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,…
    ## $ species_endemic_inverts      <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 2…
    ## $ species_endemic_mollusks     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ species_endemic_plants       <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,…
    ## $ species_endemic_birds        <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ species_endemic_mammals      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ species_vulnerable           <dbl> 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 5…
    ## $ species_listed               <dbl> 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 1…
    ## $ species_endemic_vulnerable   <dbl> 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 2…
    ## $ species_endemic_listed       <dbl> 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,…
    ## $ genus                        <dbl> 160, 160, 160, 160, 160, 160, 160, 160, 1…
    ## $ family                       <dbl> 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 8…
    ## $ tax_order                    <dbl> 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 4…
    ## $ tax_class                    <dbl> 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 1…
    ## $ phylum                       <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
    ## $ species_current              <dbl> 165, 165, 165, 165, 165, 165, 165, 165, 1…
    ## $ species_historical           <dbl> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
    ## $ species_other                <dbl> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 4…
    ## $ hyd_cat                      <fct> Mixed, Mixed, Mixed, Mixed, Mixed, Mixed,…
    ## $ hyd_cls                      <fct> "Low-volume snowmelt and rain", "Low-volu…
    ## $ nf_bfl_dry_cfs               <dbl> 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.3,…
    ## $ nf_bfl_wet_cfs               <dbl> 178, 178, 178, 178, 178, 178, 178, 178, 1…
    ## $ nf_wet_start                 <dbl> 85.2, 85.2, 85.2, 85.2, 85.2, 85.2, 85.2,…
    ## $ nf_wet_dur                   <dbl> 136, 136, 136, 136, 136, 136, 136, 136, 1…
    ## $ peak_q2_cfs                  <dbl> 3870, 3870, 3870, 3870, 3870, 3870, 3870,…
    ## $ peak_q5_cfs                  <dbl> 7430, 7430, 7430, 7430, 7430, 7430, 7430,…
    ## $ peak_q10_cfs                 <dbl> 8720, 8720, 8720, 8720, 8720, 8720, 8720,…
    ## $ merit_width_m                <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, N…
    ## $ vb_width_transect            <dbl> 31.4700, 31.4700, 31.4700, 31.4700, 31.47…
    ## $ frac_leveed_longitudinal     <dbl> 0.4390244, 0.4390244, 0.4390244, 0.439024…
    ## $ lateral_levee_confinement_ft <dbl> 1872.016, 1872.016, 1872.016, 1872.016, 1…
    ## $ bf_xarea_m                   <dbl> 54.1284, 54.1284, 54.1284, 54.1284, 54.12…
    ## $ chan_width_m                 <dbl> 31.47, 31.47, 31.47, 31.47, 31.47, 31.47,…
    ## $ velocity_m_s                 <dbl> 5.164862, 5.164862, 5.164862, 5.164862, 5…
    ## $ wetted_perimeter_m           <dbl> 34.91, 34.91, 34.91, 34.91, 34.91, 34.91,…
    ## $ hydraulic_radius_m           <dbl> 1.550513, 1.550513, 1.550513, 1.550513, 1…
    ## $ critical_shields_number      <dbl> 0.03428697, 0.03428697, 0.03428697, 0.034…
    ## $ grain_size_mobilized_mm      <dbl> 0.74819322, 0.74819322, 0.74819322, 0.748…
    ## $ shear_velocity_cm_s          <dbl> 20.377361, 20.377361, 20.377361, 20.37736…
    ## $ settling_velocity_ndim       <dbl> 522.746144, 522.746144, 522.746144, 522.7…
    ## $ grain_size_suspended_ndim    <dbl> 1746.03995, 1746.03995, 1746.03995, 1746.…
    ## $ grain_size_suspended_mm      <dbl> 0.148708751, 0.148708751, 0.148708751, 0.…
    ## $ mtpi30_min                   <dbl> -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -…
    ## $ vb_bf_w_ratio                <dbl> 1.00000, 1.00000, 1.00000, 1.00000, 1.000…
    ## $ dataset                      <chr> "Deer Creek", "Deer Creek", "Deer Creek",…
    ## $ flow_cfs                     <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000,…
    ## $ area_tot_ft2                 <dbl> 165174, 195682, 218107, 232765, 243522, 2…
    ## $ area_wua_ft2                 <dbl> 126770, 155882, 166777, 171971, 174341, 1…
    ## $ hsi_frac                     <dbl> 0.7674937, 0.7966088, 0.7646568, 0.738818…

## Exploration and PCA of training data

``` r
df <- flowline_attributes |>
  select(
         # predictors of flow (as would be found in a regional regression)
         slope, da_area_sq_km, da_elev_mean, da_ppt_mean_mm, 
         # baseflow and peakflow statistics
         nf_bfl_dry_cfs, nf_bfl_wet_cfs, erom_q_ma_cfs, peak_q2_cfs, peak_q5_cfs, peak_q10_cfs,
         # flow and channel characteristics, hydrologically predicted
         bf_depth_m, bf_w_d_ratio, erom_v_ma_fps,
         # misc characteristics of the catchment
         da_avg_slope, da_k_erodibility, mean_ndvi,
         # misc characteristics of the locality
         loc_bfi, loc_pct_clay, loc_pct_sand, loc_permeability, loc_bedrock_depth, loc_ppt_mean_mm,
         # channel confinement characteristics
         mtpi30_min, vb_width_transect, vb_bf_w_ratio, 
         # channel sinuosity
         sinuosity,
         # aquatic
         bio_aq_rank_sw, species, species_fish
         ) |> 
  na.omit() |> 
  scale() 

df |> cor() |> corrr::rplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](model-expl_files/figure-gfm/pca-1.png)<!-- -->

``` r
# uniqueness: 
# df |> cor() |> abs() |> rowMeans() |> rank() 
pca <- princomp(df) 
#summary(pca)

factoextra::fviz_eig(pca, addlabels = TRUE)
```

![](model-expl_files/figure-gfm/pca-2.png)<!-- -->

``` r
factoextra::fviz_cos2(pca, choice = "var", axes = 1:2)
```

![](model-expl_files/figure-gfm/pca-3.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(1,2), col.var = "cos2", repel = TRUE)
```

![](model-expl_files/figure-gfm/pca-4.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(1,3), col.var = "cos2", repel = TRUE)
```

![](model-expl_files/figure-gfm/pca-5.png)<!-- -->

``` r
factoextra::fviz_pca_var(pca, axes=c(2,3), col.var = "cos2", repel = TRUE)
```

![](model-expl_files/figure-gfm/pca-6.png)<!-- -->

``` r
# flowlines |> GGally::ggpairs()
```

## Set up training data

``` r
td <- train_data |>
  mutate(flow_norm_cfs = flow_cfs / erom_q_ma_cfs) |> # flow as a percent of mean annual flow
  mutate(wua_per_lf = area_wua_ft2 / (reach_length_km*3280.84),
         log_wua_per_lf = log(wua_per_lf)) |> # transect-wise habitat area per linear foot
  select(dataset, comid, 
         # suitable habitat area normalized by reach length
         hsi_frac, wua_per_lf, log_wua_per_lf,
         # flow cfs normalized by mean annual flow
         flow_cfs, flow_norm_cfs,
         # predictors of flow (as would be found in a regional regression)
         slope, da_area_sq_km, da_elev_mean, da_ppt_mean_mm, 
         # baseflow and peakflow statistics
         nf_bfl_dry_cfs, nf_bfl_wet_cfs, erom_q_ma_cfs, peak_q2_cfs, peak_q5_cfs, peak_q10_cfs,
         # flow and channel characteristics, hydrologically predicted
         bf_depth_m, bf_w_d_ratio, erom_v_ma_fps,
         # misc characteristics of the catchment
         da_avg_slope, da_k_erodibility, mean_ndvi,
         # misc characteristics of the locality
         loc_bfi, loc_pct_clay, loc_pct_sand, loc_permeability, loc_bedrock_depth, loc_ppt_mean_mm,
         # channel confinement characteristics
         mtpi30_min, vb_width_transect, vb_bf_w_ratio, frac_leveed_longitudinal,
         # channel sinuosity
         sinuosity,
         # aquatic
         bio_aq_rank_sw, species, species_fish
         ) |> 
  mutate(nf_bfl_dry_cfs_norm = nf_bfl_dry_cfs/erom_q_ma_cfs, 
         nf_bfl_wet_cfs_norm = nf_bfl_wet_cfs/erom_q_ma_cfs) |>
  drop_na() |> glimpse()
```

    ## Rows: 738
    ## Columns: 39
    ## $ dataset                  <chr> "Deer Creek", "Deer Creek", "Deer Creek", "De…
    ## $ comid                    <dbl> 12071480, 12071480, 12071480, 12071480, 12071…
    ## $ hsi_frac                 <dbl> 0.7674937, 0.7966088, 0.7646568, 0.7388181, 0…
    ## $ wua_per_lf               <dbl> 11.31795, 13.91706, 14.88976, 15.35347, 15.56…
    ## $ log_wua_per_lf           <dbl> 2.426390, 2.633115, 2.700673, 2.731342, 2.745…
    ## $ flow_cfs                 <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 500…
    ## $ flow_norm_cfs            <dbl> 0.2877963, 0.7194907, 0.8633889, 1.1511851, 1…
    ## $ slope                    <dbl> 0.00272993, 0.00272993, 0.00272993, 0.0027299…
    ## $ da_area_sq_km            <dbl> 573.3900, 573.3900, 573.3900, 573.3900, 573.3…
    ## $ da_elev_mean             <dbl> 1206.24, 1206.24, 1206.24, 1206.24, 1206.24, …
    ## $ da_ppt_mean_mm           <dbl> 1437.351, 1437.351, 1437.351, 1437.351, 1437.…
    ## $ nf_bfl_dry_cfs           <dbl> 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.…
    ## $ nf_bfl_wet_cfs           <dbl> 178, 178, 178, 178, 178, 178, 178, 178, 178, …
    ## $ erom_q_ma_cfs            <dbl> 347.468, 347.468, 347.468, 347.468, 347.468, …
    ## $ peak_q2_cfs              <dbl> 3870, 3870, 3870, 3870, 3870, 3870, 3870, 387…
    ## $ peak_q5_cfs              <dbl> 7430, 7430, 7430, 7430, 7430, 7430, 7430, 743…
    ## $ peak_q10_cfs             <dbl> 8720, 8720, 8720, 8720, 8720, 8720, 8720, 872…
    ## $ bf_depth_m               <dbl> 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.7…
    ## $ bf_w_d_ratio             <dbl> 18.29651, 18.29651, 18.29651, 18.29651, 18.29…
    ## $ erom_v_ma_fps            <dbl> 1.57425, 1.57425, 1.57425, 1.57425, 1.57425, …
    ## $ da_avg_slope             <dbl> 25.97, 25.97, 25.97, 25.97, 25.97, 25.97, 25.…
    ## $ da_k_erodibility         <dbl> 0.2964, 0.2964, 0.2964, 0.2964, 0.2964, 0.296…
    ## $ mean_ndvi                <dbl> 0.4869165, 0.4869165, 0.4869165, 0.4869165, 0…
    ## $ loc_bfi                  <dbl> 42.5746, 42.5746, 42.5746, 42.5746, 42.5746, …
    ## $ loc_pct_clay             <dbl> 13.3988, 13.3988, 13.3988, 13.3988, 13.3988, …
    ## $ loc_pct_sand             <dbl> 48.7830, 48.7830, 48.7830, 48.7830, 48.7830, …
    ## $ loc_permeability         <dbl> 7.306, 7.306, 7.306, 7.306, 7.306, 7.306, 7.3…
    ## $ loc_bedrock_depth        <dbl> 151.2547, 151.2547, 151.2547, 151.2547, 151.2…
    ## $ loc_ppt_mean_mm          <dbl> 584.914, 584.914, 584.914, 584.914, 584.914, …
    ## $ mtpi30_min               <dbl> -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -…
    ## $ vb_width_transect        <dbl> 31.4700, 31.4700, 31.4700, 31.4700, 31.4700, …
    ## $ vb_bf_w_ratio            <dbl> 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, …
    ## $ frac_leveed_longitudinal <dbl> 0.4390244, 0.4390244, 0.4390244, 0.4390244, 0…
    ## $ sinuosity                <dbl> 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.2…
    ## $ bio_aq_rank_sw           <int> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
    ## $ species                  <dbl> 199, 199, 199, 199, 199, 199, 199, 199, 199, …
    ## $ species_fish             <dbl> 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 1…
    ## $ nf_bfl_dry_cfs_norm      <dbl> 0.2138326, 0.2138326, 0.2138326, 0.2138326, 0…
    ## $ nf_bfl_wet_cfs_norm      <dbl> 0.5122774, 0.5122774, 0.5122774, 0.5122774, 0…

``` r
td |> 
  mutate(across(slope:last_col(), function(x) factor(cut(x, breaks=5, labels=F)))) |>
  pivot_longer(cols=slope:last_col(), names_to="varname", values_to="qtile") |>
  group_by(flow_cfs, varname, qtile) |> summarize(wua_per_lf=mean(wua_per_lf)) |> ungroup() |>
  ggplot() + geom_line(aes(x=flow_cfs, y=wua_per_lf, color=qtile)) + facet_wrap(~varname) + scale_color_viridis_d() + scale_y_log10() + scale_x_log10()
```

    ## `summarise()` has grouped output by 'flow_cfs', 'varname'. You can override
    ## using the `.groups` argument.

![](model-expl_files/figure-gfm/td-gridplot-1.png)<!-- -->

## Split into training and testing datasets

``` r
set.seed(47)
td_split <- initial_split(td, strata=dataset)
```

## Modeling

### Linear regression

``` r
lm_spec <- linear_reg() 
```

First set up recipes (sets of predictor and response variables), then
train and predict

#### Scale-independent (%HSI versus dimensionless flow)

Scale-independent model version predicts the percent suitable habitat
area using normalized flow, and without using any direct correlates of
watershed size.

``` r
si_rec <- recipe(data=training(td_split), 
      formula=hsi_frac ~ 
        # flow as percent of mean annual flow
        flow_norm_cfs + 
        # channel characteristics: gradient and sinuosity
        slope + sinuosity + 
        # flow and channel characteristics, hydrologically predicted
        erom_v_ma_fps + bf_depth_m + bf_w_d_ratio + 
        # misc characteristics of the catchment
        da_k_erodibility + da_avg_slope + mean_ndvi +
        # misc characteristics of the locality
        loc_bfi + loc_pct_clay + loc_pct_sand + loc_permeability + loc_bedrock_depth + loc_ppt_mean_mm +
        # channel confinement
        mtpi30_min + vb_bf_w_ratio + frac_leveed_longitudinal +
        # baseflow as percent of mean annual flow
        nf_bfl_dry_cfs_norm + nf_bfl_wet_cfs_norm + 
        # aquatic
        bio_aq_rank_sw + species + species_fish
        ) |>
  step_log(all_numeric_predictors(), -mtpi30_min, -mean_ndvi, -bio_aq_rank_sw, -vb_bf_w_ratio, -frac_leveed_longitudinal) |>
  # let effect of flow vary with gradient
  step_interact(terms = ~ flow_norm_cfs:slope)
  # try interacting with all numeric predictors
  #step_interact(terms = ~ flow_norm_cfs:all_numeric_predictors()) 

lm_si <- workflow() |>
  add_recipe(si_rec |> step_normalize(all_numeric_predictors())) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_si |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic  p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.636         0.620 0.129      38.5 1.05e-99    24   359. -667. -555.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_si |> tidy()
```

    ## # A tibble: 25 × 5
    ##    term             estimate std.error statistic   p.value
    ##    <chr>               <dbl>     <dbl>     <dbl>     <dbl>
    ##  1 (Intercept)       0.439     0.00550    79.8   7.66e-297
    ##  2 flow_norm_cfs    -0.170     0.0440     -3.87  1.24e-  4
    ##  3 slope             0.105     0.114       0.917 3.60e-  1
    ##  4 sinuosity        -0.0169    0.00760    -2.23  2.64e-  2
    ##  5 erom_v_ma_fps    -0.133     0.123      -1.08  2.81e-  1
    ##  6 bf_depth_m        0.0283    0.231       0.123 9.03e-  1
    ##  7 bf_w_d_ratio      0.00731   0.0463      0.158 8.75e-  1
    ##  8 da_k_erodibility  0.427     0.246       1.73  8.34e-  2
    ##  9 da_avg_slope      0.0491    0.0238      2.07  3.92e-  2
    ## 10 mean_ndvi         0.0110    0.0116      0.949 3.43e-  1
    ## # ℹ 15 more rows

``` r
lm_si$fit$fit |> dotwhisker::dwplot()
```

![](model-expl_files/figure-gfm/lm-si-1.png)<!-- -->

``` r
lm_si_res <- testing(td_split) |>
  mutate(hsi_frac_pred = predict(lm_si, testing(td_split))[[".pred"]]) |>
  select(comid, hsi_frac, hsi_frac_pred, flow_norm_cfs) 

lm_si_res |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-independent model, linear regression")
```

![](model-expl_files/figure-gfm/lm-si-2.png)<!-- -->

``` r
lm_si_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs)) |> #, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*100, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*100, color="actual", group=1)) + 
  xlab("Flow (cfs)") + ylab("Suitable Habitat Area (% of Total Area)")
```

![](model-expl_files/figure-gfm/lm-si-3.png)<!-- -->

#### Scale-dependent (WUA-per-linear-ft versus flow)

Scale-dependent model version predicts the WUA per linear foot, using
absolute flow values as well as watershed scale parameters.

- interact flow and mean annual flow; flow and area/elev/ppt; to add
  scale dependence to the flow effect

``` r
sd_rec <- recipe(data=training(td_split), 
      formula = log_wua_per_lf ~ 
        # flow (cfs)
        flow_cfs + 
        # predictors of flow (catchment area, elevation, and MAP) -- attributes for gradient and upstream drainage area are interrelated
        slope + da_area_sq_km + da_elev_mean + da_ppt_mean_mm +
        # baseflow/peak flow statistics
        nf_bfl_dry_cfs + nf_bfl_wet_cfs + erom_q_ma_cfs + #log(peak_q2_cfs)  + log(peak_q5_cfs) + log(peak_q10_cfs) +
        # flow and channel characteristics, hydrologically predicted
        erom_v_ma_fps + bf_depth_m + bf_w_d_ratio + 
        # misc characteristics of the catchment
        da_k_erodibility + da_avg_slope + mean_ndvi +
        # misc characteristics of the locality
        loc_bfi + loc_pct_clay + loc_pct_sand + loc_permeability + loc_bedrock_depth + loc_ppt_mean_mm +
        # channel confinement characteristics
        mtpi30_min + vb_bf_w_ratio + frac_leveed_longitudinal + vb_width_transect +
        # channel sinuosity
        sinuosity + 
        # aquatic
        bio_aq_rank_sw + species + species_fish
        ) |>
  step_log(all_numeric_predictors(), -mtpi30_min, -mean_ndvi, -vb_bf_w_ratio, -frac_leveed_longitudinal) |>
  step_interact(terms = ~ flow_cfs:erom_q_ma_cfs) |>
  step_interact(terms = ~ slope:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_elev_mean) |>
  step_interact(terms = ~ flow_cfs:da_ppt_mean_mm) 

lm_sd <- workflow() |>
  add_recipe(sd_rec |> step_normalize(all_numeric_predictors())) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_sd |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.965         0.962 0.450      428.       0    33  -325.  720.  871.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_sd |> tidy()
```

    ## # A tibble: 34 × 5
    ##    term           estimate std.error statistic   p.value
    ##    <chr>             <dbl>     <dbl>     <dbl>     <dbl>
    ##  1 (Intercept)      0.833     0.0191   43.6    1.33e-175
    ##  2 flow_cfs       -22.9       9.02     -2.54   1.13e-  2
    ##  3 slope            2.39      0.896     2.66   7.99e-  3
    ##  4 da_area_sq_km  -25.9       6.89     -3.76   1.88e-  4
    ##  5 da_elev_mean    41.9      10.2       4.11   4.53e-  5
    ##  6 da_ppt_mean_mm  55.6       7.60      7.32   9.54e- 13
    ##  7 nf_bfl_dry_cfs  33.6       9.68      3.48   5.52e-  4
    ##  8 nf_bfl_wet_cfs  27.7       4.58      6.06   2.65e-  9
    ##  9 erom_q_ma_cfs   20.9       5.31      3.92   9.88e-  5
    ## 10 erom_v_ma_fps    0.0474    0.859     0.0552 9.56e-  1
    ## # ℹ 24 more rows

``` r
lm_sd$fit$fit |> dotwhisker::dwplot()
```

![](model-expl_files/figure-gfm/lm-sd-1.png)<!-- -->

``` r
lm_sd_res <-
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(lm_sd, testing(td_split))[[".pred"]]) |>
  transmute(comid, wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) 

lm_sd_res |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-dependent model, linear regression")
```

![](model-expl_files/figure-gfm/lm-sd-2.png)<!-- -->

``` r
lm_sd_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() + xlab("Flow (cfs)") + ylab("Suitable Habitat Area per LF Channel (ft)")
```

![](model-expl_files/figure-gfm/lm-sd-3.png)<!-- -->

### Linear regression with regularized (lasso) feature selection

``` r
#install.packages("glmnet")
lasso_spec <- linear_reg(penalty = 0.01, mixture = 1) |> set_engine("glmnet")
```

#### Scale-independent (%HSI versus dimensionless flow)

``` r
drop_vars <- workflow() |>
  add_recipe(si_rec |> step_normalize(all_numeric_predictors())) |>
  add_model(lasso_spec) |>
  fit(data=training(td_split)) |>
  tidy() |> 
  filter(estimate==0) |> 
  pull(term)
```

    ## Warning: package 'glmnet' was built under R version 4.3.2

    ## Warning: package 'Matrix' was built under R version 4.3.2

``` r
lasso_si <- workflow() |>
  add_recipe(si_rec |> step_rm(all_of(drop_vars))) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lasso_si |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic   p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>     <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.597         0.592 0.134      115. 2.65e-103     7   331. -644. -605.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_si |> tidy()
```

    ## # A tibble: 8 × 5
    ##   term                  estimate std.error statistic   p.value
    ##   <chr>                    <dbl>     <dbl>     <dbl>     <dbl>
    ## 1 (Intercept)             0.581   0.949        0.612 5.41e-  1
    ## 2 sinuosity              -0.135   0.0300      -4.51  7.80e-  6
    ## 3 erom_v_ma_fps          -0.134   0.0342      -3.91  1.04e-  4
    ## 4 bf_depth_m             -0.258   0.0252     -10.2   1.27e- 22
    ## 5 bf_w_d_ratio           -0.288   0.113       -2.55  1.11e-  2
    ## 6 da_avg_slope            0.223   0.339        0.658 5.11e-  1
    ## 7 bio_aq_rank_sw          0.0940  0.0185       5.07  5.37e-  7
    ## 8 flow_norm_cfs_x_slope   0.0199  0.000753    26.5   5.88e-100

``` r
lasso_si$fit$fit |> dotwhisker::dwplot()
```

![](model-expl_files/figure-gfm/lasso-si-1.png)<!-- -->

``` r
lasso_si_res <- 
testing(td_split) |>
  mutate(hsi_frac_pred = predict(lasso_si, testing(td_split))[[".pred"]]) |>
  select(comid, hsi_frac, hsi_frac_pred, flow_norm_cfs)

lasso_si_res |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c()  + 
  ggtitle("Scale-independent model, linear regression w/ lasso feature selection")
```

![](model-expl_files/figure-gfm/lasso-si-2.png)<!-- -->

``` r
lasso_si_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs)) |> #, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*100, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*100, color="actual", group=1)) + 
  xlab("Flow (cfs)") + ylab("Suitable Habitat Area (% of Total Area)")
```

![](model-expl_files/figure-gfm/lasso-si-3.png)<!-- -->

#### Scale-dependent (WUA-per-linear-ft versus flow)

``` r
drop_vars <- workflow() |>
  add_recipe(sd_rec |> step_normalize(all_numeric_predictors())) |>
  add_model(lasso_spec) |>
  fit(data=training(td_split)) |>
  tidy() |> 
  filter(estimate==0) |> 
  pull(term)

lasso_sd <- workflow() |>
  add_recipe(sd_rec |>
               step_rm(all_of(drop_vars))) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lasso_sd |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic   p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>     <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.944         0.943 0.554      569. 4.94e-324    16  -449.  934. 1012.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_sd |> tidy()
```

    ## # A tibble: 17 × 5
    ##    term                      estimate std.error statistic  p.value
    ##    <chr>                        <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)               -19.8     18.4        -1.08  2.81e- 1
    ##  2 da_elev_mean               -3.95     1.55       -2.55  1.11e- 2
    ##  3 da_ppt_mean_mm              9.94     1.29        7.69  7.06e-14
    ##  4 erom_v_ma_fps               0.825    0.193       4.26  2.36e- 5
    ##  5 bf_depth_m                  0.618    0.133       4.66  4.00e- 6
    ##  6 bf_w_d_ratio               -4.90     0.696      -7.04  5.97e-12
    ##  7 mean_ndvi                  -0.147    0.608      -0.242 8.09e- 1
    ##  8 loc_bfi                     1.15     0.474       2.42  1.60e- 2
    ##  9 loc_pct_sand               -0.332    0.219      -1.52  1.30e- 1
    ## 10 loc_bedrock_depth          -0.380    0.334      -1.14  2.55e- 1
    ## 11 mtpi30_min                  0.0520   0.00909     5.72  1.80e- 8
    ## 12 frac_leveed_longitudinal   -0.334    0.161      -2.08  3.81e- 2
    ## 13 vb_width_transect           0.0776   0.0378      2.05  4.05e- 2
    ## 14 sinuosity                  -0.621    0.138      -4.49  8.76e- 6
    ## 15 bio_aq_rank_sw             -0.0874   0.814      -0.107 9.15e- 1
    ## 16 species_fish               -3.34     0.676      -4.94  1.04e- 6
    ## 17 flow_cfs_x_da_ppt_mean_mm   0.0266   0.00298     8.92  7.07e-18

``` r
lasso_sd$fit$fit |> dotwhisker::dwplot()
```

![](model-expl_files/figure-gfm/lasso-sd-1.png)<!-- -->

``` r
lasso_sd_res <- 
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(lasso_sd, testing(td_split))[[".pred"]]) |>
  transmute(comid, wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) 

lasso_sd_res |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c()   + 
  ggtitle("Scale-dependent model, linear regression w/ lasso feature selection")
```

![](model-expl_files/figure-gfm/lasso-sd-2.png)<!-- -->

``` r
lasso_sd_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() + xlab("Flow (cfs)") + ylab("Suitable Habitat Area per LF Channel (ft)")
```

![](model-expl_files/figure-gfm/lasso-sd-3.png)<!-- -->

### Random Forest Regrssion

``` r
#install.packages("ranger")
#rfr_spec <- rand_forest(mode = "regression", trees = 2^12)
```

#### Scale-independent (%HSI versus dimensionless flow)

``` r
# set up cross validation
folds <- vfold_cv(training(td_split), v = 10, strata = dataset)
val_metrics <- metric_set(rmse, rsq, ccc)

# tuned random forest: how many trees to use?
rfr_spec <- rand_forest(mode = "regression", trees = tune())
rfr_grid <- tibble(trees = 2^seq(0, 12, 1))

rfr_si_tune <- 
  workflow() |>
  add_recipe(si_rec) |>
  add_model(rfr_spec) |>
  tune_grid(resamples = folds,
            grid = rfr_grid,
            metrics = val_metrics)
```

    ## Warning: package 'ranger' was built under R version 4.3.2

``` r
rfr_si_tune |>
  collect_metrics() |>
  ggplot(aes(x = trees, y = mean, color = .metric)) +
  geom_errorbar(aes(ymin = mean - std_err,
                    ymax = mean + std_err),
                width = 0.05) +
  geom_line() +
  facet_wrap(. ~ .metric, ncol = 1, scales="free") +
  theme(legend.position = "none") + scale_x_log10()
```

![](model-expl_files/figure-gfm/rfr-si-tune-1.png)<!-- -->

``` r
rfr_spec <- rand_forest(mode = "regression", trees = 2^8)

rfr_si <- workflow() |>
  add_recipe(si_rec) |>
  add_model(rfr_spec) |>
  fit(data=training(td_split))

rfr_si$fit$fit |> print()
```

    ## parsnip model object
    ## 
    ## Ranger result
    ## 
    ## Call:
    ##  ranger::ranger(x = maybe_data_frame(x), y = y, num.trees = ~2^8,      num.threads = 1, verbose = FALSE, seed = sample.int(10^5,          1)) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  256 
    ## Sample size:                      553 
    ## Number of independent variables:  24 
    ## Mtry:                             4 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.007243483 
    ## R squared (OOB):                  0.8352352

``` r
rfr_si_res <- 
testing(td_split) |>
  mutate(hsi_frac_pred = predict(rfr_si, testing(td_split))[[".pred"]]) |>
  select(comid, hsi_frac, hsi_frac_pred, flow_norm_cfs) 

rfr_si_res |> 
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-independent model, random forest regression")
```

![](model-expl_files/figure-gfm/rfr-si-1.png)<!-- -->

``` r
rfr_si_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs)) |> #, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*100, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*100, color="actual", group=1)) + 
  xlab("Flow (cfs)") + ylab("Suitable Habitat Area (% of Total Area)")
```

![](model-expl_files/figure-gfm/rfr-si-2.png)<!-- -->

#### Scale-dependent (WUA-per-linear-ft versus flow)

``` r
rfr_spec <- rand_forest(mode = "regression", trees = tune())
rfr_grid <- tibble(trees = 2^seq(0, 12, 1))

rfr_sd_tune <- 
  workflow() |>
  add_recipe(sd_rec) |>
  add_model(rfr_spec) |>
  tune_grid(resamples = folds,
            grid = rfr_grid,
            metrics = val_metrics)
rfr_sd_tune |>
  collect_metrics() |>
  ggplot(aes(x = trees, y = mean, color = .metric)) +
  geom_errorbar(aes(ymin = mean - std_err,
                    ymax = mean + std_err),
                width = 0.05) +
  geom_line() +
  facet_wrap(. ~ .metric, ncol = 1, scales="free") +
  theme(legend.position = "none") + scale_x_log10()
```

![](model-expl_files/figure-gfm/rfr-sd-tune-1.png)<!-- -->

``` r
rfr_spec <- rand_forest(mode = "regression", trees = 2^8)

rfr_sd <- workflow() |>
  add_recipe(sd_rec) |>
  add_model(rfr_spec) |>
  fit(data=training(td_split))

rfr_sd$fit$fit |> print()
```

    ## parsnip model object
    ## 
    ## Ranger result
    ## 
    ## Call:
    ##  ranger::ranger(x = maybe_data_frame(x), y = y, num.trees = ~2^8,      num.threads = 1, verbose = FALSE, seed = sample.int(10^5,          1)) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  256 
    ## Sample size:                      553 
    ## Number of independent variables:  33 
    ## Mtry:                             5 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.0575351 
    ## R squared (OOB):                  0.9892562

``` r
rfr_sd_res <-
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(rfr_sd, testing(td_split))[[".pred"]]) |>
  transmute(comid, wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) 

rfr_sd_res |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c()   + 
  ggtitle("Scale-dependent model, random forest regression")
```

![](model-expl_files/figure-gfm/rfr-sd-1.png)<!-- -->

``` r
rfr_sd_res |> filter(comid %in% head(unique(testing(td_split)$comid),12)) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() + xlab("Flow (cfs)") + ylab("Suitable Habitat Area per LF Channel (ft)") 
```

![](model-expl_files/figure-gfm/rfr-sd-2.png)<!-- -->

## Random Forest model prediction output

``` r
# assemble prediction dataset
pd <- flowline_attributes |> 
  # filter for just the areas that are hydrologically similar to the training data
  filter(hyd_cls %in% c("High-volume snowmelt and rain",
                        "Low-volume snowmelt and rain")) |>
  # select the variables that are used in the model and drop NAs
  mutate(nf_bfl_dry_cfs_norm = coalesce(nf_bfl_dry_cfs/erom_q_ma_cfs,0), # normalized baseflow values
         nf_bfl_wet_cfs_norm = coalesce(nf_bfl_wet_cfs/erom_q_ma_cfs,0)) |>
  select(comid, any_of(c(sd_rec$var_info$variable, si_rec$var_info$variable))) |> drop_na() |>
  # create series of flow prediction points
  expand_grid(flow_cfs = c(0,signif(100*2^seq(-2,7,0.5),2))) |>
  mutate(flow_cfs = if_else(flow_cfs>0, flow_cfs, erom_q_ma_cfs)) |> # also evaluate at mean annual flow
  mutate(flow_norm_cfs = coalesce(flow_cfs/erom_q_ma_cfs,0)) |> # flow as a percent of mean annual flow
  arrange(comid, flow_cfs) |>
  glimpse()
```

    ## Rows: 244,180
    ## Columns: 32
    ## $ comid                    <dbl> 342439, 342439, 342439, 342439, 342439, 34243…
    ## $ slope                    <dbl> 0.00931547, 0.00931547, 0.00931547, 0.0093154…
    ## $ da_area_sq_km            <dbl> 8.3835, 8.3835, 8.3835, 8.3835, 8.3835, 8.383…
    ## $ da_elev_mean             <dbl> 2284.09, 2284.09, 2284.09, 2284.09, 2284.09, …
    ## $ da_ppt_mean_mm           <dbl> 1248.088, 1248.088, 1248.088, 1248.088, 1248.…
    ## $ nf_bfl_dry_cfs           <dbl> 0.864, 0.864, 0.864, 0.864, 0.864, 0.864, 0.8…
    ## $ nf_bfl_wet_cfs           <dbl> 4.33, 4.33, 4.33, 4.33, 4.33, 4.33, 4.33, 4.3…
    ## $ erom_q_ma_cfs            <dbl> 4.109, 4.109, 4.109, 4.109, 4.109, 4.109, 4.1…
    ## $ erom_v_ma_fps            <dbl> 0.88129, 0.88129, 0.88129, 0.88129, 0.88129, …
    ## $ bf_depth_m               <dbl> 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.7…
    ## $ bf_w_d_ratio             <dbl> 8.460526, 8.460526, 8.460526, 8.460526, 8.460…
    ## $ da_k_erodibility         <dbl> 0.1827, 0.1827, 0.1827, 0.1827, 0.1827, 0.182…
    ## $ da_avg_slope             <dbl> 19.95, 19.95, 19.95, 19.95, 19.95, 19.95, 19.…
    ## $ mean_ndvi                <dbl> 0.2641405, 0.2641405, 0.2641405, 0.2641405, 0…
    ## $ loc_bfi                  <dbl> 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 5…
    ## $ loc_pct_clay             <dbl> 13.60, 13.60, 13.60, 13.60, 13.60, 13.60, 13.…
    ## $ loc_pct_sand             <dbl> 57.30, 57.30, 57.30, 57.30, 57.30, 57.30, 57.…
    ## $ loc_permeability         <dbl> 10.16, 10.16, 10.16, 10.16, 10.16, 10.16, 10.…
    ## $ loc_bedrock_depth        <dbl> 123.50, 123.50, 123.50, 123.50, 123.50, 123.5…
    ## $ loc_ppt_mean_mm          <dbl> 1197.80, 1197.80, 1197.80, 1197.80, 1197.80, …
    ## $ mtpi30_min               <dbl> -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -…
    ## $ vb_bf_w_ratio            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ frac_leveed_longitudinal <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ vb_width_transect        <dbl> 6.43, 6.43, 6.43, 6.43, 6.43, 6.43, 6.43, 6.4…
    ## $ sinuosity                <dbl> 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.3…
    ## $ bio_aq_rank_sw           <int> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, …
    ## $ species                  <dbl> 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 6…
    ## $ species_fish             <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ nf_bfl_dry_cfs_norm      <dbl> 0.2102701, 0.2102701, 0.2102701, 0.2102701, 0…
    ## $ nf_bfl_wet_cfs_norm      <dbl> 1.0537844, 1.0537844, 1.0537844, 1.0537844, 1…
    ## $ flow_cfs                 <dbl> 4.109, 25.000, 35.000, 50.000, 71.000, 100.00…
    ## $ flow_norm_cfs            <dbl> 1.000000, 6.084205, 8.517888, 12.168411, 17.2…

``` r
# # get min/max ranges from the training data
# ranges <- td |> select(all_of(names(pd))) |> 
#   reframe(across(everything(), function(x) (range(x)))) 
# # filter the prediction data for just those falling within the training data
# pd_idx <- pd |> 
#   select(comid, all_of(names(ranges))) |>
#   mutate(across(-comid, function(x) x>ranges[[cur_column()]][[1]] & x<=ranges[[cur_column()]][[2]])) |>
#   filter(if_all(-comid)) |>
#   pull(comid)
# # update pd to exclude the out-of-range data
# pd <- pd |> filter(comid %in% pd_idx) |> glimpse()
```

``` r
pd_rf <- pd |> #head(100) |>
  mutate(log_wua_per_lf_pred = predict(rfr_sd, pd)[[".pred"]],
         wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs,
         hsi_frac_pred = predict(rfr_si, pd)[[".pred"]]) 

pd_rf |> select(comid, 
                flow_cfs, wua_per_lf_pred, 
                flow_norm_cfs, hsi_frac_pred) |> glimpse()
```

    ## Rows: 244,180
    ## Columns: 5
    ## $ comid           <dbl> 342439, 342439, 342439, 342439, 342439, 342439, 342439…
    ## $ flow_cfs        <dbl> 4.109, 25.000, 35.000, 50.000, 71.000, 100.000, 140.00…
    ## $ wua_per_lf_pred <dbl> 3.849073, 3.904366, 3.904366, 3.904366, 3.946281, 3.94…
    ## $ flow_norm_cfs   <dbl> 1.000000, 6.084205, 8.517888, 12.168411, 17.279143, 24…
    ## $ hsi_frac_pred   <dbl> 0.6146102, 0.4218199, 0.4151587, 0.4129182, 0.4116372,…

``` r
pd_rf |> saveRDS("../data/rf_predictions.Rds")
```

Example plot of habitat area per linear foot at 800 cfs (scale-dependent
model)

``` r
flowlines |> st_zm() |>
  inner_join(pd_rf |> filter(flow_cfs==800), by=join_by(comid)) |>
  ggplot() + geom_sf(aes(color = wua_per_lf_pred)) + 
  scale_color_viridis_c(trans="log")
```

![](model-expl_files/figure-gfm/map-output-sd-1.png)<!-- -->

Example plot of percent habitat area at mean annual flow
(scale-independent model)

``` r
flowlines |> st_zm() |>
  inner_join(pd_rf |> filter(flow_norm_cfs==1), by=join_by(comid)) |>
  ggplot() + geom_sf(aes(color = hsi_frac_pred)) + 
  scale_color_viridis_c(trans="log")
```

![](model-expl_files/figure-gfm/map-output-si-1.png)<!-- -->

## Predictions summarized by DSMHabitat reach

``` r
# this is a rough join method, just taking any COMIDs within 50 ft buffer of flowline
mainstems <- 
  st_read("/vsizip/rearing_spatial_data/habitat_extents_combined_gradients_v3.shp.zip", as_tibble=T) |>
  janitor::clean_names() |>
  group_by(habitat, species, river) |>
  summarize() |>
  st_buffer(dist=50) |>
  st_union(by_feature=T) |>
  st_transform(st_crs(flowlines)) |>
  filter(habitat=="rearing" & species=="Fall Run Chinook")
```

    ## Reading layer `habitat_extents_combined_gradients_v3' from data source 
    ##   `/vsizip/rearing_spatial_data/habitat_extents_combined_gradients_v3.shp.zip' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 323 features and 12 fields
    ## Geometry type: MULTILINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: -251634.9 ymin: -139292.4 xmax: 61284.41 ymax: 364940.4
    ## Projected CRS: NAD83 / California Albers

    ## `summarise()` has grouped output by 'habitat', 'species'. You can override
    ## using the `.groups` argument.

``` r
mainstems_comid <- 
  flowlines |> st_zm() |>
  st_join(mainstems, join=st_intersects, left=F) |>
  mutate(length_ft = st_length(geometry) |> units::set_units("ft") |> units::drop_units())

mainstems_comid |> ggplot() + geom_sf(aes(color=river))
```

![](model-expl_files/figure-gfm/dsmhabitat-prep-1.png)<!-- -->

``` r
pd_dsmhabitat <- flowline_attributes |> 
  # select the variables that are used in the model and drop NAs
  select(comid, any_of(sd_rec$var_info$variable)) |> 
  # create series of flow prediction points
  expand_grid(flow_cfs = c(100,250,300,400,500,600,1000,3000,5000,6000,7000,9000,10000,11000,12000,13000,14000,15000)) |>
  arrange(comid, flow_cfs) |>
  filter(comid %in% mainstems_comid$comid) |>
  drop_na() 

pd_dsmhabitat_pred <-
  pd_dsmhabitat|>
  mutate(log_wua_per_lf_pred = predict(rfr_sd, pd_dsmhabitat)[[".pred"]],
         wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) |>
  select(comid, 
         flow_cfs, wua_per_lf_pred) |> glimpse()
```

    ## Rows: 32,544
    ## Columns: 3
    ## $ comid           <dbl> 348543, 348543, 348543, 348543, 348543, 348543, 348543…
    ## $ flow_cfs        <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 5000, 6000, …
    ## $ wua_per_lf_pred <dbl> 3.9255960, 3.9818639, 3.9818639, 4.0561232, 4.1002234,…

``` r
pd_rf_by_mainstem <-
  mainstems_comid |> 
  st_drop_geometry() |>
  inner_join(pd_dsmhabitat_pred, by=join_by(comid), relationship="many-to-many") |>
  group_by(habitat, species, river, flow_cfs) |>
  summarize(tot_wua_per_lf_pred = sum(coalesce(wua_per_lf_pred,0)),
            tot_length_ft = sum(coalesce(length_ft,0))) |>
  mutate(tot_wua_ft2 = tot_wua_per_lf_pred * tot_length_ft,
         tot_wua_ac = tot_wua_ft2 / 43560,
         avg_wua_ft2_per_1000ft = 1000 * tot_wua_ft2 / tot_length_ft) |>
  glimpse()
```

    ## `summarise()` has grouped output by 'habitat', 'species', 'river'. You can
    ## override using the `.groups` argument.

    ## Rows: 540
    ## Columns: 9
    ## Groups: habitat, species, river [30]
    ## $ habitat                <chr> "rearing", "rearing", "rearing", "rearing", "re…
    ## $ species                <chr> "Fall Run Chinook", "Fall Run Chinook", "Fall R…
    ## $ river                  <chr> "American River", "American River", "American R…
    ## $ flow_cfs               <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 5000,…
    ## $ tot_wua_per_lf_pred    <dbl> 280.8356, 285.6920, 285.7659, 287.7898, 291.378…
    ## $ tot_length_ft          <dbl> 142202.8, 142202.8, 142202.8, 142202.8, 142202.…
    ## $ tot_wua_ft2            <dbl> 39935603, 40626194, 40636707, 40924506, 4143476…
    ## $ tot_wua_ac             <dbl> 916.7953, 932.6491, 932.8904, 939.4974, 951.211…
    ## $ avg_wua_ft2_per_1000ft <dbl> 280835.6, 285692.0, 285765.9, 287789.8, 291378.…

``` r
pd_rf_by_mainstem |>
  ggplot() + 
  geom_line(aes(x = flow_cfs, y = tot_wua_ac, color=species, linetype=habitat)) + 
  facet_wrap(~river, scales="free_y") + 
  scale_x_log10(breaks=c(100,300,1000,3000,10000)) + 
  scale_y_log10() + theme(legend.position="top", panel.grid.minor = element_blank()) 
```

![](model-expl_files/figure-gfm/dsmhabitat-val-1.png)<!-- -->

``` r
watersheds <- pd_rf_by_mainstem |> 
  filter(river != "South Delta") |> 
  pull(river) |> unique()

watershed_name <- tolower(gsub(pattern = "-| ", replacement = "_", x = watersheds))
watershed_rda_name <- paste(watershed_name, "floodplain", sep = "_")

dsm_habitat_floodplain <- map_df(watershed_rda_name, function(watershed) {
  df <- as.data.frame(do.call(`::`, list(pkg = "DSMhabitat", name = watershed)))
}) |> 
  rename(river = watershed,
         flow_cfs_dsm = flow_cfs) |>
  mutate(FR_floodplain_m2 = FR_floodplain_acres * 4046.86,
         FR_floodplain_m2_suitable = DSMhabitat::apply_suitability(FR_floodplain_m2),
         FR_floodplain_acres_suitable = FR_floodplain_m2_suitable / 4046.86)

pd_rf_by_mainstem |> 
  left_join(dsm_habitat_floodplain) |> 
  #filter(river == "American River") |> 
  ggplot() +
  geom_line(aes(x = flow_cfs, y = tot_wua_ac, color='modeled')) + 
  geom_line(aes(x = flow_cfs_dsm, y = FR_floodplain_acres_suitable, color = 'DSMhabitat')) + 
  facet_wrap(~river) +#, scales="free_y") + 
  scale_x_log10() + #breaks=c(100,300,1000,3000,10000)) + 
  scale_y_log10() + theme(legend.position="top", panel.grid.minor = element_blank()) 
```

    ## Joining with `by = join_by(river)`

    ## Warning in left_join(pd_rf_by_mainstem, dsm_habitat_floodplain): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 1 of `x` matches multiple rows in `y`.
    ## ℹ Row 1 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](model-expl_files/figure-gfm/dsmhabitat-watershed-comp-1.png)<!-- -->

``` r
watersheds <- pd_rf_by_mainstem |> 
  filter(river != "South Delta") |> 
  pull(river) |> unique()

watershed_name <- tolower(gsub(pattern = "-| ", replacement = "_", x = watersheds))
watershed_rda_name <- paste(watershed_name, "instream", sep = "_")

dsm_habitat_instream <- map_df(paste(watershed_name, "instream", sep = "_"), 
                               possibly(function(watershed) {
                                 df <- as.data.frame(do.call(`::`, list(pkg = "DSMhabitat", name = watershed)))
                                 }, otherwise = NULL)) |> rename(river = watershed) 

ggplot() +
  geom_line(data=pd_rf_by_mainstem, aes(x = flow_cfs, y = avg_wua_ft2_per_1000ft, color='modeled')) + 
  geom_line(data=dsm_habitat_instream, aes(x = flow_cfs, y = FR_juv_wua, color = 'DSMhabitat')) + 
  facet_wrap(~river) +#, scales="free_y") + 
  scale_x_log10() + 
  scale_y_log10() + theme(legend.position="top", panel.grid.minor = element_blank()) 
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 5 rows containing missing values (`geom_line()`).

![](model-expl_files/figure-gfm/dsmhabitat-watershed-comp-instream-1.png)<!-- -->

===

``` r
knitr::knit_exit()
```
