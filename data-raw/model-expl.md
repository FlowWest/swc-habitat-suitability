Exploratory Modeling
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-02-13

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
- [TESTS](#tests)
  - [Series interpolated (panel)
    regression](#series-interpolated-panel-regression)
  - [Use regression to predict FSA curve
    shape](#use-regression-to-predict-fsa-curve-shape)

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
    ## Columns: 86

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

    ## $ comid                     <dbl> 12071480, 12071480, 12071480, 12071480, 1207…
    ## $ reachcode                 <fct> 18020157005039, 18020157005039, 180201570050…
    ## $ gnis_id                   <fct> 1655075, 1655075, 1655075, 1655075, 1655075,…
    ## $ gnis_name                 <fct> "Deer Creek", "Deer Creek", "Deer Creek", "D…
    ## $ lengthkm                  <dbl> 3.414, 3.414, 3.414, 3.414, 3.414, 3.414, 3.…
    ## $ ftype                     <fct> StreamRiver, StreamRiver, StreamRiver, Strea…
    ## $ fcode                     <int> 46006, 46006, 46006, 46006, 46006, 46006, 46…
    ## $ huc_8                     <chr> "18020157", "18020157", "18020157", "1802015…
    ## $ huc_10                    <chr> "1802015700", "1802015700", "1802015700", "1…
    ## $ huc_12                    <chr> "180201570050", "180201570050", "18020157005…
    ## $ ftype_desc                <chr> "Stream/River", "Stream/River", "Stream/Rive…
    ## $ hydro_seq                 <dbl> 10012820, 10012820, 10012820, 10012820, 1001…
    ## $ reach_code                <fct> 18020157005039, 18020157005039, 180201570050…
    ## $ stream_level              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ stream_order              <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
    ## $ us_length_km              <dbl> 389.906, 389.906, 389.906, 389.906, 389.906,…
    ## $ ds_length_km              <dbl> 346.603, 346.603, 346.603, 346.603, 346.603,…
    ## $ da_area_sq_km             <dbl> 573.3900, 573.3900, 573.3900, 573.3900, 573.…
    ## $ reach_length_km           <dbl> 3.414, 3.414, 3.414, 3.414, 3.414, 3.414, 3.…
    ## $ slope                     <dbl> 0.00272993, 0.00272993, 0.00272993, 0.002729…
    ## $ elev_min                  <dbl> 51.80, 51.80, 51.80, 51.80, 51.80, 51.80, 51…
    ## $ elev_max                  <dbl> 61.12, 61.12, 61.12, 61.12, 61.12, 61.12, 61…
    ## $ stream_power              <dbl> 1.565314563, 1.565314563, 1.565314563, 1.565…
    ## $ da_ppt_mean_mm            <dbl> 1437.351, 1437.351, 1437.351, 1437.351, 1437…
    ## $ loc_ppt_mean_mm           <dbl> 584.914, 584.914, 584.914, 584.914, 584.914,…
    ## $ vogel_q_ma_cfs            <dbl> 593.8086, 593.8086, 593.8086, 593.8086, 593.…
    ## $ vogel_v_ma_fps            <dbl> 1.63041, 1.63041, 1.63041, 1.63041, 1.63041,…
    ## $ erom_q_ma_cfs             <dbl> 347.468, 347.468, 347.468, 347.468, 347.468,…
    ## $ erom_v_ma_fps             <dbl> 1.57425, 1.57425, 1.57425, 1.57425, 1.57425,…
    ## $ sinuosity                 <dbl> 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.…
    ## $ da_avg_slope              <dbl> 25.97, 25.97, 25.97, 25.97, 25.97, 25.97, 25…
    ## $ da_elev_mean              <dbl> 1206.24, 1206.24, 1206.24, 1206.24, 1206.24,…
    ## $ da_elev_min               <dbl> 51.80, 51.80, 51.80, 51.80, 51.80, 51.80, 51…
    ## $ da_elev_max               <dbl> 2393.52, 2393.52, 2393.52, 2393.52, 2393.52,…
    ## $ da_elev_rel               <dbl> 2341.72, 2341.72, 2341.72, 2341.72, 2341.72,…
    ## $ loc_bfi                   <dbl> 42.5746, 42.5746, 42.5746, 42.5746, 42.5746,…
    ## $ da_bfi                    <dbl> 63.3563, 63.3563, 63.3563, 63.3563, 63.3563,…
    ## $ loc_twi                   <dbl> 977.304, 977.304, 977.304, 977.304, 977.304,…
    ## $ da_twi                    <dbl> 676.9103, 676.9103, 676.9103, 676.9103, 676.…
    ## $ loc_pct_clay              <dbl> 13.3988, 13.3988, 13.3988, 13.3988, 13.3988,…
    ## $ da_pct_clay               <dbl> 21.6023, 21.6023, 21.6023, 21.6023, 21.6023,…
    ## $ loc_pct_sand              <dbl> 48.7830, 48.7830, 48.7830, 48.7830, 48.7830,…
    ## $ da_pct_sand               <dbl> 38.4830, 38.4830, 38.4830, 38.4830, 38.4830,…
    ## $ loc_permeability          <dbl> 7.306, 7.306, 7.306, 7.306, 7.306, 7.306, 7.…
    ## $ loc_bedrock_depth         <dbl> 151.2547, 151.2547, 151.2547, 151.2547, 151.…
    ## $ da_permeability           <dbl> 6.1970, 6.1970, 6.1970, 6.1970, 6.1970, 6.19…
    ## $ da_bedrock_depth          <dbl> 105.5002, 105.5002, 105.5002, 105.5002, 105.…
    ## $ loc_k_erodibility         <dbl> 0.2617, 0.2617, 0.2617, 0.2617, 0.2617, 0.26…
    ## $ da_k_erodibility          <dbl> 0.2964, 0.2964, 0.2964, 0.2964, 0.2964, 0.29…
    ## $ loc_precip                <dbl> 582.4351, 582.4351, 582.4351, 582.4351, 582.…
    ## $ da_precip                 <dbl> 1417.543, 1417.543, 1417.543, 1417.543, 1417…
    ## $ loc_runoff                <dbl> 521, 521, 521, 521, 521, 521, 521, 521, 521,…
    ## $ da_runoff                 <dbl> 822.3441, 822.3441, 822.3441, 822.3441, 822.…
    ## $ bf_width_m                <dbl> 31.47, 31.47, 31.47, 31.47, 31.47, 31.47, 31…
    ## $ bf_depth_m                <dbl> 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.…
    ## $ bf_w_d_ratio              <dbl> 18.29651, 18.29651, 18.29651, 18.29651, 18.2…
    ## $ mean_ndvi                 <dbl> 0.4869165, 0.4869165, 0.4869165, 0.4869165, …
    ## $ hyd_cat                   <fct> Mixed, Mixed, Mixed, Mixed, Mixed, Mixed, Mi…
    ## $ hyd_cls                   <fct> "Low-volume snowmelt and rain", "Low-volume …
    ## $ nf_bfl_dry_cfs            <dbl> 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74…
    ## $ nf_bfl_wet_cfs            <dbl> 178, 178, 178, 178, 178, 178, 178, 178, 178,…
    ## $ nf_wet_start              <dbl> 85.2, 85.2, 85.2, 85.2, 85.2, 85.2, 85.2, 85…
    ## $ nf_wet_dur                <dbl> 136, 136, 136, 136, 136, 136, 136, 136, 136,…
    ## $ peak_q2_cfs               <dbl> 3870, 3870, 3870, 3870, 3870, 3870, 3870, 38…
    ## $ peak_q5_cfs               <dbl> 7430, 7430, 7430, 7430, 7430, 7430, 7430, 74…
    ## $ peak_q10_cfs              <dbl> 8720, 8720, 8720, 8720, 8720, 8720, 8720, 87…
    ## $ merit_width_m             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ vb_width_transect         <dbl> 31.47, 31.47, 31.47, 31.47, 31.47, 31.47, 31…
    ## $ bf_xarea_m                <dbl> 54.1284, 54.1284, 54.1284, 54.1284, 54.1284,…
    ## $ chan_width_m              <dbl> 31.47, 31.47, 31.47, 31.47, 31.47, 31.47, 31…
    ## $ velocity_m_s              <dbl> 5.164862, 5.164862, 5.164862, 5.164862, 5.16…
    ## $ wetted_perimeter_m        <dbl> 34.91, 34.91, 34.91, 34.91, 34.91, 34.91, 34…
    ## $ hydraulic_radius_m        <dbl> 1.550513, 1.550513, 1.550513, 1.550513, 1.55…
    ## $ critical_shields_number   <dbl> 0.03428697, 0.03428697, 0.03428697, 0.034286…
    ## $ grain_size_mobilized_mm   <dbl> 0.74819322, 0.74819322, 0.74819322, 0.748193…
    ## $ shear_velocity_cm_s       <dbl> 20.377361, 20.377361, 20.377361, 20.377361, …
    ## $ settling_velocity_ndim    <dbl> 522.746144, 522.746144, 522.746144, 522.7461…
    ## $ grain_size_suspended_ndim <dbl> 1746.03995, 1746.03995, 1746.03995, 1746.039…
    ## $ grain_size_suspended_mm   <dbl> 0.148708751, 0.148708751, 0.148708751, 0.148…
    ## $ mtpi30_min                <dbl> -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, …
    ## $ vb_bf_w_ratio             <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ dataset                   <chr> "Deer Creek", "Deer Creek", "Deer Creek", "D…
    ## $ flow_cfs                  <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 50…
    ## $ area_tot_ft2              <dbl> 165174, 195682, 218107, 232765, 243522, 2551…
    ## $ area_wua_ft2              <dbl> 126770, 155882, 166777, 171971, 174341, 1787…
    ## $ hsi_frac                  <dbl> 0.7674937, 0.7966088, 0.7646568, 0.7388181, …

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
         mtpi30_min, #vb_width_transect, vb_bf_w_ratio, 
         # channel sinuosity
         sinuosity,
         ) |> 
  mutate(nf_bfl_dry_cfs_norm = nf_bfl_dry_cfs/erom_q_ma_cfs, 
         nf_bfl_wet_cfs_norm = nf_bfl_wet_cfs/erom_q_ma_cfs) |>
  drop_na() |> glimpse()
```

    ## Rows: 738
    ## Columns: 33
    ## $ dataset             <chr> "Deer Creek", "Deer Creek", "Deer Creek", "Deer Cr…
    ## $ comid               <dbl> 12071480, 12071480, 12071480, 12071480, 12071480, …
    ## $ hsi_frac            <dbl> 0.7674937, 0.7966088, 0.7646568, 0.7388181, 0.7159…
    ## $ wua_per_lf          <dbl> 11.31795, 13.91706, 14.88976, 15.35347, 15.56507, …
    ## $ log_wua_per_lf      <dbl> 2.426390, 2.633115, 2.700673, 2.731342, 2.745029, …
    ## $ flow_cfs            <dbl> 100, 250, 300, 400, 500, 600, 1000, 3000, 5000, 60…
    ## $ flow_norm_cfs       <dbl> 0.2877963, 0.7194907, 0.8633889, 1.1511851, 1.4389…
    ## $ slope               <dbl> 0.00272993, 0.00272993, 0.00272993, 0.00272993, 0.…
    ## $ da_area_sq_km       <dbl> 573.3900, 573.3900, 573.3900, 573.3900, 573.3900, …
    ## $ da_elev_mean        <dbl> 1206.24, 1206.24, 1206.24, 1206.24, 1206.24, 1206.…
    ## $ da_ppt_mean_mm      <dbl> 1437.351, 1437.351, 1437.351, 1437.351, 1437.351, …
    ## $ nf_bfl_dry_cfs      <dbl> 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74.3, 74…
    ## $ nf_bfl_wet_cfs      <dbl> 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, …
    ## $ erom_q_ma_cfs       <dbl> 347.468, 347.468, 347.468, 347.468, 347.468, 347.4…
    ## $ peak_q2_cfs         <dbl> 3870, 3870, 3870, 3870, 3870, 3870, 3870, 3870, 38…
    ## $ peak_q5_cfs         <dbl> 7430, 7430, 7430, 7430, 7430, 7430, 7430, 7430, 74…
    ## $ peak_q10_cfs        <dbl> 8720, 8720, 8720, 8720, 8720, 8720, 8720, 8720, 87…
    ## $ bf_depth_m          <dbl> 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72, 1.…
    ## $ bf_w_d_ratio        <dbl> 18.29651, 18.29651, 18.29651, 18.29651, 18.29651, …
    ## $ erom_v_ma_fps       <dbl> 1.57425, 1.57425, 1.57425, 1.57425, 1.57425, 1.574…
    ## $ da_avg_slope        <dbl> 25.97, 25.97, 25.97, 25.97, 25.97, 25.97, 25.97, 2…
    ## $ da_k_erodibility    <dbl> 0.2964, 0.2964, 0.2964, 0.2964, 0.2964, 0.2964, 0.…
    ## $ mean_ndvi           <dbl> 0.4869165, 0.4869165, 0.4869165, 0.4869165, 0.4869…
    ## $ loc_bfi             <dbl> 42.5746, 42.5746, 42.5746, 42.5746, 42.5746, 42.57…
    ## $ loc_pct_clay        <dbl> 13.3988, 13.3988, 13.3988, 13.3988, 13.3988, 13.39…
    ## $ loc_pct_sand        <dbl> 48.7830, 48.7830, 48.7830, 48.7830, 48.7830, 48.78…
    ## $ loc_permeability    <dbl> 7.306, 7.306, 7.306, 7.306, 7.306, 7.306, 7.306, 7…
    ## $ loc_bedrock_depth   <dbl> 151.2547, 151.2547, 151.2547, 151.2547, 151.2547, …
    ## $ loc_ppt_mean_mm     <dbl> 584.914, 584.914, 584.914, 584.914, 584.914, 584.9…
    ## $ mtpi30_min          <dbl> -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5…
    ## $ sinuosity           <dbl> 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.27, 1.…
    ## $ nf_bfl_dry_cfs_norm <dbl> 0.2138326, 0.2138326, 0.2138326, 0.2138326, 0.2138…
    ## $ nf_bfl_wet_cfs_norm <dbl> 0.5122774, 0.5122774, 0.5122774, 0.5122774, 0.5122…

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
        mtpi30_min + 
        # baseflow as percent of mean annual flow
        nf_bfl_dry_cfs_norm + nf_bfl_wet_cfs_norm 
        ) |>
  step_log(all_numeric_predictors(), -mtpi30_min) |>
  # let effect of flow vary with gradient
  # step_interact(terms = ~ flow_norm_cfs:slope)
  # try interacting with all numeric predictors
  step_interact(terms = ~ flow_norm_cfs:all_numeric_predictors())

lm_si <- workflow() |>
  add_recipe(si_rec) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_si |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic   p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>     <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.727         0.709 0.113      39.4 1.13e-122    35   439. -804. -644.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_si |> tidy()
```

    ## # A tibble: 36 × 5
    ##    term             estimate std.error statistic p.value
    ##    <chr>               <dbl>     <dbl>     <dbl>   <dbl>
    ##  1 (Intercept)       -1.95      2.78      -0.703 0.482  
    ##  2 flow_norm_cfs     -0.383     2.22      -0.173 0.863  
    ##  3 slope             -0.0314    0.0993    -0.316 0.752  
    ##  4 sinuosity         -0.0830    0.0401    -2.07  0.0388 
    ##  5 erom_v_ma_fps      0.351     0.834      0.421 0.674  
    ##  6 bf_depth_m        -0.596     0.317     -1.88  0.0604 
    ##  7 bf_w_d_ratio      -0.751     0.273     -2.76  0.00607
    ##  8 da_k_erodibility  -0.737     1.31      -0.564 0.573  
    ##  9 da_avg_slope       0.528     1.15       0.459 0.646  
    ## 10 mean_ndvi          0.268     0.118      2.26  0.0240 
    ## # ℹ 26 more rows

``` r
lm_si_res <- testing(td_split) |>
  mutate(hsi_frac_pred = predict(lm_si, testing(td_split))[[".pred"]]) |>
  select(comid, hsi_frac, hsi_frac_pred, flow_norm_cfs) 

lm_si_res |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-independent model, linear regression")
```

![](model-expl_files/figure-gfm/lm-si-1.png)<!-- -->

``` r
lm_si_res |> filter(comid<=2821716) |> left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*chan_width_ft, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*chan_width_ft, color="actual", group=1)) + 
  scale_y_log10()
```

![](model-expl_files/figure-gfm/lm-si-2.png)<!-- -->

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
        mtpi30_min + #log(vb_width_transect) + log(vb_bf_w_ratio) +
        # channel sinuosity
        sinuosity
        ) |>
  step_log(all_numeric_predictors(), -mtpi30_min) |>
  step_interact(terms = ~ flow_cfs:erom_q_ma_cfs) |>
  step_interact(terms = ~ slope:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_elev_mean) |>
  step_interact(terms = ~ flow_cfs:da_ppt_mean_mm)

lm_sd <- workflow() |>
  add_recipe(sd_rec) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_sd |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.961         0.959 0.469      477.       0    27  -352.  762.  887.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_sd |> tidy()
```

    ## # A tibble: 28 × 5
    ##    term            estimate std.error statistic  p.value
    ##    <chr>              <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)    -1116.      519.     -2.15    3.19e- 2
    ##  2 flow_cfs         -18.9       8.38   -2.25    2.47e- 2
    ##  3 slope              0.821     0.307   2.68    7.67e- 3
    ##  4 da_area_sq_km    -14.9       2.83   -5.25    2.18e- 7
    ##  5 da_elev_mean       0.700    83.6     0.00837 9.93e- 1
    ##  6 da_ppt_mean_mm   164.       37.6     4.37    1.52e- 5
    ##  7 nf_bfl_dry_cfs    36.5      10.0     3.63    3.05e- 4
    ##  8 nf_bfl_wet_cfs    43.1       6.07    7.09    4.27e-12
    ##  9 erom_q_ma_cfs      8.45      2.06    4.11    4.57e- 5
    ## 10 erom_v_ma_fps      1.72      1.99    0.862   3.89e- 1
    ## # ℹ 18 more rows

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

![](model-expl_files/figure-gfm/lm-sd-1.png)<!-- -->

``` r
lm_sd_res |> filter(comid<=2821716) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() 
```

![](model-expl_files/figure-gfm/lm-sd-2.png)<!-- -->

### Linear regression with regularized (lasso) feature selection

``` r
#install.packages("glmnet")
lasso_spec <- linear_reg(penalty = 0.01, mixture = 1) |> set_engine("glmnet")
```

#### Scale-independent (%HSI versus dimensionless flow)

``` r
drop_vars <- workflow() |>
  add_recipe(si_rec |> 
               step_normalize(all_numeric_predictors())) |>
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
    ## 1     0.609         0.602 0.132      93.9 1.15e-104     9   339. -656. -609.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_si |> tidy()
```

    ## # A tibble: 10 × 5
    ##    term                              estimate std.error statistic  p.value
    ##    <chr>                                <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)                      -0.0567    0.374       -0.152 8.80e- 1
    ##  2 sinuosity                        -0.0895    0.0331      -2.71  7.00e- 3
    ##  3 erom_v_ma_fps                    -0.105     0.0364      -2.87  4.20e- 3
    ##  4 bf_depth_m                       -0.211     0.0387      -5.44  7.92e- 8
    ##  5 bf_w_d_ratio                     -0.132     0.0974      -1.35  1.76e- 1
    ##  6 loc_ppt_mean_mm                   0.210     0.0346       6.06  2.51e- 9
    ##  7 flow_norm_cfs_x_slope             0.0175    0.000891    19.7   1.69e-65
    ##  8 flow_norm_cfs_x_sinuosity        -0.0953    0.0244      -3.90  1.08e- 4
    ##  9 flow_norm_cfs_x_erom_v_ma_fps    -0.0275    0.00768     -3.58  3.73e- 4
    ## 10 flow_norm_cfs_x_loc_permeability -0.000394  0.00266     -0.148 8.83e- 1

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

![](model-expl_files/figure-gfm/lasso-si-1.png)<!-- -->

``` r
lasso_si_res |> filter(comid<=2821716) |> left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*chan_width_ft, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*chan_width_ft, color="actual", group=1)) + 
  scale_y_log10()
```

![](model-expl_files/figure-gfm/lasso-si-2.png)<!-- -->

#### Scale-dependent (WUA-per-linear-ft versus flow)

``` r
drop_vars <- workflow() |>
  add_recipe(sd_rec |> 
               step_normalize(all_numeric_predictors())) |>
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
    ## 1     0.939         0.937 0.580      688. 6.50e-318    12  -477.  981. 1042.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_sd |> tidy()
```

    ## # A tibble: 13 × 5
    ##    term                      estimate std.error statistic  p.value
    ##    <chr>                        <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)               -60.3     16.7        -3.61  3.41e- 4
    ##  2 da_elev_mean               -1.13     1.43       -0.788 4.31e- 1
    ##  3 da_ppt_mean_mm             11.3      1.10       10.2   1.33e-22
    ##  4 erom_v_ma_fps               0.536    0.182       2.94  3.41e- 3
    ##  5 bf_depth_m                  0.516    0.116       4.43  1.13e- 5
    ##  6 bf_w_d_ratio               -3.52     0.623      -5.65  2.55e- 8
    ##  7 loc_bfi                     0.455    0.428       1.06  2.88e- 1
    ##  8 loc_pct_sand                0.665    0.470       1.41  1.58e- 1
    ##  9 loc_permeability           -0.180    0.0925     -1.94  5.23e- 2
    ## 10 loc_bedrock_depth          -1.19     0.346      -3.44  6.20e- 4
    ## 11 mtpi30_min                  0.0479   0.00897     5.34  1.40e- 7
    ## 12 sinuosity                  -0.661    0.143      -4.63  4.56e- 6
    ## 13 flow_cfs_x_da_ppt_mean_mm   0.0264   0.00312     8.47  2.39e-16

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

![](model-expl_files/figure-gfm/lasso-sd-1.png)<!-- -->

``` r
lasso_sd_res |> filter(comid<=2821716) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() 
```

![](model-expl_files/figure-gfm/lasso-sd-2.png)<!-- -->

### Random Forest Regrssion

``` r
#install.packages("ranger")
rfr_spec <- rand_forest(mode = "regression", trees = 2^12)
```

#### Scale-independent (%HSI versus dimensionless flow)

``` r
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
    ##  ranger::ranger(x = maybe_data_frame(x), y = y, num.trees = ~2^12,      num.threads = 1, verbose = FALSE, seed = sample.int(10^5,          1)) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  4096 
    ## Sample size:                      553 
    ## Number of independent variables:  35 
    ## Mtry:                             5 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.005088926 
    ## R squared (OOB):                  0.8842441

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
rfr_si_res |> filter(comid<=2821716) |> left_join(transmute(flowline_attributes, comid, erom_q_ma_cfs, chan_width_ft=chan_width_m/0.3048)) |>
  ggplot(aes(x=flow_norm_cfs*erom_q_ma_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=hsi_frac_pred*chan_width_ft, color="predicted", group=1)) + geom_line(aes(y=hsi_frac*chan_width_ft, color="actual", group=1)) + 
  scale_y_log10()
```

![](model-expl_files/figure-gfm/rfr-si-2.png)<!-- -->

#### Scale-dependent (WUA-per-linear-ft versus flow)

``` r
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
    ##  ranger::ranger(x = maybe_data_frame(x), y = y, num.trees = ~2^12,      num.threads = 1, verbose = FALSE, seed = sample.int(10^5,          1)) 
    ## 
    ## Type:                             Regression 
    ## Number of trees:                  4096 
    ## Sample size:                      553 
    ## Number of independent variables:  27 
    ## Mtry:                             5 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.0597608 
    ## R squared (OOB):                  0.9888406

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
rfr_sd_res |> filter(comid<=2821716) |>
  ggplot(aes(x=flow_cfs, group=1)) + facet_wrap(~comid) + 
  geom_line(aes(y=wua_per_lf_pred, color="predicted", group=1)) + geom_line(aes(y=wua_per_lf, color="actual", group=1)) + 
  scale_y_log10() 
```

![](model-expl_files/figure-gfm/rfr-sd-2.png)<!-- -->

## TESTS

### Series interpolated (panel) regression

``` r
# interpolate series so there is a value fo every 100 flows
td_series <- td |>
  group_by(comid) |>
  arrange(comid, flow_cfs) |>
  complete(flow_cfs = seq(from = min(flow_cfs), to = max(flow_cfs), by = 100)) |>
  mutate(across(hsi_frac:last_col(), function(var) zoo::na.approx(var, x = flow_cfs))) |>
  filter(flow_cfs%%100==0) |> ungroup()

sd_series_rec <- recipe(data=td_series, 
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
        mtpi30_min + #log(vb_width_transect) + log(vb_bf_w_ratio) +
        # channel sinuosity
        sinuosity
        ) |>
  step_log(all_numeric_predictors(), -mtpi30_min) |>
  step_interact(terms = ~ flow_cfs:erom_q_ma_cfs) |>
  step_interact(terms = ~ slope:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_area_sq_km) |>
  step_interact(terms = ~ flow_cfs:da_elev_mean) |>
  step_interact(terms = ~ flow_cfs:da_ppt_mean_mm) |>
  step_lag(flow_cfs) |>
  step_interact(terms = ~ flow_cfs:lag_1_flow_cfs) 
  #step_interact(terms = ~ lag_1_flow_cfs:all_numeric_predictors())

lm_sd_series <- workflow() |>
  add_recipe(sd_series_rec) |>
  add_model(lm_spec) |>
  fit(data=td_series)

lm_sd_series |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.973         0.972 0.422     5557.       0    29 -2530. 5122. 5322.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_sd_series |> tidy()
```

    ## # A tibble: 30 × 5
    ##    term            estimate std.error statistic   p.value
    ##    <chr>              <dbl>     <dbl>     <dbl>     <dbl>
    ##  1 (Intercept)    -1707.     152.        -11.2  7.20e- 29
    ##  2 flow_cfs         -19.9      3.34       -5.96 2.65e-  9
    ##  3 slope              0.794    0.0931      8.53 2.01e- 17
    ##  4 da_area_sq_km    -23.0      1.29      -17.9  5.92e- 69
    ##  5 da_elev_mean      83.4     25.0         3.33 8.80e-  4
    ##  6 da_ppt_mean_mm   207.      11.1        18.6  2.38e- 74
    ##  7 nf_bfl_dry_cfs    41.8      3.19       13.1  1.09e- 38
    ##  8 nf_bfl_wet_cfs    50.5      1.83       27.7  6.58e-156
    ##  9 erom_q_ma_cfs     14.2      0.977      14.5  9.59e- 47
    ## 10 erom_v_ma_fps      3.80     0.635       5.98 2.33e-  9
    ## # ℹ 20 more rows

``` r
td_series |>
  mutate(log_wua_per_lf_pred = predict(lm_sd_series, td_series)[[".pred"]]) |>
  transmute(wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Series interpolated scale-dependent model, linear regression")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](model-expl_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Use regression to predict FSA curve shape

- model flow-to-suitable area curve shape with exponential function
- use spatial predictors to predict the value of this exponent

``` r
exponents <- 
  td |> 
  group_by(comid) |>
  group_modify(~ broom::tidy(lm(log(wua_per_lf) ~ 0 + log(flow_cfs), data = .x))) 
                                    
td_exp <- 
  flowlines |> st_drop_geometry() |>
  left_join(readRDS("../data/flowline_attributes.Rds"), by=join_by("comid"), relationship="one-to-one") |>
  inner_join(exponents |> select(comid, fsa_exponent=estimate)) |>
  select(comid, fsa_exponent,
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
       mtpi30_min, #vb_width_transect, vb_bf_w_ratio, 
       # channel sinuosity
       sinuosity,
       )
```

    ## Joining with `by = join_by(comid)`
