Exploratory Modeling
================
[Skyler Lewis](mailto:slewis@flowwest.com)
2024-02-12

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
    ## $ vb_width_transect         <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ vb_bf_w_ratio             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
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
  select(comid, 
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
    ## Columns: 32
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
  ggplot() + geom_point(aes(x=flow_cfs, y=wua_per_lf, color=qtile)) + facet_wrap(~varname) + scale_color_viridis_d()
```

    ## `summarise()` has grouped output by 'flow_cfs', 'varname'. You can override
    ## using the `.groups` argument.

![](model-expl_files/figure-gfm/td-gridplot-1.png)<!-- -->

## Split into training and testing datasets

``` r
set.seed(47)
td_split <- initial_split(td)
```

## Modeling

### Linear regression

``` r
lm_spec <- linear_reg() 
```

First set up recipes (sets of predictor and response variables), then
train and predict

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
  step_log(all_numeric_predictors(), -mtpi30_min)

lm_si <- workflow() |>
  add_recipe(si_rec) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_si |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic  p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>    <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.615         0.602 0.130      47.4 2.81e-98    18   353. -666. -579.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_si |> tidy()
```

    ## # A tibble: 19 × 5
    ##    term                estimate std.error statistic   p.value
    ##    <chr>                  <dbl>     <dbl>     <dbl>     <dbl>
    ##  1 (Intercept)         -5.01      2.00       -2.51  1.25e-  2
    ##  2 flow_norm_cfs       -0.140     0.00524   -26.6   4.08e-100
    ##  3 slope                0.0573    0.0410      1.40  1.62e-  1
    ##  4 sinuosity           -0.0931    0.0367     -2.54  1.14e-  2
    ##  5 erom_v_ma_fps       -0.501     0.330      -1.52  1.30e-  1
    ##  6 bf_depth_m          -0.349     0.263      -1.33  1.85e-  1
    ##  7 bf_w_d_ratio        -0.283     0.190      -1.49  1.37e-  1
    ##  8 da_k_erodibility     1.27      0.756       1.68  9.28e-  2
    ##  9 da_avg_slope         3.02      0.762       3.97  8.29e-  5
    ## 10 mean_ndvi            0.0679    0.0970      0.700 4.84e-  1
    ## 11 loc_bfi             -0.242     0.125      -1.93  5.43e-  2
    ## 12 loc_pct_clay         0.149     0.0667      2.23  2.62e-  2
    ## 13 loc_pct_sand         0.127     0.161       0.788 4.31e-  1
    ## 14 loc_permeability     0.0384    0.0261      1.47  1.42e-  1
    ## 15 loc_bedrock_depth   -0.0930    0.114      -0.818 4.14e-  1
    ## 16 loc_ppt_mean_mm      0.0846    0.253       0.335 7.38e-  1
    ## 17 mtpi30_min           0.00494   0.00250     1.98  4.85e-  2
    ## 18 nf_bfl_dry_cfs_norm  1.10      0.452       2.42  1.57e-  2
    ## 19 nf_bfl_wet_cfs_norm -1.11      0.478      -2.33  2.00e-  2

``` r
testing(td_split) |>
  mutate(hsi_frac_pred = predict(lm_si, testing(td_split))[[".pred"]]) |>
  select(hsi_frac, hsi_frac_pred, flow_norm_cfs) |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-independent model, linear regression")
```

![](model-expl_files/figure-gfm/lm-si-1.png)<!-- -->

Scale-dependent model version predicts the WUA per linear foot, using
absolute flow values as well as watershed scale parameters.

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
  step_interact(terms = ~ slope:da_area_sq_km) 

lm_sd <- workflow() |>
  add_recipe(sd_rec) |>
  add_model(lm_spec) |>
  fit(data=training(td_split))

lm_sd |> glance()
```

    ## # A tibble: 1 × 12
    ##   r.squared adj.r.squared sigma statistic p.value    df logLik   AIC   BIC
    ##       <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>
    ## 1     0.957         0.955 0.486      512.       0    23  -373.  796.  904.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lm_sd |> tidy()
```

    ## # A tibble: 24 × 5
    ##    term            estimate std.error statistic  p.value
    ##    <chr>              <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)    -1417.     512.        -2.77  5.87e- 3
    ##  2 flow_cfs           0.177    0.0196     9.03  3.14e-18
    ##  3 slope              1.16     0.316      3.68  2.55e- 4
    ##  4 da_area_sq_km     -3.34     1.34      -2.50  1.26e- 2
    ##  5 da_elev_mean     -21.8     83.4       -0.262 7.94e- 1
    ##  6 da_ppt_mean_mm   224.      37.0        6.07  2.48e- 9
    ##  7 nf_bfl_dry_cfs    42.2     10.2        4.13  4.17e- 5
    ##  8 nf_bfl_wet_cfs    48.7      6.25       7.80  3.31e-14
    ##  9 erom_q_ma_cfs     -1.19     0.556     -2.13  3.34e- 2
    ## 10 erom_v_ma_fps      0.439    2.05       0.214 8.31e- 1
    ## # ℹ 14 more rows

``` r
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(lm_sd, testing(td_split))[[".pred"]]) |>
  transmute(wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-dependent model, linear regression")
```

![](model-expl_files/figure-gfm/lm-sd-1.png)<!-- -->

### Linear regression with regularized (lasso) feature selection

``` r
#install.packages("glmnet")
lasso_spec <- linear_reg(penalty = 0.01, mixture = 1) |> set_engine("glmnet")
```

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

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.3.2

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loaded glmnet 4.1-8

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
    ## 1     0.591         0.586 0.133      113. 1.42e-101     7   336. -654. -615.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_si |> tidy()
```

    ## # A tibble: 8 × 5
    ##   term            estimate std.error statistic  p.value
    ##   <chr>              <dbl>     <dbl>     <dbl>    <dbl>
    ## 1 (Intercept)     -4.00      0.834     -4.79   2.14e- 6
    ## 2 flow_norm_cfs   -0.139     0.00529  -26.2    2.29e-98
    ## 3 sinuosity       -0.137     0.0316    -4.33   1.78e- 5
    ## 4 bf_depth_m      -0.430     0.0281   -15.3    4.40e-44
    ## 5 bf_w_d_ratio    -0.543     0.0792    -6.85   1.98e-11
    ## 6 da_avg_slope     1.99      0.299      6.66   6.54e-11
    ## 7 loc_pct_clay     0.0112    0.0219     0.513  6.08e- 1
    ## 8 loc_ppt_mean_mm  0.00233   0.0289     0.0806 9.36e- 1

``` r
testing(td_split) |>
  mutate(hsi_frac_pred = predict(lasso_si, testing(td_split))[[".pred"]]) |>
  select(hsi_frac, hsi_frac_pred, flow_norm_cfs) |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c()  + 
  ggtitle("Scale-independent model, linear regression w/ lasso feature selection")
```

![](model-expl_files/figure-gfm/lasso-si-1.png)<!-- -->

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
    ## 1     0.941         0.940 0.563      784. 9.88e-324    11  -461.  947. 1003.
    ## # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

``` r
lasso_sd |> tidy()
```

    ## # A tibble: 12 × 5
    ##    term              estimate std.error statistic  p.value
    ##    <chr>                <dbl>     <dbl>     <dbl>    <dbl>
    ##  1 (Intercept)       -73.4     14.8        -4.96  9.51e- 7
    ##  2 flow_cfs            0.172    0.0226      7.63  1.08e-13
    ##  3 da_elev_mean       -0.851    1.22       -0.698 4.86e- 1
    ##  4 da_ppt_mean_mm     12.4      1.01       12.3   1.12e-30
    ##  5 bf_depth_m          0.368    0.106       3.46  5.92e- 4
    ##  6 bf_w_d_ratio       -2.06     0.368      -5.61  3.29e- 8
    ##  7 loc_bfi             0.363    0.415       0.876 3.82e- 1
    ##  8 loc_pct_sand        0.486    0.439       1.11  2.69e- 1
    ##  9 loc_permeability   -0.126    0.0853     -1.47  1.42e- 1
    ## 10 loc_bedrock_depth  -1.08     0.336      -3.21  1.39e- 3
    ## 11 mtpi30_min          0.0522   0.00898     5.81  1.09e- 8
    ## 12 sinuosity          -0.681    0.144      -4.75  2.65e- 6

``` r
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(lasso_sd, testing(td_split))[[".pred"]]) |>
  transmute(wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c()   + 
  ggtitle("Scale-dependent model, linear regression w/ lasso feature selection")
```

![](model-expl_files/figure-gfm/lasso-sd-1.png)<!-- -->

### Random Forest Regrssion

``` r
#install.packages("ranger")
rfr_spec <- rand_forest(mode = "regression", trees = 2^12)
```

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
    ## Number of independent variables:  18 
    ## Mtry:                             4 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.00792259 
    ## R squared (OOB):                  0.8137704

``` r
testing(td_split) |>
  mutate(hsi_frac_pred = predict(rfr_si, testing(td_split))[[".pred"]]) |>
  select(hsi_frac, hsi_frac_pred, flow_norm_cfs) |>
  ggplot() + geom_point(aes(x=hsi_frac, y=hsi_frac_pred, color=flow_norm_cfs)) + 
  coord_fixed() + geom_abline() + scale_color_viridis_c() + 
  ggtitle("Scale-independent model, random forest regression")
```

![](model-expl_files/figure-gfm/rfr-si-1.png)<!-- -->

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
    ## Number of independent variables:  23 
    ## Mtry:                             4 
    ## Target node size:                 5 
    ## Variable importance mode:         none 
    ## Splitrule:                        variance 
    ## OOB prediction error (MSE):       0.1131795 
    ## R squared (OOB):                  0.9784653

``` r
testing(td_split) |>
  mutate(log_wua_per_lf_pred = predict(rfr_sd, testing(td_split))[[".pred"]]) |>
  transmute(wua_per_lf=exp(log_wua_per_lf), wua_per_lf_pred=exp(log_wua_per_lf_pred), flow_cfs) |>
  ggplot() + geom_point(aes(x=wua_per_lf, y=wua_per_lf_pred, color=flow_cfs)) + 
  scale_y_log10() + scale_x_log10() +  coord_fixed() + geom_abline() + scale_color_viridis_c()   + 
  ggtitle("Scale-dependent model, random forest regression")
```

![](model-expl_files/figure-gfm/rfr-sd-1.png)<!-- -->
