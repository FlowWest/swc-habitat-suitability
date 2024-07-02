Spawning Criteria (HSC) Comparison
================
[Maddee Rubenson](mailto:mrubenson@flowwest.com)
2024-07-02

- [HSC](#hsc)
  - [Read in data](#read-in-data)
  - [Depth HSI](#depth-hsi)
  - [Velocity HSI](#velocity-hsi)

## HSC

### Read in data

Data is from various sources and was compiled by Mark Gard for use in
CVPIA DSMHabitat.

``` r
hsc_raw <- readxl::read_excel(here::here('data-raw', 'source', 'hsc', 'HSC.xlsx'), sheet = "flat") |> 
  janitor::clean_names() |> 
  glimpse()
```

    ## Rows: 2,818
    ## Columns: 9
    ## $ citation           <chr> "USFWS.1997a", "USFWS.1997a", "USFWS.1997a", "USFWS…
    ## $ river              <chr> "American", "American", "American", "American", "Am…
    ## $ species            <chr> "Chinook", "Chinook", "Chinook", "Chinook", "Chinoo…
    ## $ race               <chr> "Fall", "Fall", "Fall", "Fall", "Fall", "Fall", "Fa…
    ## $ life_stage         <chr> "Spawning", "Spawning", "Spawning", "Spawning", "Sp…
    ## $ suitability_metric <chr> "Substrate", "Substrate", "Substrate", "Substrate",…
    ## $ units              <chr> "Code", "Code", "Code", "Code", "Code", "Code", "Co…
    ## $ units_si           <dbl> 0.1000000, 1.0000000, 1.2000000, 1.3000000, 1.40000…
    ## $ suitability_index  <dbl> 0.00, 0.00, 0.36, 1.00, 0.97, 0.97, 0.53, 0.28, 0.0…

``` r
spawning_hsc <- hsc_raw |> 
  filter(life_stage == "Spawning") |> glimpse()
```

    ## Rows: 1,259
    ## Columns: 9
    ## $ citation           <chr> "USFWS.1997a", "USFWS.1997a", "USFWS.1997a", "USFWS…
    ## $ river              <chr> "American", "American", "American", "American", "Am…
    ## $ species            <chr> "Chinook", "Chinook", "Chinook", "Chinook", "Chinoo…
    ## $ race               <chr> "Fall", "Fall", "Fall", "Fall", "Fall", "Fall", "Fa…
    ## $ life_stage         <chr> "Spawning", "Spawning", "Spawning", "Spawning", "Sp…
    ## $ suitability_metric <chr> "Substrate", "Substrate", "Substrate", "Substrate",…
    ## $ units              <chr> "Code", "Code", "Code", "Code", "Code", "Code", "Co…
    ## $ units_si           <dbl> 0.1000000, 1.0000000, 1.2000000, 1.3000000, 1.40000…
    ## $ suitability_index  <dbl> 0.00, 0.00, 0.36, 1.00, 0.97, 0.97, 0.53, 0.28, 0.0…

``` r
unique(spawning_hsc$suitability_metric)
```

    ## [1] "Substrate" "Depth"     "Velocity"

### Depth HSI

``` r
spawning_hsc |> 
  filter(suitability_metric == "Depth") |> 
  filter(units_si < 3,
         species == "Chinook",
         race == "Fall") |>   
  ggplot() + 
  geom_point(aes(x = units_si, y = suitability_index, color = paste0(river, " (", citation, ")"))) + 
  geom_line(aes(x = units_si, y = suitability_index, color = paste0(river, " (", citation, ")"))) + 
  facet_wrap(~race, scales = "free_x") + 
  xlab('depth (meters)') + 
  theme(legend.title=element_blank())
```

![](spawning-criteria-comparison_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Take an average across all reaches
spawning_hsc |> 
  filter(suitability_metric == "Depth") |> 
  filter(units_si < 3,
         species == "Chinook",
         race == "Fall") |>   
  mutate(units_si = round(units_si, 1)) |> 
  group_by(units_si) |> 
  summarise(mean_suit_index = mean(suitability_index)) |> 
  ggplot() + 
  geom_line(aes(x = units_si, y = mean_suit_index), size = 1) + 
  xlab('depth (meters)') + 
  ylab('mean suitability index') 
```

![](spawning-criteria-comparison_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

### Velocity HSI

``` r
spawning_hsc |> 
  filter(suitability_metric == "Velocity") |> 
  filter(units_si < 30,
         species == "Chinook",
         race == "Fall") |>   
  ggplot() + 
  geom_point(aes(x = units_si, y = suitability_index, color = paste0(river, " (", citation, ")"))) + 
  geom_line(aes(x = units_si, y = suitability_index, color = paste0(river, " (", citation, ")"))) + 
  facet_wrap(~race, scales = "free_x") + 
  xlab('velocity (meters/second)') + 
  theme(legend.title=element_blank())
```

![](spawning-criteria-comparison_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Take an average
spawning_hsc |> 
  filter(suitability_metric == "Velocity") |> 
  filter(units_si < 30,
         species == "Chinook",
         race == "Fall") |>   
  mutate(units_si = round(units_si, 1)) |> 
  group_by(units_si) |> 
  summarise(mean_suit_index = mean(suitability_index)) |> 
  ggplot() + 
  geom_line(aes(x = units_si, y = mean_suit_index), size = 1) + 
  xlab('velocity (meters/second)')  + 
  ylab('mean suitability index')
```

![](spawning-criteria-comparison_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->
