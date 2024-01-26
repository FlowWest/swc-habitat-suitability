Aquatic Species
================
Maddee Rubenson (FlowWest)
2024-01-26

``` r
selected_huc_8 <- c("18020107", "18020125")

# HUC-12 watersheds and higher level heirarchies
watersheds <- st_read("/vsizip/nhdplusv2/WBD_Subwatershed.zip") |> 
  janitor::clean_names() |>
  filter(huc_8 %in% selected_huc_8) |>
  st_transform(project_crs)
```

    ## Reading layer `WBD_Subwatershed' from data source 
    ##   `/vsizip/nhdplusv2/WBD_Subwatershed.zip' using driver `ESRI Shapefile'
    ## Simple feature collection with 4564 features and 21 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -124.5351 ymin: 32.133 xmax: -114.6198 ymax: 43.34273
    ## Geodetic CRS:  NAD83

##### Aquatic Biodiveristy Summary

The aquatic biodiversity summary combines the three measures of
biodiversity developed for ACE into a single metric: 1) aquatic native
species richness, which represents overall native diversity of all
species in the state, both common and rare; 2) aquatic rare species
richness, which represents diversity of rare species; and, 3) aquatic
irreplaceability, which is a weighted measure of rarity and endemism.

<https://data-cdfw.opendata.arcgis.com/datasets/CDFW>::aquatic-species-list-ace-ds2740-2/explore

``` r
# CDFW:
aquatic_sf <- read_sf('aquatic_species_list/Aquatic_Species_List_-_ACE_[ds2740].shp') |> 
  janitor::clean_names() |> 
  rename(huc_12 = huc12) |> glimpse()
```

    ## Rows: 4,473
    ## Columns: 6
    ## $ objectid      <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1…
    ## $ huc_12        <chr> "180201590400", "180600080503", "180400010703", "1807030…
    ## $ bio_aq_rank_s <int> 5, 5, 4, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,…
    ## $ shape_are     <dbl> 1166747899, 235531615, 250458905, 69133067, 6926633, 137…
    ## $ shape_len     <dbl> 226741.89, 110700.72, 209016.53, 167295.31, 21538.59, 84…
    ## $ geometry      <MULTIPOLYGON [°]> MULTIPOLYGON (((-121.6377 3..., MULTIPOLYGO…

``` r
aquatic_sf|>
  st_transform(project_crs) |> 
  #st_intersection(watersheds) |>
  filter(huc_12 %in% unique(watersheds$huc_12)) |>
  ggplot() +
  geom_sf(aes(fill = bio_aq_rank_s))
```

![](aquatic-species_files/figure-gfm/aqu-spec-1.png)<!-- -->

``` r
aquatic <- aquatic_sf |> st_drop_geometry()
```

##### California Freshwater Species Database V2

The California Freshwater Species Database is the first comprehensive
geospatial database of California’s freshwater species compiled and
standardized into single format from nearly 500 sources. It provides a
single source for geodata covering the plants and animals that rely on
California’s freshwater resources to survive.

<https://www.scienceforconservation.org/products/california-freshwater-species-database>

``` r
## TNC: 

richness <- read_csv('ca_freshwater_species/RichnessSummary.csv') |> 
  janitor::clean_names() |> 
  rename(huc_12 = au_id) |> 
  mutate(huc_12 = as.character(huc_12)) 
```

    ## Rows: 4450 Columns: 33
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (32): AU_ID, species, species_fish, species_crust, species_herps, specie...
    ## num  (1): OBJECTID
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
analysis_units <- read_sf('ca_freshwater_species/AnalysisUnits.shp') |> 
  janitor::clean_names() |> 
  rename(huc_12 = au_id) |> 
  left_join(richness) |> 
  glimpse()
```

    ## Joining with `by = join_by(huc_12)`

    ## Rows: 4,450
    ## Columns: 35
    ## $ huc_12                     <chr> "150301010305", "150301010307", "1503010104…
    ## $ au_name                    <chr> "Montana Wash-Colorado River", "Red Spring-…
    ## $ geometry                   <MULTIPOLYGON [m]> MULTIPOLYGON (((486938.1 -3...…
    ## $ objectid                   <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, …
    ## $ species                    <dbl> 8, 27, 8, 9, 56, 9, 13, 43, 38, 2, 2, 2, 3,…
    ## $ species_fish               <dbl> 2, 3, 1, 3, 3, 1, 2, 3, 3, 0, 0, 0, 0, 0, 0…
    ## $ species_crust              <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_herps              <dbl> 1, 1, 1, 2, 4, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2…
    ## $ species_inverts            <dbl> 0, 0, 0, 0, 0, 0, 0, 6, 2, 0, 0, 0, 0, 0, 0…
    ## $ species_mollusks           <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0…
    ## $ species_plants             <dbl> 2, 1, 0, 0, 9, 0, 1, 4, 5, 0, 0, 0, 0, 1, 1…
    ## $ species_birds              <dbl> 1, 19, 4, 2, 37, 5, 7, 27, 23, 0, 0, 0, 1, …
    ## $ species_mammals            <dbl> 2, 3, 2, 2, 3, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0…
    ## $ species_mollusks_crust     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic            <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_fish       <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_crust      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_herps      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_inverts    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_mollusks   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_plants     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_birds      <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_mammals    <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_vulnerable         <dbl> 3, 14, 2, 4, 14, 6, 7, 11, 15, 1, 1, 1, 2, …
    ## $ species_listed             <dbl> 2, 5, 1, 2, 4, 3, 4, 5, 5, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_vulnerable <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ species_endemic_listed     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ genus                      <dbl> 8, 26, 8, 8, 43, 9, 13, 36, 35, 2, 2, 2, 3,…
    ## $ family                     <dbl> 8, 19, 8, 7, 29, 8, 12, 24, 27, 2, 2, 2, 3,…
    ## $ tax_order                  <dbl> 6, 11, 7, 4, 19, 5, 8, 15, 19, 1, 1, 1, 2, …
    ## $ tax_class                  <dbl> 6, 5, 4, 4, 6, 4, 5, 7, 8, 1, 1, 1, 2, 3, 3…
    ## $ phylum                     <dbl> 2, 2, 1, 1, 2, 1, 2, 3, 4, 1, 1, 1, 1, 2, 2…
    ## $ species_current            <dbl> 2, 20, 3, 3, 41, 5, 8, 31, 32, 0, 0, 0, 0, …
    ## $ species_historical         <dbl> 2, 8, 1, 3, 15, 1, 2, 9, 3, 0, 0, 0, 0, 0, …
    ## $ species_other              <dbl> 5, 5, 4, 5, 7, 4, 5, 6, 7, 2, 2, 2, 3, 2, 3…

``` r
analysis_units |> 
  st_transform(project_crs) |> 
  filter(huc_12 %in% unique(watersheds$huc_12)) |>
  #st_intersection(watersheds) |> 
  ggplot() +
  geom_sf(aes(fill = species))
```

![](aquatic-species_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
