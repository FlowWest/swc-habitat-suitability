#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:11:25 2024

@author: maddeerubenson
"""

import geopandas as gpd
import rasterio

all_rivers_path = 'VBET/data-raw/salmonid_habitat_extents.shp'

# Read the shapefile
all_rivers = gpd.read_file(all_rivers_path)

yuba_river = all_rivers[(all_rivers['River'] == 'Yuba River') & 
                        (all_rivers['Species'] == 'Fall Run Chinook') & 
                        (all_rivers['Habitat'] == "rearing")]

yuba_river.to_file('VBET/data-raw/yuba_river_extent.shp')


# dem

dem_path = 'data-raw/dem_nhdplushr_yuba_meters.tif'

with rasterio.open(dem_path, 'r') as src:
    # Now you can work with the opened DEM dataset
    # For example, you can access its metadata or read specific raster bands
    print(src.profile)  # Print metadata
    dem_data = src.read(1) 
    
    dem_data.crs
    
# processing test

yuba_nhd = gpd.read_file('data-raw/yuba_thalweg_nhdplushr_vaa_prj.shp')
yuba_nhd.crs


yuba_nhd.crs.to_string() 

dem_data.crs.to_string()

rasterio.open(dem_path).crs.to_string()
