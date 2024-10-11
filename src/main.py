from rasterio.plot import show
#from gpd import *
import rasterio
import geopandas as gpd
import pandas as pd
from rasterio.mask import mask


def main():
    path = r"C:\Users\Gray\OneDrive - Texas A&M University\Fall 2024\GEOG 392 500\Project\NDVI_2023\NDVI_45_20231030_20231112_clip_20241009160501_2125730837.tif"
    cities_gdf = gpd.read_file(r"C:\Users\Gray\Downloads\tl_2017_48_place\tl_2017_48_place.shp")
    
        
    
    with rasterio.open(path) as src:
        raster_crs = src.crs
        
        if cities_gdf.crs != raster_crs:
            cities_gdf = cities_gdf.to_crs(src.crs)
            
        out_image, out_transform = mask(src, cities_gdf.geometry)
        show(out_image)

    

if __name__ == "__main__":
    main()