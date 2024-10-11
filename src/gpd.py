import rasterio
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

# Takes in a GeoTIFF and converts it to a GeoDataFrame
def tiffToGDF(path):
    with rasterio.open(path) as src:
        band1 = src.read(1)
        height, width = band1.shape
        rows, cols = np.mgrid[0:height, 0:width]
        xs, ys = rasterio.transform.xy(src.transform, rows, cols)

    # Flatten the 2D arrays
    lons = xs.ravel()
    lats = ys.ravel()
    ndvi = band1.ravel()

    points = [Point(lon, lat) for lon, lat in zip(lons, lats)]

    # Create DataFrame and GeoDataFrame
    gdf = gpd.GeoDataFrame({
        'ndvi': ndvi,
        'geometry': points
    }, crs="EPSG:5070")

    return gdf