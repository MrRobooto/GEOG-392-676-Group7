import rasterio
from rasterio.plot import show
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

# Open the raster file
with rasterio.open(r"C:\Users\gmay4\Downloads\polygonclip_2005749709.tar\polygonclip_2005749709\CACHE\VEGSCAPE\NDVI_2023\NDVI_50_20231204_20231217_clip_20241009160513_935735382.tif") as src:
    #show(src)
    band1 = src.read(1)
    #print('Band1 has shape', band1.shape)
    height = band1.shape[0]
    width = band1.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rasterio.transform.xy(src.transform, rows, cols)
    lons = np.array(xs)
    lats = np.array(ys)
    #print('lons shape', lons.shape)
    #print(lons)

    ndvi = []

    for row in range(height):
        for col in range(width):
            ndvi.append(band1[row, col])

    vals = []
    for i in range(10):
        vals.append(ndvi[i])
    #coords = tuple(zip(lons, lats))
    Points = []
    for i in range(10):
        Points.append(Point(lons[i], lats[i]))
    print(Points)
    
    
    data = {
        'ndvi': vals
    }

    df = pd.DataFrame(data)
    df['geometry'] = Points
    print(df.head())
    
    
    
    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry='geometry')

    # Set the coordinate reference system (CRS) to WGS 84 (EPSG:4326)
    gdf.set_crs(epsg=5070, inplace=True)
    print(gdf)