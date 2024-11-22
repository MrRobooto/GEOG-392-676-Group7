import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.warp import reproject, Resampling
import numpy as np
from shapely.geometry import mapping, Polygon
import pandas as pd
from scipy.ndimage import binary_fill_holes
import fiona
import os
from tqdm import tqdm

# Define base path
BASE_PATH = r"C:\Users\gmay4\Desktop\392 Project Python\Data"

# Define all file paths
PATHS = {
    'ndvi_2021': os.path.join(BASE_PATH, 'NDVI_2021.tif'),
    'ndvi_2023': os.path.join(BASE_PATH, 'NDVI_2023.tif'),
    'cropland': os.path.join(BASE_PATH, 'Cropland_Clip.tif'),
    'wildfire': os.path.join(BASE_PATH, r'Wildfire\bf8e3e3e-092e-4175-9148-49b35704d306.gdb'),
    # Output paths
    'output_cropland': os.path.join(BASE_PATH, 'processed_cropland.shp'),
    'output_fire_buffers': os.path.join(BASE_PATH, 'fire_buffers.shp'),
    'output_cropland_fire': os.path.join(BASE_PATH, 'cropland_near_fire_buffer.shp'),
    'output_ndvi': os.path.join(BASE_PATH, 'ndvi_comparison.csv'),
    'texas_boundary': os.path.join(BASE_PATH, 'texas_boundary.shp'),
}

# Register the GDB driver
fiona.supported_drivers['FileGDB'] = 'r'  # Read-only driver for geodatabases

def process_cropland(cropland_path=PATHS['cropland']):
    """
    Process cropland raster that's already clipped to contain only farmland
    """
    print(f"Processing cropland from: {cropland_path}")
    with rasterio.open(cropland_path) as src:
        print("Step 1/5: Reading cropland data...")
        cropland = src.read(1)
        transform = src.transform
        crs = src.crs
        
        # 2. Resample to 100x100m resolution
        print("Step 2/5: Resampling to 100x100m resolution...")
        scale_factor = 100/30
        new_height = int(cropland.shape[0] / scale_factor)
        new_width = int(cropland.shape[1] / scale_factor)
        
        print(f"Resampling from {cropland.shape} to {(new_height, new_width)}")
        
        out_profile = src.profile.copy()
        out_profile.update({
            'height': new_height,
            'width': new_width,
            'transform': rasterio.Affine(
                transform.a * scale_factor,
                transform.b,
                transform.c,
                transform.d,
                transform.e * scale_factor,
                transform.f
            )
        })
        
        # Show progress during resampling
        print("Step 3/5: Performing resampling...")
        resampled = np.empty((new_height, new_width), dtype=cropland.dtype)
        with tqdm(total=100, desc="Resampling") as pbar:
            reproject(
                source=cropland,
                destination=resampled,
                src_transform=transform,
                src_crs=crs,
                dst_transform=out_profile['transform'],
                dst_crs=crs,
                resampling=Resampling.mode,
                callback=lambda *args: pbar.update(1)
            )
        
        print("Step 4/5: Processing binary data and filling holes...")
        cropland_binary = (resampled > 0).astype(np.uint8)
        filled_cropland = binary_fill_holes(cropland_binary)
        
        print("Step 5/5: Converting to polygons...")
        shapes = list(rasterio.features.shapes(
            filled_cropland.astype('uint8'),
            transform=out_profile['transform']
        ))
        
        print("Creating final polygons...")
        polygons = []
        for shape in tqdm(shapes, desc="Creating polygons"):
            if shape[1] == 1:
                polygons.append(Polygon(shape[0]['coordinates'][0]))
        
        cropland_gdf = gpd.GeoDataFrame(
            {'geometry': polygons}, 
            crs=crs
        )
        
        print(f"Created {len(polygons)} cropland polygons")
        return cropland_gdf

def process_fire_data(gdb_path=PATHS['wildfire'], texas_boundary=PATHS['texas_boundary']):
    """
    Process fire polygons from geodatabase with proper CRS handling
    """
    print(f"Processing fire data from: {gdb_path}")
    try:
        print("Step 1/3: Reading geodatabase layers...")
        layers = fiona.listlayers(gdb_path)
        print(f"Available layers in geodatabase: {layers}")
        
        print("Step 2/3: Reading fire polygons...")
        fire_layer = layers[0]
        fires = gpd.read_file(gdb_path, layer=fire_layer)
        
        # Convert to a projected CRS (Texas State Plane - EPSG:3857)
        fires = fires.to_crs(epsg=3857)
        print(f"Read {len(fires)} fire polygons")
        
        if texas_boundary is not None:
            print("Step 3/3: Clipping to Texas boundary...")
            texas = gpd.read_file(texas_boundary)
            # Ensure Texas boundary is in the same CRS
            texas = texas.to_crs(epsg=3857)
            fires = gpd.overlay(fires, texas, how='intersection')
            print(f"Clipped to {len(fires)} fire polygons within Texas")
        
        # Create buffers with progress bar (now in projected CRS)
        print("Creating fire buffers...")
        fire_buffers = []
        for idx in tqdm(range(len(fires)), desc="Creating fire buffers"):
            buffer = fires.iloc[idx:idx+1].geometry.buffer(1000)  # Buffer in meters
            fire_buffers.append(buffer.iloc[0])
            
        fire_buffer_gdf = gpd.GeoDataFrame(
            geometry=fire_buffers, 
            crs=fires.crs
        )
        return fire_buffer_gdf, fires
        
    except Exception as e:
        print(f"Error processing fire data: {str(e)}")
        raise

def process_ndvi_comparison(ndvi_2021_path=PATHS['ndvi_2021'],
                          ndvi_2023_path=PATHS['ndvi_2023'],
                          cropland_buffer=None,
                          all_cropland=None):
    """
    Process and compare NDVI data for different areas
    """
    print(f"Processing NDVI data from:\n{ndvi_2021_path}\n{ndvi_2023_path}")
    
    def extract_ndvi_stats(ndvi_path, mask_geometry=None, category='all'):
        with rasterio.open(ndvi_path) as src:
            if mask_geometry is not None:
                print(f"Processing {category} for {os.path.basename(ndvi_path)}...")
                if mask_geometry.crs != src.crs:
                    print("Reprojecting geometry...")
                    mask_geometry = mask_geometry.to_crs(src.crs)
                
                geometries = [mapping(geom) for geom in tqdm(mask_geometry.geometry, desc="Preparing geometries")]
                clipped_ndvi, _ = mask(src, geometries, crop=True)
                ndvi_values = clipped_ndvi[0].flatten()
            else:
                print(f"Processing all land for {os.path.basename(ndvi_path)}...")
                ndvi_values = src.read(1).flatten()
            
            print("Filtering no-data values...")
            ndvi_values = ndvi_values[ndvi_values != src.nodata]
            
            return pd.DataFrame({
                'ndvi_value': ndvi_values,
                'category': category,
                'year': '2021' if '2021' in ndvi_path else '2023'
            })
    
    results = []
    total_steps = (2 if cropland_buffer is not None else 0) + \
                  (2 if all_cropland is not None else 0) + 2
    current_step = 1
    
    if cropland_buffer is not None:
        for ndvi_path in [ndvi_2021_path, ndvi_2023_path]:
            print(f"\nStep {current_step}/{total_steps}: Processing cropland near fire for {os.path.basename(ndvi_path)}")
            df = extract_ndvi_stats(ndvi_path, cropland_buffer, 'cropland_near_fire')
            results.append(df)
            current_step += 1
    
    if all_cropland is not None:
        for ndvi_path in [ndvi_2021_path, ndvi_2023_path]:
            print(f"\nStep {current_step}/{total_steps}: Processing all cropland for {os.path.basename(ndvi_path)}")
            df = extract_ndvi_stats(ndvi_path, all_cropland, 'all_cropland')
            results.append(df)
            current_step += 1
    
    for ndvi_path in [ndvi_2021_path, ndvi_2023_path]:
        print(f"\nStep {current_step}/{total_steps}: Processing all land for {os.path.basename(ndvi_path)}")
        df = extract_ndvi_stats(ndvi_path, category='all_land')
        results.append(df)
        current_step += 1
    
    print("\nCombining results...")
    ndvi_comparison = pd.concat(results, ignore_index=True)
    
    return ndvi_comparison

def main():
    try:
        # 1. Process cropland
        print("\nStep 1: Processing cropland...")
        cropland_gdf = process_cropland()
        # Ensure cropland is in the correct CRS
        cropland_gdf = cropland_gdf.to_crs(epsg=3857)
        cropland_gdf.to_file(PATHS['output_cropland'])
        print(f"Saved processed cropland to: {PATHS['output_cropland']}")
        
        # 2. Process fires
        print("\nStep 2: Processing fire data...")
        fire_buffer_gdf, fires = process_fire_data()
        fire_buffer_gdf.to_file(PATHS['output_fire_buffers'])
        print(f"Saved fire buffers to: {PATHS['output_fire_buffers']}")
        
        # 3. Select cropland near fires and create buffer
        print("\nStep 3: Processing cropland near fires...")
        # Ensure both geometries are in the same CRS before spatial join
        fire_buffer_gdf = fire_buffer_gdf.to_crs(cropland_gdf.crs)
        
        cropland_near_fire = gpd.sjoin(
            cropland_gdf,
            fire_buffer_gdf,
            how='inner',
            predicate='intersects'
        )
        
        # Check if we got any intersections
        if len(cropland_near_fire) == 0:
            raise ValueError("No intersections found between cropland and fire buffers. Check your data extent and CRS.")
            
        cropland_fire_buffer = cropland_near_fire.geometry.buffer(1000)
        cropland_fire_buffer_gdf = gpd.GeoDataFrame(
            geometry=cropland_fire_buffer,
            crs=cropland_near_fire.crs
        )
        cropland_fire_buffer_gdf.to_file(PATHS['output_cropland_fire'])
        print(f"Saved cropland near fire to: {PATHS['output_cropland_fire']}")
        
        # 4. Process NDVI comparison
        print("\nStep 4: Processing NDVI comparison...")
        # Ensure the buffer is in the same CRS as the NDVI data
        with rasterio.open(PATHS['ndvi_2021']) as src:
            ndvi_crs = src.crs
            cropland_fire_buffer_gdf = cropland_fire_buffer_gdf.to_crs(ndvi_crs)
            cropland_gdf = cropland_gdf.to_crs(ndvi_crs)
        
        ndvi_comparison = process_ndvi_comparison(
            cropland_buffer=cropland_fire_buffer_gdf,
            all_cropland=cropland_gdf
        )
        
        # 5. Save results
        ndvi_comparison.to_csv(PATHS['output_ndvi'])
        print(f"Saved NDVI comparison to: {PATHS['output_ndvi']}")
        
        # 6. Calculate and display statistics
        stats = ndvi_comparison.groupby(['year', 'category']).agg({
            'ndvi_value': ['count', 'mean', 'std', 'min', 'max']
        }).round(4)
        
        print("\nNDVI Statistics:")
        print(stats)
        
        # 7. Calculate NDVI change
        print("\nCalculating NDVI change between 2021 and 2023...")
        ndvi_change = ndvi_comparison.pivot_table(
            index='category',
            columns='year',
            values='ndvi_value',
            aggfunc='mean'
        )
        ndvi_change['change'] = ndvi_change['2023'] - ndvi_change['2021']
        
        print("\nNDVI Change by Category:")
        print(ndvi_change)
        
    except Exception as e:
        print(f"\nError in main execution: {str(e)}")
        raise

if __name__ == "__main__":
    main()