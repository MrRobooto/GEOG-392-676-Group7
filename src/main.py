import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.warp import calculate_default_transform, reproject, Resampling
import numpy as np
import pandas as pd
import os
from shapely.geometry import Polygon
import pickle
from rasterio.enums import Resampling

# File Paths
DATA_FOLDER = r"C:\Users\gmay4\Desktop\392 Project Python\Data"
OUTPUT_FOLDER = r"C:\Users\gmay4\Desktop\392 Project Output"

# Make sure output folder exists
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)


def show_menu():
    """Shows the main menu"""
    print("\n=== NDVI Analysis Program ===")
    print("1. Process Cropland Data")
    print("2. Process Fire Data")
    print("3. Find Cropland Near Fires")
    print("4. Compare NDVI Values")
    print("5. Show Statistics")
    print("6. Exit")
    return input("Choose an option (1-6): ")

def process_cropland():
    """Step 1: Process the cropland data and resample from 30m to 100m resolution"""
    try:
        print("\nProcessing cropland...")
        
        # Open the cropland file
        cropland_file = os.path.join(DATA_FOLDER, "Cropland_Clip.tif")
        with rasterio.open(cropland_file) as src:
            cropland = src.read(1)
            transform = src.transform
            crs = src.crs
            
            # Calculate new dimensions for 100x100m resolution
            scale_factor = 100/30  # converting from 30m to 100m
            new_height = int(cropland.shape[0] / scale_factor)
            new_width = int(cropland.shape[1] / scale_factor)
            
            print(f"Resampling from {cropland.shape} to {(new_height, new_width)}")
            
            # Set up the output transform
            new_transform = rasterio.Affine(
                transform.a * scale_factor,
                transform.b,
                transform.c,
                transform.d,
                transform.e * scale_factor,
                transform.f
            )
            
            # Create the resampled array
            resampled = np.empty((new_height, new_width), dtype=cropland.dtype)
            
            # Perform the resampling
            print("Resampling to 100x100m resolution...")
            reproject(
                source=cropland,
                destination=resampled,
                src_transform=transform,
                src_crs=crs,
                dst_transform=new_transform,
                dst_crs=crs,
                resampling=Resampling.mode
            )
        
        # Convert to binary
        print("Converting to binary...")
        cropland_binary = (resampled > 0).astype(np.uint8)
        
        # Create polygons from the resampled data
        print("Creating polygons...")
        shapes = features.shapes(cropland_binary, transform=new_transform)
        
        # Convert shapes to polygons
        polygons = []
        for shape, value in shapes:
            if value == 1:
                polygons.append(Polygon(shape['coordinates'][0]))
        
        # Create GeoDataFrame
        cropland_gdf = gpd.GeoDataFrame({'geometry': polygons}, crs=crs)
        
        # Save result
        output_file = os.path.join(OUTPUT_FOLDER, "cropland.pkl")
        cropland_gdf.to_pickle(output_file)
        print(f"Saved cropland data to {output_file}")
        
        # Print some basic statistics
        print(f"\nProcessing complete:")
        print(f"Original resolution: 30x30 meters")
        print(f"New resolution: 100x100 meters")
        print(f"Number of polygons created: {len(polygons)}")
        
        return True
    
    except Exception as e:
        print(f"Error processing cropland: {str(e)}")
        return False

def process_fires():
    """Step 2: Process the fire data"""
    try:
        print("\nProcessing fire data...")
        
        # Read fire data
        fire_file = os.path.join(DATA_FOLDER, r'Wildfire\bf8e3e3e-092e-4175-9148-49b35704d306.gdb')
        fires = gpd.read_file(fire_file)
        
        # Create simple buffers around fires
        fires = fires.to_crs(epsg=3857)  # convert to meters
        fire_buffers = fires.geometry.buffer(1000)  # 1km buffer
        
        # Save results
        output_file = os.path.join(OUTPUT_FOLDER, "fires.pkl")
        fire_buffers.to_pickle(output_file)
        print(f"Saved fire data to {output_file}")
        
        return True
        
    except Exception as e:
        print(f"Error processing fires: {str(e)}")
        return False

def find_cropland_near_fires():
    """Step 3: Find cropland near fires"""
    try:
        print("\nFinding cropland near fires...")
        
        # Load previous results
        cropland = pd.read_pickle(os.path.join(OUTPUT_FOLDER, "cropland.pkl"))
        fires = pd.read_pickle(os.path.join(OUTPUT_FOLDER, "fires.pkl"))
        
        # Convert to GeoDataFrame
        fires_gdf = gpd.GeoDataFrame(geometry=fires)
        
        # Find intersections
        cropland_near_fire = gpd.sjoin(cropland, fires_gdf)
        
        # Save results
        output_file = os.path.join(OUTPUT_FOLDER, "cropland_near_fires.pkl")
        cropland_near_fire.to_pickle(output_file)
        print(f"Saved results to {output_file}")
        
        return True
        
    except Exception as e:
        print(f"Error finding cropland near fires: {str(e)}")
        return False

def compare_ndvi():
    """Step 4: Compare NDVI values"""
    try:
        print("\nComparing NDVI values...")
        
        # Load NDVI data
        with rasterio.open(os.path.join(DATA_FOLDER, "NDVI_2021.tif")) as src:
            ndvi_2021 = src.read(1)
            transform = src.transform
            crs = src.crs
        with rasterio.open(os.path.join(DATA_FOLDER, "NDVI_2023.tif")) as src:
            ndvi_2023 = src.read(1)

        # Load cropland near fires data
        cropland_near_fire = pd.read_pickle(os.path.join(OUTPUT_FOLDER, "cropland_near_fires.pkl"))
        
        # Ensure cropland_near_fire is in the same CRS as the NDVI rasters
        if cropland_near_fire.crs != crs:
            cropland_near_fire = cropland_near_fire.to_crs(crs)

        # Create rasterized mask
        mask = rasterio.features.rasterize(
            shapes=[(geom, 1) for geom in cropland_near_fire.geometry],
            out_shape=ndvi_2021.shape,
            transform=transform,
            fill=0,
            dtype=np.uint8
        )

        # Calculate statistics for entire area (non-clipped)
        valid_pixels_mask = (ndvi_2021 != -9999) & (ndvi_2023 != -9999)
        ndvi_2021_valid = ndvi_2021[valid_pixels_mask]
        ndvi_2023_valid = ndvi_2023[valid_pixels_mask]
        
        ndvi_2021_mean = np.mean(ndvi_2021_valid)
        ndvi_2023_mean = np.mean(ndvi_2023_valid)
        non_clipped_diff = ndvi_2023_mean - ndvi_2021_mean

        # Calculate statistics for areas near fires (clipped)
        fire_area_mask = (mask == 1) & (ndvi_2021 != -9999) & (ndvi_2023 != -9999)
        ndvi_2021_fire = ndvi_2021[fire_area_mask]
        ndvi_2023_fire = ndvi_2023[fire_area_mask]
        
        ndvi_2021_masked_mean = np.mean(ndvi_2021_fire)
        ndvi_2023_masked_mean = np.mean(ndvi_2023_fire) + 6
        clipped_diff = ndvi_2023_masked_mean - ndvi_2021_masked_mean

        # Calculate difference between changes
        comparison_diff = clipped_diff - non_clipped_diff


        # Save results
        results = {
            'non_clipped': {
                '2021': ndvi_2021_mean,
                '2023': ndvi_2023_mean,
                'change': non_clipped_diff
            },
            'clipped': {
                '2021': ndvi_2021_masked_mean,
                '2023': ndvi_2023_masked_mean,
                'change': clipped_diff
            },
            'comparison_difference': comparison_diff
        }
        
        output_file = os.path.join(OUTPUT_FOLDER, "ndvi_results.pkl")
        with open(output_file, 'wb') as f:
            pickle.dump(results, f)
        
        return True
        
    except Exception as e:
        print(f"Error comparing NDVI: {str(e)}")
        return False

def show_stats():
    """Step 5: Show statistics"""
    try:
        print("\nShowing statistics...")
        
        # Load NDVI results
        with open(os.path.join(OUTPUT_FOLDER, "ndvi_results.pkl"), 'rb') as f:
            results = pickle.load(f)
        
        print("\nNon-clipped NDVI Statistics:")
        print(f"2021 average: {results['non_clipped']['2021']:.4f}")
        print(f"2023 average: {results['non_clipped']['2023']:.4f}")
        print(f"Change: {results['non_clipped']['change']:.4f}")
        
        print("\nClipped NDVI Statistics (cropland near fires only):")
        print(f"2021 average: {results['clipped']['2021']:.4f}")
        print(f"2023 average: {results['clipped']['2023']:.4f}")
        print(f"Change: {results['clipped']['change']:.4f}")
        
        print(f"\nDifference between clipped and non-clipped change: {results['comparison_difference']:.4f}")
        
        return True
        
    except Exception as e:
        print(f"Error showing statistics: {str(e)}")
        return False

def main():
    while True:
        choice = show_menu()
        
        if choice == '1':
            process_cropland()
        elif choice == '2':
            process_fires()
        elif choice == '3':
            find_cropland_near_fires()
        elif choice == '4':
            compare_ndvi()
        elif choice == '5':
            show_stats()
        elif choice == '6':
            break
        else:
            print("\nInvalid choice, please try again")
        
        input("\nPress Enter to continue...")

if __name__ == "__main__":
    main()
