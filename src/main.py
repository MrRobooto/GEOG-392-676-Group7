import rasterio
from rasterio.plot import show

'''
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
INPUT_DB_PATH = BASE_DIR + "\data"
CSV_PATH = BASE_DIR + "\data"
OUTPUT_DB_PATH = BASE_DIR + "\scr\output"'''

# Open and display the raster file
with rasterio.open(r"C:\Users\gmay4\Downloads\polygonclip_2005749709.tar\polygonclip_2005749709\CACHE\VEGSCAPE\NDVI_2023\NDVI_50_20231204_20231217_clip_20241009160513_935735382.tif") as src:
    show(src)
