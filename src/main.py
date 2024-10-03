import arcpy
import os

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
INPUT_DB_PATH = BASE_DIR + "\data"
CSV_PATH = BASE_DIR + "\data"
OUTPUT_DB_PATH = BASE_DIR + "\scr\output"

# make input GDB path as the base working path 
arcpy.env.workspace = INPUT_DB_PATH