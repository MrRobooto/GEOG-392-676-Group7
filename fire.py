import arcpy

PROJ_PATH = r"C:\Users\caitl\OneDrive\Documents\Group7 Project\Group7WildfireProject\Group7WildfireProject.aprx"

# Open the project
project = arcpy.mp.ArcGISProject(PROJ_PATH)

# Debug: List all map names
print("Maps in the project:")
for map_obj in project.listMaps():
    print(map_obj.name)

# Access the map (assuming there is only 1 map object named 'Map')
mapObj = project.listMaps('Map')[0] if 'Map' in [m.name for m in project.listMaps()] else None
if not mapObj:
    print("No map named 'Map' found.")
    exit()

# Debug: List all layers and check if they are feature layers
print("Layers in the map:")
for layer in mapObj.listLayers():
    print(f"Layer: {layer.name}, IsFeatureLayer: {layer.isFeatureLayer}")

# Loop through layers in the map
for layer in mapObj.listLayers():
    if layer.isFeatureLayer:
        symbology = layer.symbology

        if hasattr(symbology, 'renderer') and layer.name == "FirePolygon":
            print(f"Updating symbology for {layer.name}")
            symbology.updateRenderer('UniqueValueRenderer')
            symbology.renderer.fields = ["Type"]
            layer.symbology = symbology
        else:
            print(f"Layer {layer.name} is not a feature layer or does not meet the symbology condition.")
    else:
        print(f"Skipping non-feature layer: {layer.name}")

# Save the updated project
output_path = r"C:\Users\caitl\OneDrive\Documents\Group7 Project\FiremapChanged\WildfireChanged\WildfireChanged.aprx"
project.saveACopy(output_path)
print(f"Project saved to: {output_path}")
