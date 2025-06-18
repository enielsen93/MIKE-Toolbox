import arcpy
import numpy as np
from mikeio import Grid2D
import math
import mikeio

# Input raster path
tif_path = r"C:\Users\elnn\Downloads\Arealdække_Viby_J (2).tif"
land_use_to_manning = {
    1: 30,  # Bare land
    2: 50,  # Water
    3: 40,  # Other paved
    6: 20,  # Shallow vegetation
    7: 5,  # Dense vegetation
    8: 20,  # Field
    9: 70,  # Paved road
    10: 40,  # Unpaved road
    12: 20,  # Railroad
    15: 30,  # Bare rock
    16: 50,  # Building
}

# Describe raster to get extent and projection
desc = arcpy.Describe(tif_path)
extent = desc.extent
spatial_ref = desc.spatialReference
raster = arcpy.Raster(tif_path)

# Get cell size from raster metadata
cellsize = raster.meanCellWidth

# Get number of columns and rows from raster
ncols = raster.width
nrows = raster.height

# Origin (center of lower-left cell)
xllcorner = extent.XMin
yllcorner = extent.YMin
origin = (xllcorner + cellsize / 2, yllcorner + cellsize / 2)

# Export projection
projection_wkt = spatial_ref.exportToString()

# Create Grid2D
geometry = Grid2D(
nx=ncols,
ny=nrows,
dx=cellsize,
origin=origin,
projection=projection_wkt
)

# Convert raster to NumPy array (original resolution)
array_raw = arcpy.RasterToNumPyArray(tif_path)

# Flip Y-axis if needed (arcpy gives origin top-left, mikeio expects bottom-left)
array = np.flipud(array_raw)


# Create a fast vectorized lookup table
lut = np.full(max(land_use_to_manning.keys()) + 1, np.nan)
for code, manning in land_use_to_manning.items():
    if manning is not None:
        lut[code] = manning

dtype = np.int8

# Apply LUT to array (values not in LUT will remain nan)
manning_array = lut[array].astype(dtype)

# Create DataArray with Manning values
da = mikeio.DataArray(data=manning_array,
               item=mikeio.ItemInfo("Roughness", mikeio.EUMType.Mannings_M),
               geometry=geometry,
               dims=("y","x") # No time dimension,
               )

da.to_dfs(tif_path.replace(".tif",".dfs2"), dtype = dtype)

# Optional: check shape
print(f"Array shape: {array.shape}")
print(f"Expected shape: ({nrows}, {ncols})")
print(f"Grid: {geometry}")
