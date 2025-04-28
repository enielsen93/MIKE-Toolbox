import sys
import time
import mikeio
import bisect
from scipy.spatial import cKDTree
from datetime import timedelta
from scipy.interpolate import RegularGridInterpolator
import arcpy
import os
import numpy as np

tic = time.time()
DFSUFile = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Hasle Torv\MIKE_URBAN\HAT_062\HAT_062_m21fm - Result Files\HAT_062_CDS_5_132_overlandBaseDefault_2D_Flood_statistics.dfsu"
DFSUField = "Maximum water depth"
RasterFileOutput = DFSUFile.replace(".dfsu",".tif")
clip_layers = None
searchDistance = 1
raster_cell_size = 0.4

def statusUpdate(a, b=2, c =2):
    print(a)

def getAvailableFilename(filepath):
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath

arcpy.env.overwriteOutput = True

clip_shapes = []
if clip_layers:
    clip_layers_list = clip_layers.split(";")
    statusUpdate("Reading Clip Layers")
    for clip_layer in clip_layers:
        clip_layer_dissolved = arcpy.Dissolve_management(clip_layer, os.path.join("in_memory", os.path.splitext(
            os.path.basename(clip_layer))[0]))[0]
        with arcpy.da.SearchCursor(clip_layer_dissolved, ["SHAPE@"]) as cursor:
            for row in cursor:
                clip_shapes.append(row[0])

statusUpdate("Reading DFSU file", tic)
dfs = mikeio.dfsu.Dfsu2DH(DFSUFile)

statusUpdate("Retrieving element coordinates from DFSU file", tic)
element_coordinates = dfs.geometry.element_coordinates

print("Reading DFSU File")
dfs_read = dfs.read(items=[i for i, a in enumerate(dfs.items) if DFSUField == a.name])
dfs_read_data = dfs_read.to_numpy()
np.nan_to_num(dfs_read_data, copy=False)

arcpy.CreateFeatureclass_management("in_memory", "ClipPolygon", "POLYGON")
boundary_xy_table = dfs.geometry.boundary_polylines[1][0].xy
polygons = arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table]))
with arcpy.da.InsertCursor("in_memory\ClipPolygon", "SHAPE@") as cursor:
    cursor.insertRow([polygons])

arcpy.CreateFeatureclass_management("in_memory", "CutPolygon", "POLYGON")
cut_polygons = dfs.geometry.boundary_polylines[3]
polygons = []
for cut_polygon in cut_polygons:
    boundary_xy_table = cut_polygon.xy
    polygons.append(arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table])))
with arcpy.da.InsertCursor("in_memory\CutPolygon", "SHAPE@") as cursor:
    for polygon in polygons:
        cursor.insertRow([polygon])

union = arcpy.Union_analysis(in_features="in_memory\CutPolygon #;in_memory\ClipPolygon #", out_feature_class=r"in_memory\Union", join_attributes="ALL", cluster_tolerance="", gaps="GAPS")[0]
total_clip_polygon = arcpy.FeatureClassToFeatureClass_conversion(union, "in_memory", "TotalClipPolygon", 'FID_CutPolygon = -1')

x_limit = [np.min(element_coordinates[:, 0]) - searchDistance,
           np.max(element_coordinates[:, 0]) + searchDistance]
y_limit = [np.min(element_coordinates[:, 1]) - searchDistance,
           np.max(element_coordinates[:, 1]) + searchDistance]
raster_xs_vector = np.arange(x_limit[0], x_limit[1], raster_cell_size)
raster_ys_vector = np.arange(y_limit[0], y_limit[1], raster_cell_size)

raster_x, raster_y = np.meshgrid(raster_xs_vector, raster_ys_vector)
raster_depth = np.zeros(raster_x.shape)
raster_x_flat = raster_x.flatten()
raster_y_flat = raster_y.flatten()
raster_depth_flat = np.zeros(raster_x.flatten().shape + tuple([1]))-100


statusUpdate("Interpolating", tic)
from scipy.interpolate import NearestNDInterpolator

interp = NearestNDInterpolator(element_coordinates[:, :2, ], dfs_read_data[0,0,:])


raster_depth_flat = interp(raster_x_flat, raster_y_flat)

total_clip_polygon = arcpy.management.CopyFeatures(total_clip_polygon, getAvailableFilename(os.path.join(r"C:\Users\ELNN\OneDrive - Ramboll\Documents\ArcGIS\scratch.gdb", "total_clip_polygon")))[0]

statusUpdate("Saving Raster", tic)
raster_depth = raster_depth_flat.reshape(raster_depth.shape[0:2] + tuple([1]))

raster_depth_compressed = raster_depth
raster = arcpy.NumPyArrayToRaster(np.flip(raster_depth_compressed[:, :, 0], axis=0),
                                  lower_left_corner=arcpy.Point(x_limit[0], y_limit[0]),
                                  x_cell_size=raster_cell_size,
                                  y_cell_size=raster_cell_size,
                                  value_to_nodata=-100)

# raster.save(RasterFileOutput)
arcpy.env.scratchWorkspace = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\ArcGIS\scratch.gdb"
raster_in_memory = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, "DFSURaster"))
raster.save(raster_in_memory)

arcpy.Clip_management(in_raster=raster_in_memory,
    rectangle="%d %d %d %d" % (x_limit[0], y_limit[1], x_limit[1], y_limit[0]),
    out_raster=RasterFileOutput,
    in_template_dataset=total_clip_polygon, nodata_value="-3.402823e+38",
    clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT")
