
# -*- coding: utf-8 -*-
import os
import arcpy
import numpy as np
import re
from arcpy import env
arcpy.env.addOutputsToMap = False
import sys
import time
import mikeio
import bisect
from scipy.spatial import cKDTree
from datetime import timedelta
from scipy.interpolate import RegularGridInterpolator
from alive_progress import alive_bar

arcpy.env.overwriteOutput = True
tic = time.time()
DFSUFile = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\05_RESULTS\03_FLOOD\Status\2024-06-10\VBT_219_mGA_CDS_20_138_FMBaseDefault_2D_overland_reduced.dfsu"
output_shape_file = r"C:\Papirkurv\VBLastTimestepPlume.shp"
max_polygons_filepath = r"C:\Papirkurv\VBLastTimestepPolygons"
minimum_water_depth = float(0.005)

dfs = mikeio.dfsu.Dfsu2DH(DFSUFile)
#DFSUField = "Total water depth"
DFSUField = "Total water depth"

# element_coordinates = dfs.element_coordinates
print(dfs.items)
dfs_read = dfs.read(items=[i for i, item in enumerate(dfs.items) if item.name in DFSUField])

dfs_read_data = dfs_read.to_numpy()
np.nan_to_num(dfs_read_data, copy=False)

elements_with_water = np.where(dfs_read_data[0, -1, :] > minimum_water_depth)[0]

element_table = dfs.element_table

plumes_class = []
plumes_nodes = []
plumes_elements = []
arcpy.SetProgressor("step", "Doing whatever...",
                    0, len(elements_with_water), 1)
arcpy.AddMessage(len(elements_with_water))
print("Classifying flooded cells")
with alive_bar(len(elements_with_water), force_tty=True) as bar:
    for element_i, element in enumerate(elements_with_water):
        arcpy.SetProgressorPosition(element_i)
        matches_found = []
        for element_node in element_table[element]:
            matches = [i for i, nodes in enumerate(plumes_nodes) if element_node in nodes]
            for match in matches:
                if match not in matches_found:
                    matches_found.append(match)

        if len(matches_found) > 0:
            plumes_elements[min(matches_found)].add(element)
            for node in element_table[element]:
                plumes_nodes[min(matches_found)].add(node)
            for match in np.flip(np.sort(matches_found))[:-1]:
                for node in plumes_nodes[match]:
                    plumes_nodes[min(matches_found)].add(node)
                for element in plumes_elements[match]:
                    plumes_elements[min(matches_found)].add(element)
                plumes_elements[min(matches_found)].add(element)
                del plumes_class[match]
                del plumes_nodes[match]
                del plumes_elements[match]
        else:
            plume_class = np.max(plumes_class) + 1 if len(plumes_class) > 0 else 0
            plumes_class.append(plume_class)
            plumes_nodes.append(set(element_table[element]))
            plumes_elements.append(set([element]))
        bar()

output_file = arcpy.CreateFeatureclass_management(r"in_memory", "boundary", "POLYGON", spatial_reference = dfs.projection_string)[0]
arcpy.AddField_management(output_file, "Class", "LONG")
arcpy.AddField_management(output_file, "Volume", "FLOAT", field_precision=15, field_scale=5)
arcpy.AddField_management(output_file, "Area", "FLOAT", field_precision=15, field_scale=5)

nodes_coordinates = dfs.node_coordinates

print("Inserting geometry from DFSU")
# arcpy.AddMessage(plumes_class)
# arcpy.AddMessage(plumes_elements[0])
with alive_bar(sum(len(sublist) for sublist in plumes_elements), force_tty=True) as bar:
    with arcpy.da.InsertCursor(output_file, ["SHAPE@", "Class", "Volume"]) as cursor:
        for plume_i in range(len(plumes_class)):
            for element in plumes_elements[plume_i]:
                node_coordinates = nodes_coordinates[element_table[element], :-1]
                triangle = arcpy.Polygon(arcpy.Array([arcpy.Point(coords[0], coords[1]) for coords in node_coordinates]))
                # arcpy.AddMessage(dfs_read_data.shape)
                # arcpy.AddMessage(dfs_read_data[0, -1, element])
                cursor.insertRow([triangle, plume_i, float(dfs_read_data[0, -1, element])])
                bar()

arcpy.CopyFeatures_management(output_file, max_polygons_filepath)

print("Calculating area of cells")
with alive_bar(int(arcpy.GetCount_management(output_file).getOutput(0)), force_tty=True) as bar:
    with arcpy.da.UpdateCursor(output_file, ["SHAPE@", "Volume", "Area"]) as cursor:
        for row in cursor:
            area = row[0].area
            row[1] = row[1] * area
            row[2] = area
            cursor.updateRow(row)
            bar()


print("Dissolving")
arcpy.Dissolve_management(output_file, output_shape_file, dissolve_field=["Class"],
                          statistics_fields=[["Volume", "SUM"], ["Area", "SUM"], ["Volume", "MAX"]])


def renameField(shape_file, old_field_name, new_field_name):
    arcpy.AddField_management(output_shape_file, new_field_name, "FLOAT", field_precision=15, field_scale=5)
    arcpy.CalculateField_management(output_shape_file, new_field_name, "!%s!" % (old_field_name), "PYTHON_9.3", "")
    arcpy.DeleteField_management(output_shape_file, old_field_name)


field_name_rename_dictionairy = {"SUM_Volume": "Volume", "MAX_Volume": "Max_depth", "SUM_Area": "Area"}

for old_field_name, new_field_name in field_name_rename_dictionairy.items():
    renameField(output_shape_file, old_field_name, new_field_name)