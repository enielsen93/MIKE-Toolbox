# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:17:40 2020

@author: eni
"""


#import arcpy
import os
import numpy as np
import mikeio
import arcpy
arcpy.env.overwriteOutput = True


arcpy.CheckOutExtension("3D")

mesh_file = r"C:\Papirkurv\VBT_219_mGA_CDS_20_138_FMBaseDefault_2D_overland_reduced.dfsu"
dfs = mikeio.dfsu.Dfsu2DH(mesh_file)

element_coordinates = dfs.element_coordinates

DFSUField = "Surface elevation"
dfs_read = dfs.read(items=[i for i, item in enumerate(dfs.items) if item.name in DFSUField])
dfs_read_data = dfs_read[0]
dfs_read_data.nanmax().to_numpy()

element_coordinates_with_water = element_coordinates[:,:]
element_coordinates_with_water[:,2] = dfs_read_data.nanmax().to_numpy()
clip = True

nodesPath = r"in_memory"
nodesName = "nodes_Z"
arcpy.CreateFeatureclass_management(nodesPath, nodesName, "POINT", has_z = "Enabled", spatial_reference = dfs.projection_string)
# arcpy.AddField_management(os.path.join(nodesPath, nodesName), "Z", "DOUBLE")

with arcpy.da.InsertCursor(os.path.join(nodesPath,nodesName), ["SHAPE@"]) as cursor:
    for node in element_coordinates:
        cursor.insertRow([arcpy.Point(node[0], node[1], node[2])])

arcpy.CreateFeatureclass_management("in_memory", "nodes_WE", "POINT", has_z = "Enabled")
# arcpy.AddField_management(os.path.join(nodesPath, nodesName), "Z", "DOUBLE")

with arcpy.da.InsertCursor(r"in_memory\nodes_WE", ["SHAPE@"]) as cursor:
    for node in element_coordinates_with_water:
        cursor.insertRow([arcpy.Point(node[0], node[1], node[2])])


if clip:
    arcpy.CreateFeatureclass_management("in_memory", "ClipPolygon", "POLYGON")
    boundary_xy_table = dfs.boundary_polylines[1][0].xy
    polygons = arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table]))
    with arcpy.da.InsertCursor("in_memory\ClipPolygon", "SHAPE@") as cursor:
        cursor.insertRow([polygons])

    arcpy.CreateFeatureclass_management("in_memory", "CutPolygon", "POLYGON")
    cut_polygons = dfs.boundary_polylines[3]
    polygons = []
    for cut_polygon in cut_polygons:
        boundary_xy_table = cut_polygon.xy
        polygons.append(arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table])))
    with arcpy.da.InsertCursor("in_memory\CutPolygon", "SHAPE@") as cursor:
        for polygon in polygons:
            cursor.insertRow([polygon])

    arcpy.ddd.CreateTin(r"C:\Papirkurv\NewTIN", dfs.projection_string,
                        r"in_memory\nodes_Z Shape.Z Mass_Points; in_memory\ClipPolygon <None> Hard_Clip; in_memory\CutPolygon <None> Hard_Erase")

    arcpy.ddd.CreateTin(r"C:\Papirkurv\NewTINWithWater", dfs.projection_string,
                        r"in_memory\nodes_WE Shape.Z Mass_Points; in_memory\ClipPolygon <None> Hard_Clip; in_memory\CutPolygon <None> Hard_Erase")
else:
    arcpy.ddd.CreateTin(r"C:\Papirkurv\NewTINWithWater", dfs.projection_string,
                   r"in_memory\nodes_WE Shape.Z Mass_Points")
    arcpy.ddd.CreateTin(r"C:\Papirkurv\NewTINWithWater", dfs.projection_string,
                        r"in_memory\nodes_WE Shape.Z Mass_Points")

