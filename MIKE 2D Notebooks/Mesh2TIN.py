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

mesh_files = [r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jens Juuls Vej\MIKE_URBAN\07_2D\Mesh_v2-00_Z.mesh"]

for mesh_file in mesh_files:
# mesh_file = r"C:\Papirkurv\MIWM\Herfolge_Flood_endelig_mesh_m_veje_haevet_grunde_interp.mesh"
    print(mesh_file)
    dfs = mikeio.dfsu.Mesh(mesh_file)
    node_coordinates = dfs.node_coordinates
    arcpy.env.overwriteOutput = True
    clip = True

    nodesPath = r"in_memory"
    nodesName = "nodes_Z"
    arcpy.CreateFeatureclass_management(nodesPath, nodesName, "POINT", has_z = "Enabled", spatial_reference = dfs.geometry.projection_string)
    # arcpy.AddField_management(os.path.join(nodesPath, nodesName), "Z", "DOUBLE")
    print("Creating point class")
    with arcpy.da.InsertCursor(os.path.join(nodesPath,nodesName), ["SHAPE@"]) as cursor:
        for node in node_coordinates:
            cursor.insertRow([arcpy.Point(node[0], node[1], node[2])])

    if clip:
        print("Preparing clip")
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

        print("Creating TIN")
        arcpy.ddd.CreateTin(mesh_file.replace(".mesh","TIN"), dfs.geometry.projection_string,
                            r"in_memory\nodes_Z Shape.Z Mass_Points; in_memory\ClipPolygon <None> Hard_Clip; in_memory\CutPolygon <None> Hard_Erase")
    else:
        print("Creating TIN")
        arcpy.ddd.CreateTin(mesh_file.replace(".mesh","TIN"), dfs.geometry.projection_string,
                       "%s Shape.Z Mass_Points")

