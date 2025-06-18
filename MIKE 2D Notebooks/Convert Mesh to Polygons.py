# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:25:43 2021

@author: mu
"""

import mikeio
import numpy as np
import arcpy
import os
from alive_progress import alive_bar

filenames = [
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_2023_10_17_interp.mesh"
]
filenames.reverse()

arcpy.env.overwriteOutput = True

for filename in filenames:
    print(filename)
    dfsu = mikeio.dfsu.Mesh(filename)

    node_coordinates = dfsu.node_coordinates
    element_coordinates = dfsu.element_coordinates

    element_table = dfsu.element_table

    shapefile_filepath = filename.replace(".mesh",".shp")

    arcpy.management.CreateFeatureclass(os.path.dirname(shapefile_filepath), os.path.basename(shapefile_filepath), "POLYGON", has_z = "ENABLED")
    arcpy.AddField_management(shapefile_filepath, "Z", "FLOAT")
    with arcpy.da.InsertCursor(shapefile_filepath, ["SHAPE@",'Z']) as cursor:
        with alive_bar(len(element_table), force_tty=True) as bar:
            for element_i, element in enumerate(element_table):
                cursor.insertRow([arcpy.Polygon(arcpy.Array([arcpy.Point(node[0], node[1])
                                                             for node in node_coordinates[
                                                                 np.array([element[0], element[1], element[2]])]])), element_coordinates[element_i,2]])
                bar()
            # cursor.insertRow([arcpy.Polygon(arcpy.Array([arcpy.Point(node[0], node[1])
            #                            for node in node_coordinates[np.array([element[0], element[1], element[2]])]])),
            #                   np.mean(node_coordinates[np.array([element[0], element[1], element[2]]),2])])
            # print(element)

    arcpy.AddSpatialIndex_management(shapefile_filepath)