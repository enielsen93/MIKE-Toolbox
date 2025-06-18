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

mesh_file = r"\\Aarhusfil\userdata\1101189_Forsyning_og_Klimatilpasning\Bruger\TBBS\Tvis\Grid_med_terrain_groeft_udfyldt_new.mesh"
dfs = mikeio.dfsu.Mesh(mesh_file)
element_coordinates = dfs.element_coordinates

# dfs_read = dfs.read()

# elements_with_water = np.where(np.max(dfs_read.data[1],axis=0)>0)[0]
# raster_search_radius = 50

# x_limit = [np.min(element_coordinates[elements_with_water,0]) - raster_search_radius, 
#            np.max(element_coordinates[elements_with_water,0]) + raster_search_radius]
# y_limit = [np.min(element_coordinates[elements_with_water,1]) - raster_search_radius, 
#            np.max(element_coordinates[elements_with_water,1]) + raster_search_radius]


nodesPath = r"C:\Papirkurv"
nodesName = "nodes_Z.shp"
arcpy.CreateFeatureclass_management(nodesPath, nodesName, "POINT")
arcpy.AddField_management(os.path.join(nodesPath, nodesName), "Z", "DOUBLE")

with arcpy.da.InsertCursor(os.path.join(nodesPath,nodesName), ["SHAPE@XY", "Z"]) as cursor:
    for element in element_coordinates:
        cursor.insertRow([arcpy.Point(element[0],element[1]),element[2]])


# nodes = np.load(r"K:\Hydrauliske modeller\Modeller\Kapacitestanalyse Silkeborg Nord\Mesh\nodes.npy")
# idxx = np.intersect1d(np.where(nodes[:,0]<np.sort(nodes[:,0])[int(len(nodes)/10)*4]),
#                       np.where(nodes[:,0]>np.sort(nodes[:,0])[int(len(nodes)/10)*3]))
# idxy = np.intersect1d(np.where(nodes[:,1]<np.sort(nodes[:,1])[int(len(nodes)/10)*4]),
#                       np.where(nodes[:,1]>np.sort(nodes[:,1])[int(len(nodes)/10)*3]))

# idxXY = np.intersect1d(idxx,idxy)

# for idx in idxXY:
#     mindist = np.sort(np.sqrt((nodes[idxXY,0]-nodes[idx,0])**2 + (nodes[idxXY,1]-nodes[idx,1])**2))[1]
#     print mindist

#plt.scatter(nodes[idxXY,0],nodes[idxXY,1])
#nodesName = "nodes2.shp"
#arcpy.CreateFeatureclass_management("K:\Hydrauliske modeller\Modeller\Kapacitestanalyse Silkeborg Nord\Mesh", nodesName, "POINT")
#arcpy.AddField_management(os.path.join("K:\Hydrauliske modeller\Modeller\Kapacitestanalyse Silkeborg Nord\Mesh", nodesName), "mindist", "DOUBLE")
#with arcpy.da.InsertCursor(os.path.join("K:\Hydrauliske modeller\Modeller\Kapacitestanalyse Silkeborg Nord\Mesh", nodesName), ["SHAPE@XY","mindist"]) as cursor:
#    for idx in idxXY:
#        mindist = np.sort(np.sqrt((nodes[idxXY,0]-nodes[idx,0])**2 + (nodes[idxXY,1]-nodes[idx,1])**2))[1]
#        cursor.insertRow([tuple(nodes[idx][0:2]), mindist])
#        