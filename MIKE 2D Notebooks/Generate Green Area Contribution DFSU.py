2# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 11:09:22 2021

@author: eni
"""


import mikeio
import datetime
import numpy as np
import arcpy
import os
from copy import deepcopy
import matplotlib.pyplot as plt

dfsu = mikeio.dfsu.Dfsu("K:\Hydrauliske modeller\Modeller\A265_1 Kapacitetsanalyse Sydbyen\Mesh\A265_1_Mesh_GA_v2_1_cropped_interp_E13.mesh")

# extent = arcpy.Describe("K:\Hydrauliske modeller\Papirkurv\GreenContribution.shp").extent
output_folder = r"K:\Hydrauliske modeller\Modeller\A265_1 Kapacitetsanalyse Sydbyen\Mesh"

initial_loss = 30 # mm
hydrological_reduction_factor = 1
sikkerhedsfaktor = [1.11, 1.14, 1.17, 1.21, 1.23]

sikkerhedsfaktor_2050 = [0.96,	0.97,	0.99,	1,	1.02]
sikkerhedsfaktor_2100 = [1.08,	1.1,	1.12,	1.16,	1.19]


cds_files = ["K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS\CDS5aar.dfs0",
            "K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS\CDS10aar.dfs0",
            "K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS\CDS20aar.dfs0",
            "K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS\CDS50aar.dfs0",
            "K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS\CDS100aar.dfs0"]

green_area_shapefile = r"K:\Hydrauliske modeller\Modeller\A265_1 Kapacitetsanalyse Sydbyen\Mesh\ImperviousArea.shp"
green_area_shape = [row[0] for row in arcpy.da.SearchCursor(green_area_shapefile, "SHAPE@")][0]

element_coordinates = dfsu.element_coordinates[:,0:2]
try:
    arcpy.management.Delete("in_memory\Points")
except Exception as e:
    pass

points_filepath = arcpy.CreateFeatureclass_management("in_memory", "Points", "POINT")[0]
with arcpy.da.InsertCursor(points_filepath, "SHAPE@") as cursor:
    for point_shape in [arcpy.Point(element[0], element[1]) for element in element_coordinates]:
        cursor.insertRow([point_shape])
        
points_within = arcpy.management.SelectLayerByLocation(points_filepath, overlap_type = "WITHIN", select_features = green_area_shapefile,
                                                       invert_spatial_relationship = "INVERT")[0]
element_IDs = [ID-1 for ID in list(points_within.getSelectionSet())]

item = mikeio.eum.ItemInfo("Precipitation rate", itemtype = mikeio.eum.EUMType.Precipitation_Rate, unit = mikeio.eum.EUMUnit.mm_per_hour)
    
# multipoint = arcpy.Multipoint([arcpy.Point(element[0], element[1]) for element in element_coordinates])


#dx = 2
#dy = 2

#cells_horizontal = int(extent.width / dx)
#cells_vertical = int(extent.height / dy)
#
#
#cells_x = np.arange(cells_horizontal) * dx + extent.lowerLeft.X
#cells_y = np.arange(cells_vertical) * dy + extent.lowerLeft.Y
#
#cells_x_meshgrid, cells_y_meshgrid = np.meshgrid(cells_x, cells_y)
#

# for row in range(element_coordinates.shape[0]):
#     data[row] = 1 if green_area_shape.contains(arcpy.Point(element_coordinates[row,0], element_coordinates[row,1])) else 0

data = np.zeros((240,len(element_coordinates)))
fig1, ax1 = plt.subplots(figsize=(8.27,3))
fig2, ax2 = plt.subplots(figsize=(8.27,3))

for year, sikkerhedsfaktor in zip((2050, 2100), (sikkerhedsfaktor_2050, sikkerhedsfaktor_2100)):
    for cds_i, cds_file in enumerate(cds_files):
        dfs0_file = mikeio.dfs0.Dfs0(cds_file)
        cds_intensity = dfs0_file.read().data[0] / 1e3*60 # mm / min
        
        # subtract initial_loss
        cds_intensity_reduced = deepcopy(cds_intensity)
        initial_loss_left = initial_loss
        for i in range(len(cds_intensity)):
            cds_intensity_reduced[i] = max(cds_intensity[i] - initial_loss_left, 0)
            initial_loss_left = initial_loss_left - (cds_intensity[i] - cds_intensity_reduced[i])
            data[i,element_IDs] = cds_intensity_reduced[i] * hydrological_reduction_factor * sikkerhedsfaktor[cds_i] * 60
        ax = ax1 if sikkerhedsfaktor == sikkerhedsfaktor_2050 else ax2
        ax.plot(cds_intensity_reduced * sikkerhedsfaktor[cds_i] / 60 * 1e3 * hydrological_reduction_factor, 
                label = u"CDS %s år, hyd. r.factor: %1.2f, sik.faktor: %1.2f, %d mm" % 
                ('{0: <3}'.format([5,10,20,50,100][cds_i]), hydrological_reduction_factor, sikkerhedsfaktor[cds_i],
                np.sum(cds_intensity_reduced * sikkerhedsfaktor[cds_i] / 60 * 1e3 * hydrological_reduction_factor)/1e3*60))
        dfsu_filename = "Green_area_contribution_Y%d_%s" % (year, os.path.splitext(os.path.basename(cds_file))[0] )
        output_file = os.path.join(output_folder, dfsu_filename + ".dfsu")
        dfsu.write(output_file, [data], 
                          start_time = datetime.datetime.strptime("01-01-2012 00:00", "%d-%m-%Y %H:%M"),
                          dt = 60,
                          items = [item])

for ax in (ax1, ax2):
    ax.set_ylabel("Reduceret intensitet [μm/s]")
    ax.set_xlabel("Tid [min]")
    ax.set_xticks(np.arange(0, 240+1, 30))
    ax.set_xlim(0, 240)
    ax.set_ylim(0, 60)
    ax.legend(loc = "upper left", fontsize = 'x-small')
    ax.grid()
    # ax.text(237, 57, "yo mama", horizontalalignment = "right", verticalalignment = "top", bbox = dict(boxstyle='round', facecolor='white', alpha=0.5))
ax1.set_title("År 2050")
ax2.set_title("År 2100")

for year,fig in zip((2050,2100), (fig1, fig2)):
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(r'K:\Hydrauliske modeller\Modeller\A265_1 Kapacitetsanalyse Sydbyen\Mesh\Grønt bidrag År %d.jpg' % (year))


# ax1.set_ylim(0, 100)
#     dfs2_filename = os.path.splitext(os.path.basename(cds_file))[0]
#     output_file = os.path.join(output_folder, dfs2_filename + ".dfs2")
#     dfs2_file = mikeio.dfs2.Dfs2()
#     item = mikeio.eum.ItemInfo("Precipitation rate", itemtype = mikeio.eum.EUMType.Precipitation_Rate, unit = mikeio.eum.EUMUnit.mm_per_hour)
#     dfs2_file.write(filename = output_file,
#                       data = [data], 
#                       start_time = datetime.datetime.strptime("01-01-2012 00:00", "%d-%m-%Y %H:%M"),
#                       dt = 60,
#                       items = [item],
#                       dx = dx,
#                       dy = dy,
#                       coordinate = ["UTM-32", extent.lowerLeft.X, extent.lowerLeft.Y, 0])
