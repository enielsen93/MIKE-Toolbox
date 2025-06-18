# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:25:43 2021

@author: mu
"""

import mikeio
import numpy as np
filename = r"F:\Mike Urban\Kapacitetsanalyse Sydbyen\2021-12-01\Sydbyen_SA1-SA4_v3_1_sc1Y2100T100Default_2D_overland.dfsu"

dfsu = mikeio.dfsu.Dfsu(filename)

dfsu_read = dfsu.read()

courant = dfsu_read.data[6]

elements = np.where(np.nanmax(courant,axis=0)>1.2)
print(elements)

element_coordinates = dfsu.element_coordinates
for element in element_coordinates[elements,0:2][0]:
    print("%d, %d" % (element[0], element[1]))
# element_coordinates[np.where(np.nanmax(courant,axis=0)>1.5),0:2]

# import arcpy
# fc = arcpy.management.CreateFeatureclass(r"K:\Hydrauliske modeller\Papirkurv", "Courant points", "POINT")


