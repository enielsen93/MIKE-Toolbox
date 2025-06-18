# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy import interpolate
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import cKDTree
import time
import sys
sys.path.append("C:/Program Files/ArcGIS/Pro/Resources/ArcPy")
from matplotlib import cm
import mikeio
tic = time.time()
import arcpy
import os
print("%d seconds: Reading Mesh" % (time.time()-tic))

filename = r"\\files\Projects\RWA2023N006XX\RWA2023N00606\Vesterbro Torv\MIKE_URBAN\07_MESH\VBT_Mesh_2023_11_29.mesh"
#extent_shape_file = r"K:\Hydrauliske modeller\Modeller\Kapacitetsanalyse Aarhusbakken\GIS\Code_2.shp"
dfs = mikeio.dfsu.Mesh(filename)

element_area = dfs.get_element_area()
plt.figure()

plt.hist(element_area, bins=200)

plt.show()
print("BOB")