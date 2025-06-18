import mikeio
from pyproj import Proj, transform
import arcpy
import numpy as np
from scipy.sparse import coo_matrix, vstack

dfsu_filepath = r"C:\Offline\VOR_Status\VOR_Status_CDS10_2dCDS10.m21fm - Result Files\VOR_Status_CDS10_2D.dfsu"
dfsu_output = dfsu_filepath.replace(".dfsu","derived.dfsu")

dfs = mikeio.dfsu.Dfsu2DH(dfsu_filepath)
element_coordinates = dfs.element_coordinates
dfs_read = dfs.read(items=[i for i, a in enumerate(dfs.items) if "Surface elevation" == a.name])
dfs_read_data = dfs_read.to_numpy()
# del dfs_read
dfs_read_data = np.vstack((dfs_read_data, np.zeros((0, dfs_read_data.shape[1], dfs_read_data.shape[2]))))
dfs_read_data[-1,:,:] = dfs_read_data[-1,:,:] - element_coordinates[:, -1]
# dfs_read_data - element_coordinates[:, -1]
# dfs.write(, dfs_read_data)

items = dfs_read.items
items.append(mikeio.ItemInfo("Total water depth", mikeio.EUMType.Water_Level))
dfs.write_header(filename = dfsu_output, start_time = dfs.start_time, items = items)
# with dfs.write_header(filename = dfsu_output, start_time = dfs.start_time, items = items) as f:
#     for i in range(len())

print("BREAK")