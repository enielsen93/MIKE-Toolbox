import mikeio
import arcpy
import numpy as np

export_shapefile = arcpy.management.CreateFeatureclass(out_path=r"C:\Papirkurv", out_name="Overlapping Cells 2",
                                                       geometry_type="POINT")
arcpy.AddField_management(export_shapefile, "num", "LONG")
dfs2_filepath = r"C:\Offline\VOR_Status\VOR_Status_006_CDS20__CDS20_Coarse.dfs2"

dfs2 = mikeio.dfs2.Dfs2(dfs2_filepath)

cells = np.array([[814, 1244], [654, 1484], [655, 1484], [656, 1484], [654, 1485], [655, 1485], [656, 1485], [512, 1680], [512, 1681], [2595, 1862], [1435, 2110], [1436, 2110], [1435, 2111], [1436, 2111]])

x, y = dfs2.geometry.x, dfs2.geometry.y

with arcpy.da.InsertCursor(export_shapefile, ["SHAPE@", "num"]) as cursor:
    for cell_row in range(cells.shape[0]):
        # print(x[cells[cell_row, 0]], y[cells[cell_row, 1]])
        cursor.insertRow((arcpy.Point(x[cells[cell_row, 0]], y[cells[cell_row, 1]]), cell_row))

print("Done")
