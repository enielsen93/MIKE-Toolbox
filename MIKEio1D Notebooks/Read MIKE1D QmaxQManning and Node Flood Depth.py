from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach
import arcpy
import os
import numpy as np
import mousereader
from mikegraph import calculate_full_flow

# res1d_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Kongelund og Marselistunnel\MIKE\KOM_Plan_027\KOM_CDS_5_sc3Base.res1d"
res1d_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\MOL\MOL_055R_opdim_LangebroBase.res1d"
MU_model = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\MOL\MOL_058R_opdim_Langebro.mdb"
# res1d_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Kongelund og Marselistunnel\MIKE\KOM_STATUS_001\KOM_STATUS_CDS20Base.res1d"
# MU_model = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Kongelund og Marselistunnel\MIKE\KOM_STATUS_001\KOM_STATUS_001.mdb"

nodes = {}
links = {}
class Node:
    def __init__(self, muid):
        self.diameter = None
        self.net_type_no = 0
        self.ground_level = 0
        self.max_level = 0
        self.id = muid

    @property
    def flood_depth(self):
        return self.max_level - self.ground_level

    @property
    def flood_volume(self):
        return self.diameter*1000*min(1, self.flood_depth)*self.flood_depth if self.diameter and self.flood_depth else 0

class Link:
    def __init__(self, MUID, diameter = 0, slope = 0, materialid = "Plastic", nettypeno = 0):
        self.MUID = MUID
        self.slope = max(slope, 1e-2) if slope else 1e-2
        self.diameter = diameter
        self.materialid = materialid
        self.net_type_no = nettypeno
        self.discharge = 0

    @property
    def QFull(self):
        if self.diameter:
            return ColebrookWhite.QFull(self.diameter, self.slope / 1e2, self.materialid)
        else:
            return 0

    @property
    def filling_degree(self):
        if self.QFull:
            return self.discharge/self.QFull
        else:
            return 0

if MU_model:
    import pyodbc
    conn = pyodbc.connect(
        r'Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;' % (MU_model))
    cursor = conn.cursor()
    cursor.execute('select MUID, Diameter, NetTypeNo from msm_Node')
    rows = cursor.fetchall()
    for row in rows:
        nodes[row[0]] = Node(row[0])
        nodes[row[0]].diameter = row[1]
        nodes[row[0]].net_type_no = row[2]

    cursor.execute('SELECT MUID, Diameter, slope_c, materialid, nettypeno FROM msm_Link')
    rows = cursor.fetchall()
    for row in rows:
        links[row[0]] = Link(row[0], row[1], row[2], row[3], row[4])


# res1d_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Kongelund og Marselistunnel\MIKE\KOM_Plan_017_sc2\KOM_Plan_017_sc2_CDS_5Base.res1d"
print(res1d_file)
node_queries = []
res1d = Res1D(res1d_file)
for node in res1d.data.Nodes:
    node_queries.append(QueryDataNode("WaterLevel", node, validate = False))

link_queries = []
for link in res1d.data.Reaches:
    link_queries.append(QueryDataReach("Discharge", link, validate = False))

nodes_DF = Res1D(res1d_file, node_queries)
links_DF = Res1D(res1d_file, link_queries)
arcpy.env.overwriteOutput = True
output_folder = r"C:\Papirkurv\Resultater"# r"C:\Users\ELNN\OneDrive - Ramboll\Documents\ArcGIS\scratch.gdb"
new_filename = os.path.basename(res1d_file).replace(".res1d",".shp")
try:
    output_filepath = arcpy.CreateFeatureclass_management(output_folder, new_filename, "POINT")[0]
    arcpy.management.AddField(output_filepath, "MUID", "TEXT")
    arcpy.management.AddField(output_filepath, "NetTypeNo", "SHORT")
    arcpy.management.AddField(output_filepath, "Flood_dep", "FLOAT", 8, 2)
    arcpy.management.AddField(output_filepath, "Flood_vol", "FLOAT", 8, 2)
    arcpy.management.AddField(output_filepath, "Result", "FLOAT", 8, 2)
except Exception as e:
    print(e)
    output_filepath = os.path.join(output_folder, new_filename)
    with arcpy.da.UpdateCursor(output_filepath, ["MUID"]) as cursor:
        for row in cursor:
            cursor.deleteRow()

with arcpy.da.InsertCursor(output_filepath, ["SHAPE@", "MUID", "Flood_dep", "Flood_vol", "NetTypeNo", "Result"]) as cursor:
    for query_node in nodes_DF.data.Nodes:
        muid = query_node.ID
        if muid not in nodes:
            nodes[muid] = Node(muid)
        node = nodes[muid]
        node.ground_level = query_node.GroundLevel
        node.max_level = np.max(nodes_DF.get_node_values(muid, "WaterLevel"))
        # print((node.max_level, node.ground_level))
        cursor.insertRow([arcpy.Point(query_node.XCoordinate, query_node.YCoordinate), muid, node.flood_depth, node.flood_volume, node.net_type_no if node.net_type_no is not None else 0, node.flood_depth])

new_filename = os.path.basename(res1d_file).replace(".res1d","QMaxQManning.shp")
try:
    output_filepath = arcpy.CreateFeatureclass_management(output_folder, new_filename, "Polyline")[0]
    arcpy.management.AddField(output_filepath, "MUID", "TEXT")
    arcpy.management.AddField(output_filepath, "NetTypeNo", "SHORT")
    arcpy.management.AddField(output_filepath, "QMax", "FLOAT", 8, 2)
    arcpy.management.AddField(output_filepath, "QManning", "FLOAT", 8, 2)
    arcpy.management.AddField(output_filepath, "Result", "FLOAT", 8, 2)
except Exception as e:
    print(e)
    output_filepath = os.path.join(output_folder, new_filename)
    with arcpy.da.UpdateCursor(output_filepath, ["MUID"]) as cursor:
        for row in cursor:
            cursor.deleteRow()

with arcpy.da.InsertCursor(output_filepath, ["SHAPE@", "MUID", "QMax", "QManning", "Result", "NetTypeNo"]) as cursor:
    for query_reach in links_DF.data.Reaches:
        points = arcpy.Array()
        for point in query_reach.LocationSpan.Coordinates:
            points.append(arcpy.Point(point.X, point.Y))
        shape = arcpy.Polyline(points)
        muid = query_reach.Name
        if muid not in links:
            links[muid] = Link(muid)
        link = links[muid]
        try:
            link.discharge = np.max(links_DF.get_reach_values(muid, 0, "Discharge"))
        except Exception as e:
            print(e)
            continue
        # print((node.max_level, node.ground_level))
        cursor.insertRow([shape, muid, link.discharge, link.QFull if link.QFull else 0, link.filling_degree, link.net_type_no if link.net_type_no is not None else 0])

print("BReak")