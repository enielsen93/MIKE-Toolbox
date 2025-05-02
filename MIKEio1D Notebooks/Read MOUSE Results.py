import arcpy
import os
import numpy as np
import traceback
import cProfile
from alive_progress import alive_bar
import math
import warnings
from scipy.optimize import bisect
import sys
import mousereader

if len(sys.argv)>1:
    res1d_file = sys.argv[1]
else:
    extension = ""
    # MU_model = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Fredensvang\MIKE\FRE_005\FRE_005.sqlite"
    # res1d_file = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Fredensvang\MIKE\FRE_005\FRE_005_m1d - Result Files\FRE_005_CDS5_156_240_valideringBaseDefault_Network_HD.res1d"
    MU_model = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Gjellerup\Mike Urban\GJE_003\GJE_003.mdb"
    res1d_file = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Gjellerup\Mike Urban\GJE_003\GJE_003_CDS_5_150Base.PRF"

if not 'MU_model' in locals(): # Guessing where MIKE+ database is. Write path of MIKE-model above if not correct
    model_folder = os.path.dirname(os.path.dirname(res1d_file))
    MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".sqlite"
    if os.path.exists(MU_model):
        print("Assuming MIKE+ database is %s" % (MU_model))
    else:
        model_folder = os.path.dirname(res1d_file)
        MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".mdb"
        if os.path.exists(MU_model):
            print("Assuming MIKE+ database is %s" % (MU_model))
        else:
            raise Exception("Did not find MIKE+ Database %s." % MU_model.replace(".sqlite", ".(sqlite|mdb)"))

ms_Catchment = os.path.join(MU_model, "ms_Catchment" if ".mdb" in MU_model else "msm_Catchment")
msm_CatchCon = os.path.join(MU_model, "msm_CatchCon")

filter_to_extent = None
# filter_to_extent = [571790, 6225063, 572819, 6226104] #HAT40
# filter_to_extent = [571721, 6219599, 573205, 6220639] #SON_215


if filter_to_extent:
    print("Skipping all reaches and nodes outside extent %s" % filter_to_extent)

print("Initializing")

nodes = {}
reaches = {}
class Node:
    def __init__(self, muid):
        self.diameter = None
        self.net_type_no = 0
        self.ground_level = 0
        self.invert_level = 0
        self.max_level = 0
        self.id = muid
        self.max_headloss = 0
        self.inlet_waterlevel = 0
        self.outlet_waterlevel = 0
        self.max_inlet_velocity = 0
        self.end_depth = 0
        self.skip = False
        self.shape = None

    @property
    def flood_depth(self):
        if self.max_level and self.ground_level:
            return self.max_level - self.ground_level
        else:
            return 0

    @property
    def flood_volume(self):
        reservoir_height = -0.25
        if self.diameter and self.flood_depth>0:
            node_area = self.diameter**2*np.pi/4
            integral1 = (math.exp(7*min(1,(self.flood_depth-reservoir_height)))/7-math.exp(7*0)/7)*node_area
            integral2 = (max(1, (self.flood_depth-reservoir_height))-1)*node_area*1000
            return integral1+integral2
        else:
            return 0

    @property
    def flow_area(self):
        return self.diameter * (self.max_level - self.invert_level) if self.max_level > 0 and self.diameter and self.invert_level else 0

    @property
    def flow_area_diameter(self):
        return np.sqrt(self.flow_area*4/np.pi) if self.flow_area > 0 else 0

class Reach:
    def __init__(self, muid):
        self.muid = muid
        self.net_type_no = 0
        self.diameter = 0
        # self.start_coordinate = None
        # self.end_coordinate = None
        self.shape = None
        self.length = None
        self.uplevel = None
        self.dwlevel = None
        self.max_discharge = None
        self.sum_discharge = None
        self.end_discharge = None
        self.min_discharge = None
        self.fromnode = None
        self.tonode = None
        self.type = "Link"
        self.max_flow_velocity = None
        self.min_start_water_level = None
        self.min_end_water_level = None
        self.max_start_water_level = None
        self.max_end_water_level = None
        self.material = None
        self.skip = False
        self.tau = None
        self.depth_difference = None

    @property
    def energy_line_gradient(self):
        return ((self.max_start_water_level - self.max_end_water_level) - (self.min_start_water_level-self.min_end_water_level)) / self.shape.length

    @property
    def friction_loss(self):
        return (self.max_start_water_level - self.max_end_water_level) - (self.min_start_water_level-self.min_end_water_level)

    @property
    def fill_degree(self):
        if all((self.max_start_water_level, self.uplevel, self.diameter)):
            if (self.max_start_water_level-self.uplevel)/self.diameter*1e2<0:
                print(self.max_start_water_level,self.uplevel,self.diameter)
            return (self.max_start_water_level-self.uplevel)/self.diameter*1e2

    @property
    def slope(self):
        if self.uplevel and self.dwlevel and self.length:
            return (self.uplevel-self.dwlevel)/self.length
        else:
            return 10e-3

    @property
    def QFull(self, resolution = 0.000001):
        if self.material[0].lower() == "p":
            k = 0.001 # Plastic roughness
        else:
            k = 0.0015 # Concrete roughness (used for all except plastic

        g = 9.82  # m2/s
        kinematic_viscosity = 0.0000013  # m2/s
        hydraulic_radius = self.diameter / 4.0 # for full pipes
        def colebrookWhite(v):
            Re = v * hydraulic_radius / kinematic_viscosity # Reynolds number
            f = 0.01 # initial guess for friction number
            # Iteratively solve the Colebrook-White equation for friction factor f
            for i in range(4):
                f = 2 / (6.4 - 2.45 * np.log(k / hydraulic_radius + 4.7 / (Re * np.sqrt(f)))) ** 2

            # Energy line gradient
            I = f * (v ** 2 / (2 * g * hydraulic_radius))

            # Return the difference between calculated and actual slope (should equal zero)
            return I-self.slope

        v = bisect(colebrookWhite, 1e-5, 500, xtol=2e-5, maxiter=50, disp=True)
        # Return the discharge
        if v:
            return v * (self.diameter / 2.0) ** 2 * np.pi
        else:
            return None

class Catchment:
    def __init__(self, muid):
        self.muid = muid
        self.nodeid = None
        self.nodeid_exists = None

print("Reading MIKE Database")
if MU_model and ".mdb" in MU_model:
    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Node"), ["MUID", "Diameter", "NetTypeNo", "groundlevel", "criticallevel", "invertlevel", "SHAPE@"]) as cursor:
        for row in cursor:
            nodes[row[0]] = Node(row[0])
            nodes[row[0]].diameter = row[1]
            nodes[row[0]].net_type_no = row[2]
            nodes[row[0]].ground_level = row[3]
            nodes[row[0]].critical_level = row[4]
            nodes[row[0]].invert_level = row[5]
            nodes[row[0]].shape = row[6]

    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Weir"),
                               ["MUID", "NetTypeNo"]) as cursor:
        for row in cursor:
            reaches[row[0]] = Reach(row[0])
            reaches[row[0]].net_type_no = row[1]
            reaches[row[0]].type = "Weir"

    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Link"),
                               ["MUID", "NetTypeNo", "Diameter", "uplevel", "uplevel_c", "dwlevel", "dwlevel_c", "materialid", "SHAPE@"]) as cursor:
        for row in cursor:
            reaches[row[0]] = Reach(row[0])
            reaches[row[0]].net_type_no = row[1]
            reaches[row[0]].diameter = row[2]
            reaches[row[0]].uplevel = row[3] if row[3] else row[4]
            reaches[row[0]].dwlevel = row[5] if row[5] else row[6]
            reaches[row[0]].material = row[7]
            reaches[row[0]].shape = row[8]

    check_catchment_connections = False

    catchments = {}
    if check_catchment_connections:
        cursor.execute('SELECT MUID FROM ms_Catchment')
        rows = cursor.fetchall()
        for row in rows:
            catchments[row[0]] = Catchment(row[0])

        cursor.execute('SELECT CatchID, NodeID FROM msm_CatchCon')
        rows = cursor.fetchall()
        for row in rows:
            catchments[row[0]].nodeid = row[1]

        for catchment in catchments.values():
            if catchment.nodeid in nodes:
                catchment.nodeid_exists = True
            else:
                catchment.nodeid_exists = False

elif MU_model and ".sqlite" in MU_model:
    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Node"), ["MUID", "Diameter", "NetTypeNo", "GroundLevel", "CriticalLevel", "InvertLevel"]) as cursor:
        for row in cursor:
            nodes[row[0]] = Node(row[0])
            nodes[row[0]].diameter = row[1]
            nodes[row[0]].net_type_no = row[2]
            nodes[row[0]].ground_level = row[3]
            nodes[row[0]].critical_level = row[4]
            nodes[row[0]].invert_level = row[5]

    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Weir"), ["MUID", "NetTypeNo"]) as cursor:
        for row in cursor:
            reaches[row[0]] = Reach(row[0])
            reaches[row[0]].net_type_no = row[1]
            reaches[row[0]].type = "Weir"

    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Link"), ["MUID", "NetTypeNo", "Diameter", "uplevel", "uplevel_c", "dwlevel", "dwlevel_c", "MaterialID"]) as cursor:
        for row in cursor:
            reaches[row[0]] = Reach(row[0])
            reaches[row[0]].net_type_no = row[1]
            reaches[row[0]].diameter = row[2]
            reaches[row[0]].uplevel = row[3] if row[3] else row[4]
            reaches[row[0]].dwlevel = row[5] if row[5] else row[6]
            reaches[row[0]].material = row[7]

    check_catchment_connections = True

    catchments = {}
    if check_catchment_connections:
        with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Catchment"), ["MUID"]) as cursor:
            for row in cursor:
                catchments[row[0]] = Catchment(row[0])

        with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_CatchCon"), ["CatchID", "NodeID"]) as cursor:
            for row in cursor:
                catchments[row[0]].nodeid = row[1]

        for catchment in catchments.values():
            if catchment.nodeid in nodes:
                catchment.nodeid_exists = True
            else:
                catchment.nodeid_exists = False

# res1d_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Kongelund og Marselistunnel\MIKE\KOM_Plan_017_sc2\KOM_Plan_017_sc2_CDS_5Base.res1d"
# dataframe = df.read()
print("Creating Shapefiles")
arcpy.env.overwriteOutput = True
output_folder = r"C:\Papirkurv\Resultater"

def getAvailableFilename(filepath):
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath

extension = extension if 'extension' in locals() else ""

def sanitizeFilename(filename):
    parts = filename.rsplit('.', 1)
    name = parts[0].replace('-', '_').replace('.', '_').replace(",","_")
    return name + '.' + parts[1] if len(parts) == 2 else name

print("Creating Nodes and Reaches")
# output_folder = r"C:\path\to\output"  # Replace with your path
gdb_name = sanitizeFilename(os.path.basename(res1d_file)).replace(".res1d","_results%s" % extension).replace(".PRF","_results%s" % extension) + ".gdb"
nodes_new_filename = "Nodes"
links_new_filename = "Reaches"

gdb_path = os.path.join(output_folder, gdb_name)
nodes_output_filepath = os.path.join(gdb_path, nodes_new_filename)
links_output_filepath = os.path.join(gdb_path, links_new_filename)

while True:
    try:
        if not arcpy.Exists(gdb_path):
            arcpy.CreateFileGDB_management(output_folder, gdb_name)

        fields = ["Diameter", "Ground_lev", "Invert_lev", "Max_elev", "Flood_dep", "Flood_vol", "max_hl",
                  "max_I_V", "flow_area", "flow_diam", "end_depth", "Surcha", "SurchaBal", "MaxSurcha"]
        if arcpy.Exists(nodes_output_filepath):
            arcpy.DeleteFeatures_management(nodes_output_filepath)
            existing_fields = [f.name for f in arcpy.ListFields(nodes_output_filepath)]
            for field in fields:
                if field not in existing_fields:
                    arcpy.management.AddField(nodes_output_filepath, field, "FLOAT")

        else:
            nodes_output_filepath = arcpy.CreateFeatureclass_management(gdb_path, nodes_new_filename, "POINT")[0]
            arcpy.management.AddField(nodes_output_filepath, "MUID", "TEXT")
            arcpy.management.AddField(nodes_output_filepath, "NetTypeNo", "SHORT")

            # for field in ["Diameter", "Ground_lev", "Invert_lev", "Max_elev", "Flood_dep", "Flood_vol", "max_hl", "max_I_V", "flow_area", "flow_diam", "end_depth", "Surcha", "SurchaBal", "MaxSurcha"]:
            # arcpy.management.AddField(nodes_output_filepath, field, "FLOAT", 8, 2)
            for field in fields:
                arcpy.management.AddField(nodes_output_filepath, field, "FLOAT")

        fields = ["Diameter", "MaxQ", "SumQ", "Up_MaxE", "Dw_MaxE", "EndQ", "MinQ", "MaxV", "FillDeg", "EnergyGr", "FrictionLo",
                          "MaxTau", "Depthdiff"]
        if arcpy.Exists(links_output_filepath):
            arcpy.DeleteFeatures_management(links_output_filepath)

            existing_fields = [f.name for f in arcpy.ListFields(links_output_filepath)]
            for field in fields:
                if field not in existing_fields:
                    arcpy.management.AddField(links_output_filepath, field, "FLOAT")

        else:
            links_output_filepath = arcpy.CreateFeatureclass_management(gdb_path, links_new_filename, "POLYLINE")[0]
            arcpy.management.AddField(links_output_filepath, "MUID", "TEXT")
            arcpy.management.AddField(links_output_filepath, "NetTypeNo", "SHORT")
            for field in fields:
                arcpy.management.AddField(links_output_filepath, field, "FLOAT")

        break
    except arcpy.ExecuteError as e:
        if "ERROR 000464: Cannot get exclusive schema lock" in str(e):
            input("The file %s is locked. Press enter to retry, after unlocking the file..." % fc_path)
        else:
            raise


def bretting(y, max_discharge, full_discharge, di):
    q_div_qf = 0.46 - 0.5 * math.cos(np.pi * y / di) + 0.04 * math.cos(2 * np.pi * y / di)
    # return q_div_qf
    return q_div_qf - max_discharge / full_discharge


print("Reading and writing Reach Results")
mouse_result = mousereader.MouseResult(res1d_file, ["ALL"], "Link_Q")
with alive_bar(len(reaches), force_tty=True) as bar:
    with arcpy.da.InsertCursor(links_output_filepath, ["SHAPE@", "MUID", "Diameter", "MaxQ", "SumQ", "NetTypeNo", "EndQ", "MinQ", "MaxV", "EnergyGr", "FrictionLo", "FillDeg", "MaxTau", "Depthdiff"]) as cursor:
        for muid in set(reaches.keys()):
            reach = reaches[muid]
            try:
                discharge = mouse_result.query(muid)
                reach.max_discharge = np.max(abs(discharge))
                reach.sum_discharge = np.sum(discharge)*60
                reach.end_discharge = discharge[-1]
                reach.min_discharge = np.min(discharge)
                reach.max_flow_velocity = 0
                energy_line_gradient = 0
                friction_loss = 0
                reach.fill_degree = 0
                reach.tau = 0
                reach.depth_difference = 0
                cursor.insertRow([reach.shape, muid, reach.diameter or 0, reach.max_discharge or 0, reach.sum_discharge or 0,
                              reach.net_type_no or 0, reach.end_discharge or 0,
                                  reach.min_discharge or 0, reach.max_flow_velocity or 0,
                                  energy_line_gradient, friction_loss, reach.fill_degree or 0, reach.tau or 0, reach.depth_difference or 0])
            except Exception as e:
                print(e)
            bar()


print("Reading and writing Node Results")
mouse_result = mousereader.MouseResult(res1d_file, ["ALL"], "Node_WL")
with arcpy.da.InsertCursor(nodes_output_filepath, ["SHAPE@", "MUID", "Diameter", "Invert_lev", "Max_elev", "Flood_dep", "Flood_vol", "NetTypeNo", "max_hl", "max_I_V", "flow_area", "flow_diam", "end_depth", "Surcha", "SurchaBal", "MaxSurcha", "Ground_lev"]) as cursor:
    with alive_bar(len(nodes), force_tty=True) as bar:
        for muid in set(nodes.keys()):
            node = nodes[muid]
            try:
                water_level = mouse_result.query(muid)
                node.max_level = np.max(water_level)
                surcharge = 0
                surcharge_balance = 0
                max_surcharge = 0
                cursor.insertRow([node.shape, muid, node.diameter or 0,
                                      node.invert_level, node.max_level, node.flood_depth, node.flood_volume or 0,
                                      node.net_type_no or 0, node.max_headloss or 0,
                                      node.max_inlet_velocity or 0, node.flow_area, node.flow_area_diameter, node.end_depth or 0,
                                      surcharge or 0, surcharge_balance or 0, max_surcharge or 0, node.ground_level or 0])
            except Exception as e:
                print(e)
            bar()
#
#
# import winsound
# winsound.Beep(300, 200)
#
# if len([catchment for catchment in catchments.values() if not catchment.nodeid])>0:
#     print("%d catchments not connected. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]), "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid])))
#
# if len([catchment for catchment in catchments.values() if not catchment.nodeid_exists])>0:
#     print("%d catchments connected to missing node. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]),
#                                                                "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid_exists])))
#
# print((nodes_new_filename, links_new_filename))

from datetime import datetime
now = datetime.now()
print("Code run at %s - simulation run at %s" % (now.strftime("%H:%M"), datetime.fromtimestamp(os.path.getmtime(res1d_file)).strftime("%H:%M")))