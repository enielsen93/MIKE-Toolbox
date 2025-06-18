from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach, QueryDataStructure
import arcpy
import os
import numpy as np
import traceback
import cProfile
from alive_progress import alive_bar
import math
import warnings
from scipy.optimize import bisect

if len(sys.argv)>1:
    res1d_file = sys.argv[1]
else:
    extension = ""
    # MU_model = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Fredensvang\MIKE\FRE_005\FRE_005.sqlite"
    # res1d_file = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Fredensvang\MIKE\FRE_005\FRE_005_m1d - Result Files\FRE_005_CDS5_156_240_valideringBaseDefault_Network_HD.res1d"
    # MU_model = r"C:\Users\elnn\OyneDrive - Ramboll\Documents\Aarhus Vand\Hasle Torv\MIKE_URBAN\HAT_055\delete\HAT_055.sqlite"
    res1d_file = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jyllands Alle\MIKE_REGNVAND\JYL_072\JYL_072_m1d - Result Files\JYL_072_CDS_5_132_240_2080BaseResult_file.res1d"

if not 'MU_model' in locals(): # Guessing where MIKE+ database is. Write path of MIKE-model above if not correct
    model_folder = os.path.dirname(os.path.dirname(res1d_file))
    MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".sqlite"
    print(MU_model)
    if os.path.exists(MU_model):
        print("Assuming MIKE database is %s" % (MU_model))
    else:
        model_folder = os.path.dirname(res1d_file)
        MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".mdb"
        if os.path.exists(MU_model):
            print("Assuming MIKE database is %s" % (MU_model))
        else:
            raise Exception("Did not find MIKE+ Database %s." % MU_model.replace(".sqlite", ".(sqlite|mdb)"))

ms_Catchment = os.path.join(MU_model, "ms_Catchment" if ".mdb" in MU_model else "msm_Catchment")
msm_CatchCon = os.path.join(MU_model, "msm_CatchCon")

filter_to_extent = None
if "hat" in res1d_file.lower():
    filter_to_extent = [571290, 6224963, 572819, 6226104] #HAT40
# filter_to_extent = [571721, 6219599, 573205, 6220639] #SON_215
if "jjv" in res1d_file.lower():
    filter_to_extent = [570066, 6218716, 571300, 6219756] # JJV
if "son" in res1d_file.lower():
    filter_to_extent = [571590, 6219471, 572939, 6220501]

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
        self.cover_type_no = None
        self._flood_volume = None

    @property
    def flood_depth(self):
        if self.max_level and self.ground_level:
            return self.max_level - self.ground_level
        else:
            return 0

    @property
    def flood_volume(self):
        if self._flood_volume:
            return self._flood_volume
        else:
            reservoir_height = -0.25
            if self.diameter and self.flood_depth>0 and self.cover_type_no != 2:
                node_area = self.diameter**2*np.pi/4
                integral1 = (math.exp(7*min(1,(self.flood_depth-reservoir_height)))/7-math.exp(7*0)/7)*node_area
                integral2 = (max(1, (self.flood_depth-reservoir_height))-1)*node_area*1000
                return integral1+integral2
            else:
                return 0

    @flood_volume.setter
    def flood_volume(self, value):
        self._flood_volume = value

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
    import pyodbc
    if not any("Access" in item for item in pyodbc.drivers()):
        raise Exception("Error. Could not find driver for Microsoft Access! Perhaps Python is 64 bit and Access is 32 bit or vice versa? Install Microsoft Access Database Engine 2016 64 bit from Software Store.")
    with pyodbc.connect(r'Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;' % (MU_model)) as conn:
        with conn.cursor() as cursor:
            cursor.execute('select MUID, Diameter, NetTypeNo, groundlevel, criticallevel, invertlevel, covertypeno from msm_Node')
            rows = cursor.fetchall()
            for row in rows:
                nodes[row[0]] = Node(row[0])
                nodes[row[0]].diameter = row[1]
                nodes[row[0]].net_type_no = row[2]
                nodes[row[0]].ground_level = row[3]
                nodes[row[0]].critical_level = row[4]
                nodes[row[0]].invert_level = row[5]
                nodes[row[0]].cover_type_no = row[6]

            cursor.execute('select MUID, NetTypeNo from msm_Weir')
            rows = cursor.fetchall()
            for row in rows:
                reaches[row[0]] = Reach(row[0])
                reaches[row[0]].net_type_no = row[1]
                reaches[row[0]].type = "Weir"

            cursor.execute('select MUID, NetTypeNo, Diameter, uplevel, uplevel_c, dwlevel, dwlevel_c, materialid from msm_Link')
            rows = cursor.fetchall()
            for row in rows:
                reaches[row[0]] = Reach(row[0])
                reaches[row[0]].net_type_no = row[1]
                reaches[row[0]].diameter = row[2]
                reaches[row[0]].uplevel = row[3] if row[3] else row[4]
                reaches[row[0]].dwlevel = row[5] if row[5] else row[6]
                reaches[row[0]].material = row[7]

        check_catchment_connections = True

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
    with arcpy.da.SearchCursor(os.path.join(MU_model, "msm_Node"), ["MUID", "Diameter", "NetTypeNo", "GroundLevel", "CriticalLevel", "InvertLevel", "CoverTypeNo"]) as cursor:
        for row in cursor:
            nodes[row[0]] = Node(row[0])
            nodes[row[0]].diameter = row[1]
            nodes[row[0]].net_type_no = row[2]
            nodes[row[0]].ground_level = row[3]
            nodes[row[0]].critical_level = row[4]
            nodes[row[0]].invert_level = row[5]
            nodes[row[0]].cover_type_no = row[6]

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
print("Reading %s" % res1d_file)
# queries = []

extension = extension if 'extension' in locals() else ""

def read_res1d():
    res1d = Res1D(res1d_file)

    print("Reading Geometry from res1d")
    res1d_nodes = [node for node in res1d.data.Nodes]
    for reach in [r for r in res1d.data.Reaches if r.Name.replace("Weir:","") in reaches]:
        muid = reach.Name.replace("Weir:","")

        # reaches[muid].shape = arcpy.Polyline(arcpy.Array([arcpy.Point(coordinate.X, coordinate.Y) for coordinate in reach.GridPoints]))
        reaches[muid].shape = arcpy.Polyline(
            arcpy.Array([arcpy.Point(coordinate.X, coordinate.Y) for coordinate in reach.GridPoints]))
        if filter_to_extent and not (reaches[muid].shape[0][0].X > filter_to_extent[0] and reaches[muid].shape[0][0].X < filter_to_extent[2]
                and reaches[muid].shape[0][0].Y > filter_to_extent[1] and reaches[muid].shape[0][0].Y < filter_to_extent[3]):
            reaches[muid].skip = True

        reaches[muid].fromnode = res1d_nodes[reach.StartNodeIndex].ID
        reaches[muid].tonode = res1d_nodes[reach.EndNodeIndex].ID
        reaches[muid].length = reach.Length

    df = res1d

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


    # nodes_new_filename = getAvailableFilename(os.path.join(output_folder, os.path.basename(res1d_file).replace(".res1d","_nodes%s.shp" % extension)))
    #
    # links_new_filename = getAvailableFilename(os.path.join(output_folder, os.path.basename(res1d_file).replace(".res1d","_links%s.shp" % extension)))

    print("Creating Nodes and Reaches")
    # output_folder = r"C:\path\to\output"  # Replace with your path
    gdb_name = os.path.basename(res1d_file).replace(".res1d","_results%s" % extension) + ".gdb"
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
                arcpy.management.AddFields(nodes_output_filepath, [[field, "FLOAT"] for field in fields])

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
                arcpy.management.AddFields(links_output_filepath, [[field, "FLOAT"] for field in fields])

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

    timeseries = [time.timestamp() for time in df.time_index]
    print("Reading and writing Reach Results")
    with alive_bar(len(reaches), force_tty=True) as bar:
        with arcpy.da.InsertCursor(links_output_filepath, ["SHAPE@", "MUID", "Diameter", "MaxQ", "SumQ", "NetTypeNo", "EndQ", "MinQ", "MaxV", "EnergyGr", "FrictionLo", "FillDeg", "MaxTau", "Depthdiff", "Up_MaxE", "Dw_MaxE"]) as cursor:
            for muid in set(reaches.keys()):
                reach = reaches[muid]
                if not reach.skip:
                    try:
                        queries = [QueryDataReach("Discharge", muid, reach.length), QueryDataReach("FlowVelocity", muid, reach.length),
                                   QueryDataReach("WaterLevel", muid, 0), QueryDataReach("WaterLevel", muid, reach.length)]
                        query_result = res1d.read(queries)
                        reach_discharge = query_result.iloc[:,0]
                        reach.max_discharge = np.max(abs(reach_discharge))
                        reach.min_discharge = np.min(reach_discharge)
                        reach.sum_discharge = np.trapz(abs(reach_discharge), timeseries)
                        reach.end_discharge = np.round((abs(reach_discharge[-1])),2)
                        reach.max_flow_velocity = np.max(abs(query_result.iloc[:,1]))
                        reach_start_values = query_result.iloc[:,2]
                        reach_end_values = query_result.iloc[:,3]
                        reach.min_start_water_level = np.min(abs(reach_start_values))
                        reach.min_end_water_level = np.min(abs(reach_end_values))
                        reach.max_start_water_level = np.max(abs(reach_start_values))
                        reach.max_end_water_level = np.max(abs(reach_end_values))

                        # Calculate tau
                        full_discharge = reach.QFull
                        if reach.max_discharge < full_discharge:
                            water_level = bisect(bretting, 0, reach.diameter, args = (reach.max_discharge, full_discharge, reach.diameter), xtol = 0.002, maxiter = 100)
                            radius = reach.diameter/2
                            theta = 2 * math.acos((radius - water_level) / radius)
                            if water_level < radius / 2:
                                wet_perimeter = radius*theta
                                wet_area = (radius**2*(theta-math.sin(theta)))/2
                            else:
                                wet_perimeter = 2*np.pi*radius - radius*theta
                                wet_area = np.pi * radius**2 - (radius**2*(theta-math.sin(theta)))/2
                            hydraulic_radius = wet_area / wet_perimeter
                            reach.tau = 999.7 * 9.81 * reach.slope * hydraulic_radius
                        else:
                            reach.tau = 1e3

                        # Calculate Depth Difference
                        if reach.fromnode in df.nodes and reach.tonode in df.nodes:
                            reach.depth_difference = (reach.max_start_water_level - df.nodes[reach.fromnode].bottom_level) - (reach.max_end_water_level - df.nodes[
                                reach.tonode].bottom_level)

                    # print(water_level)

                    except Exception as e:
                        if muid in res1d.structures.keys():
                            try:
                                queries = [QueryDataStructure("Discharge", muid)]
                                query_result = res1d.read(queries)
                                reach_discharge = query_result.iloc[:, 0]
                                reach.max_discharge = np.max(abs(reach_discharge))
                                reach.min_discharge = np.min(reach_discharge)
                                reach.sum_discharge = np.trapz(abs(reach_discharge), timeseries)
                                reach.end_discharge = np.round(abs(reach_discharge[-1]),4)

                            except Exception as e:
                                warnings.warn("Failed to get discharge from %s" % (muid))

                    # if True:
                    if all((reach.min_start_water_level, reach.min_end_water_level, reach.max_start_water_level, reach.max_end_water_level)):
                        energy_line_gradient = reach.energy_line_gradient
                        friction_loss = reach.friction_loss
                    else:
                        energy_line_gradient = 0
                        friction_loss = 0
                    cursor.insertRow([reach.shape, muid, reach.diameter or 0, reach.max_discharge or 0, reach.sum_discharge or 0,
                                  reach.net_type_no or 0, reach.end_discharge or 0,
                                      reach.min_discharge or 0, reach.max_flow_velocity or 0,
                                      energy_line_gradient, friction_loss, reach.fill_degree or 0, reach.tau or 0, reach.depth_difference or 0, reach.max_start_water_level or 0, reach.max_end_water_level or 0])
                bar()
    res1d_quantities = res1d.quantities

    print("Reading and writing Node Results")
    with arcpy.da.InsertCursor(nodes_output_filepath, ["SHAPE@", "MUID", "Diameter", "Invert_lev", "Max_elev", "Flood_dep", "Flood_vol", "NetTypeNo", "max_hl", "max_I_V", "flow_area", "flow_diam", "end_depth", "Surcha", "SurchaBal", "MaxSurcha", "Ground_lev"]) as cursor:
        with alive_bar(len(list(df.data.Nodes)), force_tty=True) as bar:
            for query_node in df.data.Nodes:
                muid = query_node.ID

                if muid not in nodes:
                    nodes[muid] = Node(muid)

                if filter_to_extent and not (
                        query_node.XCoordinate > filter_to_extent[0] and query_node.XCoordinate < filter_to_extent[2]
                        and query_node.YCoordinate > filter_to_extent[1] and query_node.YCoordinate < filter_to_extent[3]):
                    nodes[muid].skip = True

                if not nodes[muid].skip:
                    node = nodes[muid]
                    try:
                        if not node.ground_level:
                            node.ground_level = query_node.CriticalLevel if hasattr(query_node, 'CriticalLevel') and query_node.CriticalLevel else query_node.GroundLevel
                        if not node.invert_level:
                            query_node.BottomLevel
                    except Exception as e:
                        print(e)

                    queries = [QueryDataNode("WaterLevel", muid)]
                    query_result = res1d.read(queries)
                    node.max_level = np.max(query_result.iloc[:,0])
                    node.end_depth = query_result.iloc[-1, 0] - node.invert_level
                    if "WaterVolumeAboveGround" in res1d_quantities:
                        query = QueryDataNode("WaterVolumeAboveGround", muid)
                        try:
                            query_result = res1d.read(query)
                            node.flood_volume = np.max(np.max(query_result.iloc[:, 0]))
                        except Exception as e:
                            pass

                    max_surcharge = None
                    surcharge = None
                    surcharge_balance = None
                    if "DischargeToSurface" in res1d_quantities:# and any(["DischargeToSurface" in str(dataitem.Quantity) for dataitem in res1d.nodes[muid].DataItems]):
                        query = QueryDataNode("DischargeToSurface", muid)
                        try:
                            query_result = res1d.read(query)
                            positive_surcharge = query_result.iloc[:, 0].copy()
                            positive_surcharge[positive_surcharge<0] = 0
                            surcharge = np.trapz(positive_surcharge, timeseries)
                            surcharge_balance = np.trapz(query_result.iloc[:, 0], timeseries)
                            max_surcharge = np.max(query_result.iloc[:, 0])
                        except Exception as e:
                            print(e)

                    if "DivertedRunoffToSurface" in res1d_quantities:# and any(
                        #["DivertedRunoffToSurface" in str(dataitem.Quantity) for dataitem in res1d.nodes[muid].DataItems]):
                        query = QueryDataNode("DivertedRunoffToSurface", muid)
                        try:
                            query_result = res1d.read(query)
                            diverted_runoff_to_surface = query_result.iloc[:, 0]
                            surcharge += np.trapz(diverted_runoff_to_surface, timeseries)
                            surcharge_balance += np.trapz(diverted_runoff_to_surface, timeseries)
                            if muid == "D24102R":
                                import matplotlib.pyplot as plt
                                plt.plot(diverted_runoff_to_surface)
                                print("BOB")
                            # if surcharge>0://
                            #     plt.plot(timeseries,query_result.iloc[:, 0])
                            # max_surcharge += np.trapz(query_result.iloc[:, 0], timeseries)
                        except Exception as e:
                            print(e)

                    if muid in [reach.tonode for reach in reaches.values()] and muid in [reach.fromnode for reach in reaches.values()]:
                        try:
                            water_levels = [reach.max_end_water_level for reach in reaches.values() if reach.tonode == muid and reach.type == "Link" and reach.max_end_water_level]
                            node.inlet_waterlevel = np.max(water_levels) if water_levels else 0
                            water_levels = [reach.max_start_water_level for reach in reaches.values() if reach.fromnode == muid and reach.type == "Link"]
                            node.outlet_waterlevel = np.max(water_levels) if water_levels else 0
                            node.max_headloss = node.inlet_waterlevel - node.outlet_waterlevel if all([node.inlet_waterlevel, node.outlet_waterlevel]) else 0
                            inlet_velocities = [reach.max_flow_velocity for reach in reaches.values() if reach.tonode == muid and reach.type == "Link" and reach.max_flow_velocity]
                            node.max_inlet_velocity = np.max(inlet_velocities) if inlet_velocities else 0

                        except Exception as e:
                            print(muid)
                            print(traceback.format_exc())
                            print(e)

                    cursor.insertRow([arcpy.Point(query_node.XCoordinate, query_node.YCoordinate), muid, node.diameter or 0,
                                      node.invert_level, node.max_level, node.flood_depth, node.flood_volume or 0,
                                      node.net_type_no or 0, node.max_headloss or 0,
                                      node.max_inlet_velocity or 0, node.flow_area, node.flow_area_diameter, node.end_depth or 0,
                                      surcharge or 0, surcharge_balance or 0, max_surcharge or 0, node.ground_level or 0])
                bar()

    import winsound
    winsound.Beep(300, 200)

    if len([catchment for catchment in catchments.values() if not catchment.nodeid])>0:
        print("%d catchments not connected. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]), "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid])))

    if len([catchment for catchment in catchments.values() if not catchment.nodeid_exists])>0:
        print("%d catchments connected to missing node. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]),
                                                                   "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid_exists])))

    # writing Last Modified time to the file
    timestamp_file = os.path.join(gdb_path, "last_updated.txt")
    import datetime
    with open(timestamp_file, "w") as f:
        f.write(f"Last update: {datetime.datetime.now().isoformat()}\n")

    print((nodes_new_filename, links_new_filename))

    from datetime import datetime
    now = datetime.now()
    print("Code run at %s - simulation run at %s" % (now.strftime("%H:%M"), datetime.fromtimestamp(os.path.getmtime(res1d_file)).strftime("%H:%M")))

# import cProfile
# import pstats
# import io
#
#
# def profile_res1d():
#     pr = cProfile.Profile()
#     pr.enable()
#
#     read_res1d()
#
#     pr.disable()
#     s = io.StringIO()
#     ps = pstats.Stats(pr, stream=s).sort_stats(pstats.SortKey.CUMULATIVE)
#     ps.print_stats()
#     print(s.getvalue())
#
# profile_res1d()

# Reading again loop
while True:
    read_res1d()
    user_input = input("Run again? (y/n): ").strip().lower()
    if user_input and user_input[0] == 'n':
        print("Exiting...")
        break