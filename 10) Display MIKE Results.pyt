# -*- coding: utf-8 -*-
import os
import arcpy
import numpy as np
import re
from arcpy import env
arcpy.env.addOutputsToMap = False
import codecs
import sys
import mousereader
import ColebrookWhite
from mikegraph import PipeNetwork
import codecs

from subprocess import call
from shutil import copyfile
import traceback
import scipy

from matplotlib.ticker import FuncFormatter, MaxNLocator



def readRes1D(res1d_file, MU_model = None, gdb_path = None, filter_to_extent = None, date_filter = None):
    from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach, QueryDataStructure

    if MU_model:
        ms_Catchment = os.path.join(MU_model, "ms_Catchment" if ".mdb" in MU_model else "msm_Catchment")
        msm_CatchCon = os.path.join(MU_model, "msm_CatchCon")

    if filter_to_extent:
        arcpy.AddMessage("Skipping all reaches and nodes outside extent %s" % filter_to_extent)

    arcpy.AddMessage("Initializing")

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
                    arcpy.AddMessage((self.max_start_water_level,self.uplevel,self.diameter))
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

    arcpy.AddMessage("Reading MIKE Database")
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
        if check_catchment_connections and MU_model:
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
    arcpy.AddMessage("Reading %s" % res1d_file)
    # queries = []

    extension = extension if 'extension' in locals() else ""

    if date_filter:
        def convertDate(date_str):
            date_str = date_str.strip()
            print(date_str)

            # Special case: only year
            if date_str.isdigit() and len(date_str) == 4:
                return datetime.datetime(int(date_str), 1, 1)

            formats = [
                "%d-%m-%Y",
                "%d/%m/%Y",
                "%d.%m.%Y",
                "%d-%m-%y",
                "%d/%m/%y",
                "%d.%m.%y",
                "%d-%m",  # Assume current year
                "%d/%m",
                "%d.%m",
                "%Y-%m-%d"
            ]

            for fmt in formats:
                try:
                    parsed = datetime.datetime.strptime(date_str, fmt)
                    # If no year was provided (default is 1900), use current year
                    if parsed.year == 1900:
                        parsed = parsed.replace(year=datetime.now().year)
                    return parsed
                except ValueError:
                    continue

            raise ValueError(f"Failed to interpret date: {date_str}")

        time_filter = convertDate(date_filter.split(" - ")[0]), convertDate(date_filter.split(" - ")[1])
    res1d = Res1D(res1d_file)#, time = time_filter if date_filter else None)

    arcpy.AddMessage("Reading Geometry from res1d")
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
    arcpy.AddMessage("Creating Shapefiles")
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

    arcpy.AddMessage("Creating Nodes and Reaches")
    # output_folder = r"C:\path\to\output"  # Replace with your path
    # gdb_name = os.path.basename(res1d_file).replace(".res1d","_results%s" % extension) + ".gdb"
    nodes_new_filename = "%s_Nodes" % os.path.basename(res1d_file).replace(".res1d","").replace("Base","").replace("Result_file","")
    links_new_filename = "%s_Reaches" % os.path.basename(res1d_file).replace(".res1d","").replace("Base","").replace("Result_file","")

    # gdb_path = arcpy.env.ScratchGDB
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

            fields = ["Diameter", "MaxQ", "SumQ", "EndQ", "MinQ", "MaxV", "FillDeg", "EnergyGr", "FrictionLo",
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
                arcpy.management.AddField(links_output_filepath, "FromNode", "TEXT")
                arcpy.management.AddField(links_output_filepath, "ToNode", "TEXT")
                arcpy.management.AddFields(links_output_filepath, [[field, "FLOAT"] for field in fields])

            # Adding metadata
            metadata_tablename = "Metadata"
            if arcpy.Exists(os.path.join(gdb_path, metadata_tablename)):
                try:
                    arcpy.management.DeleteRows(os.path.join(gdb_path, metadata_tablename))
                except Exception as e:
                    arcpy.AddMessage(e)
                    arcpy.AddMessage(os.path.join(gdb_path, metadata_tablename))
                    pass
            else:
                metadata_filepath = arcpy.management.CreateTable(gdb_path, metadata_tablename)[0]
                arcpy.management.AddField(metadata_filepath, "res1d_path", "TEXT", field_length=500)
                arcpy.management.AddField(metadata_filepath, "simulation_date", "DATE")
                arcpy.management.AddField(metadata_filepath, "result_analysis_date", "DATE")

            with arcpy.da.InsertCursor(os.path.join(gdb_path, metadata_tablename), ["res1d_path", "simulation_date", "result_analysis_date"]) as cursor:
                cursor.insertRow([res1d_file, datetime.datetime.fromtimestamp(os.path.getmtime(res1d_file)), datetime.datetime.now()])


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

    for reach in res1d.reaches.values():
        if reach.group == "Reach":
            reach_quantities = reach.quantities

    timeseries = [time.timestamp() for time in df.time_index]
    arcpy.AddMessage("Reading and writing Reach Results")
    with arcpy.da.InsertCursor(links_output_filepath, ["SHAPE@", "MUID", "Diameter", "FromNode", "ToNode", "MaxQ", "SumQ", "NetTypeNo", "EndQ", "MinQ", "MaxV", "EnergyGr", "FrictionLo", "FillDeg", "MaxTau", "Depthdiff"]) as cursor:
        for muid in set(reaches.keys()):
            reach = reaches[muid]
            if not reach.skip:
                if muid in res1d.reaches.keys():
                    queries = []
                    query_labels = []
                    for quantity in ["Discharge", "FlowVelocity", "WaterLevel"]:
                        if quantity in reach_quantities:
                            if quantity == "WaterLevel":
                                queries.append(QueryDataReach(quantity, muid, 0))
                                query_labels.append("WaterLevel_start")
                                queries.append(QueryDataReach(quantity, muid, reach.length))
                                query_labels.append("WaterLevel_end")
                            else:
                                queries.append(QueryDataReach(quantity, muid, reach.length))
                                query_labels.append(quantity)

                    query_result = res1d.read(queries)
                    query_result.columns = query_labels
                    # reach_discharge = query_result.iloc[:,0]
                    if "Discharge" in reach_quantities:
                        reach.max_discharge = np.max(abs(query_result["Discharge"]))
                        reach.min_discharge = np.min(query_result["Discharge"])
                        reach.sum_discharge = np.trapz(abs(query_result["Discharge"]), timeseries)
                        reach.end_discharge = np.round((abs(query_result["Discharge"][-1])), 2)
                    if "FlowVelocity" in reach_quantities:
                        reach.max_flow_velocity = np.max(abs(query_result["FlowVelocity"]))
                    if "WaterLevel" in reach_quantities:
                        reach_start_values = query_result["WaterLevel_start"]
                        reach_end_values = query_result["WaterLevel_end"]
                        reach.min_start_water_level = np.min(abs(reach_start_values))
                        reach.min_end_water_level = np.min(abs(reach_end_values))
                        reach.max_start_water_level = np.max(abs(reach_start_values))
                        reach.max_end_water_level = np.max(abs(reach_end_values))

                    # Calculate tau
                    try:
                        full_discharge = reach.QFull
                        if reach.max_discharge < full_discharge:
                            water_level = bisect(bretting, 0, reach.diameter,
                                                 args=(reach.max_discharge, full_discharge, reach.diameter), xtol=0.002,
                                                 maxiter=100)
                            radius = reach.diameter / 2
                            theta = 2 * math.acos((radius - water_level) / radius)
                            if water_level < radius / 2:
                                wet_perimeter = radius * theta
                                wet_area = (radius ** 2 * (theta - math.sin(theta))) / 2
                            else:
                                wet_perimeter = 2 * np.pi * radius - radius * theta
                                wet_area = np.pi * radius ** 2 - (radius ** 2 * (theta - math.sin(theta))) / 2
                            hydraulic_radius = wet_area / wet_perimeter
                            reach.tau = 999.7 * 9.81 * reach.slope * hydraulic_radius
                        else:
                            reach.tau = 1e3
                    except Exception as e:
                        pass

                    # Calculate Depth Difference
                    if reach.fromnode in df.nodes and reach.tonode in df.nodes:
                        reach.depth_difference = (reach.max_start_water_level - df.nodes[
                            reach.fromnode].bottom_level) - (reach.max_end_water_level - df.nodes[
                            reach.tonode].bottom_level)
                elif muid in res1d.structures.keys():
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
                cursor.insertRow([reach.shape, muid, reach.diameter if reach.diameter and reach.diameter>0 else 0, reach.fromnode, reach.tonode, reach.max_discharge or 0, reach.sum_discharge or 0,
                              reach.net_type_no or 0, reach.end_discharge or 0,
                                  reach.min_discharge or 0, reach.max_flow_velocity or 0,
                                  energy_line_gradient, friction_loss, reach.fill_degree or 0, reach.tau or 0, reach.depth_difference or 0])
    res1d_quantities = res1d.quantities

    arcpy.AddMessage("Reading and writing Node Results")
    with arcpy.da.InsertCursor(nodes_output_filepath, ["SHAPE@", "MUID", "Diameter", "Invert_lev", "Max_elev", "Flood_dep", "Flood_vol", "NetTypeNo", "max_hl", "max_I_V", "flow_area", "flow_diam", "end_depth", "Surcha", "SurchaBal", "MaxSurcha", "Ground_lev"]) as cursor:
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
                    arcpy.AddMessage(e)

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
                        arcpy.AddMessage(e)

                if "DivertedRunoffToSurface" in res1d_quantities:# and any(
                    #["DivertedRunoffToSurface" in str(dataitem.Quantity) for dataitem in res1d.nodes[muid].DataItems]):
                    query = QueryDataNode("DivertedRunoffToSurface", muid)
                    try:
                        query_result = res1d.read(query)
                        diverted_runoff_to_surface = query_result.iloc[:, 0]
                        surcharge += np.trapz(diverted_runoff_to_surface, timeseries)
                        surcharge_balance += np.trapz(diverted_runoff_to_surface, timeseries)
                        # if surcharge>0://
                        #     plt.plot(timeseries,query_result.iloc[:, 0])
                        # max_surcharge += np.trapz(query_result.iloc[:, 0], timeseries)
                    except Exception as e:
                        arcpy.AddMessage(e)

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
                        arcpy.AddMessage(muid)
                        arcpy.AddMessage(traceback.format_exc())
                        arcpy.AddMessage(e)

                cursor.insertRow([arcpy.Point(query_node.XCoordinate, query_node.YCoordinate), muid, node.diameter or 0,
                                  node.invert_level, node.max_level, node.flood_depth, node.flood_volume or 0,
                                  node.net_type_no or 0, node.max_headloss or 0,
                                  node.max_inlet_velocity or 0, node.flow_area, node.flow_area_diameter, node.end_depth or 0,
                                  surcharge or 0, surcharge_balance or 0, max_surcharge or 0, node.ground_level or 0])

    if MU_model and len([catchment for catchment in catchments.values() if not catchment.nodeid])>0:
        arcpy.AddMessage("%d catchments not connected. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]), "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid])))

    if MU_model and len([catchment for catchment in catchments.values() if not catchment.nodeid_exists])>0:
        arcpy.AddMessage("%d catchments connected to missing node. ('%s')" % (len([catchment for catchment in catchments.values() if not catchment.nodeid_exists]),
                                                                   "', '".join([catchment.muid for catchment in catchments.values() if not catchment.nodeid_exists])))

    now = datetime.datetime.now()
    arcpy.AddMessage("Code run at %s - simulation run at %s" % (now.strftime("%H:%M"), datetime.datetime.fromtimestamp(os.path.getmtime(res1d_file)).strftime("%H:%M")))
    
    return nodes_output_filepath, links_output_filepath

def m11extrapath():
    m11extraPath = r"C:\Program Files (x86)\DHI\2016\bin\m11extra.exe"
    i = 2030
    while not os.path.exists(m11extraPath.replace("2016",str(i))) and i < 20:
        i -= 1
    m11extraPath = m11extraPath.replace("2016",str(i))
    return m11extraPath

if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument
# arcpy.env.workspace = arcpy.env.scratchGDB

def getAvailableFilename(filepath, parent = None):
    parent = "F%s" % (parent) if parent and parent[0].isdigit() else None
    parent = os.path.basename(re.sub(r"\.[^\.\\]+$","", parent)).replace(".","_").replace("-","_").replace(" ","_").replace(",","_") if parent else None
    filepath = "%s\%s_%s" % (os.path.dirname(filepath), parent, os.path.basename(filepath)) if parent else filepath
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath
 
class Toolbox(object):
    def __init__(self):
        self.label =  "Display Mike Urban Results"
        self.alias  = "Display Mike Urban Results"

        # List of tool classes associated with this toolbox
        if arcgis_pro:
            self.tools = [DisplayMIKE1DResults, ReadMIKE1DResults, PlotRes1D]
        else:
            self.tools = [DisplayFloodReturnPeriodFun, DisplayWeirStatistics, DisplayFlowStatistics, DisplayQFullQMax, DisplayWeirReturnPeriod, DisplayMIKE1DResults, DisplayExtent]
class DisplayFloodReturnPeriodFun(object):
    def __init__(self):
        self.label       = "Display Flood Return Period"
        self.description = "Display Flood Return Period"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        critical_return_period = arcpy.Parameter(
            displayName="Critical Return Period (5 years, 10 years or 20 years)",
            name="critical_return_period",
            datatype="Double",
            parameterType="Optional",
            direction="Input")

        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        # exportShape = arcpy.Parameter(
        #     displayName="Output manholes with Flood Return Period",
        #     name="exportShape",
        #     datatype="DEFeatureClass",
        #     parameterType="Required",
        #     direction="Output")

        #
        # exportBasins = arcpy.Parameter(
        #     displayName="Output basins with Flood Return Period",
        #     name="exportBasins",
        #     datatype="DEFeatureClass",
        #     parameterType="Required",
        #     direction="Output")
            
        flowFile = arcpy.Parameter(
            displayName="Include PRF file in order to show maximum discharge from basins and permanent water volume",
            name="flowFile",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        flowFile.filter.list=["prf"]
        
        traceNetwork = arcpy.Parameter(
            displayName="Calculate connected catchment area for each basin",
            name="traceNetwork",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        traceNetwork.category = "Get Catchment Area"
        traceNetwork.enabled = False
        
        reaches = arcpy.Parameter(
            displayName="Trace network through:",
            name="reaches",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        reaches.filter.type = "ValueList"  
        reaches.filter.list = ["Orifice","Weir","Pump", "Basin"]
        reaches.value = ["Orifice","Weir","Pump", "Basin"]
        reaches.category = "Get Catchment Area"
        reaches.enabled = False
        
        breakChainOnNodes = arcpy.Parameter(
            displayName="End each trace at following node MUIDs (each node should by delimited by a comma: node1, node2)",
            name="breakChainOnNodes",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        breakChainOnNodes.category = "Additional settings"
        breakChainOnNodes.category = "Get Catchment Area"
        breakChainOnNodes.enabled = False

        use_critical_level = arcpy.Parameter(
            displayName="Use Critical Level from nodes",
            name="use_critical_level",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        use_critical_level.enabled = True

        parameters = [erfFile, observationPeriod, critical_return_period, mike_urban_database, flowFile, traceNetwork, reaches, breakChainOnNodes, use_critical_level]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        erfFile = parameters[0].ValueAsText
        observationPeriod = float(parameters[1].ValueAsText)
        critical_return_period = float(parameters[2].ValueAsText)
        mike_urban_database = parameters[3].ValueAsText
        msm_Node = mike_urban_database + "\msm_Node"
        msm_Link = mike_urban_database + "\msm_Link"
        msm_Weir = mike_urban_database + "\msm_Weir"
        msm_Orifice = mike_urban_database + "\msm_Orifice"
        msm_Pump = os.path.join(mike_urban_database,"msm_Pump")
        msm_CatchCon = os.path.join(mike_urban_database,"msm_CatchCon")
        ms_Catchment = os.path.join(mike_urban_database,"ms_Catchment")
        msm_HModA = os.path.join(mike_urban_database,"msm_HModA")
        msm_HParA = os.path.join(mike_urban_database,"msm_HParA")

        flowFile = parameters[4].ValueAsText
        traceNetwork = parameters[5].ValueAsText
        reaches = parameters[6].ValueAsText
        # break_chain_on_nodes = parameters[7].ValueAsText
        use_critical_level = parameters[8].Value

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)
        MIKE_gdb = os.path.join(MIKE_folder, os.path.splitext(os.path.basename(mike_urban_database))[0])
        no_dir = True
        dir_ext = 0
        while no_dir:
            try:
                if arcpy.Exists(MIKE_gdb):
                    os.rmdir(MIKE_gdb)
                os.mkdir(MIKE_gdb)
                no_dir = False
            except Exception as e:
                dir_ext += 1
                MIKE_gdb = os.path.join(MIKE_folder, "%s_%d" % (os.path.splitext(os.path.basename(mike_urban_database))[0], dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        arcpy.env.addOutputsToMap = False

        
        arcpy.SetProgressorLabel("Getting critical level for manholes")
        MUIDsHCrit = {}
        with arcpy.da.SearchCursor(msm_Node,["MUID","GroundLevel","CriticalLevel"]) as cursor:
            for row in cursor:
                if use_critical_level and row[2]:
                    MUIDsHCrit[row[0]] = row[2]
                else:
                    MUIDsHCrit[row[0]] = row[1]

        MUIDs = MUIDsHCrit.keys()
        if flowFile:
            arcpy.SetProgressorLabel("Getting Discharge from Basins")
            workingDir = os.path.dirname(__file__)
            M11OUT = workingDir + "\M11.OUT"
            M11IN = workingDir + "\M11.IN"

            prfFileCopy = workingDir + "\prffile.PRF"
            resultFile = workingDir + "\ResultFile.txt"
            
            copyfile(flowFile,prfFileCopy)
            call([m11extrapath(), prfFileCopy])

            lines = ""
            links = []
            with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
                for linei,line in enumerate(M11OUTFile):
                    try:
                        if ("Link_Q" in line or "Weir_Q" in line or "Orifice_Q" in line) and re.findall("<([^>]+)>",line)[0] not in links:		
                            links.append(re.findall("<([^>]+)>",line)[0])
                            line = re.sub("^0","1",line)
                        lines += line
                    except UnicodeDecodeError:
                        lines += line
                        arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                    except Exception as e:
                        raise(e)
                           
            with codecs.open(M11IN,'w','cp1252') as M11INFile:
                M11INFile.write(lines)

            call([m11extrapath(), prfFileCopy, resultFile])

            linksFlow = []
            with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
                txt = M11OUTFile.readlines()

            linksFlow = {}
            for link in links:
                linksFlow[link] = []
            lines = range([i for i,line in enumerate(txt) if "M11 DATA" in line][0]+1,len(txt))
            for i in lines:
                line = txt[i]
                values = line.split("  ")
                for valuei in range(1,len(values)):
                    linksFlow[links[valuei-1]].append(float(values[valuei]))
            os.remove(resultFile)
            os.remove(M11IN)
            
            lines = ""
            nodes = []
            basins = [row[0] for row in arcpy.da.SearchCursor(msm_Node, "MUID", where_clause = "TypeNo = 2")]
            with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
                for linei,line in enumerate(M11OUTFile):
                    try:
                        if "Node_WL" in line and re.findall("<([^>]+)>",line)[0] in basins:		
                            nodes.append(re.findall("<([^>]+)>",line)[0])
                            line = re.sub("^0","1",line)
                        lines += line
                    except UnicodeDecodeError:
                        lines += line
                        arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                    except Exception as e:
                        raise(e)
                           
            with codecs.open(M11IN,'w','cp1252') as M11INFile:
                M11INFile.write(lines)
            
            call([m11extrapath(), prfFileCopy, resultFile])
            nodesMinWL = {}
            with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
                txt = M11OUTFile.readlines()
            
            findValues = re.compile(r" ([\d\.]+)")
            for vali,val in enumerate(findValues.findall(txt[[i for i,line in enumerate(txt) if "M11 DATA" in line][0]+1])[1:]):
                nodesMinWL[nodes[vali]] = float(val)
            
            os.remove(prfFileCopy)
            os.remove(M11OUT)
            os.remove(resultFile)
            os.remove(M11IN)
            
            msm_LinkFromNode = {}
            msm_LinkToNode = {}
            with arcpy.da.SearchCursor(msm_Link,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]
            with arcpy.da.SearchCursor(msm_Weir,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]
            with arcpy.da.SearchCursor(msm_Orifice,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]

        arcpy.SetProgressorLabel("Reading ERF-file")
        dataTables = mousereader.readERF(erfFile,"MaxLevel_Ranked",MUIDs)
        # arcpy.AddMessage(dataTables)
        arcpy.SetProgressorLabel("Getting return period of flooding")
        MUIDsTCrit = {}
        MUIDs_elevation_at_crit = {}
        for i,MUID in enumerate(MUIDs):
            try:
                # arcpy.AddMessage(dataTables[i])
                nodeH = dataTables[i]["col2"]
                MUIDs_elevation_at_crit[MUID] = -1

                if MUIDsHCrit[MUID]>np.max(nodeH):
                    MUIDsTCrit[MUID] = 1000
                elif MUIDsHCrit[MUID]<np.min(nodeH):
                    MUIDsTCrit[MUID] = 0
                else:
                    MUIDsTCrit[MUID] = np.interp(MUIDsHCrit[MUID],
                          np.flipud(nodeH),
                          np.flipud(float(observationPeriod)/(np.array(range(len(nodeH)))+1)))
                    # arcpy.AddMessage((MUIDsHCrit[MUID], nodeH))
                if critical_return_period:
                    try:
                        MUIDs_elevation_at_crit[MUID] = np.interp(critical_return_period,
                                                           np.flipud(float(observationPeriod) /
                                                                     (np.array(range(len(nodeH))) + 1)),
                                                           np.flipud(nodeH))
                    except Exception as e:
                        pass
            except Exception as e:
                arcpy.AddError(MUID)
                arcpy.AddError(traceback.format_exc())
        
        arcpy.SetProgressorLabel("Creating manhole result file")
        msm_NodeNew = arcpy.CopyFeatures_management(msm_Node, getAvailableFilename(os.path.join(arcpy.env.scratchGDB, "NodeFloodReturn"),
                                                                          parent=mike_urban_database)).getOutput(0)
        # arcpy.AddMessage(msm_NodeNew)
        arcpy.AddField_management (msm_NodeNew, "TCrit", "DOUBLE", 10, 5)
        with arcpy.da.UpdateCursor(msm_NodeNew,["MUID","TCrit"]) as cursor:
            for row in cursor:
                if row[0] in MUIDsTCrit:
                    row[1] = MUIDsTCrit[row[0]]
                cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False

        def addLayer(layer_source, source, group=None, workspace_type = None, new_name=None,
                     definition_query=None):
            if not workspace_type:
                if ".mdb" in source:
                    workspace_type = "ACCESS_WORKSPACE"
                elif ".shp" in source:
                    workspace_type = "SHAPEFILE_WORKSPACE"
                elif ".gdb" in source:
                    workspace_type = "FILEGDB_WORKSPACE"
            if ".sqlite" in source:
                source_layer = arcpy.mapping.Layer(source)

                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, source_layer, "TOP")
                update_layer = arcpy.mapping.ListLayers(mxd, source_layer.name, df)[0]

                layer_source_mike_plus = layer_source.replace("MOUSE",
                                                              "MIKE+") if "MOUSE" in layer_source and os.path.exists(
                    layer_source.replace("MOUSE", "MIKE+")) else None
                layer_source = layer_source_mike_plus if layer_source_mike_plus else layer_source
                layer = arcpy.mapping.Layer(layer_source)
                update_layer.visible = layer.visible
                update_layer.labelClasses = layer.labelClasses
                update_layer.showLabels = layer.showLabels
                update_layer.name = layer.name
                update_layer.definitionQuery = definition_query

                try:
                    arcpy.mapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                except Exception as e:
                    arcpy.AddWarning(source)
                    pass
            else:
                layer = arcpy.mapping.Layer(layer_source)
                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, layer, "TOP")
                update_layer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name
                # arcpy.AddMessage((unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                #                                workspace_type, unicode(os.path.basename(source))))
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                               workspace_type, unicode(os.path.basename(source)))

        # arcpy.AddMessage((os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Node.lyr", msm_NodeNew, "FILEGDB_WORKSPACE"))
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Node.lyr", msm_NodeNew, workspace_type="FILEGDB_WORKSPACE")
        return
        arcpy.SetProgressorLabel("Creating basin result file")
        exportBasins = arcpy.Select_analysis(mike_urban_database + r"\mu_Geometry\msm_Node", getAvailableFilename(arcpy.env.scratchGDB + "\BasinFloodReturn",
                                                                         parent=mike_urban_database), where_clause = "TypeNo = 2").getOutput(0)
        fields = ["MUID", "GeometryID","CriticalLe" if ".shp" in exportBasins else "CriticalLevel", "Volume", "TCrit", "CatchArea", "CatchImpA", "Discharge", "VolCrit", 'LvlTCrit']

        arcpy.AddField_management(exportBasins, "Discharge", "FLOAT")
        arcpy.AddField_management(exportBasins, "Volume", "FLOAT")
        arcpy.AddField_management(exportBasins, "TCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "VolCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "CritLevel", "FLOAT")
        arcpy.AddField_management(exportBasins, "LvlTCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "CatchArea", "FLOAT")
        arcpy.AddField_management(exportBasins, "CatchImpA", "FLOAT")


        class Basin:
            def __init__(self, geometryID):
                self.geometryID = geometryID if geometryID else ""
                self.value1 = []
                self.value3 = []
                self.edges = []
                self.permanent_level = None
                self.invert_level = None

            class Edge:
                def __init__(self, name, uplevel):
                    self.name = name
                    self.uplevel = uplevel

            @property
            def critical_level(self):  # overwritten if critical level in msm_Node
                return np.max(self.elevations)

            @property
            def elevations(self):
                if np.min(self.value1) < self.invert_level:
                    return self.value1 + (self.invert_level - np.min(self.value1))
                else:
                    return self.value1

            @property
            def edges_sort(self):
                if len(self.edges) > 1:
                    idx_sort = np.argsort([edge.uplevel for edge in self.edges])
                    return [self.edges[i] for i in idx_sort]
                else:
                    return [self.edges[0]]

            @property
            def terrain_elevation(self):
                return [elevation for elevation in self.elevations if elevation < self.critical_level] + [
                    self.critical_level]

            @property
            def max_volume(self):
                idxSort = np.argsort(self.elevations)
                elevations = np.array(self.elevations)[idxSort]
                surface_areas = np.array(self.value3)[idxSort]
                elevations = [elevation for elevation in elevations if elevation < self.critical_level] + [
                    self.critical_level]
                surface_areas = np.interp(elevations, np.sort(self.elevations), surface_areas)
                if self.permanent_level:
                    return np.interp(np.max(elevations), elevations,
                                     [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations))) - np.interp(
                        self.permanent_level, elevations,
                        [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))
                else:
                    return np.interp(np.max(elevations), elevations, [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))

            @property
            def max_area(self):
                return np.max(self.value3)

            def get_volume(self, level):
                idxSort = np.argsort(self.elevations)
                elevations = np.array(self.elevations)[idxSort]
                surface_areas = np.array(self.value3)[idxSort]
                elevations = [elevation for elevation in elevations if elevation < self.critical_level] + [
                    self.critical_level]
                surface_areas = np.interp(elevations, np.sort(self.elevations), surface_areas)
                if self.permanent_level:
                    return np.interp(level, elevations,
                                     [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations))) - np.interp(
                        self.permanent_level, elevations,
                        [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))
                else:
                    return np.interp(level, elevations, [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))

        basins = {}
        with arcpy.da.SearchCursor(msm_Node, ["MUID", "GeometryID", "InvertLevel"], where_clause="TypeNo = 2") as cursor:
            for row in cursor:
                basins[row[0]] = Basin(row[1])
                basins[row[0]].invert_level = row[2]
        
        # if traceNetwork:
        #     arcpy.SetProgressorLabel("Creating basin result file - analyzing basin catchment area")
        #     import networkx as nx
        #     breakChainOnNodes = parameters[8].ValueAsText
        #     networkLinks = [msm_Link]
        #     networkLinks.append(msm_Orifice) if "Orifice" in reaches else None
        #     networkLinks.append(msm_Weir) if "Weir" in reaches else None
        #     networkLinks.append(msm_Pump) if "Pump" in reaches else None
        #
        #     network = nx.DiGraph()
        #
        #     hParADict = {}
        #     with arcpy.da.SearchCursor(msm_HParA,["MUID","RedFactor"]) as cursor:
        #         for row in cursor:
        #             hParADict[row[0]] = row[1]
        #
        #     catchmentImperviousnessDict = {}
        #     catchmentReductionFactor = {}
        #     with arcpy.da.SearchCursor(msm_HModA, ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor"]) as cursor:
        #         for row in cursor:
        #             catchmentImperviousnessDict[row[0]] = row[1]
        #             catchmentReductionFactor[row[0]] = hParADict[row[2]] if row[3] == 0 else row[4]
        #
        #     catchmentConnectionDict = {}
        #     with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
        #         for row in cursor:
        #             if row[1] not in catchmentConnectionDict:
        #                 catchmentConnectionDict[row[1]] = [row[0]]
        #             else:
        #                 catchmentConnectionDict[row[1]].append(row[0])
        #
        #     catchmentAreaDict = {}
        #     catchmentPersonsDict = {}
        #     catchmentImperviousAreaDict = {}
        #     catchmentReducedAreaDict = {}
        #     catchmentNetTypeNoDict = {}
        #     catchmentStatus = {}
        #     with arcpy.da.SearchCursor(ms_Catchment, ['MUID','SHAPE@AREA','Area','Persons',"NetTypeNo", 'Element_S']) as cursor:
        #         for row in cursor:
        #             catchmentPersonsDict[row[0]] = row[3] if row[3] is not None else 0
        #             catchmentAreaDict[row[0]] = row[2]*1e4 if row[2] is not None else row[1]
        #             catchmentNetTypeNoDict[row[0]] = row[4]
        #             catchmentStatus[row[0]] = row[5]
        #             if row[0] in catchmentImperviousnessDict:
        #                 catchmentImperviousAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]]
        #                 catchmentReducedAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]] * catchmentReductionFactor[row[0]]
        #             else:
        #                 arcpy.AddWarning("Warning: Could not find model record for Catchment %s" % (row[0]))
        #
        #     nodeTypeDict = {}
        #     nodeTypes = {1:u"Brønd",2:"Bassin",3:u"Udløb"}
        #     with arcpy.da.SearchCursor(msm_Node,["MUID","TypeNo"]) as cursor:
        #         for row in cursor:
        #             network.add_node(row[0])
        #             nodeTypeDict[row[0]] = nodeTypes[row[1]]
        #
        #     weights = {"msm_Link":1, "msm_Pump":1e4, "msm_Orifice":1e4, "msm_Weir":1e4}
        #     for networkLink in networkLinks:
        #         weight = weights[os.path.basename(networkLink)]
        #         with arcpy.da.SearchCursor(networkLink,["FromNode","ToNode","SHAPE@LENGTH"]) as cursor:
        #             for row in cursor:
        #                 network.add_edge(row[0],row[1], weight = weight*row[2])
        #
        #
        #     if not "Basin" in reaches:
        #         for basin in basins.values():
        #             for edge in network.out_edges(basin.MUID):
        #                 network.remove_edge(basin.MUID,edge[1])
        #                 arcpy.AddMessage("Removed edge %s-%s because tracing through basins is disabled" % (basin.MUID,edge[1]))
        #
        #     if breakChainOnNodes:
        #         breakEdges = [edge for edge in network.edges if edge[0] in re.findall("([^'^(),; \n]+)",breakChainOnNodes)]
        #         network.remove_edges_from(breakEdges)
        #         for edge in breakEdges:
        #             arcpy.AddMessage("Removed edge %s-%s because %s is included in list of nodes to end trace at" % (edge[0],edge[1]))
        #
        #     outlets = []
        #     junctions = []
        #     for node in list(network.nodes):
        #         if node:
        #             if not network.out_edges(node):
        #                 outlets.append(node)
        #             if len(network.out_edges(node))>1:
        #                 junctions.append(node)
        #
        #     for source in junctions:
        #         if source != None:
        #             lengths = np.ones((len(outlets),1))*1e9
        #             for i,target in enumerate(outlets):
        #                 try:
        #                     if nx.has_path(network, source, target):
        #                         lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
        #                 except:
        #                     arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
        #             toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[1]
        #             for edge in network.out_edges(source):
        #                 if not edge[1] == toNode:
        #                     network.remove_edge(source,edge[1])
        #                     arcpy.AddMessage("Removed edge %s-%s so that node %s exclusively leads to outlet %s" % (source,edge[1],source,outlets[np.argmin(lengths)]))
        #
        #     for basin in [row[0] for row in arcpy.da.SearchCursor(exportBasins, ["MUID"])]:
        #         nodesUpstream = nx.ancestors(network,basin)
        #         nodesUpstream.add(basin)
        #         catchIDs = [catchmentConnectionDict[n] for n in nodesUpstream if n in catchmentConnectionDict]
        #         catchIDs = [catchID for sublist in catchIDs for catchID in sublist]
        #
        #         basins[basin].total_catchment_area = round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100
        #         basins[basin].total_catchment_area_impervious = round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100

        with arcpy.da.SearchCursor(os.path.join(mike_urban_database, r"ms_TabD"), ["TabID", "Value1", "Value3"],
                                   where_clause="TabID IN ('%s')" % ("', '".join(
                                       [basin.geometryID for basin in basins.values()]))) as cursor:
            for row in cursor:
                basin = [basin for basin in basins.values() if basin.geometryID == row[0]][0]
                basin.value1.append(row[1])
                basin.value3.append(row[2])

        # arcpy.AddMessage(fields)
        arcpy.SetProgressorLabel("Creating basin result file - analyzing basin volume")
        with arcpy.da.UpdateCursor(exportBasins, fields) as cursor:
            for row in cursor:
                if flowFile:
                    basins[row[0]].permanent_level = nodesMinWL[row[0]]

                row[3] = basins[row[0]].max_volume
                try:
                    row[4] = MUIDsTCrit[row[0]]

                except Exception as e:
                    arcpy.AddWarning("Failed to find %s in General Report" % row[0])
                row[7] = 0
                if flowFile:
                    flow = np.array(())
                    for link in [MUID for MUID,FromNode in msm_LinkFromNode.iteritems() if FromNode==row[0]]:
                        if len(flow) == 0:
                            flow = np.array(linksFlow[link])
                        else:
                            flow = flow + np.array(linksFlow[link])
                    if len(flow)==0:
                        arcpy.AddMessage("Could not find link that discharges from basin %s (%s). Assuming intlet is also outlet." % (row[0],row[1]))
                        for link in [MUID for MUID,ToNode in msm_LinkToNode.iteritems() if ToNode==row[0]]:
                            if len(flow) == 0:
                                flow = np.array([a*(-1) for a in linksFlow[link]])
                            else:
                                flow = flow + np.array([a*(-1) for a in linksFlow[link]])
                    # arcpy.AddMessage(flow)
                    row[7] += np.max(flow)
                    
                # row[5] = basins[row[0]].total_catchment_area if basins[row[0]].total_catchment_area else 0
                # row[6] = basins[row[0]].total_catchment_area_impervious if basins[row[0]].total_catchment_area_impervious else 0
                if critical_return_period:
                    row[8] = basins[row[0]].get_volume(MUIDs_elevation_at_crit[row[0]])
                    row[9] = MUIDs_elevation_at_crit[row[0]]
                cursor.updateRow(row)
        
        # if flowFile:
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin Discharge.lyr", exportBasins, workspace_type="FILEGDB_WORKSPACE")
        # basinLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin Discharge.lyr")
        # # else:
        # #     basinLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin.lyr")
        # basinLayer = arcpy.mapping.AddLayer(df, basinLayer)
        # basinLayer = arcpy.mapping.ListLayers(mxd, basinLayer, df)[0]
        # basinLayer.replaceDataSource(os.path.dirname(exportBasins), "SHAPEFILE_WORKSPACE", os.path.basename(exportBasins).split(".")[0])
        #arcpy.RefreshTOC()
        #arcpy.RefreshActiveView()
        return

class DisplayWeirStatistics(object):
    def __init__(self):
        self.label       = "Display Weir Statistics"
        self.description = "Display Weir Statistics"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        htmlFile = arcpy.Parameter(
            displayName="Input LTS ERF File",
            name="htmlFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        htmlFile.filter.list=["html"]
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
            
        parameters = [htmlFile, mike_urban_database]
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        htmlFile = parameters[0].ValueAsText
        mike_urban_database = parameters[1].ValueAsText
        arcpy.env.addOutputsToMap = False
        msm_weir = arcpy.CopyFeatures_management(mike_urban_database + "\msm_Weir", getAvailableFilename(arcpy.env.scratchGDB + "\msm_Weir"))

        weirs = {}
        arcpy.SetProgressorLabel("Getting FROMNODE and TONODE for pipes")
        link_network = PipeNetwork(mike_urban_database, map_only = "weir")
        with arcpy.da.SearchCursor(mike_urban_database + "\msm_Weir",["MUID"]) as cursor:
            try:
                for row in cursor:
                    arcpy.AddMessage(row[0])
                    fromnode = link_network.weirs[row[0]].fromnode
                    tonode = link_network.weirs[row[0]].tonode if link_network.weirs[row[0]].tonode is not None else "0"
                    weirs["%s-%s" % (fromnode,tonode)] = row[0]
                    arcpy.AddMessage("%s-%s" % (fromnode,tonode))
            except Exception as e:
                arcpy.AddError(traceback.format_exc())
                raise(e)
                
        arcpy.SetProgressorLabel("Reading HTML-file")
        with codecs.open(htmlFile,'r',encoding='mbcs') as f:
            htmlFileTxt = f.read().split("\r\n")
            htmlFileTxt = np.array(htmlFileTxt)
            for line in htmlFileTxt:
                line = unicode(line)

        weirStart = [i for i,a in enumerate(htmlFileTxt) if "<H3>General statistics in weirs" in a][0]
        weirEnd = [i for i,a in enumerate(htmlFileTxt) if "TABLE" in a and i > weirStart][0]
        #nodeEnd = [i for i,a in enumerate(htmlFileTxt) if "TABLE" in a][1]
        htmlFileTxtWeirs = htmlFileTxt[weirStart:weirEnd]

        arcpy.SetProgressorLabel("Writing Shapefile")
        if "fromnode" not in [f.name for f in arcpy.ListFields(msm_weir[0])]:
            arcpy.AddField_management(msm_weir[0],"fromnode","TEXT")
        if "tonode" not in [f.name for f in arcpy.ListFields(msm_weir[0])]:
            arcpy.AddField_management(msm_weir[0],"tonode","TEXT")
        arcpy.AddField_management(msm_weir[0],"QVol","FLOAT")
        arcpy.AddField_management(msm_weir[0],"QNo","FLOAT")
        arcpy.AddField_management(msm_weir[0],"QHour","FLOAT")
        
        MUIDsQVol = {}
        MUIDsQNo = {}
        MUIDsQHours = {}
        # arcpy.AddMessage(weirs)
        getMUIDRe = re.compile(r"ALIGN=LEFT>([^<]+)")
        getQs = re.compile(r"<TD>([-0-9 <>\.]+)<\/TD><TD>([-0-9 <>\.]+)<\/TD><TD>([-0-9 <>\.]+)<\/TD><\/TR>$")
        for line in htmlFileTxtWeirs:
            if "ALIGN=LEFT" in line:
                try:
                    # arcpy.AddMessage(weirs["%s-%s" % (getMUIDRe.findall(line)[0],getMUIDRe.findall(line)[1])])
                    # arcpy.AddMessage("%s-%s" % (getMUIDRe.findall(line)[0],getMUIDRe.findall(line)[1]))
                    fromnode = getMUIDRe.findall(line)[0]
                    tonode = getMUIDRe.findall(line)[1] if not "Weir Outlet" in line else "0"
                    MUIDsQVol[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][0])
                    MUIDsQNo[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][1])
                    MUIDsQHours[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][2])
                except Exception as e:
                    arcpy.AddError("Error on line %s" % (line))
                    arcpy.AddError(traceback.format_exc())

        
        with arcpy.da.UpdateCursor(msm_weir[0],["MUID", "QVol", "QNo", "QHour", "fromnode", "tonode"]) as cursor:
            for row in cursor:
                if row[0] not in MUIDsQVol:
                    arcpy.AddError("Error: Could not find results of weir %s in result file" % (row[0]))
                else:
                    row[1] = MUIDsQVol[row[0]]
                    row[2] = MUIDsQNo[row[0]]
                    row[3] = MUIDsQHours[row[0]]
                    row[4] = link_network.weirs[row[0]].fromnode
                    row[5] = link_network.weirs[row[0]].tonode
                    cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        weirLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Weir.lyr")
        weirLayer = arcpy.mapping.AddLayer(df, weirLayer,'TOP')
        weirLayer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        weirLayer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE", os.path.basename(msm_weir[0]).split(".")[0])
        
        weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        #arcpy.RefreshTOC()
        
        return
        
class DisplayFlowStatistics(object):
    def __init__(self):
        self.label       = "Display Flow Statistics"
        self.description = "Display Flow Statistics"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")      
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database with links",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        exportShape = arcpy.Parameter(
            displayName="Output shapefile with link flow statistics",
            name="exportShape",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
            
        use_networker = arcpy.Parameter(
            displayName="Get fromnode and tonode attributes of pipes from spatial analysis?",
            name="use_networker",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="Optional",
            direction="Input")
        use_networker.value = False
        
        export_cad = arcpy.Parameter(
            displayName="Export CAD File with results and outlet flows",
            name="export_cad",
            datatype="File",
            parameterType="Optional",
            category="Additional Settings",
            direction="Output")
        # export_cad.filter.list=["dxf"]
            
        parameters = [erfFile, observationPeriod, mike_urban_database, exportShape, use_networker, export_cad]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        arcpy.env.addOutputsToMap = False
        erfFile = parameters[0].ValueAsText
        mike_urban_database = parameters[2].ValueAsText
        geometryFile = mike_urban_database + "\mu_Geometry\msm_Link"
        msm_Weir = mike_urban_database + "\mu_Geometry\msm_Weir"
        exportShape = parameters[3].ValueAsText
        observationPeriod = float(parameters[1].ValueAsText)
        use_networker = parameters[4].ValueAsText
        export_cad = parameters[5].ValueAsText
        
        MUIDs = {}
        if use_networker == "False":
            with arcpy.da.SearchCursor(geometryFile,["MUID","FROMNODE","TONODE"]) as cursor:
                for row in cursor:
                    MUIDs[row[0]] = "'%s', '%s'" % (row[1],row[2])
            with arcpy.da.SearchCursor(msm_Weir,["MUID","FROMNODE","TONODE"]) as cursor:
                for row in cursor:
                    MUIDs[row[0]] = "'%s', '%s'" % (row[1],row[2])
        else:
            link_network = PipeNetwork(mike_urban_database, map_only = "link weir").links
            for link in link_network.values():
                MUIDs[link.MUID] = "'%s', '%s'" % (link.fromnode, link.tonode)
            
        dataTables = mousereader.readERF(erfFile, "MaxFlow_Ranked", MUIDs.values(), ignore = True)

        msmLinkNew = arcpy.CopyFeatures_management(geometryFile,exportShape)
        
        if "FROMNODE" not in [field.name for field in arcpy.ListFields(msmLinkNew)] and "TONODE" not in [field.name for field in arcpy.ListFields(msmLinkNew)]:
            arcpy.AddField_management(msmLinkNew, "FROMNODE", "TEXT")
            arcpy.AddField_management(msmLinkNew, "TONODE", "TEXT")
            
        RPs = [1, 2, 5, 10, 20, 1000]
        fields = ["T1", "T2", "T5", "T10", "T20", "TMax"]
        for field in fields:
            arcpy.AddField_management (msmLinkNew, field, "FLOAT", field_precision = 8, field_scale = 5)
            
        if use_networker:
            fromnodes = [link.fromnode for link in link_network.values()]
            outlet_MUIDs = [link.MUID for link in link_network.values() if link.tonode not in fromnodes]
        else:
            fromnodes = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["FROMNODE"])]
            outlet_MUIDs = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["MUID","TONODE"]) if row[1] not in fromnodes]
        
        with arcpy.da.UpdateCursor(msmLinkNew,["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"]) as cursor:
            for row in cursor:
                try:
                    if row[0] in MUIDs.keys() and not dataTables[MUIDs.keys().index(row[0])] == None:
                            for i in range(len(RPs)):
                                try:
                                    flows = dataTables[MUIDs.keys().index(row[0])]["col2"]
                                    erfRPs = np.flipud(np.array(observationPeriod)/np.arange(1,len(flows)+1))
                                    row[i+1] = np.interp(RPs[i],erfRPs,np.flipud(flows))*1e3
                                    if use_networker:
                                        row[-3] = link_network[row[0]].fromnode
                                        row[-2] = link_network[row[0]].tonode
                                except Exception as e:
                                    arcpy.AddWarning(e)
                            
                # row[-3] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.X
                # row[-2] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.Y
                    cursor.updateRow(row)   
                except Exception as e:
                    arcpy.AddWarning("Failed on MUID %s: %s" % (row[0], e))
        
        with arcpy.da.InsertCursor(msmLinkNew, ["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"]) as cursor:
                extra_fields = ["FROMNODE","TONODE"] if not use_networker else []                
                row = [None] * len(["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"])
                with arcpy.da.SearchCursor(msm_Weir, ["MUID", "SHAPE@"] + extra_fields) as weir_cursor:
                    for weir_row in weir_cursor:
                        try:
                            if weir_row[0] in MUIDs.keys() and not dataTables[MUIDs.keys().index(weir_row[0])] == None:
                                row[0] = weir_row[0]
                                for i in range(len(RPs)):
                                    try:
                                        flows = dataTables[MUIDs.keys().index(weir_row[0])]["col2"]
                                        erfRPs = np.flipud(np.array(observationPeriod)/np.arange(1,len(flows)+1))
                                        row[i+1] = np.interp(RPs[i],erfRPs,np.flipud(flows))*1e3
                                        row[-1] = weir_row[1]
                                        if use_networker:
                                            row[-3] = link_network[weir_row[0]].fromnode
                                            row[-2] = link_network[weir_row[0]].tonode
                                    except Exception as e:
                                        arcpy.AddWarning(e)
                                    
                        # row[-3] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.X
                        # row[-2] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.Y
                                cursor.insertRow(row)   
                        except Exception as e:
                            arcpy.AddWarning("Failed on MUID %s: %s" % (weir_row[0], e))
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        add_layer = arcpy.mapping.Layer(exportShape)
        if (arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database))) 
            and arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database)))[0].isGroupLayer):
            group_layer = arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database)))[0]
            arcpy.mapping.AddLayerToGroup(df, group_layer, add_layer, "TOP")
            update_layer = [layer for layer in arcpy.mapping.ListLayers(mxd, add_layer.name, df) if layer.longName == os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring"]
            
            arcpy.mapping.AddLayerToGroup(df, group_layer, add_layer, "TOP")
            update_layer_outlet = [layer for layer in arcpy.mapping.ListLayers(mxd, add_layer.name, df) if layer.longName == os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring udløb"]
            
            arcpy.AddMessage(os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring")
        else:
            arcpy.mapping.AddLayer(df, add_layer, "TOP")
            update_layer = arcpy.mapping.ListLayers(mxd, add_layer.name, df)[0]
            
            arcpy.mapping.AddLayer(df, add_layer, "TOP")
            update_layer_outlet = arcpy.mapping.ListLayers(mxd, add_layer.name, df)[0]
        source_layer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_LinkQMax.lyr")
        arcpy.mapping.UpdateLayer(df,update_layer,source_layer,False)
        update_layer.replaceDataSource(unicode(add_layer.workspacePath), 'SHAPEFILE_WORKSPACE', unicode(add_layer.datasetName))
        update_layer.definitionQuery = "MUID NOT IN ('%s')" % ("', '".join(outlet_MUIDs))
        
        source_layer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_LinkQMaxOutlet.lyr")
        arcpy.mapping.UpdateLayer(df,update_layer_outlet,source_layer,False)
        update_layer_outlet.replaceDataSource(unicode(add_layer.workspacePath), 'SHAPEFILE_WORKSPACE', unicode(add_layer.datasetName))
        update_layer_outlet.definitionQuery = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))
        
        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        
        # if export_cad:
            # import ezdxf
            
            # tonodes = {}
            # if use_networker:
                # fromnodes = [link.fromnode for link in link_network.values()]
                # outlet_MUIDs = [link.MUID for link in link_network.values() if link.tonode not in fromnodes]
                # for link in link_network.values():
                    # tonodes[link.MUID] = link.tonode
            # else:
                # fromnodes = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["FROMNODE"])]
                # outlet_MUIDs = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["MUID","TONODE"]) if row[1] not in fromnodes]
                # with arcpy.da.SearchCursor(exportShape, ["FROMNODE", "MUID"], where_clause = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))):
                    # tonodes[row[1]] = row[0]
                
            
            # doc = ezdxf.new(dxfversion='R2010')
            # doc.layers.new(u'Udloeb', dxfattribs={'color': 255})
            # msp = doc.modelspace()
            # arcpy.AddMessage(outlet_MUIDs)
            # arcpy.AddMessage( "MUID IN ('%s')" % ("', '".join(outlet_MUIDs)))
            # with arcpy.da.SearchCursor(exportShape, ["SHAPE@", "MUID", "T1", "T2", "T5", "T10", "T20", "TMax"], where_clause = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))) as cursor:
                # for row in cursor:
                    # a = msp.add_mtext(
                    # 'ID: %s%sT1: %1.1f' % ("bay", "\n", ("%1.1f" % (max(row[2],0))).replace(".","."))) #tonodes[row[1]]
                    # a.set_location((row[0].lastPoint.X, row[0].lastPoint.Y))

        # # Save DXF document.
            # doc.saveas(export_cad)
        return
        
class DisplayQFullQMax(object):
    def __init__(self):
        self.label       = "Display Filling Degree of Pipe"
        self.description = "Display Filling Degree of Pipe"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        prfFile = arcpy.Parameter(
            displayName="PRF file",
            name="prfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        prfFile.filter.list=["prf"]
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database with links",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        exportShape = arcpy.Parameter(
            displayName="Output shapefile with Filling Degree of Pipe",
            name="exportShape",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        
        minimumSlope = arcpy.Parameter(
            displayName="Replace low slopes [o/oo]",
            name="minimumSlope",
            datatype="double",
            parameterType="Optional",
            direction="Input")
        minimumSlope.filter.list = [5, 10]
            
        parameters = [prfFile, mike_urban_database, exportShape, minimumSlope]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        arcpy.env.addOutputsToMap = False
        prfFile = parameters[0].ValueAsText
        mike_urban_database = parameters[1].ValueAsText
        geometryFile = mike_urban_database + "\mu_Geometry\msm_Link"
        exportShape = parameters[2].ValueAsText
        minimumSlope = float(parameters[3].ValueAsText)
        
        workingDir = os.path.dirname(__file__)
        M11OUT = workingDir + "\M11.OUT"
        M11IN = workingDir + "\M11.IN"

        prfFileCopy = workingDir + "\prffile.PRF"
        resultFile = workingDir + "\ResultFile.txt"
        replaceLowSlopes = 5

        if not minimumSlope:
            minimumSlope = 0
        copyfile(prfFile,prfFileCopy)
        call([m11extrapath(), prfFileCopy])

        lines = ""
        links = []
        with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
            for linei,line in enumerate(M11OUTFile):
                try:
                    if "Link_Q" in line:		
                        links.append(re.findall("<([^>]+)>",line)[0])
                        line = re.sub("^0","1",line)
                    lines += line
                except UnicodeDecodeError:
                    lines += line
                    arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                except Exception as e:
                    raise(e)
                       
        with codecs.open(M11IN,'w','cp1252') as M11INFile:
            M11INFile.write(lines)

        call([m11extrapath(), prfFileCopy, resultFile, "/MAX"])
        os.remove(prfFileCopy)
        os.remove(M11IN)
        os.remove(M11OUT)

        linksFlow = []
        with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
            for linei,line in enumerate(M11OUTFile):
                linksFlow.append(float(re.findall(" +([\-\d\.]+)",line)[0]))
        os.remove(resultFile)
        
        arcpy.CopyFeatures_management(os.path.join(mike_urban_database,"msm_Link"), exportShape)
        arcpy.AddField_management(exportShape, "QMax", "DOUBLE", field_scale = 5, field_precision = 10)
        arcpy.AddField_management(exportShape, "QFull", "DOUBLE", field_scale = 5, field_precision = 10)
        arcpy.AddField_management(exportShape, "Filldeg", "DOUBLE", field_scale = 5, field_precision = 10)
        with arcpy.da.UpdateCursor(exportShape, ["MUID","Diameter","Slope_C","MaterialID","QMax","QFull","FillDeg"]) as cursor:
            for row in cursor:
                try:
                    row[4] = linksFlow[[i for i,a in enumerate(links) if a == row[0]][0]]
                except Exception as e:
                    arcpy.AddMessage("Failed on row %s" % row)
                    raise(e)
                row[5] = colebrookWhite.QFull(row[1],max(row[2]*1e-2,minimumSlope*1e-3),row[3])
                row[6] = max(0,row[4]/row[5])
                cursor.updateRow(row)
            
        
                
        # mxd = arcpy.mapping.MapDocument("CURRENT")
        # df = arcpy.mapping.ListDataFrames(mxd)[0]
        # arcpy.env.addOutputsToMap = False
        # linkLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\QmaxQFull.lyr.lyr")
        # linkLayer = arcpy.mapping.AddLayer(df, linkLayer)
        # linkLayer = arcpy.mapping.ListLayers(mxd, linkLayer, df)[0]
        # arcpy.AddMessage(os.path.dirname(exportShape))
        # arcpy.AddMessage(os.path.dirname(os.path.basename(exportShape).split(".")[0]))
        # linkLayer.replaceDataSource(os.path.dirname(exportShape), "FILEGDB_WORKSPACE", os.path.basename(exportShape).split(".")[0])
        return
        
class DisplayWeirReturnPeriod(object):
    def __init__(self):
        self.label       = "Display Weir Return Period"
        self.description = "Display Weir Return Period"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        critical_return_period = arcpy.Parameter(
            displayName="Critical Return Period (5 years, 10 years or 20 years)",
            name="critical_return_period",
            datatype="Double",
            parameterType="Optional",
            direction="Input")

        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        parameters = [erfFile, observationPeriod, critical_return_period, mike_urban_database]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        erf_file = parameters[0].ValueAsText
        mike_urban_database = parameters[3].ValueAsText

        return_period = 38
        critical_return_period = 10
        
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        msm_Weir = os.path.join(mike_urban_database, "msm_Weir")

        class Weir:
            def __init__(self, muid, fromnode, tonode):
                self.fromnode = fromnode
                self.tonode = tonode if tonode else "0"
                self.muid  = muid

            result = {}

            @property
            def name(self):
                return "'%s', '%s'" % (self.fromnode, self.tonode)

            @property
            def events_count(self):
                return float(len(self.result['col2'])) if self.result['col2'] else 999

            @property
            def critical_discharge(self):
                # arcpy.AddMessage((self.name, critical_return_period, return_period / (1.0+np.arange(self.events_count)), self.result['col2'], np.min((1.0+np.arange(self.events_count)))))
                return np.interp(critical_return_period, np.flipud(return_period / (1.0+np.arange(self.events_count))), np.flipud(self.result['col2'])) if np.min(return_period / (1.0+np.arange(self.events_count))) <= critical_return_period else 0


        weirs = {}
        link_network = PipeNetwork(mike_urban_database, map_only = "weir")
        for weir in link_network.weirs:
            weirs[weir] = Weir(weir, link_network.weirs[weir].fromnode, link_network.weirs[weir].tonode)


        results = mousereader.readERF(erf_file, "Total_Discharge_Ranked", [weir.name for weir in weirs.values()])
        arcpy.AddMessage((erf_file, "Total_Discharge_Ranked", [weir.name for weir in weirs.values()]))

        for result, weir in zip(results, weirs.values()):
            arcpy.AddMessage((weir,result))
            if result is not None:
                weir.result = result
                print(weir.critical_discharge)
        
        empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
        empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(mike_urban_database))[0]
        
        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)
        MIKE_gdb = os.path.join(MIKE_folder, empty_group_layer.name)
        no_dir = True
        dir_ext = 0
        while no_dir:
            try:
                if arcpy.Exists(MIKE_gdb):
                    os.rmdir(MIKE_gdb)
                os.mkdir(MIKE_gdb)
                no_dir = False                
            except Exception as e:
                dir_ext += 1
                MIKE_gdb = os.path.join(MIKE_folder, "%s_%d" % (empty_group_layer.name, dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb
            
        msm_WeirNew = arcpy.CopyFeatures_management(msm_Weir, os.path.join(arcpy.env.scratchGDB, "Weir_return_period"))[0]

        arcpy.AddField_management (msm_WeirNew, "TCrit", "DOUBLE", 10, 5)
        arcpy.AddField_management (msm_WeirNew, "QCrit", "DOUBLE", 10, 5)
        
        with arcpy.da.UpdateCursor(msm_WeirNew, ["MUID", "TCrit", "QCrit"]) as cursor:
            for row in cursor:
                if weir.result:
                    # arcpy.AddMessage((row[0], weirs[row[0]].events_count, weirs[row[0]].result["col2"]))
                    row[1] = return_period / weirs[row[0]].events_count
                    row[2] = weirs[row[0]].critical_discharge
                    cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False

        def addLayer(layer_source, source, group = None, workspace_type = "ACCESS_WORKSPACE", new_name = None, definition_query = None):
            if ".sqlite" in source:
                source_layer = arcpy.mapping.Layer(source)
                
                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, source_layer, "BOTTOM")
                update_layer = arcpy.mapping.ListLayers(mxd, source_layer.name, df)[0]
                
                layer_source_mike_plus = layer_source.replace("MOUSE", "MIKE+") if "MOUSE" in layer_source and os.path.exists(layer_source.replace("MOUSE", "MIKE+")) else None
                layer_source = layer_source_mike_plus if layer_source_mike_plus else layer_source
                layer = arcpy.mapping.Layer(layer_source)
                update_layer.visible = layer.visible
                update_layer.labelClasses = layer.labelClasses
                update_layer.showLabels = layer.showLabels
                update_layer.name = layer.name
                update_layer.definitionQuery = definition_query
                
                try:
                    arcpy.mapping.UpdateLayer(df, update_layer, layer, symbology_only = True)
                except Exception as e:
                    arcpy.AddWarning(source)
                    pass
            else:
                layer = arcpy.mapping.Layer(layer_source)
                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, layer, "BOTTOM")
                update_layer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name
                arcpy.AddMessage((unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source))))
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source)))
                

        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Weir_return_period.lyr", msm_WeirNew, group = empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

class DisplayMIKE1DResults(object):
    def __init__(self):
        self.label = "a) Display MIKE1D Results"
        self.description = "a) Display MIKE1D Results"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        folder = arcpy.Parameter(
            displayName="Folder",
            name="folder",
            datatype="Folder",
            parameterType="Optional",
            direction="Input")
        if os.path.exists(r"C:\Papirkurv\Resultater"):
            folder.value = r"C:\Papirkurv\Resultater"

        node_featureclass = arcpy.Parameter(
            displayName="Nodes Result File",
            name="node_featureclass",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        node_featureclass.filter.list = ["POINT"]

        reach_featureclass = arcpy.Parameter(
            displayName="Reaches Result File",
            name="reach_featureclass",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        reach_featureclass.filter.list = ["LINE", "POLYLINE"]

        display_type = arcpy.Parameter(
            displayName="Display with fitting symbology",
            name="display_type",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        display_type.filter.list = ["Flood Volume", "Flood Depth", "Max Elevation / Headloss", "Surcharge Balance", "Peak Discharge", "Link Depth Difference", "Total Discharge"]
        display_type.value = ["Flood Volume", "Peak Discharge"]

        # res1d_filepath  = arcpy.Parameter(
        #     displayName="res1d Filepath",
        #     name="res1d_filepath",
        #     datatype="File",
        #     parameterType="Optional",
        #     direction="Input")
        # res1d_filepath.filter.list = ["res1d"]
        #
        # python3_path = arcpy.Parameter(
        #     displayName="Python 3 Filepath",
        #     name="python3_path",
        #     datatype="File",
        #     parameterType="Optional",
        #     category="IF RES1D: Python 3 Path with MIKEIO and Arcpy",
        #     direction="Input")
        # python3_path.filter.list = ["exe"]
        #
        # python3_path.value = r"C:\Users\elnn\AppData\Local\anaconda3\envs\myenv_py3_v2\python.exe"

        parameters = [folder, node_featureclass, reach_featureclass, display_type]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        # res1d_filepath = parameters[4]
        # python3_path = parameters[5]
        if parameters[0].altered and parameters[0].value and not parameters[1].value and not parameters[2].value:
            folder = parameters[0].ValueAsText
            import glob, os

            # Gather both shapefiles and file geodatabases
            candidate_files = glob.glob(os.path.join(folder, "*.shp"))
            candidate_files.extend(glob.glob(os.path.join(folder, "*.gdb")))

            # Get modification times and sort the candidates (newest first)
            file_mod_times = [(file, os.path.getmtime(file)) for file in candidate_files]

            for i, (file, _) in enumerate(file_mod_times):
                if ".gdb" in file and os.path.exists(os.path.join(file, "last_updated.txt")):
                    file_mod_times[i] = (file, os.path.getmtime(os.path.join(file, "last_updated.txt")))

            sorted_files = sorted(file_mod_times, key=lambda x: x[1], reverse=True)

            # Assign the newest file containing "nodes" to parameters[1]
            for file, _ in sorted_files:
                if "gdb" in file:
                    files_inside_gdb = [os.path.join(dirpath, f) for dirpath, _, files in
                     arcpy.da.Walk(file,
                                   datatype="FeatureClass") for f in files]
                    for subfile in files_inside_gdb:
                        if "nodes" in subfile.lower():
                            parameters[1].value = os.path.join(file, subfile)
                        elif "reaches" in subfile.lower():
                            parameters[2].value = os.path.join(file, subfile)
                    break
                if "nodes" in file:
                    parameters[1].Value = file
                    break

            # Assign the newest file containing "links" to parameters[2]
            if not parameters[2].Value:
                for file, _ in sorted_files:
                    if "links" in file:
                        parameters[2].Value = file
                        break

        # if res1d_filepath.value:
        #     parameters[0].enabled = False
        #     parameters[1].enabled = False
        #     parameters[2].enabled = False
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        nodes_featureclass = parameters[1].ValueAsText
        reaches_featureclass = parameters[2].ValueAsText
        display_type = parameters[3].ValueAsText
        # res1d_filepath = parameters[4].ValueAsText
        # python3_path = parameters[5].ValueAsText
        arcpy.AddMessage(display_type)

        if arcgis_pro:
            mxd = arcpy.mp.ArcGISProject("CURRENT")
            df = mxd.listMaps()[0]
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]

        empty_group_mapped = arcpymapping.LayerFile(os.path.dirname(
            os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr") if arcgis_pro else arcpy.mapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = df.addLayer(empty_group_mapped) if arcgis_pro else arcpymapping.AddLayer(df, empty_group_mapped,
                                                                                               "TOP")
        empty_group_layer = df.listLayers('Empty Group')[0] if arcgis_pro else \
            arcpymapping.ListLayers(mxd, "Empty Group", df)[0]
        if (nodes_featureclass and ".gdb" in nodes_featureclass) or (reaches_featureclass and ".gdb" in reaches_featureclass):
            empty_group_layer.name = os.path.basename(os.path.dirname(nodes_featureclass)).replace(".gdb","") if nodes_featureclass else os.path.dirname(reaches_featureclass)

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
                     definition_query=None):
            if arcgis_pro and not ".lyrx" in layer_source and os.path.exists(layer_source.replace(".lyr", ".lyrx")):
                layer_source = layer_source.replace(".lyr", ".lyrx")
            arcpy.AddMessage(layer_source)
            if source and ".sqlite" in source:
                source_layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpy.mapping.Layer(source)

                if group:
                    if arcgis_pro:
                        update_layer = df.addLayerToGroup(group, source_layer, "BOTTOM")
                    else:
                        arcpymapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    if arcgis_pro:
                        update_layer = df.addLayer(source_layer, "TOP")
                    else:
                        arcpymapping.AddLayer(df, source_layer, "TOP")

                if not arcgis_pro: update_layer = df.listLayers(mxd, source_layer.name, df)[0] if arcgis_pro else \
                    arcpy.mapping.ListLayers(mxd, source_layer.name, df)[0]

                if arcgis_pro:
                    new_connection_properties = update_layer.connectionProperties
                    new_connection_properties["workspace_factory"] = 'Sql'
                    new_connection_properties["connection_info"]["database"] = os.path.dirname(source)
                    update_layer.updateConnectionProperties()
                else:
                    if ".sqlite" in source:
                        layer = arcpymapping.Layer(layer_source)
                        update_layer.visible = layer.visible
                        update_layer.labelClasses = layer.labelClasses
                        update_layer.showLabels = layer.showLabels
                        update_layer.name = layer.name
                        update_layer.definitionQuery = definition_query

                        try:
                            arcpymapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                        except Exception as e:
                            arcpy.AddWarning(source)
                            pass
                    else:
                        update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                                       workspace_type, os.path.basename(source))

                try:
                    arcpymapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                except Exception as e:
                    arcpy.AddWarning(source)
                    pass
            else:
                layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpymapping.Layer(layer_source)
                if group:
                    if arcgis_pro:
                        df.addLayerToGroup(group, layer, "TOP")
                    else:
                        arcpymapping.AddLayerToGroup(df, group, layer, "TOP")
                else:
                    if arcgis_pro:
                        df.addLayer(layer, "TOP")
                    else:
                        arcpymapping.AddLayer(df, layer, "TOP")
                update_layer = df.listLayers(layer.listLayers()[0].name)[0] if arcgis_pro else \
                    arcpymapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name

                if source:
                    if arcgis_pro:
                        # CONFIRMED WORKING FOR SHAPEFILE -> FILEGDB
                        arcpy.AddMessage(update_layer)
                        cp = update_layer.connectionProperties
                        if workspace_type == "FILEGDB_WORKSPACE":
                            workspace_type = "File Geodatabase"
                        arcpy.AddMessage(workspace_type)
                        cp["connection_info"]['database'] = os.path.dirname(
                            source.replace(r"\mu_Geometry", ""))  # output db path+name
                        cp['dataset'] = os.path.basename(source)
                        cp['workspace_factory'] = workspace_type
                        update_layer.updateConnectionProperties(update_layer.connectionProperties, cp)
                    else:
                        update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                                       workspace_type, os.path.basename(source))
            return update_layer

        if nodes_featureclass:
            if "_spill.shp" in nodes_featureclass:
                arcpy.AddMessage(nodes_featureclass)

                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_spill.lyr",
                                 nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                 workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                 new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                layer.showLabels = True
            else:
                # layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                #          nodes_featureclass.replace(".shp",""), group=None, workspace_type = "SHAPEFILE_WORKSPACE", new_name = os.path.basename(nodes_featureclass).replace(".shp",""))

                if "flood volume" in display_type.lower():
                    layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_floodvol.lyr",
                                     nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                     workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                                     new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = True

                if "Flood Depth".lower() in display_type.lower():
                    arcpy.AddMessage(nodes_featureclass)
                    layer = addLayer(
                        os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_depthdiff.lyr",
                        nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                        workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                        new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = False

                if "Surcharge Balance".lower() in display_type.lower():
                    arcpy.AddMessage(nodes_featureclass)
                    layer = addLayer(
                        os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_Surcharge_balance.lyr",
                        nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                        workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE")
                    layer.showLabels = True

                if "headloss" in display_type.lower():
                    layer = addLayer(
                        os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                        nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                        workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                        new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = False

        if reaches_featureclass:
            if "depth difference" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_depthdiff.lyr",
                                 reaches_featureclass.replace(".shp",""), group=empty_group_layer, workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE", new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                layer.showLabels = False
            if "peak discharge" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links.lyr",
                                 reaches_featureclass.replace(".shp",""), group=empty_group_layer,
                                 workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                 new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                layer.showLabels = False
            if "total discharge" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_SumQ.lyr",
                                 reaches_featureclass.replace(".shp", ""), group=empty_group_layer,
                                 workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                 new_name=os.path.basename(reaches_featureclass).replace(".shp", ""))
                layer.showLabels = True
        if not arcgis_pro:
            arcpy.RefreshTOC()
        # def addLayer(layer_source, source):
        #     layer = arcpy.mapping.Layer(layer_source)
        #     layer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
        #     layer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        #     layer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
        #                                 os.path.basename(msm_weir[0]).split(".")[0])
        #
        #     weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        #
        # weirLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Weir.lyr")
        # weirLayer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
        # weirLayer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        # weirLayer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
        #                             os.path.basename(msm_weir[0]).split(".")[0])
        #
        # weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        return

class DisplayExtent(object):
    def __init__(self):
        self.label = "Display Extent of Dataframe"
        self.description = "Display Extent of Dataframe"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions
        parameters = []

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional

        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        import subprocess
        def copy2clip(txt):
            cmd = 'echo ' + txt.strip() + '|clip'
            return subprocess.check_call(cmd, shell=True)

        copy2clip("[%d, %d, %d, %d]" % (df.extent.lowerLeft.X, df.extent.lowerLeft.Y, df.extent.upperRight.X, df.extent.upperRight.Y))
        arcpy.AddMessage("[%d, %d, %d, %d]" % (df.extent.lowerLeft.X, df.extent.lowerLeft.Y, df.extent.upperRight.X, df.extent.upperRight.Y))
        return


class DrawSelection(object):
    def __init__(self):
        self.label = "2) Draw Results for Selection"
        self.description = "2) Draw Results for Selection"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")

        result_files = arcpy.Parameter(
            displayName="RES1D Network Result Files",
            name="result_files",
            datatype="File",
            multiValue=True,
            parameterType="optional",
            direction="Input")
        result_files.filter.list = ["res1d"]
        #
        # # # new_names (auto-filled)
        # # new_names = arcpy.Parameter(
        # #     displayName="New Names",
        # #     name="Rename of Result Files?",
        # #     datatype="String",
        # #     multiValue=True,
        # #     parameterType="Optional",
        # #     direction="Input"
        # # )
        #
        # pdf_output = arcpy.Parameter(
        #     displayName="PDF Output",
        #     name="pdf_output",
        #     datatype="File",
        #     parameterType="Optional",
        #     direction="Output")
        # # pdf_output.filter.list = ["pdf"]
        #
        # default_name = "Longitudinal Profiles.pdf"
        # default_path = os.path.join(arcpy.env.scratchFolder, default_name)
        # pdf_output.value = default_path
        #
        # overwrite_or_append = arcpy.Parameter(
        #     displayName="Append to PDF (Default is overwrite)",
        #     name="overwrite_or_append",
        #     datatype="Boolean",
        #     parameterType="optional",
        #     direction="Input")
        #
        # backup_tempfile = arcpy.Parameter(
        #     displayName="Backup Temp File Path",
        #     name="backup_tempfile",
        #     datatype="String",
        #     parameterType="Derived",  # hidden and not user editable
        #     direction="Output")


        parameters = [pipe_layer, result_files, ]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        backup_tempfile = parameters[4]
        if arcgis_pro:
            # Reference the active map in the current project
            aprx = arcpymapping.ArcGISProject("CURRENT")
            map_view = aprx.activeMap

            # List layers with selected features
            layers = None
            for layer in map_view.listLayers():
                try:
                    if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                        layers = layer.longName
                        break
                except:
                    pass
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
            layers = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                      lyr.getSelectionSet() if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'][0]

        if layers and not parameters[0].ValueAsText:
            parameters[0].value = layers

        pdf_output = parameters[2]
        overwrite_or_append = parameters[3]
        if pdf_output.ValueAsText and os.path.exists(pdf_output.ValueAsText):
            overwrite_or_append.enabled = True
            if overwrite_or_append.enabled:
                import tempfile
                import shutil
                with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as temp_file:
                    temp_backup_path = temp_file.name

                # Copy the existing file to this temp file
                shutil.copy2(pdf_output.ValueAsText, temp_backup_path)
                backup_tempfile.value = temp_backup_path


        else:
            overwrite_or_append.enabled = False

        # if parameters[1].altered and parameters[1].Values and not parameters[2].altered:  # result_files
        #     input_files = parameters[1].values
        #     cleaned_names = []
        #
        #     for path in input_files:
        #         basename = os.path.basename(str(path))
        #         basename = os.path.splitext(basename)[0]
        #         # Remove unwanted substrings
        #         cleaned = basename.replace("Default_Network_HD", "") \
        #             .replace("Default_Network", "") \
        #             .replace("Default", "") \
        #             .replace("HD", "")
        #         cleaned = cleaned.strip("_- ")  # Clean up any leftover junk
        #         cleaned_names.append(cleaned)
        #
        #     parameters[2].values = cleaned_names
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value

        if arcgis_pro:
            from mikeio1d import open as open_res1d
            from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach, QueryDataStructure
            from shapely import wkb
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            map_obj = aprx.activeMap
            view = aprx.activeView
            old_scale = map_obj.referenceScale
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]

        result_files = [f.replace("'", "") for f in parameters[1].ValueAsText.split(";")] if parameters[1].ValueAsText else None
        output_pdf = parameters[2].ValueAsText
        overwrite_or_append = parameters[3].Value
        backup_tempfile = parameters[4].Value

        pipe_layer_reference = pipe_layer
        for lyr in (map_obj.listLayers() if arcgis_pro else arcpy.mapping.ListLayers(mxd)):
            if lyr.name == pipe_layer.name and lyr.getSelectionSet():
                pipe_layer_reference = lyr
                break

        selection_set = pipe_layer_reference.getSelectionSet()


class ReadMIKE1DResults(object):
    def __init__(self):
        self.label = "1) Read MIKE1D Results"
        self.description = "1) Read MIKE1D Results"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        res1d_filepath = arcpy.Parameter(
            displayName="Res1D Filepath",
            name="res1d_filepath",
            datatype="File",
            multiValue=True,
            parameterType="Required",
            direction="Input")
        res1d_filepath.filter.list = ["res1d"]

        mike_database = arcpy.Parameter(
            displayName="MIKE+ database",
            name="mike_database",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        mike_database.filter.list = ["sqlite"]

        display_type = arcpy.Parameter(
            displayName="Display with fitting symbology",
            name="display_type",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        display_type.filter.list = ["Flood Volume", "Flood Depth", "Max Elevation / Headloss", "Surcharge Balance", "Peak Discharge", "Link Depth Difference", "Total Discharge"]
        display_type.value = ["Flood Volume", "Peak Discharge"]

        display_results = arcpy.Parameter(
            displayName="Display Results",
            name="display_results",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input"
        )
        display_results.value = True

        read_only_extent = arcpy.Parameter(
            displayName="Read only results in ArcGIS Pro Extent",
            name="read_only_extent",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input"
        )

        date_filter = arcpy.Parameter(
            displayName="Filter results to these dates (StartDate - EndDate)",
            name="date_filter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        date_filter.category = "Filter results"

        parameters = [res1d_filepath, mike_database, display_type, display_results, read_only_extent, date_filter]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        res1d_filepaths = [f.replace("'", "") for f in parameters[0].ValueAsText.split(";")] if parameters[0].ValueAsText else None
        if res1d_filepaths:
            res1d_filepath = res1d_filepaths[0]
        else:
            res1d_filepath = None

        display_results = parameters[3]
        display_type = parameters[2]
        if display_results.Value:
            display_type.enabled = True
        else:
            display_type.enabled = False

        if res1d_filepath and not parameters[1].Value:
            model_folder = os.path.dirname(os.path.dirname(res1d_filepath))
            MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".sqlite"

            if os.path.exists(MU_model):
                parameters[1].Value = MU_model
            else:
                model_folder = os.path.dirname(res1d_filepath)
                MU_model = os.path.join(model_folder, os.path.basename(model_folder)) + ".mdb"
                if os.path.exists(MU_model):
                    parameters[1].Value = MU_model
                    # print("Assuming MIKE database is %s" % (MU_model))
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach, QueryDataStructure
        res1d_filepaths = [f.replace("'", "") for f in parameters[0].ValueAsText.split(";")] if parameters[0].ValueAsText else None
        mike_database = parameters[1].ValueAsText
        display_type = parameters[2].ValueAsText
        display_results = parameters[3].Value
        read_only_extent = parameters[4].Value
        date_filter = parameters[5].Value

        if read_only_extent:
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            extent = aprx.activeView.camera.getExtent()

        for res1d_filepath in res1d_filepaths:
            nodes_featureclass, reaches_featureclass = readRes1D(res1d_filepath, mike_database, gdb_path = arcpy.env.scratchGDB, filter_to_extent = [extent.lowerLeft.X-50, extent.lowerLeft.Y-50, extent.upperRight.X+50, extent.upperRight.Y+50] if read_only_extent else None, date_filter = None)
            arcpy.AddMessage("BOBOBOBOB")
            arcpy.AddMessage(date_filter is None)
            nodes_featureclass = arcpy.Describe(nodes_featureclass).catalogPath
            reaches_featureclass = arcpy.Describe(reaches_featureclass).catalogPath

            if display_results:
                if arcgis_pro:
                    mxd = arcpy.mp.ArcGISProject("CURRENT")
                    df = mxd.listMaps()[0]
                else:
                    mxd = arcpy.mapping.MapDocument("CURRENT")
                    df = arcpy.mapping.ListDataFrames(mxd)[0]

                empty_group_mapped = arcpymapping.LayerFile(os.path.dirname(
                    os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr") if arcgis_pro else arcpy.mapping.Layer(
                    os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
                empty_group = df.addLayer(empty_group_mapped) if arcgis_pro else arcpymapping.AddLayer(df, empty_group_mapped,
                                                                                                       "TOP")
                empty_group_layer = df.listLayers('Empty Group')[0] if arcgis_pro else \
                    arcpymapping.ListLayers(mxd, "Empty Group", df)[0]
                # if (nodes_featureclass and ".gdb" in nodes_featureclass) or (reaches_featureclass and ".gdb" in reaches_featureclass):
                empty_group_layer.name = os.path.basename(res1d_filepath).replace(".res1d","").replace("Base","").replace("Result_file","")

                def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
                             definition_query=None):
                    if arcgis_pro and not ".lyrx" in layer_source and os.path.exists(layer_source.replace(".lyr", ".lyrx")):
                        layer_source = layer_source.replace(".lyr", ".lyrx")
                    arcpy.AddMessage(layer_source)
                    if source and ".sqlite" in source:
                        source_layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpy.mapping.Layer(source)

                        if group:
                            if arcgis_pro:
                                update_layer = df.addLayerToGroup(group, source_layer, "BOTTOM")
                            else:
                                arcpymapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                        else:
                            if arcgis_pro:
                                update_layer = df.addLayer(source_layer, "TOP")
                            else:
                                arcpymapping.AddLayer(df, source_layer, "TOP")

                        if not arcgis_pro: update_layer = df.listLayers(mxd, source_layer.name, df)[0] if arcgis_pro else \
                            arcpy.mapping.ListLayers(mxd, source_layer.name, df)[0]

                        if arcgis_pro:
                            new_connection_properties = update_layer.connectionProperties
                            new_connection_properties["workspace_factory"] = 'Sql'
                            new_connection_properties["connection_info"]["database"] = os.path.dirname(source)
                            update_layer.updateConnectionProperties()
                        else:
                            if ".sqlite" in source:
                                layer = arcpymapping.Layer(layer_source)
                                update_layer.visible = layer.visible
                                update_layer.labelClasses = layer.labelClasses
                                update_layer.showLabels = layer.showLabels
                                update_layer.name = layer.name
                                update_layer.definitionQuery = definition_query

                                try:
                                    arcpymapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                                except Exception as e:
                                    arcpy.AddWarning(source)
                                    pass
                            else:
                                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                                               workspace_type, os.path.basename(source))

                        try:
                            arcpymapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                        except Exception as e:
                            arcpy.AddWarning(source)
                            pass
                    else:
                        layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpymapping.Layer(layer_source)
                        if group:
                            if arcgis_pro:
                                df.addLayerToGroup(group, layer, "TOP")
                            else:
                                arcpymapping.AddLayerToGroup(df, group, layer, "TOP")
                        else:
                            if arcgis_pro:
                                df.addLayer(layer, "TOP")
                            else:
                                arcpymapping.AddLayer(df, layer, "TOP")
                        update_layer = df.listLayers(layer.listLayers()[0].name)[0] if arcgis_pro else \
                            arcpymapping.ListLayers(mxd, layer.name, df)[0]
                        if definition_query:
                            update_layer.definitionQuery = definition_query
                        if new_name:
                            update_layer.name = new_name

                        if source:
                            if arcgis_pro:
                                # CONFIRMED WORKING FOR SHAPEFILE -> FILEGDB
                                arcpy.AddMessage(update_layer)
                                cp = update_layer.connectionProperties
                                if workspace_type == "FILEGDB_WORKSPACE":
                                    workspace_type = "File Geodatabase"
                                arcpy.AddMessage(workspace_type)
                                cp["connection_info"]['database'] = os.path.dirname(
                                    source.replace(r"\mu_Geometry", ""))  # output db path+name
                                cp['dataset'] = os.path.basename(source)
                                cp['workspace_factory'] = workspace_type
                                update_layer.updateConnectionProperties(update_layer.connectionProperties, cp)
                            else:
                                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                                               workspace_type, os.path.basename(source))
                    return update_layer

                if nodes_featureclass:
                    if "_spill.shp" in nodes_featureclass:
                        arcpy.AddMessage(nodes_featureclass)

                        layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_spill.lyr",
                                         nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                         workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                         new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                        layer.showLabels = True
                    else:
                        # layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                        #          nodes_featureclass.replace(".shp",""), group=None, workspace_type = "SHAPEFILE_WORKSPACE", new_name = os.path.basename(nodes_featureclass).replace(".shp",""))

                        if "flood volume" in display_type.lower():
                            layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_floodvol.lyr",
                                             nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                             workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                                             new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                            layer.showLabels = True

                        if "Flood Depth".lower() in display_type.lower():
                            arcpy.AddMessage(nodes_featureclass)
                            layer = addLayer(
                                os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_depthdiff.lyr",
                                nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                                new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                            layer.showLabels = False

                        if "Surcharge Balance".lower() in display_type.lower():
                            arcpy.AddMessage(nodes_featureclass)
                            layer = addLayer(
                                os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_Surcharge_balance.lyr",
                                nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE")
                            layer.showLabels = True

                        if "headloss" in display_type.lower():
                            layer = addLayer(
                                os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                                nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                                new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                            layer.showLabels = False

                if reaches_featureclass:
                    if "depth difference" in display_type.lower():
                        layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_depthdiff.lyr",
                                         reaches_featureclass.replace(".shp",""), group=empty_group_layer, workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE", new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                        layer.showLabels = False
                    if "peak discharge" in display_type.lower():
                        layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links.lyr",
                                         reaches_featureclass.replace(".shp",""), group=empty_group_layer,
                                         workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                         new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                        layer.showLabels = False
                    if "total discharge" in display_type.lower():
                        layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_SumQ.lyr",
                                         reaches_featureclass.replace(".shp", ""), group=empty_group_layer,
                                         workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                         new_name=os.path.basename(reaches_featureclass).replace(".shp", ""))
                        layer.showLabels = True
                if not arcgis_pro:
                    arcpy.RefreshTOC()
            # def addLayer(layer_source, source):
            #     layer = arcpy.mapping.Layer(layer_source)
            #     layer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
            #     layer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
            #     layer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
            #                                 os.path.basename(msm_weir[0]).split(".")[0])
            #
            #     weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
            #
            # weirLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Weir.lyr")
            # weirLayer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
            # weirLayer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
            # weirLayer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
            #                             os.path.basename(msm_weir[0]).split(".")[0])
            #
            # weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        return


class PlotRes1D(object):
    def __init__(self):
        self.label = "2) Plot Res1D Results"
        self.description = "2) Plot Res1D Results"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        manhole_layer = arcpy.Parameter(
            displayName="Manhole feature layer",
            name="manhole_layer",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        manhole_layer.filter.list = ["Point"]

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")

        result_files = arcpy.Parameter(
            displayName="RES1D Network Result Files or DFS0 Rain Series",
            name="result_files",
            datatype="File",
            multiValue=True,
            parameterType="optional",
            direction="Input")
        result_files.filter.list = ["res1d", "dfs0"]

        stop_updating = arcpy.Parameter(
            displayName="Stop Updating Parameters",
            name="stop_updating",
            datatype="Boolean",
            parameterType="Derived",  # It's not an input!
            direction="Output"
        )
        stop_updating.enabled = False  # Hides it from the UI
        
        date_filter = arcpy.Parameter(
            displayName="Filter results to these dates (StartDate - EndDate)",
            name="date_filter",
            datatype="String",
            parameterType="Optional",
            direction="Input"
        )
        date_filter.category = "Filter results"

        step_every = arcpy.Parameter(
            displayName="Step every",
            name="step_every",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input"
        )
        step_every.value = 1
        step_every.category = "Filter results"

        font_size = arcpy.Parameter(
            displayName="Font Size",
            name="font_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input"
        )
        font_size.value = 7.0  # Default value
        font_size.category = "Additional Settings"

        parameters = [manhole_layer, pipe_layer, result_files, stop_updating, date_filter, step_every, font_size]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if not parameters[3].Value:
            parameters[3].Value = True
            if arcgis_pro:
                # Reference the active map in the current project
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap

                # List layers with selected features
                layers = None
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Point":
                            layers = layer.longName
                            break
                    except:
                        pass
            else:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd)[0]
                layers = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                          lyr.getSelectionSet() if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Point'][0]

            if layers and not parameters[0].ValueAsText and not parameters[0].altered:
                parameters[0].value = layers

            if arcgis_pro:
                # Reference the active map in the current project
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap

                # List layers with selected features
                layers = None
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            layers = layer.longName
                            break
                    except:
                        pass
            else:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                df = arcpy.mapping.ListDataFrames(mxd)[0]
                layers = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                          lyr.getSelectionSet() if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'][0]

            if layers and not parameters[1].ValueAsText:
                parameters[1].value = layers
            for parameter in [parameters[0], parameters[1]]:
                if parameter.ValueAsText:
                    if not parameters[2].value and ".gdb" in parameter.value.dataSource:
                        metadata_filepath = os.path.join(os.path.dirname(parameter.value.dataSource), "metadata")
                        # parameters[1].Value = [metadata_filepath]
                        if arcpy.Exists(metadata_filepath):
                            res1d_filepath = [row[0] for row in arcpy.da.SearchCursor(metadata_filepath, ["res1d_path"])][0]
                            if arcpy.Exists(res1d_filepath):
                                parameters[2].Value = [res1d_filepath]

        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        manhole_layer = parameters[0].ValueAsText
        pipe_layer = parameters[1].ValueAsText
        # dfs0_file = parameters[3].ValueAsText
        date_filter = parameters[4].ValueAsText
        step_every = parameters[5].Value
        font_size = parameters[6].Value

        if date_filter:
            # import dateparser

            def convertDate(date_str):
                date_str = date_str.strip()
                print(date_str)

                # Special case: only year
                if date_str.isdigit() and len(date_str) == 4:
                    return datetime.datetime(int(date_str), 1, 1)

                formats = [
                    "%d-%m-%Y",
                    "%d/%m/%Y",
                    "%d.%m.%Y",
                    "%d-%m-%y",
                    "%d/%m/%y",
                    "%d.%m.%y",
                    "%d-%m",  # Assume current year
                    "%d/%m",
                    "%d.%m",
                    "%Y-%m-%d",
                ]

                for fmt in formats:
                    try:
                        parsed = datetime.datetime.strptime(date_str, fmt)
                        # If no year was provided (default is 1900), use current year
                        if parsed.year == 1900:
                            parsed = parsed.replace(year=datetime.now().year)
                        return parsed
                    except ValueError:
                        continue

                raise ValueError(f"Failed to interpret date: {date_str}")

            time_filter = convertDate(date_filter.split(" - ")[0]), convertDate(date_filter.split(" - ")[1])

            arcpy.AddMessage(
                f"Filtering to Start: {time_filter[0].strftime('%Y-%m-%d %H:%M:%S')}, End: {time_filter[1].strftime('%Y-%m-%d %H:%M:%S')}")
        else:
            time_filter = None

        result_files = [f.replace("'", "") for f in parameters[2].ValueAsText.split(";")] if parameters[2].ValueAsText else None

        from mikeio1d.res1d import Res1D, QueryDataNode, QueryDataReach, QueryDataStructure

        manholes_selected = []
        pipes_selected = []
        if manhole_layer:
            manholes_selected = [row[0] for row in arcpy.da.SearchCursor(manhole_layer, ["MUID"])]

        if pipe_layer:
            pipes_selected = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]

        import matplotlib.pyplot as plt

        # Parameters
        subplots_count = 2 if manholes_selected and pipes_selected else 1
        fig_width_cm = 15.7  # Total figure width in centimeters
        fig_width_in = fig_width_cm / 2.54  # Convert to inches
        aspect_ratio = subplots_count  # Adjust for desired height (e.g., 0.6 for landscape-like)

        # Calculate figure height based on aspect ratio and number of rows (1 row here)
        fig_height_in = fig_width_in * aspect_ratio

        # Create subplots
        fig, axs = plt.subplots(subplots_count, 1, figsize=(fig_width_in, fig_height_in), dpi=300, sharex = True)

        # If only one subplot, make axs iterable
        if subplots_count == 1:
            axs = [axs]

        manhole_queries = {muid: QueryDataNode("WaterLevel", muid,0) for muid in manholes_selected}
        pipe_queries = {muid: QueryDataReach("Discharge", muid,0) for muid in pipes_selected}
        queries = {**manhole_queries, **pipe_queries}

        arcpy.SetProgressor("default", "Reading res1d")

        cmap = plt.cm.get_cmap('tab10' if arcgis_pro else "Set1")
        linestyles = ['-', '--', '-.', ':',
                      (0, (1, 1)),
                      (0, (5, 5)),
                      (0, (3, 5, 1, 5)),
                      (0, (3, 1, 1, 1))]

        plt.rcParams.update({'font.size': float(font_size)})
        plt.rcParams['font.family'] = 'Verdana'
        plt.rcParams['svg.fonttype'] = 'none'

        arcpy.AddMessage(font_size)
        import matplotlib.dates as dates

        for result_file_i, result_file in enumerate(result_files):
            if ".res1d" in result_file:
                res1d = Res1D(result_file, time = time_filter, step_every = step_every if step_every > 1 else None)
                try:
                    result_df = res1d.read([queries[key] for key in queries])
                except Exception as e:
                    import pandas as pd
                    arcpy.AddWarning(e)
                    dfs = []
                    for query in queries:
                        try:
                            dfs.append(res1d.read(queries[query]))
                        except:
                            pass
                    result_df = pd.concat(dfs, axis=1)

                columns = result_df.columns
                # arcpy.AddMessage(queries)
                col_i_WL = 0
                col_i_D = 0

                discharge_row = 1 if manholes_selected else 0
                axs[0].set_ylabel("Stuvningsniveau [m]")
                if pipes_selected:
                    axs[discharge_row].set_ylabel(u"Vandføring [L/s]")

                arcpy.SetProgressor("step", "Processing queries...", 0, len(columns), 1)

                linestyle = linestyles[result_file_i % len(linestyles)]

                def insert_gaps(series, gap_threshold: float):
                    """
                    From a pandas Series with a DatetimeIndex, return x and y arrays where
                    NaNs are inserted after large time gaps.

                    Parameters:
                        series (pd.Series): Time series with DatetimeIndex.
                        gap_threshold_seconds (float): Max allowed gap in seconds.

                    Returns:
                        x (np.ndarray): Index values with NaNs inserted (as datetime64).
                        y (np.ndarray): Series values with NaNs inserted.
                    """
                    idx = series.index
                    values = series.values

                    # Convert datetime to seconds since epoch for diffing
                    idx_seconds = idx.view(np.int64) / 1e9
                    diffs = np.diff(idx_seconds)
                    gap_locs = np.where(diffs > gap_threshold)[0]

                    x_out = []
                    y_out = []

                    for i in range(len(series)):
                        x_out.append(idx[i])
                        y_out.append(values[i])
                        if i in gap_locs:
                            x_out.append(np.datetime64('NaT'))
                            y_out.append(np.nan)

                    return np.array(x_out), np.array(y_out)

                i = 0
                for col in columns:
                    arcpy.SetProgressorLabel(f"Plotting result {i}/{len(columns)}")
                    arcpy.SetProgressorPosition(i)
                    i += 1

                    if len(result_files) > 1:
                        label = "%s (%s)" % (col.split(":")[1], os.path.basename(result_file).replace("Base", "").replace(".res1d",""))
                    else:
                        label = "%s" % (col.split(":")[1])

                    x, y = insert_gaps(result_df[col],
                                       gap_threshold=30 * 60)  # for datetime index, 30 min gap
                    arcpy.AddMessage(col)
                    if "waterlevel" in col.split(":")[0].lower():
                        axs[0].plot(x, y,  label = label, color = cmap(col_i_WL % 10), linestyle = linestyle, linewidth=0.8)
                        col_i_WL += 1

                    if "discharge" in col.split(":")[0].lower():
                        axs[discharge_row].plot(x, y*1e3, label = label, color=cmap(col_i_D % 10), linestyle = linestyle, linewidth=0.8)
                        col_i_D += 1
                        arcpy.AddMessage(col_i_D)

            elif ".dfs0" in result_file:
                for ax in axs:
                    import mikeio
                    xlim = ax.get_xlim()
                    dfs0 = mikeio.read(result_file).to_dataframe()
                    ax2 = ax.twinx()
                    ax2.step(dfs0.index, dfs0.values, 'k', linewidth = 0.5, label = u"Nedbør")
                    ax2.set_ylabel(r"Regnintensitet [µm/s]")
                    ax.set_xlim(xlim)

        arcpy.SetProgressor("default", "Showing Plot")


        for subplot_i in range(subplots_count):
            if len(queries)<9:
                axs[subplot_i].legend()
            locator = dates.AutoDateLocator(interval_multiples=True)
            axs[subplot_i].xaxis.set_major_locator(locator)



            def y_formatter(x, pos):
                # Format float with max 3 significant digits, no scientific notation
                # Use format specifier '.3f' but strip trailing zeros smartly
                s = f"{x:.3f}".rstrip('0').rstrip('.')
                # For very small numbers close to zero, just show 0
                if s == '-0':
                    s = '0'
                return s

            axs[subplot_i].yaxis.set_major_formatter(FuncFormatter(y_formatter))
            axs[subplot_i].yaxis.set_major_locator(MaxNLocator(nbins=6))  # max 6 ticks on y-axis

            for label in axs[subplot_i].get_xticklabels():
                label.set_rotation(45)
                label.set_horizontalalignment('right')

            axs[subplot_i].grid(which='major', linestyle='--', alpha=0.7)


        # axs[1].legend()
        plt.tight_layout()
        fig.autofmt_xdate()
        plt.show()
        # time.sleep(10)


