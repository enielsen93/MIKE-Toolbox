# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:21:31 2018

@author: eni
"""
import os
import arcpy
import numpy as np
import re
import hashlib
import math
import os
import sys
import networkx as nx
import ColebrookWhite
import networker
import pandas as pd
import traceback
import datetime
if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
    import pythonaddins
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument

import importlib.util
mikeio1d_installed = importlib.util.find_spec("mikeio1d") is not None

import mikegraph
import sqlite3
import timearea
from copy import deepcopy
import configparser


if "mapping" in dir(arcpy):
    import arcpy.mapping as apmapping
    from arcpy.mapping import MapDocument as MapDocument
else:
    import arcpy.mp as apmapping
    from arcpy.mp import ArcGISProject as MapDocument

diameters_plastic = [188, 235, 297, 377, 493, 588, 781, 985, 1185, 1385, 1485, 1585, 2000, 2200, 2400, 2600, 2800, 3000]
diameters_concrete = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2250, 2500, 3000, 3500]

def getAvailableFilename(filepath, parent = None):
    parent = "F%s" % (parent) if parent and parent[0].isdigit() else None
    parent = os.path.basename(re.sub(r"\.[^\.\\]+$","", parent)).replace(".","_").replace("-","_").replace(" ","_").replace(",","_") if parent else None
    filepath = "%s\%s_%s" % (os.path.dirname(filepath), parent, os.path.basename(filepath)) if parent else filepath
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
        # try:
            # arcpy.Delete_management(filepath)
            # return filepath
        # except:
            # i = 1
            # while arcpy.Exists(filepath + "%d" % i):
                # try:
                    # arcpy.Delete_management(filepath + "%d" % i)
                    # return filepath + "%d" % i
                # except:
                    # i += 1
            # return filepath + "%d" % i
    else:
        return filepath

class Config:
    def __init__(self, config_file):
        self.config_file = config_file

    def write(self, parameters):
        config = configparser.ConfigParser(allow_no_value=True)
        config.add_section("ArcGIS input parameters")
        for par_i, parameter in enumerate(parameters):
            if par_i == 0:
                config.set("ArcGIS input parameters", "# " + parameter.displayName)
            else:
                config.set(
                    "ArcGIS input parameters",
                    "\r\n# " + parameter.displayName)
            config.set("ArcGIS input parameters", str(parameter.name), str(parameter.value))

        with open(self.config_file,
                  "w") as file_write:
            config.write(file_write)

    def read(self):
        config = configparser.ConfigParser()
        config.read(self.config_file)
        parameters_dict = {}
        for option in config.options("ArcGIS input parameters"):
            parameters_dict[option] = config.get("ArcGIS input parameters", option)

        return parameters_dict

def confirm_assignment(txt, title = "Confirm", messagebox_typeno = 4, yes_return=None, no_return = None):
    import tkinter as tk
    from tkinter import messagebox

    root = tk.Tk()
    root.withdraw()  # Hide the main window
    if messagebox_typeno == 1:
        result = messagebox.askokcancel(title,txt)
        no_return = "Cancel" if no_return is None else no_return
        yes_return = "OK" if yes_return is None else yes_return
    elif messagebox_typeno == 4:
        result = messagebox.askyesno(title, txt)
        no_return = "No" if no_return is None else no_return
        yes_return = "Yes" if yes_return is None else yes_return
    else:
        result = messagebox.askyesno(title, txt)
        no_return = "No" if no_return is None else no_return
        yes_return = "Yes" if yes_return is None else yes_return
    root.destroy()
    if result:
        return yes_return
    else:
        return no_return

class Toolbox(object):
    def __init__(self):
        self.label =  "Pipe Dimension Tool"
        self.alias  = "Pipe Dimension Tool"
        self.canRunInBackground = True
        # List of tool classes associated with this toolbox
        self.tools = [PipeDimensionToolTAPro, upgradeDimensions, downgradeDimensions, setOutletLoss, reverseChange, InterpolateInvertLevels, GetMinimumSlope, CopyDiameter, CalculateSlopeOfPipe, ResetUpLevelDwlevel, SetDischargeRegulation, PipeDimensionToolResultFile, AnalyzeCatchmentArea, DrawLongitudinalProfiles, DrawClash]

class PipeDimensionToolTAPro(object):
    def __init__(self):
        self.label       = "1a) Automated Pipe Sizing | Time-Area Method"
        self.description = "1a) Automated Pipe Sizing | Time-Area Method"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        # pipe_layer.filter.list = ["Polyline"]

        reaches = arcpy.Parameter(
            displayName="Trace network through:",
            name="reaches",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        reaches.filter.type = "ValueList"
        reaches.filter.list = ["Orifice","Weir","Pump", "Basin"]
        
        result_field = arcpy.Parameter(
			displayName= "Field to assign diameter to",
			name="result_field",
			datatype="GPString",
			parameterType="Optional",
			direction="Input")
        # result_field.value = "Diameter"

        breakChainOnNodes = arcpy.Parameter(
            displayName="End each trace at following node MUIDs (each node should by delimited by a comma: node1, node2)",
            name="breakChainOnNodes",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        breakChainOnNodes.category = "Additional settings"

        runoff_file = arcpy.Parameter(
            displayName="Runoff rain event in ASCII Format",
            name="runoff",
            datatype="file",
            parameterType="Required",
            direction="Input")
        runoff_file.filter.list = ["txt","km2","kmd","csv"]

        scaling_factor = arcpy.Parameter(
			displayName= "Scaling Factor for runoff rain event",
			name="scaling_factor",
			datatype="GPString",
			parameterType="Required",
			direction="Input")
        scaling_factor.value = "1"

        useMaxInflow = arcpy.Parameter(
            displayName="Use Max. Inflow field as additional discharge to system (checkbox on msm_Node must be left unticked) - unit m3/s",
            name="useMaxInflow",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        useMaxInflow.category = "Additional settings"
        useMaxInflow.value = True

        slopeOverwrite = arcpy.Parameter(
            displayName="Overwrite slope of pipes to:",
            name="slopeOverwrite",
            datatype="double",
            parameterType="optional",
            direction="Input")
        slopeOverwrite.category = "Additional settings"

        writeDFS0 = arcpy.Parameter(
            displayName="Write csv file of runoff for every selected pipe to txt file:",
            name="writeDFS0",
            datatype="file",
            parameterType="optional",
            direction="Output")
        writeDFS0.category = "Additional settings"

        runoff_file = arcpy.Parameter(
            displayName="Runoff rain event in ASCII or DFS0" if mikeio1d_installed else "Runoff rain event in ASCII",
            name="runoff",
            datatype="file",
            parameterType="Required",
            direction="Input")
        if mikeio1d_installed:
            runoff_file.filter.list = ["txt","km2","kmd","csv", "dfs0"]
        else:
            runoff_file.filter.list = ["txt", "km2", "kmd", "csv"]

        keep_largest_diameter = arcpy.Parameter(
            displayName="Only change diameter if it's greater than existing",
            name="keep_largest_diameter",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        keep_largest_diameter.category = "Additional settings"
        keep_largest_diameter.value = False
        
        change_material = arcpy.Parameter(
            displayName="Change material (Plastic if less than 500 mm, concrete if greater than or equal to 500 mm",
            name="change_material",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        change_material.category = "Additional settings"
        change_material.value = True

        debug_output = arcpy.Parameter(
            displayName="Export Layer with attributes instead",
            name="debug_output",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")

        show_graphs = arcpy.Parameter(
            displayName="Show Discharge Graphs",
            name="show_graphs",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        show_graphs.category = "Additional settings"

        parameters = [pipe_layer, reaches, result_field, runoff_file, scaling_factor, breakChainOnNodes, useMaxInflow, slopeOverwrite, writeDFS0, keep_largest_diameter, change_material, debug_output, show_graphs]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        pipe_layer = parameters[0]
        runoff_file = parameters[3]
        scaling_factor = parameters[4]

        if arcgis_pro:
            # Reference the active map in the current project
            aprx = arcpymapping.ArcGISProject("CURRENT")
            map_view = aprx.activeMap

            if not pipe_layer.value:
                # List layers with selected features
                layers = None
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            if not pipe_layer.value:
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                     lyr.getSelectionSet() and "uplevel" in [field.name.lower() for field in arcpy.ListFields(lyr)]
                     and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and (
                                 "sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                     and lyr.visible]

                if links:
                    parameters[0].value = links[0]

        if pipe_layer.value:
            parameters[2].filter.list = [f.name for f in arcpy.ListFields(pipe_layer.ValueAsText)]
            parameters[2].Value = "Diameter" if "diameter" in [field.lower() for field in parameters[2].filter.list] else parameters[2].Value


        if pipe_layer.ValueAsText and not runoff_file.value:
            MU_database = os.path.dirname(arcpy.Describe(pipe_layer.ValueAsText).catalogPath).replace("\mu_Geometry", "")
            MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
            config_folder = os.path.join(MIKE_folder, "Config")
            config_file = os.path.join(config_folder, os.path.splitext(os.path.basename(MU_database))[0] + ".ini")
            if os.path.exists(config_file):
                config = Config(config_file)
                parameters_dict = config.read()
                for config_parameter in parameters_dict:
                    # raise(Exception([parameter.name for par_i, parameter in enumerate(parameters)]))
                    try:
                        i = [par_i for par_i, parameter in enumerate(parameters) if parameter.name == config_parameter][0]
                        parameters[i] = parameters_dict[config_parameter]
                    except Exception as e:
                        pass

            if arcgis_pro and mikeio1d_installed:
                table_path = MU_database + r"\msm_BBoundary"
                fields = ["muid", "applyboundaryno", "fraction", "variationno", "tsconnection"]

                with arcpy.da.SearchCursor(table_path, fields) as cursor:
                    for row in cursor:
                        id, applyboundaryno, fraction, variationno, tsconnection = row
                        # Do your thing with the variables
                        if applyboundaryno and variationno == 3:
                            scaling_factor.value = fraction

                            runoff_file.Value = os.path.abspath(os.path.join(os.path.dirname(MU_database), tsconnection))
                            break




        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].ValueAsText
        MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "")
        MU_database = MU_database.replace("!delete!", "")
        arcpy.AddMessage(MU_database)
        reaches = parameters[1].ValueAsText
        reaches = reaches + ";Link" if reaches else "Link"
        result_field = parameters[2].ValueAsText
        runoff_file = parameters[3].ValueAsText
        scaling_factor = float(parameters[4].Value)
        breakChainOnNodes = parameters[5].ValueAsText
        useMaxInflow = parameters[6].Value
        slopeOverwrite = parameters[7].Value
        writeDFS0 = parameters[8].ValueAsText
        keep_largest_diameter = parameters[9].Value
        change_material = parameters[10].Value
        debug_output = parameters[11].Value
        show_graph = parameters[12].Value

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)

        MIKE_gdb = os.path.join(MIKE_folder, os.path.splitext(os.path.basename(MU_database))[0])
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
                MIKE_gdb = os.path.join(MIKE_folder,
                                        "%s_%d" % (os.path.splitext(os.path.basename(MU_database))[0], dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        arcpy.SetProgressorLabel("Preparing")
        selected_pipes = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["muid"])]

        if arcgis_pro:
            mxd = arcpy.mp.ArcGISProject("CURRENT")
            df = mxd.listMaps()[0]
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]


        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE"):
            layer = apmapping.Layer(layer_source)
            if group:
                apmapping.AddLayerToGroup(df, group, layer, "BOTTOM")
            else:
                apmapping.AddLayer(df, layer, "TOP")
            updatelayer = apmapping.ListLayers(mxd, layer.name, df)[0]
            updatelayer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))), workspace_type,
                                          unicode(os.path.basename(source)))

        is_sqlite = True if ".sqlite" in MU_database else False

        msm_Link = os.path.join(MU_database, "msm_Link")

        arcpy.SetProgressorLabel("Mapping Network")
        graph = mikegraph.Graph(MU_database, useMaxInFlow = useMaxInflow)
        graph.map_network()

        if breakChainOnNodes:
            breakEdges = [edge for edge in graph.graph.edges if
                          edge[0] in re.findall("([^'^(),; \n]+)", breakChainOnNodes)]
            graph.network.remove_edges_from(breakEdges)
            for edge in breakEdges:
                arcpy.AddMessage(
                    "Removed edge %s-%s because %s is included in list of nodes to end trace at" % (edge[0], edge[1]))


        arcpy.AddMessage(graph.maxInflow)

        arcpy.SetProgressorLabel("Reading Rain Series")
        rainseries = timearea.TimeArea(runoff_file)

        rainseries.additional_discharge = graph.maxInflow
        rainseries.scaling_factor = scaling_factor

        target_manholes = [graph.network.links[link].fromnode for link in selected_pipes]
        arcpy.SetProgressor("step", "Tracing to every pipe selected", 0, len(target_manholes), 1)

        timearea_curves = {}
        peak_discharge = {}
        peak_discharge_time = {}

        graphs_count = 0

        for target_i, target_manhole in enumerate(target_manholes):
            arcpy.SetProgressorPosition(target_i)
            timearea_curves[target_manhole] = rainseries.timeareaCurve(target_manhole, graph)

            if show_graph and graphs_count < 15:
                if graphs_count > 15:
                    arcpy.AddMessage("Displaying no more than 15 graphs!")
                else:
                    graphs_count += 1
                    old_setting = arcpy.env.addOutputsToMap
                    arcpy.env.addOutputsToMap = False
                    table = arcpy.management.CreateTable(arcpy.env.scratchGDB, "Tab" + target_manhole, os.path.dirname(
                        os.path.realpath(__file__)) + "\Data\PipeDimensionTool\Template.dbf")[0]

                    rationel_curves = rainseries.rationelCurve(target_manhole, graph)

                    with arcpy.da.InsertCursor(table, ["Disch_ta", "Disch_rat"]) as cursor:
                        for discharge_ta, discharge_rat in zip(timearea_curves[target_manhole], rationel_curves):
                            cursor.insertRow([discharge_ta, discharge_rat])

                    table_view = arcpy.MakeTableView_management(table, r"%s_tv" % os.path.basename(table))

                    arcpy.env.addOutputsToMap = True
                    gr = arcpy.Graph()

                    oid_fieldname = arcpy.Describe(table).OIDFieldName

                    disch_ta_plot = gr.addSeriesAreaVertical(table_view, 'Disch_ta', oid_fieldname)
                    # arcpy.AddMessage(dir(gr.graphPropsGeneral))
                    gr.addSeriesLineVertical(table_view, 'Disch_rat', oid_fieldname)
                    gr.graphPropsGeneral.title = target_manhole

                    graph_template = os.path.dirname(
                        os.path.realpath(__file__)) + "\Data\PipeDimensionTool\graph_template_rev1.grf"
                    arcpy.MakeGraph_management(graph_template, gr, target_manhole)
                    arcpy.env.addOutputsToMap = old_setting
                # return

            peak_discharge[target_manhole] = np.max(timearea_curves[target_manhole])
            peak_discharge_time[target_manhole] = np.argmax(timearea_curves[target_manhole])


        if writeDFS0:
            dfs0_text = np.empty((2 + len(timearea_curves[timearea_curves.keys()[0]])), dtype=object)
            dfs0_text[0] = "\t".join(["Discharge[meter^3/sec]:Instantaneous"] * len(timearea_curves.keys()))
            dfs0_text[1] = "Time"
            for i in range(len(dfs0_text) - 2):
                dfs0_text[i + 2] = str(rainseries.series.index[0] + datetime.timedelta(minutes=i))

            for target_manhole in sorted(timearea_curves.keys()):
                dfs0_text[1] += "\t" + target_manhole
                for discharge_i, discharge in enumerate(timearea_curves[target_manhole] / 1e3):
                    dfs0_text[discharge_i + 2] += "\t%1.6f" % (discharge)
            filepath = writeDFS0
            with open(filepath, 'w') as f:
                for i in range(len(dfs0_text)):
                    f.write(dfs0_text[i] + "\n")

        arcpy.SetProgressorLabel("Calculating Pipe Dimensions")
        arcpy.AddMessage(peak_discharge)
        arcpy.AddMessage(peak_discharge_time)

        def addField(shapefile, field_name, datatype):
            i = 1
            while field_name in [f.name for f in arcpy.Describe(shapefile).fields]:
                field_name = "%s_%d" % (field_name, i)
            arcpy.AddField_management(shapefile, field_name, datatype)
            return field_name

        if debug_output:
            try:
                result_layer = getAvailableFilename(arcpy.env.scratchGDB + "\Pipe_Dimensions", parent=MU_database)
                arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, os.path.basename(result_layer), "POLYLINE")
                fields = ["muid", "diameter", "materialid", "slope", "nettypeno", "enabled", "ImpArea", "RedArea",
                          "MaxFlow", "CritTime"]

                addField(result_layer, "muid", "TEXT")
                addField(result_layer, "diameter", "FLOAT")
                addField(result_layer, "materialid", "TEXT")
                addField(result_layer, "slope", "FLOAT")
                addField(result_layer, "nettypeno", "SHORT")
                addField(result_layer, "enabled", "SHORT")
                addField(result_layer, "ImpArea", "FLOAT")
                addField(result_layer, "RedArea", "FLOAT")
                addField(result_layer, "MaxFlow", "FLOAT")
                addField(result_layer, "CritTime", "SHORT")
                with arcpy.da.InsertCursor(result_layer, ["SHAPE@"] + fields) as ins_cursor:
                    with arcpy.da.SearchCursor(msm_Link,
                                               ["Slope" if is_sqlite else "Slope_C", "Diameter", "MaterialID", "MUID",
                                                "SHAPE@", "NetTypeNo",
                                                "enabled"],
                                               where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                        for row_i, row in enumerate(cursor):
                            diameter_old = row[4]
                            arcpy.SetProgressorPosition(row_i)
                            if change_material:
                                D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                    diameter for diameter in diameters_concrete if diameter > 450]
                            else:
                                D = diameters_plastic if not "concrete" in row[1].lower() and not "beton" in row[
                                    1].lower() else diameters_concrete
                            slope = slopeOverwrite if slopeOverwrite else row[0] * 1e-2
                            # if writeDischargeInstead:
                            #     diameter = # ins_row = (row[4], row[3], peak_discharge[msm_Link_Network.links[row[3]].fromnode], row[2], row[0], row[5], row[6])
                            #     ins_cursor.insertRow(ins_row)
                            QFull = 0
                            Di = -1

                            if peak_discharge[graph.network.links[row[3]].fromnode] == 0:
                                Di = 0
                            else:
                                while QFull is not None and QFull * 1e3 < peak_discharge[
                                    graph.network.links[row[3]].fromnode] and Di + 1 < len(D):
                                    Di += 1
                                    QFull = ColebrookWhite.QFull(D[Di] / 1e3, slope, row[2])
                                    # arcpy.AddMessage((D[Di], QFull))
                                    # arcpy.AddMessage(QFull)
                            diameter = D[Di] / 1.0e3

                            material = row[2]
                            if keep_largest_diameter and diameter < row[1]:
                                diameter = row[1]
                            else:
                                if change_material:
                                    material = "Concrete (Normal)" if diameter > 0.45 else "Plastic"
                                arcpy.AddMessage("Changed %s from %d to %d" % (row[3], row[1] * 1e3, D[Di]))

                            upstream_nodes = graph.find_upstream_nodes(graph.network.links[row[3]].fromnode)
                            catchments = [graph.find_connected_catchments(node) for node in upstream_nodes][0]

                            ins_row = (row[4], row[3], diameter, material, row[0], row[5], row[6],
                                       np.sum([catchment.impervious_area for catchment in catchments]),
                                       np.sum([catchment.reduced_area for catchment in catchments]),
                                       peak_discharge[graph.network.links[row[3]].fromnode],
                                       peak_discharge_time[graph.network.links[row[3]].fromnode])
                            ins_cursor.insertRow(ins_row)
                            # if diameter_old != row[1]:
                            # arcpy.AddMessage("Changed diameter from %1.2f to %1.2f for pipe %s" % (diameter_old, row[1], row[3]))
            except Exception as e:
                if row[0] not in peak_discharge:
                    arcpy.AddError(traceback.format_exc())
                    arcpy.AddError("Failed to analyze catchments connected to node %s on pipe %s" % (
                    graph.network.links[row[3]].fromnode, graph.network.links[row[3]].MUID))
                    raise (e)
                arcpy.AddError(row)
                arcpy.AddError(traceback.format_exc())
                raise (e)
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Links Dimensioned.lyr", result_layer,
                     workspace_type="FILEGDB_WORKSPACE")
        elif is_sqlite:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                with arcpy.da.SearchCursor(msm_Link,
                                           ["Slope" if is_sqlite else "Slope_C", "Diameter", "MaterialID", "MUID",
                                            "SHAPE@",
                                            "NetTypeNo",
                                            "enabled"],
                                           where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                    for row_i, row in enumerate(cursor):
                        arcpy.SetProgressorPosition(row_i)
                        old_material = row[2]
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[1].lower() and not "beton" in row[
                                1].lower() else diameters_concrete
                        slope = slopeOverwrite if slopeOverwrite else row[0] * 1e-2
                        # if writeDischargeInstead:
                        #     diameter = # ins_row = (row[4], row[3], peak_discharge[msm_Link_Network.links[row[3]].fromnode], row[2], row[0], row[5], row[6])
                        #     ins_cursor.insertRow(ins_row)
                        QFull = 0
                        Di = -1

                        if peak_discharge[graph.network.links[row[3]].fromnode] == 0:
                            Di = 0
                        else:
                            while QFull is not None and QFull * 1e3 < peak_discharge[
                                graph.network.links[row[3]].fromnode] and Di + 1 < len(D):
                                Di += 1
                                QFull = ColebrookWhite.QFull(D[Di] / 1e3, slope, row[2])
                        diameter = D[Di] / 1.0e3

                        material = row[2]
                        if keep_largest_diameter and diameter <= row[1]:
                            diameter = row[1]
                        else:
                            if change_material:
                                material = "Concrete (Normal)" if diameter > 0.45 else "Plastic"
                            def materialName(material_id):
                                if "pl" in material_id.lower():
                                    material_name = "pl"
                                elif "concrete" in material_id.lower():
                                    material_name = "bt"
                                else:
                                    material_name = ""
                                return material_name

                            arcpy.AddMessage("Changed link '%s': %d%s → %d%s" % (row[3], row[1] * 1e3, materialName(old_material), D[Di], materialName(material)))
                            # arcpy.AddMessage("UPDATE msm_Link SET Diameter = %1.3f, SET MaterialID = %s WHERE MUID = %s" % (diameter, material, row[3]))
                            update_cursor.execute(
                                "UPDATE msm_Link SET Diameter = %1.3f, MaterialID = '%s' WHERE MUID = '%s'" % (
                                diameter, material, row[3]))

        else:
            try:
                edit = arcpy.da.Editor(MU_database)
                edit.startEditing(False, True)
                edit.startOperation()

                fields = ["Slope_C", "MaterialID", "MUID", "Diameter"]
                if result_field and result_field not in fields:
                    result_field_index = len(fields) - 1
                    fields.append(result_field)

                arcpy.AddMessage(fields)
                with arcpy.da.UpdateCursor(msm_Link, fields,
                                           where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                    for row_i, row in enumerate(cursor):
                        # diameter_old = row[4]
                        arcpy.SetProgressorPosition(row_i)
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[1].lower() and not "beton" in row[
                                1].lower() else diameters_concrete
                        slope = slopeOverwrite if slopeOverwrite else row[0] * 1e-2
                        QFull = 0
                        Di = -1

                        if peak_discharge[graph.network.links[row[2]].fromnode] == 0:
                            Di = 0
                        else:
                            while QFull is not None and QFull * 1e3 < peak_discharge[
                                graph.network.links[row[2]].fromnode] and Di + 1 < len(D):
                                Di += 1
                                QFull = ColebrookWhite.QFull(D[Di] / 1e3, slope, row[1])

                        diameter = D[Di]
                        if keep_largest_diameter and diameter / 1.0e3 <= row[3]:
                            diameter = row[3]
                        else:
                            arcpy.AddMessage(
                                "Changed %s from %d to %d" % (row[2], row[3] * 1e3 if row[3] else 0, D[Di]))
                            row[-1] = D[Di] / 1.0e3
                            if change_material:
                                row[1] = "Concrete (Normal)" if row[-1] > 0.45 else "Plastic"
                            # if diameter_old != row[1]:
                            # arcpy.AddMessage("Changed diameter from %1.2f to %1.2f for pipe %s" % (diameter_old, row[1], row[3]))
                        try:
                            cursor.updateRow(row)
                        except Exception as e:
                            arcpy.AddWarning("Could not update row:")
                            arcpy.AddWarning(row)
                edit.stopOperation()
                edit.stopEditing(True)
            except Exception as e:
                arcpy.AddError(row)
                arcpy.AddError(traceback.format_exc())
                raise (e)
        return


class PipeDimensionToolResultFile(object):
    def __init__(self):
        self.label = "1b) Automated Pipe Sizing | MIKE1D Results"
        self.description = "1b) Automated Pipe Sizing | MIKE1D Results"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        result_field = arcpy.Parameter(
            displayName="Field to assign diameter to",
            name="result_field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        result_layer = arcpy.Parameter(
            displayName="Result Layer with Max Q",
            name="result_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        result_layer_field = arcpy.Parameter(
            displayName="Field with Max_Q",
            name="result_layer_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        slopeOverwrite = arcpy.Parameter(
            displayName="Overwrite slope of pipes to:",
            name="slopeOverwrite",
            datatype="double",
            parameterType="optional",
            direction="Input")
        slopeOverwrite.category = "Additional settings"

        keep_largest_diameter = arcpy.Parameter(
            displayName="Only change diameter if it's greater than existing",
            name="keep_largest_diameter",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        keep_largest_diameter.category = "Additional settings"
        keep_largest_diameter.value = False

        change_material = arcpy.Parameter(
            displayName="Change material (Plastic if less than 500 mm, concrete if greater than or equal to 500 mm",
            name="change_material",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        change_material.category = "Additional settings"
        change_material.value = True

        parameters = [pipe_layer, result_field, result_layer, result_layer_field, slopeOverwrite, keep_largest_diameter, change_material]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  # optional
        pipe_layer = parameters[0].ValueAsText
        fields = [f.name for f in arcpy.Describe(pipe_layer).fields]
        if pipe_layer:
            parameters[1].filter.list = [f.name for f in arcpy.Describe(pipe_layer).fields]
        if "diameter" in [f.lower() for f in fields]:
            parameters[1].Value = "Diameter"

        result_layer = parameters[2].ValueAsText
        if result_layer:
            parameters[3].filter.list = [f.name for f in arcpy.Describe(result_layer).fields]
            if "maxq" in [f.name.lower() for f in arcpy.Describe(result_layer).fields]:
                parameters[3].Value = "MaxQ"
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].ValueAsText
        MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "").replace("!delete!","")
        result_field = parameters[1].ValueAsText
        result_layer = parameters[2].ValueAsText
        result_layer_field = parameters[3].ValueAsText
        slopeOverwrite = parameters[4].Value
        keep_largest_diameter = parameters[5].Value
        change_material = parameters[6].Value

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)
        #
        # config_folder = os.path.join(MIKE_folder, "Config")
        # if not os.path.exists(config_folder):
        #     os.mkdir(config_folder)
        # config_file = os.path.join(config_folder, os.path.splitext(os.path.basename(MU_database))[0] + ".ini")
        #
        # config = Config(config_file)
        # config.write(parameters)

        MIKE_gdb = os.path.join(MIKE_folder, os.path.splitext(os.path.basename(MU_database))[0])
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
                MIKE_gdb = os.path.join(MIKE_folder,
                                        "%s_%d" % (os.path.splitext(os.path.basename(MU_database))[0], dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        arcpy.SetProgressorLabel("Preparing")
        selected_pipes = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["muid"])]

        is_sqlite = True if ".sqlite" in MU_database else False

        msm_Link = os.path.join(MU_database, "msm_Link")

        arcpy.SetProgressorLabel("Calculating Pipe Dimensions")

        def addField(shapefile, field_name, datatype):
            i = 1
            while field_name in [f.name for f in arcpy.Describe(shapefile).fields]:
                field_name = "%s_%d" % (field_name, i)
            arcpy.AddField_management(shapefile, field_name, datatype)
            return field_name

        peak_discharge = {row[0]:row[1]*1e3 for row in arcpy.da.SearchCursor(result_layer, ["MUID", result_layer_field])}
        # arcpy.AddMessage(peak_discharge)

        if is_sqlite:
            # arcpy.AddMessage(MU_database)
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                with arcpy.da.SearchCursor(msm_Link,
                                           ["Slope" if is_sqlite else "Slope_C", "Diameter", "MaterialID", "MUID",
                                            "SHAPE@",
                                            "NetTypeNo",
                                            "enabled"],
                                           where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                    for row_i, row in enumerate(cursor):
                        arcpy.SetProgressorPosition(row_i)
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[2].lower() and not "beton" in row[
                                2].lower() else diameters_concrete
                        slope = slopeOverwrite if slopeOverwrite else row[0] * 1e-2
                        # if writeDischargeInstead:
                        #     diameter = # ins_row = (row[4], row[3], peak_discharge[msm_Link_Network.links[row[3]].fromnode], row[2], row[0], row[5], row[6])
                        #     ins_cursor.insertRow(ins_row)
                        QFull = 0
                        Di = -1
                        # arcpy.AddMessage(peak_discharge)
                        if peak_discharge[row[3]] == 0:
                            Di = 0
                        else:
                            while QFull is not None and QFull * 1e3 < peak_discharge[row[3]] and Di + 1 < len(D):
                                Di += 1
                                QFull = ColebrookWhite.QFull(D[Di] / 1e3, slope, row[2])
                        diameter = D[Di] / 1.0e3

                        material = row[2]
                        if keep_largest_diameter and diameter <= row[1]:
                            diameter = row[1]
                        else:
                            old_material = row[2]
                            if change_material:
                                material = "Concrete (Normal)" if diameter > 0.45 else "Plastic"
                            def materialName(material_id):
                                if "pl" in material_id.lower():
                                    material_name = "pl"
                                elif "concrete" in material_id.lower():
                                    material_name = "bt"
                                else:
                                    material_name = ""
                                return material_name

                            arcpy.AddMessage("Changed link '%s': %d%s → %d%s" % (row[3], row[1] * 1e3, materialName(old_material), D[Di], materialName(material)))
                            # arcpy.AddMessage("UPDATE msm_Link SET Diameter = %1.3f, SET MaterialID = %s WHERE MUID = %s" % (diameter, material, row[3]))
                            update_cursor.execute(
                                "UPDATE msm_Link SET Diameter = %1.3f, MaterialID = '%s' WHERE MUID = '%s'" % (
                                    diameter, material, row[3]))

        else:
            try:
                edit = arcpy.da.Editor(MU_database)
                edit.startEditing(False, True)
                edit.startOperation()

                fields = ["Slope_C", "MaterialID", "MUID", "Diameter"]
                if result_field and result_field not in fields:
                    result_field_index = len(fields) - 1
                    fields.append(result_field)

                arcpy.AddMessage(fields)
                with arcpy.da.UpdateCursor(msm_Link, fields,
                                           where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                    for row_i, row in enumerate(cursor):
                        # diameter_old = row[4]
                        arcpy.SetProgressorPosition(row_i)
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[1].lower() and not "beton" in row[
                                1].lower() else diameters_concrete
                        slope = slopeOverwrite if slopeOverwrite else max(row[0] * 1e-2, 5e-3)
                        QFull = 0
                        Di = -1

                        if peak_discharge[row[2]] == 0:
                            Di = 0
                        else:
                            while QFull is not None and QFull * 1e3 < peak_discharge[row[2]] and Di + 1 < len(D):
                                Di += 1
                                QFull = ColebrookWhite.QFull(D[Di] / 1e3, slope, row[1])

                        diameter = D[Di]
                        if keep_largest_diameter and diameter / 1.0e3 <= row[3]:
                            diameter = row[3]
                        else:
                            arcpy.AddMessage(
                                "Changed %s from %d to %d" % (row[2], row[3] * 1e3 if row[3] else 0, D[Di]))
                            row[-1] = D[Di] / 1.0e3
                            if change_material:
                                row[1] = "Concrete (Normal)" if row[-1] > 0.45 else "Plastic"
                            # if diameter_old != row[1]:
                            # arcpy.AddMessage("Changed diameter from %1.2f to %1.2f for pipe %s" % (diameter_old, row[1], row[3]))
                        try:
                            cursor.updateRow(row)
                        except Exception as e:
                            arcpy.AddWarning("Could not update row:")
                            arcpy.AddWarning(row)
                edit.stopOperation()
                edit.stopEditing(True)
            except Exception as e:
                arcpy.AddError(row)
                arcpy.AddError(traceback.format_exc())
                raise (e)
        # edit.stopOperation()
        # edit.stopEditing(True)
        # import pickle
        # import saveParameters
        # save_parameters = saveParameters.Parameters(parameters)
        # arcpy.AddMessage(save_parameters.parameters)
        # for parameter in save_parameters.parameters:
        #     arcpy.AddMessage((type(parameter.Value), type(parameter.ValueAsText)))
        # with open(r"C:\Papirkurv\Parameters",'w') as f:
        #     pickle.dump(save_parameters.parameters, f)
        return

class upgradeDimensions(object):
    def __init__(self):
        self.label       = "2b) Upgrade dimensions"
        self.description = "2b) Upgrade dimensions"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        change_material = arcpy.Parameter(
            displayName="Change material (Plastic if less than 500 mm, concrete if greater than or equal to 500 mm",
            name="change_material",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        change_material.value = "true"
        
        parameters = [pipe_layer, change_material]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            if not parameters[0].value:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                        and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and "diameter" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                         and lyr.visible][0]
                if links:
                    parameters[0].value = links

        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        change_material = parameters[1].Value
        MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "")
        MU_database = MU_database.replace("!delete!", "")
        is_sqlite = True if ".sqlite" in MU_database else False

        MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer,["MUID"])]
        if len(MUIDs) == len([row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).CatalogPath,["MUID"])]):
            userquery = pythonaddins.MessageBox("Change dimension of %d pipes?" % (len(MUIDs)), "Confirm Assignment", 4)
            if not userquery == "Yes":
                return

        if is_sqlite:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                with arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath.replace("!delete!",""), ["MUID", "Diameter", "MaterialID"],
                                           where_clause="MUID IN ('%s')" % ("', '".join(MUIDs))) as cursor:
                    for row in cursor:
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[2].lower() and not "beton" in row[
                                2].lower() else diameters_concrete
                        D = np.array(D)
                        oldDiameter = row[1] * 1e3
                        diameter = D[np.where(row[1] * 1e3 < D)[0][0]] / 1e3
                        material = row[2]
                        if change_material:
                            material = "Concrete (Normal)" if diameter > 0.45 else "Plastic"
                        update_cursor.execute("UPDATE msm_Link SET Diameter = %1.3f, MaterialID = '%s' WHERE MUID = '%s'" % (diameter, material, row[0]))
                        arcpy.AddMessage("Upgraded pipe %s from %d to %d" % (row[0], oldDiameter, diameter*1e3))
                connection.commit()
        else:
            edit = arcpy.da.Editor(os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)))
            edit.startEditing(False, True)
            edit.startOperation()

            with arcpy.da.UpdateCursor(arcpy.Describe(pipe_layer).catalogPath,["MUID", "Diameter", "MaterialID"], where_clause = "MUID IN ('%s')" % ("', '".join(MUIDs))) as cursor:
                for row in cursor:
                    if change_material:
                        D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                            diameter for diameter in diameters_concrete if diameter > 450]
                    else:
                        D = diameters_plastic if not "concrete" in row[2].lower() and not "beton" in row[
                            2].lower() else diameters_concrete
                    D = np.array(D)
                    oldDiameter = row[1]*1e3
                    row[1] = D[np.where(row[1]*1e3<D)[0][0]]/1e3
                    if change_material:
                        row[2] = "Concrete (Normal)" if row[1]>0.45 else "Plastic"
                    cursor.updateRow(row)
                    arcpy.AddMessage("Upgraded pipe %s from %d to %d" % (row[0],oldDiameter,row[1]*1e3))

            edit.stopOperation()
            edit.stopEditing(True)
        return

class downgradeDimensions(object):
    def __init__(self):
        self.label       = "2a) Downgrade dimensions"
        self.description = "2a) Downgrade dimensions"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        change_material = arcpy.Parameter(
            displayName="Change material (Plastic if less than 500 mm, concrete if greater than or equal to 500 mm",
            name="change_material",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        change_material.value = "true"
        # pipe_layer.filter.list = ["Polyline"] # does not work for some reason

        parameters = [pipe_layer, change_material]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        # mxd = arcpy.mapping.MapDocument("CURRENT")  
        # df = arcpy.mapping.ListDataFrames(mxd)[0]  
        # workspaces = set()
        # for lyr in arcpy.mapping.ListLayers(mxd, df):
            # if lyr.supports("workspacepath"):
                # workspaces.add(lyr.workspacePath)
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            if not parameters[0].value:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                        and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and "diameter" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                        and lyr.visible][0]
                if links:
                    parameters[0].value = links
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        change_material = parameters[1].Value
        MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "")
        MU_database = MU_database.replace("!delete!", "")
        is_sqlite = True if ".sqlite" in MU_database else False

        MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer,["MUID"])]
        if len(MUIDs) == len([row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).CatalogPath,["MUID"])]):
            userquery = pythonaddins.MessageBox("Change dimension of %d pipes?" % (len(MUIDs)), "Confirm Assignment", 4)
            if not userquery == "Yes":
                return

        if is_sqlite:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                D_plastic = np.array(diameters_plastic)
                D_concrete = np.array(diameters_concrete)
                with arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath.replace("!delete!",""), ["MUID", "Diameter", "MaterialID"],
                                           where_clause="MUID IN ('%s')" % ("', '".join(MUIDs))) as cursor:
                    for row in cursor:
                        if change_material:
                            D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                                diameter for diameter in diameters_concrete if diameter > 450]
                        else:
                            D = diameters_plastic if not "concrete" in row[2].lower() and not "beton" in row[
                                2].lower() else diameters_concrete
                        D = np.array(D)
                        oldDiameter = row[1] * 1e3
                        diameter = D[np.where(row[1]*1e3>D)[0][-1]]/1e3
                        material = row[2]
                        if change_material:
                            material = "Concrete (Normal)" if diameter > 0.45 else "Plastic"
                        arcpy.AddMessage("UPDATE msm_Link SET Diameter = %1.3f, MaterialID = '%s' WHERE MUID = '%s'" % (diameter, material, row[0]))
                        update_cursor.execute(
                            "UPDATE msm_Link SET Diameter = %1.3f, MaterialID = '%s' WHERE MUID = '%s'" % (diameter, material, row[0]))
                        arcpy.AddMessage("Downgraded pipe %s from %d to %d" % (row[0], oldDiameter, diameter*1e3))
                connection.commit()
        else:
            edit = arcpy.da.Editor(os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)))
            edit.startEditing(False, True)
            edit.startOperation()
            D_plastic = np.array(diameters_plastic)
            D_concrete = np.array(diameters_concrete)
            with arcpy.da.UpdateCursor(arcpy.Describe(pipe_layer).catalogPath,["MUID", "Diameter", "MaterialID"], where_clause = "MUID IN ('%s')" % ("', '".join(MUIDs))) as cursor:
                for row in cursor:
                    if change_material:
                        D = [diameter for diameter in diameters_plastic if diameter < 450] + [
                            diameter for diameter in diameters_concrete if diameter > 450]
                    else:
                        D = diameters_plastic if not "concrete" in row[2].lower() and not "beton" in row[
                            2].lower() else diameters_concrete
                    D = np.array(D)
                    oldDiameter = row[1]*1e3
                    row[1] = D[np.where(row[1]*1e3>D)[0][-1]]/1e3
                    if change_material:
                        row[2] = "Concrete (Normal)" if row[1]>0.45 else "Plastic"
                    cursor.updateRow(row)
                    arcpy.AddMessage("Downgraded pipe %s from %d to %d" % (row[0],oldDiameter,row[1]*1e3))

            edit.stopOperation()
            edit.stopEditing(True)
        return

class CopyDiameter(object):
    def __init__(self):
        self.label       = "6) Copy field value from layer to layer"
        self.description = "6) Copy field value from layer to layer"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        reference_feature_layer = arcpy.Parameter(
            displayName="Reference feature layer",
            name="reference_feature_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        target_feature_layer = arcpy.Parameter(
            displayName="Target feature layer",
            name="target_feature_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        copy_field = arcpy.Parameter(
            displayName="Field to copy",
            name="result_field",
            datatype="GPString",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        copy_field.filter.type = "ValueList"

        match_by = arcpy.Parameter(
            displayName="Match Feature Classes by",
            name="match_by",
            datatype="GPString",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        match_by.filter.type = "ValueList"
        match_by.filter.list = ["SHAPE@", "OBJECTID", "MUID", "FROMNODE-TONODE"]
        match_by.value = "SHAPE@"

        parameters = [reference_feature_layer, target_feature_layer, copy_field, match_by]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        reference_feature_layer = parameters[0].Value
        target_feature_layer = parameters[1].Value
        copy_field = parameters[2].ValueAsText
        match_by = parameters[3].ValueAsText
    
        def changeShapeFieldname(fields):
                for i in range(len(fields)):
                    if fields[i] == "SHAPE":
                        fields[i] = "SHAPE@"
                return fields

        if parameters[0].altered and not parameters[0].hasBeenValidated:
            reference_feature_layer = parameters[0].ValueAsText

            parameters[2].filter.list = changeShapeFieldname([f.name for f in arcpy.Describe(reference_feature_layer).fields])

        # if "diameter" in [field.lower() for field in parameters[2].filter.list] and not parameters[2].hasBeenValidated:
        #     parameters[2].value = "Diameter"
        #     if "MaterialID" in parameters[2].filter.list:
        #         parameters[2].value = ["Diameter", "MaterialID"]
        
        if parameters[1].ValueAsText:
            MU_database = os.path.dirname(os.path.dirname(arcpy.Describe(target_feature_layer).catalogPath))
            
            is_sqlite = True if ".sqlite" in MU_database else False
            if is_sqlite:
                parameters[3].value = "MUID"

        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        reference_feature_layer = parameters[0].Value
        target_feature_layer = parameters[1].Value
        copy_field = parameters[2].ValueAsText.split(";")
        match_by = parameters[3].ValueAsText
        MU_database = os.path.dirname(arcpy.Describe(target_feature_layer).catalogPath).replace("mu_Geometry","").replace("!delete!","")
        reference_MU_database = os.path.dirname(arcpy.Describe(reference_feature_layer).catalogPath).replace("mu_Geometry","").replace("!delete!","")
        target_MU_database = os.path.dirname(arcpy.Describe(target_feature_layer).catalogPath).replace("mu_Geometry",
                                                                                                          "").replace("!delete!","")
        
        # arcpy.AddMessage(MU_database)
        
        is_sqlite = True if ".sqlite" in MU_database else False
        arcpy.AddMessage(match_by)
        match_by = match_by.lower().replace("shape@", "shape")
        if is_sqlite:
            if "shape" in [field.lower() for field in copy_field]:
                arcpy.AddError("Copy field %s is not supported for sqlite" % (copy_field))
            if "shape" in match_by.lower():
                arcpy.AddError("Match field %s is not supported for sqlite" % (match_by))
        
        field_for_where_clause = "objectid" if not is_sqlite else "muid"
        check_other_layer = True
        arcpy.AddMessage("Confirm Query - Might be hidden behind window!")
        try:
            if arcpy.Describe(reference_feature_layer).fidSet:
                count = int(arcpy.GetCount_management(reference_feature_layer).getOutput(0))
                if arcgis_pro:
                    userquery = confirm_assignment(
                        "Change for %d features (selected in reference layer)?" % (count), "Confirm Assignment", 4)
                else:
                    userquery = pythonaddins.MessageBox("Change for %d features (selected in reference layer)?" % (count), "Confirm Assignment", 4)
                if not userquery == "Yes":
                    pass
                else:
                    if field_for_where_clause.lower() == "MUID".lower():
                        reference_where_clause = "%s IN ('%s')" % (field_for_where_clause, "', '".join([str(row[0]) for row in arcpy.da.SearchCursor(reference_feature_layer, [field_for_where_clause])])) 
                    else:
                        reference_where_clause = "%s IN (%s)" % (field_for_where_clause, ", ".join([str(row[0]) for row in arcpy.da.SearchCursor(reference_feature_layer, [field_for_where_clause])])) 
                    target_where_clause = ""
                    check_other_layer = False
        except Exception as e:
            arcpy.AddError(e)
            pass
       
        # arcpy.AddMessage((check_other_layer and arcpy.Describe(target_feature_layer).fidSet))
        try:
            if check_other_layer and arcpy.Describe(target_feature_layer).fidSet:
                count = int(arcpy.GetCount_management(target_feature_layer).getOutput(0))
                if arcgis_pro:
                    userquery = confirm_assignment("Change for %d features (selected in target layer)?" % (count),
                                                   "Confirm Assignment", 4)
                else:
                    userquery = pythonaddins.MessageBox("Change for %d features (selected in target layer)?" % (count), "Confirm Assignment", 4)
                if not userquery == "Yes": 
                    arcpy.AddMessage("Cancelled both user queries")
                    return
                else:
                    if field_for_where_clause.lower() == "MUID".lower():
                        target_where_clause = "%s IN ('%s')" % (field_for_where_clause, "', '".join([str(row[0]) for row in arcpy.da.SearchCursor(target_feature_layer, [field_for_where_clause])]))
                    else:
                        target_where_clause = "%s IN (%s)" % (field_for_where_clause, ", ".join([str(row[0]) for row in arcpy.da.SearchCursor(target_feature_layer, [field_for_where_clause])])) 
                    reference_where_clause = ""
            else:
                if check_other_layer:
                    arcpy.AddMessage("Nothing selected")
                    return
        except Exception as e:
            arcpy.AddError(e)
            return  
            
        # arcpy.AddMessage(target_where_clause)
        # arcpy.AddMessage(reference_where_clause)
                
        class Reference():
            pass
            
        def changeShapeFieldname(fields):
            for i in range(len(fields)):
                if fields[i].lower() == "shape":
                    fields[i] = "shape@"

            if "shape@" not in fields:
                fields.append("shape@")
            return fields

        references = []
        fields = changeShapeFieldname([field.name.lower() for field in arcpy.ListFields(reference_feature_layer)])
        if "esri_oid" in fields:
            fields.remove("esri_oid")
        # arcpy.AddMessage(fields)

        arcpy.AddMessage(reference_where_clause)
        arcpy.AddMessage(arcpy.Describe(reference_feature_layer).catalogPath)
        arcpy.AddMessage(fields)

        with arcpy.da.SearchCursor(arcpy.Describe(reference_feature_layer).catalogPath, fields, where_clause = reference_where_clause) as cursor:
            for row in cursor:
                reference = Reference()
                for field_i, field in enumerate(fields):
                    setattr(reference, field.lower().replace("shape@","shape"), row[field_i])
                
                references.append(reference)
        # arcpy.AddMessage([reference.muid for reference in references])
        # arcpy.AddMessage([reference.diameter for reference in references])
        
        fields = changeShapeFieldname([field.name.lower() for field in arcpy.ListFields(target_feature_layer)])
        if "ESRI_OID" in fields:
            fields.remove("ESRI_OID")

        # arcpy.AddMessage(fields)
        if match_by.lower() != "FROMNODE-TONODE".lower():
            # arcpy
            match_by_field_i = [field_i for field_i, field in enumerate(fields) if field.lower() == match_by.lower().replace("shape","shape@")][0]
        # arcpy.AddMessage(fields)
        # arcpy.AddMessage(match_by_field_i)
        
        target_where_clause = target_where_clause if target_where_clause else ""
        # arcpy.AddMessage(reference_where_clause)

        if match_by.lower() == "FROMNODE-TONODE".lower():
            import networker
            reference_network = networker.NetworkLinks(reference_MU_database, map_only="link")
            target_network = networker.NetworkLinks(target_MU_database, map_only="link")

        if is_sqlite:
            # arcpy.AddMessage(MU_database)
            with sqlite3.connect(
                        MU_database) as connection:
                update_cursor = connection.cursor()
                layer_name = os.path.basename(arcpy.Describe(target_feature_layer).catalogPath).replace("main.","")
                arcpy.AddMessage(copy_field)
                # copy_field = copy_field if type(copy_field) is str else copy_field.split(";")
                # arcpy.AddMessage(copy_field)
                # arcpy.AddMessage(target_where_clause)
                # arcpy.AddMessage([row[0] for row in arcpy.da.SearchCursor(target_feature_layer, ["MUID"], target_where_clause)])

                # arcpy.AddMessage("BOB")
                # arcpy.AddMessage(references[0].muid)
                for MUID in [row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(target_feature_layer).catalogPath, ["MUID"], target_where_clause)]:
                    if match_by.lower() == "FROMNODE-TONODE".lower():
                        if MUID in target_network.links:
                            muid_field_i = [i for i, field in enumerate(fields) if field.lower() == "muid"][0]
                            # arcpy.AddMessage(row[muid_field_i])
                            # arcpy.AddMessage(MUID in target_network.links)
                            # arcpy.AddMessage(reference.muid in reference_network.links)
                            # arcpy.AddMessage((reference_network.links[reference.muid].fromnode, reference_network.links[reference.muid].tonode))
                            # arcpy.AddMessage((reference_network.links[reference.muid].fromnode,
                            #                   reference_network.links[reference.muid].tonode))
                            try:
                                match = [reference for reference in references if reference_network.links[reference.muid].fromnode.lower() == target_network.links[MUID].fromnode.lower() and reference_network.links[reference.muid].tonode.lower() == target_network.links[MUID].tonode.lower()
                                         if reference.muid in reference_network.links and MUID in target_network.links]
                            except Exception as e:
                                for reference in references:
                                    if not reference_network.links[reference.muid].fromnode or not reference_network.links[reference.muid].tonode:
                                        arcpy.AddError("Could not find fromnode or tonode for link %s" % reference.muid)
                                arcpy.AddError("Could not match fromnode / tonode for link %s" % (MUID))

                            # arcpy.AddMessage(match)
                        else:
                            continue
                        # match = [reference for reference in references if getattr(reference, "muid") == match_MUID][0]
                    else:
                        pass
                        # arcpy.AddMessage(MUID)

                        match = [reference for reference in references if getattr(reference, "muid") == MUID]
                    if match:
                        for field in copy_field:
                            field_value = getattr(match[0], field.lower())
                            # arcpy.AddMessage(type(field_value) is str or type(field_value) is unicode)

                            old_field_value = update_cursor.execute("SELECT %s FROM %s WHERE MUID = '%s'" % (field, layer_name, MUID)).fetchone()[0]
                            sql_expression = "UPDATE %s SET %s = %s WHERE MUID = '%s'" % (layer_name, field,
                                                "'%s'" % (field_value) if type(field_value) is str or (arcgis_pro or type(field_value) is unicode) else "%s" % (field_value),
                                                MUID)
                            arcpy.AddMessage(sql_expression)
                            arcpy.AddMessage(
                                        "Changed %s field %s from %s to %s" % (MUID, field, old_field_value, field_value))
                            try:
                                update_cursor.execute(sql_expression)
                            except Exception as e:
                                arcpy.AddMessage(sql_expression)
                                raise(e)
        else:
            edit = arcpy.da.Editor(MU_database)
            edit.startEditing(False, True)
            edit.startOperation()
            with arcpy.da.UpdateCursor(arcpy.Describe(target_feature_layer).catalogPath, fields, where_clause = target_where_clause) as cursor:
                for row in cursor:
                    # arcpy.AddMessage((fields, row))
                    if match_by.lower() == "FROMNODE-TONODE".lower():
                        muid_field_i = [i for i,field in enumerate(fields) if field.lower() == "muid"]
                        match_MUID = [link.muid for link in reference_links if link.fromnode.lower() == target_network.links[row[muid_field_i]].fromnode and link.tonode.lower() == target_network.links[row[muid_field_i]].tonode]
                        arcpy.AddMessage("FEATURE NOT SUPPORTED YET")
                        return
                        # match = [reference for reference in references if reference_network.links[

                    else:
                        match = [reference for reference in references if getattr(reference, match_by.lower()) == row[match_by_field_i]]

                    # arcpy.AddMessage((getattr(reference, match_by.lower()), row[match_by_field_i]))
                    # arcpy.AddMessage([getattr(reference, match_by.lower()) for reference in references])
                    if match:
                        reference = match[0]
                        for field_i, field in enumerate(fields):
                            # arcpy.AddMessage(copy_field)
                            if field.lower() in [f.lower() for f in copy_field]:
                                if row[field_i] != getattr(reference, field.lower()):
                                    arcpy.AddMessage(
                                        "Changed %s field %s from %s to %s" % (row[match_by_field_i], field, row[field_i], getattr(reference, field.lower())))
                                    row[field_i] = getattr(reference, field.lower())
                            elif field == "SHAPE":
                                shape = deepcopy(row[field_i])
                                row[field_i] = shape
                        # arcpy.AddMessage(row)
                        cursor.updateRow(row)


            edit.stopOperation()
            edit.stopEditing(True)
        return

class InterpolateInvertLevels(object):
    def __init__(self):
        self.label       = "3) Invert Level Interpolation"
        self.description = "3) Invert Level Interpolation"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        usePipeElevations = arcpy.Parameter(
            displayName="Include UpLevel and DwLevel",
            name="usePipeElevations",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
            
        fixed_slope = arcpy.Parameter(
            displayName="Use Fixed Slope instead of Interpolation:",
            name="fixed_slope",
            datatype="double",
            parameterType="optional",
            direction="Input")
        fixed_slope.category = "Use Fixed Slope"
        
        use_slope_from_upstream = arcpy.Parameter(
            displayName="Use Upstream Node as Fix Point for setting Fixed Slope",
            name="use_slope_from_upstream",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        use_slope_from_upstream.category = "Use Fixed Slope"
        use_slope_from_upstream.value = True
            
        use_slope_from_downstream = arcpy.Parameter(
            displayName="Use Downstream Node as Fix Point for setting Fixed Slope",
            name="use_slope_from_downstream",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")
        use_slope_from_downstream.category = "Use Fixed Slope"

        parameters = [pipe_layer, usePipeElevations, fixed_slope, use_slope_from_upstream, use_slope_from_downstream]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        if not parameters[0].value:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                    and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                    and lyr.visible][0]
            if links:
                parameters[0].value = links
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        usePipeElevations = parameters[1].Value
        fixed_slope = parameters[2].Value
        use_slope_from_upstream = parameters[3].Value
        use_slope_from_downstream = parameters[4].Value
        
        MU_database = (os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)) if ".mdb" in arcpy.Describe(pipe_layer).catalogPath else
                        os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)).replace("!delete!","")
        is_sqlite = True if ".sqlite" in MU_database else False        

        links_MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer,["MUID"])]
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")
        end_node_critical = None

        if len(links_MUIDs) == len(
                [row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).CatalogPath, ["MUID"])]):
            userquery = pythonaddins.MessageBox("Interpolate invert slope for %d pipes?" % (len(links_MUIDs)), "Confirm Assignment", 4)
            if not userquery == "Yes":
                return
        
        if not len(links_MUIDs) == len(set(links_MUIDs)):
            arcpy.AddError("Error: There's two or more pipes with identical names.")
            duplicates = []
            for MUID in links_MUIDs:
                if links_MUIDs.count(MUID)>1 and MUID not in duplicates:
                    arcpy.AddError("Pipe %s" % (MUID))
                    duplicates.append(MUID)
            raise(Exception("Cancelling toolbox"))

        msm_Link_Network = networker.NetworkLinks(MU_database, map_only="link", filter_sql_query = "MUID IN ('%s')" % ("', '".join(links_MUIDs)))
        
        tonodes = [msm_Link_Network.links[MUID].tonode for MUID in links_MUIDs]
        fromnodes = [msm_Link_Network.links[MUID].fromnode for MUID in links_MUIDs]
        try:
            start_node = [fromnode for fromnode in fromnodes if fromnode not in tonodes][0]
            end_node = [tonode for tonode in tonodes if tonode not in fromnodes][0]
        except Exception as e:
            for MUID in links_MUIDs:
                arcpy.AddMessage("%s: %s-%s" % (MUID, msm_Link_Network.links[MUID].fromnode, msm_Link_Network.links[MUID].tonode))
            arcpy.AddWarning(("fromnodes:", fromnodes))
            arcpy.AddWarning(("tonodes:", tonodes))
            raise(e)

        invert_levels = {row[0]: row[1] for row in
                         arcpy.da.SearchCursor(msm_Node,
                                               ["MUID", "InvertLevel"],
                                               where_clause="MUID IN ('%s', '%s')" % (start_node, end_node))}
        
        if usePipeElevations:
            last_link = [MUID for MUID in links_MUIDs if msm_Link_Network.links[MUID].tonode == end_node][0]
            dwlevel = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["dwlevel"], where_clause = "MUID = '%s'" % (last_link))][0]
            
            first_link = [MUID for MUID in links_MUIDs if msm_Link_Network.links[MUID].fromnode == start_node][0]
            uplevel = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["dwlevel"], where_clause = "MUID = '%s'" % (first_link))][0]
            if uplevel:
                invert_levels[start_node] = uplevel
            if dwlevel:
                invert_levels[end_node] = dwlevel

        network = nx.DiGraph()
        for link in msm_Link_Network.links.values():
            network.add_edge(link.fromnode, link.tonode,
                             weight=link.length)

        path_nodes = nx.bellman_ford_path(network, start_node, end_node, weight="weight")
        lengths = np.zeros(len(path_nodes) - 1, dtype=np.float)
        for i in range(1, len(path_nodes)):
            lengths[i - 1] = network.edges[path_nodes[i - 1], path_nodes[i]]["weight"]

        try:
            if fixed_slope is not None:
                slope = fixed_slope
                invert_levels[end_node] = invert_levels[start_node] - fixed_slope * np.sum(lengths) if use_slope_from_upstream else invert_levels[end_node]
            else:
                slope = (invert_levels[start_node] - invert_levels[end_node]) / np.sum(lengths)
        except Exception as e:
            arcpy.AddError(("start_node", start_node))
            arcpy.AddError(("end_node", end_node))
            arcpy.AddError(("tonodes", tonodes))
            arcpy.AddError(("fromnodes", fromnodes))
            arcpy.AddError(("invert_levels[start_node]", invert_levels[start_node]))
            arcpy.AddError(("invert_levels[end_node]", invert_levels[end_node]))
            arcpy.AddError(("np.sum(lengths)", np.sum(lengths)))
            raise(e)
        arcpy.AddMessage("Assuming slope %d o/oo" % (slope*1e3))

        if is_sqlite:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                arcpy.AddMessage(MU_database)
                with arcpy.da.SearchCursor(msm_Node, ["MUID", "InvertLevel", "GroundLevel"],
                                           where_clause="MUID IN ('%s')" % ("', '".join(path_nodes))) as cursor:
                    for row in cursor:
                        if (row[0] != end_node 
                            or (fixed_slope is not None and use_slope_from_downstream and row[0] != end_node)
                            or (fixed_slope is not None and use_slope_from_upstream and row[0] != start_node)):
                            total_length = nx.bellman_ford_path_length(network, row[0], end_node, weight="weight")
                            
                            new_invert_level = round(invert_levels[end_node] + total_length * slope,2)
                            if new_invert_level != row[1]:
                                update_cursor.execute("UPDATE msm_Node SET InvertLevel = %1.2f WHERE MUID = '%s'" % (new_invert_level, row[0]))
                                arcpy.AddMessage(
                                    "Changed invert level of %s from %1.2f to %1.2f" % (row[0], row[1] if row[1] else 0, new_invert_level))
        else:
            edit = arcpy.da.Editor(MU_database)
            edit.startEditing(False, True)
            edit.startOperation()

            with arcpy.da.UpdateCursor(msm_Node, ["MUID", "InvertLevel", "GroundLevel"],
                                       where_clause="MUID IN ('%s')" % ("', '".join(path_nodes))) as cursor:
                for row in cursor:
                    if row[0] != end_node or (fixed_slope is not None and use_slope_from_downstream and row[0] != end_node) or (fixed_slope is not None and use_slope_from_upstream and row[0] != start_node):
                        total_length = nx.bellman_ford_path_length(network, row[0], end_node, weight="weight")

                        try:
                            new_invert_level = round(invert_levels[end_node] + total_length * slope,2)
                        except Exception as e:
                            arcpy.AddError((invert_levels[end_node], total_length, slope))
                            raise(e)
                        if new_invert_level != row[1]:
                            arcpy.AddMessage(
                                "Changed invert level of %s from %1.2f to %1.2f" % (row[0], row[1] if row[1] else 0, new_invert_level))
                            row[1] = new_invert_level
                            cursor.updateRow(row)

            edit.stopOperation()
            edit.stopEditing(True)
        return

class GetMinimumSlope(object):
    def __init__(self):
        self.label = "1c) Calculate Minimum Slope of Energy Gradient"
        self.description = "1c) Calculate Minimum Slope of Energy Gradient"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        end_node_critical = arcpy.Parameter(
            displayName="Critical Level at End Node",
            name="end_node_critical",
            datatype="double",
            parameterType="Optional",
            direction="Input")

        parameters = [pipe_layer, end_node_critical]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            if not parameters[0].value:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                        and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and "diameter" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                        and lyr.visible][0]
                if links:
                    parameters[0].value = links
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        end_node_critical = parameters[1].Value
        MU_database = (os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)) if ".mdb" in arcpy.Describe(pipe_layer).catalogPath else
                        os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)).replace("!delete!","")
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")

        if arcgis_pro:
            mxd = arcpy.mp.ArcGISProject("CURRENT")
            df = mxd.listMaps()[0]
        else:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        def addLayer(layer_source, source, group = None, workspace_type = "ACCESS_WORKSPACE"):
            layer = apmapping.Layer(layer_source)
            if group:
                apmapping.AddLayerToGroup(df, group, layer, "BOTTOM")
            else:
                apmapping.AddLayer(df, layer, "TOP")
            updatelayer = apmapping.ListLayers(mxd, layer.name, df)[0]
            updatelayer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source)))

        links_MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]
        where_clause = "MUID IN ('%s')" % ("', '".join(links_MUIDs))
        msm_Link_Network = networker.NetworkLinks(MU_database, map_only="link", filter_sql_query = where_clause)

        tonodes = [msm_Link_Network.links[MUID].tonode for MUID in links_MUIDs]
        fromnodes = [msm_Link_Network.links[MUID].fromnode for MUID in links_MUIDs]
        start_nodes = [fromnode for fromnode in fromnodes if fromnode not in tonodes]
        end_node = [tonode for tonode in tonodes if tonode not in fromnodes][0]

        levels = {row[0]: [row[1], row[2]] for row in
                  arcpy.da.SearchCursor(msm_Node,
                                        ["MUID", "InvertLevel", "GroundLevel"],
                                        where_clause="MUID IN ('%s')" % ("', '".join(set(fromnodes + tonodes))))}

        with arcpy.da.SearchCursor(msm_Link, ["MUID", "Diameter"], where_clause = where_clause) as cursor:
            for row in cursor:
                msm_Link_Network.links[row[0]].diameter = row[1]

        # end_node_critical = (end_node_critical if
        #                      end_node_critical else levels[end_node])

        network = nx.DiGraph()
        for link in msm_Link_Network.links.values():
            network.add_edge(link.fromnode, link.tonode,
                             weight=link.length)

        critical_energy_gradient = {}
        for start_node in start_nodes:
            path_nodes = nx.bellman_ford_path(network, start_node, end_node, weight="weight")

            minimum_ground_slope = {}
            last_critical_energy_gradient = None
            debug = []
            for fromnode_i, fromnode in enumerate(path_nodes[:-1]):
                total_gradient = ((levels[fromnode][1] -
                                   (end_node_critical if end_node_critical else levels[end_node][0] +
                                                                                [link.diameter for link in
                                                                                 msm_Link_Network.links.values() if
                                                                                 link.tonode == end_node and
                                                                                 link.fromnode in path_nodes][0])) /
                                  nx.bellman_ford_path_length(network, fromnode, end_node, weight="weight"))

                for tonode in path_nodes[fromnode_i + 1:]:
                    gradient = ((levels[fromnode][1] - levels[tonode][0] -
                                 [link.diameter for link in msm_Link_Network.links.values() if link.tonode == tonode and
                                  link.fromnode in path_nodes][0]) /
                                nx.bellman_ford_path_length(network, fromnode, tonode, weight="weight"))

                    gradient = gradient if gradient < total_gradient else total_gradient
                    gradient = gradient if not last_critical_energy_gradient or gradient < last_critical_energy_gradient else last_critical_energy_gradient
                    critical_energy_gradient[
                        fromnode] = gradient if not fromnode in critical_energy_gradient or gradient < \
                                                critical_energy_gradient[fromnode] else critical_energy_gradient[
                        fromnode]
                last_critical_energy_gradient = critical_energy_gradient[fromnode]
        
        if ".sqlite" in MU_database:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                with arcpy.da.SearchCursor(msm_Link, ["MUID"],
                                       where_clause="MUID IN ('%s')" % ("', '".join(links_MUIDs))) as cursor:
                    for row in cursor:
                        fromnode = msm_Link_Network.links[row[0]].fromnode
                        if fromnode in critical_energy_gradient:
                            update_cursor.execute(
                                "UPDATE msm_Link SET Slope = %1.2f WHERE MUID = '%s'" % (
                                    critical_energy_gradient[fromnode] * 1e2, row[0]))
                        else:
                            arcpy.AddWarning("%s not in critical_energy_gradient" % (fromnode))
            # msm_Link_result = getAvailableFilename(arcpy.env.scratchGDB + "\msm_Link")
            # arcpy.Select_analysis(msm_Link, msm_Link_result, where_clause="MUID IN ('%s')" % ("', '".join(links_MUIDs)))
        else:
            msm_Link_result = msm_Link
            edit = arcpy.da.Editor(MU_database.replace("!delete!",""))
            edit.startEditing(False, True)
            edit.startOperation()
            with arcpy.da.UpdateCursor(msm_Link.replace("!delete!",""), ["MUID", "Slope_C" if not ".sqlite" in MU_database else "Slope"],
                                       where_clause="MUID IN ('%s')" % ("', '".join(links_MUIDs))) as cursor:
                for row in cursor:
                    fromnode = msm_Link_Network.links[row[0]].fromnode
                    if fromnode in critical_energy_gradient:
                        row[1] = critical_energy_gradient[fromnode] * 1e2
                        cursor.updateRow(row)
                    else:
                        arcpy.AddWarning("%s not in critical_energy_gradient" % (fromnode))
            edit.stopOperation()
            edit.stopEditing(True)
        # if ".sqlite" in MU_database:
        #     addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Links Dimensioned.lyr", msm_Link_result, workspace_type = "FILEGDB_WORKSPACE")
        # else:
        #     edit.stopOperation()
        #     edit.stopEditing(True)

        return

class setOutletLoss(object):
    def __init__(self):
        self.label       = "a) Change outlet loss of nodes"
        self.description = "a) Change outlet loss of nodes"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")


        parameters = [MU_database]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        # mxd = arcpy.mapping.MapDocument("CURRENT")  
        # df = arcpy.mapping.ListDataFrames(mxd)[0]  
        # workspaces = set()
        # for lyr in arcpy.mapping.ListLayers(mxd, df):
            # if lyr.supports("workspacepath"):
                # workspaces.add(lyr.workspacePath)

        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):

        MU_database = parameters[0].ValueAsText
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")

        msm_Node_saddle = set()
        try:
            with arcpy.da.SearchCursor(msm_Link, ["FromNode","ToNode","Diameter"]) as cursor:
                for row in cursor:
                    for node in [row[0],row[1]]:
                        if row[2] > 0.950:
                            msm_Node_saddle.add(node)
        except Exception as e:
            arcpy.AddError("Error: Try running Project Check Tool.")
            raise(e)


        edit = arcpy.da.Editor(MU_database)
        edit.startEditing(False, True)
        edit.startOperation()
        with arcpy.da.UpdateCursor(msm_Node, ["MUID", "LossParID"]) as cursor:
            for row in cursor:
                if row[0] in msm_Node_saddle and not row[1] == "No Cross Section Changes":
                    old_loss_par_id = row[1]
                    row[1] = "No Cross Section Changes"
                    cursor.updateRow(row)
                    arcpy.AddMessage("Changed Loss Par ID of node %s from %s to %s because connected pipe has diameter above 1000 mm" % (row[0], old_loss_par_id, row[1]))
                elif row[0] not in msm_Node_saddle and not row[1] == "Weighted Inlet Energy":
                    old_loss_par_id = row[1]
                    row[1] = "Weighted Inlet Energy"
                    cursor.updateRow(row)
                    arcpy.AddMessage("Changed Loss Par ID of node %s from %s to %s" % (row[0], old_loss_par_id, row[1]))

        edit.stopOperation()
        edit.stopEditing(True)
        # MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer,["MUID"])]

        # edit = arcpy.da.Editor(os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)))
        # edit.startEditing(False, True)
        # edit.startOperation()
        # D = np.array([102.4, 149, 188, 235, 278, 396, 470.8, 497, 588, 781, 885, 985, 1085, 1185, 1285, 1385, 1485, 1585, 2000, 2200, 2400, 2600])
        # with arcpy.da.UpdateCursor(arcpy.Describe(pipe_layer).catalogPath,["MUID", "Diameter"], where_clause = "MUID IN ('%s')" % ("', '".join(MUIDs))) as cursor:
            # for row in cursor:
                # oldDiameter = row[1]*1e3
                # row[1] = D[np.where(row[1]*1e3>D)[0][-1]]/1e3
                # cursor.updateRow(row)
                # arcpy.AddMessage("Downgraded pipe %s from %d to %d" % (row[0],oldDiameter,row[1]*1e3))

        # edit.stopOperation()
        # edit.stopEditing(True)
        return

# class internalDimensions(object):
    # def __init__(self):
        # self.label       = "Change diameter to internal dimensions for plastic pipes"
        # self.description = "Change diameter to internal dimensions for plastic pipes"
        # self.canRunInBackground = False

    # def getParameterInfo(self):
        # #Define parameter definitions

        # MU_database = arcpy.Parameter(
            # displayName="Mike Urban database",
            # name="database",
            # datatype="DEWorkspace",
            # parameterType="Required",
            # direction="Input")


        # parameters = [MU_database]
        # return parameters

    # def isLicensed(self):
        # return True

    # def updateParameters(self, parameters): #optional
        # # mxd = arcpy.mapping.MapDocument("CURRENT")
        # # df = arcpy.mapping.ListDataFrames(mxd)[0]
        # # workspaces = set()
        # # for lyr in arcpy.mapping.ListLayers(mxd, df):
            # # if lyr.supports("workspacepath"):
                # # workspaces.add(lyr.workspacePath)

        # return

    # def updateMessages(self, parameters): #optional
        # return

    # def execute(self, parameters, messages):

        # MU_database = parameters[0].ValueAsText
        # msm_Link = os.path.join(MU_database, "msm_Link")

        # externalDiameters = {0.200: 0.180, 0.250: 0.233, 0.300: 0.276, 0.400: 0.392, 0.500: 0.493, 0.600: 0.588, 0.800: 0.781, 1.000: 0.985, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600}
        # [180, 233, 276, 392, 493, 588, 781, 985, 1185, 1285, 1385, 1485, 1585, 2000, 2200, 2400, 2600]

        # with arcpy.da.UpdateCursor(msm_Link,["MUID","Diameter","MaterialID"]), where_clause = "MaterialID = 'Plastic' AND Diameter IN (0.160, 0.200, 0.250, 0.300, 0.400, 0.500, 0.600, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600)") as cursor:
            # for row in cursor:

        # return

class reverseChange(object):
    def __init__(self):
        self.label       = "2c) Revert changed dimensions"
        self.description = "2c) Revert changed dimensions"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        reverse_text = arcpy.Parameter(
            displayName="Output log of upgrade/downgrade dimensions result",
            name="reverse_text",
            datatype="String",
            parameterType="Required",
            direction="Input")

        parameters = [pipe_layer, reverse_text]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        # mxd = arcpy.mapping.MapDocument("CURRENT")  
        # df = arcpy.mapping.ListDataFrames(mxd)[0]  
        # workspaces = set()
        # for lyr in arcpy.mapping.ListLayers(mxd, df):
            # if lyr.supports("workspacepath"):
                # workspaces.add(lyr.workspacePath)

        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):

        pipe_layer = parameters[0].Value
        reverse_text = parameters[1].Value

        find_pipe_dimension = re.compile("pipe ([^ ]+) from ([\d\.]+)")
        results = find_pipe_dimension.findall(reverse_text)

        pipe_dimension_dictionairy = {}
        for result in results:
            pipe_dimension_dictionairy[result[0]] = float(result[1])/1e3


        edit = arcpy.da.Editor(os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)))
        edit.startEditing(False, True)
        edit.startOperation()
        with arcpy.da.UpdateCursor(pipe_layer, ["MUID","Diameter"], where_clause = "MUID IN ('%s')" % "', '".join(pipe_dimension_dictionairy.keys())) as cursor:
            for row in cursor:
                arcpy.AddMessage("Changing pipe %s from %d to %d" % (row[1], pipe_dimension_dictionairy[row[0]]))
                row[1] = pipe_dimension_dictionairy[row[0]]
                cursor.updateRow(row)
        edit.stopOperation()
        edit.stopEditing(True)
        return
        

class CalculateSlopeOfPipe(object):
    def __init__(self):
        self.label       = "b) Calculate Slope of Pipe"
        self.description = "b) Calculate Slope of Pipe"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
            multiValue = True)

        parameters = [pipe_layer]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            if not parameters[0].value:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                        and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                         and lyr.visible]
                if links:
                    parameters[0].value = ";".join(links)
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        pipe_layers = parameters[0].ValueAsText.split(";")

        for pipe_layer in pipe_layers:
            is_sqlite = ".sqlite" in arcpy.Describe(pipe_layer).catalogPath
            if is_sqlite:
                links_MUID = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]
                MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "")

                msm_Node = os.path.join(MU_database, "msm_Node")

                nodes_invert_level = {row[0]: row[1] for row in
                                      arcpy.da.SearchCursor(msm_Node, ["MUID", "InvertLevel"])}

                arcpy.SetProgressorLabel("Networking Database")
                network = networker.NetworkLinks(MU_database, map_only = "link", filter_sql_query = "MUID IN ('%s')" % ("', '".join(links_MUID)))
                arcpy.SetProgressor("step", "Calculating slope of pipes", 0, len(links_MUID), 1)

                with sqlite3.connect(MU_database) as connection:
                    update_cursor = connection.cursor()
                    with arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath, ["MUID", "UpLevel", "DwLevel"], where_clause = "MUID IN ('%s')" % ("', '".join(links_MUID))) as cursor:
                        for row_i, row in enumerate(cursor):
                            try:

                                arcpy.SetProgressorPosition(row_i)
                                uplevel = nodes_invert_level[network.links[row[0]].fromnode] if not row[1] else row[1]
                                dwlevel = nodes_invert_level[network.links[row[0]].tonode] if not row[2] else row[2]

                                length = network.links[row[0]].length
                                slope = (uplevel - dwlevel) / length * 1e2
                                update_cursor.execute("UPDATE msm_Link SET Slope = %1.6f WHERE MUID = '%s'" % (slope, row[0]))
                                # arcpy.AddMessage("UPDATE msm_Link SET Slope = %1.6f WHERE MUID = '%s'" % (slope, row[0]))
                            except Exception as e:
                                arcpy.AddError(traceback.format_exc())
                                arcpy.AddError(row)

                    connection.commit()
            else:
                links_OID = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["OID@"])]
                OID_fieldname = arcpy.Describe(pipe_layer).OIDFieldName

                MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry","")

                msm_Node = os.path.join(MU_database, "msm_Node")

                nodes_invert_level = {row[0]: row[1] for row in arcpy.da.SearchCursor(msm_Node, ["MUID", "InvertLevel"])}

                arcpy.SetProgressorLabel("Networking Database")
                network = networker.NetworkLinks(MU_database, map_only = "link", filter_sql_query = "%s IN (%s)" % (OID_fieldname, ', '.join([str(OID) for OID in links_OID])))

                edit = arcpy.da.Editor(MU_database)
                edit.startEditing(False, True)
                edit.startOperation()

                arcpy.SetProgressor("step", "Calculating slope of pipes", 0, len(links_OID), 1)
                with arcpy.da.UpdateCursor(arcpy.Describe(pipe_layer).catalogPath, ["MUID", "Slope_C", "UpLevel", "DwLevel"], where_clause = "%s IN (%s)" % (OID_fieldname, ', '.join([str(OID) for OID in links_OID]))) as cursor:
                    for row_i, row in enumerate(cursor):
                        try:
                            arcpy.SetProgressorPosition(row_i)
                            uplevel = nodes_invert_level[network.links[row[0]].fromnode] if not row[2] else row[2]
                            dwlevel = nodes_invert_level[network.links[row[0]].tonode] if not row[3] else row[3]
                            length = network.links[row[0]].length
                            slope = (uplevel-dwlevel)/length*1e2
                            row[1] = slope
                            # arcpy.AddMessage(row)
                            # arcpy.AddMessage((uplevel, dwlevel, length, slope))
                            cursor.updateRow(row)
                        except Exception as e:
                            arcpy.AddError(traceback.format_exc())
                            arcpy.AddError(row)
                edit.stopOperation()
                edit.stopEditing(True)
        return


class ResetUpLevelDwlevel(object):
    def __init__(self):
        self.label = "5) Set Uplevel and Dwlevel to NULL"
        self.description = "5) Set Uplevel and Dwlevel to NULL"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        if_equal_to_invert_level = arcpy.Parameter(
            displayName="Only do it for links where uplevel and dwlevel are already equal to Invert Level?",
            name="if_equal_to_invert_level",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")

        parameters = [pipe_layer, if_equal_to_invert_level]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  # optional
        if arcgis_pro:
            if not parameters[0].value:
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet() and arcpy.Describe(layer).shapeType == "Polyline":
                            parameters[0].value = layer.longName
                            break
                    except:
                        pass
        else:
            if not parameters[0].value:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                         lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                         and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]
                         and lyr.visible][0]
                if links:
                    parameters[0].value = links
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        if_equal_to_invert_level = parameters[1].Value

        is_sqlite = True if ".sqlite" in arcpy.Describe(pipe_layer).catalogPath else False

        if not is_sqlite:
            links_OID = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["OID@"])]
            if len(links_OID) == 0:
                links_OID = [row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath, ["OID@"])]
            OID_fieldname = arcpy.Describe(pipe_layer).OIDFieldName
            where_clause = "%s IN (%s)" % (OID_fieldname, ", ".join([str(OID) for OID in links_OID]))
        else:
            links_OID = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]
            if len(links_OID) == 0:
                links_OID = [row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath, ["MUID"])]
            OID_fieldname = "MUID"
            if if_equal_to_invert_level:
                where_clause = "%s IN ('%s') AND (UpLevel IS NOT NULL OR DwLevel IS NOT NULL)" % (
                    OID_fieldname, "', '".join([str(OID) for OID in links_OID]))
            else:
                where_clause = "%s IN ('%s')" % (OID_fieldname, "', '".join([str(OID) for OID in links_OID]))

        MU_database = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath).replace("\mu_Geometry", "")

        msm_Node = os.path.join(MU_database, "msm_Node")

        nodes_invert_level = {row[0]: row[1] for row in arcpy.da.SearchCursor(msm_Node, ["MUID", "InvertLevel"])}

        network = networker.NetworkLinks(MU_database, map_only="link", filter_sql_query = where_clause)

        if is_sqlite:
            with sqlite3.connect(
                    MU_database) as connection:
                update_cursor = connection.cursor()
                arcpy.AddMessage(where_clause)
                with arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath.replace("!delete!",""), ["MUID", "uplevel", "dwlevel"],
                                           where_clause=where_clause) as cursor:
                    for row in cursor:
                        uplevel = row[1]
                        dwlevel = row[2]
                        if not if_equal_to_invert_level or row[1] == nodes_invert_level[network.links[row[0]].fromnode]:
                            update_cursor.execute("UPDATE msm_Link SET uplevel = NULL WHERE MUID = '%s'" % (row[0]))
                            arcpy.AddMessage("Set pipe %s uplevel to NULL" % (row[0]))
                        if not if_equal_to_invert_level or row[2] == nodes_invert_level[
                            network.links[row[0]].tonode]:
                            update_cursor.execute("UPDATE msm_Link SET dwlevel = NULL WHERE MUID = '%s'" % (row[0]))
                            arcpy.AddMessage("Set pipe %s dwlevel to NULL" % (row[0]))
                connection.commit()
        else:
            edit = arcpy.da.Editor(MU_database)
            edit.startEditing(False, True)
            edit.startOperation()

            links_count = np.sum([1 for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).catalogPath, ["MUID"],
                                                                   where_clause = where_clause)])
            arcpy.SetProgressor("step", "Setting UpLevel and DwLevel", 0, links_count, 1)
            with arcpy.da.UpdateCursor(arcpy.Describe(pipe_layer).catalogPath, ["MUID", "UpLevel", "DwLevel"],
                                       where_clause=where_clause) as cursor:
                for row_i, row in enumerate(cursor):
                    arcpy.SetProgressorPosition(row_i)
                    try:
                        uplevel = row[1]
                        dwlevel = row[2]
                        if row[1] == nodes_invert_level[network.links[row[0]].fromnode]:
                            row[1] = None
                            arcpy.AddMessage("Set UpLevel for %s to Null" % (row[0]))
                        if row[2] == nodes_invert_level[network.links[row[0]].tonode]:
                            row[2] = None
                            arcpy.AddMessage("Set DwLevel for %s to Null" % (row[0]))
                        cursor.updateRow(row)
                    except Exception as e:
                        arcpy.AddError(traceback.format_exc())
                        arcpy.AddError(row)
                        raise(e)
            edit.stopOperation()
            edit.stopEditing(True)

        return


class SetDischargeRegulation(object):
    def __init__(self):
        self.label = "4) Set Discharge Regulation"
        self.description = "4) Set Discharge Regulation"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        discharge = arcpy.Parameter(
            displayName="Discharge [m3/s]",
            name="discharge",
            datatype="double",
            parameterType="optional",
            direction="Input")

        parameters = [pipe_layer, discharge]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if not parameters[0].value:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                     lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                     and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                     and lyr.visible][0]
            if links:
                parameters[0].value = links
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        pipe_layer = parameters[0].Value
        discharge = parameters[1].Value
        links_MUID = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]

        MU_database = (
            os.path.dirname(os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)) if ".mdb" in arcpy.Describe(
                pipe_layer).catalogPath else
            os.path.dirname(arcpy.Describe(pipe_layer).catalogPath))
        is_sqlite = True if ".sqlite" in MU_database else False

        msm_PasReg = os.path.join(MU_database, "msm_PasReg")
        ms_Tab = os.path.join(MU_database, "ms_Tab")
        ms_TabD = os.path.join(MU_database, "ms_TabD")

        links_MUIDs = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")
        end_node_critical = None

        if len(links_MUIDs) == len(
                [row[0] for row in arcpy.da.SearchCursor(arcpy.Describe(pipe_layer).CatalogPath, ["MUID"])]):
            userquery = pythonaddins.MessageBox("Set Discharge Regulation for %d pipes?" % (len(links_MUIDs)),
                                                "Confirm Assignment", 4)
            if not userquery == "Yes":
                return

        passive_regulations = {row[0]:row[1] for row in arcpy.da.SearchCursor(msm_PasReg, ["LinkID", "FunctionID"], where_clause = "LinkID IN ('%s')" % ("', '".join(links_MUIDs)))}
        missing_links = [MUID for MUID in links_MUID if MUID not in passive_regulations]
        msm_Link_Network = networker.NetworkLinks(MU_database, map_only="link", filter_sql_query = "MUID IN ('%s')" % ("', '".join(links_MUIDs)))
        arcpy.AddMessage("TabID IN ('%s')" % ("', '".join(passive_regulations.values())))

        with arcpy.da.UpdateCursor(ms_TabD, ["TabID", "Value2"], where_clause = "TabID IN ('%s')" % ("', '".join(passive_regulations.values()))) as cursor:
            for row in cursor:
                old_discharge = row[1]
                row[1] = discharge
                cursor.updateRow(row)
                arcpy.AddMessage("Changed %s from %d L/s to %d L/s " % (row[0], old_discharge*1e3 if old_discharge else 0, row[1]*1e3))

        with arcpy.da.InsertCursor(msm_PasReg, ["LinkID", "TypeNo", "FunctionID", "ControlNodeAID"]) as cursor:
            for link in missing_links:
                row = [link, 1, "Reg_%s" % (link), msm_Link_Network.links[link].fromnode]
                cursor.insertRow(row)

        with arcpy.da.InsertCursor(ms_Tab, ["MUID", "TypeNo"]) as cursor:
            for link in missing_links:
                row = ["Reg_%s" % (link), 4]
                cursor.insertRow(row)
                arcpy.AddMessage("Inserted regulation for link %s at %d L/s " % (row[0], discharge*1e3))

        with arcpy.da.InsertCursor(ms_TabD, ["TabID", "Sqn", "Value1", "Value2"]) as cursor:
            for link in missing_links:
                cursor.insertRow(["Reg_%s" % (link), 1, -100, discharge])
                cursor.insertRow(["Reg_%s" % (link), 2, 100, discharge])

        return
        
class IncreaseBasinSize(object):
    def __init__(self):
        self.label = "Increase volume of Basin"
        self.description = "Increase volume of Basin"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        basin_layer = arcpy.Parameter(
            displayName="Basin feature layer",
            name="basin_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
    
        depth = arcpy.Parameter(
            displayName="Depth to increase volume at",
            name="depth",
            datatype="double",
            parameterType="Required",
            direction="Input")

        size_at_depth = arcpy.Parameter(
            displayName="Volume",
            name="size_at_depth",
            datatype="double",
            parameterType="Required",
            direction="Input")

        parameters = [basin_layer, depth, size_at_depth]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        basin_layer = parameters[0].Value
        depth = parameters[1].Value
        size_at_depth = parameters[2].Value
        
        
        return


class AnalyzeCatchmentArea(object):
    def __init__(self):
        self.label = "c) Analyze Catchment Area to Link"
        self.description = "c) Analyze Catchment Area to Link"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        link_layer = arcpy.Parameter(
            displayName="Link Layer",
            name="link_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        sql_query_catchments = arcpy.Parameter(
            displayName="Set as Definition Query for Catchments",
            name="sql_query_catchments",
            datatype="GPString",
            category="Additional Settings",
            parameterType="optional",
            direction="Input")

        parameters = [link_layer, sql_query_catchments]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        link_layer = parameters[0].Value
        sql_query_catchments = parameters[1].Value

        MU_database = os.path.dirname(arcpy.Describe(link_layer).catalogPath).replace("\mu_Geometry", "")
        graph = mikegraph.Graph(MU_database, ignore_regulations = True)

        graph.map_network()
        graph._read_catchments(where_clause=sql_query_catchments)

        selected_pipes = [row[0] for row in arcpy.da.SearchCursor(link_layer, ["muid"])]

        total_area = {}
        impervious_area = {}
        reduced_area = {}

        target_manholes = [graph.network.links[link].fromnode for link in selected_pipes]
        for manhole in target_manholes:
            upstream_nodes = graph.find_upstream_nodes(manhole)
            upstream_catchments = graph.find_connected_catchments(upstream_nodes[0])
            total_area[manhole] = np.sum([catchment.area for catchment in upstream_catchments if catchment.area])/1e4
            impervious_area[manhole] = np.sum([catchment.impervious_area for catchment in upstream_catchments if catchment.area])/1e4
            reduced_area[manhole] = np.sum([catchment.reduced_area for catchment in upstream_catchments if catchment.area])/1e4

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)

        MIKE_gdb = os.path.join(MIKE_folder, os.path.splitext(os.path.basename(MU_database))[0])
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
                MIKE_gdb = os.path.join(MIKE_folder,
                                        "%s_%d" % (os.path.splitext(os.path.basename(MU_database))[0], dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        result_layer = getAvailableFilename(arcpy.env.scratchGDB + "\Pipe_Area", parent=MU_database)
        arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, os.path.basename(result_layer), "POLYLINE")
        fields = ["muid", "diameter", "materialid", "slope", "nettypeno", "enabled", "Area", "ImpArea", "RedArea"]

        def addField(shapefile, field_name, datatype):
            i = 1
            while field_name in [f.name for f in arcpy.Describe(shapefile).fields]:
                field_name = "%s_%d" % (field_name, i)
            arcpy.AddField_management(shapefile, field_name, datatype)
            return field_name

        addField(result_layer, "muid", "TEXT")
        addField(result_layer, "diameter", "FLOAT")
        addField(result_layer, "materialid", "TEXT")
        addField(result_layer, "slope", "FLOAT")
        addField(result_layer, "nettypeno", "SHORT")
        addField(result_layer, "enabled", "SHORT")
        addField(result_layer, "Area", "FLOAT")
        addField(result_layer, "ImpArea", "FLOAT")
        addField(result_layer, "RedArea", "FLOAT")

        is_sqlite = True if ".sqlite" in MU_database else False
        with arcpy.da.InsertCursor(result_layer, ["SHAPE@"] + fields) as ins_cursor:
            with arcpy.da.SearchCursor(link_layer,
                                       ["Slope" if is_sqlite else "Slope_C", "Diameter", "MaterialID", "MUID",
                                        "SHAPE@", "NetTypeNo",
                                        "enabled"],
                                       where_clause="MUID IN ('%s')" % ("','".join(selected_pipes))) as cursor:
                for row_i, row in enumerate(cursor):
                    # diameter_old = row[4]
                    slope = row[0] * 1e-2 if row[0] else 0
                    # if writeDischargeInstead:
                    #     diameter = # ins_row = (row[4], row[3], peak_discharge[msm_Link_Network.links[row[3]].fromnode], row[2], row[0], row[5], row[6])
                    #     ins_cursor.insertRow(ins_row)
                    diameter = row[1]
                    fromnode, tonode = graph.network.links[row[3]].fromnode, graph.network.links[row[3]].tonode
                    material = row[2]
                    ins_row = (row[4], row[3], diameter, material, row[0], row[5], row[6],
                               total_area[fromnode],
                               impervious_area[fromnode],
                               reduced_area[fromnode])
                    ins_cursor.insertRow(ins_row)

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE"):
            layer = apmapping.Layer(layer_source)
            if group:
                apmapping.AddLayerToGroup(df, group, layer, "BOTTOM")
            else:
                apmapping.AddLayer(df, layer, "TOP")
            updatelayer = apmapping.ListLayers(mxd, layer.name, df)[0]
            updatelayer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))), workspace_type,
                                          unicode(os.path.basename(source)))

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]

        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Links Dimensioned.lyr", result_layer,
                 workspace_type="FILEGDB_WORKSPACE")
        return


class DrawLongitudinalProfiles(object):
    def __init__(self):
        self.label = "7) Draw Longitudinal Profiles"
        self.description = "2b) Draw Longitudinal Profiles"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        result_files = arcpy.Parameter(
            displayName="RES1D Network Result Files",
            name="result_files",
            datatype="File",
            multiValue=True,
            parameterType="optional",
            direction="Input")
        result_files.filter.list = ["res1d"]

        # # new_names (auto-filled)
        # new_names = arcpy.Parameter(
        #     displayName="New Names",
        #     name="Rename of Result Files?",
        #     datatype="String",
        #     multiValue=True,
        #     parameterType="Optional",
        #     direction="Input"
        # )

        pdf_output = arcpy.Parameter(
            displayName="PDF Output",
            name="pdf_output",
            datatype="File",
            parameterType="Optional",
            direction="Output")
        # pdf_output.filter.list = ["pdf"]

        default_name = "Longitudinal Profiles.pdf"
        default_path = os.path.join(arcpy.env.scratchFolder, default_name)
        pdf_output.value = default_path

        overwrite_or_append = arcpy.Parameter(
            displayName="Append to PDF (Default is overwrite)",
            name="overwrite_or_append",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")

        backup_tempfile = arcpy.Parameter(
            displayName="Backup Temp File Path",
            name="backup_tempfile",
            datatype="String",
            parameterType="Derived",  # hidden and not user editable
            direction="Output")


        parameters = [pipe_layer, result_files, pdf_output, overwrite_or_append, backup_tempfile]
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
            arcpy.AddMessage(old_scale)
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


        import os
        import sqlite3
        import pandas as pd
        import numpy as np
        import networkx as nx
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.patches import Rectangle
        from matplotlib import transforms
        import matplotlib.image as mpimg
        import tempfile

        links_selected = [row[0] for row in arcpy.da.SearchCursor(pipe_layer, ["MUID"])]

        # -----------------------
        # User Inputs
        # -----------------------
        # RES1D_FOLDER = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jens Juuls Vej\MIKE_URBAN\JJV_032\JJV_032_m1d - Result Files\JJV_032_CDS_5_133_240BaseDefault_Network_HD.res1d"
        SQLITE_FILE = os.path.dirname(arcpy.Describe(pipe_layer).catalogPath)
        # SQLITE_FILE = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Hasle Torv\MIKE_URBAN\HAT_063\HAT_063.sqlite"
        OUTPUT_PDF = output_pdf
        # output_pdf_temporary = OUTPUT_PDF
        # if overwrite_or_append:
        #     tmp_pdf
        #     with tempfile.NamedTemporaryFile(suffix='.pdf', delete=True) as tmp_pdf:
        #         tmp_path = tmp_pdf.name
        #     output_pdf_temporary =

        import time

        class TimerLogger:
            def __init__(self):
                import time
                self.time = time
                self.last = self.time.time()

            def log(self, msg):
                now = self.time.time()
                text = "{} (Elapsed: {:.2f} seconds)".format(msg, now - self.last)
                self.last = now
                arcpy.SetProgressor("default", msg)

                return text

        tlog = TimerLogger()

        # -----------------------
        # 1) Load Nodes & Links + Geoms
        # -----------------------
        arcpy.AddMessage(tlog.log("Loading Nodes & Links"))
        node_table = os.path.join(SQLITE_FILE, "msm_Node")
        link_table = os.path.join(SQLITE_FILE, "msm_Link")
        # --- Data container classes ---
        class Link:
            """Holds Link attributes and native geometry"""

            def __init__(self, muid, fromnodeid, tonodeid,
                         length, diameter, uplevel, dwlevel, geometry):
                self.muid = muid
                self.fromnodeid = fromnodeid
                self.tonodeid = tonodeid
                self.length = length
                self.diameter = diameter
                self.uplevel = uplevel
                self.dwlevel = dwlevel
                self.geometry = geometry
                if arcgis_pro:
                    self.shapely_geom = wkb.loads(bytes(geometry.WKB))

        class Node:
            """Holds Node attributes and native geometry"""

            def __init__(self, muid, invertlevel, groundlevel, geometry):
                self.muid = muid
                self.invertlevel = invertlevel if invertlevel != -99 else None
                self.groundlevel = groundlevel if groundlevel and groundlevel != -99 else self.invertlevel
                self.geometry = geometry
                if arcgis_pro:
                    self.shapely_geom = wkb.loads(bytes(geometry.WKB))

        # -----------------------
        # 1) Read Links
        # -----------------------
        # Build WHERE clause for selected links
        link_placeholders = ",".join("'{}'".format(m) for m in links_selected)
        link_where = "muid IN ({})".format(link_placeholders)

        # Determine actual field names for 'from' and 'to' nodes (handles naming variations)
        all_link_fields = [f.name for f in arcpy.ListFields(link_table)]
        fromnode_field = "fromnodeid" if "fromnodeid" in all_link_fields else "fromnode"
        tonode_field = "tonodeid" if "tonodeid" in all_link_fields else "tonode"

        # Fields to retrieve, including the SHAPE token
        link_fields = [
            "muid",
            fromnode_field,
            tonode_field,
            "length",
            "diameter",
            "UpLevel",
            "DwLevel",
            "SHAPE@"
        ]

        links = {}
        # Use arcpy.da.SearchCursor to fetch rows and attach geometries
        with arcpy.da.SearchCursor(link_table, link_fields, link_where) as cursor:
            for muid, frm, to, length, diam, up, dw, shape in cursor:
                link = Link(
                    muid=muid,
                    fromnodeid=frm,
                    tonodeid=to,
                    length=length,
                    diameter=diam,
                    uplevel=up,
                    dwlevel=dw,
                    geometry=shape
                )
                link.length = link.length if link.length else shape.length
                links[muid] = link

        # -----------------------
        # 2) Read Nodes
        # -----------------------
        # Extract unique node IDs from loaded links
        node_ids = {link.fromnodeid for link in links.values()} | {link.tonodeid for link in links.values()}

        # Build WHERE clause for selected nodes
        node_placeholders = ",".join("'{}'".format(nid) for nid in node_ids)
        node_where = "muid IN ({})".format(node_placeholders)

        # Fields to retrieve, including the SHAPE token
        node_fields = [
            "muid",
            "invertlevel",
            "groundlevel",
            "SHAPE@"
        ]

        nodes = {}
        # Use arcpy.da.SearchCursor to fetch node rows and attach geometries
        with arcpy.da.SearchCursor(node_table, node_fields, node_where) as cursor:
            for muid, invertlevel, groundlevel, shape in cursor:
                node = Node(
                    muid=muid,
                    invertlevel=invertlevel,
                    groundlevel=groundlevel if groundlevel else invertlevel,
                    geometry=shape
                )
                nodes[muid] = node

        # Log results
        arcpy.AddMessage("Loaded {} links and {} nodes".format(len(links), len(nodes)))

        # -----------------------
        # 2) Read Scenarios
        # -----------------------
        arcpy.AddMessage(tlog.log("Reading Result Files"))
        scenario_data = {}
        if result_files and arcgis_pro:
            # for f in result_files:
            #     name = os.path.basename(os.path.splitext(f)[0])
            #     res = open_res1d(f, nodes = (list(de_nodes.index)))
            #     gdf = res.network.nodes.to_geopandas(agg='max', include_derived=False)
            #     gdf = gdf.rename(columns={'max_WaterLevel': 'water_level'}).set_index('name')
            #     scenario_data[name] = gdf['water_level'].to_dict()

            class Pipe:
                def __init__(self, muid, start_node, end_node):
                    self.muid = muid
                    self.start_node = start_node
                    self.end_node = end_node
                    self.water_level_start = None
                    self.water_level_end = None

            for f in result_files:
                name = os.path.basename(os.path.splitext(f)[0])
                scenario_data[name] = []
                res1d = open_res1d(f, reaches = links_selected)
                for pipe in links_selected:
                    if pipe in res1d.reaches:
                        queries = [QueryDataReach("WaterLevel", pipe, 0),
                                   QueryDataReach("WaterLevel", pipe, res1d.reaches[pipe].length)]
                        query_result = res1d.read(queries).max()
                        pipe_result = Pipe(pipe, res1d.reaches[pipe].start_node, res1d.reaches[pipe].end_node)

                        pipe_result.water_level_start = query_result.iloc[0]
                        pipe_result.water_level_end = query_result.iloc[1]
                        scenario_data[name].append(pipe_result)

                # res = open_res1d(f, nodes = (list(de_nodes.index)))
                # gdf = res.network.nodes.to_geopandas(agg='max', include_derived=False)
                # gdf = gdf.rename(columns={'max_WaterLevel': 'water_level'}).set_index('name')
                # scenario_data[name] = gdf['water_level'].to_dict()

        arcpy.AddMessage(tlog.log("Graphing"))
        # -----------------------
        # 3) Graph & Paths
        # -----------------------
        G = nx.DiGraph()
        for link in links.values():
            G.add_edge(link.fromnodeid, link.tonodeid)

        sources = [n for n, d in G.in_degree() if d == 0]
        sinks = [n for n, d in G.out_degree() if d == 0]
        paths = []
        for s in sources:
            for t in sinks:
                for p in nx.all_simple_paths(G, s, t):
                    if p not in paths:
                        paths.append(p)

        arcpy.AddMessage(tlog.log("Plotting"))

        if arcgis_pro:
            import geopandas as gpd
            # Creating Geoseries
            import contextily as ctx

            origin_CRS = "EPSG:32632"
            webmercator_crs = "EPSG:3857"

            links_geoseries = gpd.GeoSeries([link.shapely_geom for link in links.values()], crs = origin_CRS).to_crs(webmercator_crs)
        # nodes_geoseries = gpd.GeoSeries([node.shapely_geom for node in nodes.values()], crs = origin_CRS).to_crs(webmercator_crs)

        # -----------------------
        # 4) Plot & Save
        # -----------------------
        plt.rcParams['pdf.fonttype'] = 42

        arcpy.AddMessage(tlog.log("Plotting"))
        arcpy.SetProgressor("step", "Plotting...", 0, len(paths), 1)
        with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as tmp_pdf:
            with PdfPages(tmp_pdf) as pdf:
                for path_i, path in enumerate(paths):
                    arcpy.SetProgressorLabel("Plotting Path %s-%s (%d of %d)" % (path[0], path[-1], path_i + 1, len(paths)))
                    arcpy.SetProgressorPosition(path_i)
                    # Calculate Chainage
                    chain = [0.0]
                    for fromnodeid, tonodeid in zip(path, path[1:]):
                        link = [link for link in links.values() if link.fromnodeid == fromnodeid and link.tonodeid == tonodeid][0]

                        chain.append(chain[-1] + link.length)

                    groundlevels = [nodes[muid].groundlevel for muid in path]
                    invertlevels = [nodes[muid].invertlevel for muid in path]

                    # Set up Figures
                    fw = max(8, 8 + 0.5 * (len(path) - 2))
                    fig = plt.figure(figsize=(fw, 5), dpi = 300)
                    ax = fig.add_axes([0.05, 0.15, 0.65, 0.80])  # left, bottom, width, height

                    # bar width
                    diffs = np.diff(chain) if len(chain) > 1 else [1]
                    w = min(diffs) * 0.2;
                    half = w / 2

                    # Draw pipes
                    for fromnodeid, tonodeid in zip(path, path[1:]):
                        i = path.index(fromnodeid)
                        x0, x1 = chain[i], chain[i + 1]
                        sx, ex = x0 + half, x1 - half
                        link = [link for link in links.values() if link.fromnodeid == fromnodeid and link.tonodeid == tonodeid][0]

                        # get uplevel, fallback to upstream node invert
                        bu = link.uplevel if link.uplevel else nodes[fromnodeid].invertlevel

                        # get dwlevel, fallback to downstream node invert
                        bd = link.dwlevel if link.dwlevel else nodes[tonodeid].invertlevel

                        d = link.diameter*1e3 # mm

                        if sx and ex and bu and bd:
                            ax.plot([sx, ex], [bu, bd], 'k-', lw=1)
                            ax.plot([sx, ex], [bu + d / 1000, bd + d / 1000], 'k-', lw=1)

                        # diameter label at bottom of main axes
                        mid = 0.5 * (sx + ex)
                        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                        ax.text(mid, 0, u'ø{}'.format(int(d)) if not np.isnan(d) else "",
                                transform=trans,
                                ha='center', va='bottom', fontsize=8)

                    # Draw Manholes
                    ax.plot(chain, groundlevels, 'g-', lw=0.8, label='Ground Level')
                    for x, g, inv in zip(chain, groundlevels, invertlevels):
                        if g and inv:
                            rect = Rectangle((x - half, inv), w, g - inv, fill=False, edgecolor='black', lw=1.0)
                            ax.add_patch(rect)

                    cmap = plt.cm.get_cmap('tab10' if arcgis_pro else "Set1")

                    # for idx, (nm, wl) in enumerate(scenario_data.items()):
                    #     vals = [wl.get(n, np.nan) for n in path]
                    #     ax.plot(chain, vals, '--', color=cmap(idx), label=nm)

                    # Draw results
                    for idx, (nm, pipe_result) in enumerate(scenario_data.items()):
                        modified_chainage = [chain[0]] + [c for c in chain[1:-1] for _ in (0, 1)] + [chain[-1]]
                        link_water_level = []
                        for fromnodeid, tonodeid in zip(path, path[1:]):
                            water_level = next(
                                ((pipe.water_level_start, pipe.water_level_end) for pipe in pipe_result
                                 if pipe.start_node == fromnodeid and pipe.end_node == tonodeid),
                                None
                            )
                            if water_level:
                                link_water_level.extend(water_level)
                            else:
                                link_water_level.extend([np.nan, np.nan])

                        ax.plot(modified_chainage, link_water_level, '--', lw = 0.8, color=cmap(idx), label=nm)

                    # Write MUIDs
                    y_max = ax.get_ylim()[1]
                    for x, n in zip(chain, path):
                        ax.text(x, y_max * 0.95, str(n), ha='center', va='bottom', rotation=90, fontsize=8)

                    # ticks, labels, legend
                    ticks = np.arange(0, np.ceil(chain[-1] / 50) * 50 + 1, 50)
                    ax.set_xticks(ticks)
                    ax.set_xticklabels([str(int(t)) for t in ticks], fontsize=8)
                    ax.set_xlabel('Chainage (m)')
                    ax.set_ylabel('Elevation (m)')
                    ax.set_title(u"Profile: {} → {}".format(path[0], path[-1]))
                    ax.grid(True, linestyle='--', lw=0.5)
                    legend = ax.legend(fontsize='small', loc='upper left', bbox_to_anchor=(1.02, 1.0))

                    # Create square inset map: e.g., 0.22 x 0.22 in figure coords
                    inset_size = 0.55
                    map_ax = fig.add_axes([0.6, 0.15, inset_size, inset_size])  # x0, y0, width, height

                    if False:
                        # Plot full network in light grey
                        links_geoseries.plot(ax=map_ax, color='black', linewidth=0.5)
                        # Highlight selected path
                        segs = []
                        for fromnodeid, tonodeid in zip(path, path[1:]):
                            link = [link for link in links.values() if link.fromnodeid == fromnodeid and link.tonodeid == link.tonodeid][0]
                            segs.append(link.shapely_geom)
                        geo_segs_web = gpd.GeoSeries(segs, crs=origin_CRS).to_crs(webmercator_crs)
                        geo_segs_web.plot(ax=map_ax, color='red', linewidth=1)

                        # Get bounds of the red path
                        xmin, ymin, xmax, ymax = geo_segs_web.total_bounds
                        xmid = (xmin + xmax) / 2
                        ymid = (ymin + ymax) / 2

                        # Ensure square extent
                        half_extent = max((xmax - xmin), (ymax - ymin)) / 2
                        expansion_factor = 2
                        half_extent *= expansion_factor

                        map_ax.set_xlim(xmid - half_extent, xmid + half_extent)
                        map_ax.set_ylim(ymid - half_extent, ymid + half_extent)

                        # Add basemap
                        ctx.add_basemap(map_ax, zoom=19, source=ctx.providers.OpenStreetMap.Mapnik, interpolation='bilinear')

                        map_ax.set_axis_off()
                        map_ax.set_axis_off()
                    elif arcgis_pro:
                        muids = []
                        for fromnodeid, tonodeid in zip(path, path[1:]):
                            muids.append([link.muid for link in links.values() if
                                    link.fromnodeid == fromnodeid and link.tonodeid == tonodeid][0])

                        where_clause = "MUID IN ('%s')" % ("' , '".join(muids))

                        def getExtents(extents):

                            XMin = min([ext.XMin for ext in extents])
                            XMax = max([ext.XMax for ext in extents])
                            YMin = min([ext.YMin for ext in extents])
                            YMax = max([ext.YMax for ext in extents])

                            return arcpy.Extent(XMin, YMin, XMax, YMax)

                        for layer in map_obj.listLayers():
                            try:
                                layer.setSelectionSet()
                            except:
                                pass
                        arcpy.management.SelectLayerByAttribute(pipe_layer_reference, "NEW_SELECTION", where_clause)

                        view.zoomToAllLayers(True)

                        current_scale = view.camera.scale
                        arcpy.AddMessage(current_scale)
                        view.camera.scale = current_scale * 1.5
                        map_obj.referenceScale = 2000

                        temp_dir = tempfile.gettempdir()
                        temp_png = os.path.join(temp_dir, "arcgis_map_snapshot")

                        # Export active map to PNG

                        view.exportToPNG(temp_png, height = 1000, width = 1000)
                        from PIL import Image

                        # Convert PNG to JPG
                        with Image.open(temp_png + ".png") as img:
                            rgb_img = img.convert('RGB')  # Remove alpha channel if present
                            temp_jpg = temp_png + ".jpg"
                            rgb_img.save(temp_jpg, 'JPEG', quality=65, optimize=True)

                        arcpy.AddMessage(temp_jpg)
                        img = mpimg.imread(temp_jpg)
                        map_ax.imshow(img)
                        map_ax.axis('off')

                    # save & close
                    pdf.savefig(fig)
                    # if path_i>2:
                    #     break
                    # plt.close(fig)

                # if overwrite or append is true (meaning it should append), it merges the temporary pdf with the output pdf
            if overwrite_or_append:
                from pypdf import PdfMerger
                merger = PdfMerger()
                merger.append(backup_tempfile)
                merger.append(tmp_pdf)
                merger.write(OUTPUT_PDF)
            else:
                import shutil
                shutil.copy(tmp_pdf.name, OUTPUT_PDF)

        print("Saved profiles to: {}".format(OUTPUT_PDF))

        os.startfile(os.path.dirname(OUTPUT_PDF))  # Opens folder in Explorer
        arcpy.AddMessage(selection_set)
        pipe_layer_reference.setSelectionSet(list(selection_set))
        map_obj.referenceScale = old_scale
        arcpy.AddMessage(old_scale)

        return


class DrawClash(object):
    def __init__(self):
        self.label       = "8) Draw Profile of Pipe Clash"
        self.description = "8) Draw Profile of Pipe Clash"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        pipe_layer = arcpy.Parameter(
            displayName="Pipe feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input",
            multiValue = True)

        parameters = [pipe_layer]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        if not parameters[0].value:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                    and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)] and ("sqlite" in arcpy.Describe(lyr).catalogPath or "mdb" in arcpy.Describe(lyr).catalogPath)
                     and lyr.visible]
            if links:
                parameters[0].value = ";".join(links)
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        import networker
        from matplotlib.patches import Rectangle
        import matplotlib.pyplot as plt

        pipelayers = parameters[0].ValueAsText.split(";")

        pipe_catalogue = pd.DataFrame({
            "Material": ["Concrete"] * 13,
            "Intern_dia": [300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000],
            "Thick_side": [54, 72, 62, 71, 82, 94, 105, 116, 138, 165, 188, 188, 188],
            "Thick_top": [54, 72, 116, 110, 128, 146, 164, 183, 219, 256, 290, 290, 290],
            "Thick_bot": [54, 72, 141, 169, 197, 225, 253, 281, 338, 394, 450, 450, 450],
            "Total_wid": [408, 544, 624, 742, 864, 988, 1110, 1232, 1476, 1730, 1976, 1976 + 200, 1976 + 400],
            "Total_hei": [408, 544, 757, 879, 1025, 1171, 1317, 1464, 1757, 2050, 2340, 2340 + 200, 2340 + 400]
        })

        def get_thick_bot(material, diameter):
            """
            Returns the Thick_bot for the given material and diameter.
            - If material is not in pipe_catalogue, return None.
            - If diameter matches an Intern_dia exactly, return that Thick_bot.
            - If not, pick the first Intern_dia > diameter and return its Thick_bot.
            - If no Intern_dia is larger (i.e., diameter > max in catalogue), return None.
            """
            # 1. Filter rows by material
            df_mat = pipe_catalogue[pipe_catalogue["Material"].apply(lambda m: m.lower() in material.lower())]

            if df_mat.empty:
                return None

            # 2. Sort by internal diameter
            df_sorted = df_mat.sort_values("Intern_dia")

            # 3. Find the first row where Intern_dia >= requested diameter
            larger = df_sorted[df_sorted["Intern_dia"] >= diameter]
            if not larger.empty:
                return larger.iloc[0]["Thick_bot"]
            else:
                return None

        # ─────────────────────────────────────────────────────────────────────────────
        # STEP 1: Read all pipes (geometry + attributes) into a Python list
        # ─────────────────────────────────────────────────────────────────────────────

        # We'll store each pipe as a dict with keys: geometry, fromnode, tonode,
        # diameter, material, uplevel, dwlevel, source_layer
        pipes = []

        networks = {}
        for layer in pipelayers:
            filter_sql_query = "MUID IN ('%s')" % ("', '".join([row[0] for row in arcpy.da.SearchCursor(layer, ["MUID"])]))
            MU_database = os.path.dirname(arcpy.Describe(layer).catalogPath).replace("\mu_Geometry", "")
            MU_database = MU_database.replace("!delete!", "")
            networks[layer] = networker.NetworkLinks(MU_database, filter_sql_query = filter_sql_query)
            # Make sure the layer actually exists:
            if not arcpy.Exists(layer):
                arcpy.AddWarning("Layer not found: {}".format(layer))
                continue

            # Define fields to pull from each pipe feature
            fields = ["SHAPE@", "MUID", "Diameter", "MaterialID"]
            with arcpy.da.SearchCursor(layer, fields) as cursor:
                for row in cursor:
                    geom, muid, dia, matid = row
                    # Skip null geometries or missing fields
                    if geom is None:
                        continue
                    thickness = get_thick_bot(matid, dia * 1000)
                    pipes.append({
                        "geometry": geom,
                        "muid": muid,
                        "diameter": dia,
                        "fromnode": networks[layer].links[muid].fromnode,
                        "tonode": networks[layer].links[muid].tonode,
                        "material": matid,
                        "uplevel": networks[layer].links[muid].uplevel,
                        "dwlevel": networks[layer].links[muid].dwlevel,
                        "source_layer": layer,
                        "thickness": thickness / 1e3 if thickness else 0
                    })

        pipe_sets_checked = []
        for pipe1 in pipes:
            for pipe2 in pipes:
                skip = False
                for pipe_set in pipe_sets_checked:
                    if pipe1["muid"] in pipe_set and pipe2["muid"] in pipe_set:
                        skip = True
                        break

                if not skip:
                    pipe_sets_checked.append([pipe1["muid"], pipe2["muid"]])

                    if not pipe1["muid"] == pipe2["muid"]:
                        intersection_points = pipe1["geometry"].intersect(pipe2["geometry"], 1)

                        if intersection_points:
                            for intersection_point in intersection_points:
                                plt.figure()
                                plt.title("%s - %s" % (pipe1["muid"], pipe2["muid"]))

                                for pipe in [pipe1, pipe2]:
                                    pipe["clash chainage"] = pipe["geometry"].queryPointAndDistance(intersection_point)[
                                        1]
                                    if pipe["clash chainage"] == 0 or pipe["clash chainage"] == pipe["geometry"].length:
                                        plt.close()
                                        continue
                                    for offset in [0, pipe["diameter"], pipe["diameter"] + pipe["thickness"],
                                                   -pipe["thickness"]]:
                                        plt.plot([0 - pipe["clash chainage"],
                                                  pipe["geometry"].length - pipe["clash chainage"]],
                                                 [pipe["uplevel"] + offset, pipe["dwlevel"] + offset], 'k-', lw=1)

                                # plt.axvspan(-3, 3, color='red', alpha=0.1, zorder=-1)

                                # Vertical line at x=0
                                plt.axvline(x=0, color='black', linewidth=1, alpha=0.5)

                                # Add tick marks every 0.1 along the y-axis
                                ymin, ymax = plt.ylim()
                                for y in np.arange(ymin, ymax, 0.1):
                                    plt.plot([-1, 1], [y, y], color='black', linewidth=0.5, alpha=0.5)
                                # plt.grid(axis='y', which='both', linestyle='--', alpha=0.5)
                                # plt.gca().yaxis.set_major_locator(plt.MultipleLocator(0.2))

                                plt.show()
                                # plt.plot(xs, ys + pipe["diameter"] / 1000, 'k-', lw=1)
        return