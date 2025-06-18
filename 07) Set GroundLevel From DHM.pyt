# Tool for retrieving ground levels for points (manholes)
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com
import os
import arcpy
import numpy as np
import re
from arcpy import env
import requests
if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
    from arcpy._mapping import Layer
    import pythonaddins
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument
import json
import sqlite3

callname = "AUEZHLODEF"
hallpass = "Berlin_1989"


class Toolbox(object):
    def __init__(self):
        self.label = "Get Terrain Elevation of Points"
        self.alias = "Get Terrain Elevation of Points"

        # List of tool classes associated with this toolbox
        self.tools = [GetTerrainElevation, GetTerrainElevationDandas]


class GetTerrainElevation(object):
    def __init__(self):
        self.label = "1) Assign terrain data MIKE Manholes or Points"
        self.description = "1) Assign terrain data to MIKE Manholes or Points"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        point_features = arcpy.Parameter(
            displayName="Input point file",
            name="point_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        terrain_layer = arcpy.Parameter(
            displayName="Terrain data (raster)",
            name="terrain_layer",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")

        fields = arcpy.Parameter(
            displayName="Field to assign terrain data to",
            name="Fields",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        decimals = arcpy.Parameter(
            displayName="Number of decimal places",
            name="decimals",
            datatype="GPLong",
            category="Additional Settings",
            parameterType="Optional",
            direction="Input")
        decimals.value = 2

        user = arcpy.Parameter(
            displayName="Username for datafordeleren",
            name="user",
            datatype="GPString",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        user.value = callname

        passw = arcpy.Parameter(
            displayName="Password for datafordeleren",
            name="passw",
            datatype="GPString",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        passw.value = hallpass

        coordinate_system = arcpy.Parameter(
            displayName="Coordinate System of input point file",
            name="coordinate_system",
            category="Additional Settings",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")
        coordinate_system.value = arcpy.SpatialReference("ETRS 1989 UTM Zone 32N").exportToString()

        parameters = [point_features, terrain_layer, fields, decimals, user, passw, coordinate_system]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        # mxd = arcpy.mapping.MapDocument("CURRENT")
        # nodes = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Point'
        # and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
        # if nodes:
        # parameters[0].value = nodes[0]

        if not parameters[0].value:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            nodes = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Point'
                    and "groundlevel" in [field.name.lower() for field in arcpy.ListFields(lyr)]][0]
            if nodes:
                parameters[0].value = nodes

        if parameters[0].altered and not parameters[2].value:
            parameters[2].filter.list = [f.name.lower() for f in arcpy.Describe(parameters[0].value).fields]
            if "GroundLevel".lower() in parameters[2].filter.list:
                parameters[2].value = "GroundLevel"

        if parameters[0].altered:
            parameters[6].value = (arcpy.Describe(parameters[0].value).spatialReference.exportToString()
                                   if arcpy.Describe(parameters[0].value).spatialReference.PCSName
                                   else parameters[6].ValueAsText)
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        point_layer = parameters[0].ValueAsText
        DHM_file = parameters[1].ValueAsText
        field = parameters[2].ValueAsText
        decimals = int(parameters[3].ValueAsText) if parameters[3].ValueAsText else 2
        user = parameters[4].ValueAsText
        passw = parameters[5].ValueAsText
        coordinate_system = arcpy.SpatialReference(text=parameters[6].Value)
        is_sqlite_database = True if ".sqlite" in arcpy.Describe(point_layer).catalogPath else False

        MU_database = os.path.dirname(arcpy.Describe(point_layer).catalogPath).replace("\mu_Geometry", "")

        OID_field = arcpy.Describe(point_layer).OIDFieldName if not is_sqlite_database else "muid"
        point_layer_OIDs = [row[0] for row in arcpy.da.SearchCursor(point_layer, OID_field)]

        if not is_sqlite_database:
            edit = arcpy.da.Editor(os.path.dirname(os.path.dirname(arcpy.Describe(point_layer).catalogPath)))
            edit.startEditing(False, True)
            edit.startOperation()

        point_layer_shapes = {}
        points_elevation = {}
        with arcpy.da.SearchCursor(point_layer, [OID_field, "SHAPE@XY"]) as cursor:
            for row in cursor:
                point_layer_shapes[row[0]] = row[1]

        def getTerrainElevation(x, y):
            if DHM_file:
                arcpy.AddMessage(str(arcpy.GetCellValue_management(DHM_file, "%1.2f %1.2f" % (
                    x, y), "1").getOutput(0).replace(".", ",")))
                row[1] = str(
                    arcpy.GetCellValue_management(DHM_file, "%1.2f %1.2f" % (x, y), "1").getOutput(0).replace(".", ","))
            else:
                if not coordinate_system.PCSCode == 25832:
                    point_reprojected = arcpy.PointGeometry(
                        arcpy.Point(x, y),
                        coordinate_system).projectAs(arcpy.SpatialReference("ETRS 1989 UTM Zone 32N"))
                    x, y = [point_reprojected.firstPoint.X, point_reprojected.firstPoint.Y]
                url = (
                        r"https://services.datafordeler.dk/DHMTerraen/DHMKoter/1.0.0/GEOREST/HentKoter?format=json&username=%s&password=%s&geop=POINT(%1.2f%s%1.2f)" %
                        (user, passw, x, " ", y))
                arcpy.AddMessage(url)
                try:
                    return json.loads(requests.get(url, verify=False).text)['HentKoterRespons']["data"][0]["kote"]
                except Exception as e:
                    arcpy.AddError(row[0])
                    arcpy.AddError(url)
                    arcpy.AddError(requests.get(url).text)
                    arcpy.AddError(e.message)

        arcpy.SetProgressor("step", "Answer messagebox (might be hidden behind window)", 0, len(point_layer_shapes), 1)
        if not arcgis_pro:
            userquery = pythonaddins.MessageBox("Assign terrain elevation to %d points?" % (len(point_layer_shapes)),
                                                "Confirm Assignment", 4)
        else:
            import tkinter as tk
            from tkinter import messagebox

            def confirm_assignment(num_points):
                root = tk.Tk()
                root.withdraw()  # Hide the main window
                result = messagebox.askyesno("Confirm Assignment", f"Assign terrain elevation to {num_points} points?")
                root.destroy()
                return result

            userquery = confirm_assignment(len(point_layer_shapes))

        if userquery == "Yes" or (isinstance(userquery, bool) and userquery is True):
            arcpy.SetProgressor("step", "Getting terrain elevation of points", 0, len(point_layer_shapes), 1)
            if is_sqlite_database and 'msm_Node'.lower() in arcpy.Describe(point_layer).catalogPath.lower():
                with sqlite3.connect(
                        MU_database.replace("!delete!", "")) as connection:
                    update_cursor = connection.cursor()
                    for muid in point_layer_shapes.keys():
                        x, y = point_layer_shapes[muid][0], point_layer_shapes[muid][1]
                        terrain_elevation = round(getTerrainElevation(x,y), decimals)
                        update_cursor.execute(
                            "UPDATE msm_Node SET %s = %s WHERE MUID = '%s'" % (
                            field, terrain_elevation, muid))
                        arcpy.AddMessage("Set %s to %1.2f" % (muid, terrain_elevation))

            else:
                with arcpy.da.UpdateCursor(arcpy.Describe(point_layer).catalogPath, [OID_field, field],
                                           where_clause="%s IN (%s)" % (OID_field, ", ".join(
                                                   map(str, point_layer_OIDs)))) as pointcursor:  # 'SHAPE@XY'
                    for i, row in enumerate(pointcursor):
                        arcpy.SetProgressorPosition(i)
                        x, y = point_layer_shapes[row[0]][0], point_layer_shapes[row[0]][1]
                        terrain_elevation = round(getTerrainElevation(x, y), decimals)
                        # arcpy.AddMessage(terrain_elevation)
                        row[1] = terrain_elevation
                        try:
                            if terrain_elevation == "NoData":
                                arcpy.AddWarning("Warning: Found NoData on location %s" % ("%1.2f %1.2f" % (x, y)))
                            else:
                                pointcursor.updateRow(row)
                                arcpy.AddMessage("Set GroundLevel to %1.2f" % terrain_elevation)
                        except Exception as e:
                            arcpy.AddError("Error on row %s" % row[1])
                            arcpy.AddError(e.message)

        if not is_sqlite_database:
            edit.stopOperation()
            try:
                edit.stopEditing(True)
            except RuntimeError as e:
                if "GDB_Release" in e.message:
                    pass
                else:
                    raise (e)
        return


class GetTerrainElevationDandas(object):
    def __init__(self):
        self.label = "2) Assign terrain data to Dandas manholes"
        self.description = "2) Assign terrain data to Dandas manholes"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter

        DDS_database = arcpy.Parameter(
            displayName="Dandas database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        terrain_layer = arcpy.Parameter(
            displayName="Terrain data (raster)",
            name="terrain_layer",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")

        point_features = arcpy.Parameter(
            displayName="Shapefile with nodes to assign terrain level to",
            name="point_features",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        point_features.filter.list = ["Point"]

        point_features = arcpy.Parameter(
            displayName="Shapefile with nodes to assign terrain level to",
            name="point_features",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        point_features.filter.list = ["Point"]

        onlywhenmissing = arcpy.Parameter(
            displayName="Assign only terrain elevation to nodes that are missing it",
            name="onlywhenmissing",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")
        onlywhenmissing.value = "false"

        user = arcpy.Parameter(
            displayName="User for datafordeleren",
            name="user",
            datatype="GPString",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        user.value = callname

        passw = arcpy.Parameter(
            displayName="Pass for datafordeleren",
            name="passw",
            datatype="GPString",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        passw.value = hallpass

        coordinate_system = arcpy.Parameter(
            displayName="Coordinate System of input point file:",
            name="coordinate_system",
            category="Additional Settings",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")
        coordinate_system.value = arcpy.SpatialReference("ETRS 1989 UTM Zone 32N").exportToString()

        parameters = [DDS_database, terrain_layer, point_features, onlywhenmissing, user, passw, coordinate_system]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        if parameters[0].altered:
            parameters[6].value = (arcpy.Describe(parameters[0].value).spatialReference.exportToString()
                                   if arcpy.Describe(parameters[0].value).spatialReference.PCSName
                                   else parameters[6].ValueAsText)
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        terrain_layer = parameters[1].ValueAsText
        DDS_database = parameters[0].ValueAsText
        NodeTable = DDS_database + r"\Knude"
        DaekselTable = DDS_database + r"\Daeksel"
        point_features = parameters[2].ValueAsText
        onlywhenmissing = parameters[3].Value
        user = parameters[4].ValueAsText
        passw = parameters[5].ValueAsText

        arcpy.env.workspace = "in_memory/TerrainData"
        Xs = [0]
        Ys = [0]
        if point_features:
            with arcpy.da.SearchCursor(point_features, ["SHAPE@XY"]) as cursor:
                for row in cursor:
                    Xs.append(row[0][0])
                    Ys.append(row[0][1])
        arcpy.SetProgressor("step", "Assigning terrain data to manhole", 0, int(len(Xs)), 1)

        nodes = {}
        with arcpy.da.UpdateCursor(NodeTable,
                                   ["XKoordinat", "YKoordinat", "Terraenkote", "Knudenavn", "ID"]) as pointcursor:
            for pointi, pointrow in enumerate(pointcursor):
                for i in range(len(Xs)):
                    if not point_features or np.sqrt((Xs[i] - pointrow[0]) ** 2 + (Ys[i] - pointrow[1]) ** 2) < 0.1:
                        # arcpy.AddMessage(pointrow[0])
                        if not onlywhenmissing or not pointrow[2]:
                            if DHM_file:
                                pointrow[2] = str(arcpy.GetCellValue_management(terrain_layer, "%1.2f %1.2f" % (
                                pointrow[0], pointrow[1]), "1").getOutput(0).replace(",", "."))
                            else:
                                url = (
                                            r"https://services.datafordeler.dk/DHMTerraen/DHMKoter/1.0.0/GEOREST/HentKoter?format=json&username=%s&password=%s!&geop=POINT(%1.2f%s%1.2f)" %
                                            (user, passw, point_layer_shapes[row[0]][0], "%",
                                             point_layer_shapes[row[0]][1]))
                                try:
                                    row[1] = \
                                    json.loads(requests.get(url, verify=False).text)['HentKoterRespons']["data"][0][
                                        "kote"]
                                except Exception as e:
                                    arcpy.AddError(row[2])
                                    arcpy.AddError(url)
                                    arcpy.AddError(requests.get(url).text)
                                    arcpy.AddError(e.message)
                            pointcursor.updateRow(pointrow)
                            arcpy.AddMessage("Changed manhole %s" % (pointrow[3]))
                            nodes[pointrow[4]] = pointrow[2]
                            # with arcpy.da.UpdateCursor(DaekselTable, ["KnudeID", "Daekselkote"]) as daekselcursor: #, where_clause = "KnudeID = %d" % (pointrow[4])
                            # for daekselrow in daekselcursor:
                            # daekselrow[1] = pointrow[2]
                            # daekselcursor.updateRow(daekselrow)
                arcpy.SetProgressorPosition(pointi)

        # arcpy.AddMessage("(%s)" % (",".join([str(a) for a in list(nodes.keys())])))
        with arcpy.da.UpdateCursor(DaekselTable, ["KnudeID", "Daekselkote"], where_clause="KnudeID IN (%s)" % (",".join(
                [str(a) for a in
                 list(nodes.keys())]))) as daekselcursor:  # , where_clause = "KnudeID = %d" % (pointrow[4])
            for daekselrow in daekselcursor:
                daekselrow[1] = nodes[daekselrow[0]]
                daekselcursor.updateRow(daekselrow)
                # arcpy.AddMessage("Did not change manhole %s" %(pointrow[3]))

        return