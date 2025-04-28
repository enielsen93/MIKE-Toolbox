# -*- coding: utf-8 -*-
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy

if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument
import numpy as np
import csv
import os
import traceback
import re
import scipy.integrate
from collections import namedtuple
import warnings

def getAvailableFilename(filepath, parent=None):
    parent = "F%s" % (parent) if parent and parent[0].isdigit() else None
    parent = os.path.basename(re.sub(r"\.[^\.\\]+$", "", parent)).replace(".", "_").replace("-", "_").replace(" ",
                                                                                                              "_").replace(
        ",", "_") if parent else None
    filepath = "%s\%s_%s" % (os.path.dirname(filepath), parent, os.path.basename(filepath)) if parent else filepath
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath


from arcpy import env


class Toolbox(object):
    def __init__(self):
        self.label = "Display Mike Urban Model"
        self.alias = "Display Mike Urban Model"

        # List of tool classes associated with this toolbox
        self.tools = [DisplayMikeUrban, CopyMUPTemplate]  # DimensionAnalysis, DisplayPipeElevation


class DisplayMikeUrban(object):
    def __init__(self):
        self.label = "1) Display Mike Urban Model"
        self.description = "1) Display Mike Urban Model"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        MU_database.filter.list = ["mdb", "sqlite"]

        join_catchments = arcpy.Parameter(
            displayName="Join catchments with imperviousness from msm_HModA table",
            name="join_catchments",
            datatype="Boolean",
            parameterType="optional",
            direction="Output")
        join_catchments.value = False

        features_to_display = arcpy.Parameter(
            displayName="Display the following features",
            name="features_to_display",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            category="Additional Settings",
            direction="Input")
        features_to_display.filter.type = "ValueList"
        features_to_display.filter.list = ["Manholes", "Pipes", "Basins", "Weirs", "Orifices", "Pumps", "Network Loads",
                                           "Boundary Water Levels", "Catchment Connections", "Catchments",
                                           "Control Actions", "Sensors", "Annotations"]
        features_to_display.value = ["Manholes", "Pipes", "Weirs", "Orifices", "Pumps", "Network Loads",
                                     "Boundary Water Levels", "Catchment Connections", "Catchments", "Control Actions",
                                     "Sensors"]

        show_loss_par = arcpy.Parameter(
            displayName="Show Loss Parameters on Manhole Text",
            name="show_loss_par",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="optional",
            direction="Output")
        show_loss_par.value = False

        show_depth = arcpy.Parameter(
            displayName="Display Depth",
            name="show_depth",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="optional",
            direction="Output")
        show_depth.value = True

        sql_query = arcpy.Parameter(
            displayName="Set as Definition Query for layers",
            name="sql_query",
            datatype="GPString",
            category="Additional Settings",
            parameterType="optional",
            direction="Output")

        save_layer = arcpy.Parameter(
            displayName="Save it to layer",
            name="save_layer",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="optional",
            direction="Output")

        parameters = [MU_database, join_catchments, features_to_display, show_loss_par, show_depth, sql_query, save_layer]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        if parameters[0].ValueAsText and '"' in parameters[0].ValueAsText:
            parameters[0].Value = parameters[0].ValueAsText.replace('"','')
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        join_catchments = parameters[1].Value
        features_to_display = [feature_name.replace("'", "").replace('"', '') for feature_name in
                               parameters[2].ValueAsText.split(";")]
        show_loss_par = parameters[3].Value
        show_depth = parameters[4].Value
        sql_query = parameters[5].Value
        save_layer = parameters[6].Value

        msm_Node = MU_database + "\msm_Node"
        msm_Link = MU_database + "\mu_Geometry\msm_Link"
        # arcpy.AddMessage(MU_database)
        ms_Catchment = MU_database + "\mu_Geometry\ms_Catchment" if not ".sqlite" in MU_database else MU_database + "\msm_Catchment"
        msm_CatchCon = MU_database + "\mu_Geometry\msm_CatchConLink" if not ".sqlite" in MU_database else MU_database + "\msm_CatchCon"
        msm_Weir = MU_database + "\mu_Geometry\msm_Weir" if not ".sqlite" in MU_database else MU_database + "\msm_Weir"
        msm_Orifice = MU_database + "\msm_Orifice"
        msm_BBoundary = MU_database + "\msm_BBoundary"
        msm_Pump = MU_database + "\msm_Pump"
        msm_BItem = MU_database + "\msm_BItem"
        msm_PasReg = MU_database + "\msm_PasReg"
        ms_TabD = MU_database + "\ms_TabD"
        msm_RTC = MU_database + "\msm_RTCDevice"
        msm_RTCSensor = MU_database + "\msm_RTCSensor"
        mxd = arcpyMapDocument("CURRENT")
        df = mxd.listMaps()[0] if arcgis_pro else arcpymapping.ListDataFrames(mxd)[0]

        is_sqlite_database = True if ".sqlite" in MU_database else False

        import time
        start_time = time.time()

        def printStepAndTime(txt):
            arcpy.AddMessage("%s - %d" % (txt, time.time() - start_time))

        printStepAndTime("Adding Empty Group")
        empty_group_mapped = arcpymapping.LayerFile(os.path.dirname(
            os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr") if arcgis_pro else arcpy.mapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = df.addLayer(empty_group_mapped) if arcgis_pro else arcpymapping.AddLayer(df, empty_group_mapped,
                                                                                               "TOP")
        empty_group_layer = df.listLayers('Empty Group')[0] if arcgis_pro else \
        arcpymapping.ListLayers(mxd, "Empty Group", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(MU_database))[0]

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

        arcpy.env.addOutputsToMap = False
        # addLayer = arcpymapping.Layer(os.path.dirname(os.path.realpath(__file__)) + ("\Data\MOUSE Manholes with LossPar.lyr" if show_loss_par else "\Data\MOUSE Manholes.lyr"))

        templates_folder = os.path.dirname(os.path.realpath(
            __file__)) + "\Data\Templates\MIKE Plus" if arcgis_pro and is_sqlite_database else os.path.dirname(
            os.path.realpath(__file__)) + "\Data"
        templates_extension = ".lyrx" if arcgis_pro and is_sqlite_database else ".lyr"


#         import timeit
#         arcpy.AddMessage(timeit.timeit(lambda: , number=5))
#         timeit.timeit('arcpy.CopyFeatures_management(msm_Node, getAvailableFilename(arcpy.env.scratchGDB + "\spatial_reference", parent=MU_database)).getOutput(0)', number=5)
#         spatial_reference_fc = arcpy.CopyFeatures_management(msm_Node, getAvailableFilename(
#                     arcpy.env.scratchGDB + "\spatial_reference", parent=MU_database)).getOutput(0)
        spatial_reference_fc = arcpy.CopyFeatures_management(msm_Node, getAvailableFilename(arcpy.env.scratchGDB + "\spatial_reference", parent=MU_database)).getOutput(0)
        spatial_reference = arcpy.Describe(spatial_reference_fc).spatialReference

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
                     definition_query=None):
            if ".sqlite" in source:
                source_layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpy.mapping.Layer(source)
                # if not "objectid" in [field.name.lower() for field in arcpy.ListFields(source)]:
                #     import sqlite3
                #     with sqlite3.connect(MU_database) as connection:
                #         update_cursor = connection.cursor()
                #         sql_expression = """
                #         PRAGMA foreign_keys=off;
                #         BEGIN TRANSACTION;
                #
                #         ALTER TABLE %s RENAME TO delete_table;
                #
                #         CREATE TABLE %s
                #         (
                #           column1 datatype [ NULL | NOT NULL ],
                #           column2 datatype [ NULL | NOT NULL ],
                #           ...
                #           CONSTRAINT constraint_name UNIQUE (uc_col1, uc_col2, ... uc_col_n)
                #         );
                #
                #         INSERT INTO table_name SELECT * FROM old_table;
                #
                #         COMMIT;
                #
                #         PRAGMA foreign_keys=on;
                #         """
                #         try:
                #             update_cursor.execute("ALTER TABLE %s ADD COLUMN OBJECTID INTEGER" % os.path.basename(source))
                #             sql_expression = "CREATE INDEX OBJECTID ON %s(OBJECTID)" % os.path.basename(source)
                #             update_cursor.execute(sql_expression)
                #         except Exception as e:
                #             arcpy.AddMessage(source)
                #             raise(e)

                if group:
                    if arcgis_pro:
                        update_layer = df.addLayerToGroup(group, source_layer, "BOTTOM")
                    else:
                        arcpymapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    if arcgis_pro:
                        update_layer = df.addLayer(source_layer, "BOTTOM")
                    else:
                        arcpymapping.AddLayer(df, source_layer, "BOTTOM")

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

                # layer_source_mike_plus = layer_source.replace("MOUSE", "MIKE+") if "MOUSE" in layer_source and os.path.exists(layer_source.replace("MOUSE", "MIKE+")) else None
                # layer_source = layer_source_mike_plus if layer_source_mike_plus else layer_source
                # layer = arcpymapping.Layer(layer_source)
                # update_layer.visible = layer.visible
                # update_layer.labelClasses = layer.labelClasses
                # update_layer.showLabels = layer.showLabels
                # update_layer.name = layer.name
                # update_layer.definitionQuery = definition_query

                try:
                    arcpymapping.UpdateLayer(df, update_layer, layer, symbology_only=True)
                except Exception as e:
                    arcpy.AddWarning(source)
                    pass
            else:
                # arcpy.AddMessage(layer_source)
                layer = arcpymapping.LayerFile(layer_source) if arcgis_pro else arcpymapping.Layer(layer_source)
                if group:
                    if arcgis_pro:
                        df.addLayerToGroup(group, layer, "BOTTOM")
                    else:
                        arcpymapping.AddLayerToGroup(df, group, layer, "BOTTOM")
                else:
                    if arcgis_pro:
                        df.addLayer(layer, "BOTTOM")
                    else:
                        arcpymapping.AddLayer(df, layer, "BOTTOM")
                update_layer = df.listLayers(layer.listLayers()[0].name)[0] if arcgis_pro else \
                arcpymapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name

                if arcgis_pro:
                    df.updateConnectionProperties(update_layer.connectionProperties['connection_info']['database'],
                                                  os.path.dirname(source.replace(r"\mu_Geometry", "")))
                else:
                    update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                                   workspace_type, os.path.basename(source))

            if "msm_Node" in source:
                for label_class in (update_layer.listLabelClasses() if arcgis_pro else update_layer.labelClasses):
                    if show_depth:
                        label_class.expression = label_class.expression.replace("return labelstr",
                                                                                'if [GroundLevel] and [InvertLevel]: labelstr += "\\nD:%1.2f" % ( convertToFloat([GroundLevel]) - convertToFloat([InvertLevel]) )\r\n  return labelstr')

        # empty_group_layer.visible = False
        if "Manholes" in features_to_display:
            printStepAndTime("Adding nodes to map")
            arcpy.SetProgressor("default", "Adding nodes to map")

            if is_sqlite_database:
                layer = addLayer(os.path.join(templates_folder,
                                              "msm_Node" + templates_extension if not show_loss_par else "msm_Node+ with LossPar" + templates_extension),
                                 msm_Node, group=empty_group_layer, definition_query=sql_query)
            else:
                layer = addLayer(os.path.join(templates_folder,
                                              "msm_Node" + templates_extension if not show_loss_par else "msm_Node with LossPar" + templates_extension),
                                 msm_Node, group=empty_group_layer, definition_query=sql_query)

        class Basin:
            def __init__(self, geometryID):
                self.geometryID = geometryID if geometryID else ""
                self.value1 = []
                self.value3 = []
                self.edges = []
                self.permanent_level = None

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
                return np.trapz(surface_areas, elevations)

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

        if "Basins" in features_to_display:
            basins = {}
            if True:
                arcpy.SetProgressor("default", "Calculating volume of basins")
                printStepAndTime("Calculating volume of basins")
                # Import basins
                with arcpy.da.SearchCursor(msm_Node, ["MUID", "GeometryID"], where_clause="TypeNo = 2") as cursor:
                    for row in cursor:
                        basins[row[0]] = Basin(row[1])

                if len(basins) > 0:
                    fromnode_fieldname = "fromnode" if not is_sqlite_database else "fromnodeid"
                    tonode_fieldname = "tonode" if not is_sqlite_database else "tonodeid"
                    outlet_feature_classes = {msm_Link: ["UpLevel", "DwLevel"], msm_Orifice: ["InvertLevel"] * 2,
                                              msm_Weir: ["CrestLevel"] * 2}
                    tonodes = {}

                    for feature_class, edgelevel_fieldname in zip(outlet_feature_classes.keys(),
                                                                  outlet_feature_classes.values()):
                        if fromnode_fieldname in [field.name.lower() for field in arcpy.ListFields(feature_class)]:
                            with arcpy.da.SearchCursor(feature_class, [fromnode_fieldname, edgelevel_fieldname[0],
                                                                       edgelevel_fieldname[1], "MUID",
                                                                       tonode_fieldname],
                                                       where_clause="%s IN ('%s')" % (fromnode_fieldname,
                                                                                      "', '".join(
                                                                                          basins.keys()))) as cursor:
                                for row in cursor:
                                    # arcpy.AddMessage(row)
                                    if row[4]:
                                        if row[0] in tonodes:
                                            tonodes[row[0]].append(row[4])
                                        else:
                                            tonodes[row[0]] = [row[4]]
                                    if row[1]:
                                        basins[row[0]].edges.append(basins[row[0]].Edge(row[3], max(row[1], row[2])))

                    for feature_class, edgelevel_fieldname in zip(outlet_feature_classes.keys(),
                                                                  outlet_feature_classes.values()):
                        if fromnode_fieldname in [field.name.lower() for field in arcpy.ListFields(feature_class)]:
                            with arcpy.da.SearchCursor(feature_class,
                                                       [fromnode_fieldname, edgelevel_fieldname[0],
                                                        edgelevel_fieldname[1], "MUID", tonode_fieldname],
                                                       where_clause="%s IN ('%s')" % (fromnode_fieldname,
                                                                                      "', '".join(
                                                                                          [item for sublist in
                                                                                           tonodes.values() for item in
                                                                                           sublist]))) as cursor:
                                for row in cursor:
                                    basin_MUIDs = [a for a in tonodes if row[0] in tonodes[a]]
                                    for basin_MUID in basin_MUIDs:
                                        basins[basin_MUID].edges.append(
                                            basins[basin_MUID].Edge(row[3], max(row[1], row[2])))

                    with arcpy.da.SearchCursor(os.path.join(MU_database, r"ms_TabD"), ["TabID", "Value1", "Value3"],
                                               where_clause="TabID IN ('%s')" % ("', '".join(
                                                   [basin.geometryID for basin in basins.values()]))) as cursor:
                        for row in cursor:
                            basin = [basin for basin in basins.values() if basin.geometryID == row[0]][0]
                            basin.value1.append(row[1])
                            basin.value3.append(row[2])

                    exportBasins = getAvailableFilename(arcpy.env.scratchGDB + r"\basins", parent=MU_database)
                    arcpy.Select_analysis(msm_Node, exportBasins, where_clause="TypeNo = 2")
                    arcpy.management.AddField(exportBasins, "Volume", "FLOAT")
                    arcpy.management.AddField(exportBasins, "MaxArea", "FLOAT")
                    with arcpy.da.UpdateCursor(exportBasins,
                                               ["MUID", "Volume", "CriticalLevel", "GeometryID", "Description",
                                                "GroundLevel", "InvertLevel", "MaxArea"]) as cursor:
                        for row in cursor:
                            if row[0] in basins:
                                basin = basins[row[0]]
                                basin.invert_level = row[6]
                                if row[2]:
                                    basin.critical_level = row[2]
                                try:
                                    row[1] = basin.max_volume
                                    row[7] = basin.max_area

                                    description = ""
#                                     arcpy.AddMessage(description)
                                    elevation_discrepancy = basin.invert_level - basin.value1[0]
                                    for elev in [e for e in sorted(basin.value1) if e < basin.critical_level]:
                                        description += "%1.2f: %d m3\n" % (
                                            elev+elevation_discrepancy, basin.get_volume(
                                                elev+elevation_discrepancy))

#                                     arcpy.AddMessage(description)

                                    # Deprecated
                                    if False and basin.edges:
                                        for edge in basin.edges_sort:
                                            description += "%s (%1.2f): %d m3\n" % (
                                            edge.name, edge.uplevel, basin.get_volume(
                                                edge.uplevel)) if edge.uplevel and edge.uplevel < row[
                                                5] and edge.uplevel > row[
                                                                      6] else ""

#                                     if len(description) > 255 - 30:
#                                         description = ""
#                                         if basin.edges:
#                                             for edge in basin.edges_sort:
#                                                 description += "%s. (%1.2f): %d m3\n" % (
#                                                     edge.name[0:5], edge.uplevel, basin.get_volume(
#                                                         edge.uplevel)) if edge.uplevel and edge.uplevel < row[
#                                                     5] and edge.uplevel > row[
#                                                                               6] else ""

#                                     if len(description) > 255 - 30:
#                                         description = ""
#                                         if basin.edges:
#                                             for edge in basin.edges_sort:
#                                                 description += "%1.2f: %d m3\n" % (edge.uplevel, basin.get_volume(
#                                                     edge.uplevel)) if edge.uplevel and edge.uplevel < row[
#                                                     5] and edge.uplevel > row[
#                                                                           6] else ""

                                    if basin.critical_level and basin.critical_level < row[5]:
                                        description += "Maks. (%1.2f): %d m3\n" % (
                                            basin.critical_level, basin.get_volume(basin.critical_level))
                                    else:
                                        description += "Maks. (%1.2f): %d m3\n" % (
                                        row[5], min(basin.max_volume, basin.get_volume(row[5])))

                                    row[4] = description[:220]
                                    # cursor.updateRow(row)
                                except Exception as e:
                                    arcpy.AddWarning("Error: Could not calculate volume of basin %s" % (row[0]))
                                    arcpy.AddWarning(basin.value1)
                                    arcpy.AddWarning(basin.value3)
                                    arcpy.AddWarning(basin.edges)
                                    arcpy.AddWarning(traceback.format_exc())

                                # arcpy.AddMessage(row)
                                # arcpy.AddMessage(len(description))
                                cursor.updateRow(row)

                    printStepAndTime("Adding basins to map")
                    arcpy.SetProgressor("default", "Adding basins to map")
                    addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basins.lyr",
                             exportBasins, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE",
                             definition_query=sql_query)

        arcpy.SetProgressor("default", "Adding links, weirs and pumps to map")

        if is_sqlite_database:
            links_sql_query = sql_query + " AND Enabled = True" if sql_query else "Enabled = True"
        else:
            links_sql_query = sql_query

        if "Pipes" in features_to_display:
            addLayer(os.path.dirname(
                os.path.realpath(__file__)) + "\Data\msm_Link.lyr" if not is_sqlite_database else os.path.dirname(
                os.path.realpath(__file__)) + "\Data\MIKE+ Links.lyr",
                     msm_Link, group=empty_group_layer, definition_query=links_sql_query)

        if "Weirs" in features_to_display:
            # arcpy.AddMessage(weirs)
            # if len([row[0] for row in arcpy.da.SearchCursor(weirs,["MUID"])])>0:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Weir.lyr",
                     msm_Weir, group=empty_group_layer, definition_query=sql_query)

        if "Orifices" in features_to_display:
            # arcpy.AddMessage(weirs)
            # if len([row[0] for row in arcpy.da.SearchCursor(weirs,["MUID"])])>0:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Orifice.lyr",
                     msm_Orifice, group=empty_group_layer, definition_query=sql_query)

        if "Pumps" in features_to_display:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Pumps.lyr",
                     msm_Pump, group=empty_group_layer, definition_query=sql_query)

        if "Network Loads" in features_to_display:
            printStepAndTime("Adding network loads to map")
            # Create Network Load
            arcpy.SetProgressor("default", "Adding network loads to map")
            networkShape = getAvailableFilename(arcpy.env.scratchGDB + r"\NetworkLoads", parent=MU_database)

            class NetworkLoad():
                Geometry = None
                net_type_no = None

                @property
                def title(self):
                    t = os.path.basename(self.TSConnection) + " " + self.TimeseriesName if self.TSConnection else None
                    return t

            network_loads = []
            networkLoadInProject = False
            fields = ["MUID", "ApplyBoundaryNo", "ConnectionTypeNo", "CatchLoadNo",
                      "IndividualConnectionNo" if not is_sqlite_database else "TypeNo", "NodeID", "LinkID", "CatchID"]
            with arcpy.da.SearchCursor(msm_BBoundary, fields,
                                       where_clause="ApplyBoundaryNo = 1 AND GroupNo = 2 AND TypeNo IN (9,10)" if is_sqlite_database else "ApplyBoundaryNo = 1 AND GroupNo = 2 AND IndividualConnectionNo = 1") as cursor:
                for row in cursor:
                    networkLoadInProject = True
                    network_load = NetworkLoad()

                    for i, val in enumerate(row):
                        setattr(network_load, fields[i], val)
                    network_loads.append(network_load)

            if True:
                # arcpy.AddMessage("MUID IN ('%s')" % "', '".join([network_load.MUID for network_load in network_loads]) if is_sqlite_database else "BoundaryID IN ('%s')" % "', '".join([network_load.MUID for network_load in network_loads]))
                fields = ["MUID" if is_sqlite_database else "BoundaryID", "VariationNo", "ConstantValue",
                          "TSConnection", "TimeseriesName"]
                with arcpy.da.SearchCursor(msm_BBoundary if is_sqlite_database else msm_BItem, fields,
                                           where_clause="MUID IN ('%s')" % "', '".join(
                                               [network_load.MUID for network_load in
                                                network_loads]) if is_sqlite_database else "BoundaryID IN ('%s')" % "', '".join(
                                               [network_load.MUID for network_load in network_loads])) as cursor:
                    for row in cursor:
                        network_load = network_loads[
                            [i for i, network_load in enumerate(network_loads) if network_load.MUID == row[0]][0]]
                        for i, val in enumerate(row):
                            setattr(network_load, fields[i], val)

                        if network_load.VariationNo != 1:
                            network_load.ConstantValue = None
                        elif network_load.VariationNo != 3:
                            network_load.TSConnection = None
                            network_load.TimeseriesName = None

                        if network_load.CatchLoadNo:
                            with arcpy.da.SearchCursor(ms_Catchment, ["MUID", "SHAPE@XY"],
                                                       where_clause="MUID = '%s'" % (network_load.CatchID)) as xycursor:
                                for xyrow in xycursor:
                                    network_load.Geometry = xyrow[1]
                        elif (is_sqlite_database and network_load.TypeNo == 9) or (
                                not is_sqlite_database and network_load.IndividualConnectionNo == 1):
                            with arcpy.da.SearchCursor(msm_Node, ["MUID", "SHAPE@XY", "NetTypeNo"],
                                                       where_clause="MUID = '%s'" % (network_load.NodeID)) as xycursor:
                                for xyrow in xycursor:
                                    network_load.Geometry = xyrow[1]
                                    network_load.net_type_no = xyrow[2]
                        elif (is_sqlite_database and network_load.TypeNo == 10) or (
                                not is_sqlite_database and network_load.IndividualConnectionNo == 2):
                            with arcpy.da.SearchCursor(msm_Link, ["MUID", "SHAPE@XY", "NetTypeNo"],
                                                       where_clause="MUID = '%s'" % (network_load.LinkID)) as xycursor:
                                for xyrow in xycursor:
                                    network_load.Geometry = xyrow[1]
                                    network_load.net_type_no = xyrow[2]
                        else:
                            arcpy.AddError("Unknown error")

            if networkLoadInProject:
                arcpy.CreateFeatureclass_management(os.path.dirname(networkShape), os.path.basename(networkShape),
                                                    "POINT", spatial_reference=spatial_reference)
                arcpy.AddField_management(networkShape, "MUID", "TEXT")
                arcpy.AddField_management(networkShape, "NetTypeNo", "SHORT")
                arcpy.AddField_management(networkShape, "Discharge", "DOUBLE")
                arcpy.AddField_management(networkShape, "Title", "STRING")
                arcpy.AddField_management(networkShape, "NodeID", "STRING")
                with arcpy.da.InsertCursor(networkShape,
                                           ["MUID", "SHAPE@XY", "Discharge", "Title", "NodeID", "NetTypeNo"]) as cursor:
                    for network_load in network_loads:
                        if network_load.Geometry:
                            cursor.insertRow([network_load.MUID, network_load.Geometry, network_load.ConstantValue,
                                              network_load.title, network_load.NodeID, network_load.net_type_no])

                addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Network Load.lyr",
                         networkShape, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE",
                         definition_query=sql_query)

            # if not is_sqlite_database:
            printStepAndTime("Adding passive regulations to map")
            # Adding passive regulations to project
            arcpy.SetProgressor("default", "Adding passive regulations to map")

            class Regulation:
                q_max = 0
                linkID = None
                tabID = None
                net_type_no = None
                shape = None

            passive_regulations = {}
            if is_sqlite_database:
                with arcpy.da.SearchCursor(msm_Link, ["MUID", "FunctionID"], where_clause="FLowRegNo = 1") as cursor:
                    for row in cursor:
                        passive_regulations[row[0]] = Regulation()
                        passive_regulations[row[0]].tabID = row[1]
                        passive_regulations[row[0]].linkID = row[0]
            else:
                with arcpy.da.SearchCursor(msm_PasReg, ["LinkID", "FunctionID"], where_clause="TypeNo = 1") as cursor:
                    for row in cursor:
                        passive_regulations[row[0]] = Regulation()
                        passive_regulations[row[0]].tabID = row[1]
                        passive_regulations[row[0]].linkID = row[0]

            regulationsShape = getAvailableFilename(arcpy.env.scratchGDB + r"\passive_regulations",
                                                    parent=MU_database)
            if len(passive_regulations) > 0:
                with arcpy.da.SearchCursor(ms_TabD, ["TabID", "Sqn", "Value2"], where_clause="TabID IN ('%s')" % (
                "', '".join([regulation.tabID for regulation in passive_regulations.values()]))) as cursor:
                    for row in cursor:
                        linkIDs = [regulation.linkID for regulation in passive_regulations.values() if
                                   regulation.tabID == row[0]]
                        for linkID in linkIDs:
                            passive_regulations[linkID].q_max = row[2] if passive_regulations[linkID].q_max < row[
                                2] else passive_regulations[linkID].q_max

            printStepAndTime("Adding pump regulations to map")
            # Adding passive regulations to project
            arcpy.SetProgressor("default", "Adding pump regulations to map")

            pumps = {}
            with arcpy.da.SearchCursor(msm_Pump, ["MUID", "QMaxSetID", "NetTypeNo", "SHAPE@"],
                                       where_clause="CapTypeNo = 1") as cursor:
                for row in cursor:
                    pumps[row[0]] = Regulation()
                    pumps[row[0]].linkID = row[0]
                    pumps[row[0]].tabID = row[1]
                    pumps[row[0]].shape = row[3]
                    pumps[row[0]].net_type_no = row[2]

            with arcpy.da.SearchCursor(ms_TabD, ["TabID", "Sqn", "Value2"], where_clause="TabID IN ('%s')" % (
                    "', '".join([regulation.tabID for regulation in pumps.values()]))) as cursor:
                for row in cursor:
                    regulations = [regulation for regulation in pumps.values() if regulation.tabID == row[0]]
                    for regulation in regulations:
                        regulation.q_max = row[2] if regulation.q_max < row[2] else regulation.q_max

            arcpy.CreateFeatureclass_management(os.path.dirname(regulationsShape), os.path.basename(regulationsShape),
                                                "POLYLINE", spatial_reference=spatial_reference)
            arcpy.AddField_management(regulationsShape, "LinkID", "TEXT")
            arcpy.AddField_management(regulationsShape, "NetTypeNo", "SHORT")
            arcpy.AddField_management(regulationsShape, "FunctionID", "TEXT")
            arcpy.AddField_management(regulationsShape, "QMax", "FLOAT")
            with arcpy.da.InsertCursor(regulationsShape,
                                       ["SHAPE@", "LinkID", "FunctionID", "QMax", "NetTypeNo"]) as regulation_cursor:
                with arcpy.da.SearchCursor(msm_Link, ["SHAPE@", "MUID", "NetTypeNo"], where_clause="MUID IN ('%s')" % (
                        "', '".join(passive_regulations.keys()))) as link_cursor:
                    for row in link_cursor:
                        passive_regulations[row[1]].net_type_no = row[2]
                        regulation_cursor.insertRow(
                            [row[0], row[1], passive_regulations[row[1]].tabID, passive_regulations[row[1]].q_max,
                             passive_regulations[row[1]].net_type_no])

                for regulation in pumps.values():
                    row = [regulation.shape, regulation.linkID, regulation.tabID, regulation.q_max,
                           regulation.net_type_no]
                    regulation_cursor.insertRow(row)

            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Passive Regulations.lyr",
                     regulationsShape, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE",
                     definition_query=sql_query)

        if "Boundary Water Levels" in features_to_display:
            printStepAndTime("Adding outlets to map")

            # Adding outlets to project
            class Outlet:
                boundary_MUID = None
                boundary_item_MUID = None
                boundary_water_level = None
                geometry = None
                TSConnection = None
                TimeseriesName = None
                net_type_no = None

                def __init__(self, nodeID):
                    self.nodeID = nodeID

                @property
                def title(self):
                    t = os.path.basename(self.TSConnection) + " " + self.TimeseriesName if self.TSConnection else None
                    return t

            outlets = {}
            with arcpy.da.SearchCursor(msm_BBoundary, ["MUID", "NodeID"],
                                       where_clause="ApplyBoundaryNo = 1 AND TypeNo = 12") as cursor:
                for row in cursor:
                    outlets[row[1]] = Outlet(row[1])
                    outlets[row[1]].boundary_MUID = row[0]

            if not is_sqlite_database:
                with arcpy.da.SearchCursor(msm_BItem,
                                           ["MUID", "BoundaryID", "ConstantValue", "VariationNo", "TSConnection",
                                            "TimeseriesName"],
                                           where_clause="BoundaryType = 12 AND TypeNo = 1") as cursor:
                    for row in cursor:
                        nodeIDs = [outlet.nodeID for outlet in outlets.values() if outlet.boundary_MUID == row[1]]
                        if nodeIDs:
                            outlets[nodeIDs[0]].boundary_item_MUID = row[0]
                            if row[3] == 1:
                                outlets[nodeIDs[0]].boundary_water_level = row[2]
                            elif row[3] == 3:
                                outlets[nodeIDs[0]].TSConnection = row[4]
                                outlets[nodeIDs[0]].TimeseriesName = row[5]
            else:
                with arcpy.da.SearchCursor(msm_BBoundary, ["NodeID", "ConstantValue"],
                                           where_clause="ApplyBoundaryNo = 1 AND TypeNo = 12 AND ConnectionTypeNo = 3") as cursor:
                    for row in cursor:
                        outlets[row[0]].boundary_water_level = row[1]

            with arcpy.da.SearchCursor(msm_Node, ["SHAPE@", "MUID", "NetTypeNo"],
                                       where_clause="TypeNo = 3 AND MUID IN ('%s')" % (
                                       "', '".join(outlets.keys()))) as cursor:
                for row in cursor:
                    outlets[row[1]].geometry = row[0]
                    outlets[row[1]].net_type_no = row[2]

            if len([outlet for outlet in outlets.values()]) > 0:  # if any outlets with water level exist
                boundariesShape = getAvailableFilename(arcpy.env.scratchGDB + r"\BoundaryWaterLevel",
                                                       parent=MU_database)
                try:
                    arcpy.CreateFeatureclass_management(os.path.dirname(boundariesShape),
                                                        os.path.basename(boundariesShape), "POINT",
                                                        spatial_reference=spatial_reference)
                    arcpy.AddField_management(boundariesShape, "NodeID", "TEXT")
                    arcpy.AddField_management(boundariesShape, "NetTypeNo", "SHORT")
                    arcpy.AddField_management(boundariesShape, "B_MUID", "TEXT")
                    arcpy.AddField_management(boundariesShape, "BI_MUID", "TEXT")
                    arcpy.AddField_management(boundariesShape, "Wat_Lev", "FLOAT")
                    arcpy.AddField_management(boundariesShape, "Title", "TEXT")

                    with arcpy.da.InsertCursor(boundariesShape,
                                               ["SHAPE@", "NodeID", "B_MUID", "BI_MUID", "Wat_Lev", "Title",
                                                "NetTypeNo"]) as cursor:
                        for outlet in outlets.values():
                            cursor.insertRow(
                                [outlet.geometry, outlet.nodeID, outlet.boundary_MUID, outlet.boundary_item_MUID,
                                 outlet.boundary_water_level, outlet.title, outlet.net_type_no])
                    addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Boundary Water Level.lyr",
                             boundariesShape, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE",
                             definition_query=sql_query)
                except Exception as e:
                    arcpy.AddError(traceback.format_exc())

        if "Catchment Connections" in features_to_display:
            printStepAndTime("Adding catchment connections to map")
            arcpy.SetProgressor("default", "Adding catchment connections to map")
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchment Connections.lyr",
                     msm_CatchCon, group=empty_group_layer)

        if "Catchments" in features_to_display:
            printStepAndTime("Adding catchments to map")
            if join_catchments:
                arcpy.SetProgressor("default", "Joining ms_Catchment and msm_HModA and adding catchments to map")
                ms_Catchment = arcpy.CopyFeatures_management(ms_Catchment, getAvailableFilename(
                    arcpy.env.scratchGDB + "\ms_CatchmentImp", parent=MU_database)).getOutput(0)
                # arcpy.JoinField_management(in_data=ms_Catchment, in_field="MUID", join_table=MU_database + r"\msm_HModA", join_field="CatchID", fields="ImpArea")

                arcpy.management.AddField(ms_Catchment, "ImpArea", "FLOAT")
                arcpy.management.AddField(ms_Catchment, "ParAID", "TEXT")
                arcpy.management.AddField(ms_Catchment, "RedFactor", "FLOAT")
                arcpy.management.AddField(ms_Catchment, "ConcTime", "FLOAT")
                arcpy.management.AddField(ms_Catchment, "InitLoss", "FLOAT")
                arcpy.management.AddField(ms_Catchment, "NodeID", "TEXT")
                arcpy.management.AddField(ms_Catchment, "NodeNT", "SHORT")

                class HParA:
                    reduction_factor = None
                    concentration_time = None
                    initial_loss = None

                hParA_dict = {}
                with arcpy.da.SearchCursor(MU_database + r"\msm_HParA",
                                           ["MUID", "RedFactor", "ConcTime", "InitLoss"]) as cursor:
                    for row in cursor:
                        hParA_dict[row[0]] = HParA()
                        hParA_dict[row[0]].reduction_factor = row[1]
                        hParA_dict[row[0]].concentration_time = row[2]
                        hParA_dict[row[0]].initial_loss = row[3]

                catchments_dict = {}

                class Catchment:
                    imperviousness = None
                    local_parameters = None
                    ParAID = None
                    reduction_factor = None
                    concentration_time = None
                    initial_loss = None
                    node_ID = None
                    node_id_net_type_no = None

                if "mdb" in MU_database:
                    cursor = arcpy.da.SearchCursor(MU_database + r"\msm_HModA",
                                                   ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor", "ConcTime",
                                                    "ILoss"])
                else:
                    cursor = arcpy.da.SearchCursor(ms_Catchment,
                                                   ["muid", "modelaimparea", "modelaparaid", "modelalocalno",
                                                    "modelarfactor", "modelaconctime",
                                                    "modelailoss"])

                for row in cursor:
                    # try:
                    catchments_dict[row[0]] = Catchment()
                    catchments_dict[row[0]].local_parameters = row[3] if "mdb" in MU_database else 1 - row[3]
                    catchments_dict[row[0]].imperviousness = row[1] if "mdb" in MU_database else row[1] * 1e2
                    catchments_dict[row[0]].ParAID = row[2] if not catchments_dict[row[0]].local_parameters else ""
                    if catchments_dict[row[0]].local_parameters or row[2] in hParA_dict:
                        catchments_dict[row[0]].reduction_factor = (hParA_dict[row[2]].reduction_factor
                                                                    if not catchments_dict[row[0]].local_parameters else
                                                                    row[4])
                        catchments_dict[row[0]].concentration_time = (hParA_dict[row[2]].concentration_time
                                                                      if not catchments_dict[
                            row[0]].local_parameters else row[5])
                        catchments_dict[row[0]].initial_loss = (hParA_dict[row[2]].initial_loss
                                                                if not catchments_dict[row[0]].local_parameters else
                                                                row[6])
                    else:
                        arcpy.AddWarning("%s not in msm_HParA" % row[2])
                    # arcpy.AddMessage((catchments_dict[row[0]].ParAID, catchments_dict[row[0]], catchments_dict[row[0]].local_parameters))
                # except Exception as e:
                #     catchments_dict[row[0]].concentration_time = 7
                #     catchments_dict[row[0]].reduction_factor = 0
                #     arcpy.AddWarning("%s not found in msm_HParA" % (row[2]))
                #     arcpy.AddWarning(e)

                del cursor

                nodes_net_type_no = {row[0]: row[1] for row in arcpy.da.SearchCursor(msm_Node, ["MUID", "NetTypeNo"])}

                with arcpy.da.SearchCursor(os.path.join(MU_database, "msm_CatchCon"), ["CatchID", "NodeID"]) as cursor:
                    for row in cursor:
                        if row[0] in catchments_dict:
                            catchments_dict[row[0]].node_ID = row[1]

                            if row[1] in nodes_net_type_no:
                                catchments_dict[row[0]].node_id_net_type_no = nodes_net_type_no[row[1]]
                        else:
                            arcpy.AddWarning(
                                "Catchment connection registers connection to nonexisting catchment %s" % row[0])

                with arcpy.da.UpdateCursor(ms_Catchment,
                                           ["MUID", "ImpArea", "ParAID", "RedFactor", "ConcTime", "InitLoss", "NodeID",
                                            "NodeNT"]) as cursor:
                    for row in cursor:
                        try:
                            catchment = catchments_dict[row[0]]
                        except Exception as e:
                            arcpy.AddWarning("Could not find catchment %s in msm_HModA" % (row[0]))
                        else:
                            row[1] = catchment.imperviousness
                            row[2] = catchment.ParAID
                            row[3] = catchment.reduction_factor
                            row[4] = catchment.concentration_time
                            row[5] = catchment.initial_loss
                            row[6] = catchment.node_ID
                            row[7] = catchment.node_id_net_type_no
                            cursor.updateRow(row)

                # arcpy.JoinField_management(in_data=ms_Catchment, in_field="MUID", join_table=MU_database + r"\msm_CatchCon", join_field="CatchID", fields="NodeID")
                # arcpy.management.AddField(ms_Catchment, "RedFactor", field_type = "FLOAT")

                addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments W Imp Area.lyr",
                         ms_Catchment, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE")

            else:
                # arcpy.AddMessage((os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments WO Imp Area.lyr",
                # catchments, empty_group_layer))
                addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments WO Imp Area.lyr",
                         ms_Catchment, group=empty_group_layer)

        if "Control Actions" in features_to_display:
            control_actions = {}

            class ControlAction:
                def __init__(self, applyno, structuretypeno, structure_id):
                    self.applyno = applyno
                    self.structuretypeno = structuretypeno
                    self.structure_id = structure_id
                    self.shape = None

            with arcpy.da.SearchCursor(msm_RTC,
                                       ["muid", "DeviceNo", "pumpid", "gateid", "weirid", "valveid"]) as cursor:
                # ["muid", "applyno", "structuretypeno", "pumpid", "weirid", "valveid", "orifricegateid", "orificeweirid"], where_clause = "active == 1") as cursor:
                for row in cursor:
                    if row[1] == 1:  # Pump
                        control_actions[row[0]] = ControlAction(1, row[1], row[2])
                    elif row[1] == 3:  # Gate
                        control_actions[row[0]] = ControlAction(1, row[1], row[3])
                    elif any(row[2:]):
                        for i, value in enumerate(row[2:]):
                            control_actions[row[0]] = ControlAction(1, row[1], row[i])

            pumps = [control_action.structure_id for control_action in control_actions.values() if
                     control_action.structuretypeno == 1]
            with arcpy.da.SearchCursor(msm_Pump, ["MUID", "SHAPE@"],
                                       where_clause="MUID IN ('%s')" % "', '".join(pumps)) as cursor:
                for row in cursor:
                    for control_action in [control_action for control_action in control_actions.values() if
                                           control_action.structure_id == row[
                                               0] and control_action.structuretypeno == 1]:
                        control_action.shape = row[1]

            orifices = [control_action.structure_id for control_action in control_actions.values() if
                        control_action.structuretypeno == 3]

            with arcpy.da.SearchCursor(msm_Orifice, ["MUID", "SHAPE@"],
                                       where_clause="MUID IN ('%s')" % "', '".join(orifices)) as cursor:
                for row in cursor:
                    for control_action in [control_action for control_action in control_actions.values() if
                                           control_action.structure_id == row[
                                               0] and control_action.structuretypeno == 3]:
                        control_action.shape = row[1]

                        # arcpy.AddMessage(control_action.shape.firstPoint.X)

            control_actions_fc = getAvailableFilename(arcpy.env.scratchGDB + r"\RTC",
                                                      parent=MU_database)
            arcpy.CreateFeatureclass_management(os.path.dirname(control_actions_fc),
                                                os.path.basename(control_actions_fc), "POLYLINE", spatial_reference=spatial_reference)

            arcpy.AddField_management(control_actions_fc, "muid", "TEXT")
            arcpy.AddField_management(control_actions_fc, "rtc_muid", "TEXT")
            arcpy.AddField_management(control_actions_fc, "struct_type_no", "SHORT")
            # arcpy.AddField_management(control_actions_fc, "muid", "TEXT")

            with arcpy.da.InsertCursor(control_actions_fc, ["SHAPE@", "muid", "rtc_muid", "struct_type_no"]) as cursor:
                for muid in control_actions:
                    cursor.insertRow([control_actions[muid].shape, control_actions[muid].structure_id, muid,
                                      control_actions[muid].structuretypeno])

            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\RTC.lyr",
                     control_actions_fc, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE", new_name="RTC")

        if "Sensors" in features_to_display:
            sensors = {}

            class Sensor:
                def __init__(self, typeno, locationno, locationid):
                    self.typeno = typeno
                    self.locationno = locationno
                    self.locationid = locationid
                    self.shape = None

            with arcpy.da.SearchCursor(msm_RTCSensor,
                                       ["MUID", "TypeNo", "LocationNo", "NodeLocationID", "LinkLocationID",
                                        "WeirLocationID", "OrificeLocationID", "ValveLocationID",
                                        "PumpLocationID"]) as cursor:
                for row in cursor:
                    if row[2] == 1:  # Node
                        sensors[row[0]] = Sensor(row[1], row[2], row[3])
                    elif row[2] == 2:  # Link
                        sensors[row[0]] = Sensor(row[1], row[2], row[4])
                    elif row[2] == 7:  # Pump
                        sensors[row[0]] = Sensor(row[1], row[2], row[8])
                    else:
                        for i, value in enumerate(row[3:]):
                            sensors[row[0]] = Sensor(row[1], row[2], row[i])

            nodes = [sensor.locationid for sensor in sensors.values() if
                     sensor.locationno == 1]

            links = [sensor.locationid for sensor in sensors.values() if
                     sensor.locationno == 2]

            pumps = [sensor.locationid for sensor in sensors.values() if
                     sensor.locationno == 7]

            with arcpy.da.SearchCursor(msm_Node, ["MUID", "SHAPE@"],
                                       where_clause="MUID IN ('%s')" % "', '".join(nodes)) as cursor:
                for row in cursor:
                    for sensor in [sensor for sensor in sensors.values() if
                                   sensor.locationid == row[
                                       0] and sensor.locationno == 1]:
                        sensor.shape = row[1]

            with arcpy.da.SearchCursor(msm_Link, ["MUID", "SHAPE@XY"],
                                       where_clause="MUID IN ('%s')" % "', '".join(links)) as cursor:
                for row in cursor:
                    for sensor in [sensor for sensor in sensors.values() if
                                   sensor.locationid == row[
                                       0] and sensor.locationno == 1]:
                        sensor.shape = row[1]

                with arcpy.da.SearchCursor(msm_Pump, ["MUID", "SHAPE@XY"],
                                           where_clause="MUID IN ('%s')" % "', '".join(pumps)) as cursor:
                    for row in cursor:
                        for sensor in [sensor for sensor in sensors.values() if
                                       sensor.locationid == row[
                                           0] and sensor.locationno == 1]:
                            sensor.shape = row[1]

            sensors_fc = getAvailableFilename(arcpy.env.scratchGDB + r"\RTCSensors",
                                              parent=MU_database)
            arcpy.CreateFeatureclass_management(os.path.dirname(sensors_fc),
                                                os.path.basename(sensors_fc), "POINT", spatial_reference=spatial_reference)

            arcpy.AddField_management(sensors_fc, "muid", "TEXT")
            arcpy.AddField_management(sensors_fc, "typeno", "SHORT")
            arcpy.AddField_management(sensors_fc, "LocationNo", "SHORT")
            arcpy.AddField_management(sensors_fc, "locationid", "TEXT")
            # arcpy.AddField_management(control_actions_fc, "muid", "TEXT")

            with arcpy.da.InsertCursor(sensors_fc, ["SHAPE@", "muid", "typeno", "locationno", "locationid"]) as cursor:
                for muid in sensors:
                    # arcpy.AddMessage((sensors[muid].locationno, sensors[muid].shape))
                    cursor.insertRow([sensors[muid].shape, muid, sensors[muid].typeno, sensors[muid].locationno,
                                      sensors[muid].locationid])

#             arcpy.AddMessage(sensors_fc)
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\RTCSensors.lyr",
                     sensors_fc, group=empty_group_layer, workspace_type="FILEGDB_WORKSPACE", new_name="RTC Sensors")

        old_workspace = arcpy.env.workspace
        if "Annotations" in features_to_display and not is_sqlite_database:
            arcpy.env.workspace = os.path.join(MU_database, "mu_Geometry")
            annotation_classes = [os.path.join(MU_database, fc) for fc in
                                  arcpy.ListFeatureClasses(feature_type="Annotation")]
            for annotation_class in annotation_classes:
                try:
                    # arcpymapping.AddLayerToGroup(df, empty_group_layer, arcpymapping.Layer(annotation_class), "BOTTOM")
                    addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Annotation.lyr",
                             annotation_class, group=empty_group_layer,
                             new_name=empty_group_layer.name + " " + os.path.basename(annotation_class))
                except Exception as e:
                    arcpy.AddWarning(e)
            arcpy.env.workspace = old_workspace

        if save_layer:
            layout_layer_path = MU_database.replace(".mdb", ".lyr")
            arcpy.management.SaveToLayerFile(empty_group_layer, layout_layer_path)
        return


class DimensionAnalysis(object):
    def __init__(self):
        self.label = "Display Mike Urban Model with Q Full"
        self.description = "Display Mike Urban Model with Q Full"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        parameters = [MU_database]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        msm_Node = MU_database + "\mu_Geometry\msm_Node"
        msm_Link = MU_database + "\mu_Geometry\msm_Link"
        msm_CatchCon = MU_database + "\mu_Geometry\msm_CatchConLink"
        mxd = arcpyMapDocument("CURRENT")
        df = arcpymapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        addLayer = arcpymapping.Layer(msm_Node)
        arcpymapping.AddLayer(df, addLayer)
        updatelayer = arcpymapping.ListLayers(mxd, "msm_Node", df)[0]
        sourcelayer = arcpymapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Manholes.lyr")
        arcpymapping.UpdateLayer(df, updatelayer, sourcelayer, False)
        updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'ACCESS_WORKSPACE',
                                      unicode(addLayer.datasetName))

        addLayer = arcpymapping.Layer(msm_Link)
        arcpymapping.AddLayer(df, addLayer)
        updatelayer = arcpymapping.ListLayers(mxd, "msm_Link", df)[0]
        sourcelayer = arcpymapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Links W QFull.lyr")
        arcpymapping.UpdateLayer(df, updatelayer, sourcelayer, False)
        updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'ACCESS_WORKSPACE',
                                      unicode(addLayer.datasetName))
        # arcpy.ApplySymbologyFromLayer_management (addLayer, "Template.lyr")

        ms_CatchmentImp = arcpy.CopyFeatures_management(MU_database + r"\mu_Geometry\ms_Catchment",
                                                        getAvailableFilename(
                                                            arcpy.env.scratchGDB + "\ms_CatchmentImp")).getOutput(0)
        ms_CatchmentImpLayer = arcpy.MakeFeatureLayer_management(ms_CatchmentImp, getAvailableFilename(
            arcpy.env.scratchGDB + "\ms_CatchmentImpLayer")).getOutput(0).name
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID",
                                   join_table=parameters[0].ValueAsText + r"\msm_HModA", join_field="CatchID",
                                   fields="ImpArea")
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID",
                                   join_table=parameters[0].ValueAsText + r"\msm_CatchCon", join_field="CatchID",
                                   fields="NodeID")
        addLayer = arcpymapping.Layer(ms_CatchmentImpLayer)
        arcpymapping.AddLayer(df, addLayer)
        updatelayer = arcpymapping.ListLayers(mxd, ms_CatchmentImpLayer, df)[0]
        sourcelayer = arcpymapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments W Imp Area.lyr")
        arcpymapping.UpdateLayer(df, updatelayer, sourcelayer, False)
        updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE',
                                      unicode(addLayer.datasetName))

        # addLayer = arcpymapping.Layer(MU_database + r"\mu_Geometry\ms_Catchment")
        # arcpymapping.AddLayer(df, addLayer)
        # updatelayer = arcpymapping.ListLayers(mxd, "ms_Catchment", df)[0]
        # sourcelayer = arcpymapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments W Imp Area.lyr")
        # arcpymapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        # arcpy.AddJoin_management(arcpymapping.ListLayers(mxd, "Delopland", df)[0], "MUID", MU_database + r"\msm_HModA", "CatchID")
        # updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'ACCESS_WORKSPACE', unicode(addLayer.datasetName))
        # for lblClass in arcpymapping.ListLayers(mxd, "Delopland", df)[0].labelClasses:
        # lblClass.expression = '"%1.2f bef. ha (%1.0f%s)" % (float( [ms_Catchment.SHAPE_Area])/1e4*float([msm_HModA.ImpArea])/1e2,round(float([msm_HModA.ImpArea])/5,0)*5,"%")'

        addLayer = arcpymapping.Layer(msm_CatchCon)
        arcpymapping.AddLayer(df, addLayer)
        updatelayer = arcpymapping.ListLayers(mxd, "msm_CatchConLink", df)[0]
        sourcelayer = arcpymapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchment Connections.lyr")
        arcpymapping.UpdateLayer(df, updatelayer, sourcelayer, False)
        updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'ACCESS_WORKSPACE',
                                      unicode(addLayer.datasetName))

        # arcpy.RefreshTOC()
        # arcpy.RefreshActiveView()
        return


class DisplayPipeElevation(object):

    def __init__(self):
        self.label = "Show Pipes with Erronous Elevation"
        self.description = "Show Pipes with Erronous Elevation"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        parameters = [MU_database]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        msm_Link = MU_database + "\mu_Geometry\msm_Link"
        mxd = arcpyMapDocument("CURRENT")
        df = arcpymapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        groupLayer = arcpymapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Links Elevation Error.lyr")
        groupLayer = arcpymapping.AddLayer(df, groupLayer)
        groupLayer = arcpymapping.ListLayers(mxd, groupLayer, df)[0]
        # arcpy.AddMessage(groupLayer)
        for layer in groupLayer:
            arcpy.AddMessage(layer.datasetName)
            layer.replaceDataSource(unicode(MU_database), 'ACCESS_WORKSPACE', unicode(layer.datasetName))

        # arcpy.RefreshTOC()
        # arcpy.RefreshActiveView()
        return


class CopyMUPTemplate(object):
    def __init__(self):
        self.label = "a) Use MUP Template for MIKE Urban Database with proper labelling"
        self.description = "a) Use MUP Template for MIKE Urban Database with proper labelling"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        MU_database.filter.list = ["mdb"]

        parameters = [MU_database]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        import shutil
        shutil.copy2(os.path.dirname(os.path.realpath(__file__)) + "\Data\Template.mup",
                     MU_database.replace(".mdb", ".mup"))

        with open(MU_database.replace(".mdb", ".mup"), 'r') as fopen:
            txt = fopen.read()

        pattern = r"(?<=<[aA]ctive[dD]atabase>)[^<+]+"

        txt = re.sub(pattern, os.path.basename(MU_database), txt)

        with open(MU_database.replace(".mdb", ".mup"), 'w') as fopen:
            fopen.write(txt)

        return
