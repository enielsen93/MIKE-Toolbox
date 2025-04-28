# -*- coding: utf-8 -*-

# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy
import os
import traceback
import networker
import xml.etree.ElementTree as ET
import xml.dom.minidom
import numpy as np
import warnings

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

class Catchment:
    muid = None
    area = None
    persons = None

    imperviousness = None
    reduction_factor = None
    concentration_time = None
    connection = None
    nettypeno = None

    nodeID = None

    use_local_parameters = None

    def __init__(self, MUID):
        self.MUID = MUID

    @property
    def impervious_area(self):
        return self.area*self.imperviousness/1e2

    @property
    def reduced_area(self):
        return self.area*self.imperviousness/1e2*self.reduction_factor

class HParA:
    reduction_factor = None
    concentration_time = None

    def __init__(self, MUID):
        self.MUID = MUID

from arcpy import env
class Toolbox(object):
    def __init__(self):
        self.label =  "Export to CAD"
        self.alias  = "Export to CAD"

        # List of tool classes associated with this toolbox
        self.tools = [DisplayMikeUrbanAsCAD, ExportToCAD, ExportToDUModelBuilder, ExportToDDS]

class DisplayMikeUrbanAsCAD(object):
    def __init__(self):
        self.label = "Display Model as CAD"
        self.description = "Display Model as CAD"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Optional",
            direction="Input")
            
        label_scale = arcpy.Parameter(
            displayName="Label Scale",
            name="label_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Output")
        label_scale.value = 500

        draw_3d = arcpy.Parameter(
            displayName="Draw in 3D",
            name="draw_3d",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        draw_3d.value = True

        separate_diameter = arcpy.Parameter(
            displayName="Separate pipe dimension and material in CAD Levels",
            name="separate_diameter",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        separate_diameter.value = True

        extent_feature = arcpy.Parameter(
            displayName="Crop data to Polygon",
            name="extent_feature",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            category="Additonal Settings",
            direction="Input")
        extent_feature.filter.list = ["Polygon"]

        parameters = [MU_database, label_scale, draw_3d, separate_diameter, extent_feature]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        label_scale = parameters[1].Value
        draw_3d = parameters[2].Value
        separate_diameter = parameters[3].Value
        extent_feature = parameters[4].Value

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        df.label_scale = label_scale

        empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
        empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(MU_database))[0]

        show_depth = True

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None):
            if ".sqlite" in source:
                source_layer = arcpy.mapping.Layer(source)

                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, source_layer, "BOTTOM")
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
                    arcpy.mapping.AddLayer(df, layer, "BOTTOM")
                update_layer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if new_name:
                    update_layer.name = new_name
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                               workspace_type, unicode(os.path.basename(source)))

            if "msm_Node" in source:
                for label_class in update_layer.labelClasses:
                    if show_depth:
                        label_class.expression = label_class.expression.replace("return labelstr",
                                                                                'if [GroundLevel] and [InvertLevel]: labelstr += "\\nD:%1.2f" % ( convertToFloat([GroundLevel]) - convertToFloat([InvertLevel]) )\r\n  return labelstr')

        arcpy.SetProgressorLabel("Creating Workspace")
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

        arcpy.CreateFileGDB_management(MIKE_gdb, "scratch.gdb")
        arcpy.AddMessage(MIKE_gdb)
        arcpy.env.scratchWorkspace = os.path.join(MIKE_gdb, "scratch.gdb")

        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")
        msm_Weir = MU_database + "\mu_Geometry\msm_Weir"
        msm_Pump = MU_database + "\msm_Pump"
        # net_types = {"Regnvand":2, "Spildvand":1, u"Fællesvand":3}

        # for net_type in net_types.keys():
            # addLayer(os.path.dirname(os.path.realpath(__file__)) + ur"\Data\ExportToCAD\",
                     # msm_Node, group=empty_group_layer)

        arcpy.CopyFeatures_management(msm_Node, os.path.join(arcpy.env.scratchWorkspace,
                                                                        os.path.basename(
                                                                            arcpy.Describe(msm_Node).catalogPath)))[0]
        msm_Node = os.path.join(arcpy.env.scratchWorkspace, os.path.join(os.path.basename(
                                                                            arcpy.Describe(msm_Node).catalogPath)))
                                                                            
        arcpy.CopyFeatures_management(msm_Link, os.path.join(arcpy.env.scratchWorkspace,
                                                                        os.path.basename(
                                                                            arcpy.Describe(msm_Link).catalogPath)))[0]
        msm_Link = os.path.join(arcpy.env.scratchWorkspace, os.path.join(os.path.basename(
                                                                            arcpy.Describe(msm_Link).catalogPath)))

        if extent_feature:
            arcpy.SetProgressorLabel("Clipping Features")
            extent_shapes = [row[0] for row in arcpy.da.SearchCursor(extent_feature, ["SHAPE@"])]
            with arcpy.da.UpdateCursor(msm_Node, "SHAPE@") as cursor:
                for row in cursor:
                    touches = False
                    for extent_shape in extent_shapes:
                        if row[0].within(extent_shape):
                            touches = True
                    if not touches:
                        cursor.deleteRow()

            with arcpy.da.UpdateCursor(msm_Link, "SHAPE@") as cursor:
                for row in cursor:
                    touches = False
                    for extent_shape in extent_shapes:
                        if row[0].within(extent_shape):
                            touches = True
                    if not touches:
                        cursor.deleteRow()


        if separate_diameter:
            arcpy.SetProgressorLabel("Filtering msm_Node and msm_Link by polygon")
            arcpy.management.AddField(msm_Link, "Layer", "TEXT")

            diameters = {0.1024: 0.11, 0.149: 0.16, 0.188: 0.2, 0.233: 0.25, 0.276: 0.3, 0.392: 0.4, 0.493: 0.5,
                         0.588: 0.6,
                         0.781: 0.8, 0.885: 0.9, 0.985: 1, 1.085: 1.1, 1.185: 1.2, 1.285: 1.3, 1.385: 1.4,
                         1.485: 1.5,
                         1.585: 1.6}
            diameters_by_index = [outer_diameter for outer_diameter in diameters.values()]

            def pipe_layer(diameter, materialid, typeno, width, height):
                label = ""
                diameter = float(diameter) if diameter else 0
                material = materialid.lower() if materialid else ""
                if "concrete" in material or "beton" in material:
                    materialAbb = "bt"
                elif "plast" in material:
                    materialAbb = "pl"
                    if diameter in diameters:
                        diameter = diameters[diameter]
                    else:
                        dist = [abs(inner_diameter - diameter) for inner_diameter in diameters]
                        diameter = diameters_by_index[dist.index(min(dist))]
                else:
                    materialAbb = ""

                if typeno == 'Rectangular' or typeno == '3':
                    label += "%dx%d" % (float(width) * 1e3, float(height) * 1e3)
                elif typeno == 'CRS' or typeno == '2':
                    label += "CRS"
                else:
                    label += "%d" % (diameter * 1000)

                label += materialAbb
                return label

            with arcpy.da.UpdateCursor(msm_Link, ["MUID", "Diameter", "MaterialID", "Layer", "TypeNo", "Width", "Height"]) as cursor:
                for row_i, row in enumerate(cursor):
                    arcpy.SetProgressorPosition(row_i)
                    row[3] = pipe_layer(row[1], row[2], row[4], row[5], row[6])
                    cursor.updateRow(row)

        if draw_3d:
            arcpy.SetProgressorLabel("Copying msm_Node")
            msm_Node_z = arcpy.CreateFeatureclass_management(arcpy.env.scratchWorkspace,
                                                             os.path.basename(arcpy.Describe(msm_Node).catalogPath) + "_z",
                                                             "POINT", template=msm_Node, has_z="ENABLED")[0]

            arcpy.management.Append(msm_Node, msm_Node_z)

            arcpy.SetProgressor("step", "Setting Z Coordinate of msm_Node", 0, int(arcpy.GetCount_management(msm_Node).getOutput(0)), 1)
            with arcpy.da.UpdateCursor(msm_Node_z, ["SHAPE@Z", "InvertLevel"],
                                       where_clause="InvertLevel IS NOT NULL") as cursor:
                for row_i, row in enumerate(cursor):
                    arcpy.SetProgressorPosition(row_i)
                    row[0] = row[1]
                    cursor.updateRow(row)

            arcpy.SetProgressorLabel("Copying msm_Link")
            msm_Link_z = arcpy.CreateFeatureclass_management(arcpy.env.scratchWorkspace,
                                                             os.path.basename(arcpy.Describe(msm_Link).catalogPath) + "_z",
                                                             "POLYLINE", template=msm_Link, has_z="ENABLED")[0]

            arcpy.management.Append(msm_Link, msm_Link_z)

            arcpy.SetProgressorLabel("Reading Invert Levels")
            nodes_invert_level = {row[0]: row[1] for row in arcpy.da.SearchCursor(os.path.join(MU_database, "msm_Node"), ["MUID", "InvertLevel"])}

            arcpy.SetProgressorLabel("Networking MIKE Urban Database")
            network = networker.NetworkLinks(MU_database)

            arcpy.SetProgressor("step", "Setting Z Coordinate of msm_Link", 0, int(arcpy.GetCount_management(msm_Link).getOutput(0)), 1)
            with arcpy.da.UpdateCursor(msm_Link_z, ["SHAPE@", "MUID", "UpLevel", "DwLevel", "length"]) as cursor:
                for row_i, row in enumerate(cursor):
                    try:
                        arcpy.SetProgressorPosition(row_i)
                        uplevel = nodes_invert_level[network.links[row[1]].fromnode] if not row[2] else row[2]
                        dwlevel = nodes_invert_level[network.links[row[1]].tonode] if not row[3] else row[3]
                        length = network.links[row[1]].length
                        slope = (uplevel - dwlevel) / length

                        linelist = []
                        for part in row[0]:
                            parts = []
                            for part_i, point in enumerate(part):
                                if part_i == 0:
                                    z = uplevel
                                elif part_i == len(part)-1:
                                    z = dwlevel
                                else:
                                    total_distance = 0
                                    point_geometries = [arcpy.PointGeometry(p) for p in part]
                                    for i in range(1, part_i + 1):
                                        total_distance += point_geometries[i - 1].distanceTo(point_geometries[i])
                                    z = uplevel - total_distance * slope

                                parts.append(arcpy.Point(point.X, point.Y, z))
                            linelist.append(parts)

                        row[0] = arcpy.Polyline(arcpy.Array(linelist), arcpy.Describe(msm_Link_z).spatialReference, True)
                        cursor.updateRow(row)
                    except Exception as e:
                        arcpy.AddError(traceback.format_exc())
                        arcpy.AddError(row)
                        raise(e)
        else:
            msm_Node_z = msm_Node
            msm_Link_z = msm_Link


        arcpy.SetProgressorLabel("Adding Layers")
        for manhole_layer in [u"Wastewater Manhole.lyr", u"Rainwater Manhole.lyr", u"Combined Manhole.lyr"]:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\%s" % manhole_layer,
                     msm_Node_z, group=empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

        for pipe_layer in [u"Wastewater Pipe.lyr", "Rainwater Pipe.lyr", "Combined Pipe.lyr"]:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\%s" % pipe_layer,
                     msm_Link_z, group=empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

        # for pipe_layer in [u"Wastewater Pipe.lyr", "Rainwater Pipe.lyr", "Combined Pipe.lyr"]:
        #     addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\%s" % pipe_layer,
        #              msm_Link, group=empty_group_layer)

        addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\Pump.lyr",
                 msm_Pump, group=empty_group_layer)
        addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\CSO.lyr",
                 msm_Weir, group=empty_group_layer)

        # for link in arcpy.da.SearchCursor(msm_Link, ["MUID", "UpLevel", "DwLevel", "SHAPE@"], where_clause = "UpLevel IS NOT NULL OR DwLevel IS NOT NULL"):
        #     for row in cursor:
        #         if row[1] is not None:


class ExportToCAD(object):
    def __init__(self):
        self.label = "Export Model to CAD"
        self.description = "Export Model to CAD"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Optional",
            direction="Input")

        dgn_file = arcpy.Parameter(
            displayName="Output dwg file",
            name="dgn_file",
            datatype="File",
            parameterType="Required",
            direction="Output")
        dgn_file.filter.list = ["dwg"]
        
        label_scale = arcpy.Parameter(
            displayName="Label Scale",
            name="label_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Output")
        label_scale.value = 500
        
        

        parameters = [MU_database, dgn_file, label_scale]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText if parameters[0].ValueAsText else "CAD export"
        dgn_file = parameters[1].ValueAsText
        label_scale = parameters[2].Value

        msm_Node = os.path.join(MU_database, "msm_Node")

        arcpy.env.overwriteOutput = True

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        df.scale = label_scale

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
                MIKE_gdb = os.path.join(MIKE_folder, "%s_%d" % (os.path.splitext(os.path.basename(MU_database))[0], dir_ext))
        mxd.saveACopy(os.path.join(MIKE_gdb, "CAD.mxd"))
        arcpy.CreateFileGDB_management(MIKE_gdb, "scratch.gdb")
        arcpy.AddMessage(MIKE_gdb)
        arcpy.env.scratchWorkspace = os.path.join(MIKE_gdb, "scratch.gdb")

        arcpy.CreateFeatureclass_management(arcpy.env.scratchWorkspace, "Extent", geometry_type="Polygon")
        if os.path.exists(msm_Node):
            extent = arcpy.Describe(msm_Node).extent
        else:
            extent = df.extent

        with arcpy.da.InsertCursor(os.path.join(arcpy.env.scratchWorkspace, "Extent"), ["SHAPE@"]) as cursor:
            # row = [arcpy.Polygon(arcpy.Array([arcpy.Point(extent.XMin, extent.YMin),
            #                                   arcpy.Point(extent.XMax, extent.YMin),
            #                                   arcpy.Point(extent.XMax, extent.YMax),
            #                                   arcpy.Point(extent.XMin, extent.YMax)]))]
            row = [arcpy.Polygon(arcpy.Array([arcpy.Point(-1e9, -1e9),
                                              arcpy.Point(-1e9, 1e9),
                                              arcpy.Point(1e9, 1e9),
                                              arcpy.Point(1e9, -1e9)]))]
            cursor.insertRow(row)

        inPolygonIndexLayer = os.path.join(arcpy.env.scratchWorkspace, "Extent")

        labels = arcpy.TiledLabelsToAnnotation_cartography(
            os.path.join(MIKE_gdb, "CAD.mxd"), "Layers", inPolygonIndexLayer, arcpy.env.scratchWorkspace,
            "anno", "_", label_scale, generate_unplaced_annotation="NOT_GENERATE_UNPLACED_ANNOTATION")

        arcpy.AddMessage(labels)

        arcpy.env.workspace = arcpy.env.scratchWorkspace
        anno_classes = [os.path.join(arcpy.env.scratchWorkspace, fc) for fc in arcpy.ListFeatureClasses(feature_type = "Annotation")]
        feature_layers = [layer.longName for layer in arcpy.mapping.ListLayers(mxd, "", df) if not layer.isGroupLayer and layer.isFeatureLayer and layer.visible]
        layers = anno_classes + feature_layers
        arcpy.AddMessage(layers)
        arcpy.AddMessage(extent)
        arcpy.ExportCAD_conversion(anno_classes + feature_layers,
                                   Output_Type = "DWG_R2010" if ".dwg" in dgn_file else "DGN_V8", Output_File = dgn_file,
                                    Seed_File = os.path.dirname(os.path.realpath(__file__)) +
                                                (r"\Data\ExportToCAD\Seedfile_2D.dwg" if ".dwg" in dgn_file
                                                else r"\Data\ExportToCAD\Seedfile_2D.dgn"))
        return
        
class ExportToDUModelBuilder(object):
    def __init__(self):
        self.label = "Export MIKE Urban Model to D&U Model Builder"
        self.description = "Export MIKE Urban Model to D&U Model Builder"
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
            
        extent_feature = arcpy.Parameter(
            displayName="Crop data to Polygon",
            name="extent_feature",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            category="Additonal Settings",
            direction="Input")
        extent_feature.filter.list = ["Polygon"]

        parameters = [MU_database, extent_feature]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        extent_feature = parameters[1].Value

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]

        empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
        empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(MU_database))[0]

        show_depth = True

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
                update_layer2 = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer2.definitionQuery = definition_query
                if new_name:
                    update_layer2.name = new_name
                update_layer2.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source)))
                
            # if "msm_Node" in source:
                # for label_class in update_layer2.labelClasses:
                    # if show_depth:
                        # label_class.expression = label_class.expression.replace("return labelstr", 'if [GroundLevel] and [InvertLevel]: labelstr += "\\nD:%1.2f" % ( convertToFloat([GroundLevel]) - convertToFloat([InvertLevel]) )\r\n  return labelstr')

        arcpy.SetProgressorLabel("Creating Workspace")
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

        arcpy.CreateFileGDB_management(MIKE_gdb, "scratch.gdb")
        arcpy.AddMessage(MIKE_gdb)
        arcpy.env.scratchWorkspace = os.path.join(MIKE_gdb, "scratch.gdb")

        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Link = os.path.join(MU_database, "msm_Link")
        msm_Weir = MU_database + "\mu_Geometry\msm_Weir"
        msm_Pump = MU_database + "\msm_Pump"
        # net_types = {"Regnvand":2, "Spildvand":1, u"Fællesvand":3}

        # for net_type in net_types.keys():
            # addLayer(os.path.dirname(os.path.realpath(__file__)) + ur"\Data\ExportToCAD\",
                     # msm_Node, group=empty_group_layer)
        
        arcpy.CopyFeatures_management(msm_Node, os.path.join(arcpy.env.scratchWorkspace,
                                                                        os.path.basename(
                                                                            arcpy.Describe(msm_Node).catalogPath)))[0]
        msm_Node = os.path.join(arcpy.env.scratchWorkspace, os.path.join(os.path.basename(
                                                                            arcpy.Describe(msm_Node).catalogPath)))
                                                                            
        arcpy.CopyFeatures_management(msm_Link, os.path.join(arcpy.env.scratchWorkspace,
                                                                        os.path.basename(
                                                                            arcpy.Describe(msm_Link).catalogPath)))[0]
        msm_Link = os.path.join(arcpy.env.scratchWorkspace, os.path.join(os.path.basename(
                                                                            arcpy.Describe(msm_Link).catalogPath)))
        if extent_feature:
            arcpy.SetProgressorLabel("Clipping Features")
            extent_shapes = [row[0] for row in arcpy.da.SearchCursor(extent_feature, ["SHAPE@"])]

            links = set()
            with arcpy.da.UpdateCursor(msm_Link, ["SHAPE@", "MUID"]) as cursor:
                for row in cursor:
                    touches = False
                    for extent_shape in extent_shapes:
                        if row[0].intersect(extent_shape, 2):
                            touches = True
                            links.add(row[1])
                    if not touches:
                        cursor.deleteRow()

            arcpy.SetProgressorLabel("Networking MIKE Urban Database")
            network = networker.NetworkLinks(MU_database, filter_sql_query = "MUID IN ('%s')" % ("', '".join(links)))
            arcpy.SetProgressorLabel("Clipping Features")

            nodes_in_links = set()
            for link in [network.links[MUID] for MUID in links]:
                nodes_in_links.add(link.fromnode)
                nodes_in_links.add(link.tonode)

            with arcpy.da.UpdateCursor(msm_Node, ["SHAPE@", "MUID"]) as cursor:
                for row in cursor:
                    touches = False
                    if row[1] in nodes_in_links:
                        touches = True
                    else:
                        for extent_shape in extent_shapes:
                            if row[0].within(extent_shape):
                                touches = True
                    if not touches:
                        cursor.deleteRow()
        else:
            arcpy.SetProgressorLabel("Networking MIKE Urban Database")
            network = networker.NetworkLinks(MU_database)

        arcpy.SetProgressorLabel("...")
        arcpy.management.AddField(msm_Node, "Label", "TEXT")
        arcpy.management.AddField(msm_Node, "X_Coordinate", "TEXT")
        arcpy.management.AddField(msm_Node, "Y_Coordinate", "TEXT")

        draw_3d = True
        if draw_3d:
            arcpy.SetProgressorLabel("Copying msm_Node")
            msm_Node_z = arcpy.CreateFeatureclass_management(arcpy.env.scratchWorkspace,
                                                             os.path.basename(arcpy.Describe(msm_Node).catalogPath) + "_z",
                                                             "POINT", template=msm_Node, has_z="ENABLED")[0]

            arcpy.management.Append(msm_Node, msm_Node_z)

            arcpy.SetProgressor("step", "Setting Z Coordinate of msm_Node", 0, int(arcpy.GetCount_management(msm_Node).getOutput(0)), 1)
            with arcpy.da.UpdateCursor(msm_Node_z, ["SHAPE@Z", "InvertLevel"],
                                       where_clause="InvertLevel IS NOT NULL") as cursor:
                for row_i, row in enumerate(cursor):
                    arcpy.SetProgressorPosition(row_i)
                    row[0] = row[1]
                    cursor.updateRow(row)

            arcpy.SetProgressorLabel("Copying msm_Link")
            msm_Link_z = arcpy.CreateFeatureclass_management(arcpy.env.scratchWorkspace,
                                                             os.path.basename(arcpy.Describe(msm_Link).catalogPath) + "_z",
                                                             "POLYLINE", template=msm_Link, has_z="ENABLED")[0]

            arcpy.management.Append(msm_Link, msm_Link_z)

            arcpy.SetProgressorLabel("Reading Invert Levels")
            nodes_invert_level = {row[0]: row[1] for row in arcpy.da.SearchCursor(os.path.join(MU_database, "msm_Node"), ["MUID", "InvertLevel"])}

            arcpy.SetProgressor("step", "Setting Z Coordinate of msm_Link", 0, int(arcpy.GetCount_management(msm_Link).getOutput(0)), 1)
            with arcpy.da.UpdateCursor(msm_Link_z, ["SHAPE@", "MUID", "UpLevel", "DwLevel", "length"]) as cursor:
                for row_i, row in enumerate(cursor):
                    try:
                        arcpy.SetProgressorPosition(row_i)
                        uplevel = nodes_invert_level[network.links[row[1]].fromnode] if not row[2] else row[2]
                        dwlevel = nodes_invert_level[network.links[row[1]].tonode] if not row[3] else row[3]
                        length = network.links[row[1]].length
                        slope = (uplevel - dwlevel) / length

                        linelist = []
                        for part in row[0]:
                            parts = []
                            for part_i, point in enumerate(part):
                                if part_i == 0:
                                    z = uplevel
                                elif part_i == len(part)-1:
                                    z = dwlevel
                                else:
                                    total_distance = 0
                                    point_geometries = [arcpy.PointGeometry(p) for p in part]
                                    for i in range(1, part_i + 1):
                                        total_distance += point_geometries[i - 1].distanceTo(point_geometries[i])
                                    z = uplevel - total_distance * slope

                                parts.append(arcpy.Point(point.X, point.Y, z))
                            linelist.append(parts)

                        row[0] = arcpy.Polyline(arcpy.Array(linelist), arcpy.Describe(msm_Link_z).spatialReference, True)
                        cursor.updateRow(row)
                    except Exception as e:
                        arcpy.AddError(traceback.format_exc())
                        arcpy.AddError(row)
                        raise(e)

            msm_Node = msm_Node_z
            msm_Link = msm_Link_z

        
        with arcpy.da.UpdateCursor(msm_Node, ["Label", "X_Coordinate", "Y_Coordinate", "MUID", "SHAPE@X", "SHAPE@Y"]) as cursor:
            for row in cursor:
                row[0] = row[3]
                row[1] = row[4]
                row[2] = row[5]
                cursor.updateRow(row)
        
        arcpy.management.AddField(msm_Link, "Label", "TEXT")
        arcpy.management.AddField(msm_Link, "Descript", "TEXT")
        arcpy.management.AddField(msm_Link, "Start_Node", "TEXT")
        arcpy.management.AddField(msm_Link, "Stop_Node", "TEXT")
        arcpy.management.AddField(msm_Link, "Invert_Start_m", "TEXT")
        arcpy.management.AddField(msm_Link, "Invert_Stop_m", "TEXT")
        arcpy.management.AddField(msm_Link, "WallThick", "DOUBLE", 8, 4)

        diameters = {0.1024: 0.11, 0.149: 0.16, 0.188: 0.2, 0.233: 0.25, 0.276: 0.3, 0.392: 0.4, 0.493: 0.5,
                     0.588: 0.6,
                     0.781: 0.8, 0.885: 0.9, 0.985: 1, 1.085: 1.1, 1.185: 1.2, 1.285: 1.3, 1.385: 1.4,
                     1.485: 1.5,
                     1.585: 1.6}
        diameters_by_index = [outer_diameter for outer_diameter in diameters.values()]

        def pipe_layer(diameter, materialid, typeno, width, height):
            label = u""
            diameter = float(diameter) if diameter else 0
            material = materialid.lower() if materialid else ""
            if "concrete" in material or "beton" in material:
                materialAbb = "bt"
            elif "plast" in material:
                materialAbb = "pl"
                if diameter in diameters:
                    diameter = diameters[diameter]
                else:
                    dist = [abs(inner_diameter - diameter) for inner_diameter in diameters]
                    diameter = diameters_by_index[dist.index(min(dist))]
            else:
                materialAbb = ""

            if typeno == 'Rectangular' or typeno == '3':
                label += "%dx%d" % (float(width) * 1e3, float(height) * 1e3)
            elif typeno == 'CRS' or typeno == '2':
                label += "CRS"
            else:
                label += "%d" % (diameter * 1000)

            label += materialAbb
            return label

        class Pipe_type:
            def __init__(self, diameter, wallthickness):
                self.diameter = diameter
                self.wallthickness = wallthickness

        pipe_catalogue = {}
        with arcpy.da.SearchCursor(os.path.dirname(os.path.realpath(__file__)) + ur"\Data\ExportToCAD\Pipe Catalogue.dbf", ["Pipe_type", "Diameter", "WallThi"]) as cursor:
            for row in cursor:
                pipe_catalogue[row[0]] = Pipe_type(row[1], row[2])


        with arcpy.da.UpdateCursor(msm_Link, ["MUID", "Diameter", "MaterialID", "Descript", "TypeNo", "Width", "Height", "Label", "WallThick"]) as cursor:
            for row_i, row in enumerate(cursor):
                arcpy.SetProgressorPosition(row_i)
                row[3] = pipe_layer(row[1], row[2], row[4], row[5], row[6])
                row[7] = "l" + row[0]
                if row[3] in pipe_catalogue:
                    row[1] = pipe_catalogue[row[3]].diameter
                    row[8] = pipe_catalogue[row[3]].wallthickness
                else:
                    try:
                        row[1] = row[1]*1000 if row[1] else 0
                    except Exception as e:
                        pass
                cursor.updateRow(row)

        arcpy.SetProgressorLabel("Reading Invert Levels")
        nodes_invert_level = {row[0]: row[1] for row in arcpy.da.SearchCursor(os.path.join(MU_database, "msm_Node"), ["MUID", "InvertLevel"])}

        with arcpy.da.UpdateCursor(msm_Link, ["MUID", "UpLevel", "DwLevel", "length", "Start_Node", "Stop_Node", "Invert_Start_m", "Invert_Stop_m"]) as cursor:
            for row_i, row in enumerate(cursor):
                arcpy.SetProgressorPosition(row_i)
                row[3] = network.links[row[0]].length
                row[4] = network.links[row[0]].fromnode
                row[5] = network.links[row[0]].tonode
                row[6] = nodes_invert_level[network.links[row[0]].fromnode] if not row[1] else row[1]
                row[7] = nodes_invert_level[network.links[row[0]].tonode] if not row[2] else row[2]           
                cursor.updateRow(row)
        
        msm_Node_req_fields = ["Label", "X_Coordinate", "Y_Coordinate", "Diameter", "InvertLevel", "GroundLevel", "NetTypeNo", "Length", "SHAPE", "OBJECTID", "SHAPE_Length"]
        arcpy.DeleteField_management(msm_Node, [field.name for field in arcpy.ListFields(msm_Node) if field.name not in msm_Node_req_fields])
        
        msm_Link_req_fields = ["Label", "Start_Node", "Stop_Node", "Diameter", "Descript", "NetTypeNo", "Invert_Start_m", "Length", "Invert_Stop_m", "WallThick", "SHAPE", "OBJECTID", "SHAPE_Length"]
        arcpy.DeleteField_management(msm_Link, [field.name for field in arcpy.ListFields(msm_Link) if field.name not in msm_Link_req_fields])

        arcpy.SetProgressorLabel("Adding Layers")
        # arcpy.AddMessage(msm_Node)
        # layer = arcpy.mapping.Layer(msm_Node)    
        # arcpy.mapping.AddLayerToGroup(df, empty_group_layer, layer, "BOTTOM")
        # layer = arcpy.mapping.Layer(msm_Link)    
        # arcpy.mapping.AddLayerToGroup(df, empty_group_layer, layer, "BOTTOM")
        for manhole_layer in [u"Wastewater Manhole.lyr", u"Rainwater Manhole.lyr", u"Combined Manhole.lyr"]:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + ur"\Data\ExportToCAD\%s" % manhole_layer,
                     msm_Node, group=empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

        for pipe_layer in [u"Wastewater Pipe.lyr", "Rainwater Pipe.lyr", "Combined Pipe.lyr"]:
            addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\%s" % pipe_layer,
                     msm_Link, group=empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

        # for pipe_layer in [u"Wastewater Pipe.lyr", "Rainwater Pipe.lyr", "Combined Pipe.lyr"]:
        #     addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\%s" % pipe_layer,
        #              msm_Link, group=empty_group_layer)

        addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\Pump.lyr",
                 msm_Pump, group=empty_group_layer)
        addLayer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCAD\CSO.lyr",
                 msm_Weir, group=empty_group_layer)

        # for link in arcpy.da.SearchCursor(msm_Link, ["MUID", "UpLevel", "DwLevel", "SHAPE@"], where_clause = "UpLevel IS NOT NULL OR DwLevel IS NOT NULL"):
        #     for row in cursor:
        #         if row[1] is not None:


class ExportToDDS(object):
    def __init__(self):
        self.label = "Export MIKE Urban Model to DanDasGraf"
        self.description = "Export MIKE Urban Model to DanDasGraf"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        msm_Node = arcpy.Parameter(
            displayName="Manhole Layers:",
            name="msm_Node",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        
        msm_Link = arcpy.Parameter(
            displayName="Link Layers:",
            name="msm_Link",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")

        dandas_knuder = arcpy.Parameter(
            displayName="Export Manholes to DanDasGraf XML-file:",
            name="dandas_knuder",
            datatype="File",
            parameterType="Required",
            direction="Output")
        dandas_knuder.filter.list = ["xml"]

        dandas_ledninger = arcpy.Parameter(
            displayName="Export Pipes to DanDasGraf XML-file:",
            name="dandas_ledninger",
            datatype="File",
            parameterType="Required",
            direction="Output")
        dandas_ledninger.filter.list = ["xml"]

        dandas_deloplande = arcpy.Parameter(
            displayName="Export Catchments to DanDasGraf XML-file:",
            name="dandas_deloplande",
            datatype="File",
            parameterType="Optional",
            direction="Output")

        use_pipe_catalogue = arcpy.Parameter(
            displayName="Use Pipe Catalogue:",
            name="use_pipe_catalogue",
            datatype="DEDbaseTable",
            parameterType="Optional",
            direction="Input")
        use_pipe_catalogue.value = os.path.dirname(os.path.realpath(__file__)) + r"\Data\ExportToCad\Pipe_Catalogue.dbf"

        use_thickness = arcpy.Parameter(
            displayName="For concrete pipes assume wall thickness corresponding to:",
            name="use_thickness",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        use_thickness.filter.list = ["Top", "Bottom", "Side"]

        parameters = [msm_Node, msm_Link, dandas_knuder, dandas_ledninger, dandas_deloplande, use_pipe_catalogue, use_thickness]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        msm_Node = parameters[0].value
        msm_Link = parameters[1].value
        mike_urban_database = os.path.dirname(arcpy.Describe(msm_Node).catalogPath).replace("\mu_Geometry", "")
        if os.path.exists(r"C:\Papirkurv"):
            basefolder = r"C:\Papirkurv"
        else:
            basefolder = os.path.dirname(mike_urban_database
                                         )
        if msm_Node and not parameters[2].value:
            parameters[2].value = os.path.join(basefolder, os.path.basename(mike_urban_database).replace(".mdb","_knuder.xml").replace(".sqlite","_knuder.xml"))

        if msm_Link and not parameters[3].value:
            parameters[3].value = os.path.join(basefolder, os.path.basename(mike_urban_database).replace(".mdb","_ledninger.xml").replace(".sqlite","_ledninger.xml"))


        # if not parameters[0].value and not parameters[1].value:
        #     mxd = arcpy.mapping.MapDocument("CURRENT")
        #
        #     nodes = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Point'
        #             and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
        #     if nodes:
        #         parameters[0].value = nodes
        #
        #     links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
        #             and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
        #     if links:
        #         parameters[1].value = links
                    
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        msm_Node = parameters[0].Value
        msm_Link = parameters[1].Value
        dandas_knuder = parameters[2].ValueAsText
        dandas_ledninger = parameters[3].ValueAsText
        dandas_deloplande = parameters[4].ValueAsText
        use_pipe_catalogue = parameters[5].Value
        use_thickness = parameters[6].Value

        pipe_catalogue = []
        class PipeSize:
            def __init__(self, material, internal_diameter, wall_thickness):
                self.material = material
                self.internal_diameter = internal_diameter
                self.wall_thickness = wall_thickness

        if use_pipe_catalogue:
            thickness_field = {"Top": "Thick_top", "Bottom": "Thick_bot", "Side": "Thick_wid"}[use_thickness]
            with arcpy.da.SearchCursor(use_pipe_catalogue, ["Material", "Intern_dia", thickness_field]) as cursor:
                for row in cursor:
                    pipe_catalogue.append(PipeSize(row[0], row[1], row[2]))

        mike_urban_database = os.path.dirname(arcpy.Describe(msm_Node).catalogPath).replace("\mu_Geometry", "")
        is_mike_plus = True if ".sqlite" in mike_urban_database else False
        network = networker.NetworkLinks(mike_urban_database)
        msm_HParA = os.path.join(mike_urban_database, "msm_HParA")
        ms_Catchment = os.path.join(mike_urban_database, "msm_Catchment") if is_mike_plus else os.path.join(mike_urban_database, "ms_Catchment")
        msm_HModA = os.path.join(mike_urban_database, "msm_HModA")
        msm_CatchCon = os.path.join(mike_urban_database, "msm_CatchCon")

        node_root = ET.Element("KnudeGroup")
        node_root.set("xmlns", "http://www.danva.dk/xml/schemas/dandas/20120102")
        Referencesys = ET.SubElement(node_root, "Referencesys")
        ET.SubElement(Referencesys, "KoordinatsysKode").text = "9"
        ET.SubElement(Referencesys, "KotesysKode").text = "1"

        link_root = ET.Element("LedningGroup")
        link_root.set("xmlns", "http://www.danva.dk/xml/schemas/dandas/20120102")
        Referencesys = ET.SubElement(node_root, "Referencesys")
        ET.SubElement(Referencesys, "KoordinatsysKode").text = "9"
        ET.SubElement(Referencesys, "KotesysKode").text = "1"

        nodes = {}

        class Node:
            def __init__(self, muid, x, y):
                self.muid = muid
                self.x = x
                self.y = y
                self.diameter = 1
                self.type_no = 1
                self.invert_level = 0
                self.ground_level = 0
                self.net_type_no = 2
                self.description = ""
                self.shape = None

        outer_diameters_plastic = np.array([0.110, 0.160, 0.200, 0.250, 0.300, 0.400, 0.500, 0.600, 0.800, 0.900, 1.000, 1.200,
                                   1.400, 1.500, 1.600])
        links = {}

        class Link:
            def __init__(self, muid, shape):
                self.muid = muid
                self.shape = shape
                self.diameter = 0.2
                self.material_id = "Concrete (Normal)"
                self.uplevel = None
                self.dwlevel = None
                self.net_type_no = 2
                self.description = ""

            def outer_diameter(self, pipe_catalogue):
                if self.material_id == "Plastic":# and self.diameter not in outer_diameters_plastic and self.diameter < 1.6:
                    wall_thickness = self.wall_thickness(pipe_catalogue)
                    return self.diameter*1e3 + (wall_thickness if wall_thickness else 0)*2
                    # i = np.where(outer_diameters_plastic - self.diameter >= 0)[0][0]
                    # return outer_diameters_plastic[i]
                else:
                    return self.diameter*1e3

            @property
            def material_code(self):
                if "concrete" in self.material_id.lower() or "beton" in self.material_id.lower():
                    return 1
                elif "plastic" in self.material_id.lower():
                    return 4
                elif "gap" in self.material_id.lower():
                    return 8
                else:
                    return 5

            def wall_thickness(self, pipe_catalogue):
                if "concrete" in self.material_id.lower() or "beton" in self.material_id.lower():
                    material = "Concrete"
                elif "plastic" in self.material_id.lower():
                    material = "Plastic"
                elif "gap" in self.material_id.lower():
                    material = "GAP"
                else:
                    material = "None"

                wall_thickness = [pipe_size.wall_thickness for pipe_size in pipe_catalogue if pipe_size.material.lower() == material.lower() and pipe_size.internal_diameter == self.diameter*1e3]
                if wall_thickness:
                    return wall_thickness[0]
                elif material.lower() in [pipe_size.material.lower() for pipe_size in pipe_catalogue]:
                    try:
                        return np.interp(self.diameter*1e3, [pipe_size.internal_diameter for pipe_size in pipe_catalogue if pipe_size.material.lower() == material.lower()], [pipe_size.wall_thickness for pipe_size in pipe_catalogue if pipe_size.material.lower() == material.lower()])
                    except Exception as e:
                        arcpy.AddWarning(e)
                        return None
                else:
                    return None

        catchments_dict = {}
        if True:
            hParA_dict = {}
            with arcpy.da.SearchCursor(msm_HParA, ["MUID", "RedFactor", "ConcTime"]) as cursor:
                for row in cursor:
                    hParA_dict[row[0]] = HParA(row[0])
                    hParA_dict[row[0]].reduction_factor = row[1]
                    hParA_dict[row[0]].concentration_time = row[2]

            msm_HModA_without_ms_Catchment = []
            if is_mike_plus:
                with arcpy.da.SearchCursor(ms_Catchment,
                                           ['MUID', 'SHAPE@AREA', 'Area', 'Persons', "NetTypeNo", "ModelAImpArea",
                                            "ModelAParAID", "ModelALocalNo", "ModelARFactor",
                                            "ModelAConcTime", "SHAPE@"]) as cursor:
                    for row in cursor:
                        catchments_dict[row[0]] = Catchment(row[0])
                        catchments_dict[row[0]].persons = row[3] if row[3] is not None else 0
                        catchments_dict[row[0]].area = row[2] if row[2] is not None else abs(row[1])
                        catchments_dict[row[0]].nettypeno = row[4]
                        catchments_dict[row[0]].use_local_parameters = not bool(row[7])
                        catchments_dict[row[0]].imperviousness = row[5] * 1e2
                        catchments_dict[row[0]].shape = row[10]

                        try:
                            catchments_dict[row[0]].reduction_factor = (hParA_dict[row[6]].reduction_factor
                                                                             if not row[7] == 0 else row[8])
                            catchments_dict[row[0]].concentration_time = (hParA_dict[row[6]].concentration_time / 60.0
                                                                               if not row[7] == 0 else row[9] / 60.0)

                        except Exception as e:
                            catchments_dict[row[0]].concentration_time = 7
                            catchments_dict[row[0]].reduction_factor = 0
                            warnings.warn("%s not found in msm_HParA" % (row[6]))

                with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
                    for row in cursor:
                        catchments_dict[row[0]].nodeID = row[1]

            else:
                with arcpy.da.SearchCursor(msm_HModA,
                                           ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor", "ConcTime"]) as cursor:
                    for row in cursor:
                        catchments_dict[row[0]] = Catchment(row[0])
                        catchments_dict[row[0]].imperviousness = row[1]

                        try:
                            catchments_dict[row[0]].reduction_factor = (hParA_dict[row[2]].reduction_factor
                                                                             if row[3] == 0 else row[4])
                            catchments_dict[row[0]].concentration_time = (hParA_dict[row[2]].concentration_time
                                                                               if row[3] == 0 else row[5])
                        except Exception as e:
                            catchments_dict[row[0]].concentration_time = 7
                            catchments_dict[row[0]].reduction_factor = 0
                            warnings.warn("%s not found in msm_HParA" % (row[2]))

                with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
                    for row in cursor:
                        if row[0] in catchments_dict:
                            catchments_dict[row[0]].nodeID = row[1]

                with arcpy.da.SearchCursor(ms_Catchment,
                                           ['MUID', 'SHAPE@AREA', 'Area', 'Persons', "NetTypeNo", "SHAPE@"]) as cursor:
                    for row in cursor:
                        if row[0] not in catchments_dict:
                            catchments_dict[row[0]] = Catchment(row[0])
                        catchments_dict[row[0]].muid = row[0]
                        catchments_dict[row[0]].persons = row[3] if row[3] is not None else 0
                        catchments_dict[row[0]].area = row[2] * 1e4 if row[2] is not None else row[1]
                        catchments_dict[row[0]].nettypeno = row[4]
                        catchments_dict[row[0]].shape = row[5]

        with arcpy.da.SearchCursor(msm_Node,
                                   ["MUID", "SHAPE@XY", "Diameter", "TypeNo", "InvertLevel", "GroundLevel", "NetTypeNo",
                                    "Description"]) as cursor:
            for row in cursor:
                if row[1] is not None:
                    nodes[row[0]] = Node(row[0], row[1][0], row[1][1])
                    node = nodes[row[0]]
                    if row[2]: node.diameter = row[2]
                    if row[3]: node.type_no = row[3]
                    if row[4]: node.invert_level = row[4]
                    if row[5]: node.ground_level = row[5]
                    if row[6]: node.net_type_no = row[6]
                    if row[7]: node.description = row[7]

        if msm_Link:
            with arcpy.da.SearchCursor(msm_Link,
                                       ["MUID", "SHAPE@", "Diameter", "MaterialID", "UpLevel", "DwLevel", "NetTypeNo",
                                        "Description"]) as cursor:
                for row in cursor:
                    links[row[0]] = Link(row[0], row[1])
                    link = links[row[0]]
                    if row[2]: link.diameter = row[2]
                    if row[3]: link.material_id = row[3]
                    if row[4]: link.uplevel = row[4]
                    if row[5]: link.dwlevel = row[5]
                    if row[6] and row[6] in [1,2,3]: link.net_type_no = row[6]
                    if row[7]: link.description = row[7]

        for node in nodes.values():
            node_dds = ET.SubElement(node_root, "Knude")
            node_dds.set("Knudenavn", node.muid)
            ET.SubElement(node_dds, "Bundkote").text = "%1.2f" % node.invert_level
            ET.SubElement(node_dds, "DiameterBredde").text = "%d" % (node.diameter * 1e3)
            ET.SubElement(node_dds, "KategoriAfloebKode").text = "1"
            ET.SubElement(node_dds, "KnudeKode").text = "1"
            ET.SubElement(node_dds, "Bemaerkning").text = node.description
            ET.SubElement(node_dds, "Terraenkote").text = "%1.2f" % node.ground_level
            ET.SubElement(node_dds, "TypeAfloebKode").text = "%d" % node.net_type_no
            ET.SubElement(node_dds, "XKoordinat").text = "%1.2f" % node.x
            ET.SubElement(node_dds, "YKoordinat").text = "%1.2f" % node.y
            ET.SubElement(node_dds, "XLabel").text = "%1.2f" % (node.x + 5)
            ET.SubElement(node_dds, "YLabel").text = "%1.2f" % (node.y + 5)

            if catchments_dict and dandas_deloplande:
                catchment_items = ET.SubElement(node_dds, "DeloplandItems")
                arcpy.AddMessage(([catchment.MUID for catchment in catchments_dict.values() if catchment.nodeID == node.muid], node.muid))
                for catchment_no, catchment in enumerate([catchment for catchment in catchments_dict.values() if catchment.nodeID == node.muid and catchment.shape.getPart()]):
                    try:
                        catchment_item = ET.SubElement(catchment_items, "Delopland")
                        catchment_item.set("Deloplandnr", "%d" % (catchment_no + 1))

                        ET.SubElement(catchment_item, "Areal").text = "%1.4f" % (catchment.area / 10000)
                        ET.SubElement(catchment_item, "BefaestelsePct").text = "%1.4f" % catchment.imperviousness
                        ET.SubElement(catchment_item, "TekstjusteringKode").text = "4"
                        ET.SubElement(catchment_item, "Tekstvinkel").text = "0"
                        if not catchment.persons == None:
                            ET.SubElement(catchment_item, "PE").text = "%1.1f" % catchment.persons
                        ET.SubElement(catchment_item, "XLabel").text = "%1.3f" % catchment.shape.centroid.X
                        ET.SubElement(catchment_item, "YLabel").text = "%1.3f" % catchment.shape.centroid.Y
                        ET.SubElement(catchment_item, "TekstFaktor").text = "1.00"

                        catchment_item_coord_items = ET.SubElement(catchment_item, "DeloplandKoordItems")
                        for coordi, coord in enumerate(catchment.shape.getPart()[0]):
                            catchment_item_coord_item = ET.SubElement(catchment_item_coord_items, "DeloplandKoord")
                            catchment_item_coord_item.set("Sortering", "%d" % (coordi + 1))
                            ET.SubElement(catchment_item_coord_item, "Xkoordinat").text = "%1.2f" % coord.X
                            ET.SubElement(catchment_item_coord_item, "Ykoordinat").text = "%1.2f" % coord.Y
                    except Exception as e:
                        arcpy.AddMessage("Error on catchment %s" % (catchment.muid))
                        arcpy.AddError(traceback.format_exc())
                    

        if msm_Link:
            for link in links.values():
                # arcpy.AddMessage(links.values())
                if link.muid in network.links and network.links[link.muid].fromnode in nodes and network.links[link.muid].tonode in nodes:
                    link_dds = ET.SubElement(link_root, "Ledning")
                    link_dds.set("NedstroemKnudenavn", network.links[link.muid].tonode)
                    link_dds.set("OpstroemKnudenavn", network.links[link.muid].fromnode)
                    ET.SubElement(link_dds, "KategoriAfloebKode").text = "1"
                    ET.SubElement(link_dds, "TransportKode").text = "1"
                    ET.SubElement(link_dds, "TypeAfloebKode").text = "%d" % link.net_type_no

                    link_dds_parts = ET.SubElement(link_dds, 'DelLedningItems')
                    link_dds_part = ET.SubElement(link_dds_parts, 'DelLedning')
                    link_dds_part.set("NedstroemKnudenavn", network.links[link.muid].tonode)
                    link_dds_part.set("OpstroemKnudenavn", network.links[link.muid].fromnode)
                    if link.uplevel:
                        ET.SubElement(link_dds_part,
                                  'BundloebskoteOpst').text = "%1.2f" % link.uplevel if link.uplevel else "%1.2f" % nodes[
                                    network.links[link.muid].fromnode].invert_level

                    if link.dwlevel:
                        ET.SubElement(link_dds_part,
                                      'BundloebskoteNedst').text = "%1.2f" % link.dwlevel if link.dwlevel else "%1.2f" % nodes[
                            network.links[link.muid].tonode].invert_level
                    ET.SubElement(link_dds_part, "MaterialeKode").text = "%d" % link.material_code
                    ET.SubElement(link_dds_part, 'Handelsmaal').text = "%d" % (link.outer_diameter(pipe_catalogue))
                    ET.SubElement(link_dds_part, 'DiameterIndv').text = "%d" % (link.diameter * 1e3)
                    wall_thickness = link.wall_thickness(pipe_catalogue)
                    if wall_thickness:
                        ET.SubElement(link_dds_part, 'Godstykkelse').text = "%1.1f" % wall_thickness

                    if len(link.shape.getPart()[0]) > 2:
                        bends = ET.SubElement(link_dds_part, 'KnaekpunktItems')
                        for part_i, part in enumerate(link.shape.getPart()[0][1:-1]):
                            bend = ET.SubElement(bends, 'Knaekpunkt')
                            bend.set('Sortering', '%d' % (part_i + 1))
                            ET.SubElement(bend, 'XKoordinat').text = "%1.2f" % (part.X)
                            ET.SubElement(bend, 'YKoordinat').text = "%1.2f" % (part.Y)

                    link_dds_part_labels = ET.SubElement(link_dds_part, 'LabelDelledningItems')
                    link_dds_part_label = ET.SubElement(link_dds_part_labels, 'LabelDelledning')
                    centroid = link.shape.positionAlongLine(link.shape.length / 2)
                    ET.SubElement(link_dds_part_label, 'PunktPaaLednKode').text = "0"
                    ET.SubElement(link_dds_part_label, 'TekstjusteringKode').text = "4"
                    ET.SubElement(link_dds_part_label, 'XLabel').text = "%1.2f" % centroid.firstPoint.X
                    ET.SubElement(link_dds_part_label, 'YLabel').text = "%1.2f" % centroid.firstPoint.Y

        # Skriv XML-fil
        with open(dandas_knuder, "w+") as f:
            f.write(xml.dom.minidom.parseString(ET.tostring(node_root, encoding="UTF-8")).toprettyxml().encode("utf-8"))

        if msm_Link:
            with open(dandas_ledninger, "w+") as f:
                f.write(xml.dom.minidom.parseString(ET.tostring(link_root, encoding="UTF-8")).toprettyxml().encode("utf-8"))

        