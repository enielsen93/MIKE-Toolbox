# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com
import os
import arcpy
import numpy as np
import re
import pythonaddins
import hashlib
import xlwt
import csv

if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument

def getAvailableFilename(filepath):
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath

def splitWhereclause(iterator):
    iterator = np.array(iterator)
    lengths = 3000-(np.cumsum([len(name)+3 for name in iterator])+50)
    length = [i for i, a in enumerate(lengths) if a < 0][0]

    idx = np.arange(0, len(iterator), length)
    if len(iterator) not in idx:
        idx = np.concatenate((idx, [len(iterator)]), axis= 0)

    where_clauses = []
    for i in range(len(idx)-1):
        where_clauses.append("('%s')" % ("','".join([str(name) for name in iterator[idx[i]:idx[i+1]]])))
    return where_clauses

from arcpy import env
class Toolbox(object):
    def __init__(self):
        self.label =  "Catchments"
        self.alias  = "Catchments"

        # List of tool classes associated with this toolbox
        self.tools = [CatchmentProcessing, CatchmentProcessingAlt, CheckCatchments, FAS2Deloplande, TransferCatchments, DuplicateCatchments, GenerateCatchmentConnections, SetImperviousness, CatchmentProcessingScalgo, CatchmentSlopeAnalysis, CatchmentConnectionDescription]

class CatchmentProcessing(object):
    def __init__(self):
        self.label       = "1a) Calculate imperviousness"
        self.description = "1a) Calculate imperviousness of catchments from shapefiles delineating impervious areas"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        imp_layer1 = arcpy.Parameter(
            displayName="Shapefile no. 1 delineating impervious areas",
            name="imp_layer1",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        imp_layer1.filter.list = ["Polygon"]
        
        imp_layer2 = arcpy.Parameter(
            displayName="Shapefile no. 2 delineating impervious areas",
            name="imp_layer2",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        imp_layer2.filter.list = ["Polygon"]

        imp_layer3 = arcpy.Parameter(
            displayName="Shapefile no. 3 delineating impervious areas",
            name="imp_layer3",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        imp_layer3.filter.list = ["Polygon"]

        catchments = arcpy.Parameter(
            displayName="Catchment layer",
            name="catchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchments.filter.list = ["Polygon"]

        parameters = [imp_layer1, imp_layer2, imp_layer3, catchments]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        # imparea = os.path.join(arcpy.env.scratchGDB,"ImpArea")
        #catchintsec = os.path.join(arcpy.env.scratchGDB,"catchintsec")
        #catchspatialjoin = os.path.join(arcpy.env.scratchGDB,"CatchSpatialJoin")
        #ms_catchments = os.path.join(arcpy.env.scratchGDB, "ms_catchments")
        imparea = os.path.join("in_memory","ImpArea")
        catchintsec = os.path.join("in_memory","catchintsec")
        catchspatialjoin = os.path.join("in_memory","CatchSpatialJoin")
        ms_catchments = os.path.join("in_memory", "ms_catchments")
        catchments = parameters[3].ValueAsText
        arcpy.SetProgressor ("default", "Preparing catchments")
        arcpy.CopyFeatures_management(catchments,ms_catchments)
        
        arcpy.SetProgressor ("default", "Answer user promt (may be behind window)")
        userquery = pythonaddins.MessageBox("Calculate and assign imperviousness to %s catchments?" % (arcpy.GetCount_management(ms_catchments)[0]), "Confirm Assignment", 4)
        if userquery == "Yes":
            arcpy.AddMessage ("Repairing catchments layer")
            arcpy.RepairGeometry_management(ms_catchments)
            arcpy.DefineProjection_management(in_dataset=ms_catchments, coor_system="PROJCS['ETRS_1989_UTM_Zone_32N',GEOGCS['GCS_ETRS_1989',DATUM['D_ETRS_1989',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',9.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]")

            arcpy.AddMessage ("Preparing and combining shapefiles delineating impervious areas")
            #arcpy.SetProgressor ("default", "Clipping shapefiles delineating impervious areas to speed up geoprocessing calculation")
            extentPolygon = arcpy.CreateFeatureclass_management("in_memory", "extentPolygon", geometry_type = "POLYGON")
            extentCoordinates = arcpy.Describe(ms_catchments).extent
            array = arcpy.Array([arcpy.Point(extentCoordinates.XMin, extentCoordinates.YMin),
                         arcpy.Point(extentCoordinates.XMax, extentCoordinates.YMin),
                         arcpy.Point(extentCoordinates.XMax, extentCoordinates.YMax),
                         arcpy.Point(extentCoordinates.XMin, extentCoordinates.YMax),])
            polygon = arcpy.Polygon(array)
            with arcpy.da.InsertCursor(extentPolygon, ['SHAPE@']) as cursor:
                cursor.insertRow([polygon])
            
            layersBasename = []
            layers = []
            for i in range(2):
                if parameters[i].value:
                    layersBasename.append(arcpy.Describe(parameters[i].value).baseName)
                    if arcpy.Describe(parameters[i].value).dataElementType == "DEShapeFile":
                        layers.append(arcpy.MakeFeatureLayer_management(parameters[i].value,"in_memory\layer%d.lyr" % i))
                    else:
                        layers.append(parameters[i].value)

            layersClipped = []
            for i,layer in enumerate(layers):
                try:
                    # layersClipped.append(os.path.join("in_memory","layersClipped%d" % i))
                    arcpy.AddMessage ("Clipping layer %s" % layersBasename[i])
                    layersClipped.append(arcpy.SelectLayerByLocation_management(in_layer=layer, overlap_type="INTERSECT", select_features=extentPolygon, search_distance="", selection_type="NEW_SELECTION", invert_spatial_relationship="NOT_INVERT"))
                    # arcpy.Clip_analysis(layer, extentPolygon, os.path.join("in_memory","layersClipped%d" % i))
                except Exception as e:
                    arcpy.AddWarning("Could not clip layer %s to extent of catchment shapefile" % layersBasename[i])
                    raise(e)
                try:
                    arcpy.AddMessage ("Repairing layer %s" % layersBasename[i])
                    arcpy.RepairGeometry_management(layersClipped[-1])
                except Exception as e:
                    arcpy.AddWarning("Could not repair layer %s" % layersBasename[i])
            if len(layersClipped)>1:
                arcpy.AddMessage ("Processing union between %s and %s" % (layersBasename[0],layersBasename[1]))
                arcpy.Union_analysis (layersClipped[0:2], imparea, join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
                for i in range(len(layersClipped)-2):
                    arcpy.SetProgressor ("default", "Processing union for %s" % (layersBasename[i+2]))
                    arcpy.Union_analysis([imparea, layersClipped[i+2]], imparea + str(i), join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
                    imparea = imparea + str(i)
            else:
                arcpy.SetProgressor ("default", "Processing union for %s" % (layersBasename[0]))
                arcpy.Union_analysis(layersClipped[0], imparea, join_attributes="ALL", cluster_tolerance="", gaps="GAPS")

            arcpy.SetProgressor ("default", "Processing intersect analysis")
            arcpy.Intersect_analysis ([[imparea,2],[ms_catchments,1]], catchintsec)
            
            arcpy.CopyFeatures_management(catchintsec, "C:\Papirkurv\CatchIntSec.shp")
            
            hashes = []
            try:
                with arcpy.da.SearchCursor(catchintsec, ['SHAPE@WKT']) as catchmentcursor:
                    for row in catchmentcursor:
                        hashes.append(hash(row[0]))
            except RuntimeError as e:
                arcpy.AddError(imparea)
                arcpy.AddError(arcpy.Exists(imparea))
                arcpy.AddError(ms_catchments)
                arcpy.AddError(arcpy.Exists(ms_catchments))
                arcpy.AddError(arcpy.Exists(catchintsec))
                raise(e)
            duplicates = set([x for x in hashes if hashes.count(x) > 1])

            duplicatesDict = {}
            for dup in duplicates:
                duplicatesDict[dup] = 0

            with arcpy.da.UpdateCursor(catchintsec, ['SHAPE@WKT']) as catchmentcursor:
                for row in catchmentcursor:
                    if hash(row[0]) in duplicatesDict.keys():
                        if duplicatesDict[hash(row[0])] > 0:
                            catchmentcursor.deleteRow()
                        duplicatesDict[hash(row[0])] = duplicatesDict[hash(row[0])] + 1


            arcpy.AddField_management(in_table=catchintsec, field_name="shapearea", field_type="FLOAT", field_precision="20", field_scale="10", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
            arcpy.DefineProjection_management(in_dataset=catchintsec, coor_system="PROJCS['ETRS_1989_UTM_Zone_32N',GEOGCS['GCS_ETRS_1989',DATUM['D_ETRS_1989',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',9.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]")
            arcpy.CalculateField_management(catchintsec, "shapearea", "!SHAPE.AREA@SQUAREMETERS!", "PYTHON_9.3")
            arcpy.SpatialJoin_analysis(target_features=ms_catchments, join_features=catchintsec, out_feature_class=catchspatialjoin, join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_ALL", field_mapping="""MUID "MUID" true true false 40 Text 0 0 ,First,#,%s,MUID,-1,-1;shapearea "shapearea" true true false 4 Float 0 0 ,Sum,#,%s,shapearea,-1,-1""" % (ms_catchments,catchintsec), match_option="CONTAINS", search_radius="", distance_field_name="")
            arcpy.AddField_management(in_table=catchspatialjoin, field_name="ImpArea", field_type="FLOAT", field_precision="20", field_scale="10", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
            arcpy.CalculateField_management(catchspatialjoin, "ImpArea", "!shapearea!/!SHAPE.AREA@SQUAREMETERS!*100", "PYTHON_9.3")

            MUIDs = []
            ImpAreas = []
            with arcpy.da.SearchCursor(catchspatialjoin, ['MUID','ImpArea']) as catchmentcursor:
                for row in catchmentcursor:
                    MUIDs.append(row[0])
                    ImpAreas.append(row[1])
                    if False and type(row[1])==float:
                        arcpy.AddMessage("%s: %f" %(row[0],row[1]))
            catchmentsWithoutModelParameters = np.ones(len(MUIDs),dtype=bool)

            if ".sqlite" in arcpy.Describe(catchments).catalogPath:
                import sqlite3
                with sqlite3.connect((os.path.join(os.path.dirname(arcpy.Describe(catchments).catalogPath))).replace("!delete!","")) as connection:
                    update_cursor = connection.cursor()
                    for muid, imparea in zip(MUIDs,ImpAreas):
                        arcpy.AddMessage(imparea)
                        if imparea is None:
                            imparea = 0
                        try:
                            sql_execute = "UPDATE msm_Catchment SET ModelAImpArea = %1.2f WHERE MUID = '%s'" % (imparea/1e2, muid)
                            arcpy.AddMessage(sql_execute)
                            update_cursor.execute(
                                sql_execute)
                        except Exception as e:
                            arcpy.AddError("Can't change %s to %s" % (muid, imparea))
                            raise(e)
            else:
                if "OplandData_GDB" in arcpy.Describe(catchments).file:
                    catchmentcursor = arcpy.da.UpdateCursor(catchments, ["MUID","BEF_GRAD"])
                else:
                    catchmentcursor = arcpy.da.UpdateCursor(os.path.join(os.path.dirname(arcpy.Describe(catchments).path),"msm_HModA"), ["CatchID","ImpArea"])
                for row in catchmentcursor:
                    j = [k for k,a in enumerate(MUIDs) if a==row[0]]
                    if len(j)==1:
                        oldValue = row[1]
                        catchmentsWithoutModelParameters[j] = False
                        if type(ImpAreas[j[0]])==float:
                            row[1] = ImpAreas[j[0]]
                        else:
                            row[1] = 0
                        if abs(oldValue-row[1])>1:
                            arcpy.AddMessage("Changed catchment %s from %1.0f to %1.0f" % (row[0],oldValue,row[1]))
                            try:
                                catchmentcursor.updateRow(row)
                            except Exception as e:
                                arcpy.AddError("Can't change %s from %d to %d" % (row[0],oldValue,row[1]))
                                raise(e)

            del catchmentcursor
            if True in catchmentsWithoutModelParameters:
                arcpy.AddWarning("Catchments without Model Records [msm_HModA] (imperviousness not set): '" + "', '".join(np.array(MUIDs)[np.where(catchmentsWithoutModelParameters)]) + "'")
            #arcpy.CopyFeatures_management(catchspatialjoin,r"C:/Dokumenter/catchspatialjoin")
        return

class CatchmentProcessingAlt(object):
    def __init__(self):
        self.label       = "1b) Calculate imperviousness (set to field of shapefile)"
        self.description = "1b) Calculate imperviousness of catchments from shapefiles delineating impervious areas"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        # Input Features parameter
        imp_layer1 = arcpy.Parameter(
            displayName="Shapefile no. 1 delineating impervious areas",
            name="imp_layer1",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        imp_layer1.filter.list = ["Polygon"]
        
        imp_layer2 = arcpy.Parameter(
            displayName="Shapefile no. 2 delineating impervious areas",
            name="imp_layer2",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        imp_layer2.filter.list = ["Polygon"]

        imp_layer3 = arcpy.Parameter(
            displayName="Shapefile no. 3 delineating impervious areas",
            name="imp_layer3",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        imp_layer3.filter.list = ["Polygon"]

        catchments = arcpy.Parameter(
            displayName="Catchment layer",
            name="catchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchments.filter.list = ["Polygon"]

        field = arcpy.Parameter(
            displayName= "Field to assign imperviousness to",
            name="Field",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        parameters = [imp_layer1, imp_layer2, imp_layer3, catchments, field]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        if parameters[3].altered:
            parameters[4].filter.list = [f.name for f in arcpy.Describe(parameters[3].value).fields]
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        # imparea = os.path.join(arcpy.env.scratchGDB,"ImpArea")
        #catchintsec = os.path.join(arcpy.env.scratchGDB,"catchintsec")
        #catchspatialjoin = os.path.join(arcpy.env.scratchGDB,"CatchSpatialJoin")
        #ms_catchments = os.path.join(arcpy.env.scratchGDB, "ms_catchments")
        imparea = os.path.join("in_memory","ImpArea")
        catchintsec = os.path.join("in_memory","catchintsec")
        catchspatialjoin = os.path.join("in_memory","CatchSpatialJoin")
        ms_catchments = os.path.join("in_memory", "ms_catchments")
        field = parameters[4].ValueAsText
        arcpy.SetProgressor ("default", "Preparing catchments")
        tempField = "temp"
        arcpy.AddField_management(parameters[3].ValueAsText, "temp","SHORT")
        with arcpy.da.UpdateCursor(parameters[3].ValueAsText,"temp") as cursor:
            for i,row in enumerate(cursor):
                row[0] = int(i)
                cursor.updateRow(row)
        arcpy.CopyFeatures_management(parameters[3].ValueAsText,ms_catchments)

        arcpy.SetProgressor ("default", "Answer user promt (may be behind window)")
        userquery = pythonaddins.MessageBox("Calculate and assign imperviousness to %s catchments?" % (arcpy.GetCount_management(ms_catchments)[0]), "Confirm Assignment", 4)
        if userquery == "Yes":
            arcpy.AddMessage ("Repairing catchments layer")
            arcpy.RepairGeometry_management(ms_catchments)
            arcpy.DefineProjection_management(in_dataset=ms_catchments, coor_system="PROJCS['ETRS_1989_UTM_Zone_32N',GEOGCS['GCS_ETRS_1989',DATUM['D_ETRS_1989',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',9.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]")

            arcpy.AddMessage ("Preparing and combining shapefiles delineating impervious areas")
            #arcpy.SetProgressor ("default", "Clipping shapefiles delineating impervious areas to speed up geoprocessing calculation")
            extentPolygon = arcpy.CreateFeatureclass_management("in_memory", "extentPolygon", geometry_type = "POLYGON")
            extentCoordinates = arcpy.Describe(ms_catchments).extent
            array = arcpy.Array([arcpy.Point(extentCoordinates.XMin, extentCoordinates.YMin),
                         arcpy.Point(extentCoordinates.XMax, extentCoordinates.YMin),
                         arcpy.Point(extentCoordinates.XMax, extentCoordinates.YMax),
                         arcpy.Point(extentCoordinates.XMin, extentCoordinates.YMax),])
            polygon = arcpy.Polygon(array)
            with arcpy.da.InsertCursor(extentPolygon, ['SHAPE@']) as cursor:
                cursor.insertRow([polygon])

            layersBasename = []
            layers = []
            for i in range(2):
                if parameters[i].value:
                    layersBasename.append(arcpy.Describe(parameters[i].value).baseName)
                    if arcpy.Describe(parameters[i].value).dataElementType == "DEShapeFile":
                        layers.append(arcpy.MakeFeatureLayer_management(parameters[i].value,"in_memory\layer%d.lyr" % i))
                    else:
                        layers.append(parameters[i].value)

            layersClipped = []
            for i,layer in enumerate(layers):
                try:
                    # layersClipped.append(os.path.join("in_memory","layersClipped%d" % i))
                    arcpy.AddMessage ("Clipping layer %s" % layersBasename[i])
                    layersClipped.append(arcpy.SelectLayerByLocation_management(in_layer=layer, overlap_type="INTERSECT", select_features=extentPolygon, search_distance="", selection_type="NEW_SELECTION", invert_spatial_relationship="NOT_INVERT"))
                    # arcpy.Clip_analysis(layer, extentPolygon, os.path.join("in_memory","layersClipped%d" % i))
                except Exception as e:
                    arcpy.AddWarning("Could not clip layer %s to extent of catchment shapefile" % layersBasename[i])
                    raise(e)
                try:
                    arcpy.AddMessage ("Repairing layer %s" % layersBasename[i])
                    arcpy.RepairGeometry_management(layersClipped[-1])
                except Exception as e:
                    arcpy.AddWarning("Could not repair layer %s" % layersBasename[i])
            if len(layersClipped)>1:
                arcpy.AddMessage ("Processing union between %s and %s" % (layersBasename[0],layersBasename[1]))
                arcpy.Union_analysis (layersClipped[0:2], imparea, join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
                for i in range(len(layersClipped)-2):
                    arcpy.SetProgressor ("default", "Processing union for %s" % (layersBasename[i+2]))
                    arcpy.Union_analysis([imparea, layersClipped[i+2]], imparea + str(i), join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
                    imparea = imparea + str(i)
            else:
                arcpy.SetProgressor ("default", "Processing union for %s" % (layersBasename[0]))
                arcpy.Union_analysis(layersClipped[0], imparea, join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
            
            arcpy.SetProgressor ("default", "Processing intersect analysis")
            arcpy.Intersect_analysis ([[imparea,2],[ms_catchments,1]], catchintsec)

            hashes = []
            try:
                with arcpy.da.SearchCursor(catchintsec, ['SHAPE@WKT']) as catchmentcursor:
                    for row in catchmentcursor:
                        hashes.append(hash(row[0]))
            except RuntimeError as e:
                arcpy.AddError(imparea)
                arcpy.AddError(arcpy.Exists(imparea))
                arcpy.AddError(ms_catchments)
                arcpy.AddError(arcpy.Exists(ms_catchments))
                arcpy.AddError(arcpy.Exists(catchintsec))
                raise(e)
            duplicates = set([x for x in hashes if hashes.count(x) > 1])

            duplicatesDict = {}
            for dup in duplicates:
                duplicatesDict[dup] = 0

            with arcpy.da.UpdateCursor(catchintsec, ['SHAPE@WKT']) as catchmentcursor:
                for row in catchmentcursor:
                    if hash(row[0]) in duplicatesDict.keys():
                        if duplicatesDict[hash(row[0])] > 0:
                            catchmentcursor.deleteRow()
                        duplicatesDict[hash(row[0])] = duplicatesDict[hash(row[0])] + 1


            arcpy.AddField_management(in_table=catchintsec, field_name="shapearea", field_type="FLOAT", field_precision="20", field_scale="10", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
            arcpy.DefineProjection_management(in_dataset=catchintsec, coor_system="PROJCS['ETRS_1989_UTM_Zone_32N',GEOGCS['GCS_ETRS_1989',DATUM['D_ETRS_1989',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',0.0],PARAMETER['central_meridian',9.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]")
            arcpy.CalculateField_management(catchintsec, "shapearea", "!SHAPE.AREA@SQUAREMETERS!", "PYTHON_9.3")
            # arcpy.CopyFeatures_management(catchintsec,"P:\Temp\catchintsec.shp")
            # arcpy.CopyFeatures_management(ms_catchments,"P:\Temp\ms_catchments.shp")
            # arcpy.CopyFeatures_management(ms_catchments,"P:\Temp\ms_catchments2.shp")
            arcpy.SpatialJoin_analysis(target_features=ms_catchments, join_features=catchintsec, out_feature_class=catchspatialjoin, join_operation="JOIN_ONE_TO_ONE", join_type="KEEP_ALL", field_mapping="""temp "temp" true true false 40 Text 0 0 ,First,#,%s,temp,-1,-1;shapearea "shapearea" true true false 4 Float 0 0 ,Sum,#,%s,shapearea,-1,-1""" % (ms_catchments,catchintsec), match_option="CONTAINS", search_radius="", distance_field_name="")
            arcpy.AddField_management(in_table=catchspatialjoin, field_name="ImpArea", field_type="FLOAT", field_precision="20", field_scale="10", field_length="", field_alias="", field_is_nullable="NULLABLE", field_is_required="NON_REQUIRED", field_domain="")
            arcpy.CalculateField_management(catchspatialjoin, "ImpArea", "!shapearea!/!SHAPE.AREA@SQUAREMETERS!*100", "PYTHON_9.3")
            # arcpy.CopyFeatures_management(catchspatialjoin,"P:\Temp\ms_catchments.shp")
            
            IDColumn = "temp"
            if IDColumn not in [f.name for f in arcpy.ListFields(parameters[3].ValueAsText)]:
                arcpy.AddError("No index column named %s found on shapefile" % (IDColumn))
            OBJECTID = []
            ImpAreas = []
            with arcpy.da.SearchCursor(catchspatialjoin, [IDColumn,'ImpArea']) as catchmentcursor:
                for row in catchmentcursor:
                    OBJECTID.append(row[0])
                    ImpAreas.append(row[1])
                    if False and type(row[1])==float:
                        arcpy.AddMessage("%s: %f" %(row[0],row[1]))
            arcpy.AddMessage(OBJECTID)
            catchmentsWithoutModelParameters = np.ones(len(OBJECTID),dtype=bool)
            arcpy.AddMessage(ImpAreas)
            catchmentcursor = arcpy.da.UpdateCursor(parameters[3].ValueAsText, [IDColumn, field])
            for row in catchmentcursor:
                j = [k for k,a in enumerate(OBJECTID) if int(a)==row[0]]
                if len(j)==1:
                    oldValue = row[1]
                    catchmentsWithoutModelParameters[j] = False
                    if type(ImpAreas[j[0]])==float:
                        row[1] = ImpAreas[j[0]]
                    else:
                        row[1] = 0
                    if oldValue and abs(oldValue-row[1])>1:
                        arcpy.AddMessage("Changed catchment %s from %1.0f to %1.0f" % (row[0],oldValue,row[1]))
                    try:
                        catchmentcursor.updateRow(row)
                    except Exception as e:
                        arcpy.AddError("Can't change %s from %d to %d" % (row[0],oldValue,row[1]))
                        raise(e)
            arcpy.DeleteField_management(parameters[3].ValueAsText,"temp")
            del catchmentcursor

            #arcpy.CopyFeatures_management(catchspatialjoin,r"C:/Dokumenter/catchspatialjoin")
        return


class CatchmentProcessingScalgo(object):
    def __init__(self):
        self.label = "1c) Calculate imperviousness based on areal cover from SCALGO LIVE"
        self.description = "1c) Calculate imperviousness based on areal cover from SCALGO LIVE"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        catchments = arcpy.Parameter(
            displayName="Catchment layer",
            name="catchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchments.filter.list = ["Polygon"]

        scalgo_raster = arcpy.Parameter(
            displayName="Raster with areal cover",
            name="DHMFile",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        impervious_table = arcpy.Parameter(
            displayName='Imperviousness of type of areal cover',
            name='impervious_table',
            datatype='GPValueTable',
            parameterType='Required',
            multiValue="True",
            direction='Input')

        gridcode_dict = {1: 'Bar jord', 2: 'Vand', 3: 'Andet befaestet', 6: 'Lav vegetation', 7: 'Hoej vegetation',
                         8: 'Mark', 9: 'Befaestet vej', 10: 'Ubefaestet vej', 16: 'Bygning'}

        imperviousness_dict = {'Bar jord': 0, 'Vand': 0, 'Andet befaestet': 100, 'Lav vegetation': 0,
                               'Mark': 0, 'Hoej vegetation': 0, 'Befaestet vej': 100, 'Ubefaestet vej': 0, 'Bygning': 100}

        impervious_table.columns = [['GPLong', 'GridCode'], ['String', 'Feature'], ['Double', 'Imperviousness']]
        impervious_table.values = [[key, value, imperviousness_dict[value]] for key,value in gridcode_dict.iteritems()]

        field = arcpy.Parameter(
            displayName="Field to assign imperviousness to",
            name="Field",
            datatype="GPString",
            parameterType="Optional",
            category="Set to field instead",
            direction="Input")

        parameters = [catchments, scalgo_raster, impervious_table, field]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  # optional
        if parameters[0].altered:
            parameters[3].filter.list = [f.name for f in arcpy.Describe(parameters[0].value).fields]
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        ms_Catchment = parameters[0].ValueAsText
        raster = parameters[1].ValueAsText
        impervious_table = parameters[2].Value
        set_to_field = parameters[3].ValueAsText
        arcpy.env.overwriteOutput = True
        arcpy.AddMessage(impervious_table)
        mike_urban_database = os.path.dirname(arcpy.Describe(ms_Catchment).catalogPath)
        is_sqlite = True if ".sqlite" in mike_urban_database else False

        class Catchment:
            def __init__(self, area):
                self.area = abs(area)
                self.impervious_area = 0
                self.old_imperviousness = 0

            @property
            def imperviousness(self):
                return self.impervious_area / self.area * 100

        catchments = {row[0]: Catchment(row[1]) for row in arcpy.da.SearchCursor(ms_Catchment, ["MUID", "SHAPE@Area"])}
        arcpy.AddMessage([field.name for field in arcpy.ListFields(ms_Catchment)])
        if is_sqlite:
            with arcpy.da.SearchCursor(ms_Catchment, ["MUID", set_to_field if set_to_field else "modelaimparea"]) as cursor:
                for row in cursor:
                    catchments[row[0]].old_imperviousness = row[1]

        # ms_Catchment_copy = arcpy.management.MinimumBoundingGeometry(ms_Catchment, r"C:\Papirkurv\ms_Catchment2", geometry_type = "RECTANGLE_BY_AREA", group_option = "ALL")[0]
        # ms_Catchment_copy = arcpy.management.CopyFeatures(ms_Catchment, r"in_memory\ms_Catchment")[0]
        arcpy.SetProgressor("default", "Answer user promt (may be behind window)")
        arcpy.AddMessage("Answer user promt (may be behind window)")
        userquery = pythonaddins.MessageBox(
            "Assign imperviousness to %s catchments?" % (len(catchments)),
            "Confirm Assignment", 4)

        arcpy.AddMessage("MUIDs IN ('%s')" % ("', '".join(catchments.keys())))
        if userquery == "Yes":
            gridcode_dict = {row[0]: row[1] for row in impervious_table}
            imperviousness_dict = {row[1]: row[2] for row in impervious_table}

            gridcodes_impervious = [gridcode for gridcode, area_cover in gridcode_dict.iteritems() if
                                    imperviousness_dict[area_cover] > 0]

            # raster_class = arcpy.Raster(raster)
            #
            # raster_np_array = arcpy.RasterToNumPyArray(raster_class)
            #
            # lowerLeft = arcpy.Point(raster_class.extent.XMin, raster_class.extent.YMin)
            # cellSize = raster_class.meanCellWidth
            #
            # raster_np_array[np.isin(raster_np_array, gridcodes_impervious, invert = True)] = 0
            # raster_edited = arcpy.NumPyArrayToRaster(raster_np_array,lowerLeft,cellSize,
            #                          value_to_nodata=0)

            # raster_clipped = arcpy.management.Clip(raster, None, r"in_memory\raster_clipped", in_template_dataset = ms_Catchment_copy)[0]
            raster_polygon = arcpy.RasterToPolygon_conversion(in_raster=raster, out_polygon_features=r"in_memory\raster_to_polygon",
                                             simplify="SIMPLIFY", raster_field="Value",
                                             create_multipart_features="SINGLE_OUTER_PART", max_vertices_per_feature="")[0]

            # from alive_progress import alive_it
            with arcpy.da.UpdateCursor(raster_polygon, ["gridcode"], where_clause = "gridcode NOT IN (%s)" % (", ".join([str(gridcode) for gridcode in gridcodes_impervious]))) as cursor:
                for row in cursor:
                    cursor.deleteRow()
            arcpy.RepairGeometry_management(raster_polygon)

            arcpy.AddMessage((raster_polygon, arcpy.Describe(ms_Catchment).catalogPath))

            # arcpy.CopyFeatures_management(arcpy.Describe(ms_Catchment).catalogPath, "in_memory\ms_Catchment")
            arcpy.Select_analysis(arcpy.Describe(ms_Catchment).catalogPath,
                                  "in_memory\ms_Catchment", where_clause="MUID IN ('%s')" % "', '".join(
                    [muid for muid in catchments.keys()]))

            # Repair the geometry features, if needed
            try:
                arcpy.RepairGeometry_management("in_memory\ms_Catchment")
            except Exception as e:
                arcpy.AddWarning(e.message)

            intersect_polygon = arcpy.Intersect_analysis(in_features="%s #;%s #" % (raster_polygon, "in_memory\ms_Catchment"),
                                                         out_feature_class=r"in_memory\intersect", join_attributes="ALL",
                                                         cluster_tolerance="-1 Unknown", output_type="INPUT")[0]

            # arcpy.CopyFeatures_management(intersect_polygon, r"C:\Papirkurv\Bob.shp")

            with arcpy.da.SearchCursor(intersect_polygon, ["MUID", "gridcode", "SHAPE@Area"]) as cursor:
                for row in cursor:
                    if row[1] in gridcode_dict and imperviousness_dict[gridcode_dict[row[1]]] > 0:
                        try:
                            catchments[row[0]].impervious_area += imperviousness_dict[gridcode_dict[row[1]]] * row[2] / 100.0
                        except Exception as e:
                            arcpy.AddError(catchments.keys())
                            arcpy.AddError(imperviousness_dict.keys())
                            arcpy.AddError(row)
                            raise(e)

            arcpy.AddMessage(arcpy.Describe(ms_Catchment).catalogPath)
            in_decimal = False
            if "sqlite" in arcpy.Describe(ms_Catchment).catalogPath:
                import sqlite3
                arcpy.AddMessage(mike_urban_database.replace("!delete!",""))
                connection = None
                field = set_to_field if set_to_field else "modelaimparea"
                try:
                    connection =sqlite3.connect(mike_urban_database.replace("!delete!", ""))
                    update_cursor = connection.cursor()
                    for muid in catchments.keys():
                        if abs(catchments[muid].old_imperviousness - catchments[muid].imperviousness/100.0) > 0.01:
                            arcpy.AddMessage("UPDATE msm_Catchment SET %s = %1.2f WHERE MUID = '%s'" % (field, catchments[muid].imperviousness/100.0, muid))
                            update_cursor.execute(
                                                "UPDATE msm_Catchment SET %s = %1.2f WHERE MUID = '%s'" % (field, catchments[muid].imperviousness/100.0, muid))
                            arcpy.AddMessage("Changed catchment %s from %d to %d" % (muid, catchments[muid].old_imperviousness*1e2, catchments[muid].imperviousness))
                    connection.commit()
                    connection.close()
                except Exception as e:
                    import traceback
                    arcpy.AddWarning(traceback.format_exc())
                    raise(e)
                finally:
                    if connection:
                        connection.close()
            else:
                if "OplandData_GDB" in arcpy.Describe(ms_Catchment).file:
                    catchmentcursor = arcpy.da.UpdateCursor(ms_Catchment, ["MUID", "BEF_GRAD"],
                                                            where_clause="MUIDs IN ('%s')" % ("', '".join(catchments.keys())))
                elif (set_to_field.lower() if set_to_field else "ModelAImpArea".lower()) in [field.name.lower() for field in arcpy.ListFields(ms_Catchment)]:
                    in_decimal = False
                    catchmentcursor = arcpy.da.UpdateCursor(arcpy.Describe(ms_Catchment).catalogPath, ["muid", set_to_field if set_to_field else "modelaimparea"],
                                                            where_clause="MUID IN ('%s')" % ("', '".join(catchments.keys())))
                    # for catchment in catchments.values():
                    #     catchment.imperviousness = catchment.imperviousness/100.0
                else:
                    catchmentcursor = arcpy.da.UpdateCursor(
                        os.path.join(os.path.dirname(arcpy.Describe(ms_Catchment).path), "msm_HModA"), ["CatchID", set_to_field if set_to_field else "ImpArea"],
                        where_clause="CatchID IN ('%s')" % ("', '".join(catchments.keys())))

                for row in catchmentcursor:
                    oldValue = row[1]
                    row[1] = catchments[row[0]].imperviousness/100.0 if in_decimal else catchments[row[0]].imperviousness
                    if abs(oldValue - row[1]) > (0.01 if in_decimal else 1):
                        if in_decimal:
                            arcpy.AddMessage("Changed catchment %s from %1.0f to %1.0f" % (row[0], oldValue*1e2, row[1]*1e2))
                        else:
                            arcpy.AddMessage("Changed catchment %s from %1.0f to %1.0f" % (row[0], oldValue, row[1]))
                        try:
                            catchmentcursor.updateRow(row)
                        except Exception as e:
                            arcpy.AddError("Can't change %s from %1.2f to %1.2f" % (row[0], oldValue, row[1]))
                            raise (e)

                del catchmentcursor

        return

class AlteredCatchments(object):
    def __init__(self):
        self.label       = "Find altered catchments"
        self.description = "Compare two catchment databases and select the catchments in the first database that differ from the second"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        newcatchments = arcpy.Parameter(
            displayName="New catchments",
            name="newcatchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        newcatchments.filter.list = ["Polygon"]

        oldcatchments = arcpy.Parameter(
            displayName="Old catchments",
            name="oldcatchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        oldcatchments.filter.list = ["Polygon"]

        selectiontypeparameter = arcpy.Parameter(
           displayName = "Selection type",
           name = "selection_type",
           datatype = "GPString",
           parameterType = "Required",
           direction = "Input")
        selectiontypeparameter.value = "NEW_SELECTION"
        selectiontypeparameter.filter.list = ["NEW_SELECTION", "ADD_TO_SELECTION", "REMOVE_FROM_SELECTION", "SUBSET_SELECTION"]


        parameters = [newcatchments, oldcatchments, selectiontypeparameter]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        oldCatchmentsFile = parameters[1].ValueAsText
        newCatchmentsFile = parameters[0].ValueAsText
        selectiontype = parameters[2].ValueAsText

        oldCatchments = {}
        with arcpy.da.SearchCursor(oldCatchmentsFile, ['MUID','SHAPE@AREA']) as catchmentcursor:
            for row in catchmentcursor:
                oldCatchments[row[0]] = row[1]
        arcpy.AddMessage(oldCatchments)
        changedCatchments = []
        with arcpy.da.SearchCursor(newCatchmentsFile, ['MUID','SHAPE@AREA']) as catchmentcursor:
            for row in catchmentcursor:
                matchMUID = [row[0] for a in oldCatchments.keys() if row[0]==a]
                if not matchMUID:
                    changedCatchments.append(row[0])
                elif oldCatchments[matchMUID[0]] != row[1]:
                    changedCatchments.append(row[0])
        where_clause=r"""MUID IN (%s)""" % ("'" + "', '".join(changedCatchments) + "'")

        arcpy.AddMessage(where_clause)
        arcpy.AddMessage("%1.0f catchments changed (%1.1f%s)" % (len(changedCatchments),float(len(changedCatchments))/len(oldCatchments)*100,"%"))

        arcpy.SelectLayerByAttribute_management(in_layer_or_view="ms_Catchment", selection_type="NEW_SELECTION", where_clause = where_clause   )
        # arcpy.SelectLayerByAttribute_management(in_layer_or_view="ReduktionKnude", selection_type="NEW_SELECTION", where_clause = "NodeID IN ('000846R','000866R')")
        return

class CheckCatchments(object):
    def __init__(self):
        self.label       = "2) Check catchments"
        self.description = "2) Check catchments"

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
             # try:
        # MU_database.filter.list = ["hey"] 
        # mxd = arcpy.mapping.MapDocument("CURRENT")
        # df = arcpy.mapping.ListDataFrames(mxd)[0]
        # databases = set()
        # for layer in arcpy.mapping.ListLayers(df):
            # desc = arcpy.Describe(layer)
            # if ".mdb" in desc.catalogPath:
                # databases.add(desc.catalogPath.split(".mdb")[0] + ".mdb") 
        # MU_database.filter.list = list(databases)
        # except:
            # pass

        parameters = [MU_database]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        catchments = MU_database + "\mu_Geometry\ms_Catchment"
        catchCon = MU_database + "\msm_CatchCon"
        hModA = MU_database + "\msm_HModA"
        nodesShape = MU_database + "\msm_Node"

        MUIDs = []
        MUIDsBadName = []
        badNameRegex = re.compile("^Catchment")
        with arcpy.da.SearchCursor(catchments,["MUID"]) as cursor:
            for row in cursor:
                MUIDs.append(row[0])
                if badNameRegex.findall(row[0]):
                    MUIDsBadName.append(row[0])

        duplicates = set([x for x in MUIDs if MUIDs.count(x) > 1])

        if len(duplicates)>0:
            arcpy.AddError("Critical error: Duplicate catchments in table ms_Catchment found: '" + "', '".join(duplicates) + "'")

        if len(MUIDsBadName)>0:
            arcpy.AddWarning("Warning: Catchments with unoriginal names found (not suitable for FLORA): '" + "', '".join(MUIDsBadName) + "'")

        catchmentsSinglePart = arcpy.MultipartToSinglepart_management(catchments,
                                               r"in_memory\singlepartCatchments")

        MUIDsSinglePart = []
        with arcpy.da.SearchCursor(catchmentsSinglePart,["MUID"]) as cursor:
            for row in cursor:
                MUIDsSinglePart.append(row[0])

        singlePartDuplicates = set([x for x in MUIDsSinglePart if MUIDsSinglePart.count(x) > 1])
        multiPartFeatures = [a for a in singlePartDuplicates if not a in duplicates]

        if len(multiPartFeatures)>0:
            arcpy.AddWarning("Warning: Multipart catchments in table ms_Catchment found (unable to import catchments into Dandas): '" + "', '".join(multiPartFeatures) + "'")

        MUIDsWithoutModelRecords = set(MUIDs)
        catchIDs = []
        with arcpy.da.SearchCursor(hModA,["CatchID"]) as cursor:
            for row in cursor:
                catchIDs.append(row[0])
                if row[0] in MUIDsWithoutModelRecords:
                    MUIDsWithoutModelRecords.remove(row[0])

        if len(MUIDsWithoutModelRecords)>0:
            arcpy.AddError("Critical error: Catchments without model records found: '" + "', '".join(MUIDsWithoutModelRecords) + "'")

        catchIDsWithoutCatchments = []
        for catchID in set(catchIDs):
            if not catchID in MUIDs:
                catchIDsWithoutCatchments.append(catchID)

        if len(catchIDsWithoutCatchments):
            arcpy.AddWarning("Warning: Model records without catchments found: '" + "', '".join(catchIDsWithoutCatchments) + "'")
        catchIDsDuplicates = set([x for x in catchIDs if catchIDs.count(x) > 1])
        if len(catchIDsDuplicates)>0:
            arcpy.AddError("Critical error: Duplicate entries in table msm_HModA found: '" + "', '".join(catchIDsDuplicates) + "'")

        nodes = [row[0] for row in arcpy.da.SearchCursor(nodesShape,"MUID")]

        catchConIDs = []
        MUIDsWithoutConnection = set(MUIDs)
        MUIDsWithConnectionToDeletedNode = set(MUIDs)
        with arcpy.da.SearchCursor(catchCon,["CatchID","NodeID"]) as cursor:
            for row in cursor:
                catchConIDs.append(row[0])
                if row[0] in MUIDsWithoutConnection:
                    MUIDsWithoutConnection.remove(row[0])
                if row[1] in nodes:
                    try:
                        MUIDsWithConnectionToDeletedNode.remove(row[0])
                    except Exception as e:
                        pass

        if len(MUIDsWithoutConnection)>0:
            arcpy.AddError("Critical error: Catchments without connection to node found: '" + "', '".join(MUIDsWithoutConnection) + "'")

        if len(MUIDsWithConnectionToDeletedNode)>0:
            arcpy.AddError("Critical error: Catchments connected to non-existant node found: '" + "', '".join(MUIDsWithConnectionToDeletedNode) + "'")

        catchConIDsDuplicates = set([x for x in catchConIDs if catchConIDs.count(x) > 1])
        if len(catchConIDsDuplicates)>0:
            arcpy.AddError("Critical error: Duplicate entries in table msm_CatchCon found: '" + "', '".join(catchConIDsDuplicates) + "'")

        catchConsWithoutCatchments = []
        for catchConID in set(catchConIDs):
            if not catchConID in MUIDs:
                catchConsWithoutCatchments.append(catchConID)
        if len(catchConsWithoutCatchments)>0:
            arcpy.AddWarning("Warning: Connections in msm_CatchCon without catchments found: '" + "', '".join(catchConsWithoutCatchments) + "'")


        arcpy.AddMessage("Catchment check done.")

        return

class FAS2Deloplande(object):
    def __init__(self):
        self.label       = "b) Assign FAS Data to Catchments"
        self.description = "b) Assign FAS Data to Catchments"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        fas_features = arcpy.Parameter(
            displayName="FAS shape file",
            name="fas_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        fas_features.filter.list = ["Point"]

        catchment_features = arcpy.Parameter(
            displayName="Catchment shape file",
            name="catchment_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        catchment_features.filter.list = ["Polygon"]

        fields = arcpy.Parameter(
            displayName= "Field to assign FAS data to",
            name="Fields",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        matchtype = arcpy.Parameter(
            displayName= "Only match if FAS address is inside catchment",
            name="matchtype",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        matchtype.value = "false"

        taksttype = arcpy.Parameter(
            displayName= "Include the following taksttyper:",
            name="taksttype",
            datatype="GPString",
            multiValue="true",
            parameterType="Required",
            direction="Input")

        searchDistance = arcpy.Parameter(
            displayName= "Search distance for matching FAS-point to catchment:",
            name="searchDistance",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        searchDistance.value = 100.0

        taksttype.filter.list = [u'1/2 vandafledning - samletanke', u'Anden form for nedsivning', u'Biologisk', u'Biologisk rensning', u'Godkendt', u'Godkendt nedsivning', u'Off. kloak', u'Off. kloak - erhverv', u'Pumpestation- el ref.', u'Samletank', u'Samletank bidrag, Kjellerup', u'Samletank toilet', u'Samletank uden m\xe5ler', u'Samletank, m. m\xe5ler', u'Statsafg. 50% Aqua', u'Statsafgift nit.biologisk', u'Statsafgift spildevand', u'Statsafgift spildevand biologisk fa', u'Trappemodel', u'Trappemodel - samlet over 500 m\xb3', u'Trappemodel %', u'Trappemodel >20.000', u'Trappemodel over 20.000 m\xb3', u'Trappemodel over 20.000 m\xb3 - Fast enhed', u'Trappemodel over 500 m\xb3', u'Trappemodel uden rabat', u'Vandafl.', u'vandafl. bl. bolig & erhverv', u'vandafl. bl. bolig og erhverv', u'Vandafl.bidrag', u'Vandafledning']
        taksttype.value = [u'Off. kloak', u'Off. kloak - erhverv', u'Pumpestation- el ref.', u'Trappemodel', u'Trappemodel - samlet over 500 m\xb3', u'Trappemodel %', u'Trappemodel >20.000', u'Trappemodel over 20.000 m\xb3', u'Trappemodel over 20.000 m\xb3 - Fast enhed', u'Trappemodel over 500 m\xb3', u'Trappemodel uden rabat', u'Vandafl.', u'vandafl. bl. bolig & erhverv', u'vandafl. bl. bolig og erhverv', u'Vandafl.bidrag', u'Vandafledning', u'Vandafledning']
        parameters = [fas_features, catchment_features, fields, matchtype, taksttype, searchDistance]

        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[1].altered:
            parameters[2].filter.list = [f.name for f in arcpy.Describe(parameters[1].value).fields]
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        searchDistance = parameters[5].ValueAsText
        FASFile = parameters[0].ValueAsText
        arcpy.env.workspace = "in_memory/FASJoin"
        OIDs = []
        with arcpy.da.SearchCursor(FASFile,"FID") as fascursor:
            for row in fascursor:
                OIDs.append("%d" % row[0])

        arcpy.CopyFeatures_management (arcpy.Describe(parameters[1].ValueAsText).catalogPath, "catchmentFile")
        arcpy.DefineProjection_management("catchmentFile",arcpy.Describe(parameters[0].ValueAsText).spatialReference)

        catchments = 0
        with arcpy.da.SearchCursor("catchmentFile", ['SHAPE@']) as catchmentcursor:
            for row in catchmentcursor:
                catchments = catchments+1
        catchmentForbrug = np.zeros((catchments),dtype=np.float)
        arcpy.SetProgressor("step","Calculating FAS per catchment", 0, int(arcpy.GetCount_management(FASFile)[0]), 1)
        selectedFAS = 0
        nodesOutsideSearchDistance = []
        nodesOutsideSearchDistanceValue = []

        with arcpy.da.SearchCursor(FASFile, ['SHAPE@',"Forbrug", "taksttype","FID"]) as fascursor:
            for i, fasrow in enumerate(fascursor):
                arcpy.SetProgressorPosition(i)
                if fasrow[2] in parameters[4].ValueAsText:
                    geometry = fascursor[0]
                    selectedFAS = selectedFAS +1
                    dist = []
                    with arcpy.da.SearchCursor("catchmentFile", ['SHAPE@']) as catchmentcursor:
                        for row in catchmentcursor:
                            newgeometry = catchmentcursor[0]
                            dist.append(geometry.distanceTo(newgeometry))
                            if dist[-1] == 0:
                                break
                    idx = np.argmin(dist)
                    if parameters[3].ValueAsText == "false" or dist[idx]==0:
                        catchmentForbrug[idx] = catchmentForbrug[idx] + fasrow[1]
                        if searchDistance and float(searchDistance)<np.min(dist):
                            nodesOutsideSearchDistance.append("%d" % fasrow[3])
                            nodesOutsideSearchDistanceValue.append(np.min(dist))

        idxSort = np.flipud(np.argsort(nodesOutsideSearchDistanceValue))
        if len(nodesOutsideSearchDistance)>0:
            nodesOutsideSearchDistanceValue = np.array(nodesOutsideSearchDistanceValue)
            nodesOutsideSearchDistance = np.array(nodesOutsideSearchDistance)
            messageStr = "Do you want to assign FAS to catchments from these nodes:\n([OBJECTID]: distance to nearest catchment)\n"
            for nodei in idxSort:
                messageStr += "%s: %dm\n" % (nodesOutsideSearchDistance[nodei], nodesOutsideSearchDistanceValue[nodei])
            arcpy.SetProgressor ("default", "Confirm assignment")
            arcpy.SelectLayerByAttribute_management(in_layer_or_view=FASFile, selection_type="NEW_SELECTION", where_clause = "FID IN (" + ",".join(nodesOutsideSearchDistance) + ")")
            userquery = pythonaddins.MessageBox(messageStr, "Confirm Assignment", 3)

            if userquery == "Yes":
                arcpy.SetProgressor("step","Assigning FAS to catchments", 0, catchments, 1)
                with arcpy.da.UpdateCursor(parameters[1].ValueAsText, [parameters[2].ValueAsText]) as catchmentcursor:
                    for i,row in enumerate(catchmentcursor):
                        arcpy.SetProgressorPosition(i)
                        if catchmentForbrug[i] <> 0:
                            row[0] = catchmentForbrug[i]
                            catchmentcursor.updateRow(row)

            elif userquery == "No":
                for OID in nodesOutsideSearchDistance:
                    OIDs.remove(OID)
                catchmentForbrug = np.zeros((catchments),dtype=np.float)
                arcpy.SelectLayerByAttribute_management(in_layer_or_view=FASFile, selection_type="NEW_SELECTION", where_clause = "FID IN (" + ",".join(OIDs) + ")"  )
                with arcpy.da.SearchCursor(FASFile, ['SHAPE@',"Forbrug", "taksttype","FID"]) as fascursor:
                    for i, fasrow in enumerate(fascursor):
                        arcpy.SetProgressorPosition(i)
                        if fasrow[2] in parameters[4].ValueAsText:
                            geometry = fascursor[0]
                            dist = []
                            with arcpy.da.SearchCursor("catchmentFile", ['SHAPE@']) as catchmentcursor:
                                for row in catchmentcursor:
                                    newgeometry = catchmentcursor[0]
                                    dist.append(geometry.distanceTo(newgeometry))
                            idx = np.argmin(dist)
                            if parameters[3].ValueAsText == "false" or dist[idx]==0:
                                if searchDistance and float(searchDistance)>np.min(dist):
                                    catchmentForbrug[idx] = catchmentForbrug[idx] + fasrow[1]
                arcpy.SetProgressor("step","Assigning FAS to catchments", 0, catchments, 1)
                with arcpy.da.UpdateCursor(parameters[1].ValueAsText, [parameters[2].ValueAsText]) as catchmentcursor:
                    for i,row in enumerate(catchmentcursor):
                        arcpy.SetProgressorPosition(i)
                        if catchmentForbrug[i] <> 0:
                            row[0] = catchmentForbrug[i]
                            catchmentcursor.updateRow(row)
            else:
                arcpy.AddMessage("Cancelled FAS assignment")
        else:
            arcpy.SetProgressor("step","Assigning FAS to catchments", 0, catchments, 1)
            with arcpy.da.UpdateCursor(parameters[1].ValueAsText, [parameters[2].ValueAsText]) as catchmentcursor:
                for i,row in enumerate(catchmentcursor):
                    arcpy.SetProgressorPosition(i)
                    if catchmentForbrug[i] <> 0:
                        row[0] = catchmentForbrug[i]
                        catchmentcursor.updateRow(row)
        return

class TransferCatchments(object):
    def __init__(self):
        self.label       = "3) Transfer Catchments from shapefile to Mike Urban Model"
        self.description = "3) Transfer Catchments from shapefile to Mike Urban Model"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        ms_Catchment = arcpy.Parameter(
            displayName="Catchments Layer",
            name="ms_Catchment",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        ms_Catchment.filter.list = ["Polygon"]
        
        muidField = arcpy.Parameter(
            displayName="MUID field",
            name="muidField",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        
        impField = arcpy.Parameter(
            displayName="Imperviousness field",
            name="impField",
            datatype="String",
            parameterType="Optional",
            direction="Input")

        nodeID = arcpy.Parameter(
            displayName="Field indicating connection to node",
            name="nodeID",
            datatype="String",
            parameterType="Optional",
            direction="Input")

        MU_database = arcpy.Parameter(
            displayName="Mike Urban Model",
            name="mikeUrbanModel",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        parameters = [ms_Catchment, muidField, impField, nodeID, MU_database]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        # if parameters[0].altered:
            # parameters[1].filter.list = [f.name for f in arcpy.Describe(parameters[0].value).fields]
            # muidi = [i for i,a in enumerate(parameters[1].filter.list) if "muid" in a.lower()]
            # parameters[1].value = parameters[1].filter.list[muidi[0]] if muidi else ""
            
            # parameters[2].filter.list = [f.name for f in arcpy.Describe(parameters[0].value).fields]
            # impi = [i for i,a in enumerate(parameters[2].filter.list) if "imp" in a.lower()]
            # parameters[2].value = parameters[2].filter.list[impi[0]] if impi else ""
            
            # parameters[3].filter.list = [f.name for f in arcpy.Describe(parameters[0].value).fields]
            # nodei = [i for i,a in enumerate(parameters[3].filter.list) if "node" in a.lower()]
            # parameters[3].value = parameters[3].filter.list[nodei[0]] if nodei else ""
            
        if not parameters[4].valueAsText:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            database = [lyr.dataSource for lyr in arcpy.mapping.ListLayers(mxd) if not lyr.getSelectionSet() and lyr.isFeatureLayer and ".mdb" in lyr.dataSource]
            if database:
                database = re.findall(r"(.+)(?=\.mdb)", database[0])[0] + ".mdb"
                parameters[4].value = database
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        ms_Catchment = parameters[0].ValueAsText
        muidField = parameters[1].ValueAsText
        impField = parameters[2].ValueAsText
        nodeID = parameters[3].ValueAsText
        MU_database = parameters[4].ValueAsText
        
        # ms_CatchmentCopy = os.path.join("P:\Temp", "ms_catchments56435")
        ms_CatchmentCopy = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, "ms_catchments"))
        arcpy.CopyFeatures_management(ms_Catchment,ms_CatchmentCopy)
        arcpy.AddMessage(ms_CatchmentCopy)
        if not muidField == "MUID" and "MUID" in [f.name for f in arcpy.ListFields(ms_CatchmentCopy)]:
            arcpy.DeleteField_management(ms_CatchmentCopy, "MUID")
            arcpy.AlterField_management(ms_CatchmentCopy, muidField, "MUID")
        
        errorMessage = ""
        if not ".sqlite" in MU_database:
            edit = arcpy.da.Editor(MU_database)
            edit.startEditing(True,False)
            edit.startOperation()
            ms_CatchmentMUIDs = [row[0] for row in arcpy.da.SearchCursor(ms_CatchmentCopy,"MUID")]
            duplicateMUIDs = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"ms_Catchment"),"MUID",where_clause = "MUID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
            duplicateHModA = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_HModA"),"CatchID",where_clause = "CatchID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
            duplicateCatchCon = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_CatchCon"),"CatchID",where_clause = "CatchID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
            if duplicateMUIDs:
                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Layer in Mike Urban Database" % ("(', '".join(duplicateMUIDs))
            if duplicateHModA:
                errorMessage = errorMessage + "\n" if errorMessage else ""
                errorMessage += "Catchments with MUID ('%s') already exist in Model Records (msm_HModA) in Mike Urban Database" % ("(', '".join(duplicateHModA))
            if duplicateCatchCon:
                errorMessage = errorMessage + "\n" if errorMessage else ""
                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Connections (msm_CatchCon) in Mike Urban Database" % ("(', '".join(duplicateCatchCon))
            # pythonaddins.MessageBox("prutskid", "Confirm Assignment", 4)
            if errorMessage:
                arcpy.AddWarning(errorMessage)
            if not errorMessage or pythonaddins.MessageBox("%s\nTransfer catchments anyway?" % (errorMessage), "Confirm Transfer", 1) == "OK":
                arcpy.Append_management(ms_CatchmentCopy, os.path.join(MU_database,"ms_Catchment"),"NO_TEST")
                
                with arcpy.da.SearchCursor(ms_CatchmentCopy,["MUID",impField,nodeID]) as catchmentCursor:
                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_CatchCon"),["CatchID","NodeID","TypeNo"]) as cursor:
                        for row in catchmentCursor:
                            nID = 0 if not row[2] else row[2]
                            cursor.insertRow((row[0],nID,1))
                    catchmentCursor.reset()
                    del cursor
                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_HModA"),["CatchID","ImpArea","ParAID","LocalNo","ConcTime","RFactor","ILoss","CoeffNo","TACoeff", "TACurveID"]) as cursor:
                        for row in catchmentCursor:
                            iArea = 0.0 if not row[1] else row[1]
                            cursor.insertRow((row[0], iArea,"-DEFAULT-",0,7,0.9,0.0006,0,0.33, "TACurve1"))
        else:
            ms_CatchmentMUIDs = [row[0] for row in arcpy.da.SearchCursor(ms_CatchmentCopy,"MUID")]
            duplicateMUIDs = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_Catchment"),"muid",where_clause = "muid IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
            duplicateCatchCon = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_CatchCon"),"catchid",where_clause = "catchid IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
            if duplicateMUIDs:
                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Layer in Mike Urban Database" % ("(', '".join(duplicateMUIDs))
            if duplicateCatchCon:
                errorMessage = errorMessage + "\n" if errorMessage else ""
                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Connections (msm_CatchCon) in Mike Urban Database" % ("(', '".join(duplicateCatchCon))
            # pythonaddins.MessageBox("prutskid", "Confirm Assignment", 4)
            if errorMessage:
                arcpy.AddWarning(errorMessage)
            if not errorMessage or pythonaddins.MessageBox("%s\nTransfer catchments anyway?" % (errorMessage), "Confirm Transfer", 1) == "OK":
                with arcpy.da.SearchCursor(ms_CatchmentCopy, field_names = ["SHAPE@", "MUID",impField,nodeID]) as catchmentCursor:
                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_Catchment"), field_names = ["SHAPE@", "muid", "hydrologicalmodelno", "modelaimparea", "modelaparaid"]) as cursor:
                        for row in catchmentCursor:
                            cursor.insertRow((row[0], row[1], 1, row[2], "-DEFAULT-", ))
                    catchmentCursor.reset()
                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_CatchCon"),["catchid","nodeid","typeno"]) as cursor:
                        for row in catchmentCursor:
                            nID = 0 if not row[2] else row[2]
                            cursor.insertRow((row[0],nID,1))
                    catchmentCursor.reset()
                    # with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_HModA"),["CatchID","ImpArea","ParAID","LocalNo","ConcTime","RFactor","ILoss","CoeffNo","TACoeff"]) as cursor:
                        # for row in catchmentCursor:
                            # iArea = 0 if not row[1] else row[1]
                            # cursor.insertRow((row[0],iArea,"-DEFAULT-",0,7,0.9,0.0006,0,0.33))
        edit.stopOperation()
        edit.stopEditing(True)
        return
        
class GetImperviousness(object):
    def __init__(self):
        self.label       = "Display imperviousness from Mike Urban database"
        self.description = "Display imperviousness from Mike Urban database"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MU_database.filter.list = ["mdb","sqlite"]
        
        includeWasteWater = arcpy.Parameter(
            displayName="Include wastewater?",
            name="includeWasteWater",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        
        parameters = [MU_database, includeWasteWater]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):        
        includeWasteWater = parameters[1].ValueAsText
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        ms_CatchmentImp = arcpy.CopyFeatures_management (parameters[0].ValueAsText + r"\mu_Geometry\ms_Catchment", getAvailableFilename(arcpy.env.scratchGDB + "\ms_CatchmentImp")).getOutput(0)
        ms_CatchmentImpLayer = arcpy.MakeFeatureLayer_management(ms_CatchmentImp, arcpy.env.scratchGDB + "\ms_CatchmentImpLayer").getOutput(0).name
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID", join_table=parameters[0].ValueAsText + r"\msm_HModA", join_field="CatchID", fields="ImpArea")
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID", join_table=parameters[0].ValueAsText + r"\msm_CatchCon", join_field="CatchID", fields="NodeID")
        addLayer = arcpy.mapping.Layer(ms_CatchmentImpLayer)
        if not includeWasteWater:
            updatelayer.definitionQuery = 'NetTypeNo IS NULL OR NetTypeNo IN (2,3)'
        else: 
            updatelayer.definitionQuery = ''
        #arcpy.ApplySymbologyFromLayer_management (addLayer, "Template.lyr")
        #arcpy.RefreshTOC()
        #arcpy.RefreshActiveView()
      
        
        return
                
class GetImperviousnessFlora(object):
    def __init__(self):
        self.label       = "Display imperviousness from Flora Catchments"
        self.description = "Display imperviousness from Flora Catchments"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        FLORA_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        includeWasteWater = arcpy.Parameter(
            displayName="Include wastewater?",
            name="includeWasteWater",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        
        parameters = [FLORA_database, includeWasteWater]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):        
        includeWasteWater = parameters[1].ValueAsText
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        addLayer = arcpy.mapping.Layer(parameters[0].ValueAsText + r"\Oplandsrelaterede_data\OplandData_GDB")
        arcpy.mapping.AddLayer(df, addLayer)
        updatelayer = arcpy.mapping.ListLayers(mxd, "OplandData_GDB", df)[0]
        sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\OplandData_GDB_Symbology.lyr")
        arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        arcpy.AddMessage(str(addLayer.workspacePath))
        arcpy.AddMessage(str(addLayer.datasetName))
        updatelayer.replaceDataSource(str(addLayer.workspacePath), 'ACCESS_WORKSPACE', str(addLayer.datasetName))
        if not includeWasteWater:
            updatelayer.definitionQuery = 'NetTypeNo = 3'
        #arcpy.RefreshTOC()
        #arcpy.RefreshActiveView()
      
        
        return
 
class GetImperviousness(object):
    def __init__(self):
        self.label       = "Display imperviousness from Mike Urban database"
        self.description = "Display imperviousness from Mike Urban database"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MU_database.filter.list = ["mdb","sqlite"]
        
        includeWasteWater = arcpy.Parameter(
            displayName="Include wastewater?",
            name="includeWasteWater",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        
        parameters = [MU_database, includeWasteWater]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):        
        includeWasteWater = parameters[1].ValueAsText
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        ms_CatchmentImp = arcpy.CopyFeatures_management (parameters[0].ValueAsText + r"\mu_Geometry\ms_Catchment", getAvailableFilename(arcpy.env.scratchGDB + "\ms_CatchmentImp")).getOutput(0)
        ms_CatchmentImpLayer = arcpy.MakeFeatureLayer_management(ms_CatchmentImp, arcpy.env.scratchGDB + "\ms_CatchmentImpLayer").getOutput(0).name
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID", join_table=parameters[0].ValueAsText + r"\msm_HModA", join_field="CatchID", fields="ImpArea")
        arcpy.JoinField_management(in_data=ms_CatchmentImpLayer, in_field="MUID", join_table=parameters[0].ValueAsText + r"\msm_CatchCon", join_field="CatchID", fields="NodeID")
        addLayer = arcpy.mapping.Layer(ms_CatchmentImpLayer)
        arcpy.mapping.AddLayer(df, addLayer)
        updatelayer = arcpy.mapping.ListLayers(mxd, ms_CatchmentImpLayer, df)[0]
        sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchments W Imp Area.lyr")
        arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        updatelayer.replaceDataSource(str(addLayer.workspacePath), 'FILEGDB_WORKSPACE', str(addLayer.datasetName))
        if not includeWasteWater:
            updatelayer.definitionQuery = 'NetTypeNo IS NULL OR NetTypeNo IN (2,3)'
        else: 
            updatelayer.definitionQuery = ''
        #arcpy.ApplySymbologyFromLayer_management (addLayer, "Template.lyr")
        #arcpy.RefreshTOC()
        #arcpy.RefreshActiveView()
      
        
        return
                
class DuplicateCatchments(object):
    def __init__(self):
        self.label       = "4) Solve Catchments with identical MUIDs"
        self.description = "4) Solve Catchments with identical MUIDs"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        catchments = arcpy.Parameter(
            displayName="Catchment Layer",
            name="catchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchments.filter.list = ["Polygon"]

        parameters = [catchments]
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):        
        catchments = parameters[0].ValueAsText
        MU_database = os.path.dirname(arcpy.Describe(catchments).catalogPath).replace("\mu_Geometry", "")
        ms_Catchment = os.path.join(MU_database, "ms_Catchment")
        msm_HModA = os.path.join(MU_database, "msm_HModA")
        msm_CatchCon = os.path.join(MU_database, "msm_CatchCon")

        arcpy.AddMessage("WARNING: EXPLODE MULTIPART REMOVES CORRESPONDING MSM_HMODA ROWS.")

        if not arcpy.Exists(msm_CatchCon):
            arcpy.AddMessage("Assuming it is not part of a MIKE database as %s does not exist" % msm_CatchCon)
            MUIDs = [row[0] for row in arcpy.da.SearchCursor(catchments, ["MUID"])]
            MUIDs_duplicate = set()
            for MUID in MUIDs:
                if MUIDs.count(MUID) > 1:
                    MUIDs_duplicate.add(MUID)

            MUIDs_reassigned = {MUID: [] for MUID in MUIDs_duplicate}
            MUIDs_used = MUIDs

            for MUID_duplicate in MUIDs_duplicate:
                idx_of_duplicates = [i for i, MUID in enumerate(MUIDs) if MUID == MUID_duplicate][1:]

                for idx in idx_of_duplicates:
                    i = 2
                    new_MUID = MUID_duplicate + "s%d" % i
                    while new_MUID in MUIDs_used:
                        i += 1
                        new_MUID = MUID_duplicate + "s%d" % i

                    MUIDs_reassigned[MUID_duplicate].append(new_MUID)
                    MUIDs_used.append(new_MUID)

            MUIDs_reassigned_loop = {MUID: 0 for MUID in
                                     MUIDs_reassigned.keys()}  # dictionairy to check how many times the loop has found this MUID
            with arcpy.da.UpdateCursor(catchments, ['MUID'], "MUID IN ('%s')" % ("', '".join(MUIDs_duplicate))) as cursor:
                for row in cursor:
                    MUID = row[0]
                    if row[0] in MUIDs_reassigned:
                        if MUIDs_reassigned_loop[row[0]] > 0:
                            row[0] = MUIDs_reassigned[row[0]][MUIDs_reassigned_loop[row[0]]-1]
                            arcpy.AddMessage(row)
                            cursor.updateRow(row)
                        MUIDs_reassigned_loop[MUID] += 1

        else:
            MUIDs = [row[0] for row in arcpy.da.SearchCursor(catchments, ["MUID"])]
            # max_catchcon_MUID = np.max([row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database, 'ms_Catchment'), ["MUID"])])

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

            arcpy.SetProgressorLabel("Creating backup in %s" % (arcpy.env.scratchWorkspace))
            for feature_class in [ms_Catchment, msm_HModA, msm_CatchCon]:
                arcpy.Copy_management(feature_class, os.path.join(arcpy.env.scratchWorkspace, os.path.basename(feature_class)))

            # Initialize an empty set to store MUIDs that appear more than once
            MUIDs_duplicate = set()
            for MUID in MUIDs:
                # If the MUID appears more than once in the list, add it to the duplicates set
                if MUIDs.count(MUID) > 1:
                    MUIDs_duplicate.add(MUID)

            # Create a dictionary where keys are duplicate MUIDs and values are empty lists (to hold reassigned MUIDs)
            MUIDs_reassigned = {MUID: [] for MUID in MUIDs_duplicate}

            # Start with the same list of MUIDs, which will later include newly reassigned ones
            MUIDs_used = MUIDs

            # Loop through each duplicate MUID
            for MUID_duplicate in MUIDs_duplicate:
                # Find the index of all occurrences of the duplicate MUID (except the first one)
                idx_of_duplicates = [i for i, MUID in enumerate(MUIDs) if MUID == MUID_duplicate][1:]

                # For each duplicate occurrence, generate a new unique MUID
                for idx in idx_of_duplicates:
                    i = 2  # Start the suffix at 2 for the first duplicate
                    new_MUID = MUID_duplicate + "s%d" % i  # Create the first new MUID by appending 's2'

                    # Increment the suffix until the new MUID is unique in the list of used MUIDs
                    while new_MUID in MUIDs_used:
                        i += 1
                        new_MUID = MUID_duplicate + "s%d" % i

                    # Store the newly generated MUID in the reassignment dictionary
                    MUIDs_reassigned[MUID_duplicate].append(new_MUID)
                    # Add the new MUID to the list of used MUIDs to ensure uniqueness
                    MUIDs_used.append(new_MUID)

            # Dictionary to store rows from the msm_HModA table where MUIDs are duplicated
            msm_HModA_table = {}
            # Use an ArcPy SearchCursor to fetch all rows with CatchID matching the duplicate MUIDs
            with arcpy.da.SearchCursor(msm_HModA, '*') as cursor:
                for row in cursor:
                    msm_HModA_table[row[1]] = row  # Store rows with the MUID as the key
            arcpy.AddMessage((msm_HModA, "CatchID IN ('%s')" % ("', '".join(MUIDs_duplicate)), msm_HModA_table))

            # Dictionary to store rows from the msm_CatchCon table where MUIDs are duplicated
            msm_CatchCon_table = {}
            # Use a SearchCursor to fetch rows from msm_CatchCon where CatchID matches duplicate MUIDs
            with arcpy.da.SearchCursor(msm_CatchCon, '*',
                                       "CatchID IN ('%s')" % ("', '".join(MUIDs_duplicate))) as cursor:
                for row in cursor:
                    msm_CatchCon_table[row[2]] = row  # Store rows with the CatchID as the key

            # Dictionary to track how many times each duplicate MUID has been reassigned
            MUIDs_reassigned_loop = {MUID: 0 for MUID in MUIDs_reassigned.keys()}

            # Use an UpdateCursor to update the catchments table, replacing duplicate MUIDs with reassigned ones
            with arcpy.da.UpdateCursor(catchments, ['MUID'],
                                       "MUID IN ('%s')" % ("', '".join(MUIDs_duplicate))) as cursor:
                for row in cursor:
                    MUID = row[0]  # Current MUID in the row
                    if row[0] in MUIDs_reassigned:  # If the MUID is one of the duplicates
                        if MUIDs_reassigned_loop[row[0]] > 0:  # If it's been reassigned before
                            old_MUID = row[0]
                            new_MUID = MUIDs_reassigned[row[0]][MUIDs_reassigned_loop[row[0]] - 1]
                            row[0] = new_MUID  # Replace with the next reassigned MUID in the list
                            cursor.updateRow(row)  # Update the row with the new MUID
                            arcpy.AddMessage("Updated row from %s to %s in catchments table." % (old_MUID, new_MUID))
                        # Increment the reassignment counter for this MUID
                        MUIDs_reassigned_loop[MUID] += 1

            # For each duplicate MUID, assign the original MUID as the first element in its reassignment list
            for MUID in MUIDs_reassigned:
                MUIDs_reassigned[MUID][0] = MUID
            arcpy.AddMessage(msm_HModA_table)

            # Use an InsertCursor to add new rows to the msm_HModA table with the reassigned MUIDs
            with arcpy.da.InsertCursor(msm_HModA, '*') as cursor:
                for original_MUID in MUIDs_reassigned.keys():
                    for new_MUID in MUIDs_reassigned[original_MUID]:
                        row = list(msm_HModA_table[original_MUID])  # Copy the original row
                        row[1] = new_MUID  # Replace the MUID with the reassigned one
                        cursor.insertRow(row)  # Insert the new row with the updated MUID
                    arcpy.AddMessage("Inserted new row %s in msm_HModA table." % row)


            # Variable to track how many new MUIDs are created
            new_MUID_i = 0
            # Find the maximum MUID in the msm_CatchCon table (to generate new unique MUIDs)
            catch_con_max_MUID = np.max([row[0] for row in arcpy.da.SearchCursor(msm_CatchCon, ["MUID"])])

            # Use an InsertCursor to add new rows to msm_CatchCon with the reassigned MUIDs
            with arcpy.da.InsertCursor(msm_CatchCon, '*') as cursor:
                for original_MUID in MUIDs_reassigned.keys():
                    for new_MUID in MUIDs_reassigned[original_MUID]:
                        new_MUID_i += 1
                        temp_row = list(msm_CatchCon_table[original_MUID])  # Copy the original row
                        temp_row[2] = new_MUID  # Update the MUID to the reassigned one
                        temp_row[1] = catch_con_max_MUID + new_MUID_i  # Assign a new unique MUID value
                        cursor.insertRow(temp_row)  # Insert the updated row into the table
                        arcpy.AddMessage("Inserted new row %s into msm_CatchCon table." % temp_row)
      
        
        return

class GenerateCatchmentConnections(object):
    def __init__(self):
        self.label       = "a) Generate Catchment Connections"
        self.description = "a) Generate Catchment Connections"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        catchment_layer = arcpy.Parameter(
            displayName="Catchment feature layer",
            name="pipe_layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchment_layer.filter.list = ["Polygon"]

        delete_old_connections = arcpy.Parameter(
            displayName="Delete All Old Catchment Connections",
            name="delete_old_connections",
            datatype="Boolean",
            parameterType="optional",
            direction="Output")
        delete_old_connections.value = True

        catchcon_layer = arcpy.Parameter(
            displayName="Generate catchment connections in this layer (default is msm_CatchConLink)",
            name="catchment_layer",
            datatype="GPFeatureLayer",
            parameterType="optional",
            category="Additional Settings",
            direction="Input")

        parameters = [catchment_layer, delete_old_connections, catchcon_layer]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):


        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        catchment_layer = parameters[0].Value
        delete_old_connections = parameters[1].Value
        catchcon_layer = parameters[2].Value
        mu_database = os.path.dirname(arcpy.Describe(catchment_layer).catalogPath).replace("mu_Geometry","")
        
        MUIDs = [row[0] for row in arcpy.da.SearchCursor(catchment_layer,["MUID"])]
        arcpy.AddMessage("%d selected catchments" % (len(MUIDs)))
        arcpy.AddMessage(mu_database)
        msm_Node = os.path.join(mu_database, "msm_Node")
        ms_Catchment = os.path.join(mu_database, "ms_Catchment") if not ".sqlite" in os.path.dirname(arcpy.Describe(catchment_layer).catalogPath) else os.path.join(mu_database, "msm_Catchment")
        msm_CatchCon = os.path.join(mu_database, "msm_CatchCon")
        msm_CatchConLink = os.path.join(mu_database, "msm_CatchConLink") if not catchcon_layer else catchcon_layer

        where_clause = "MUID IN ('%s')" % "', '".join(MUIDs)
        if len(where_clause) > 2900:
            where_clause = ""
            arcpy.AddMessage("Query exceeds 2900 chars")

        muids_split = splitWhereclause(MUIDs)
        for muids in muids_split:
            where_clause = "MUID IN %s" % muids
            # arcpy.AddMessage(where_clause)
            catchments_coordinates = {row[0]:row[1] for row in arcpy.da.SearchCursor(ms_Catchment, ["MUID", "SHAPE@XY"],
                                        where_clause = where_clause)}

            arcpy.AddMessage(where_clause.replace("MUID","CatchID"))

            links = {row[0]:[row[1], row[2]] for row in arcpy.da.SearchCursor(msm_CatchCon, ["MUID", "CatchID", "NodeID"],
                                        where_clause = where_clause.replace("MUID","CatchID")) if row[1] in MUIDs}
            # arcpy.AddMessage(where_clause.replace("MUID","CatchID"))
            # arcpy.AddMessage(links)

            if not catchcon_layer:
                with arcpy.da.UpdateCursor(msm_CatchConLink, ["CatchConID"], where_clause = "CatchConID IN (%s)" % ", ".join([str(key) for key in links.keys()])) as cursor:
                    for row in cursor:
                        if row[0] in links.keys():
                            # arcpy.AddMessage(row)
                            cursor.deleteRow()

            nodes_MUIDs = [node for catchment, node in links.values()]

            nodes_coordinates = {row[0]:row[1] for row in arcpy.da.SearchCursor(msm_Node, ["MUID", "SHAPE@XY"],
                                        where_clause = "MUID IN ('%s')" % "', '".join(nodes_MUIDs) if where_clause else "")}

            # arcpy.AddMessage("Generating %d links" % (len(links)))
            arcpy.SetProgressor("step","Generating Catchment Connections", 0, len(links), 1)
            fields = ["shape@"]
            for field in ["MUID", "CatchConID","SHAPE_Length", "CatchID"]:
                if field.lower() in [f.name.lower() for f in arcpy.ListFields(msm_CatchConLink)]:
                    fields.append(field.lower())
            # arcpy.AddMessage(fields)
            arcpy.SetProgressor("step", "Inserting links", 0, len(links.keys()), 1)
            with arcpy.da.InsertCursor(msm_CatchConLink, fields) as cursor:
                for link_i, link in enumerate(links.keys()):
                    arcpy.SetProgressorPosition(link_i)
                    # arcpy.AddMessage(links[link])
                    catchment, node = links[link]
                    catchment_coordinate = catchments_coordinates[catchment]
                    if node not in nodes_coordinates:
                        arcpy.AddWarning("Could not find node %s for catchment %s" % (node, catchment))
                    else:
                        node_coordinate = nodes_coordinates[node]
                        shape = arcpy.Polyline(arcpy.Array([arcpy.Point(catchment_coordinate[0], catchment_coordinate[1]),
                                                     arcpy.Point(node_coordinate[0], node_coordinate[1])]))
                        # row = [shape, link_i, link, shape.length]
                        row = [None]*len(fields)
                        # arcpy.AddMessage([i for i,field in enumerate(fields) if field=="shape@"][0])
                        row[[i for i,field in enumerate(fields) if field=="shape@"][0]] = shape
                        if "muid" in fields:
                            row[[i for i, field in enumerate(fields) if field == "muid"][0]] = link_i
                        if "catchconid" in fields:
                            row[[i for i, field in enumerate(fields) if field == "catchconid"][0]] = link
                        if "shape_length" in fields:
                            row[[i for i, field in enumerate(fields) if field == "shape_length"][0]] = shape.length
                        if "catchid" in fields:
                            row[[i for i, field in enumerate(fields) if field == "catchid"][0]] = catchment
                        cursor.insertRow(row)

                
        return


class SetImperviousness(object):
    def __init__(self):
        self.label = "1d) Set Imperviousness"
        self.description = "1d) Set Imperviousness"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions
        catchments = arcpy.Parameter(
            displayName="Catchment layer",
            name="catchments",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchments.filter.list = ["Polygon"]

        imperviousness = arcpy.Parameter(
            displayName="Imperviousness:",
            name="imperviousness",
            datatype="double",
            parameterType="Required",
            direction="Input")

        parameters = [catchments, imperviousness]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        # imparea = os.path.join(arcpy.env.scratchGDB,"ImpArea")
        # catchintsec = os.path.join(arcpy.env.scratchGDB,"catchintsec")
        # catchspatialjoin = os.path.join(arcpy.env.scratchGDB,"CatchSpatialJoin")
        # ms_catchments = os.path.join(arcpy.env.scratchGDB, "ms_catchments")
        catchments = parameters[0].ValueAsText
        imperviousness = parameters[1].Value

        MUIDs = [row[0] for row in arcpy.da.SearchCursor(catchments, ["MUID"])]

        userquery = pythonaddins.MessageBox(
            "Assign imperviousness to %s catchments?" % (len(MUIDs)),
            "Confirm Assignment", 4)
        arcpy.AddMessage("MUIDs IN ('%s')" % ("', '".join(MUIDs)))
        if userquery == "Yes":
            if "OplandData_GDB" in arcpy.Describe(catchments).file:
                catchmentcursor = arcpy.da.UpdateCursor(catchments, ["MUID", "BEF_GRAD"], where_clause = "MUIDs IN ('%s')" % ("', '".join(MUIDs)))
            else:
                catchmentcursor = arcpy.da.UpdateCursor(
                    os.path.join(os.path.dirname(arcpy.Describe(catchments).path), "msm_HModA"), ["CatchID", "ImpArea"], where_clause = "CatchID IN ('%s')" % ("', '".join(MUIDs)))
            for row in catchmentcursor:
                oldValue = row[1]
                row[1] = imperviousness
                if abs(oldValue - row[1]) > 1:
                    arcpy.AddMessage("Changed catchment %s from %1.0f to %1.0f" % (row[0], oldValue, row[1]))
                    try:
                        catchmentcursor.updateRow(row)
                    except Exception as e:
                        arcpy.AddError("Can't change %s from %d to %d" % (row[0], oldValue, row[1]))
                        raise (e)

            del catchmentcursor
            # if True in catchmentsWithoutModelParameters:
            #     arcpy.AddWarning(
            #         "Catchments without Model Records [msm_HModA] (imperviousness not set): '" + "', '".join(
            #             np.array(MUIDs)[np.where(catchmentsWithoutModelParameters)]) + "'")
                # arcpy.CopyFeatures_management(catchspatialjoin,r"C:/Dokumenter/catchspatialjoin")
        return

class CatchmentSlopeAnalysis(object):
    def __init__(self):
        self.label = "5) Calculate Slope from Catchment to Sewage System"
        self.description = "5) Calculate Slope from Catchment to Sewage System"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        catchment_layer = arcpy.Parameter(
            displayName="Catchments Layer:",
            name="ms_Catchment",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchment_layer.filter.list = ["Polygon"]

        msm_Link_field = arcpy.Parameter(
            displayName="Field from msm_Link to use:",
            name="msm_Link_field",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        msm_Link_field.filter.list = ["InvertLevel", "GroundLevel", "CriticalLevel"]
        msm_Link_field.value = "InvertLevel"

        raster = arcpy.Parameter(
            displayName="Terrain Raster Path:",
            name="raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        initial_depth = arcpy.Parameter(
            displayName="Initial depth (e.g. below frost line) [m]:",
            name="initial_depth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        initial_depth.value = 1.0

        resolution = arcpy.Parameter(
            displayName="Resolution of analysis (cell size) [m]:",
            name="resolution",
            datatype="GPLong",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        resolution.value = 8

        required_slope = arcpy.Parameter(
            displayName="Required Slope [m/m]:",
            name="required_slope",
            datatype="GPDouble",
            category="Calculate Vertical Difference",
            parameterType="Optional",
            direction="Input")
        required_slope.value = 20e-3

        add_debug_output = arcpy.Parameter(
            displayName="Create feature class with connections:",
            name="add_debug_output",
            datatype="GPBoolean",
            category="Additional Settings",
            parameterType="Optional",
            direction="Input")

        result_layer = arcpy.Parameter(
            displayName="Get Elevation from Result Layer:",
            name="result_layer",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            category="Get Elevation from Result Layer",
            direction="Input")
        result_layer.filter.list = ["Point", "Polyline"]

        result_layer_field = arcpy.Parameter(
            displayName="Field from Result Layer to use:",
            name="result_layer_field",
            datatype="GPString",
            parameterType="Required",
            category="Get Elevation from Result Layer",
            direction="Input")
        # msm_Link_field.filter.list = ["InvertLevel", "GroundLevel", "CriticalLevel"]
        # msm_Link_field.value = "InvertLevel"

        parameters = [catchment_layer, msm_Link_field, raster, initial_depth, resolution, required_slope, add_debug_output, result_layer, result_layer_field]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        result_layer = parameters[7].Value
        if result_layer:
            parameters[1].enabled = False
            parameters[8].enabled = True
            fields = [f.name for f in arcpy.Describe(parameters[7].value).fields]
            if "up_maxe" in [field.lower() for field in fields]:
                parameters[8].Value = "Up_MaxE & Dw_MaxE"
            else:
                parameters[8].filter.list = fields
        else:
            parameters[1].enabled = True
            parameters[8].enabled = False
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        catchment_layer = parameters[0].ValueAsText
        msm_Link_field = parameters[1].ValueAsText
        raster = parameters[2].ValueAsText
        initial_depth = parameters[3].Value
        resolution = parameters[4].Value
        required_slope = parameters[5].Value
        add_debug_output = parameters[6].Value
        result_layer = parameters[7].ValueAsText
        result_layer_field = parameters[8].ValueAsText

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]

        arcpy.env.overwriteOutput = True

        # Define input data
        MU_database = os.path.dirname(arcpy.Describe(catchment_layer).catalogPath).replace("\mu_Geometry", "")
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_Catchment = catchment_layer


        catchments_output = r"in_memory\catchments_export"
        clipped_raster = "in_memory\DTM_clipped"
        # clipped_raster = "in_memory\DTM"
        clipped_raster_resampled = "in_memory\DTM_c"
        raster_points = "in_memory\DTM_clipped_resampled_points"
        raster_points_mapped = "in_memory\DTM_clipped_resampled_points_mapped"
        # raster_points = "C:\Papirkurv\DTM_clipped_resampled_points.shp"
        # catchment_loop = "in_memory\Catchment_loop"

        arcpy.management.CopyFeatures(msm_Catchment, catchments_output)
        arcgis_pro = False

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
                        arcpymapping.AddLayer(df, layer, "TOP")
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

        arcpy.SetProgressor("default", "Creating MIKE Folder")
        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchWorkspace), "MIKE URBAN")
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
        arcpy.env.scratchWorkspace = os.path.join(MIKE_gdb, "scratch.gdb")

        output_raster_path = getAvailableFilename(os.path.join(arcpy.env.scratchWorkspace, "Slope"))
        output_vertical_add_raster_path = getAvailableFilename(os.path.join(arcpy.env.scratchWorkspace, "VerticalAdd"))

        arcpy.SetProgressor("default", "Graphing MIKE Database")
        # graphing MU_database
        import mikegraph
        graph = mikegraph.Graph(MU_database)
        graph.map_network()

        arcpy.SetProgressor("default", "Catchment Analysis")
        catchments_MUID = [row[0] for row in arcpy.da.SearchCursor(msm_Catchment, ["MUID"])]
        arcpy.AddMessage("Analyzing %d Catchments" % len(catchments_MUID))
        catchments = [graph.catchments_dict[muid] for muid in catchments_MUID]
        arcpy.Select_analysis(arcpy.Describe(catchment_layer).catalogPath,
                              catchments_output, where_clause="MUID IN ('%s')" % "', '".join([catchment.MUID for catchment in catchments]))

        arcpy.RepairGeometry_management(catchments_output)
        arcpy.AddField_management(catchments_output, 'NodeID', 'TEXT')

        with arcpy.da.UpdateCursor(catchments_output, ["MUID", "NodeID"]) as cursor:
            for row in cursor:
                for catchment in catchments:
                    if catchment.MUID == row[0]:
                        row[1] = catchment.nodeID
                        cursor.updateRow(row)
                        break

        arcpy.SetProgressor("default", "Raster Preparation")
        arcpy.management.Clip(raster, "#", clipped_raster, in_template_dataset=catchments_output,
                              clipping_geometry="ClippingGeometry", nodata_value=-999)
        arcpy.management.Resample(
            in_raster=clipped_raster,
            out_raster=clipped_raster_resampled,
            cell_size="%d %d" % (resolution, resolution),
            resampling_type="NEAREST"
        )

        arcpy.SetProgressor("default", "Converting Raster to Points")
        raster = arcpy.Raster(clipped_raster_resampled)
        raster_array = arcpy.RasterToNumPyArray(raster)
        cell_size = raster.meanCellHeight  # Assuming square cells
        lower_left = arcpy.Point(raster.extent.XMin, raster.extent.YMin)

        arcpy.CreateFeatureclass_management(
            out_path=os.path.dirname(raster_points),
            out_name=os.path.basename(raster_points),
            geometry_type='POINT'
            # spatial_reference=raster.spatialReference
        )
        arcpy.AddField_management(raster_points, 'grid_code', 'DOUBLE')
        arcpy.AddField_management(raster_points, 'oldID', 'FLOAT', field_scale=0, field_precision=15)

        points_slope = np.zeros(raster_array.shape, dtype=float).flatten()
        points_slope.fill(-999)

        points_vertical_add = np.zeros(raster_array.shape, dtype=float).flatten()
        points_vertical_add.fill(-999)
        i = 0
        arcpy.SetProgressor("step", "Converting Raster to Points", 0, int(raster_array.size), 1)
        # arcpy.conversion.RasterToPoint(raster, raster_points, "VALUE")

        rows, cols = raster_array.shape
        x_coords = lower_left.X + (np.arange(cols) * cell_size) + (cell_size / 2)
        y_coords = lower_left.Y + (np.arange(rows) * cell_size) + (cell_size / 2)
        y_coords = np.flipud(y_coords)  # Flip Y to match raster indexing

        # Flatten arrays for vectorized processing
        x_grid, y_grid = np.meshgrid(x_coords, y_coords)
        valid_mask = raster_array<1e100  # Mask for valid values
        flat_x = x_grid[valid_mask]
        flat_y = y_grid[valid_mask]
        flat_values = raster_array[valid_mask]

        # Batch insert points
        point_data = [
            (arcpy.PointGeometry(arcpy.Point(x, y), raster.spatialReference), value, i)
            for i, (x, y, value) in enumerate(zip(flat_x, flat_y, flat_values)) if value != -999]

        with arcpy.da.InsertCursor(raster_points, ['SHAPE@', 'grid_code', "oldID"]) as cursor:
            for row_i, row in enumerate(point_data):
                cursor.insertRow(row)
                if row_i % 30 == 0:
                    arcpy.SetProgressorPosition(i)

        class Node:
            def __init__(self, muid, shape, critical_level):
                self.muid = muid
                self.shape = shape
                self.critical_level = critical_level


        arcpy.SetProgressor("default", "Getting Catchments and Nodes")
        nodes = {}
        with arcpy.da.SearchCursor(msm_Node, ["muid", "SHAPE@", msm_Link_field]) as cursor:
            for row in cursor:
                nodes[row[0]] = Node(row[0], row[1], row[2])

        class LinkResult:
            def __init__(self, muid, up_max_elevation, dw_max_elevation):
                self.muid = muid
                self.up_max_elevation = up_max_elevation
                self.dw_max_elevation = dw_max_elevation

        link_results = {}
        # if result layer, set it to node critical level
        if result_layer:
            if result_layer_field.lower() == "Up_MaxE & Dw_MaxE".lower():
                with arcpy.da.SearchCursor(arcpy.Describe(result_layer).catalogPath, ["MUID", "Up_MaxE", "Dw_MaxE"]) as cursor:
                    for row in cursor:
                        link_results[row[0]] = LinkResult(*row)
            else:
                with arcpy.da.SearchCursor(arcpy.Describe(result_layer).catalogPath, ["MUID", result_layer_field]) as cursor:
                    for row in cursor:
                        if row[0] in nodes.keys():
                            nodes[row[0]].critical_level = row[1]

        with arcpy.da.SearchCursor(catchments_output, ["OID@", "MUID"]) as cursor:
            for row in cursor:
                for catchment in catchments:
                    if row[1] == catchment.MUID:
                        catchment.OID = row[0]

        # arcpy.AddSpatialIndex_management(raster_points)

        arcpy.SetProgressor("step", "Converting Raster to Points", 0, int(raster_array.size), 1)
        arcpy.SetProgressorPosition(i)

        if add_debug_output:
            debug_output_path = getAvailableFilename(os.path.join(arcpy.env.scratchWorkspace, "CatchmentConnections"))
            arcpy.management.CreateFeatureclass(arcpy.env.scratchWorkspace, os.path.basename(debug_output_path), geometry_type = "POLYLINE", spatial_reference = raster.spatialReference)[0]
            debug_fields = ["Length", "UpLevel", "DwLevel", "Slope"]
            for field in debug_fields:
                arcpy.management.AddField(debug_output_path ,field, "FLOAT", field_precision = 8, field_scale = 4)
            debug_cursor = arcpy.da.InsertCursor(debug_output_path, ["SHAPE@"] + debug_fields)


        arcpy.MakeFeatureLayer_management(catchments_output, "Catchment Layer")

        # arcpy.AddSpatialIndex_management(catchments_output)
        # arcpy.AddSpatialIndex_management(raster_points)
        arcpy.SetProgressor("default", "Spatial Join of raster points and catchments")
        raster_points_name = arcpy.Describe(raster_points)
        catchments_output_name = arcpy.Describe(catchments_output)

        # FIELD MAPPINGS
        field_mappings = arcpy.FieldMappings()

        # Function to add a field with a specific merge rule
        def add_field_mapping(field_mappings, input_fc, field_name, merge_rule, output_name=None, field_type=None,
                              field_length=None):
            field_map = arcpy.FieldMap()
            field_map.addInputField(input_fc, field_name)

            # Set merge rule
            field_map.mergeRule = merge_rule

            # Modify output field properties if specified
            if output_name or field_type or field_length:
                field = field_map.outputField
                if output_name:
                    field.name = output_name
                    field.aliasName = output_name
                if field_type:
                    field.type = field_type  # 'TEXT', 'FLOAT', 'LONG'
                if field_length and field_type == "TEXT":
                    field.length = field_length  # Set text length explicitly
                field_map.outputField = field

            field_mappings.addFieldMap(field_map)

        # Add fields from raster_points (target)
        add_field_mapping(field_mappings, raster_points, "grid_code", "FIRST")
        add_field_mapping(field_mappings, raster_points, "oldID", "FIRST")

        add_field_mapping(field_mappings, catchments_output, "NodeID", "FIRST", field_type="TEXT", field_length=50)

        arcpy.analysis.SpatialJoin(raster_points, catchments_output, raster_points_mapped, join_operation = "JOIN_ONE_TO_ONE", match_option = "INTERSECT", join_type = "KEEP_COMMON", field_mapping = field_mappings)

        arcpy.SetProgressor("step", "Slope Analysis", 0, int(raster_array.size), 1)
        with arcpy.da.SearchCursor(raster_points_mapped, ["SHAPE@", "oldID", "grid_code", "NodeID"]) as cursor:
            for row_i, row in enumerate(cursor):
                node_id = row[3]
                if len(node_id)>0:
                    links = [link for link in graph.network.links.values() if
                             link.fromnode == node_id or link.tonode == node_id]

                    links_3d = []
                    for link in links:
                        if result_layer_field.lower() == "Up_MaxE & Dw_MaxE".lower():
                            links_3d.append(link.shape_3d(uplevel=link_results[link.MUID].up_max_elevation,
                                                          dwlevel=link_results[link.MUID].dw_max_elevation))
                        else:
                            links_3d.append(link.shape_3d(uplevel=nodes[link.fromnode].critical_level,
                                                          dwlevel=nodes[link.tonode].critical_level))
                    shortest_distance = 1e9
                    nearest_point = None
                    for link in links_3d:
                        point, _, distance, _ = link.queryPointAndDistance(row[0])
                        if distance < shortest_distance:
                            shortest_distance = distance
                            nearest_point = point
                    frostfri_dybde = initial_depth

                    if shortest_distance > 1 and nearest_point is not None:
                        points_slope[int(row[1])] = (row[2] - frostfri_dybde - nearest_point.centroid.Z) / shortest_distance * 1e3

                        points_vertical_add[int(row[1])] = required_slope * min(shortest_distance, 10) + 10e-3 * max(0,
                                                                                                                     shortest_distance - 10) - (
                                                                   row[
                                                                       2] - frostfri_dybde - nearest_point.centroid.Z)
                        if add_debug_output:
                            debug_cursor.insertRow([arcpy.Polyline(arcpy.Array([row[0].centroid, nearest_point.centroid])), shortest_distance, row[2], nearest_point.centroid.Z, points_slope[int(row[1])]])
                if row_i % 10 == 0:
                    arcpy.SetProgressorPosition(row_i)

        empty_group_mapped = arcpymapping.LayerFile(os.path.dirname(
            os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr") if arcgis_pro else arcpy.mapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = df.addLayer(empty_group_mapped) if arcgis_pro else arcpymapping.AddLayer(df, empty_group_mapped,
                                                                                               "TOP")
        empty_group_layer = df.listLayers('Empty Group')[0] if arcgis_pro else \
            arcpymapping.ListLayers(mxd, "Empty Group", df)[0]
        # if ".gdb" in nodes_featureclass or ".gdb" in reaches_featureclass:
        import datetime
        empty_group_layer.name = "Catchment Analysis (%s)" % (datetime.datetime.now().strftime("%H:%M"))

        if add_debug_output:
            del debug_cursor
            addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchment Connections.lyr", debug_output_path,
                     workspace_type="FILEGDB_WORKSPACE", new_name = "Catchment Link", group = empty_group_layer)

        arcpy.SetProgressor("default", "Creating Rasters")
        slope_raster = arcpy.NumPyArrayToRaster(
            in_array=points_slope.reshape(raster_array.shape),
            lower_left_corner=lower_left,
            x_cell_size=cell_size,
            y_cell_size=cell_size,
            value_to_nodata=-999  # Set this if you have a specific NoData value
        )
        slope_raster.save(output_raster_path)

        vertical_add_raster = arcpy.NumPyArrayToRaster(
            in_array=(points_vertical_add.reshape(raster_array.shape)),
            lower_left_corner=lower_left,
            x_cell_size=cell_size,
            y_cell_size=cell_size,
            value_to_nodata=-999  # Set this if you have a specific NoData value
        )
        vertical_add_raster.save(output_vertical_add_raster_path)

        arcpy.SetProgressor("default", "Adding to Layout")

        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchment_slope.lyr",
                 output_raster_path, workspace_type="FILEGDB_WORKSPACE", group = empty_group_layer)

        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Catchment_elevation_adjustment.lyr",
                 output_vertical_add_raster_path, workspace_type="FILEGDB_WORKSPACE", group = empty_group_layer)

        return

class CatchmentConnectionDescription(object):
    def __init__(self):
        self.label = "c) Assign NetTypeNo to Catchment Connections"
        self.description = "c) Assign NetTypeNo to Catchment Connections"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        catchment_layer = arcpy.Parameter(
            displayName="Catchments Layer:",
            name="ms_Catchment",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        catchment_layer.filter.list = ["Polygon"]

        assign_type = arcpy.Parameter(
            displayName="Assign using the following method:",
            name="assign_type",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        assign_type.filter.list = ["NetTypeNo of Catchment", "NetTypeNo of NodeID of Catchment", "Parameter Set of Catchment"]

        parameters = [catchment_layer, assign_type]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional

        return

    def execute(self, parameters, messages):
        catchment_layer = parameters[0].ValueAsText
        assign_type = parameters[1].ValueAsText
        MU_database = os.path.dirname(arcpy.Describe(catchment_layer).catalogPath).replace("\mu_Geometry", "")
        is_sqlite = True if ".sqlite" in MU_database else False
        if not is_sqlite:
            raise(Exception("Database Type is not supported"))
        msm_Node = os.path.join(MU_database, "msm_Node")
        msm_CatchCon = os.path.join(MU_database, "msm_CatchCon") if is_sqlite else os.path.join(MU_database, "msm_CatchConLink")

        if assign_type.lower() == "NetTypeNo of Catchment".lower():
            raise(Exception("Method not currently supported"))
            catchments_nettypeno = {row[0]:row[1] for row in arcpy.da.SearchCursor(catchment_layer, ["MUID", "NETTYPENO"])}
        elif assign_type.lower() == "NetTypeNo of NodeID of Catchment".lower():
            raise (Exception("Method not currently supported"))
        elif assign_type.lower() == "Parameter Set of Catchment".lower():
            arcpy.AddMessage([row for row in
                                    arcpy.da.SearchCursor(catchment_layer, ["MUID", "ModelAParAID"])])
            catchments_HParAID = {row[0]: row[1] for row in
                                    arcpy.da.SearchCursor(catchment_layer, ["MUID", "ModelAParAID"])}
            import sqlite3
            with sqlite3.connect(MU_database.replace("!delete!","")) as connection:
                update_cursor = connection.cursor()
                for muid, description in catchments_HParAID.iteritems():
                    sql_execute = "UPDATE msm_CatchCon SET Description = '%s' WHERE CatchID = '%s'" % (
                            description, muid)
                    update_cursor.execute(
                        sql_execute)
        else:
            raise(Exception("Method not found"))


        return