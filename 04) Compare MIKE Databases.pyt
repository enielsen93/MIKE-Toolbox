import arcpy
import os
import re
import numpy as np
if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
    from arcpy._mapping import Layer
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument

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

def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
             definition_query=None):
    if arcgis_pro:
        mxd = arcpy.mp.ArcGISProject("CURRENT")
        df = mxd.listMaps()[0]
    else:
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
    if arcgis_pro and not ".lyrx" in layer_source and os.path.exists(layer_source.replace(".lyr", ".lyrx")):
        layer_source = layer_source.replace(".lyr", ".lyrx")
    arcpy.AddMessage(source)
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

        if not arcgis_pro: update_layer = arcpymapping.listLayers(mxd, source_layer.name, df)[0] if arcgis_pro else \
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
                df.addLayerToGroup(group, layer, "BOTTOM")
            else:
                arcpymapping.AddLayerToGroup(df, group, layer, "BOTTOM")
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
                import copy
                cp = copy.deepcopy(update_layer.connectionProperties)
                arcpy.AddMessage(cp)
                if workspace_type == "FILEGDB_WORKSPACE" or workspace_type == "FILEGDB":
                    workspace_type = "File Geodatabase"
                cp["connection_info"]['database'] = os.path.dirname(
                    source.replace(r"\mu_Geometry", ""))
                arcpy.AddMessage(cp["connection_info"]['database'])
                cp['dataset'] = os.path.basename(source)
                cp['workspace_factory'] = workspace_type
                arcpy.AddMessage(("😀", update_layer.connectionProperties))
                arcpy.AddMessage(("😗", cp))
                update_layer.updateConnectionProperties(update_layer.connectionProperties, cp)
                arcpy.AddMessage(update_layer.connectionProperties)
            else:
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                               workspace_type, os.path.basename(source))
    return update_layer

class Toolbox(object):
    def __init__(self):
        self.label = "Compare MIKE Models"
        self.alias = "Compare MIKE Models"
        self.canRunInBackground = True
        # List of tool classes associated with this toolbox
        self.tools = [CompareMikeModels, FixSimulationModelName, CompareMikeModelsLabels]


class CompareMikeModels(object):
    def __init__(self):
        self.label = "1) Compare MIKE Models"
        self.description = "1) Compare MIKE Models"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        database1 = arcpy.Parameter(
            displayName="Reference MIKE Model",
            name="database1",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        database1.filter.list = ["sqlite", "mdb"]

        database2 = arcpy.Parameter(
            displayName="MIKE Model to compare to reference model",
            name="database2",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        database2.filter.list = ["sqlite", "mdb"]

        check_features = arcpy.Parameter(
            displayName="Check these features:",
            name="check_features",
            datatype="GPString",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        check_features.filter.list = ["msm_Catchment", "msm_Node", "msm_Link", "msm_Weir", "msm_Orifice", "msm_Project", "msm_HParA", "msm_BBoundary"]
        check_features.value = ["msm_Catchment", "msm_Node", "msm_Link", "msm_Weir", "msm_Orifice"]

        ignore_fields = arcpy.Parameter(
            displayName="Ignore following common fields:",
            name="ignore_fields",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        ignore_fields.filter.list = ["OBJECTID", "SHAPE", "Slope", "UpLevel_C", "DwLevel_C", "Length_C", "UpLevel", "DwLevel", "Diameter", "NetTypeNo", "GroundLevel", "InvertLevel", "CriticalLevel", "Area", "Description", "AssetName", "Fricno"]
        ignore_fields.value = ["OBJECTID", "Slope", "UpLevel_C", "DwLevel_C", "Length_C", "CriticalLevel", "Area"]
        #ignore_fields.value = ["msm_Catchment", "msm_Node", "msm_Link", "msm_Weir", "msm_Orifice"]

        parameters = [database1, database2, check_features, ignore_fields]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        database1 = parameters[0].ValueAsText
        database2 = parameters[1].ValueAsText
        check_features = parameters[2].ValueAsText.split(";")
        ignore_fields = [field.lower() for field in parameters[3].ValueAsText.split(";")]

        if "slope" in ignore_fields:
            ignore_fields.append("slope_c")
        if "shape" in ignore_fields:
            ignore_fields.append("shape@")

        #arcpy.AddMessage(ignore_fields)

        def ignore_field(fieldname):
            if fieldname.lower() in ignore_fields:
                return True
            else:
                return False

        if arcgis_pro:
            mxd = arcpy.mp.ArcGISProject("CURRENT")
            df = mxd.listMaps()[0]
        else:
            mxd = arcpymapping.MapDocument("CURRENT")
            df = arcpymapping.ArcGISProject("current").activeMap

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE"):
            layer = arcpy.mapping.Layer(layer_source)
            if group:
                arcpy.mapping.AddLayerToGroup(df, group, layer, "BOTTOM")
            else:
                arcpy.mapping.AddLayer(df, layer, "BOTTOM")
            updatelayer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
            updatelayer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))), workspace_type,
                                          unicode(os.path.basename(source)))
            return updatelayer

        if arcgis_pro:
            # Load the .lyr file
            layer_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), r"Data\EmptyGroup.lyr")

            # Add it to the map
            df.addDataFromPath(layer_path)

            # Rename it
            for lyr in df.listLayers():
                if lyr.name == "Empty Group":
                    lyr.name = "{} vs. {}".format(os.path.basename(database1).replace(".sqlite","_sqlite"), os.path.basename(database2).replace(".sqlite","_sqlite"))
                    empty_group_layer = lyr
                    break
        else:
            empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
            empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
            empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
            empty_group_layer.name = "%s vs. %s" % (os.path.basename(database1), os.path.basename(database2))
        
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

        for feature in check_features:
            arcpy.AddMessage("Checking feature %s" % (feature))
            feature_path_1 = os.path.join(database1, feature)
            if ".mdb" in feature_path_1 and "Catchment" in feature_path_1:
                feature_path_1 = feature_path_1.replace("msm_Catchment", "ms_Catchment")
            feature_path_2 = os.path.join(database2, feature)
            if ".mdb" in feature_path_2 and "Catchment" in feature_path_2:
                feature_path_2 = feature_path_2.replace("msm_Catchment", "ms_Catchment")

            result_layer = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, feature))

            if arcpy.Describe(feature_path_1).dataType == "FeatureClass":
                arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, os.path.basename(result_layer),
                                                template=feature_path_1)
            else:
                arcpy.CreateTable_management(arcpy.env.scratchGDB, os.path.basename(result_layer),
                                                template=feature_path_1)

            arcpy.AddField_management(result_layer, "fields", "TEXT", field_is_nullable="NULLABLE")

            features_1 = {}
            features_2 = {}
            fields_path_1 = [field.name.lower() if field.name != "geometry" else "SHAPE@" for field in arcpy.ListFields(feature_path_1)
                      if not ignore_field(field.name)]
            fields_path_2 = [field.name.lower() if field.name != "geometry" else "SHAPE@" for field in
                             arcpy.ListFields(feature_path_2)
                             if not ignore_field(field.name)]

            fields = list(set(fields_path_1) & set(fields_path_2))

            if "msm_Catchment" in feature_path_1:
                fields.append("SHAPE@AREA")

            MUID_field_i = [i for i, field in enumerate(fields) if field.lower() == "muid"][0]

            if "SHAPE@" not in fields and arcpy.Describe(feature_path_1).dataType == "FeatureClass":
                fields.append("SHAPE@")

            arcpy.SetProgressor("step","Checking feature %s" % (feature_path_1), 0, len([1 for row in arcpy.da.SearchCursor(feature_path_1, "muid")]), 1)
            with arcpy.da.SearchCursor(feature_path_1, fields) as cursor:
                step = 0
                for row in cursor:
                    arcpy.SetProgressorPosition(step)
                    step += 1
                    features_1[row[MUID_field_i]] = row

            arcpy.SetProgressor("step","Checking feature %s" % (feature_path_2), 0, len([1 for row in arcpy.da.SearchCursor(feature_path_2, "muid")]), 1)
            arcpy.AddMessage((feature_path_2, fields))
            with arcpy.da.SearchCursor(feature_path_2, fields) as cursor:
                step = 0
                for row in cursor:
                    arcpy.SetProgressorPosition(step)
                    step += 1
                    features_2[row[MUID_field_i]] = row
            
            arcpy.SetProgressor("default", "Getting MUIDs")
            # Check MUIDs
            MUIDs = set(features_1.keys()) | set(features_2.keys())

            missing_MUIDs = set(features_1).symmetric_difference(features_2)
            MUIDs_to_check = [MUID for MUID in MUIDs if MUID not in missing_MUIDs]
            features_changed = {}

            def compare_rows(row1, row2, fields = None):
                fields_diff = []
                try:
                    for field_i in range(len(row1)):
                        if fields and fields[field_i] == "SHAPE@AREA":
                            if abs(abs(row1[field_i]) - abs(row2[field_i])) > 1:
                                fields_diff.append(field_i)
                        elif row1[field_i] != row2[field_i]:
                            fields_diff.append(field_i)
                except Exception as e:
                    arcpy.AddMessage(row1)
                    arcpy.AddMessage(row2)
                    arcpy.AddMessage(field_i)
                return fields_diff

            arcpy.SetProgressor("step","Comparing rows for feature %s" % (feature), 0, len(MUIDs_to_check), 1)
            MUIDs_field_changed = {}
            for step, MUID in enumerate(MUIDs_to_check):
                if step % 30 == 0:
                    arcpy.SetProgressorPosition(step)
                idx = compare_rows(features_1[MUID], features_2[MUID], fields = fields)
                if idx:
                    MUIDs_field_changed[MUID] = [fields[i] for i in idx]

            if "catchment" in feature.lower():
                msm_CatchCon_1 = {}
                msm_CatchCon_fields = [field.name for field in arcpy.ListFields(os.path.join(database1, "msm_CatchCon"))
                                       if not ignore_field(field.name) and not field.name == "MUID"]
                catchID_field_i = [i for i, field in enumerate(msm_CatchCon_fields) if field.lower() == "catchid"][0]
                with arcpy.da.SearchCursor(os.path.join(database1, "msm_CatchCon"), msm_CatchCon_fields) as cursor:
                    for row in cursor:
                        msm_CatchCon_1[row[catchID_field_i]] = row


                msm_CatchCon_2 = {}
                msm_CatchCon_fields = [field.name for field in arcpy.ListFields(os.path.join(database2, "msm_CatchCon")) if not ignore_field(field.name) and not field.name == "MUID"]
                catchID_field_i = [i for i, field in enumerate(msm_CatchCon_fields) if field.lower() == "catchid"][0]
                with arcpy.da.SearchCursor(os.path.join(database2, "msm_CatchCon"), msm_CatchCon_fields) as cursor:
                    for row in cursor:
                        msm_CatchCon_2[row[catchID_field_i]] = row

                for MUID in MUIDs_to_check:
                    if MUID in msm_CatchCon_1 and MUID in msm_CatchCon_2:
                        idx = compare_rows(msm_CatchCon_1[MUID], msm_CatchCon_2[MUID])
                        if idx:
                            MUIDs_field_changed[MUID] = [msm_CatchCon_fields[i] for i in idx]

            if feature == "msm_Catchment" and ".mdb" in feature_path_1:
                msm_HModA_1 = {}
                msm_HModA_fields = [field.name for field in arcpy.ListFields(os.path.join(database1, "msm_HModA")) if not ignore_field(field.name)]
                catchID_field_i = [i for i, field in enumerate(msm_HModA_fields) if field.lower() == "catchid"][0]

                with arcpy.da.SearchCursor(os.path.join(database1, "msm_HModA"), msm_HModA_fields) as cursor:
                    for row in cursor:
                        msm_HModA_1[row[catchID_field_i]] = row

                msm_HModA_2 = {}
                msm_HModA_fields = [field.name for field in arcpy.ListFields(os.path.join(database2, "msm_HModA")) if not ignore_field(field.name)]
                catchID_field_i = [i for i, field in enumerate(msm_HModA_fields) if field.lower() == "catchid"][0]

                with arcpy.da.SearchCursor(os.path.join(database2, "msm_HModA"), msm_HModA_fields) as cursor:
                    for row in cursor:
                        msm_HModA_2[row[catchID_field_i]] = row

                for MUID in MUIDs_to_check:
                    if MUID in msm_HModA_1 and MUID in msm_HModA_2:
                        idx = compare_rows(msm_HModA_1[MUID], msm_HModA_2[MUID])
                        if idx:
                            MUIDs_field_changed[MUID] = [msm_HModA_fields[i] for i in idx]

            with arcpy.da.InsertCursor(result_layer, fields + ["fields"]) as cursor:
                for missing_MUID in missing_MUIDs:
                    if missing_MUID in features_1.keys():
                        row = list(features_1[missing_MUID]) + ["Not in DB2"]
                        cursor.insertRow(row)
                    else:
                        row= list(features_2[missing_MUID]) + ["Not in DB1"]
                        cursor.insertRow(row)

            #arcpy.AddMessage(fields)
            if arcpy.Describe(feature_path_1).dataType == "FeatureClass":
                geometry_field_i = [field_i for field_i, field in enumerate(fields) if field.lower() == "shape@"][0]
                MUID_field_i = [field_i for field_i, field in enumerate(fields) if field.lower() == "muid"][0]
                with arcpy.da.InsertCursor(result_layer, ["MUID", "SHAPE@", "fields"]) as cursor:
                    for MUID in MUIDs_field_changed.keys():
                        row = (features_1[MUID][MUID_field_i], features_1[MUID][geometry_field_i],
                               ", ".join(MUIDs_field_changed[MUID]))
                        cursor.insertRow(row)
            else:
                MUID_field_i = [field_i for field_i, field in enumerate(fields) if field.lower() == "muid"][0]
                with arcpy.da.InsertCursor(result_layer, ["MUID", "fields"]) as cursor:
                    for MUID in MUIDs_field_changed.keys():
                        row = (features_1[MUID][MUID_field_i], ", ".join(MUIDs_field_changed[MUID]))
                        cursor.insertRow(row)

            if arcpy.Describe(feature_path_1).dataType == "FeatureClass":
                newlayer = arcpy.management.MakeFeatureLayer(result_layer)[0]
                if arcgis_pro:
                    for label_class in (newlayer.listLabelClasses()):
                        label_class.expression = "$feature.fields"
                    newlayer.showLabels = True
                    newlayer.name = newlayer.name + " (%d features)" % (np.sum(
                        [1 for row in arcpy.da.SearchCursor(result_layer, ["MUID"])]))
                    update_layer = df.addLayerToGroup(empty_group_layer, newlayer, "TOP")
                else:
                    for label_class in (newlayer.listLabelClasses() if arcgis_pro else newlayer.labelClasses):
                        label_class.expression = "[fields]"
                    newlayer.showLabels = True
                    newlayer.name = newlayer.name + " (%d features)" % (np.sum(
                        [1 for row in arcpy.da.SearchCursor(result_layer, ["MUID"])]))
                    update_layer = arcpy.mapping.AddLayerToGroup(df, empty_group_layer, newlayer, "TOP")
            else:
                newlayer = arcpy.management.MakeFeatureLayer(result_layer)[0]
                newlayer.name = newlayer.name + " (%d features)" % (np.sum(
                    [1 for row in arcpy.da.SearchCursor(result_layer, ["MUID"])]))
                update_layer = df.addLayerToGroup(empty_group_layer, newlayer, "TOP")

        return

class FixSimulationModelName(object):
    def __init__(self):
        self.label = "a) Fix MIKE Urban Model Name"
        self.description = "a) Fix MIKE Urban Model Name"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions
        
        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Optional",
            direction="Input")

        folder = arcpy.Parameter(
            displayName="Folder",
            name="folder",
            datatype="Folder",
            parameterType="Optional",
            direction="Input")
            
        pattern = arcpy.Parameter(
            displayName="String to rename:",
            name="pattern",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
            
        replacement= arcpy.Parameter(
            displayName="Replacement for string:",
            name="replacement",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        
        parameters = [MU_database, folder, pattern, replacement]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if not parameters[0].Value and not parameters[1].Value:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            layers = arcpy.mapping.ListLayers(mxd)
            for layer in layers:
                try:
                    if ".mdb" in layer.workspacePath:
                        MU_database = layer.workspacePath
                        break
                except Exception as e:
                    pass
            parameters[0].Value = MU_database
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        folder = parameters[1].ValueAsText
        pattern = parameters[2].ValueAsText
        replacement = parameters[3].ValueAsText


        if MU_database:
            msm_Project = os.path.join(MU_database, "msm_Project")
            if not pattern:
                model_name, model_version = re.findall(r'\\(?:Model)*(?:MIKE)*\\([\w]+)_(\d+)\\', MU_database)[0]
                # model_name = "KOM"
                arcpy.AddMessage((model_name, model_version))
                with arcpy.da.UpdateCursor(msm_Project, ["MUID", "PRFHotStartFileName", "CRFFileName", "CatchmentDischargeFileName"]) as cursor:
                    for row in cursor:
                        for val_i, val in enumerate(row):
                            # arcpy.AddMessage((r"LIS_008", model_name + "_" + model_version, val))
                            arcpy.AddMessage((r"%s_\d+" % model_name, r"%s_%s" % (model_name, model_version) , val))
                            if val:
                                val = re.sub("%s_\d+" % model_name, r"%s_%s" % (model_name, model_version) , val)
                            row[val_i] = val
                            arcpy.AddMessage((row[val_i], val))
                        cursor.updateRow(row)
            else:
                with arcpy.da.UpdateCursor(msm_Project, ["MUID", "PRFHotStartFileName", "CRFFileName", "CatchmentDischargeFileName"]) as cursor:
                    for row in cursor:
                        for val_i, val in enumerate(row):
                            # arcpy.AddMessage((r"LIS_008", model_name + "_" + model_version, val))
                            if val:
                                val = re.sub(pattern, replacement, val)
                            row[val_i] = val
                            arcpy.AddMessage((row[val_i], val))
                        cursor.updateRow(row)

        if folder:
            for root, dirs, files in os.walk(folder):
                for file in files:
                    if re.search(pattern, file):
                        new_filename = re.sub(pattern, replacement, file)
                        arcpy.AddMessage((os.path.join(root, file), os.path.join(root, new_filename)))
                        try:
                            os.rename(os.path.join(root, file), os.path.join(root, new_filename))
                        except Exception as e:
                            arcpy.AddMessage(e)


class CompareMikeModelsLabels(object):
    def __init__(self):
        self.label = "2) Display MIKE Manhole and Pipe Changes"
        self.description = "2) Display MIKE Manhole and Pipe Changes"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        database1 = arcpy.Parameter(
            displayName="Reference MIKE Model",
            name="database1",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        database1.filter.list = ["sqlite", "mdb"]

        database2 = arcpy.Parameter(
            displayName="MIKE Model to compare to reference model",
            name="database2",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        database2.filter.list = ["sqlite", "mdb"]

        parameters = [database1, database2]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        db1 = parameters[0].ValueAsText
        db2 = parameters[1].ValueAsText

        fc_name = "msm_Link"
        id_field = "MUID"
        diam_field = "Diameter"

        if arcgis_pro:
            mxd = arcpy.mp.ArcGISProject("CURRENT")
            df = mxd.listMaps()[0]
        else:
            mxd = arcpymapping.MapDocument("CURRENT")
            df = arcpymapping.ArcGISProject("current").activeMap

        if arcgis_pro:
            # Load the .lyr file
            layer_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), r"Data\EmptyGroup.lyr")

            # Add it to the map
            df.addDataFromPath(layer_path)

            # Rename it
            for lyr in df.listLayers():
                if lyr.name == "Empty Group":
                    lyr.name = "{} vs. {}".format(os.path.basename(db1).replace(".sqlite","_sqlite"), os.path.basename(db2).replace(".sqlite","_sqlite"))
                    empty_group_layer = lyr
                    break
        else:
            empty_group_mapped = arcpy.mapping.Layer(
                os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
            empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
            empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
            empty_group_layer.name = "%s vs. %s" % (os.path.basename(db1), os.path.basename(db2))

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

        nodes_output_filepath = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, "Nodes"))
        links_output_filepath = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, "Links"))

        # Load db1 into dict: {ID: (geometry, diameter)}
        db1_fc = os.path.join(db1, fc_name)
        db2_fc = os.path.join(db2, fc_name)

        links_db1 = {}
        with arcpy.da.SearchCursor(db1_fc, [id_field, "SHAPE@", "Diameter", "MaterialID", "nettypeno"]) as cursor:
            for row in cursor:
                links_db1[row[0]] = (row[1], row[2], row[3], row[4])

        # Prepare output shapefile
        spatial_ref = arcpy.Describe(db1_fc).spatialReference
        arcpy.CreateFeatureclass_management(
            out_path=os.path.dirname(links_output_filepath),
            out_name=os.path.basename(links_output_filepath),
            geometry_type="POLYLINE",
            spatial_reference=spatial_ref
        )
        arcpy.AddField_management(links_output_filepath, "Change", "TEXT", field_length=30)
        arcpy.AddField_management(links_output_filepath, "NetTypeNo", "SHORT")

        def renameMaterial(materialid):
            if "plastic" in materialid.lower():
                return "pl"
            elif "concrete" in materialid.lower():
                return "bt"

        # Compare and write differences
        with arcpy.da.SearchCursor(db2_fc, [id_field, diam_field, "materialID", "nettypeno"]) as cursor2, \
                arcpy.da.InsertCursor(links_output_filepath, ["SHAPE@", "Change", "NetTypeNo"]) as insert:

            for row in cursor2:
                link_id, diam2, material2, nettypeno2 = row
                if link_id in links_db1:
                    geom1, diam1, material1, nettypeno1 = links_db1[link_id]
                    if diam1 != diam2 or material1 != material2:
                        text = u"ø{}{}→ø{}{}".format(int(diam1 * 1000), renameMaterial(material1), int(diam2 * 1000), renameMaterial(material2))
                        insert.insertRow((geom1, text, nettypeno1))

        # arcpy.AddMessage("✅ Done. Output saved to: %s" % links_output_filepath)
        # arcpy.AddMessage(links_output_filepath)
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\CompareMIKEModels_Links.lyrx",
                links_output_filepath, group=empty_group_layer, workspace_type = "FILEGDB")

        # --- Compare msm_Node invert levels ---
        node_fc_name = "msm_Node"
        invert_field = "InvertLevel"

        db1_nodes = os.path.join(db1, node_fc_name)
        db2_nodes = os.path.join(db2, node_fc_name)

        # Load db1 nodes into dict: {ID: (geometry, invert)}
        nodes_db1 = {}
        with arcpy.da.SearchCursor(db1_nodes, [id_field, "SHAPE@", invert_field, "Diameter", "nettypeno"]) as cursor:
            for row in cursor:
                nodes_db1[row[0]] = (row[1], row[2], row[3], row[4])

        # Prepare output shapefile
        spatial_ref_nodes = arcpy.Describe(db1_nodes).spatialReference
        arcpy.CreateFeatureclass_management(
            out_path=os.path.dirname(nodes_output_filepath),
            out_name=os.path.basename(nodes_output_filepath),
            geometry_type="POINT",
            spatial_reference=spatial_ref_nodes
        )
        arcpy.AddField_management(nodes_output_filepath, "Change", "TEXT", field_length=30)
        arcpy.AddField_management(nodes_output_filepath, "NetTypeNo", "SHORT")

        # Compare and write changed inverts
        with arcpy.da.SearchCursor(db2_nodes, [id_field, invert_field, "Diameter", "nettypeno"]) as cursor2, \
                arcpy.da.InsertCursor(nodes_output_filepath, ["SHAPE@", "Change"]) as insert:

            for row in cursor2:
                node_id, inv2, diameter2, nettypeno2 = row
                if node_id in nodes_db1:
                    geom1, inv1, diameter1, nettypeno1 = nodes_db1[node_id]
                    text = ""
                    if inv1 != inv2:
                        text = "BK: {:.2f}→{:.2f}".format(inv1, inv2)

                    if diameter2 != diameter1:
                        if text:
                            text += "\n"
                        text += "D: ø{:.0f}→ø{:.0f}".format(diameter1 * 1e3, diameter2 * 1e3)
                    if len(text)>0:
                        insert.insertRow((geom1, text))


        arcpy.AddMessage("✅ Done. Node changes saved to: %s" % nodes_output_filepath)
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\CompareMIKEModels_Nodes.lyrx",
                 nodes_output_filepath, group=empty_group_layer, workspace_type="FILEGDB")

        # def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
        #              definition_query=None):




        return