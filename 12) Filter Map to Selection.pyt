import os
import numpy as np
import re

if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
    from arcpy._mapping import Layer
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument

class Toolbox(object):
    def __init__(self):
        self.label = "Set Definition Query to Selection"
        self.alias = "Set Definition Query to Selection"
        self.canRunInBackground = True
        # List of tool classes associated with this toolbox
        self.tools = [SetDefinitionQuery]


class SetDefinitionQuery(object):
    def __init__(self):
        self.label = "Set Definition Query to Selection"
        self.description = "Set Definition Query to Selection"
        self.canRunInBackground = False



    def getParameterInfo(self):
        # Define parameter definitions

        layer = arcpy.Parameter(
            displayName="Layer",
            name="layer",
            datatype="GPFeatureLayer",
            parameterType="Required",
            multiValue = True,
            direction="Input")

        append = arcpy.Parameter(
            displayName="Append to existing Definition Query",
            name="append",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        remove_selection = arcpy.Parameter(
            displayName="Remove from filter instead",
            name="remove_selection",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        specific_definition_query = arcpy.Parameter(
            displayName="Specific Definition Query",
            name="specific_definition_query",
            datatype="GPString",
            parameterType="Optional",
            category="Additional Settings",
            direction="Input")

        apply_to_labels = arcpy.Parameter(
            displayName="Apply to labels instead of definition query",
            name="apply_to_labels",
            datatype="GPBoolean",
            parameterType="Optional",
            category="Additional Settings",
            direction="Input")

        parameters = [layer, append, remove_selection, specific_definition_query, apply_to_labels]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if not parameters[0].Values:
            if arcgis_pro:
                # Reference the active map in the current project
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap

                layers = []
                # List layers with selected features
                for layer in map_view.listLayers():
                    try:
                        if layer.getSelectionSet():
                            layers.append(layer.longName)
                    except:
                        pass
            else:
                mxd = arcpy.mapping.MapDocument("CURRENT")

                layers = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                         lyr.getSelectionSet()]
            if layers:
                parameters[0].value = layers


        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        layers = parameters[0].Values
        append = parameters[1].Value
        remove_selection = parameters[2].Value
        specific_definition_query = parameters[3].ValueAsText
        apply_to_labels = parameters[4].Value

        def setDefinitionQuery(layer, definition_query):
            if arcgis_pro:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap

                layers = m.listLayers()
                layer = [lyr for lyr in layers if lyr.name == layer.name and lyr.dataSource == layer.dataSource][0]

                layer.updateDefinitionQueries([{"name": "Query 1", "sql": definition_query, "isActive": True}])
                aprx.save()
            else:
                layer.definitionQuery = definition_query

        def setLabelQuery(layer, definition_query, append = True):
            if arcgis_pro:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap

                layers = m.listLayers()
                layer = [lyr for lyr in layers if lyr.name == layer.name and lyr.dataSource == layer.dataSource][0]

                for lbl_class in layer.listLabelClasses():
                    existing_query = lbl_class.SQLQuery.strip()

                    if existing_query and append:
                        lbl_class.SQLQuery = "{} OR {}".format(existing_query, definition_query)
                    else:
                        lbl_class.SQLQuery = definition_query
                    arcpy.AddMessage(definition_query)

                    lbl_class.showClassLabels = True  # Make sure labeling is enabled
            else:
                raise RuntimeError("This script is currently not supported in ArcMap. Use ArcGIS Pro.")

        if apply_to_labels:
            for layer in layers:
                if specific_definition_query:
                    setLabelQuery(layer, specific_definition_query, append = append)
                else:
                    oid_fieldname = arcpy.Describe(layer).OIDFieldName

                    if "muid" in [field.name.lower() for field in arcpy.ListFields(layer)] and (
                            hasattr(layer, "datasetName") and "CatchConLink" not in layer.datasetName):
                        arcpy.AddMessage((layer.getSelectionSet()))
                        new_definition_query = "muid %sIN ('%s')" % ("NOT " if remove_selection else "", "', '".join(
                            [row[0] for row in arcpy.da.SearchCursor(layer, ["muid"], where_clause="%s IN (%s)" % (
                            oid_fieldname, ", ".join([str(l) for l in layer.getSelectionSet()])))]))
                    else:
                        new_definition_query = "%s %sIN (%s)" % (oid_fieldname, "NOT " if remove_selection else "",
                                                                 ", ".join([str(g) for g in layer.getSelectionSet()]))

                    setLabelQuery(layer, new_definition_query, append = append)
        else:
            for layer in layers:
                if specific_definition_query: # specific definition query
                    old_definition_query = layer.definitionQuery
                    if old_definition_query and append:
                        old_definition_query = layer.definitionQuery
                        if arcgis_pro:
                            setDefinitionQuery(layer, old_definition_query + " AND " + specific_definition_query)
                        else:
                            layer.definitionQuery = old_definition_query + " AND " + specific_definition_query
                    else:
                        if arcgis_pro:
                            setDefinitionQuery(layer, specific_definition_query)
                        else:
                            layer.definitionQuery = specific_definition_query
                    # layer.updateDefinitionQueries([specific_definition_query])

                else:
                    oid_fieldname = arcpy.Describe(layer).OIDFieldName
                    # arcpy.AddMessage([row for row in arcpy.da.SearchCursor(layer, ["muid"], where_clause = "objectid IN (%s)" % (", ".join([str(l) for l in layer.getSelectionSet()])))])
                    # arcpy.AddMessage("objectid IN (%s)" % (", ".join([str(l) for l in layer.getSelectionSet()])))
                    if "muid" in [field.name.lower() for field in arcpy.ListFields(layer)] and (hasattr(layer, "datasetName") and "CatchConLink" not in layer.datasetName):
                        arcpy.AddMessage((layer.getSelectionSet()))
                        new_definition_query = "muid %sIN ('%s')" % ("NOT " if remove_selection else "", "', '".join([row[0] for row in arcpy.da.SearchCursor(layer, ["muid"], where_clause = "%s IN (%s)" % (oid_fieldname, ", ".join([str(l) for l in layer.getSelectionSet()])))]))
                    else:
                        new_definition_query = "%s %sIN (%s)" % (oid_fieldname, "NOT " if remove_selection else "", ", ".join([str(g) for g in layer.getSelectionSet()]))


                    old_definition_query = layer.definitionQuery
                    if old_definition_query and append:
                        if arcgis_pro:
                            setDefinitionQuery(layer, old_definition_query + " AND " + new_definition_query)
                        else:
                            layer.definitionQuery = old_definition_query + " AND " + new_definition_query

                    else:
                        if arcgis_pro:
                            setDefinitionQuery(layer, new_definition_query)
                        else:
                            layer.definitionQuery = new_definition_query

        return
