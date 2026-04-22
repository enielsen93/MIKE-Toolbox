# Tool for reading DFS0 or KM2 files and creating LTS files from it
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import os
import sys
import numpy as np
import sqlite3
import arcpy

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


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "SQLITE Field Calculator"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [FieldCalculator]
    
    
class FieldCalculator(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Calculate Field"
        self.description = "Calculate Field"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        featureclass = arcpy.Parameter(
            displayName="Feature Class",
            name="featureclass",
            datatype="GPFeatureLayer",
            parameterType="Required",
            multiValue = "True",
            direction="Input")

        field1 = arcpy.Parameter(
            displayName="Field 1 to assign value to",
            name="field1",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        value1 = arcpy.Parameter(
            displayName="Value 1",
            name="value1",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        field2 = arcpy.Parameter(
            displayName="Field 2 to assign value to",
            name="field2",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        value2 = arcpy.Parameter(
            displayName="Value 2",
            name="value2",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        field3 = arcpy.Parameter(
            displayName="Field 3 to assign value to",
            name="field3",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        value3 = arcpy.Parameter(
            displayName="Value 3",
            name="value3",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        # Combine all parameters
        params = [featureclass, field1, value1, field2, value2, field3, value3]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if not parameters[0].Value:
            if arcgis_pro:
                # Reference the active map in the current project
                aprx = arcpymapping.ArcGISProject("CURRENT")
                map_view = aprx.activeMap

                if not parameters[0].value:
                    # List layers with selected features
                    layers = []
                    for layer in map_view.listLayers():
                        try:
                            if layer.getSelectionSet():
                                layers.append(layer.longName)
                        except:
                            pass
                    parameters[0].value = "; ".join(layers)
            else:
                mxd = arcpy.mapping.MapDocument("CURRENT")
                featureclasses = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if
                         lyr.getSelectionSet() and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
                if featureclasses:
                    parameters[0].value = "; ".join([featureclass for featureclass in featureclasses])

        if parameters[0].Value:
            fields = [f.name for f in arcpy.Describe(parameters[0].ValueAsText.split(";")[0]).fields]
            for i in [1, 3, 5]:
                if not parameters[i].Value:
                    parameters[i].filter.list = fields

        if parameters[1].Value:
            parameters[3].enabled = True
            parameters[4].enabled = True
            if parameters[3].Value:
                parameters[5].enabled = True
                parameters[6].enabled = True
            else:
                parameters[5].enabled = False
                parameters[6].enabled = False
        else:
            parameters[3].enabled = False
            parameters[4].enabled = False
            parameters[5].enabled = False
            parameters[6].enabled = False




        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        featureclasses = parameters[0].Values
        field1 = parameters[1].ValueAsText
        value1 = parameters[2].ValueAsText
        field2 = parameters[3].ValueAsText
        value2 = parameters[4].ValueAsText
        field3 = parameters[5].ValueAsText
        value3 = parameters[6].ValueAsText

        for featureclass in featureclasses:
            MU_database = os.path.dirname(arcpy.Describe(featureclass).catalogPath).replace("\mu_Geometry", "").replace("!delete!","")
            print(MU_database)
            featureclass_name = arcpy.Describe(featureclass).name
            arcpy.AddMessage(featureclass)

            arcpy.AddMessage("Confirm Query - Might be hidden behind window!")
            selection = [row[0] for row in arcpy.da.SearchCursor(featureclass, ["MUID"])]

            if arcgis_pro:
                import tkinter as tk
                from tkinter import messagebox

                def confirm_assignment(num_features):
                    root = tk.Tk()
                    root.withdraw()  # Hide the main window
                    result = messagebox.askyesno("Confirm Assignment", f"Assign value to {num_features} features?")
                    root.destroy()
                    return result

                userquery = confirm_assignment(len(selection))
                if userquery:
                    userquery = "Yes"
            else:
                userquery = pythonaddins.MessageBox(
                    "Assign value to %d features?" % (len(selection)),
                    "Confirm Assignment", 4)
                if len(selection)>100:
                    userquery = pythonaddins.MessageBox(
                        "Are you sure? Assign value to %d features?" % (len(selection)),
                        "Confirm Assignment", 4)

            if userquery == "Yes":
                arcpy.AddMessage(MU_database)
                if ".sqlite" in MU_database:
                    try:
                        connection = sqlite3.connect(
                            MU_database)
                        update_cursor = connection.cursor()
                        for field, value in zip([field1, field2, field3], [value1, value2, value3]):
                            if value:

                                arcpy.AddMessage("UPDATE %s SET %s = %s WHERE MUID IN %s" % (featureclass_name.replace("main.",""), field, value,
                                                                                             "('%s')" % ("','".join(selection))))
                                update_query = "UPDATE %s SET %s = %s WHERE MUID IN %s" % (featureclass_name.replace("main.",""), field, value,
                                                                                             "('%s')" % ("','".join(selection)))
                                update_cursor.execute(update_query)
                        connection.commit()
                        connection.close()
                    except Exception as e:
                        import traceback
                        arcpy.AddWarning(traceback.format_exc())
                        raise (e)

                    finally:
                        if connection:
                            connection.close()
                elif ".mdb" in MU_database:
                    for field, value in zip([field1, field2, field3], [value1, value2, value3]):
                        if value:
                            edit = arcpy.da.Editor(MU_database)
                            edit.startEditing(False, True)
                            edit.startOperation()
                            print(featureclass)
                            with arcpy.da.UpdateCursor(arcpy.Describe(featureclass).catalogPath, [field], where_clause = "MUID IN %s" % ("('%s')" % ("','".join(selection)))) as cursor:
                                for row in cursor:
                                    row[0] = value
                                    cursor.updateRow(row)

                            edit.stopOperation()
                            edit.stopEditing(True)
        return
        