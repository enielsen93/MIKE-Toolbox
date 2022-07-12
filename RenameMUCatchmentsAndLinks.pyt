# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy
import numpy as np
from arcpy import env
from shutil import copyfile
class Toolbox(object):
    def __init__(self):
        self.label =  "Rename Mike Urban features"
        self.alias  = "Rename Mike Urban features"

        # List of tool classes associated with this toolbox
        self.tools = [RenameMUFeatures] 

class RenameMUFeatures(object):
    def __init__(self):
        self.label       = "Rename Mike Urban Features"
        self.description = "Rename Mike Urban Features"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        database = arcpy.Parameter(
            displayName="Input Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
            
        outputDatabase = arcpy.Parameter(
            displayName="Output Mike Urban database",
            name="outputDatabase",
            datatype="File",
            parameterType="Required",
            direction="Output")
        outputDatabase.filter.list=["mdb"]
            
        renameLinks = arcpy.Parameter(
            displayName="Rename pipes?",
            name="renameLinksBool",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
            
        renameCatchments = arcpy.Parameter(
            displayName="Rename catchments?",
            name="renameCatchmentsBool",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
            
        parameters = [database, outputDatabase, renameLinks, renameCatchments]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        
        return

    def execute(self, parameters, messages):
        database = parameters[0].ValueAsText
        outputDatabase = parameters[1].ValueAsText
        renameLinks = parameters[2].Value
        renameCatchments = parameters[3].Value
        
        copyfile(database, outputDatabase)
        
        linksFile = outputDatabase + r"\mu_Geometry\msm_Link"
        catchmentsFile = outputDatabase + r"\mu_Geometry\ms_Catchment"
        catchConFile = outputDatabase + r"\msm_CatchCon"
        hModAFile = outputDatabase + r"\msm_HModA"
        links = {}
        
        if renameLinks:
            arcpy.SetProgressorLabel("Renaming pipes")
            arcpy.CalculateField_management(linksFile, "MUID", "!FROMNODE! + str(!OBJECTID!)", "PYTHON")
            try:
                with arcpy.da.SearchCursor(linksFile, ["MUID","FROMNODE","TONODE"]) as cursor:
                    for row in cursor:
                        links[row[0]] = row[1]
            except Exception as e:
                raise Exception('Error retrieving link connections. Try running "project check tool" or "recompute link features".\n%s' % (e.message))
            linksNewNames = {}

            existingNames = []
            for link in links:
                try:
                    linkID = 1
                    while links[link] + "l%1.0f" % linkID in existingNames:
                        linkID += 1
                    linksNewNames[link] = links[link] + "l%1.0f" % linkID
                    existingNames.append(links[link] + "l%1.0f" % linkID)
                except Exception as e:
                    arcpy.AddMessage("Error on Link %s" % link)
                    raise Exception('Error retrieving link connections. Try running "project check tool" or "recompute link features".\n%s' % (e.message))
            
            arcpy.SetProgressor("step","Renaming pipes", 0, int(arcpy.GetCount_management(linksFile)[0]), 1)       
            edit = arcpy.da.Editor(outputDatabase)
            edit.startEditing(False, False) 
            edit.startOperation()  
            with arcpy.da.UpdateCursor(linksFile, ["MUID"]) as cursor:
                for i,row in enumerate(cursor):
                    arcpy.SetProgressorPosition(i)
                    oldname = row[0]
                    if row[0] != linksNewNames[row[0]]:
                        row[0] = linksNewNames[row[0]]
                        try:
                            cursor.updateRow(row)
                        except Exception as e:
                            arcpy.AddError(str(e.message))
                            arcpy.AddError("Failed to update pipe %s to %s" % (oldname,linksNewNames[row[0]]))
                            raise arcpy.ExecuteError
            edit.stopOperation()  
            edit.stopEditing(True)  

        if renameCatchments:
            arcpy.SetProgressor("step","Renaming HModA", 0, int(arcpy.GetCount_management(hModAFile)[0]), 1)
            catchments = {}
            catchmentsCatchID = {}
            with arcpy.da.SearchCursor(catchConFile, ["MUID","CatchID","NodeID"]) as cursor:
                for row in cursor:
                    catchments[row[0]] = row[2]
                    catchmentsCatchID[row[1]] = row[0]
            
            catchmentsNewNames = {}
            existingNames = []
            for catchment in catchments:
                catchmentID = 1
                while catchments[catchment] + "_cc%1.0f" % catchmentID in existingNames or catchments[catchment] + "_cc%1.0f" % catchmentID in catchmentsNewNames:
                    #arcpy.AddMessage(catchmentsNewNames)
                    catchmentID += 1
                catchmentsNewNames[catchment] = catchments[catchment] + "_cc%1.0f" % catchmentID
                existingNames.append(catchments[catchment] + "_cc%1.0f" % catchmentID)
            
            
            with arcpy.da.UpdateCursor(hModAFile, ["CatchID"]) as cursor:
                for i,row in enumerate(cursor):
                    try:
                        arcpy.SetProgressorPosition(i)
                        row[0] = catchmentsNewNames[catchmentsCatchID[row[0]]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in HModA" % (row[0]))
                            raise arcpy.ExecuteError
                    except:
                        arcpy.AddError("Failed to update catchment %s in HModA - perhaps because catchment is unconnected" % (row[0]))
                        raise arcpy.ExecuteError
                    
                            
            arcpy.SetProgressor("step","Renaming catchments", 0, int(arcpy.GetCount_management(catchmentsFile)[0]), 1)
            with arcpy.da.UpdateCursor(catchmentsFile, ["MUID"]) as cursor:
                for i,row in enumerate(cursor):
                    arcpy.SetProgressorPosition(i)
                    # arcpy.AddMessage(catchmentsCatchID)
                    if row[0] != catchmentsNewNames[catchmentsCatchID[row[0]]]:
                        row[0] = catchmentsNewNames[catchmentsCatchID[row[0]]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in catchments" % (row[0]))
                            raise arcpy.ExecuteError
            
            arcpy.SetProgressor("step","Renaming catchment connections", 0, int(arcpy.GetCount_management(catchConFile)[0]), 1)
            with arcpy.da.UpdateCursor(catchConFile, ["MUID","CatchID"]) as cursor:
                for i,row in enumerate(cursor):
                    arcpy.SetProgressorPosition(i)
                    if row[1] != catchmentsNewNames[row[0]]: #  or row[0] != catchmentsNewNames[row[0]]
                        row[1] = catchmentsNewNames[row[0]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in catchment connection file" % (row[0]))
                            raise arcpy.ExecuteError
            arcpy.SetProgressor("step","Renaming HModA", 0, int(arcpy.GetCount_management(hModAFile)[0]), 1)
            catchments = {}
            catchmentsCatchID = {}
            with arcpy.da.SearchCursor(catchConFile, ["MUID","CatchID","NodeID"]) as cursor:
                for row in cursor:
                    catchments[row[0]] = row[2]
                    catchmentsCatchID[row[1]] = row[0]
            
            catchmentsNewNames = {}
            existingNames = []
            for catchment in catchments:
                catchmentID = 1
                while catchments[catchment] + "_c%1.0f" % catchmentID in existingNames or catchments[catchment] + "_c%1.0f" % catchmentID in catchmentsNewNames:
                    catchmentID += 1
                catchmentsNewNames[catchment] = catchments[catchment] + "_c%1.0f" % catchmentID
                existingNames.append(catchments[catchment] + "_c%1.0f" % catchmentID)
            
            
            with arcpy.da.UpdateCursor(hModAFile, ["CatchID"]) as cursor:
                for i,row in enumerate(cursor):
                    try:
                        arcpy.SetProgressorPosition(i)
                        row[0] = catchmentsNewNames[catchmentsCatchID[row[0]]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in HModA" % (row[0]))
                            raise arcpy.ExecuteError
                    except:
                        arcpy.AddError("Failed to update catchment %s in HModA - perhaps because catchment is unconnected" % (row[0]))
                        raise arcpy.ExecuteError
                    
                            
            arcpy.SetProgressor("step","Renaming catchments", 0, int(arcpy.GetCount_management(catchmentsFile)[0]), 1)
            with arcpy.da.UpdateCursor(catchmentsFile, ["MUID"]) as cursor:
                for i,row in enumerate(cursor):
                    arcpy.SetProgressorPosition(i)
                    if row[0] != catchmentsNewNames[catchmentsCatchID[row[0]]]:
                        row[0] = catchmentsNewNames[catchmentsCatchID[row[0]]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in catchments" % (row[0]))
                            raise arcpy.ExecuteError
            
            arcpy.SetProgressor("step","Renaming catchment connections", 0, int(arcpy.GetCount_management(catchConFile)[0]), 1)
            with arcpy.da.UpdateCursor(catchConFile, ["MUID","CatchID"]) as cursor:
                for i,row in enumerate(cursor):
                    arcpy.SetProgressorPosition(i)
                    if row[1] != catchmentsNewNames[row[0]]: #  or row[0] != catchmentsNewNames[row[0]]
                        row[1] = catchmentsNewNames[row[0]]
                        try:
                            cursor.updateRow(row)
                        except:
                            arcpy.AddError("Failed to update catchment %s in catchment connection file" % (row[0]))
                            raise arcpy.ExecuteError         
        return