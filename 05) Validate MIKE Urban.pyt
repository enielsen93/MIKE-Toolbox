# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:21:31 2018

@author: eni
"""
import os
import arcpy
import numpy as np
import re
# import pythonaddins
import hashlib
import math
import sys
# functionsPath = [r"K:\Hydrauliske modeller\Makroer & Beregningsark\Functions", r"C:\Dokumenter\Makroer\Functions", os.path.join(os.path.dirname(os.path.dirname(__file__)),"Functions")]
# i = 0
# while not os.path.exists(functionsPath[i]):
#     i += 1
# sys.path.append(functionsPath[i])
# import TranslateMDBToSQLLite
thisFolder = os.path.dirname(__file__)
# scriptFolder = os.path.join(thisFolder, r"Scripts")
# sys.path.append(scriptFolder)
# import networkx as nx

def translate(field, table="", dbExt="sqlite"):
    dictionairy = {"ms_Catchment": "msm_Catchment",
                   "msm_HModA": "msm_Catchment",
                   "FromNode": "fromnodeid",
                   "ToNode": "tonodeid"}
    if dbExt == "sqlite":
        return dictionairy[field]
    else:
        return field


def getAvailableFilename(filepath):
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
        self.label =  "Validate Mike Urban Database"
        self.alias  = "Validate Mike Urban Database"

        # List of tool classes associated with this toolbox
        self.tools = [CheckMikeUrbanDatabase, CleanupMUS]

class CheckMikeUrbanDatabase(object):
    def __init__(self):
        self.label       = "1) Validate Mike Urban Database"
        self.description = "1) Validate Mike Urban Database"

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
            
        set_basin_cross_section_area = arcpy.Parameter(
            displayName="Set cross section area to basins (assuming 1:3 length width relationship)",
            name="set_basin_cross_section_area",
            datatype="Boolean",
            parameterType="optional",
            direction="Output")
        set_basin_cross_section_area.value = False
        set_basin_cross_section_area.category = "Additional settings"
        
        # check_basin_outlet_elevation = arcpy.Parameter(
            # displayName="Check basin weir level and critical level",
            # name="check_basin_outlet_elevation",
            # datatype="Boolean",
            # parameterType="optional",
            # direction="Output")
        # check_basin_outlet_elevation.value = False
        # check_basin_outlet_elevation.category = "Additional settings"
        
        # check_basin_upstream_critical_level = arcpy.Parameter(
            # displayName="Check if any o",
            # name="check_basin_outlet_elevation",
            # datatype="Boolean",
            # parameterType="optional",
            # direction="Output")
        # check_basin_outlet_elevation.value = False
        # check_basin_outlet_elevation.category = "Additional settings"
        
        run_cleanup = arcpy.Parameter(
            displayName="Run clean-up: Delete all external boundaries that belong to nodes that do not exist",
            name="run_cleanup",
            datatype="Boolean",
            parameterType="optional",
            direction="Output")
        run_cleanup.value = False
        run_cleanup.category = "Additional settings"
        
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

        parameters = [MU_database, set_basin_cross_section_area, run_cleanup]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        if not parameters[0].valueAsText:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            database = [lyr.dataSource for lyr in arcpy.mapping.ListLayers(mxd) if not lyr.getSelectionSet() and lyr.isFeatureLayer and ".mdb" in lyr.dataSource]
            if database:
                database = re.findall(r"(.+)(?=\.mdb)", database[0])[0] + ".mdb"
                parameters[0].value = database
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        set_basin_cross_section_area = parameters[1].Value
        run_cleanup = parameters[2].Value
        # check_basin_outlet_elevation = parameters[2].Value
        dbExtension = MU_database.split(".")[-1]
        catchments = os.path.join(MU_database, translate("ms_Catchment", dbExt = dbExtension))
        catchCon = MU_database + "\msm_CatchCon"
        hModA = os.path.join(MU_database, translate("msm_HModA", dbExt = dbExtension))
        hParA = MU_database + "\msm_HParA"
        msm_Node = MU_database + "\msm_Node"
        msm_Link = MU_database + "\msm_Link"
        msm_Weir = MU_database + "\msm_Weir"
        msm_Orifice = MU_database + "\msm_Orifice"
        msm_PasReg = MU_database + "\msm_PasReg"
        msm_Pump = MU_database + "\msm_Pump"
        networkLinks = [msm_Link, msm_Weir, msm_Orifice]
        msm_BItem = MU_database + "\msm_BItem"
        msm_BBoundary = MU_database + "\msm_BBoundary"       

        is_sqlite = True if ".sqlite" in MU_database else False

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

        try:
            catchmentsSinglePart = arcpy.MultipartToSinglepart_management(catchments,
                                                   r"in_memory\singlepartCatchments")

            MUIDsSinglePart = []
            with arcpy.da.SearchCursor(catchmentsSinglePart,["MUID"]) as cursor:
                for row in cursor:
                    MUIDsSinglePart.append(row[0])

            singlePartDuplicates = set([x for x in MUIDsSinglePart if MUIDsSinglePart.count(x) > 1])
            multiPartFeatures = [a for a in singlePartDuplicates if not a in duplicates]
            if len(multiPartFeatures) > 0:
                arcpy.AddWarning(
                    "Warning: Multipart catchments in table ms_Catchment found (unable to import catchments into Dandas): '" + "', '".join(
                        multiPartFeatures) + "'")
        except Exception as e:
            arcpy.AddWarning("Could not check for multipart features")

        if dbExtension == "mdb":
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

        nodes = [row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID")]

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

        if dbExtension == "mdb":
            rows = int(np.sum([1 for i in arcpy.da.SearchCursor(hModA,"OBJECTID")]))
            columns = ["CatchID","ParAID","LocalNo","RFactorHModA", "RFactorHParA"]
            catchmentsReductionFactor = np.empty((rows, len(columns)), dtype='O')
            with arcpy.da.SearchCursor(hModA,["CatchID","ParAID","LocalNo","RFactor"]) as cursor:
                for i,row in enumerate(cursor):
                    catchmentsReductionFactor[i,0:4] = row[0:4]
            with arcpy.da.SearchCursor(hParA,["MUID","RedFactor"]) as cursor:
                for i,row in enumerate(cursor):
                    catchmentsReductionFactor[catchmentsReductionFactor[:,1] == row[0], 4] = row[1]

            nonlocalParIdx = np.where(catchmentsReductionFactor[:,2] == 0)[0]
            for rFactorHModA in set(catchmentsReductionFactor[nonlocalParIdx, 4]):
                try:    
                    arcpy.AddMessage("%d/%d catchments with reduction factor set to %1.2f" % (len(np.intersect1d(nonlocalParIdx, np.where(catchmentsReductionFactor[:,4] == rFactorHModA)[0])),rows,rFactorHModA))
                except Exception as e:
                    pass
                    
            localParIdx = np.where(catchmentsReductionFactor[:,2] == 1)[0]
            for rFactorHParA in set(catchmentsReductionFactor[localParIdx, 3]):
                arcpy.AddMessage("%d/%d catchments with reduction factor set to %1.2f (Local Parameters)" % (len(np.intersect1d(localParIdx,np.where(catchmentsReductionFactor[:,3] == rFactorHParA)[0])),rows,rFactorHParA))
                
        catchmentsDrainageAreaSame = []
        catchmentsDrainageAreaDiffer = []
        with arcpy.da.SearchCursor(catchments, ["MUID", "Area", "SHAPE@AREA"], where_clause = 'Area IS NOT NULL') as cursor:
            for row in cursor:
                shapeArea = row[2]/1e4 if not is_sqlite else row[2]
                if min(row[1],shapeArea)>0 and abs(abs(shapeArea)-abs(row[1]))/min(row[1],shapeArea)*1e2>10:
                    catchmentsDrainageAreaDiffer.append(row[0])
                else:
                    catchmentsDrainageAreaSame.append(row[0])
                    
        if len(catchmentsDrainageAreaSame)>0:
            arcpy.AddMessage("%d catchments with drainage area set (approximately same as geometric area of catchment): '%s'" % (len(catchmentsDrainageAreaSame),"', '".join(catchmentsDrainageAreaSame)))

        if len(catchmentsDrainageAreaSame)>0:
            arcpy.AddMessage("%d catchments with drainage area set (different from geometric area of catchment): '%s'" % (len(catchmentsDrainageAreaDiffer),"', '".join(catchmentsDrainageAreaDiffer)))
        
        arcpy.AddMessage("Catchment check done.")
        
        # nodesOutletShape = [row[0] for row in arcpy.da.SearchCursor(msm_Node,["MUID","LossParID"]) if "Classic" in row[1]]
        # if len(nodesOutletShape)>0:
            # arcpy.AddWarning("Warning: Manholes with outlet head loss ID set to MOUSE Classic(Engelund) (should probably be Weighted Inlet Energy): '" + "', '".join(
                # nodesOutletShape
                # ) + "'")

        node_muids = [row[0] for row in arcpy.da.SearchCursor(msm_Node, "MUID")]
        duplicate_nodes = [node for node in node_muids if node_muids.count(node) > 1]

        if duplicate_nodes:
            arcpy.AddWarning("Warning: Duplicate MUIDs in manholes: ('" + "', '".join(duplicate_nodes) + "')")

        link_muids = [row[0] for row in arcpy.da.SearchCursor(msm_Link, "MUID")]
        duplicate_links = [link for link in link_muids if link_muids.count(link) > 1]

        if duplicate_links:
            arcpy.AddWarning("Warning: Duplicate MUIDs in links: ('" + "', '".join(duplicate_links) + "')")

        if not is_sqlite:
            regulated_links = [row[0] for row in arcpy.da.SearchCursor(msm_PasReg, ["LinkID"])]
            broken_regulations = [link_id for link_id in regulated_links if link_id not in link_muids]
            if broken_regulations:
                arcpy.AddWarning("Warning: Regulations that apply to missing links: ('" + "', '".join(broken_regulations) + "')")

        nodesBadName = [row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID",where_clause = "MUID LIKE 'Node_*'")]
        if len(nodesBadName)>0:
            arcpy.AddWarning("Warning: Manholes with unoriginal names found (not suitable for FLORA): '" + "', '".join(nodesBadName) + "'")
            
        linksBadName = [row[0] for row in arcpy.da.SearchCursor(msm_Link, "MUID", where_clause = "MUID LIKE 'Link_*'")]
        if len(linksBadName)>0:
            arcpy.AddWarning("Warning: Pipes with unoriginal names found (not suitable for FLORA): '" + "', '".join(linksBadName) + "'")

        links_different_length = []
        with arcpy.da.SearchCursor(msm_Link, ["MUID", "Length", "SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                if row[1] and row[2] and row[2]>10 and row[1] != 10 and abs(row[2]-row[1])/row[2]>0.1:
                    links_different_length.append(row[0])
        if links_different_length:
            arcpy.AddWarning("Pipes with length set that's different from geometric length: ('" + "', '".join(links_different_length) + "')")

        linksDiameterPlastic = [row[0] for row in arcpy.da.SearchCursor(msm_Link,"Diameter",where_clause = "MaterialID = 'Plastic'")]
        outerDiametersFound = set()
        for d in linksDiameterPlastic:
            if d in [0.16, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]:
                outerDiametersFound.add(str(d))
        if len(outerDiametersFound)>0:
            arcpy.AddWarning("Plastic Pipes found with diameter that might be the outer diameter instead of inner (diameters: %s)" %  (", ".join(outerDiametersFound)))

        if not is_sqlite:
            with arcpy.da.SearchCursor(msm_Node, ["MUID", "GeometryID","InvertLevel"],where_clause = "TypeNo = 2") as cursor:
                for row in cursor:
                    try:
                        geometryInvertLevel = np.min([a[0] for a in arcpy.da.SearchCursor(MU_database + r"\ms_TabD",["Value1"],where_clause = "TabID = '%s'" % (row[1]))])
                    except Exception as e:
                        arcpy.AddError(row)
                        raise(e)
                    if geometryInvertLevel != row[2]:
                        arcpy.AddWarning("Invert Level of basin node %s is different from invert level of basin geometry %s (%1.2f - %1.2f)" % (row[0], row[1], row[2], geometryInvertLevel))
        
        if set_basin_cross_section_area:
            basinGeometries = [row[0] for row in arcpy.da.SearchCursor(msm_Node, ["GeometryID"],where_clause = "TypeNo = 2")]
            for geometry in basinGeometries:
                ms_TabDArray = np.zeros((0,4))
                with arcpy.da.UpdateCursor(MU_database + r"\ms_TabD", ["Sqn", "Value1", "Value2", "Value3"], where_clause = "TabID = '%s'" % (geometry)) as cursor:
                     for row in cursor:
                         ms_TabDArray = np.append(ms_TabDArray, [row], axis=0)
                     idx = np.argsort(ms_TabDArray[:,1])
                     ms_TabDArray = ms_TabDArray[idx,:]
                     for i in np.arange(1,ms_TabDArray.shape[0]):
                         ms_TabDArray[i,2] = max(1,np.round((ms_TabDArray[i,1]-ms_TabDArray[0,1])*math.sqrt(ms_TabDArray[i,3]/3)))
                     cursor.reset()
                     for row in cursor:
                         if row[2]>0:
                             row[2] = ms_TabDArray[ms_TabDArray[:,0]==row[0],2][0]
                             cursor.updateRow(row)
                     arcpy.AddMessage("Updated basin geometry %s with cross section area" % (geometry))
                     
        # if check_basin_outlet_elevation:
        
        msm_Node_diameter_lossParID = {}
        with arcpy.da.SearchCursor(msm_Node, ["MUID","Diameter","LossParID"], where_clause = "TypeNo NOT IN (2,3)") as cursor:
            for row in cursor:
                msm_Node_diameter_lossParID[row[0]] = [row[1], row[2]]

        MUIDsCritLevel = {}
        MUIDsInvertLevel = {}
        with arcpy.da.SearchCursor(msm_Node,["MUID","GroundLevel","CriticalLevel","InvertLevel"]) as cursor:
            for row in cursor:
                if False:#row[2]:
                    MUIDsCritLevel[row[0]] = row[2]
                else:
                    MUIDsCritLevel[row[0]] = row[1]        
                MUIDsInvertLevel[row[0]] = row[3]
                
        msm_Node_saddle = []
        msm_Node_PipeAboveGround = []
        try:
            with arcpy.da.SearchCursor(msm_Link, [translate("FromNode", dbExt = dbExtension), translate("ToNode", dbExt = dbExtension), "Diameter"]) as cursor:
                for row in cursor:
                    for node in [row[0],row[1]]:
                        if (node not in msm_Node_saddle and
                            node in msm_Node_diameter_lossParID and 
                            row[2] > msm_Node_diameter_lossParID[node][0] 
                            and ("Weighted" in msm_Node_diameter_lossParID[node][1] or 
                                 "Classic" in msm_Node_diameter_lossParID[node][1])):
                            msm_Node_saddle.append(node)
                    
        except Exception as e:
            arcpy.AddError("Error: Try running Project Check Tool.")
            raise(e)
            
        if msm_Node_saddle:
            arcpy.AddWarning("Warning: Manholes found with diameter less than diameter of one more connected links with an outlet head loss defined. Consider removing outlet head loss: ('%s')" % ("', '".join(msm_Node_saddle)))
            
        try: 
            with arcpy.da.SearchCursor(msm_Link, [translate("FromNode", dbExt = dbExtension), translate("ToNode", dbExt = dbExtension), "Diameter", "MUID", "UpLevel", "DwLevel", "UpLevel_C", "DwLevel_C","CrsID"]) as cursor:
                    for row in cursor:
                        if (MUIDsCritLevel[row[0]] < max(row[4], row[6])+row[2] or MUIDsCritLevel[row[1]] < max(row[5], row[7])+row[2]) and not row[8]:
                            msm_Node_PipeAboveGround.append(row[3])
        except Exception as e:
            arcpy.AddError("Error on link %s: Try running Project Check Tool and Recompute on msm_Link." % (row[3]))
            arcpy.AddError(e.message)
        if len(msm_Node_PipeAboveGround)>0: 
            arcpy.AddWarning("Warning: Pipes found with top of inlet or outlet above critical level of connected nodes (Invert Level + Pipe Diameter > Ground Level): ('%s')" % ("', '".join(msm_Node_PipeAboveGround)))
                       
        # network = nx.DiGraph()
        # network.add_nodes_from([row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID")])
        # for networkLink in networkLinks:
            # with arcpy.da.SearchCursor(networkLink,["FromNode","ToNode"]) as cursor:
                # for row in cursor:
                    # network.add_edge(row[0],row[1])
                
        basins = [row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID", where_clause = "TypeNo = 2")]
            
        # for basin in basins:
            # nodesUpstream = nx.ancestors(network,basin)
            # for node, critlevel in MUIDsCritLevel.items():
                # if node in nodesUpstream and critlevel < MUIDsCritLevel[basin]:
                    # arcpy.AddWarning("Warning: Critical level of Manhole %s is %1.2f which is below critical level of downstream basin %s with critical level %1.2f. Consider changing cover type to Sealed." % (node, critlevel, basin, MUIDsCritLevel[basin]))
            
        # if dbExtension == "mdb":
        #     boundaryItemsWithScalingFactor = [row[0] for row in arcpy.da.SearchCursor(msm_BItem,["BoundaryID","Fraction"]) if row[1] != 1]
        #     if len(boundaryItemsWithScalingFactor)>0:
        #         arcpy.AddWarning("Boundary Item %s has scaling factor different from 1. Safety and hydrological reduction factors should be Time-Area parameters." % ( " and ".join(boundaryItemsWithScalingFactor)))
        
        if run_cleanup:
            nodes = [row[0] for row in arcpy.da.SearchCursor(msm_Node,["MUID"])]
            boundariesNodeID = [row[0] for row in arcpy.da.SearchCursor(msm_BBoundary,["NodeID"], where_clause = "TypeNo = 12") if row[0] not in nodes]
            arcpy.AddMessage("Deleting Water Level boundaries '%s' as they are located at non-existing manholes" % ("', '".join(boundariesNodeID)))
            with arcpy.da.UpdateCursor(msm_BBoundary,["NodeID"], where_clause = "NodeID IN ('%s')" % ("', '".join(boundariesNodeID))) as cursor:
                for row in cursor:
                    cursor.deleteRow()
                    
        
        return


class CleanupMUS(object):
    def __init__(self):
        self.label = "a) Cleanup missing elements from MUS files"
        self.description = "a) Cleanup missing elements from MUS files"

    def getParameterInfo(self):
        # Define parameter definitions

        reference_feature_class = arcpy.Parameter(
            displayName="Reference Feature Class",
            name="reference_feature_class",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        selection_files = arcpy.Parameter(
            displayName="Selection Files",
            name="selection_files",
            datatype="file",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        selection_files.filter.list = ["mus"]

        parameters = [reference_feature_class, selection_files]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):  # optional
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        reference_feature_class = parameters[0].Value
        selection_files = parameters[1].ValueAsText.split(";")

        to_ignore = ["msm_Node", "msm_Link", "ms_Catchment", "msm_Catchment", "msm_Orifice", "msm_Weir"]

        MUIDs = [row[0] for row in arcpy.da.SearchCursor(reference_feature_class, ["MUID"])]

        for selection_file in selection_files:
            with open(selection_file, 'r') as f:
                selection_file_lines = f.readlines()
                selection_file_lines_cleaned = []
                for line in selection_file_lines:
                    line_MUID = line[:-1]
                    if len(line_MUID) > 1 and not line_MUID in to_ignore and not line_MUID in MUIDs:
                        arcpy.AddMessage((selection_file, line_MUID))
                    else:
                        selection_file_lines_cleaned.append(line)
            with open(selection_file, 'w') as f:
                for line in selection_file_lines:
                    f.write(line)
        return
