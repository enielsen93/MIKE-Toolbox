# -*- coding: utf-8 -*-
import os
import arcpy
import numpy as np
import re
from arcpy import env
arcpy.env.addOutputsToMap = False
import codecs
import sys
import mousereader
import ColebrookWhite
import networker
import codecs

from subprocess import call
from shutil import copyfile
import traceback
import scipy

def m11extrapath():
    m11extraPath = r"C:\Program Files (x86)\DHI\2016\bin\m11extra.exe"
    i = 2030
    while not os.path.exists(m11extraPath.replace("2016",str(i))) and i < 20:
        i -= 1
    m11extraPath = m11extraPath.replace("2016",str(i))
    return m11extraPath

if "mapping" in dir(arcpy):
    arcgis_pro = False
    import arcpy.mapping as arcpymapping
    from arcpy.mapping import MapDocument as arcpyMapDocument
else:
    arcgis_pro = True
    import arcpy.mp as arcpymapping
    from arcpy.mp import ArcGISProject as arcpyMapDocument
# arcpy.env.workspace = arcpy.env.scratchGDB

def getAvailableFilename(filepath, parent = None):
    parent = "F%s" % (parent) if parent and parent[0].isdigit() else None
    parent = os.path.basename(re.sub(r"\.[^\.\\]+$","", parent)).replace(".","_").replace("-","_").replace(" ","_").replace(",","_") if parent else None
    filepath = "%s\%s_%s" % (os.path.dirname(filepath), parent, os.path.basename(filepath)) if parent else filepath
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath
 
class Toolbox(object):
    def __init__(self):
        self.label =  "Display Mike Urban Results"
        self.alias  = "Display Mike Urban Results"

        # List of tool classes associated with this toolbox
        self.tools = [DisplayFloodReturnPeriodFun, DisplayWeirStatistics, DisplayFlowStatistics, DisplayQFullQMax, DisplayWeirReturnPeriod, DisplayMIKE1DResults, DisplayExtent]

class DisplayFloodReturnPeriodFun(object):
    def __init__(self):
        self.label       = "Display Flood Return Period"
        self.description = "Display Flood Return Period"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        critical_return_period = arcpy.Parameter(
            displayName="Critical Return Period (5 years, 10 years or 20 years)",
            name="critical_return_period",
            datatype="Double",
            parameterType="Optional",
            direction="Input")

        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        # exportShape = arcpy.Parameter(
        #     displayName="Output manholes with Flood Return Period",
        #     name="exportShape",
        #     datatype="DEFeatureClass",
        #     parameterType="Required",
        #     direction="Output")
        #
        # exportBasins = arcpy.Parameter(
        #     displayName="Output basins with Flood Return Period",
        #     name="exportBasins",
        #     datatype="DEFeatureClass",
        #     parameterType="Required",
        #     direction="Output")
            
        flowFile = arcpy.Parameter(
            displayName="Include PRF file in order to show maximum discharge from basins and permanent water volume",
            name="flowFile",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        flowFile.filter.list=["prf"]
        
        traceNetwork = arcpy.Parameter(
            displayName="Calculate connected catchment area for each basin",
            name="traceNetwork",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        traceNetwork.category = "Get Catchment Area"
        traceNetwork.enabled = False
        
        reaches = arcpy.Parameter(
            displayName="Trace network through:",
            name="reaches",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        reaches.filter.type = "ValueList"  
        reaches.filter.list = ["Orifice","Weir","Pump", "Basin"]
        reaches.value = ["Orifice","Weir","Pump", "Basin"]
        reaches.category = "Get Catchment Area"
        reaches.enabled = False
        
        breakChainOnNodes = arcpy.Parameter(
            displayName="End each trace at following node MUIDs (each node should by delimited by a comma: node1, node2)",
            name="breakChainOnNodes",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        breakChainOnNodes.category = "Additional settings"
        breakChainOnNodes.category = "Get Catchment Area"
        breakChainOnNodes.enabled = False

        use_critical_level = arcpy.Parameter(
            displayName="Use Critical Level from nodes",
            name="use_critical_level",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        use_critical_level.enabled = True

        parameters = [erfFile, observationPeriod, critical_return_period, mike_urban_database, flowFile, traceNetwork, reaches, breakChainOnNodes, use_critical_level]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        erfFile = parameters[0].ValueAsText
        observationPeriod = float(parameters[1].ValueAsText)
        critical_return_period = float(parameters[2].ValueAsText)
        mike_urban_database = parameters[3].ValueAsText
        msm_Node = mike_urban_database + "\msm_Node"
        msm_Link = mike_urban_database + "\msm_Link"
        msm_Weir = mike_urban_database + "\msm_Weir"
        msm_Orifice = mike_urban_database + "\msm_Orifice"
        msm_Pump = os.path.join(mike_urban_database,"msm_Pump")
        msm_CatchCon = os.path.join(mike_urban_database,"msm_CatchCon")
        ms_Catchment = os.path.join(mike_urban_database,"ms_Catchment")
        msm_HModA = os.path.join(mike_urban_database,"msm_HModA")
        msm_HParA = os.path.join(mike_urban_database,"msm_HParA")

        flowFile = parameters[4].ValueAsText
        traceNetwork = parameters[5].ValueAsText
        reaches = parameters[6].ValueAsText
        # break_chain_on_nodes = parameters[7].ValueAsText
        use_critical_level = parameters[8].Value

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)
        MIKE_gdb = os.path.join(MIKE_folder, os.path.splitext(os.path.basename(mike_urban_database))[0])
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
                MIKE_gdb = os.path.join(MIKE_folder, "%s_%d" % (os.path.splitext(os.path.basename(mike_urban_database))[0], dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        arcpy.env.addOutputsToMap = False

        
        arcpy.SetProgressorLabel("Getting critical level for manholes")
        MUIDsHCrit = {}
        with arcpy.da.SearchCursor(msm_Node,["MUID","GroundLevel","CriticalLevel"]) as cursor:
            for row in cursor:
                if use_critical_level and row[2]:
                    MUIDsHCrit[row[0]] = row[2]
                else:
                    MUIDsHCrit[row[0]] = row[1]

        MUIDs = MUIDsHCrit.keys()
        if flowFile:
            arcpy.SetProgressorLabel("Getting Discharge from Basins")
            workingDir = os.path.dirname(__file__)
            M11OUT = workingDir + "\M11.OUT"
            M11IN = workingDir + "\M11.IN"

            prfFileCopy = workingDir + "\prffile.PRF"
            resultFile = workingDir + "\ResultFile.txt"
            
            copyfile(flowFile,prfFileCopy)
            call([m11extrapath(), prfFileCopy])

            lines = ""
            links = []
            with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
                for linei,line in enumerate(M11OUTFile):
                    try:
                        if ("Link_Q" in line or "Weir_Q" in line or "Orifice_Q" in line) and re.findall("<([^>]+)>",line)[0] not in links:		
                            links.append(re.findall("<([^>]+)>",line)[0])
                            line = re.sub("^0","1",line)
                        lines += line
                    except UnicodeDecodeError:
                        lines += line
                        arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                    except Exception as e:
                        raise(e)
                           
            with codecs.open(M11IN,'w','cp1252') as M11INFile:
                M11INFile.write(lines)

            call([m11extrapath(), prfFileCopy, resultFile])

            linksFlow = []
            with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
                txt = M11OUTFile.readlines()

            linksFlow = {}
            for link in links:
                linksFlow[link] = []
            lines = range([i for i,line in enumerate(txt) if "M11 DATA" in line][0]+1,len(txt))
            for i in lines:
                line = txt[i]
                values = line.split("  ")
                for valuei in range(1,len(values)):
                    linksFlow[links[valuei-1]].append(float(values[valuei]))
            os.remove(resultFile)
            os.remove(M11IN)
            
            lines = ""
            nodes = []
            basins = [row[0] for row in arcpy.da.SearchCursor(msm_Node, "MUID", where_clause = "TypeNo = 2")]
            with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
                for linei,line in enumerate(M11OUTFile):
                    try:
                        if "Node_WL" in line and re.findall("<([^>]+)>",line)[0] in basins:		
                            nodes.append(re.findall("<([^>]+)>",line)[0])
                            line = re.sub("^0","1",line)
                        lines += line
                    except UnicodeDecodeError:
                        lines += line
                        arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                    except Exception as e:
                        raise(e)
                           
            with codecs.open(M11IN,'w','cp1252') as M11INFile:
                M11INFile.write(lines)
            
            call([m11extrapath(), prfFileCopy, resultFile])
            nodesMinWL = {}
            with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
                txt = M11OUTFile.readlines()
            
            findValues = re.compile(r" ([\d\.]+)")
            for vali,val in enumerate(findValues.findall(txt[[i for i,line in enumerate(txt) if "M11 DATA" in line][0]+1])[1:]):
                nodesMinWL[nodes[vali]] = float(val)
            
            os.remove(prfFileCopy)
            os.remove(M11OUT)
            os.remove(resultFile)
            os.remove(M11IN)
            
            msm_LinkFromNode = {}
            msm_LinkToNode = {}
            with arcpy.da.SearchCursor(msm_Link,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]
            with arcpy.da.SearchCursor(msm_Weir,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]
            with arcpy.da.SearchCursor(msm_Orifice,["MUID","FromNode","ToNode"]) as cursor:
                for row in cursor:
                    msm_LinkFromNode[row[0]] = row[1]
                    msm_LinkToNode[row[0]] = row[2]

        arcpy.SetProgressorLabel("Reading ERF-file")
        dataTables = mousereader.readERF(erfFile,"MaxLevel_Ranked",MUIDs)
        # arcpy.AddMessage(dataTables)
        arcpy.SetProgressorLabel("Getting return period of flooding")
        MUIDsTCrit = {}
        MUIDs_elevation_at_crit = {}
        for i,MUID in enumerate(MUIDs):
            try:
                # arcpy.AddMessage(dataTables[i])
                nodeH = dataTables[i]["col2"]
                MUIDs_elevation_at_crit[MUID] = -1

                if MUIDsHCrit[MUID]>np.max(nodeH):
                    MUIDsTCrit[MUID] = 1000
                elif MUIDsHCrit[MUID]<np.min(nodeH):
                    MUIDsTCrit[MUID] = 0
                else:
                    MUIDsTCrit[MUID] = np.interp(MUIDsHCrit[MUID],
                          np.flipud(nodeH),
                          np.flipud(float(observationPeriod)/(np.array(range(len(nodeH)))+1)))
                    # arcpy.AddMessage((MUIDsHCrit[MUID], nodeH))
                if critical_return_period:
                    try:
                        MUIDs_elevation_at_crit[MUID] = np.interp(critical_return_period,
                                                           np.flipud(float(observationPeriod) /
                                                                     (np.array(range(len(nodeH))) + 1)),
                                                           np.flipud(nodeH))
                    except Exception as e:
                        pass
            except Exception as e:
                arcpy.AddError(MUID)
                arcpy.AddError(traceback.format_exc())
        
        arcpy.SetProgressorLabel("Creating manhole result file")
        msm_NodeNew = arcpy.CopyFeatures_management(msm_Node, getAvailableFilename(arcpy.env.scratchGDB + "\NodeFloodReturn",
                                                                          parent=mike_urban_database)).getOutput(0)
        # arcpy.AddMessage(msm_NodeNew)
        arcpy.AddField_management (msm_NodeNew, "TCrit", "DOUBLE", 10, 5)
        with arcpy.da.UpdateCursor(msm_NodeNew,["MUID","TCrit"]) as cursor:
            for row in cursor:
                if row[0] in MUIDsTCrit:
                    row[1] = MUIDsTCrit[row[0]]
                cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False

        def addLayer(layer_source, source, group=None, workspace_type = None, new_name=None,
                     definition_query=None):
            if not workspace_type:
                if ".mdb" in source:
                    workspace_type = "ACCESS_WORKSPACE"
                elif ".shp" in source:
                    workspace_type = "SHAPEFILE_WORKSPACE"
                elif ".gdb" in source:
                    workspace_type = "FILEGDB_WORKSPACE"
            if ".sqlite" in source:
                source_layer = arcpy.mapping.Layer(source)

                if group:
                    arcpy.mapping.AddLayerToGroup(df, group, source_layer, "BOTTOM")
                else:
                    arcpy.mapping.AddLayer(df, source_layer, "TOP")
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
                update_layer.definitionQuery = definition_query

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
                    arcpy.mapping.AddLayer(df, layer, "TOP")
                update_layer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name
                # arcpy.AddMessage((unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                #                                workspace_type, unicode(os.path.basename(source))))
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry", ""))),
                                               workspace_type, unicode(os.path.basename(source)))

        # arcpy.AddMessage((os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Node.lyr", msm_NodeNew, "FILEGDB_WORKSPACE"))
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Node.lyr", msm_NodeNew, workspace_type="FILEGDB_WORKSPACE")
        return
        arcpy.SetProgressorLabel("Creating basin result file")
        exportBasins = arcpy.Select_analysis(mike_urban_database + r"\mu_Geometry\msm_Node", getAvailableFilename(arcpy.env.scratchGDB + "\BasinFloodReturn",
                                                                         parent=mike_urban_database), where_clause = "TypeNo = 2").getOutput(0)
        fields = ["MUID", "GeometryID","CriticalLe" if ".shp" in exportBasins else "CriticalLevel", "Volume", "TCrit", "CatchArea", "CatchImpA", "Discharge", "VolCrit", 'LvlTCrit']

        arcpy.AddField_management(exportBasins, "Discharge", "FLOAT")
        arcpy.AddField_management(exportBasins, "Volume", "FLOAT")
        arcpy.AddField_management(exportBasins, "TCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "VolCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "CritLevel", "FLOAT")
        arcpy.AddField_management(exportBasins, "LvlTCrit", "FLOAT")
        arcpy.AddField_management(exportBasins, "CatchArea", "FLOAT")
        arcpy.AddField_management(exportBasins, "CatchImpA", "FLOAT")


        class Basin:
            def __init__(self, geometryID):
                self.geometryID = geometryID if geometryID else ""
                self.value1 = []
                self.value3 = []
                self.edges = []
                self.permanent_level = None
                self.invert_level = None

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
                if self.permanent_level:
                    return np.interp(np.max(elevations), elevations,
                                     [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations))) - np.interp(
                        self.permanent_level, elevations,
                        [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))
                else:
                    return np.interp(np.max(elevations), elevations, [0] + list(scipy.integrate.cumtrapz(surface_areas, elevations)))

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

        basins = {}
        with arcpy.da.SearchCursor(msm_Node, ["MUID", "GeometryID", "InvertLevel"], where_clause="TypeNo = 2") as cursor:
            for row in cursor:
                basins[row[0]] = Basin(row[1])
                basins[row[0]].invert_level = row[2]
        
        # if traceNetwork:
        #     arcpy.SetProgressorLabel("Creating basin result file - analyzing basin catchment area")
        #     import networkx as nx
        #     breakChainOnNodes = parameters[8].ValueAsText
        #     networkLinks = [msm_Link]
        #     networkLinks.append(msm_Orifice) if "Orifice" in reaches else None
        #     networkLinks.append(msm_Weir) if "Weir" in reaches else None
        #     networkLinks.append(msm_Pump) if "Pump" in reaches else None
        #
        #     network = nx.DiGraph()
        #
        #     hParADict = {}
        #     with arcpy.da.SearchCursor(msm_HParA,["MUID","RedFactor"]) as cursor:
        #         for row in cursor:
        #             hParADict[row[0]] = row[1]
        #
        #     catchmentImperviousnessDict = {}
        #     catchmentReductionFactor = {}
        #     with arcpy.da.SearchCursor(msm_HModA, ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor"]) as cursor:
        #         for row in cursor:
        #             catchmentImperviousnessDict[row[0]] = row[1]
        #             catchmentReductionFactor[row[0]] = hParADict[row[2]] if row[3] == 0 else row[4]
        #
        #     catchmentConnectionDict = {}
        #     with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
        #         for row in cursor:
        #             if row[1] not in catchmentConnectionDict:
        #                 catchmentConnectionDict[row[1]] = [row[0]]
        #             else:
        #                 catchmentConnectionDict[row[1]].append(row[0])
        #
        #     catchmentAreaDict = {}
        #     catchmentPersonsDict = {}
        #     catchmentImperviousAreaDict = {}
        #     catchmentReducedAreaDict = {}
        #     catchmentNetTypeNoDict = {}
        #     catchmentStatus = {}
        #     with arcpy.da.SearchCursor(ms_Catchment, ['MUID','SHAPE@AREA','Area','Persons',"NetTypeNo", 'Element_S']) as cursor:
        #         for row in cursor:
        #             catchmentPersonsDict[row[0]] = row[3] if row[3] is not None else 0
        #             catchmentAreaDict[row[0]] = row[2]*1e4 if row[2] is not None else row[1]
        #             catchmentNetTypeNoDict[row[0]] = row[4]
        #             catchmentStatus[row[0]] = row[5]
        #             if row[0] in catchmentImperviousnessDict:
        #                 catchmentImperviousAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]]
        #                 catchmentReducedAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]] * catchmentReductionFactor[row[0]]
        #             else:
        #                 arcpy.AddWarning("Warning: Could not find model record for Catchment %s" % (row[0]))
        #
        #     nodeTypeDict = {}
        #     nodeTypes = {1:u"Brønd",2:"Bassin",3:u"Udløb"}
        #     with arcpy.da.SearchCursor(msm_Node,["MUID","TypeNo"]) as cursor:
        #         for row in cursor:
        #             network.add_node(row[0])
        #             nodeTypeDict[row[0]] = nodeTypes[row[1]]
        #
        #     weights = {"msm_Link":1, "msm_Pump":1e4, "msm_Orifice":1e4, "msm_Weir":1e4}
        #     for networkLink in networkLinks:
        #         weight = weights[os.path.basename(networkLink)]
        #         with arcpy.da.SearchCursor(networkLink,["FromNode","ToNode","SHAPE@LENGTH"]) as cursor:
        #             for row in cursor:
        #                 network.add_edge(row[0],row[1], weight = weight*row[2])
        #
        #
        #     if not "Basin" in reaches:
        #         for basin in basins.values():
        #             for edge in network.out_edges(basin.MUID):
        #                 network.remove_edge(basin.MUID,edge[1])
        #                 arcpy.AddMessage("Removed edge %s-%s because tracing through basins is disabled" % (basin.MUID,edge[1]))
        #
        #     if breakChainOnNodes:
        #         breakEdges = [edge for edge in network.edges if edge[0] in re.findall("([^'^(),; \n]+)",breakChainOnNodes)]
        #         network.remove_edges_from(breakEdges)
        #         for edge in breakEdges:
        #             arcpy.AddMessage("Removed edge %s-%s because %s is included in list of nodes to end trace at" % (edge[0],edge[1]))
        #
        #     outlets = []
        #     junctions = []
        #     for node in list(network.nodes):
        #         if node:
        #             if not network.out_edges(node):
        #                 outlets.append(node)
        #             if len(network.out_edges(node))>1:
        #                 junctions.append(node)
        #
        #     for source in junctions:
        #         if source != None:
        #             lengths = np.ones((len(outlets),1))*1e9
        #             for i,target in enumerate(outlets):
        #                 try:
        #                     if nx.has_path(network, source, target):
        #                         lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
        #                 except:
        #                     arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
        #             toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[1]
        #             for edge in network.out_edges(source):
        #                 if not edge[1] == toNode:
        #                     network.remove_edge(source,edge[1])
        #                     arcpy.AddMessage("Removed edge %s-%s so that node %s exclusively leads to outlet %s" % (source,edge[1],source,outlets[np.argmin(lengths)]))
        #
        #     for basin in [row[0] for row in arcpy.da.SearchCursor(exportBasins, ["MUID"])]:
        #         nodesUpstream = nx.ancestors(network,basin)
        #         nodesUpstream.add(basin)
        #         catchIDs = [catchmentConnectionDict[n] for n in nodesUpstream if n in catchmentConnectionDict]
        #         catchIDs = [catchID for sublist in catchIDs for catchID in sublist]
        #
        #         basins[basin].total_catchment_area = round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100
        #         basins[basin].total_catchment_area_impervious = round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100

        with arcpy.da.SearchCursor(os.path.join(mike_urban_database, r"ms_TabD"), ["TabID", "Value1", "Value3"],
                                   where_clause="TabID IN ('%s')" % ("', '".join(
                                       [basin.geometryID for basin in basins.values()]))) as cursor:
            for row in cursor:
                basin = [basin for basin in basins.values() if basin.geometryID == row[0]][0]
                basin.value1.append(row[1])
                basin.value3.append(row[2])

        # arcpy.AddMessage(fields)
        arcpy.SetProgressorLabel("Creating basin result file - analyzing basin volume")
        with arcpy.da.UpdateCursor(exportBasins, fields) as cursor:
            for row in cursor:
                if flowFile:
                    basins[row[0]].permanent_level = nodesMinWL[row[0]]

                row[3] = basins[row[0]].max_volume
                try:
                    row[4] = MUIDsTCrit[row[0]]

                except Exception as e:
                    arcpy.AddWarning("Failed to find %s in General Report" % row[0])
                row[7] = 0
                if flowFile:
                    flow = np.array(())
                    for link in [MUID for MUID,FromNode in msm_LinkFromNode.iteritems() if FromNode==row[0]]:
                        if len(flow) == 0:
                            flow = np.array(linksFlow[link])
                        else:
                            flow = flow + np.array(linksFlow[link])
                    if len(flow)==0:
                        arcpy.AddMessage("Could not find link that discharges from basin %s (%s). Assuming intlet is also outlet." % (row[0],row[1]))
                        for link in [MUID for MUID,ToNode in msm_LinkToNode.iteritems() if ToNode==row[0]]:
                            if len(flow) == 0:
                                flow = np.array([a*(-1) for a in linksFlow[link]])
                            else:
                                flow = flow + np.array([a*(-1) for a in linksFlow[link]])
                    # arcpy.AddMessage(flow)
                    row[7] += np.max(flow)
                    
                # row[5] = basins[row[0]].total_catchment_area if basins[row[0]].total_catchment_area else 0
                # row[6] = basins[row[0]].total_catchment_area_impervious if basins[row[0]].total_catchment_area_impervious else 0
                if critical_return_period:
                    row[8] = basins[row[0]].get_volume(MUIDs_elevation_at_crit[row[0]])
                    row[9] = MUIDs_elevation_at_crit[row[0]]
                cursor.updateRow(row)
        
        # if flowFile:
        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin Discharge.lyr", exportBasins, workspace_type="FILEGDB_WORKSPACE")
        # basinLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin Discharge.lyr")
        # # else:
        # #     basinLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Basin.lyr")
        # basinLayer = arcpy.mapping.AddLayer(df, basinLayer)
        # basinLayer = arcpy.mapping.ListLayers(mxd, basinLayer, df)[0]
        # basinLayer.replaceDataSource(os.path.dirname(exportBasins), "SHAPEFILE_WORKSPACE", os.path.basename(exportBasins).split(".")[0])
        #arcpy.RefreshTOC()
        #arcpy.RefreshActiveView()
        return

class DisplayWeirStatistics(object):
    def __init__(self):
        self.label       = "Display Weir Statistics"
        self.description = "Display Weir Statistics"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        htmlFile = arcpy.Parameter(
            displayName="Input LTS ERF File",
            name="htmlFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        htmlFile.filter.list=["html"]
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
            
        parameters = [htmlFile, mike_urban_database]
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        htmlFile = parameters[0].ValueAsText
        mike_urban_database = parameters[1].ValueAsText
        arcpy.env.addOutputsToMap = False
        msm_weir = arcpy.CopyFeatures_management(mike_urban_database + "\msm_Weir", getAvailableFilename(arcpy.env.scratchGDB + "\msm_Weir"))

        weirs = {}
        arcpy.SetProgressorLabel("Getting FROMNODE and TONODE for pipes")
        link_network = networker.NetworkLinks(mike_urban_database, map_only = "weir")
        with arcpy.da.SearchCursor(mike_urban_database + "\msm_Weir",["MUID"]) as cursor:
            try:
                for row in cursor:
                    arcpy.AddMessage(row[0])
                    fromnode = link_network.weirs[row[0]].fromnode
                    tonode = link_network.weirs[row[0]].tonode if link_network.weirs[row[0]].tonode is not None else "0"
                    weirs["%s-%s" % (fromnode,tonode)] = row[0]
                    arcpy.AddMessage("%s-%s" % (fromnode,tonode))
            except Exception as e:
                arcpy.AddError(traceback.format_exc())
                raise(e)
                
        arcpy.SetProgressorLabel("Reading HTML-file")
        with codecs.open(htmlFile,'r',encoding='mbcs') as f:
            htmlFileTxt = f.read().split("\r\n")
            htmlFileTxt = np.array(htmlFileTxt)
            for line in htmlFileTxt:
                line = unicode(line)

        weirStart = [i for i,a in enumerate(htmlFileTxt) if "<H3>General statistics in weirs" in a][0]
        weirEnd = [i for i,a in enumerate(htmlFileTxt) if "TABLE" in a and i > weirStart][0]
        #nodeEnd = [i for i,a in enumerate(htmlFileTxt) if "TABLE" in a][1]
        htmlFileTxtWeirs = htmlFileTxt[weirStart:weirEnd]

        arcpy.SetProgressorLabel("Writing Shapefile")
        if "fromnode" not in [f.name for f in arcpy.ListFields(msm_weir[0])]:
            arcpy.AddField_management(msm_weir[0],"fromnode","TEXT")
        if "tonode" not in [f.name for f in arcpy.ListFields(msm_weir[0])]:
            arcpy.AddField_management(msm_weir[0],"tonode","TEXT")
        arcpy.AddField_management(msm_weir[0],"QVol","FLOAT")
        arcpy.AddField_management(msm_weir[0],"QNo","FLOAT")
        arcpy.AddField_management(msm_weir[0],"QHour","FLOAT")
        
        MUIDsQVol = {}
        MUIDsQNo = {}
        MUIDsQHours = {}
        # arcpy.AddMessage(weirs)
        getMUIDRe = re.compile(r"ALIGN=LEFT>([^<]+)")
        getQs = re.compile(r"<TD>([-0-9 <>\.]+)<\/TD><TD>([-0-9 <>\.]+)<\/TD><TD>([-0-9 <>\.]+)<\/TD><\/TR>$")
        for line in htmlFileTxtWeirs:
            if "ALIGN=LEFT" in line:
                try:
                    # arcpy.AddMessage(weirs["%s-%s" % (getMUIDRe.findall(line)[0],getMUIDRe.findall(line)[1])])
                    # arcpy.AddMessage("%s-%s" % (getMUIDRe.findall(line)[0],getMUIDRe.findall(line)[1]))
                    fromnode = getMUIDRe.findall(line)[0]
                    tonode = getMUIDRe.findall(line)[1] if not "Weir Outlet" in line else "0"
                    MUIDsQVol[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][0])
                    MUIDsQNo[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][1])
                    MUIDsQHours[weirs["%s-%s" % (fromnode, tonode)]] = float(getQs.findall(line)[0][2])
                except Exception as e:
                    arcpy.AddError("Error on line %s" % (line))
                    arcpy.AddError(traceback.format_exc())

        
        with arcpy.da.UpdateCursor(msm_weir[0],["MUID", "QVol", "QNo", "QHour", "fromnode", "tonode"]) as cursor:
            for row in cursor:
                if row[0] not in MUIDsQVol:
                    arcpy.AddError("Error: Could not find results of weir %s in result file" % (row[0]))
                else:
                    row[1] = MUIDsQVol[row[0]]
                    row[2] = MUIDsQNo[row[0]]
                    row[3] = MUIDsQHours[row[0]]
                    row[4] = link_network.weirs[row[0]].fromnode
                    row[5] = link_network.weirs[row[0]].tonode
                    cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        weirLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Weir.lyr")
        weirLayer = arcpy.mapping.AddLayer(df, weirLayer,'TOP')
        weirLayer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        weirLayer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE", os.path.basename(msm_weir[0]).split(".")[0])
        
        weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        #arcpy.RefreshTOC()
        
        return
        
class DisplayFlowStatistics(object):
    def __init__(self):
        self.label       = "Display Flow Statistics"
        self.description = "Display Flow Statistics"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")      
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database with links",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        exportShape = arcpy.Parameter(
            displayName="Output shapefile with link flow statistics",
            name="exportShape",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
            
        use_networker = arcpy.Parameter(
            displayName="Get fromnode and tonode attributes of pipes from spatial analysis?",
            name="use_networker",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="Optional",
            direction="Input")
        use_networker.value = False
        
        export_cad = arcpy.Parameter(
            displayName="Export CAD File with results and outlet flows",
            name="export_cad",
            datatype="File",
            parameterType="Optional",
            category="Additional Settings",
            direction="Output")
        # export_cad.filter.list=["dxf"]
            
        parameters = [erfFile, observationPeriod, mike_urban_database, exportShape, use_networker, export_cad]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        arcpy.env.addOutputsToMap = False
        erfFile = parameters[0].ValueAsText
        mike_urban_database = parameters[2].ValueAsText
        geometryFile = mike_urban_database + "\mu_Geometry\msm_Link"
        msm_Weir = mike_urban_database + "\mu_Geometry\msm_Weir"
        exportShape = parameters[3].ValueAsText
        observationPeriod = float(parameters[1].ValueAsText)
        use_networker = parameters[4].ValueAsText
        export_cad = parameters[5].ValueAsText
        
        MUIDs = {}
        if use_networker == "False":
            with arcpy.da.SearchCursor(geometryFile,["MUID","FROMNODE","TONODE"]) as cursor:
                for row in cursor:
                    MUIDs[row[0]] = "'%s', '%s'" % (row[1],row[2])
            with arcpy.da.SearchCursor(msm_Weir,["MUID","FROMNODE","TONODE"]) as cursor:
                for row in cursor:
                    MUIDs[row[0]] = "'%s', '%s'" % (row[1],row[2])
        else:
            link_network = networker.NetworkLinks(mike_urban_database, map_only = "link weir").links
            for link in link_network.values():
                MUIDs[link.MUID] = "'%s', '%s'" % (link.fromnode, link.tonode)
            
        dataTables = mousereader.readERF(erfFile, "MaxFlow_Ranked", MUIDs.values(), ignore = True)

        msmLinkNew = arcpy.CopyFeatures_management(geometryFile,exportShape)
        
        if "FROMNODE" not in [field.name for field in arcpy.ListFields(msmLinkNew)] and "TONODE" not in [field.name for field in arcpy.ListFields(msmLinkNew)]:
            arcpy.AddField_management(msmLinkNew, "FROMNODE", "TEXT")
            arcpy.AddField_management(msmLinkNew, "TONODE", "TEXT")
            
        RPs = [1, 2, 5, 10, 20, 1000]
        fields = ["T1", "T2", "T5", "T10", "T20", "TMax"]
        for field in fields:
            arcpy.AddField_management (msmLinkNew, field, "FLOAT", field_precision = 8, field_scale = 5)
            
        if use_networker:
            fromnodes = [link.fromnode for link in link_network.values()]
            outlet_MUIDs = [link.MUID for link in link_network.values() if link.tonode not in fromnodes]
        else:
            fromnodes = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["FROMNODE"])]
            outlet_MUIDs = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["MUID","TONODE"]) if row[1] not in fromnodes]
        
        with arcpy.da.UpdateCursor(msmLinkNew,["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"]) as cursor:
            for row in cursor:
                try:
                    if row[0] in MUIDs.keys() and not dataTables[MUIDs.keys().index(row[0])] == None:
                            for i in range(len(RPs)):
                                try:
                                    flows = dataTables[MUIDs.keys().index(row[0])]["col2"]
                                    erfRPs = np.flipud(np.array(observationPeriod)/np.arange(1,len(flows)+1))
                                    row[i+1] = np.interp(RPs[i],erfRPs,np.flipud(flows))*1e3
                                    if use_networker:
                                        row[-3] = link_network[row[0]].fromnode
                                        row[-2] = link_network[row[0]].tonode
                                except Exception as e:
                                    arcpy.AddWarning(e)
                            
                # row[-3] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.X
                # row[-2] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.Y
                    cursor.updateRow(row)   
                except Exception as e:
                    arcpy.AddWarning("Failed on MUID %s: %s" % (row[0], e))
        
        with arcpy.da.InsertCursor(msmLinkNew, ["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"]) as cursor:
                extra_fields = ["FROMNODE","TONODE"] if not use_networker else []                
                row = [None] * len(["MUID"] + fields + ["FROMNODE", "TONODE"] + ["SHAPE@"])
                with arcpy.da.SearchCursor(msm_Weir, ["MUID", "SHAPE@"] + extra_fields) as weir_cursor:
                    for weir_row in weir_cursor:
                        try:
                            if weir_row[0] in MUIDs.keys() and not dataTables[MUIDs.keys().index(weir_row[0])] == None:
                                row[0] = weir_row[0]
                                for i in range(len(RPs)):
                                    try:
                                        flows = dataTables[MUIDs.keys().index(weir_row[0])]["col2"]
                                        erfRPs = np.flipud(np.array(observationPeriod)/np.arange(1,len(flows)+1))
                                        row[i+1] = np.interp(RPs[i],erfRPs,np.flipud(flows))*1e3
                                        row[-1] = weir_row[1]
                                        if use_networker:
                                            row[-3] = link_network[weir_row[0]].fromnode
                                            row[-2] = link_network[weir_row[0]].tonode
                                    except Exception as e:
                                        arcpy.AddWarning(e)
                                    
                        # row[-3] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.X
                        # row[-2] = row[-1].projectAs(arcpy.SpatialReference("WGS 1984")).lastPoint.Y
                                cursor.insertRow(row)   
                        except Exception as e:
                            arcpy.AddWarning("Failed on MUID %s: %s" % (weir_row[0], e))
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        add_layer = arcpy.mapping.Layer(exportShape)
        if (arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database))) 
            and arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database)))[0].isGroupLayer):
            group_layer = arcpy.mapping.ListLayers(mxd, os.path.splitext(os.path.basename(mike_urban_database)))[0]
            arcpy.mapping.AddLayerToGroup(df, group_layer, add_layer, "TOP")
            update_layer = [layer for layer in arcpy.mapping.ListLayers(mxd, add_layer.name, df) if layer.longName == os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring"]
            
            arcpy.mapping.AddLayerToGroup(df, group_layer, add_layer, "TOP")
            update_layer_outlet = [layer for layer in arcpy.mapping.ListLayers(mxd, add_layer.name, df) if layer.longName == os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring udløb"]
            
            arcpy.AddMessage(os.path.splitext(os.path.basename(mike_urban_database))[0] + r"\\" + u"Vandføring")
        else:
            arcpy.mapping.AddLayer(df, add_layer, "TOP")
            update_layer = arcpy.mapping.ListLayers(mxd, add_layer.name, df)[0]
            
            arcpy.mapping.AddLayer(df, add_layer, "TOP")
            update_layer_outlet = arcpy.mapping.ListLayers(mxd, add_layer.name, df)[0]
        source_layer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_LinkQMax.lyr")
        arcpy.mapping.UpdateLayer(df,update_layer,source_layer,False)
        update_layer.replaceDataSource(unicode(add_layer.workspacePath), 'SHAPEFILE_WORKSPACE', unicode(add_layer.datasetName))
        update_layer.definitionQuery = "MUID NOT IN ('%s')" % ("', '".join(outlet_MUIDs))
        
        source_layer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_LinkQMaxOutlet.lyr")
        arcpy.mapping.UpdateLayer(df,update_layer_outlet,source_layer,False)
        update_layer_outlet.replaceDataSource(unicode(add_layer.workspacePath), 'SHAPEFILE_WORKSPACE', unicode(add_layer.datasetName))
        update_layer_outlet.definitionQuery = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))
        
        arcpy.RefreshTOC()
        arcpy.RefreshActiveView()
        
        # if export_cad:
            # import ezdxf
            
            # tonodes = {}
            # if use_networker:
                # fromnodes = [link.fromnode for link in link_network.values()]
                # outlet_MUIDs = [link.MUID for link in link_network.values() if link.tonode not in fromnodes]
                # for link in link_network.values():
                    # tonodes[link.MUID] = link.tonode
            # else:
                # fromnodes = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["FROMNODE"])]
                # outlet_MUIDs = [row[0] for row in arcpy.da.SearchCursor(exportShape, ["MUID","TONODE"]) if row[1] not in fromnodes]
                # with arcpy.da.SearchCursor(exportShape, ["FROMNODE", "MUID"], where_clause = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))):
                    # tonodes[row[1]] = row[0]
                
            
            # doc = ezdxf.new(dxfversion='R2010')
            # doc.layers.new(u'Udloeb', dxfattribs={'color': 255})
            # msp = doc.modelspace()
            # arcpy.AddMessage(outlet_MUIDs)
            # arcpy.AddMessage( "MUID IN ('%s')" % ("', '".join(outlet_MUIDs)))
            # with arcpy.da.SearchCursor(exportShape, ["SHAPE@", "MUID", "T1", "T2", "T5", "T10", "T20", "TMax"], where_clause = "MUID IN ('%s')" % ("', '".join(outlet_MUIDs))) as cursor:
                # for row in cursor:
                    # a = msp.add_mtext(
                    # 'ID: %s%sT1: %1.1f' % ("bay", "\n", ("%1.1f" % (max(row[2],0))).replace(".","."))) #tonodes[row[1]]
                    # a.set_location((row[0].lastPoint.X, row[0].lastPoint.Y))

        # # Save DXF document.
            # doc.saveas(export_cad)
        return
        
class DisplayQFullQMax(object):
    def __init__(self):
        self.label       = "Display Filling Degree of Pipe"
        self.description = "Display Filling Degree of Pipe"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        prfFile = arcpy.Parameter(
            displayName="PRF file",
            name="prfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        prfFile.filter.list=["prf"]
        
        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database with links",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        exportShape = arcpy.Parameter(
            displayName="Output shapefile with Filling Degree of Pipe",
            name="exportShape",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        
        minimumSlope = arcpy.Parameter(
            displayName="Replace low slopes [o/oo]",
            name="minimumSlope",
            datatype="double",
            parameterType="Optional",
            direction="Input")
        minimumSlope.filter.list = [5, 10]
            
        parameters = [prfFile, mike_urban_database, exportShape, minimumSlope]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        arcpy.env.addOutputsToMap = False
        prfFile = parameters[0].ValueAsText
        mike_urban_database = parameters[1].ValueAsText
        geometryFile = mike_urban_database + "\mu_Geometry\msm_Link"
        exportShape = parameters[2].ValueAsText
        minimumSlope = float(parameters[3].ValueAsText)
        
        workingDir = os.path.dirname(__file__)
        M11OUT = workingDir + "\M11.OUT"
        M11IN = workingDir + "\M11.IN"

        prfFileCopy = workingDir + "\prffile.PRF"
        resultFile = workingDir + "\ResultFile.txt"
        replaceLowSlopes = 5

        if not minimumSlope:
            minimumSlope = 0
        copyfile(prfFile,prfFileCopy)
        call([m11extrapath(), prfFileCopy])

        lines = ""
        links = []
        with codecs.open(M11OUT,'r','cp1252') as M11OUTFile:
            for linei,line in enumerate(M11OUTFile):
                try:
                    if "Link_Q" in line:		
                        links.append(re.findall("<([^>]+)>",line)[0])
                        line = re.sub("^0","1",line)
                    lines += line
                except UnicodeDecodeError:
                    lines += line
                    arcpy.AddWarning(u"Line no. %d in file %s contains illegal character." % (linei,M11OUT))
                except Exception as e:
                    raise(e)
                       
        with codecs.open(M11IN,'w','cp1252') as M11INFile:
            M11INFile.write(lines)

        call([m11extrapath(), prfFileCopy, resultFile, "/MAX"])
        os.remove(prfFileCopy)
        os.remove(M11IN)
        os.remove(M11OUT)

        linksFlow = []
        with codecs.open(resultFile,'r','cp1252') as M11OUTFile:
            for linei,line in enumerate(M11OUTFile):
                linksFlow.append(float(re.findall(" +([\-\d\.]+)",line)[0]))
        os.remove(resultFile)
        
        arcpy.CopyFeatures_management(os.path.join(mike_urban_database,"msm_Link"), exportShape)
        arcpy.AddField_management(exportShape, "QMax", "DOUBLE", field_scale = 5, field_precision = 10)
        arcpy.AddField_management(exportShape, "QFull", "DOUBLE", field_scale = 5, field_precision = 10)
        arcpy.AddField_management(exportShape, "Filldeg", "DOUBLE", field_scale = 5, field_precision = 10)
        with arcpy.da.UpdateCursor(exportShape, ["MUID","Diameter","Slope_C","MaterialID","QMax","QFull","FillDeg"]) as cursor:
            for row in cursor:
                try:
                    row[4] = linksFlow[[i for i,a in enumerate(links) if a == row[0]][0]]
                except Exception as e:
                    arcpy.AddMessage("Failed on row %s" % row)
                    raise(e)
                row[5] = colebrookWhite.QFull(row[1],max(row[2]*1e-2,minimumSlope*1e-3),row[3])
                row[6] = max(0,row[4]/row[5])
                cursor.updateRow(row)
            
        
                
        # mxd = arcpy.mapping.MapDocument("CURRENT")
        # df = arcpy.mapping.ListDataFrames(mxd)[0]
        # arcpy.env.addOutputsToMap = False
        # linkLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\QmaxQFull.lyr.lyr")
        # linkLayer = arcpy.mapping.AddLayer(df, linkLayer)
        # linkLayer = arcpy.mapping.ListLayers(mxd, linkLayer, df)[0]
        # arcpy.AddMessage(os.path.dirname(exportShape))
        # arcpy.AddMessage(os.path.dirname(os.path.basename(exportShape).split(".")[0]))
        # linkLayer.replaceDataSource(os.path.dirname(exportShape), "FILEGDB_WORKSPACE", os.path.basename(exportShape).split(".")[0])
        return
        
class DisplayWeirReturnPeriod(object):
    def __init__(self):
        self.label       = "Display Weir Return Period"
        self.description = "Display Weir Return Period"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        erfFile = arcpy.Parameter(
            displayName="ERF file",
            name="erfFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        erfFile.filter.list=["erf"]
        
        observationPeriod = arcpy.Parameter(
            displayName="Observation period of ERF file",
            name="observationPeriod",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        critical_return_period = arcpy.Parameter(
            displayName="Critical Return Period (5 years, 10 years or 20 years)",
            name="critical_return_period",
            datatype="Double",
            parameterType="Optional",
            direction="Input")

        mike_urban_database = arcpy.Parameter(
            displayName="Mike Urban Database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
        parameters = [erfFile, observationPeriod, critical_return_period, mike_urban_database]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered and parameters[0].value and not parameters[1].value:
            with open(parameters[0].ValueAsText,"r") as f:
                txt = f.read()
                observationPeriod = float(re.findall(r"Observation_period =[^,]+,[^,]+,([^\n]+)",txt)[0])
                parameters[1].value = observationPeriod
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        erf_file = parameters[0].ValueAsText
        mike_urban_database = parameters[3].ValueAsText

        return_period = 38
        critical_return_period = 10
        
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        msm_Weir = os.path.join(mike_urban_database, "msm_Weir")

        class Weir:
            def __init__(self, muid, fromnode, tonode):
                self.fromnode = fromnode
                self.tonode = tonode if tonode else "0"
                self.muid  = muid

            result = {}

            @property
            def name(self):
                return "'%s', '%s'" % (self.fromnode, self.tonode)

            @property
            def events_count(self):
                return float(len(self.result['col2'])) if self.result['col2'] else 999

            @property
            def critical_discharge(self):
                # arcpy.AddMessage((self.name, critical_return_period, return_period / (1.0+np.arange(self.events_count)), self.result['col2'], np.min((1.0+np.arange(self.events_count)))))
                return np.interp(critical_return_period, np.flipud(return_period / (1.0+np.arange(self.events_count))), np.flipud(self.result['col2'])) if np.min(return_period / (1.0+np.arange(self.events_count))) <= critical_return_period else 0


        weirs = {}
        link_network = networker.NetworkLinks(mike_urban_database, map_only = "weir")
        for weir in link_network.weirs:
            weirs[weir] = Weir(weir, link_network.weirs[weir].fromnode, link_network.weirs[weir].tonode)


        results = mousereader.readERF(erf_file, "Total_Discharge_Ranked", [weir.name for weir in weirs.values()])
        arcpy.AddMessage((erf_file, "Total_Discharge_Ranked", [weir.name for weir in weirs.values()]))

        for result, weir in zip(results, weirs.values()):
            arcpy.AddMessage((weir,result))
            if result is not None:
                weir.result = result
                print(weir.critical_discharge)
        
        empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
        empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(mike_urban_database))[0]
        
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
            
        msm_WeirNew = arcpy.CopyFeatures_management(msm_Weir, os.path.join(arcpy.env.scratchGDB, "Weir_return_period"))[0]

        arcpy.AddField_management (msm_WeirNew, "TCrit", "DOUBLE", 10, 5)
        arcpy.AddField_management (msm_WeirNew, "QCrit", "DOUBLE", 10, 5)
        
        with arcpy.da.UpdateCursor(msm_WeirNew, ["MUID", "TCrit", "QCrit"]) as cursor:
            for row in cursor:
                if weir.result:
                    # arcpy.AddMessage((row[0], weirs[row[0]].events_count, weirs[row[0]].result["col2"]))
                    row[1] = return_period / weirs[row[0]].events_count
                    row[2] = weirs[row[0]].critical_discharge
                    cursor.updateRow(row)
                
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False

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
                update_layer.definitionQuery = definition_query
                
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
                update_layer = arcpy.mapping.ListLayers(mxd, layer.name, df)[0]
                if definition_query:
                    update_layer.definitionQuery = definition_query
                if new_name:
                    update_layer.name = new_name
                arcpy.AddMessage((unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source))))
                update_layer.replaceDataSource(unicode(os.path.dirname(source.replace(r"\mu_Geometry",""))), workspace_type, unicode(os.path.basename(source)))
                

        addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Weir_return_period.lyr", msm_WeirNew, group = empty_group_layer, workspace_type = "FILEGDB_WORKSPACE")

class DisplayMIKE1DResults(object):
    def __init__(self):
        self.label = "Display MIKE1D Results (alpha)"
        self.description = "Display MIKE1D Results (alpha)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        folder = arcpy.Parameter(
            displayName="Folder",
            name="folder",
            datatype="Folder",
            parameterType="Optional",
            direction="Input")
        if os.path.exists(r"C:\Papirkurv\Resultater"):
            folder.value = r"C:\Papirkurv\Resultater"

        node_featureclass = arcpy.Parameter(
            displayName="Nodes Result File",
            name="node_featureclass",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        node_featureclass.filter.list = ["POINT"]

        reach_featureclass = arcpy.Parameter(
            displayName="Reaches Result File",
            name="reach_featureclass",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        reach_featureclass.filter.list = ["LINE", "POLYLINE"]

        display_type = arcpy.Parameter(
            displayName="Display with fitting symbology",
            name="display_type",
            datatype="GPString",
            parameterType="Optional",
            multiValue=True,
            direction="Input")
        display_type.filter.list = ["Flood volume", "Flood Depth", "Max Elevation / Headloss", "Peak Discharge", "Link Depth Difference", "Total Discharge"]
        display_type.value = ["Flood Volume", "Peak Discharge"]

        # res1d_filepath  = arcpy.Parameter(
        #     displayName="res1d Filepath",
        #     name="res1d_filepath",
        #     datatype="File",
        #     parameterType="Optional",
        #     direction="Input")
        # res1d_filepath.filter.list = ["res1d"]
        #
        # python3_path = arcpy.Parameter(
        #     displayName="Python 3 Filepath",
        #     name="python3_path",
        #     datatype="File",
        #     parameterType="Optional",
        #     category="IF RES1D: Python 3 Path with MIKEIO and Arcpy",
        #     direction="Input")
        # python3_path.filter.list = ["exe"]
        #
        # python3_path.value = r"C:\Users\elnn\AppData\Local\anaconda3\envs\myenv_py3_v2\python.exe"

        parameters = [folder, node_featureclass, reach_featureclass, display_type]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        # res1d_filepath = parameters[4]
        # python3_path = parameters[5]
        if parameters[0].altered and parameters[0].value and not parameters[1].value and not parameters[2].value:
            folder = parameters[0].ValueAsText
            import glob, os

            # Gather both shapefiles and file geodatabases
            candidate_files = glob.glob(os.path.join(folder, "*.shp"))
            candidate_files.extend(glob.glob(os.path.join(folder, "*.gdb")))

            # Get modification times and sort the candidates (newest first)
            file_mod_times = [(file, os.path.getmtime(file)) for file in candidate_files]
            sorted_files = sorted(file_mod_times, key=lambda x: x[1], reverse=True)

            # Assign the newest file containing "nodes" to parameters[1]
            for file, _ in sorted_files:
                if "gdb" in file:
                    files_inside_gdb = [os.path.join(dirpath, f) for dirpath, _, files in
                     arcpy.da.Walk(file,
                                   datatype="FeatureClass") for f in files]
                    for subfile in files_inside_gdb:
                        if "nodes" in subfile.lower():
                            parameters[1].value = os.path.join(file, subfile)
                        elif "reaches" in subfile.lower():
                            parameters[2].value = os.path.join(file, subfile)
                    break
                if "nodes" in file:
                    parameters[1].Value = file
                    break

            # Assign the newest file containing "links" to parameters[2]
            if not parameters[2].Value:
                for file, _ in sorted_files:
                    if "links" in file:
                        parameters[2].Value = file
                        break

        # if res1d_filepath.value:
        #     parameters[0].enabled = False
        #     parameters[1].enabled = False
        #     parameters[2].enabled = False
        #
        # if parameters[0].altered:
        #     parameters[3].value = getAvailableFilename(os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF","_NodeFlood"))
        #     parameters[4].value = getAvailableFilename(
        #         os.path.join(arcpy.env.scratchGDB, os.path.basename(parameters[0].ValueAsText)).replace(".ERF", "_BasinFlood"))
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        nodes_featureclass = parameters[1].ValueAsText
        reaches_featureclass = parameters[2].ValueAsText
        display_type = parameters[3].ValueAsText
        # res1d_filepath = parameters[4].ValueAsText
        # python3_path = parameters[5].ValueAsText
        arcpy.AddMessage(display_type)

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]

        empty_group_mapped = arcpymapping.LayerFile(os.path.dirname(
            os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr") if arcgis_pro else arcpy.mapping.Layer(
            os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
        empty_group = df.addLayer(empty_group_mapped) if arcgis_pro else arcpymapping.AddLayer(df, empty_group_mapped,
                                                                                               "TOP")
        empty_group_layer = df.listLayers('Empty Group')[0] if arcgis_pro else \
            arcpymapping.ListLayers(mxd, "Empty Group", df)[0]
        if (nodes_featureclass and ".gdb" in nodes_featureclass) or (reaches_featureclass and ".gdb" in reaches_featureclass):
            empty_group_layer.name = os.path.basename(os.path.dirname(nodes_featureclass)).replace(".gdb","") if nodes_featureclass else os.path.dirname(reaches_featureclass)

        def addLayer(layer_source, source, group=None, workspace_type="ACCESS_WORKSPACE", new_name=None,
                     definition_query=None):
            if ".sqlite" in source:
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
                        df.addLayer(layer, "TOP")
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
            return update_layer

        if nodes_featureclass:
            if "_spill.shp" in nodes_featureclass:
                arcpy.AddMessage(nodes_featureclass)

                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_spill.lyr",
                                 nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                 workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                 new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                layer.showLabels = True
            else:
                # layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                #          nodes_featureclass.replace(".shp",""), group=None, workspace_type = "SHAPEFILE_WORKSPACE", new_name = os.path.basename(nodes_featureclass).replace(".shp",""))

                if "flood volume" in display_type.lower():
                    layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_floodvol.lyr",
                                     nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                                     workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                                     new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = True

                if "Flood Depth".lower() in display_type.lower():
                    layer = addLayer(
                        os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes_depthdiff.lyr",
                        nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                        workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                        new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = False

                if "headloss" in display_type.lower():
                    layer = addLayer(
                        os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_nodes.lyr",
                        nodes_featureclass.replace(".shp", ""), group=empty_group_layer,
                        workspace_type="SHAPEFILE_WORKSPACE" if "shp" in nodes_featureclass else "FILEGDB_WORKSPACE",
                        new_name=os.path.basename(nodes_featureclass).replace(".shp", ""))
                    layer.showLabels = False

        if reaches_featureclass:
            if "depth difference" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_depthdiff.lyr",
                                 reaches_featureclass.replace(".shp","")    , group=empty_group_layer, workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE", new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                layer.showLabels = False
            if "peak discharge" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links.lyr",
                                 reaches_featureclass.replace(".shp","")    , group=empty_group_layer, workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE", new_name = os.path.basename(reaches_featureclass).replace(".shp",""))
                layer.showLabels = False
            if "total discharge" in display_type.lower():
                layer = addLayer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MIKE1D_results_links_SumQ.lyr",
                                 reaches_featureclass.replace(".shp", ""), group=empty_group_layer,
                                 workspace_type="SHAPEFILE_WORKSPACE" if "shp" in reaches_featureclass else "FILEGDB_WORKSPACE",
                                 new_name=os.path.basename(reaches_featureclass).replace(".shp", ""))
                layer.showLabels = True
        arcpy.RefreshTOC()
        # def addLayer(layer_source, source):
        #     layer = arcpy.mapping.Layer(layer_source)
        #     layer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
        #     layer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        #     layer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
        #                                 os.path.basename(msm_weir[0]).split(".")[0])
        #
        #     weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        #
        # weirLayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Weir.lyr")
        # weirLayer = arcpy.mapping.AddLayer(df, weirLayer, 'TOP')
        # weirLayer = arcpy.mapping.ListLayers(mxd, weirLayer, df)[0]
        # weirLayer.replaceDataSource(os.path.dirname(msm_weir[0]), "FILEGDB_WORKSPACE",
        #                             os.path.basename(msm_weir[0]).split(".")[0])
        #
        # weirLayer.name = os.path.splitext(os.path.basename(htmlFile))[0] + u" Weir Discharge"
        return

class DisplayExtent(object):
    def __init__(self):
        self.label = "Display Extent of Dataframe"
        self.description = "Display Extent of Dataframe"
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Define parameter definitions
        parameters = []

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional

        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        import subprocess
        def copy2clip(txt):
            cmd = 'echo ' + txt.strip() + '|clip'
            return subprocess.check_call(cmd, shell=True)

        copy2clip("[%d, %d, %d, %d]" % (df.extent.lowerLeft.X, df.extent.lowerLeft.Y, df.extent.upperRight.X, df.extent.upperRight.Y))
        arcpy.AddMessage("[%d, %d, %d, %d]" % (df.extent.lowerLeft.X, df.extent.lowerLeft.Y, df.extent.upperRight.X, df.extent.upperRight.Y))
        return