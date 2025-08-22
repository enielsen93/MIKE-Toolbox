# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 11:21:31 2018

@author: eni
"""
import os
import arcpy
import numpy as np
import re
import pythonaddins
import hashlib
import math
import os
import sys
import networkx as nx
import xlwt
functionsPath = [r"K:\Hydrauliske modeller\Makroer & Beregningsark\Functions", r"C:\Dokumenter\Makroer\Functions", os.path.join(os.path.dirname(os.path.dirname(__file__)),"Functions")]
i = 0
while not os.path.exists(functionsPath[i]):
    i += 1
sys.path.append(functionsPath[i])
from mikegraph import PipeNetwork
import ColebrookWhite

def getAvailableFilename(filepath):
    if arcpy.Exists(filepath):
        i = 1
        while arcpy.Exists(filepath + "%d" % i):
            i += 1
        return filepath + "%d" % i
    else:
        return filepath

class Toolbox(object):
    def __init__(self):
        self.label =  "Network Tracer"
        self.alias  = "Network Tracer"
        self.canRunInBackground = False
        # List of tool classes associated with this toolbox
        self.tools = [TraceNetwork, TraceNodeToOutlet, TimeAreaMethod]

class TraceNetwork(object):
    def __init__(self):
        self.label       = "Trace Network and summarize upstream catchment"
        self.description = "Trace Network and summarize upstream catchment"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
            
        # flags = arcpy.Parameter(
            # displayName="Layer to trace upstream network and analyze catchment for",
            # name="flags",
            # datatype="GPFeatureLayer",
            # parameterType="Required",
            # direction="Input")
        # flags.filter.list = ["Simple Junction"]
        
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
        
        sheet = arcpy.Parameter(
            displayName="Excel Sheet with Upstream Catchment Summary",
            name="Excel_sheet",
            datatype="File",
            parameterType="Optional",
            direction="Output")
        # sheet.filter = ["xls"]
        
        writeShapeFile = arcpy.Parameter(
            displayName="Export Node Shapefile with Upstream Catchment Summary?",
            name="includeWasteWater",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        writeShapeFile.value = True
        
        breakChainOnNodes = arcpy.Parameter(
            displayName="End each trace at following node MUIDs (each node should by delimited by a comma: node1, node2)",
            name="breakChainOnNodes",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        breakChainOnNodes.category = "Additional settings"
        
        catchment_where_clause = arcpy.Parameter(
            displayName="Include only catchments that match following SQL Query",
            name="catchment_where_clause",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        catchment_where_clause.category = "Additional settings"
        
        msm_HModA_where_clause = arcpy.Parameter(
            displayName="Include only catchments whose corresponding msm_HModA match following SQL Query",
            name="msm_HModA_where_clause",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        msm_HModA_where_clause.category = "Additional settings"
        
        shapefile = arcpy.Parameter(
            displayName="Add the fields to the following field",
            name="shapefile",
            datatype="GPFeatureLayer",
            parameterType="optional",
            direction="Input")
        shapefile.category = "Additional settings"
        shapefile.enabled = False
        
        field = arcpy.Parameter(
            displayName= "Field to connect to Manhole MUID",
            name="Field",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        field.category = "Additional settings"
        field.enabled = False

        parameters = [MU_database, reaches, sheet, writeShapeFile, breakChainOnNodes, catchment_where_clause, msm_HModA_where_clause, shapefile, field]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        if not parameters[0].valueAsText:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            database = [lyr.dataSource for lyr in arcpy.mapping.ListLayers(mxd) if not lyr.getSelectionSet() and lyr.isFeatureLayer and ".mdb" in lyr.dataSource]
            if database:
                database = re.findall(r"(.+)(?=\.mdb)", database[0])[0] + ".mdb"
                parameters[0].value = database
        if parameters[3].value == True:
            parameters[7].enabled = True
            parameters[8].enabled = True
            if parameters[7].altered and not parameters[8].value:
                parameters[8].filter.list = [f.name for f in arcpy.Describe(parameters[7].value).fields]            
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        # flags = parameters[1].ValueAsText
        reaches = parameters[1].ValueAsText if parameters[1].ValueAsText else []
        writeShapeFile = parameters[3].ValueAsText
        excel_sheet = parameters[2].ValueAsText
        breakChainOnNodes = parameters[4].ValueAsText
        catchment_where_clause = parameters[5].ValueAsText
        msm_HModA_where_clause = parameters[6].ValueAsText
        shapefile = parameters[7].ValueAsText
        field = parameters[8].ValueAsText if parameters[8].ValueAsText else "MUID"
        
        mxd = arcpy.mapping.MapDocument("CURRENT")  
        df = arcpy.mapping.ListDataFrames(mxd)[0]  
        
        msm_Link = os.path.join(MU_database,"msm_Link")
        msm_Node = os.path.join(MU_database,"msm_Node")
        msm_Orifice = os.path.join(MU_database,"msm_Orifice")
        msm_Weir = os.path.join(MU_database,"msm_Weir")
        msm_Pump = os.path.join(MU_database,"msm_Pump")
        msm_CatchCon = os.path.join(MU_database,"msm_CatchCon")
        ms_Catchment = os.path.join(MU_database,"ms_Catchment")
        msm_HModA = os.path.join(MU_database,"msm_HModA")
        msm_HParA = os.path.join(MU_database,"msm_HParA")
        
        networkLinks = [msm_Link]
        networkLinks.append(msm_Orifice) if "Orifice" in reaches else None
        networkLinks.append(msm_Weir) if "Weir" in reaches else None
        networkLinks.append(msm_Pump) if "Pump" in reaches else None

        network = nx.DiGraph()
        # network.add_nodes_from([row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID")])
       
        arcpy.SetProgressorLabel("Reading catchments")
        hParADict = {}
        with arcpy.da.SearchCursor(msm_HParA,["MUID","RedFactor"]) as cursor:
            for row in cursor:
                hParADict[row[0]] = row[1]
        
        catchmentImperviousnessDict = {}
        catchmentReductionFactor = {}
        with arcpy.da.SearchCursor(msm_HModA, ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor"], where_clause = msm_HModA_where_clause) as cursor:
            for row in cursor:
                catchmentImperviousnessDict[row[0]] = row[1]
                catchmentReductionFactor[row[0]] = hParADict[row[2]] if row[3] == 0 else row[4]
        
        catchmentAreaDict = {}
        catchmentPersonsDict = {}
        catchmentImperviousAreaDict = {}
        catchmentReducedAreaDict = {}
        catchmentNetTypeNoDict = {}
        catchmentStatus = {}
        with arcpy.da.SearchCursor(ms_Catchment, ['MUID','SHAPE@AREA','Area','Persons','NetTypeNo', 'Element_S'], where_clause = catchment_where_clause) as cursor:
            for row in cursor:
                if row[0] in catchmentImperviousnessDict:
                    catchmentPersonsDict[row[0]] = row[3] if row[3] is not None else 0
                    catchmentAreaDict[row[0]] = row[2]*1e4 if row[2] is not None else row[1]
                    catchmentNetTypeNoDict[row[0]] = row[4]
                    catchmentStatus[row[0]] = row[5]
                    catchmentImperviousAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]]
                    catchmentReducedAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]] * catchmentReductionFactor[row[0]]
                else:
                    arcpy.AddWarning("Warning: Could not find model record for Catchment %s" % (row[0]))
        
        catchmentConnectionDict = {}
        with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
            for row in cursor:
                if row[0] in catchmentImperviousnessDict and row[0] in catchmentNetTypeNoDict:
                    if row[1] not in catchmentConnectionDict:
                        catchmentConnectionDict[row[1]] = [row[0]]
                    else:
                        catchmentConnectionDict[row[1]].append(row[0])
                
        
        arcpy.SetProgressorLabel("Mapping Network")            
        nodeTypeDict = {}
        basins = []
        nodeTypes = {1:u"Brønd",2:"Bassin",3:u"Udløb"}
        with arcpy.da.SearchCursor(msm_Node,["MUID","TypeNo"]) as cursor:
            for row in cursor:
                network.add_node(row[0])
                nodeTypeDict[row[0]] = nodeTypes[row[1]]
                if nodeTypeDict[row[0]] == "Bassin":
                    basins.append(row[0])
        
        arcpy.SetProgressorLabel("Mapping Network: Getting FROMNODE and TONODE for pipes")
        link_network = PipeNetwork(MU_database)
        weights = {"msm_Link":1, "msm_Pump":1e4, "msm_Orifice":1e4, "msm_Weir":1e4}
        link_network_by_feature_file = {"msm_Link":link_network.links, "msm_Pump":link_network.pumps, "msm_Orifice":link_network.pumps, "msm_Weir":link_network.weirs}
        
        arcpy.SetProgressorLabel("Mapping Network")       
        for networkLink in networkLinks:
            weight = weights[os.path.basename(networkLink)]
            # arcpy.AddMessage(link_network_by_feature_file[os.path.basename(networkLink)])
            for muid, link in link_network_by_feature_file[os.path.basename(networkLink)].items():
                if link.tonode:
                    link_length = 0
                    try:
                        link_length = [row[1] for row in arcpy.da.SearchCursor(networkLink, ["MUID","SHAPE@LENGTH"], where_clause = "MUID = '%s'" % muid)][0]
                    except Exception as e:
                        arcpy.AddMessage(networkLink)
                        arcpy.AddMessage("MUID = '%s'" % muid)
                        
                    network.add_edge(link.fromnode, link.tonode, weight = weight*link_length)
                    
        arcpy.SetProgressorLabel("Mapping Network: Removing links")       
        if not "Basin" in reaches:
            for basin in basins:
                for edge in network.out_edges(basin):
                    network.remove_edge(basin,edge[1])
                    arcpy.AddMessage("Removed edge %s-%s because tracing through basins is disabled" % (basin,edge[1]))            
        
        if breakChainOnNodes:
            breakEdges = [edge for edge in network.edges if edge[0] in re.findall("([^'^(),;\- \n]+)",breakChainOnNodes)]
            network.remove_edges_from(breakEdges)
            for edge in breakEdges:
                if len(edge) == 2:
                    arcpy.AddMessage("Removed edge %s-%s because %s is included in list of nodes to end trace at" % (edge[0],edge[1], edge[0]))
        
        arcpy.SetProgressorLabel("Mapping Network: Making Network unidirectional")       
        outlets = []        
        junctions = []
        for node in list(network.nodes):
            if node:
                if not network.out_edges(node):
                    outlets.append(node)
                if len(network.out_edges(node))>1:
                    junctions.append(node)
        
        for source in junctions:
            if source != None:
                lengths = np.ones((len(outlets),1))*1e9
                for i,target in enumerate(outlets):
                    try:
                        if nx.has_path(network, source, target):
                            lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                    except:
                        arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
                toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[1]
                for edge in network.out_edges(source):
                    if not edge[1] == toNode:
                        network.remove_edge(source,edge[1])
                        arcpy.AddMessage("Removed edge %s-%s so that node %s exclusively leads to outlet %s" % (source,edge[1],source,outlets[np.argmin(lengths)]))
        
        nodes_to_outlet = {}
        for source in network.nodes:
            lengths = np.ones((len(outlets),1))*1e9
            for i,target in enumerate(outlets):
                try:
                    if nx.has_path(network, source, target):
                        lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                except:
                    arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
            toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[-1]
            
            nodes_to_outlet[source] = toNode

        if excel_sheet:
            arcpy.SetProgressorLabel("Writing Excel Sheet")       
            if os.path.exists(excel_sheet):
                os.remove(excel_sheet)
            wb = xlwt.Workbook()
            oplandArk = wb.add_sheet("Opland")
            oplandArk.col(0).width = int(2962*(10/9.44))
            oplandArk.col(2).width = int(2962*(10/9.44))
            oplandArk.col(3).width = int(2962*(17/9.44))
            oplandArk.col(4).width = int(2962*(17/9.44))
            oplandArk.write(0,0,u"Knude", xlwt.easyxf('font: bold on'))
            oplandArk.write(0,1,"Type", xlwt.easyxf('font: bold on'))
            oplandArk.write(0,2,"Opland [ha]", xlwt.easyxf('font: bold on'))
            oplandArk.write(0,3,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
            oplandArk.write(0,4,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))
            oplandArk.write(0,5,"PE", xlwt.easyxf('font: bold on'))
            
            oplandFordeltArk = wb.add_sheet("Opland fordelt")
            oplandFordeltArk.col(0).width = int(2962*(10/9.44))
            oplandFordeltArk.col(2).width = int(2962*(10/9.44))
            oplandFordeltArk.col(3).width = int(2962*(17/9.44))
            oplandFordeltArk.col(4).width = int(2962*(17/9.44))
            oplandFordeltArk.col(6).width = int(2962*(10/9.44))
            oplandFordeltArk.col(7).width = int(2962*(17/9.44))
            oplandFordeltArk.col(8).width = int(2962*(17/9.44))
            oplandFordeltArk.col(10).width = int(2962*(10/9.44))
            oplandFordeltArk.col(11).width = int(2962*(17/9.44))
            oplandFordeltArk.col(12).width = int(2962*(17/9.44))
            oplandFordeltArk.write(1,0,"Knude", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,1,"Type", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write_merge(0, 0, 2, 5,u"Fælles", xlwt.easyxf('align: horiz center; font: bold on'))
            oplandFordeltArk.write_merge(0, 0, 6, 9,u"Spildevand", xlwt.easyxf('align: horiz center; font: bold on'))
            oplandFordeltArk.write_merge(0, 0, 10, 13,u"Regnvand", xlwt.easyxf('align: horiz center; font: bold on'))
            oplandFordeltArk.write(1,2,"Opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,3,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,4,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,5,"PE", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,6,"Opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,7,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,8,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,9,"PE", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,10,"Opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,11,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,12,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))
            oplandFordeltArk.write(1,13,"PE", xlwt.easyxf('font: bold on'))
            
            if 6 in catchmentStatus.values():
                skrivKalibreret = True
                oplandKalibreretArk = wb.add_sheet("Opland kalibreret")
                oplandKalibreretArk.col(0).width = int(2962*(10/9.44))
                oplandKalibreretArk.col(2).width = int(2962*(10/9.44))
                oplandKalibreretArk.col(3).width = int(2962*(17/9.44))
                oplandKalibreretArk.col(4).width = int(2962*(17/9.44))
                oplandKalibreretArk.col(5).width = int(2962*(10/9.44))
                oplandKalibreretArk.col(6).width = int(2962*(17/9.44))
                oplandKalibreretArk.col(7).width = int(2962*(17/9.44))
                oplandKalibreretArk.write_merge(0, 0, 2, 4,u"Status", xlwt.easyxf('align: horiz center; font: bold on'))
                oplandKalibreretArk.write_merge(0, 0, 5, 7,u"Kalibreret", xlwt.easyxf('align: horiz center; font: bold on'))
                oplandKalibreretArk.write(1,0,u"Knude", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,1,"Type", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,2,"Opland [ha]", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,3,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,4,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,5,"Opland [ha]", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,6,u"Befæstet opland [ha]", xlwt.easyxf('font: bold on'))
                oplandKalibreretArk.write(1,7,u"Reduceret opland [ha]", xlwt.easyxf('font: bold on'))   
            else:
                skrivKalibreret = False
                
            i = 0
            nodesUsed = []
            catchmentOutlet = {}
            for node in list(network.nodes):
                if node == "":
                    continue
                nodesUpstream = nx.ancestors(network,node)
                nodesUpstream.add(node)
                catchIDs = [catchmentConnectionDict[n] for n in nodesUpstream if n in catchmentConnectionDict]
                catchIDs = [catchID for sublist in catchIDs for catchID in sublist]
                oplandArk.write(i+1,0, node)
                if node in nodeTypeDict:
                    oplandArk.write(i+1,1, nodeTypeDict[node])
                oplandArk.write(i+1,2, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100)
                oplandArk.write(i+1,3, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100)
                oplandArk.write(i+1,4, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100)
                oplandArk.write(i+1,5, round(np.sum([catchmentPersonsDict[catchID] for catchID in catchIDs])*10)/10)
                
                oplandFordeltArk.write(i+2,0, node)
                if node in nodeTypeDict:
                    oplandFordeltArk.write(i+2,1, nodeTypeDict[node])
                oplandFordeltArk.write(i+2,2, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 1])/1e4*100)/100)
                oplandFordeltArk.write(i+2,3, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 1])/1e4*100)/100)
                oplandFordeltArk.write(i+2,4, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 1])/1e4*100)/100)
                oplandFordeltArk.write(i+2,5, round(np.sum([catchmentPersonsDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 1])*10)/10)
                oplandFordeltArk.write(i+2,6, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 2])/1e4*100)/100)
                oplandFordeltArk.write(i+2,7, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 2])/1e4*100)/100)
                oplandFordeltArk.write(i+2,8, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 2])/1e4*100)/100)
                oplandFordeltArk.write(i+2,9, round(np.sum([catchmentPersonsDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 2])*10)/10)
                oplandFordeltArk.write(i+2,10, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 3])/1e4*100)/100)
                oplandFordeltArk.write(i+2,11, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 3])/1e4*100)/100)
                oplandFordeltArk.write(i+2,12, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 3])/1e4*100)/100)
                oplandFordeltArk.write(i+2,13, round(np.sum([catchmentPersonsDict[catchID] for catchID in catchIDs if catchmentNetTypeNoDict[catchID] == 3])*10)/10)
                
                if skrivKalibreret:
                    oplandKalibreretArk.write(i+2,0, node)
                    if node in nodeTypeDict:
                        oplandKalibreretArk.write(i+2,1, nodeTypeDict[node])
                    oplandKalibreretArk.write(i+2,2, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] != 6])/1e4*100)/100)
                    oplandKalibreretArk.write(i+2,3, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] != 6])/1e4*100)/100)
                    oplandKalibreretArk.write(i+2,4, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] != 6])/1e4*100)/100)
                    oplandKalibreretArk.write(i+2,5, round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] == 6])/1e4*100)/100)
                    oplandKalibreretArk.write(i+2,6, round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] == 6])/1e4*100)/100)
                    oplandKalibreretArk.write(i+2,7, round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs if catchmentStatus[catchID] == 6])/1e4*100)/100)
                i += 1
            wb.save(excel_sheet)
        
        def addField(shapefile, field_name, datatype):
            i = 1
            while field_name in [f.name for f in arcpy.Describe(shapefile).fields]:
                field_name = "%s_%d" % (field_name, i)
            arcpy.AddField_management(shapefile,field_name,datatype)
            return field_name
        
        if writeShapeFile:
            arcpy.SetProgressorLabel("Writing Shapefile")
            if shapefile:
                msm_NodeFeatureLayer = shapefile
            else:
                msm_NodeFeatureLayer = str(arcpy.CopyFeatures_management (msm_Node, getAvailableFilename(arcpy.env.scratchGDB + "\msm_Node")))
            opl_field = addField(msm_NodeFeatureLayer, "Opl", "FLOAT")
            befopl_field = addField(msm_NodeFeatureLayer, "BefOpl", "FLOAT")
            redopl_field = addField(msm_NodeFeatureLayer, "RedOpl", "FLOAT")
            pe_field = addField(msm_NodeFeatureLayer, "PE", "FLOAT")
            outlet_field = addField(msm_NodeFeatureLayer, "Outlet", "STRING")
            with arcpy.da.UpdateCursor(msm_NodeFeatureLayer,["MUID",opl_field,befopl_field, redopl_field, pe_field,outlet_field,field]) as cursor:
                for row in cursor:
                    node = row[6]
                    nodesUpstream = nx.ancestors(network,node)
                    nodesUpstream.add(node)
                    catchIDs = [catchmentConnectionDict[n] for n in nodesUpstream if n in catchmentConnectionDict]
                    catchIDs = [catchID for sublist in catchIDs for catchID in sublist]
                    row[1] = round(np.sum([catchmentAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100
                    row[2] = round(np.sum([catchmentImperviousAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100
                    row[3] = round(np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs])/1e4*100)/100
                    row[4] = round(np.sum([catchmentPersonsDict[catchID] for catchID in catchIDs])*10)/10
                    row[5] = nodes_to_outlet[node] if node in nodes_to_outlet else ""
                    cursor.updateRow(row)
            if not shapefile:
                addLayer = arcpy.mapping.Layer(msm_NodeFeatureLayer)
                arcpy.mapping.AddLayer(df, addLayer, "TOP")
                updatelayer = arcpy.mapping.ListLayers(mxd, os.path.basename(msm_NodeFeatureLayer), df)[0]
                sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\NetworkTracer_Manhole w. Catchment.lyr")
                arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
                updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))
        return

class TraceNodeToOutlet(object):
    def __init__(self):
        self.label       = "Trace Node to Outlet"
        self.description = "Trace Node to Outlet"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        
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
        
        writeShapeFile = arcpy.Parameter(
            displayName="Export Node Shapefile with outlet?",
            name="writeShapeFile",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        writeShapeFile.value = True

        parameters = [MU_database, reaches, writeShapeFile]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        # flags = parameters[1].ValueAsText
        reaches = parameters[1].ValueAsText
        writeShapeFile = parameters[2].ValueAsText
        
        mxd = arcpy.mapping.MapDocument("CURRENT")  
        df = arcpy.mapping.ListDataFrames(mxd)[0]  
        
        msm_Link = os.path.join(MU_database,"msm_Link")
        msm_Node = os.path.join(MU_database,"msm_Node")
        msm_Orifice = os.path.join(MU_database,"msm_Orifice")
        msm_Weir = os.path.join(MU_database,"msm_Weir")
        msm_Pump = os.path.join(MU_database,"msm_Pump")
        
        networkLinks = [msm_Link]
        networkLinks.append(msm_Orifice) if "Orifice" in reaches else None
        networkLinks.append(msm_Weir) if "Weir" in reaches else None
        networkLinks.append(msm_Pump) if "Pump" in reaches else None

        network = nx.DiGraph()
        # network.add_nodes_from([row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID")])

        nodeTypeDict = {}
        basins = []
        nodeTypes = {1:u"Brønd",2:"Bassin",3:u"Udløb"}
        with arcpy.da.SearchCursor(msm_Node,["MUID","TypeNo"]) as cursor:
            for row in cursor:
                network.add_node(row[0])
                nodeTypeDict[row[0]] = nodeTypes[row[1]]
                if nodeTypeDict[row[0]] == "Bassin":
                    basins.append(row[0])
        
        link_network = PipeNetwork(MU_database)
        weights = {"msm_Link":1, "msm_Pump":1e4, "msm_Orifice":1e4, "msm_Weir":1e4}
        link_network_by_feature_file = {"msm_Link":link_network.links, "msm_Pump":link_network.pumps, "msm_Orifice":link_network.pumps, "msm_Weir":link_network.weirs}
        for networkLink in networkLinks:
            weight = weights[os.path.basename(networkLink)]
            # arcpy.AddMessage(link_network_by_feature_file[os.path.basename(networkLink)])
            for muid, link in link_network_by_feature_file[os.path.basename(networkLink)].items():
                if link.tonode:
                    link_length = 0
                    try:
                        link_length = [row[1] for row in arcpy.da.SearchCursor(networkLink, ["MUID","SHAPE@LENGTH"], where_clause = "MUID = '%s'" % muid)][0]
                    except Exception as e:
                        pass    
                        # arcpy.AddMessage(networkLink)
                        # arcpy.AddMessage("MUID = '%s'" % muid)
                        # arcpy.AddWarning(e)
                        
                    network.add_edge(link.fromnode, link.tonode, weight = weight*link_length)
                        
        if not "Basin" in reaches:
            for basin in basins:
                for edge in network.out_edges(basin):
                    network.remove_edge(basin,edge[1])
                    arcpy.AddMessage("Removed edge %s-%s because tracing through basins is disabled" % (basin,edge[1]))            
        
        outlets = []        
        junctions = []
        for node in list(network.nodes):
            if node:
                if not network.out_edges(node):
                    outlets.append(node)
                if len(network.out_edges(node))>1:
                    junctions.append(node)
        
        for source in junctions:
            if source != None:
                lengths = np.ones((len(outlets),1))*1e9
                for i,target in enumerate(outlets):
                    try:
                        if nx.has_path(network, source, target):
                            lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                    except:
                        arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
                toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[1]
                for edge in network.out_edges(source):
                    if not edge[1] == toNode:
                        network.remove_edge(source,edge[1])
                        arcpy.AddMessage("Removed edge %s-%s so that node %s exclusively leads to outlet %s" % (source,edge[1],source,outlets[np.argmin(lengths)]))
        
        nodes_to_outlet = {}
        for source in network.nodes:
            
            lengths = np.ones((len(outlets),1))*1e9
            for i,target in enumerate(outlets):
                try:
                    if nx.has_path(network, source, target):
                        lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                except:
                    arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
            toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[-1]
            
            nodes_to_outlet[source] = toNode
        arcpy.AddMessage(nodes_to_outlet)
            
        
        if writeShapeFile:
            msm_NodeFeatureLayer = str(arcpy.CopyFeatures_management (msm_Node, getAvailableFilename(arcpy.env.scratchGDB + "\msm_Node")))
            arcpy.AddField_management(msm_NodeFeatureLayer,"Outlet","STRING")
            with arcpy.da.UpdateCursor(msm_NodeFeatureLayer,["MUID","OUTLET"]) as cursor:
                for row in cursor:
                    node = row[0]
                    row[1] = nodes_to_outlet[row[0]]
                    cursor.updateRow(row)
                    
            addLayer = arcpy.mapping.Layer(msm_NodeFeatureLayer)
            arcpy.mapping.AddLayer(df, addLayer, "TOP")
            # updatelayer = arcpy.mapping.ListLayers(mxd, os.path.basename(msm_NodeFeatureLayer), df)[0]
            # sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\NetworkTracer_Manhole w. Catchment.lyr")
            # arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
            # updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))
        return

class TimeAreaMethod(object):
    def __init__(self):
        self.label       = "Time Area Method Alpha"
        self.description = "Time Area Method Alpha"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        MU_database = arcpy.Parameter(
            displayName="Mike Urban database",
            name="database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        manhole_target = arcpy.Parameter(
            displayName="Node to trace network to",
            name="manhole_target",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
            
        # flags = arcpy.Parameter(
            # displayName="Layer to trace upstream network and analyze catchment for",
            # name="flags",
            # datatype="GPFeatureLayer",
            # parameterType="Required",
            # direction="Input")
        # flags.filter.list = ["Simple Junction"]
        
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
        
        breakChainOnNodes = arcpy.Parameter(
            displayName="End each trace at following node MUIDs (each node should by delimited by a comma: node1, node2)",
            name="breakChainOnNodes",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        breakChainOnNodes.category = "Additional settings"
        
        catchment_where_clause = arcpy.Parameter(
            displayName="Include only catchments that match following SQL Query",
            name="catchment_where_clause",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        catchment_where_clause.category = "Additional settings"
        
        msm_HModA_where_clause = arcpy.Parameter(
            displayName="Include only catchments whose corresponding msm_HModA match following SQL Query",
            name="msm_HModA_where_clause",
            datatype="GPString",
            parameterType="optional",
            direction="Input")
        msm_HModA_where_clause.category = "Additional settings"
        
        parameters = [MU_database, manhole_target, reaches, breakChainOnNodes, catchment_where_clause, msm_HModA_where_clause]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].valueAsText and not parameters[1].valueAsText:
            MU_database = parameters[0].ValueAsText
            msm_Node = os.path.join(MU_database,"msm_Node")
            parameters[1].Filter.List = np.sort([row[0] for row in arcpy.da.SearchCursor(msm_Node, ["MUID"])])
        return

    def updateMessages(self, parameters): #optional       
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].ValueAsText
        # flags = parameters[1].ValueAsText
        target_manhole = parameters[1].ValueAsText
        reaches = parameters[2].ValueAsText if parameters[2].ValueAsText else []
        breakChainOnNodes = parameters[3].ValueAsText
        catchment_where_clause = parameters[4].ValueAsText
        msm_HModA_where_clause = parameters[5].ValueAsText
        
        mxd = arcpy.mapping.MapDocument("CURRENT")  
        df = arcpy.mapping.ListDataFrames(mxd)[0]  
        
        msm_Link = os.path.join(MU_database,"msm_Link")
        msm_Node = os.path.join(MU_database,"msm_Node")
        msm_Orifice = os.path.join(MU_database,"msm_Orifice")
        msm_Weir = os.path.join(MU_database,"msm_Weir")
        msm_Pump = os.path.join(MU_database,"msm_Pump")
        msm_CatchCon = os.path.join(MU_database,"msm_CatchCon")
        ms_Catchment = os.path.join(MU_database,"ms_Catchment")
        msm_HModA = os.path.join(MU_database,"msm_HModA")
        msm_HParA = os.path.join(MU_database,"msm_HParA")
        
        networkLinks = [msm_Link]
        networkLinks.append(msm_Orifice) if "Orifice" in reaches else None
        networkLinks.append(msm_Weir) if "Weir" in reaches else None
        networkLinks.append(msm_Pump) if "Pump" in reaches else None

        network = nx.DiGraph()
        # network.add_nodes_from([row[0] for row in arcpy.da.SearchCursor(msm_Node,"MUID")])
       
        arcpy.SetProgressorLabel("Reading catchments")
        hParADict = {}
        with arcpy.da.SearchCursor(msm_HParA,["MUID","RedFactor"]) as cursor:
            for row in cursor:
                hParADict[row[0]] = row[1]
        
        catchmentImperviousnessDict = {}
        catchmentReductionFactor = {}
        with arcpy.da.SearchCursor(msm_HModA, ["CatchID", "ImpArea", "ParAID", "LocalNo", "RFactor"], where_clause = msm_HModA_where_clause) as cursor:
            for row in cursor:
                catchmentImperviousnessDict[row[0]] = row[1]
                catchmentReductionFactor[row[0]] = hParADict[row[2]] if row[3] == 0 else row[4]
        
        catchmentAreaDict = {}
        catchmentPersonsDict = {}
        catchmentImperviousAreaDict = {}
        catchmentReducedAreaDict = {}
        catchmentNetTypeNoDict = {}
        catchmentStatus = {}
        with arcpy.da.SearchCursor(ms_Catchment, ['MUID','SHAPE@AREA','Area','Persons','NetTypeNo', 'Element_S'], where_clause = catchment_where_clause) as cursor:
            for row in cursor:
                if row[0] in catchmentImperviousnessDict:
                    catchmentPersonsDict[row[0]] = row[3] if row[3] is not None else 0
                    catchmentAreaDict[row[0]] = row[2]*1e4 if row[2] is not None else row[1]
                    catchmentNetTypeNoDict[row[0]] = row[4]
                    catchmentStatus[row[0]] = row[5]
                    catchmentImperviousAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]]
                    catchmentReducedAreaDict[row[0]] = catchmentImperviousnessDict[row[0]]/100 * catchmentAreaDict[row[0]] * catchmentReductionFactor[row[0]]
                else:
                    arcpy.AddWarning("Warning: Could not find model record for Catchment %s" % (row[0]))
        
        catchmentConnectionDict = {}
        with arcpy.da.SearchCursor(msm_CatchCon, ["CatchID", "NodeID"]) as cursor:
            for row in cursor:
                if row[0] in catchmentImperviousnessDict and row[0] in catchmentNetTypeNoDict:
                    if row[1] not in catchmentConnectionDict:
                        catchmentConnectionDict[row[1]] = [row[0]]
                    else:
                        catchmentConnectionDict[row[1]].append(row[0])
                
        
        arcpy.SetProgressorLabel("Mapping Network")            
        nodeTypeDict = {}
        basins = []
        nodeTypes = {1:u"Brønd",2:"Bassin",3:u"Udløb"}
        with arcpy.da.SearchCursor(msm_Node,["MUID","TypeNo"]) as cursor:
            for row in cursor:
                network.add_node(row[0])
                nodeTypeDict[row[0]] = nodeTypes[row[1]]
                if nodeTypeDict[row[0]] == "Bassin":
                    basins.append(row[0])
        
        arcpy.SetProgressorLabel("Mapping Network: Getting FROMNODE and TONODE for pipes")
        link_network = PipeNetwork(MU_database)
        weights = {"msm_Link":1, "msm_Pump":1e4, "msm_Orifice":1e4, "msm_Weir":1e4}
        link_network_by_feature_file = {"msm_Link":link_network.links, "msm_Pump":link_network.pumps, "msm_Orifice":link_network.pumps, "msm_Weir":link_network.weirs}
        
        arcpy.SetProgressorLabel("Mapping Network")       
        for networkLink in networkLinks:
            weight = weights[os.path.basename(networkLink)]
            # arcpy.AddMessage(link_network_by_feature_file[os.path.basename(networkLink)])
            for muid, link in link_network_by_feature_file[os.path.basename(networkLink)].items():
                if link.tonode:
                    link_length = 0
                    try:
                        link_length = [row[1] for row in arcpy.da.SearchCursor(networkLink, ["MUID","SHAPE@LENGTH"], where_clause = "MUID = '%s'" % muid)][0]
                    except Exception as e:
                        arcpy.AddMessage(networkLink)
                        arcpy.AddMessage("MUID = '%s'" % muid)
                        
                    network.add_edge(link.fromnode, link.tonode, weight = weight*link_length)
                    
        arcpy.SetProgressorLabel("Mapping Network: Removing links")       
        if not "Basin" in reaches:
            for basin in basins:
                for edge in network.out_edges(basin):
                    network.remove_edge(basin,edge[1])
                    arcpy.AddMessage("Removed edge %s-%s because tracing through basins is disabled" % (basin,edge[1]))            
        
        if breakChainOnNodes:
            breakEdges = [edge for edge in network.edges if edge[0] in re.findall("([^'^(),;\- \n]+)",breakChainOnNodes)]
            network.remove_edges_from(breakEdges)
            for edge in breakEdges:
                if len(edge) == 2:
                    arcpy.AddMessage("Removed edge %s-%s because %s is included in list of nodes to end trace at" % (edge[0],edge[1], edge[0]))
        
        arcpy.SetProgressorLabel("Mapping Network: Making Network unidirectional")       
        outlets = []        
        junctions = []
        for node in list(network.nodes):
            if node:
                if not network.out_edges(node):
                    outlets.append(node)
                if len(network.out_edges(node))>1:
                    junctions.append(node)
        
        for source in junctions:
            if source != None:
                lengths = np.ones((len(outlets),1))*1e9
                for i,target in enumerate(outlets):
                    try:
                        if nx.has_path(network, source, target):
                            lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                    except:
                        arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
                toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[1]
                for edge in network.out_edges(source):
                    if not edge[1] == toNode:
                        network.remove_edge(source,edge[1])
                        arcpy.AddMessage("Removed edge %s-%s so that node %s exclusively leads to outlet %s" % (source,edge[1],source,outlets[np.argmin(lengths)]))
        
        nodes_to_outlet = {}
        for source in network.nodes:
            lengths = np.ones((len(outlets),1))*1e9
            for i,target in enumerate(outlets):
                try:
                    if nx.has_path(network, source, target):
                        lengths[i] = (nx.bellman_ford_path_length(network, source, target, weight="weight"))
                except:
                    arcpy.AddError("Failed upon tracing network from %s to %s" % (source,target))
            toNode = nx.bellman_ford_path(network, source, outlets[np.argmin(lengths)])[-1]
            
            nodes_to_outlet[source] = toNode
            
        arcpy.SetProgressorLabel("Calculating full velocity of pipe")
        msm_Link_TravelTime = {}
        with arcpy.da.SearchCursor(msm_Link, ["MUID", "Slope_C", "Diameter", "Length_C", "Length"]) as cursor:
            for row in cursor:
                try:
                    VFull = colebrookWhite.QFull(row[2], row[1]/1e2 if row[1]/1e2>1e-3 else 1e-3, "PL")/((row[2]/2)**2*3.1415)
                except Exception as e:
                    arcpy.AddWarning(row)
                    arcpy.AddWarning([type(field) for field in row])
                    VFull = 1
                length = row[4] if row[4] else row[3]
                length = 10 if not length else length
                msm_Link_TravelTime[row[0]] = length/VFull
        
        arcpy.SetProgressorLabel("Tracing")       
        time_delays = {}
        for source in network.nodes:
            if nx.has_path(network, source, target_manhole):
                arcpy.AddMessage("%s has path to %s" % (source, target_manhole))
                # time_delay = mnx
                path = nx.shortest_path(network, source, target_manhole)
                time_delay = 0
                for link in path:
                    time_delay += msm_Link_TravelTime[row[0]]/60
                time_delays[source] = time_delay
        arcpy.AddMessage(time_delays)
        
        
        
        def time_area(rain_event, conc_time, travel_time):
            runoff = np.zeros(np.shape(rain_event))
            for time_i , rain_intensity in enumerate(rain_event):
                time_i_adjusted = time_i - travel_time
                print(max(time_i_adjusted-conc_time,0),max(0,time_i_adjusted))
                rain = rain_event[max(time_i_adjusted-conc_time,0):max(0,time_i_adjusted)]
                runoff[time_i] = np.sum(rain)/conc_time if rain.any() else 0
            return runoff
		
		import matplotlib
		matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
        box_rains = {5 : 36.63306334, 10 : 26.8206342, 30 : 13.94609554, 60 : 8.63879167, 180 : 3.842244756, 360 : 2.269104561}
        for duration, intensity in box_rains.items():
            rain_event = np.concatenate((np.zeros(10), np.ones(duration)*intensity, np.zeros(30)),axis=0)
        #rain_event = np.concatenate((np.zeros(10), np.ones(30)*10, np.zeros(80)),axis=0)
            runoffs = []
            rain_event_instantaneous = []
            for MUID, time_delay in time_delays.items():
                if MUID in catchmentConnectionDict:
                    catchIDs = catchmentConnectionDict[MUID]
                    conc_time = 7
                    #arcpy.AddMessage(catchIDs)
                    red_area = np.sum([catchmentReducedAreaDict[catchID] for catchID in catchIDs])
                    runoff = time_area(rain_event, conc_time, time_delay)/1e6*red_area*1e3
                    #arcpy.AddMessage(runoff)
                    #arcpy.AddMessage(time_delay)
                    runoffs.append(runoff)
                    rain_event_instantaneous.append(rain_event/1e6*red_area*1e3)
            
            #plt.step(np.arange(0,len(runoff)), rain_event)
            plt.figure()
            plt.step(np.arange(0,len(runoff)), np.sum(np.array(runoffs), axis=0))
            plt.step(np.arange(0,len(runoff)), np.sum(np.array(rain_event_instantaneous), axis=0))
            plt.title("%d: %d" % (duration, np.max(np.sum(np.array(runoffs), axis=0))))
            plt.ylabel("L/s")
            plt.savefig(r"C:\Papirkurv\%s-%d.png" % (target_manhole, duration))
