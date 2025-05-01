# -*- coding: utf-8 -*-


# Tool for converting Dandas .XMLs to ArcGIS features
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy
import os
import re
from xml.dom import minidom
import xml.etree.ElementTree as ET  
from arcpy import env
import csv
import numpy as np  
import sys
from collections import OrderedDict
import traceback
import pythonaddins
import time
import sqlite3

def statusUpdate(message,tic):
    arcpy.AddMessage("%d seconds: %s" % (time.time()-tic, message))
    arcpy.SetProgressorLabel(message)

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
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Dandas To ArcGIS"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Dandas2MULinks, GPSManholes, DDS2ArcGISMU, DDS2ArcGISRK, DDS2Tilbudsliste, CopyMikeUrbanFeatures] # Dandas2ArcGIS
         
class Dandas2MULinks(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "1) Convert Dandas .XMLs to Mike Urban Features"
        self.description = "1) Convert Dandas .XMLs to Mike Urban Features"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        
        dandas_knuder = arcpy.Parameter(
            displayName="Dandas Knuder XML-file:",
            name="dandas_knuder",
            datatype="File",
            parameterType="Required",
            direction="Input")
        dandas_knuder.filter.list=["xml"]
        
        dandas_ledninger = arcpy.Parameter(
            displayName="Dandas Ledninger XML-file:",
            name="dandas_ledninger",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        dandas_ledninger.filter.list=["xml"]
                    
        afloebkodeparameter = arcpy.Parameter(
            displayName="Include only following sewer type:",
            name="afloebkodeparameter",
            datatype="GPString",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        afloebkodeparameter.filter.type = "ValueList"  
        afloebkodeparameter.filter.list = ["Waste Water", "Storm Drain", "Combined Sewer", "Drainage"]
        
        afloebkategori = arcpy.Parameter(
            displayName="Include also the following:",
            name="afloebkategori",
            datatype="GPString",
            parameterType="Optional",
            multiValue = True,
            direction="Input")
        afloebkategori.filter.type = "ValueList"  
        afloebkategori.filter.list = ["Stikledninger", "Ukendt", "Internt"]
        afloebkategori.value = ["Ukendt"]
                    
        coordinate_system = arcpy.Parameter(
            displayName="Coordinate System of Dandas database:",
            name="coordinate_system",
            category="Additional Settings",
            datatype="GPCoordinateSystem",
            parameterType="Optional",
            direction="Input")

        only_import_extent = arcpy.Parameter(
            displayName="Only import elements located within the extent of the dataframe:",
            name="only_import_extent",
            category="Additional Settings",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")

        import_catchments = arcpy.Parameter(
            displayName="Import catchments:",
            name="import_catchments",
            category="Additional Settings",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")

        ignore_elevation_difference = arcpy.Parameter(
            displayName="Remove elevation difference for manhole / link if less than [m]:",
            name="ignore_elevation_difference",
            category="Additional Settings",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        ignore_elevation_difference.value = 0.15
		
        params = [dandas_knuder, dandas_ledninger, afloebkodeparameter, afloebkategori, coordinate_system, only_import_extent, import_catchments, ignore_elevation_difference]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        if parameters[0].ValueAsText and not parameters[1].ValueAsText:
            if re.search(r"(?<=.)knude(?=.)", os.path.basename(parameters[0].ValueAsText), flags=re.IGNORECASE):
                ledninger_filename = re.sub(r"(?<=.)knude(?=.)", "Ledning", os.path.basename(parameters[0].ValueAsText), flags=re.IGNORECASE)
                ledninger_filepath = os.path.join(os.path.dirname(parameters[0].ValueAsText), ledninger_filename)
                if os.path.exists(ledninger_filepath):
                    parameters[1].Value = ledninger_filepath
        return

    def updateMessages(self, parameters):
        # if parameters[0].value:  
            # with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    # parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
        # else:  
            # parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
        arcpy.env.outputCoordinateSystem = arcpy.SpatialReference("ETRS 1989 UTM Zone 32N")
        
        tic = time.time()
        arcpy.SetProgressor("default","Reading Knuder .XML file")        
        dandas_knuder = parameters[0].ValueAsText
        dandas_ledninger = parameters[1].ValueAsText
        afloebkodeparameter = parameters[2].ValueAsText
        afloebkategori = parameters[3].ValueAsText
        coordinate_system = parameters[4].Value
        only_import_extent = parameters[5].Value
        import_catchments = parameters[6].Value
        ignore_elevation_difference = parameters[7].Value

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        extent = df.extent

        MIKE_folder = os.path.join(os.path.dirname(arcpy.env.scratchGDB), "MIKE URBAN")
        if not os.path.exists(MIKE_folder):
            os.mkdir(MIKE_folder)
        MIKE_gdb = os.path.join(MIKE_folder, "DDS_Import")
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
                MIKE_gdb = os.path.join(MIKE_folder, "%s_%d" % ("DDS_Import", dir_ext))
        arcpy.env.scratchWorkspace = MIKE_gdb

        if afloebkategori is None:
            afloebkategori = []
        statusUpdate("Reading %s" % (os.path.basename(dandas_knuder)),tic)
        with open(dandas_knuder,"r") as infile:
            txt = infile.read()
            txt = re.sub(r' xmlns="[^"]+"',"",txt)
        
        statusUpdate("Converting %s to XML Tree" % (os.path.basename(dandas_knuder)),tic)
        tree = ET.fromstring(txt)
        
        nodes = tree.findall("Knude")
        
        statusUpdate("Creating manhole shapefile", tic)
        msm_Node = getAvailableFilename(arcpy.env.scratchGDB + "\msm_Node")
        arcpy.management.CreateFeatureclass(arcpy.env.scratchGDB, os.path.basename(msm_Node), geometry_type = "POINT", 
                                            template = os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\msm_Node")
                                            
        if coordinate_system:
            arcpy.DefineProjection_management(msm_Node, coordinate_system) 
        # arcpy.CopyFeatures_management(os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\msm_Node", msm_Node)
        
        statusUpdate("Adding fields to manholes shapefile", tic)
        msm_Node_Template = OrderedDict()
        templateFile = os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\msm_Node_Template"
        msm_Node_Fields = [a.name for a in arcpy.ListFields(templateFile) if a.name not in ["FID","Shape","SHAPE","OBJECTID"]]
        with arcpy.da.SearchCursor(templateFile,msm_Node_Fields) as cursor:
            for row in cursor:
                for i in range(len(msm_Node_Fields)):
                    msm_Node_Template[msm_Node_Fields[i]] = row[i]

        msm_Node_Cursor = arcpy.da.InsertCursor(msm_Node, ["SHAPE@XY"] + msm_Node_Fields)
        ID = 0
        knudenavn = []

        manhole_table_dict = {0: u"Uoplyst", 1: u"Brønd (standard)", 2: u"Rensebrønd", 3: u"Tømme-/aftapningsbrønd", 4: u"Spulebrønd", 5: u"Ventilbrønd", 6: u"Udluftningsbrønd", 7: u"Målerbrønd", 8: u"Sivebrønd", 9: u"Nedløbsbrønd", 10: u"Samlebrønd", 11: u"Rendestensbrønd", 12: u"Nedgangsbrønd", 13: u"Tilslutningsbrønd", 14: u"Etagebrønd", 50: u"Andet"}

        statusUpdate("Reading manholes", tic)
        nodesDict = {}
        afloebkode = ["Waste Water", "Storm Drain", "Combined Sewer", "Drainage"]
        TypeAfloebKode = 3
        arcpy.SetProgressor("step","Reading Manholes (%d)" % (len(nodes)), 0, len(nodes), 1)
        for nodei in range(len(nodes)):
            arcpy.SetProgressorPosition(nodei)  
            try:
                if only_import_extent:
                    x_coordinate, y_coordinate = float(nodes[nodei].find("XKoordinat").text), float(nodes[nodei].find("YKoordinat").text)
                    # arcpy.AddMessage((x_coordinate > extent.lowerLeft.X, x_coordinate < extent.upperRight.X, y_coordinate > extent.lowerLeft.Y, y_coordinate < extent.upperRight.Y))
                    if not all((x_coordinate > extent.lowerLeft.X, x_coordinate < extent.upperRight.X, y_coordinate > extent.lowerLeft.Y, y_coordinate < extent.upperRight.Y)):
                        continue

                msm_Node_Table = msm_Node_Template.copy()   
                skip = False

                msm_Node_Table["MUID"] = nodes[nodei].attrib['Knudenavn']
                # arcpy.AddMessage(dir(nodes[nodei]))
                knudekode = int(nodes[nodei].find("KnudeKode").text) if nodes[nodei].find("KnudeKode") is not None else 0

                statuskode_transformer = {0: 'Uoplyst', 1: 'I brug/drift', 2: 'Ikke I brug', 3: 'Afproppet',
                                          4: 'Opfyldt', 5: u'Død', 6: 'Projekteret / planlagt', 7: 'Anlagt',
                                          8: 'Fjernet', 50: 'Andet'}

                
                msm_Node_Table["Status"] = statuskode_transformer[int(nodes[nodei].find("StatusKode").text)] if nodes[nodei].find(
                    "StatusKode") is not None and int(nodes[nodei].find(
                    "StatusKode").text) in statuskode_transformer else None
                
                if knudekode == 3:
                    msm_Node_Table["TypeNo"] = 2
                elif knudekode == 9:
                    msm_Node_Table["TypeNo"] = 3
                else:
                    msm_Node_Table["TypeNo"] = 1
                    
                owner = re.sub(r'[^\d\w]+', "", nodes[nodei].find("Ejerfordelingsnavn").text) if nodes[nodei].find("Ejerfordelingsnavn") is not None else ""
                msm_Node_Table["Ejer"] = owner

                msm_Node_Table["Bemaerkning"] = nodes[nodei].find("Bemaerkning").text[:50] if nodes[nodei].find(
                    "Bemaerkning") is not None else ""

                if nodes[nodei].find("Broend/BroendKode") is not None:
                    broendkode = int(nodes[nodei].find("Broend/BroendKode").text)
                    msm_Node_Table["Broendkode"] = manhole_table_dict[broendkode] if broendkode in manhole_table_dict else broendkode
                    # arcpy.AddMessage((msm_Node_Table["Broendkode"], broendkode))
                # if not nodes[nodei].find()
                # msm_Node_Table["Broendkode"] = manhole_table_dict[nodes[nodei].attrib['Broendkode']]
                
                # statuskode = (int(nodes[nodei].find("KnudeKode").text) if nodes[nodei].find("KnudeKode") is not None else 0)
                # if statuskode == 8:
                #     msm_Node_Table["Description"] = u"Sløjfet"
                
                msm_Node_Table["GroundLevel"] = float(nodes[nodei].find("DaekselItems").find("Daeksel").find("Daekselkote").text) if nodes[nodei].find("DaekselItems") is not None and nodes[nodei].find("DaekselItems").find("Daeksel") is not None and nodes[nodei].find("DaekselItems").find("Daeksel").find("Daekselkote") is not None else None
                if msm_Node_Table["GroundLevel"] is None:
                    msm_Node_Table["GroundLevel"] = float(nodes[nodei].find("Terraenkote").text) if nodes[nodei].find("Terraenkote") is not None else None
                
                
                msm_Node_Table["InvertLevel"] = float(nodes[nodei].find("Bundkote").text) if nodes[nodei].find("Bundkote") is not None else None
                
                msm_Node_Table["Diameter"] = float(nodes[nodei].find("DiameterBredde").text)/1e3 if nodes[nodei].find("DiameterBredde") is not None else 1
                
                msm_Node_Table["NetTypeNo"] = int(nodes[nodei].find("TypeAfloebKode").text) if nodes[nodei].find("TypeAfloebKode") is not None else None
                
                if msm_Node_Table["NetTypeNo"] is not None and msm_Node_Table["NetTypeNo"] in range(1,5) and afloebkode[msm_Node_Table["NetTypeNo"]-1] in afloebkodeparameter:
                    if ((nodes[nodei].find("KategoriAfloebKode") is not None and int(nodes[nodei].find("KategoriAfloebKode").text) == 4)  # hvis KategoriAfloebKode er stik
                    or (nodes[nodei].find("Broend/BroendKode") is not None and int(nodes[nodei].find("Broend/BroendKode").text) == 16) # eller br�ndkode er skelbr�nd
                    or (nodes[nodei].find("Punkt") is not None and nodes[nodei].find("Punkt").find("PunktKode") is not None and int(nodes[nodei].find("Punkt").find("PunktKode").text) == 8)): # eller punktkode er tilslutningspunkt
                        if "Stikledning" in afloebkategori: # Tilf�j, hvis stikledninger skal tilf�jes
                            msm_Node_Table["Description"] = "Stikledning" if not msm_Node_Table["Description"] else msm_Node_Table["Description"] + ", Stikledning"
                        else:
                            skip = True
                    elif not nodes[nodei].find("KategoriAfloebKode") == None and int(nodes[nodei].find("KategoriAfloebKode").text) == 1: # hvis KategoriAfloebKode er hovedledning
                        msm_Node_Table["Description"] = "Hovedledning" if not msm_Node_Table["Description"] else msm_Node_Table["Description"] + ", Hovedledning"
                    elif not nodes[nodei].find("KategoriAfloebKode") == None and int(nodes[nodei].find("KategoriAfloebKode").text) == 6: # hvis KategoriAfloebKode er internt
                        if "Internt" in afloebkategori:
                            msm_Node_Table["Description"] = "Internt" if not msm_Node_Table["Description"] else msm_Node_Table["Description"] + ", Internt"
                        else:
                            skip = True
                    elif not "Ukendt" in  afloebkategori: # hvis KategoriAfloebKode ikke er stik, skelbr�nd, hovedledning eller internt (dvs. ukendt)
                        skip = True
                    else:
                        msm_Node_Table["Description"] = "Ukendt" if not msm_Node_Table["Description"] else msm_Node_Table["Description"] + ", Ukendt"
                    if not skip:
                        
                        msm_Node_Cursor.insertRow([(float(nodes[nodei].find("XKoordinat").text),float(nodes[nodei].find("YKoordinat").text))]+msm_Node_Table.values())
                        nodesDict[nodes[nodei].attrib['Knudenavn']] = (float(nodes[nodei].find("XKoordinat").text),float(nodes[nodei].find("YKoordinat").text), msm_Node_Table["InvertLevel"])
                    ID = ID+1
            except Exception as e:
                arcpy.AddWarning(traceback.format_exc())
                arcpy.AddWarning(e.message)
        del msm_Node_Cursor

        if import_catchments:
            arcpy.SetProgressor("default","Reading Catchments")
            ms_Catchment = getAvailableFilename(arcpy.env.scratchGDB + "\ms_Catchment")
            arcpy.CopyFeatures_management(os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\ms_Catchment", ms_Catchment)
            if coordinate_system:
                arcpy.DefineProjection_management(ms_Catchment, coordinate_system)
            templateFile = os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\ms_Catchment_Template"
            ms_Catchment_Fields = [a.name for a in arcpy.ListFields(templateFile) if a.name not in ["FID","Shape","SHAPE","OBJECTID"]]
            ms_Catchment_Cursor = arcpy.da.InsertCursor(ms_Catchment, ["SHAPE@"] + ms_Catchment_Fields)
            ms_Catchment_Template = OrderedDict()
            with arcpy.da.SearchCursor(templateFile,ms_Catchment_Fields) as cursor:
                for row in cursor:
                    for i in range(len(ms_Catchment_Fields)):
                        ms_Catchment_Template[ms_Catchment_Fields[i]] = row[i]
            catchmentMUIDs = []
            for nodei in range(len(nodes)):
                if nodes[nodei].attrib['Knudenavn'] in nodesDict:
                    catchmentDictionary = ms_Catchment_Template.copy()
                    try:
                        if nodes[nodei].find("DeloplandItems"):
                            deloplandrows = nodes[nodei].find("DeloplandItems").findall("Delopland")
                            for deloplandei in range(len(deloplandrows)):
                                catchmentDictionary["NodeID"] = nodes[nodei].attrib['Knudenavn']
                                i = 1
                                MUID = catchmentDictionary["NodeID"] + "c%d" % (i)
                                while MUID in catchmentMUIDs:
                                    i += 1
                                    MUID = catchmentDictionary["NodeID"] + "c%d" % (i)
                                catchmentMUIDs.append(MUID)
                                catchmentDictionary["MUID"] = MUID
                                array = arcpy.Array()
                                deloplandekoorrows = deloplandrows[deloplandei].find("DeloplandKoordItems").findall("DeloplandKoord")
                                for deloplandekoori in range(len(deloplandekoorrows)):
                                    array.add(arcpy.Point(float(deloplandekoorrows[deloplandekoori].find("Xkoordinat").text),float(deloplandekoorrows[deloplandekoori].find("Ykoordinat").text)))
                                try:
                                    catchmentDictionary["ImpArea"] = float(deloplandrows[deloplandei].find("BefaestelsePct").text)
                                except:
                                    catchmentDictionary["ImpArea"] = 0
                                try:
                                    ms_Catchment_Cursor.insertRow([arcpy.Polygon(array)] + catchmentDictionary.values())
                                except Exception as e:
                                    arcpy.AddMessage("Failed to draw catchment for node %s (%s)" % (nodes[nodei].attrib['Knudenavn'],e))
                    except Exception as e:
                        arcpy.AddMessage("Failed to create catchment for node %s (%s)" % (nodes[nodei].attrib['Knudenavn'],e))
            del ms_Catchment_Cursor
                            
        if dandas_ledninger:
            arcpy.SetProgressor("default","Reading links")    
            with open(dandas_ledninger,"r") as infile:
                txt = infile.read()
                txt = re.sub(r' xmlns="[^"]+"',"",txt)
            
            tree = ET.fromstring(txt)
            
            links = tree.findall("Ledning")
            msm_Link = getAvailableFilename(arcpy.env.scratchGDB + "\msm_Link")
            arcpy.CopyFeatures_management(os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\msm_Link", msm_Link)
            if coordinate_system:
                arcpy.DefineProjection_management(msm_Link, coordinate_system) 
            msm_Link_Template = OrderedDict()
            templateFile = os.path.dirname(os.path.realpath(__file__)) + "\Data\Templates.mdb\msm_Link_Template"
            fields = [a.name for a in arcpy.ListFields(templateFile) if a.name not in [u"FID",u"SHAPE",u"OBJECTID","SHAPE_Length"]]

            with arcpy.da.SearchCursor(templateFile,fields) as cursor:
                for row in cursor:
                    for i in range(len(fields)):
                        msm_Link_Template[fields[i]] = row[i]
            msm_Link_Cursor = arcpy.da.InsertCursor(msm_Link, ["SHAPE@"] + fields) # 
            ID = 0
            
            materialList = {15:'Ceramics',
                            16:'Ceramics',
                            1:'Concrete (Normal)',
                            14:'Concrete (Normal)',
                            13:'Iron (cast)',
                            12:'Iron (wrought)',
                            4:'Plastic',
                            5:'Plastic',
                            8:'GAP',
                            9:'Plastic',
                            10:'Plastic',
                            18:'Plastic',
                            19:'Plastic',
                            20:'Plastic',
                            21:'Plastic',
                            24:'Plastic',
                            17:'Stone',
                            0:'Unknown',
                            50:'Unknown',
                            99:'Unknown'}
            
            MUIDsUsed = []
            
            arcpy.SetProgressor("step","Reading Links (%d)" % (len(links)), 0, len(links), 1)
            for linki, link in enumerate(links):
                arcpy.SetProgressorPosition(linki) 
                if link.attrib["OpstroemKnudenavn"] in nodesDict and link.attrib["NedstroemKnudenavn"] in nodesDict:
                    linkDictionary = msm_Link_Template.copy()
                    linkDictionary["FROMNODE"] = link.attrib["OpstroemKnudenavn"]
                    linkDictionary["TONODE"] = link.attrib["NedstroemKnudenavn"]


                    statuskode_transformer = {0: 'Uoplyst', 1: 'I brug/drift', 2: 'Ikke I brug', 3: 'Afproppet', 4: 'Opfyldt', 5: u'Død', 6: 'Projekteret / planlagt', 7: 'Anlagt', 8: 'Fjernet', 50: 'Andet'}
            
                    linkDictionary["Status"] = statuskode_transformer[int(link.find("StatusKode").text)] if link.find("StatusKode") is not None and int(link.find("StatusKode").text) in statuskode_transformer else None
                    #arcpy.AddMessage("Link %s-%s"% (link.attrib["OpstroemKnudenavn"],link.attrib["NedstroemKnudenavn"]))
                    
                    i = 1
                    MUID = linkDictionary["FROMNODE"] + "_" + linkDictionary["TONODE"] + "l%d" % (i)
                    while MUID in MUIDsUsed:
                        i += 1
                        MUID = linkDictionary["FROMNODE"] + "_" + linkDictionary["TONODE"] + "l%d" % (i)
                    MUIDsUsed.append(MUID)
                    linkDictionary["MUID"] = MUID[:40]
                    
                    linkDictionary["NetTypeNo"] = int(link.find("TypeAfloebKode").text) if link.find("TypeAfloebKode") is not None else 1

                    if link.find(
                            "DatoEtableret") is not None:
                        linkDictionary["YearEst"] = int(link.find("DatoEtableret").text[:4])

                    if link.find("DelLedningItems") is not None and link.find("DelLedningItems").find("DelLedning") is not None:
                        link_delledning = link.find("DelLedningItems").find("DelLedning")

                        if link_delledning.find("TvaersnitKode") is not None:
                            k = int(link_delledning.find("TvaersnitKode").text)
                            crosssection_catalogue = {0: 1, 1: 1, 5: 2, 6: 2, 8: 2, 9: 2, 10: 2, 50: 2, 99: 2, 3: 3, 4: 3, 7: 4, 2: 5}
                            linkDictionary["TypeNo"] = crosssection_catalogue[k] if k in crosssection_catalogue else 1


                        if link_delledning.find("DiameterIndv") is not None:
                            linkDictionary["Diameter"] = float(link_delledning.find("DiameterIndv").text)/1.0e3
                        elif link_delledning.find("Handelsmaal") is not None:
                            linkDictionary["Diameter"] = float(link_delledning.find("Handelsmaal").text)/1.0e3

                        if link_delledning.find("Handelsmaal") is not None:
                            linkDictionary["Handmaal"] = float(link_delledning.find("Handelsmaal").text) / 1.0e3

                        if not linkDictionary["YearEst"] and link_delledning.find("DatoEtableret") is not None:
                            linkDictionary["YearEst"] = int(link_delledning.find("DatoEtableret").text[:4])

                        linkDictionary["Description"] = link.find("Bemaerkning").text[:50] if link.find(
                            "Bemaerkning") is not None else ""


                        if link_delledning.find("MaterialeKode") is not None:
                            material = int(link_delledning.find("MaterialeKode").text)
                            linkDictionary["MaterialID"] = materialList[material] if material in materialList else 'Concrete (Normal)'
                        else:
                            linkDictionary["MaterialID"] = 'Concrete (Normal)'

                        if link_delledning.find("BundloebskoteOpst") is not None:
                            if (link_delledning.find("DeltaKoteOpst") is not None and 
                                not link_delledning.find("DeltaKoteOpst").text == "0.00"):
                                if abs(nodesDict[linkDictionary["FROMNODE"]][2] - float(link_delledning.find("BundloebskoteOpst").text))>ignore_elevation_difference:
                                    linkDictionary["UpLevel"] = float(link_delledning.find("BundloebskoteOpst").text)
                            
                        if link_delledning.find("BundloebskoteNedst") is not None:
                            if (link_delledning.find("DeltaKoteNedst") is not None and 
                                    not link_delledning.find("DeltaKoteNedst").text == "0.00"):
                                if abs(nodesDict[linkDictionary["TONODE"]][2] - float(link_delledning.find("BundloebskoteNedst").text))>ignore_elevation_difference:
                                    linkDictionary["DwLevel"] = float(link_delledning.find("BundloebskoteNedst").text)
                        
                        if link_delledning.find("Fald") is not None:
                            linkDictionary["Slope_C"] = float(link_delledning.find("Fald").text)/10

                    try:
                        vertices = [arcpy.Point(nodesDict[linkDictionary["FROMNODE"]][0],nodesDict[linkDictionary["FROMNODE"]][1])]
                        
                        try:
                            for vertex in link.find("DelLedningItems").find("DelLedning").find("KnaekpunktItems"):
                                vertices.append(arcpy.Point(float(vertex.find("XKoordinat").text),float(vertex.find("YKoordinat").text)))
                        except:
                            pass
                        vertices.append(arcpy.Point(nodesDict[linkDictionary["TONODE"]][0],nodesDict[linkDictionary["TONODE"]][1]))
                        polyline = arcpy.Polyline(arcpy.Array(vertices))
                        msm_Link_Cursor.insertRow([polyline]+linkDictionary.values())
                        del linkDictionary["FROMNODE"]
                        del linkDictionary["TONODE"]
                    except Exception as e:
                        arcpy.AddWarning("Error on link %s-%s" % (link.attrib["OpstroemKnudenavn"], link.attrib["NedstroemKnudenavn"]))
                        arcpy.AddWarning(traceback.format_exc())
                        pass
            del msm_Link_Cursor
        
        empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\Dandas Import.lyr")
        empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
        empty_group_layer = arcpy.mapping.ListLayers(mxd, "Dandas Import", df)[0]
        empty_group_layer.name = os.path.splitext(os.path.basename(dandas_knuder))[0]
        if dandas_ledninger:
            empty_group_layer.name = empty_group_layer.name + " " + os.path.splitext(os.path.basename(dandas_ledninger))[0]
        # msm_Node_reprojected = arcpy.management.Project(msm_Node, 
                                                        # getAvailableFilename(arcpy.env.scratchGDB + "\msm_Node_reprojected"),
                                                        # arcpy.SpatialReference("ETRS 1989 UTM Zone 32N"),
                                                        # in_coor_system = coordinate_system)[0] if coordinate_system else msm_Node
                                                        
        addLayer = arcpy.mapping.Layer(msm_Node)
        arcpy.mapping.AddLayerToGroup(df, empty_group_layer, addLayer, "TOP")
        updatelayer = arcpy.mapping.ListLayers(mxd, msm_Node.split("\\")[-1], df)[0]
        
        sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Node.lyr")
        arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))
        
        if dandas_ledninger:
            addLayer = arcpy.mapping.Layer(msm_Link)
            arcpy.mapping.AddLayerToGroup(df, empty_group_layer, addLayer, "BOTTOM")   
            updatelayer = arcpy.mapping.ListLayers(mxd, msm_Link.split("\\")[-1], df)[0]
            sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\msm_Link_DDS.lyr")
            arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
            updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))
           
                    
            links_with_faulty_uplevel = []
            links_with_faulty_dwlevel = []
            with arcpy.da.UpdateCursor(msm_Link, ["MUID","UpLevel","DwLevel", "FROMNODE", "TONODE"]) as cursor:
                for row in cursor:
                    if row[1] and nodesDict[row[3]][2] > row[1]:
                        links_with_faulty_uplevel.append(row[0])
                    if row[2] and nodesDict[row[4]][2] > row[2]:
                        links_with_faulty_dwlevel.append(row[0])
            arcpy.AddMessage(links_with_faulty_uplevel)
            arcpy.AddMessage(links_with_faulty_dwlevel)
        
        if import_catchments:
            addLayer = arcpy.mapping.Layer(ms_Catchment)
            arcpy.mapping.AddLayerToGroup(df, empty_group_layer, addLayer, "BOTTOM")
            updatelayer = arcpy.mapping.ListLayers(mxd, ms_Catchment.split("\\")[-1], df)[0]
            sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\MOUSE Catchments.lyr")
            arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
            updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))
        
        return
        
class DDS2Tilbudsliste(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "a) Convert to Tender Document"
        self.description = "a) Convert to Tender Document"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        
        manholes_shapefile = arcpy.Parameter(
            displayName="Manholes:",
            name="manholes_shapefile",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        manholes_shapefile.filter.list=["Point"]
        
        pipes_shapefile = arcpy.Parameter(
            displayName="Pipes:",
            name="pipes_shapefile",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        pipes_shapefile.filter.list=["Polyline"]
        
        write_depth = arcpy.Parameter(
            displayName="Include depth for manholes?",
            name="write_depth",
            datatype="Boolean",
            category="Additional Settings",
            parameterType="Optional",
            direction="Input")
        write_depth.value = True
        
        nodes_order = arcpy.Parameter(
            displayName="Write pipes as [From Manhole]-[To Manhole] or [To Manhole]-[From Manhole]",
            name="nodes_order",
            datatype="String",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        nodes_order.filter.list = [r"[Downstream Manhole]-[Upstream Manhole]", r"[Upstream Manhole]-[Downstream Manhole]"]
        nodes_order.value = r"[Downstream Manhole]-[Upstream Manhole]"
        
        pipes_nodes_order = arcpy.Parameter(
            displayName="List all manholes first and then all pipes, or manhole->pipe->manhole",
            name="pipes_nodes_order",
            datatype="String",
            category="Additional Settings",
            parameterType="Required",
            direction="Input")
        pipes_nodes_order.filter.list = [r"All manholes first and then all pipes", r"[Manhole]-[Pipe]-[Manhole]"]
        pipes_nodes_order.value = r"All manholes first and then all pipes"
        
        params = [manholes_shapefile, pipes_shapefile, write_depth, nodes_order, pipes_nodes_order]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        return

    def updateMessages(self, parameters):
        # if parameters[0].value:  
            # with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    # parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
        # else:  
            # parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
        manhole_shapefile = parameters[0].ValueAsText
        all_manholes_shapefile = arcpy.Describe(manhole_shapefile).catalogPath
        pipe_shapefile = parameters[1].ValueAsText
        write_depth = parameters[2].Value
        nodes_order = parameters[3].Value        
        pipes_nodes_order = parameters[4].Value
        
        class Manhole:
            def __init__(self, muid, invertlevel = None, groundlevel = None, diameter = None, lineno = -1, formatted_line = ""):
                self.muid = muid
                self.invertlevel = invertlevel
                self.groundlevel = groundlevel
                self.diameter = diameter
                self.lineno = lineno
                self.formatted_line = muid
                
            def getDepth(self):
                if self.invertlevel and self.groundlevel and write_depth:
                    return ", dybde %1.2f m" % (self.groundlevel - self.invertlevel)
                else:
                    return ""

            def depth(self):
                if self.invertlevel and self.groundlevel:
                    return self.groundlevel - self.invertlevel
        class Pipe:
            def __init__(self, muid, length, fromnode, tonode, diameter, material):
                self.muid = muid
                self.length = length
                self.fromnode = fromnode
                self.tonode = tonode
                self.diameter = diameter if diameter else 0
                self.lineno = -1
                self.formatted_line = muid
                self.material = material

            @property
            def outer_diameter(self):
                if self.material == "Plastic":
                    diameter = self.diameter
                    diameters = {0.1024: 0.11, 0.149: 0.16, 0.188: 0.2, 0.233: 0.25, 0.276: 0.3, 0.377: 0.4, 0.493: 0.5,
                                 0.588: 0.6,
                                 0.781: 0.8, 0.885: 0.9, 0.985: 1, 1.085: 1.1, 1.185: 1.2, 1.285: 1.3, 1.385: 1.4,
                                 1.485: 1.5,
                                 1.585: 1.6}
                    diameters_by_index = [outer_diameter for outer_diameter in diameters.values()]
                    if diameter in diameters:
                        outer_diameter = diameters[diameter]
                    else:
                        dist = [abs(inner_diameter - diameter) for inner_diameter in diameters]
                        outer_diameter = diameters_by_index[dist.index(min(dist))]
                else:
                    outer_diameter = self.diameter
                return outer_diameter

            def depth(self, nodes):
                return np.mean([nodes[self.tonode].depth(), nodes[self.fromnode].depth()])

            def getPipeName(self):
                if nodes_order == r"[Upstream Manhole]-[Downstream Manhole]":
                    return "%s-%s" % (pipe.fromnode, pipe.tonode)
                else:
                    return "%s-%s" % (pipe.tonode, pipe.fromnode)

        is_sqlite = True if ".sqlite" in arcpy.Describe(manhole_shapefile).catalogPath else False
        manholes = OrderedDict()
        with arcpy.da.SearchCursor(manhole_shapefile, ["MUID", "InvertLevel", "GroundLevel", "Diameter"], where_clause = "Description <> 'Stikledning' OR Description IS NULL") as cursor:
            for row in cursor:
                manholes[row[0]] = Manhole(row[0],row[1],row[2],row[3])
                
        all_manholes = OrderedDict()
        with arcpy.da.SearchCursor(arcpy.Describe(all_manholes_shapefile).catalogPath, ["MUID", "InvertLevel", "GroundLevel", "Diameter"], where_clause = "Description <> 'Stikledning' OR Description IS NULL") as cursor:
            for row in cursor:
                all_manholes[row[0]] = Manhole(row[0],row[1],row[2],row[3])
               
        pipes = {}
        with arcpy.da.SearchCursor(pipe_shapefile, ["MUID", "SHAPE@LENGTH", "FROMNODEID" if is_sqlite else "FROMNODE", "TONODEID" if is_sqlite else "TONODE", "Diameter", "MaterialID"]) as cursor:
            for row in cursor:
                pipes[row[0]] = Pipe(row[0],row[1],row[2],row[3], row[4], row[5])
        pipes = OrderedDict(sorted(pipes.items(), key = lambda item: item[1].diameter))
        
        for manhole in manholes.values():
            manhole.formatted_line = u"Broend %s, oe%1.0f mm plast%s" % (manhole.muid, manhole.diameter*1e3, manhole.getDepth())
            
        for pipe in pipes.values():
            # if pipe.fromnode in all_manholes and pipe.tonode in all_manholes:
            arcpy.AddMessage(pipe.depth(all_manholes))
            # pipe.formatted_line = u"\t\tStr. %s\toe%1.0f mm %s\tlbm\t%1.0f\t\t\t" % (pipe.getPipeName(), pipe.outer_diameter*1e3, "pl" if pipe.material == "Plastic" else "bt", pipe.length)
            pipe.formatted_line = u"Str. %s\t%1.0f\t%s\t%1.0f\t%1.1f" % (
                    pipe.getPipeName(), pipe.outer_diameter * 1e3, "pl" if pipe.material == "Plastic" else "bt", pipe.length, pipe.depth(all_manholes))
        
        # Each manhole and pipe is assigned a lineno, which determines the order that the object should appear in the tender documents

        def getMaxLineno():
            return max(np.max([manhole.lineno for manhole in manholes.values()]), np.max([pipe.lineno for pipe in pipes.values()]))

        arcpy.AddMessage(pipes)

        def listParents(manhole): # Recursive loop that travels upstream and assigns an object the next lineno 
             parent_pipes = [pipe for pipe in pipes.values() if pipe.tonode == manhole.muid if pipe.fromnode in manholes]
             for parent_pipe in parent_pipes:
                 parent_pipe.lineno = getMaxLineno() + 1
                 manholes[parent_pipe.fromnode].lineno = getMaxLineno() + 1
                 listParents(manholes[parent_pipe.fromnode])
        
        # Either start travelling from a manhole that has no pipes downstream
        # for manhole in manholes.values():
        #     if ([pipe.muid for pipe in pipes.values() if pipe.tonode == manhole.muid]  and not [pipe.muid for pipe in pipes.values() if pipe.fromnode == manhole.muid]):
        #         manhole.lineno = getMaxLineno() + 1
        #         listParents(manhole)
        
        # or start travelling from a pipe ending at a pipe that is not selected in ArcMap
        for pipe in pipes.values():
            if pipe.tonode not in manholes and pipe.fromnode in manholes and pipe.lineno == -1:
                pipe.lineno = getMaxLineno() + 1
                manholes[pipe.fromnode].lineno = getMaxLineno() + 1
                try:
                    listParents(manholes[pipe.fromnode])
                except:
                    arcpy.AddMessage(pipe.fromnode)
                    arcpy.AddMessage(pipe.tonode)
                    arcpy.AddMessage(pipe.muid)

        tender_text = u"bob" # all the text that should be copied/exported
        if pipes_nodes_order == r"[Manhole]-[Pipe]-[Manhole]":
            for lineno in range(getMaxLineno()+1):
                for line in [pipe.formatted_line for pipe in pipes.values() if pipe.lineno == lineno]:
                    tender_text += u"\n" + line
                for line in [manhole.formatted_line for manhole in manholes.values() if manhole.lineno == lineno]:
                    tender_text += u"\n" + line
        else:
            for lineno in range(getMaxLineno()+1):
                for line in [manhole.formatted_line for manhole in manholes.values() if manhole.lineno == lineno]:
                    tender_text += u"\n" + line
            for lineno in range(getMaxLineno()+1):
                for line in [pipe.formatted_line for pipe in pipes.values() if pipe.lineno == lineno]:
                    tender_text += u"\n" + line

        tender_text = tender_text[1:]
        arcpy.AddMessage(tender_text)
        
        # copying the tender text to the clipboard using the library tkinter
        # try:
            # from Tkinter import Tk
        # except ImportError:
            # from tkinter import Tk
        # r = Tk()
        # r.withdraw()
        # r.clipboard_clear()
        import subprocess
        def copy2clip(txt):
            # Convert da text to UTF-8 bytes (since Python 2.7 uses bytes natively)
            if isinstance(txt, unicode):
                txt = txt.encode('utf-8')
            # Pipe da raw text directly into clip
            process = subprocess.Popen(
                'clip',
                stdin=subprocess.PIPE,
                shell=True
            )
            process.communicate(input=txt)  # Send da UTF-8 encoded bytes

        # Example usage (with Unicode string)
        # tender_text = u"Dette er en test med æ, ø, og å"
        copy2clip(tender_text)
        arcpy.AddMessage(tender_text)

        # copy2clip("ø")
        # r.clipboard_append(tender_text)
        # r.update()
        # r.destroy()
        
        
        return        

class CopyMikeUrbanFeatures(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "2) Copy Mike Urban Features to Mike Urban Database"
        self.description = "2) Copy Mike Urban Features to Mike Urban Database"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        
        MU_database = arcpy.Parameter(
            displayName="Mike Urban Database to copy features to:",
            name="MU_database",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        MU_database.filter.list = ["mdb", "sqlite"]
        
        msm_Node = arcpy.Parameter(
            displayName="Manhole Layers:",
            name="msm_Node",
            datatype="GPFeatureLayer",
            multiValue = True,
            parameterType="Optional",
            direction="Input")
        
        msm_Link = arcpy.Parameter(
            displayName="Link Layers:",
            name="msm_Link",
            datatype="GPFeatureLayer",
            multiValue = True,
            parameterType="Optional",
            direction="Input")
        
        ms_Catchment = arcpy.Parameter(
            displayName="Catchment Layers:",
            name="ms_Catchment",
            datatype="GPFeatureLayer",
            multiValue = True,
            parameterType="Optional",
            direction="Input")
        
        params = [MU_database, msm_Node, msm_Link, ms_Catchment]
        
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        if not parameters[1].value and not parameters[2].value and not parameters[3].value:
            mxd = arcpy.mapping.MapDocument("CURRENT")

            nodes = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Point'
                    and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
            if nodes:
                parameters[1].value = nodes

            links = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polyline'
                    and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
            if links:
                parameters[2].value = links

            catchments = [lyr.longName for lyr in arcpy.mapping.ListLayers(mxd) if lyr.getSelectionSet() and arcpy.Describe(lyr).shapeType == 'Polygon'
                        and "muid" in [field.name.lower() for field in arcpy.ListFields(lyr)]]
            if catchments:
                parameters[3].value = catchments

        if not parameters[0].valueAsText:
            mxd = arcpy.mapping.MapDocument("CURRENT")
            database = [lyr.dataSource for lyr in arcpy.mapping.ListLayers(mxd) if not lyr.getSelectionSet() and lyr.isFeatureLayer and ".mdb" in lyr.dataSource]
            if database:
                database = re.findall(r"(.+)(?=\.mdb)", database[0])[0] + ".mdb"
                parameters[0].value = database

        # if parameters[0].valueAsText:
        #     MU_database = parameters[0].valueAsText
        #     is_sqlite = True if ".sqlite" in MU_database else False
        #     if is_sqlite:
        #         parameters[3].enabled = False
            
        return

    def updateMessages(self, parameters):
        # if parameters[0].value:  
            # with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    # parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
        # else:  
            # parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
        MU_database = parameters[0].valueAsText
        msm_Nodes = parameters[1].values
        msm_Links = parameters[2].values
        ms_Catchments = parameters[3].values
        # arcpy.AddMessage((msm_Links, ".sqlite" in arcpy.Describe(msm_Links[0]).catalogPath))
        if (msm_Nodes and ".sqlite" in arcpy.Describe(msm_Nodes[0]).catalogPath) or (msm_Links and ".sqlite" in arcpy.Describe(msm_Links[0]).catalogPath) or (ms_Catchments and ".sqlite" in arcpy.Describe(ms_Catchments[0]).catalogPath):
            source_type = 'MUPlusDB'
        else:
            source_type = 'Geodatabase'

        is_sqlite = True if ".sqlite" in MU_database else False
        arcpy.env.overwriteOutput = True

        import random
        nodes_count = 0
        if msm_Nodes:
            for i, msm_Node in enumerate(msm_Nodes):
                selected = arcpy.Select_analysis(msm_Node, "in_memory\msm_Node_%d" % (i))
                nodes_in_msm_Node = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                for node in nodes_in_msm_Node:
                    nodes_count += 1

        link_count = 0
        if msm_Links:
            for i, msm_Link in enumerate(msm_Links):
                selected = arcpy.Select_analysis(msm_Link, "in_memory\msm_Link_%d" % (i))
                links_in_msm_Link = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                for link in links_in_msm_Link:
                    link_count += 1

        catchment_count = 0
        if ms_Catchments:
            for i, ms_Catchment in enumerate(ms_Catchments):
                selected = arcpy.Select_analysis(ms_Catchment, "in_memory\ms_Catchment_%d" % (i))
                catchments_in_ms_Catchment = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                for catchment in catchments_in_ms_Catchment:
                    catchment_count += 1

        if is_sqlite:
            arcpy.AddMessage("SQLITE database. Using XML file in MIKE+ for import.")
            with open(os.path.dirname(os.path.realpath(__file__)) + r"\Data\XML\DDS_to_MIKE+.xml", 'r') as f:
                xml_txt = f.readlines()

            if not msm_Nodes:
                enabled_lineno = [i for i, line in enumerate(xml_txt) if 'property="msm_Node" value=' in line][0]
                xml_txt[enabled_lineno] = re.sub('value *= *"[^"]+"', 'value="False"', xml_txt[enabled_lineno])

            if not msm_Links:
                enabled_lineno = [i for i, line in enumerate(xml_txt) if 'property="msm_Link" value=' in line][0]
                xml_txt[enabled_lineno] = re.sub('value *= *"[^"]+"', 'value="False"', xml_txt[enabled_lineno])

        if pythonaddins.MessageBox("You are copying %d manholes, %d pipes, and %d catchments. Continue?" % (nodes_count, link_count, catchment_count),
                                   "Confirm copy?", 1) == "OK":
            if msm_Nodes:
                for i, msm_Node in enumerate(msm_Nodes):
                    reference_MU_database = os.path.dirname(arcpy.Describe(msm_Node).catalogPath)
                    nodes_in_database = [row[0] for row in arcpy.da.SearchCursor(MU_database + "\msm_Node", ["MUID"])]
                    selected = arcpy.Select_analysis(msm_Node, "in_memory\msm_Node_%d" % (i+1))
                    nodes_in_msm_Node = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                    duplicate_nodes = np.intersect1d(nodes_in_database, nodes_in_msm_Node)
                    if duplicate_nodes.size > 0:
                        response = pythonaddins.MessageBox("The following manholes are already in the Mike Urban Database: %s" % ", ".join(duplicate_nodes)
                                                           + " Would you like to remove the existing manholes first? (No: The tool will add those duplicate manholes anyway."
                                                           + " Cancel: Skip those manholes)",
                                                           "Remove duplicate features?", 3)
                        if is_sqlite and response == "Yes":
                            arcpy.AddError("Function not supported for .sqlite database")
                            return
                        else:
                            if response == "Yes":
                                edit = arcpy.da.Editor(MU_database)
                                edit.startEditing(False, True)
                                edit.startOperation()
                                with arcpy.da.UpdateCursor(MU_database + "\msm_Node", ["MUID"], where_clause = "MUID IN ('%s')" % "', '".join(duplicate_nodes)) as cursor:
                                    for row in cursor:
                                        cursor.deleteRow()
                                edit.stopOperation()
                                edit.stopEditing(True)
                            elif response == "Cancel":
                                selected = arcpy.Select_analysis(selected, "in_memory\msm_Node_%d_filtered" % (i), where_clause = "MUID NOT IN ('%s')" % "', '".join(duplicate_nodes))

                    if is_sqlite:
                        MUIDs = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                        # sql_expression = "ATTACH DATABASE %s AS source; SELECT * INTO main.msm_Node FROM source.msm_Node WHERE MUID IN ('%s')" % (reference_MU_database, "', '".join(MUIDs))

                        source_lineno = \
                            [i for i, line in enumerate(xml_txt) if '<JobPropertyValue property="Source"' in line][0]
                        xml_txt[source_lineno] = re.sub('value *= *"[^"]+"', 'value="%s"' % reference_MU_database,
                                                        xml_txt[source_lineno])

                        msm_Node_where_clause_lineno = [i for i, line in enumerate(xml_txt) if '"msm_Node_where_clause"' in line][0]
                        xml_txt[msm_Node_where_clause_lineno] = re.sub('value *= *"[^"]+"', "value=\"MUID IN ('%s')\"" % "', '".join(MUIDs),
                                                                       xml_txt[msm_Node_where_clause_lineno])

                        # Change source type if MIKE+
                        source_type_lineno = [i for i, line in enumerate(xml_txt) if 'sourceTypeParameter' in line][0]

                        xml_txt[source_type_lineno] = xml_txt[source_type_lineno].replace("sourceTypeParameter",
                                                                                                source_type)

                        # arcpy.AddMessage(sql_expression)
                    else:
                        arcpy.Append_management(selected, MU_database + "\msm_Node", schema_type = "NO_TEST")

            if msm_Links:
                links_in_database = [row[0] for row in arcpy.da.SearchCursor(MU_database + "\msm_Link", ["MUID"])]
                for i, msm_Link in enumerate(msm_Links):
                    reference_MU_database = os.path.dirname(arcpy.Describe(msm_Link).catalogPath)
                    selected = arcpy.Select_analysis(msm_Link, "in_memory\msm_Link_%d" % (i))
                    links_in_msm_Link = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                    duplicate_links = np.intersect1d(links_in_database, links_in_msm_Link)
                    if duplicate_links.size > 0:
                        response = pythonaddins.MessageBox("The following pipes are already in the Mike Urban Database: %s" % ", ".join(duplicate_links)
                                                           + " Would you like to remove the existing pipes first? (No: The tool will add those duplicate pipes anyway."
                                                           + " Cancel: Skip those pipes)",
                                                           "Remove duplicate features?", 3)
                        if response == "Yes":
                            if is_sqlite:
                                arcpy.AddError("Function not supported for .sqlite database")
                                return
                            edit = arcpy.da.Editor(MU_database)
                            edit.startEditing(False, True)
                            edit.startOperation()
                            with arcpy.da.UpdateCursor(MU_database + "\msm_Link", ["MUID"], where_clause = "MUID IN ('%s')" % "', '".join(duplicate_links)) as cursor:
                                for row in cursor:
                                    cursor.deleteRow()
                            edit.stopOperation()
                            edit.stopEditing(True)
                        elif response == "Cancel":
                            selected = arcpy.Select_analysis(selected, "in_memory\msm_Link_%d_filtered" % (i), where_clause = "MUID NOT IN ('%s')" % "', '".join(duplicate_links))
                    if is_sqlite:
                        MUIDs = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                        source_lineno = [i for i,line in enumerate(xml_txt) if '<JobPropertyValue property="Source"' in line][0]
                        xml_txt[source_lineno] = re.sub('value *= *"[^"]+"', 'value="%s"' % reference_MU_database, xml_txt[source_lineno])

                        msm_Link_where_clause_lineno = \
                            [i for i, line in enumerate(xml_txt) if '"msm_Link_where_clause"' in line][0]
                        xml_txt[msm_Link_where_clause_lineno] = re.sub('value *= *"[^"]+"',
                                                                       "value=\"MUID IN ('%s')\"" % "', '".join(MUIDs),
                                                                       xml_txt[msm_Link_where_clause_lineno])

                        fields = [field.name.upper() for field in arcpy.ListFields(selected)]

                        if "FROMNODEID" in fields:
                            arcpy.AddMessage("BOB")
                            FROMNODEID_where_clause_lineno = \
                                [i for i, line in enumerate(xml_txt) if 'property="FromNode' in line][0]
                            xml_txt[FROMNODEID_where_clause_lineno] = xml_txt[FROMNODEID_where_clause_lineno].replace('value="FROMNODE"','value="FromNodeID"')

                            TONODEID_where_clause_lineno = \
                                [i for i, line in enumerate(xml_txt) if 'property="ToNode' in line][0]
                            xml_txt[TONODEID_where_clause_lineno] = xml_txt[TONODEID_where_clause_lineno].replace('value="TONODE"','value="ToNodeID"')

                        try:
                            source_type_lineno = [i for i, line in enumerate(xml_txt) if 'sourceTypeParameter' in line][0]

                            xml_txt[source_type_lineno] = xml_txt[source_type_lineno].replace("sourceTypeParameter",
                                                                                          source_type)
                        except Exception as e:
                            pass
                    else:
                        try:
                            # arcpy.AddMessage([field.name.lower() for field in arcpy.ListFields(MU_database + "\msm_Link")])
                            if not "fromnode" in [field.name.lower() for field in arcpy.ListFields(MU_database + "\msm_Link")]:
                                arcpy.DeleteField_management(selected, ["FROMNODE","TONODE"])
                            elif "fromnode" in [field.name.lower() for field in arcpy.ListFields(MU_database + "\msm_Link")] and not "fromnode" in [field.name.lower() for field in arcpy.ListFields(selected)]:
                                for field in ["FROMNODE","TONODE"]:
                                    arcpy.AddField_management(selected, field, "TEXT", field_length = 50, field_is_nullable="NULLABLE")
                                # pythonaddins.MessageBox("Attempting to add FROMNODE and TONODE to %s. Close the Mike Urban model before proceeding." % (MU_database + "\msm_Link"), "Close Mike Urban", 0)
                                # for field in ["FROMNODE", "TONODE"]:
                                # arcpy.AddField_management(MU_database + "\msm_Link", field, "TEXT", field_length = 255)
                            arcpy.Append_management(selected, MU_database + "\msm_Link", schema_type = "NO_TEST")
                        except Exception as e:

                            arcpy.AddError("FromNode and ToNode not found in Mike Urban Database. Try manually copying the Pipes")
                            arcpy.AddError(traceback.format_exc())
                            arcpy.AddError(e)
                            # raise(e)


            if ms_Catchments:
                arcpy.AddMessage(is_sqlite)
                if is_sqlite:
                    i = 0
                    arcpy.AddMessage(ms_Catchments)
                    reference_MU_database = os.path.dirname(arcpy.Describe(ms_Catchments[0]).catalogPath)
                    arcpy.AddMessage(reference_MU_database)
                    source_lineno = \
                    [i for i, line in enumerate(xml_txt) if '<JobPropertyValue property="Source"' in line][0]
                    xml_txt[source_lineno] = re.sub('value *= *"[^"]+"', 'value="%s"' % reference_MU_database,
                                                    xml_txt[source_lineno])


                    if len([i for i, line in enumerate(xml_txt) if 'sourceTypeParameter' in line])>0:
                        source_type_lineno = [i for i, line in enumerate(xml_txt) if 'sourceTypeParameter' in line][0]

                        xml_txt[source_type_lineno] = xml_txt[source_type_lineno].replace("sourceTypeParameter",
                                                                                          source_type)

                    selected = arcpy.Select_analysis(ms_Catchments[0], "in_memory\msm_Catchment_%d" % (i))
                    MUIDs = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                    msm_Catchment_where_clause_lineno = \
                        [i for i, line in enumerate(xml_txt) if '"msm_Catchment_where_clause"' in line][0]
                    xml_txt[msm_Catchment_where_clause_lineno] = re.sub('value *= *"[^"]+"',
                                                                   "value=\"MUID IN ('%s')\"" % "', '".join(MUIDs),
                                                                   xml_txt[msm_Catchment_where_clause_lineno])

                    msm_Catchcon_where_clause_lineno = \
                        [i for i, line in enumerate(xml_txt) if '"msm_Catchcon_where_clause"' in line][0]
                    xml_txt[msm_Catchcon_where_clause_lineno] = re.sub('value *= *"[^"]+"',
                                                                   "value=\"CatchID IN ('%s')\"" % "', '".join(MUIDs),
                                                                   xml_txt[msm_Catchcon_where_clause_lineno])

                    xml_txt = [line.replace(r'property="msm_Catchment" value="False"','property="msm_Catchment" value="True"')  for line in xml_txt]
                    xml_txt = [line.replace(r'property="msm_CatchCon" value="False"',
                                            'property="msm_CatchCon" value="True"') for line in xml_txt]
                else:
                    # arcpy.AddMessage(ms_Catchments)
                    for i, ms_Catchment in enumerate(ms_Catchments):
                        # arcpy.AddMessage(type(ms_Catchment))
                        arcpy.AddMessage(arcpy.Describe(ms_Catchment).catalogPath)
                        if not ".mdb" in arcpy.Describe(ms_Catchment).catalogPath:
                            selected = arcpy.Select_analysis(ms_Catchment, "in_memory\ms_Catchment_%d" % (i))

                            ms_CatchmentMUIDs = [row[0] for row in arcpy.da.SearchCursor(selected,"MUID")]
                            duplicateMUIDs = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_Catchment"),"MUID",where_clause = "MUID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
                            # duplicateHModA = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_HModA"),"CatchID",where_clause = "CatchID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]
                            duplicateCatchCon = [row[0] for row in arcpy.da.SearchCursor(os.path.join(MU_database,"msm_CatchCon"),"CatchID",where_clause = "CatchID IN ('%s')" % ("', '".join(ms_CatchmentMUIDs)))]

                            errorMessage = ""
                            if duplicateMUIDs:
                                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Layer in Mike Urban Database" % ("(', '".join(duplicateMUIDs))
                            # if duplicateHModA:
                            #     errorMessage = errorMessage + "\n" if errorMessage else ""
                            #     errorMessage += "Catchments with MUID ('%s') already exist in Model Records (msm_HModA) in Mike Urban Database" % ("(', '".join(duplicateHModA))
                            if duplicateCatchCon:
                                errorMessage = errorMessage + "\n" if errorMessage else ""
                                errorMessage += "Catchments with MUID ('%s') already exist in Catchment Connections (msm_CatchCon) in Mike Urban Database" % ("(', '".join(duplicateCatchCon))
                            if errorMessage:
                                arcpy.AddWarning(errorMessage)
                            if not errorMessage or pythonaddins.MessageBox("%s\nTransfer catchments anyway?" % (errorMessage), "Confirm Transfer", 1) == "OK":
                                arcpy.Append_management(selected, os.path.join(MU_database,"ms_Catchment"), schema_type = "NO_TEST")
                                fields = [field.name.lower() for field in arcpy.ListFields(selected)]
                                def readField(field):
                                    if field.lower() in fields:
                                        return field
                                    else:
                                        return "SHAPE@XY"
                                with arcpy.da.SearchCursor(selected,["MUID","ImpArea","NodeID", readField("ParAID"), readField("RedFactor"), readField("ConcTime"), readField("InitLoss")]) as catchmentCursor:
                                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_CatchCon"),["CatchID","NodeID","TypeNo"]) as cursor:
                                        for row in catchmentCursor:
                                            nID = 0 if not row[2] else row[2]
                                            cursor.insertRow((row[0],nID,1))
                                    catchmentCursor.reset()
                                    with arcpy.da.InsertCursor(os.path.join(MU_database,"msm_HModA"),["CatchID","ImpArea","ParAID","LocalNo","ConcTime","RFactor","ILoss","CoeffNo","TACoeff"]) as cursor:
                                        for row in catchmentCursor:
                                            iArea = 0 if not row[1] else row[1]
                                            cursor.insertRow((row[0],iArea,"-DEFAULT-" if not row[3] or type(row[3]) is tuple else row[3],
                                                              0,
                                                              7 if not row[5] or type(row[5]) is tuple else row[5],
                                                              0.9 if not row[4] or type(row[4]) is tuple else row[4],
                                                              0.0006 if not row[6] or type(row[6]) is tuple else row[6],
                                                              0,0.33))
                        else:
                            input_database = arcpy.Describe(ms_Catchment).catalogPath.split(".mdb")[0] + ".mdb"
                            selected = arcpy.Select_analysis(ms_Catchment, "in_memory\ms_Catchment_%d" % (i))
                            arcpy.Append_management(selected, os.path.join(MU_database,"ms_Catchment"),"NO_TEST")
                            MUIDs = [row[0] for row in arcpy.da.SearchCursor(selected, ["MUID"])]
                            # arcpy.management.Append(selected, MU_database + "\ms_Catchment")
                            arcpy.AddMessage( "MUID IN ('%s')" % "', '".join(MUIDs))
                            selected_HModA = arcpy.TableSelect_analysis(os.path.join(input_database, "msm_HModA"), "in_memory\msm_HModA", where_clause = "CatchID IN ('%s')" % "', '".join(MUIDs))[0]
                            # fields = [field.name for field in arcpy.ListFields(MU_database + "\msm_HModA") if not field.name == "OBJECTID"]
                            # arcpy.AddMessage(fields)
                            # with arcpy.da.InsertCursor(MU_database + "\msm_HModA", ["CatchID","ImpArea","ParAID","LocalNo","ConcTime","RFactor","ILoss","CoeffNo","TACoeff"]) as target_cursor:
                            # with arcpy.da.SearchCursor(os.path.join(input_database, "msm_HModA"), fields, where_clause = "CatchID IN ('%s')" % "', '".join(MUIDs)) as reference_cursor:
                            # for row in reference_cursor:
                            # arcpy.AddMessage(row)
                            # target_cursor.insertRow(("gay",0,"-DEFAULT-",0,7,0.9,0.0006,0,0.33))
                            # # arcpy.AddMessage(selected_HModA)
                            # arcpy.CopyFeatures_management(selected_HModA, "K:\Hydrauliske modeller\Papirkurv\SonsOfKemet")
                            arcpy.management.Append(selected_HModA, MU_database + "\msm_HModA")
                            selected_CatchCon = arcpy.TableSelect_analysis(os.path.join(input_database, "msm_CatchCon"), "in_memory\msm_CatchCon", where_clause = "CatchID IN ('%s')" % "', '".join(MUIDs))[0]
                            arcpy.AddMessage("CatchID IN ('%s')" % "', '".join(MUIDs))
                            arcpy.AddMessage([row[0] for row in arcpy.da.SearchCursor(selected_CatchCon, ["CatchID"])])
                            arcpy.management.Append(selected_CatchCon, MU_database + "\msm_CatchCon")
            if is_sqlite:
                import tempfile
                output_path = os.path.join(tempfile.gettempdir(), '%s.xml' % os.path.basename(reference_MU_database).split(".")[0])
                with open(output_path,'w') as f:
                    f.writelines(xml_txt)
                import subprocess
                subprocess.call('explorer /select,"%s"' % output_path)
        return


class GPSManholes(object):
    def __init__(self):
        self.label       = "b) Find GPS-measured manholes from Dandas XML"
        self.description = "b) Find GPS-measured manholes from Dandas XML"

    def getParameterInfo(self):
        #Define parameter definitions
        
        knuderparameter = arcpy.Parameter(
            displayName="Dandas Knuder XML-file:",
            name="dandas_knuder",
            datatype="File",
            parameterType="Required",
            direction="Input")
        knuderparameter.filter.list=["xml"]
                    
        afloebkodeparameter = arcpy.Parameter(
            displayName="Include only following sewer type:",
            name="afloebkodeparameter",
            datatype="GPString",
            parameterType="Required",
            multiValue=True,
            direction="Input")
        afloebkodeparameter.filter.type = "ValueList"  
        afloebkodeparameter.filter.list = ["Waste Water", "Storm Drain", "Combined Sewer"]
        
        afloebkategori = arcpy.Parameter(
            displayName="Include also the following:",
            name="afloebkategori",
            datatype="GPString",
            parameterType="Optional",
            multiValue = True,
            direction="Input")
        afloebkategori.filter.type = "ValueList"  
        afloebkategori.filter.list = ["Stikledninger", "Ukendt"]
        
        outputfil = arcpy.Parameter(
            displayName="Output shape-file:",
            name="outputfil",
            datatype="DEShapefile",
            parameterType="Optional",
            direction="Ouptut")
                    
        params = [knuderparameter, afloebkodeparameter, afloebkategori, outputfil] 
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        return

    def updateMessages(self, parameters):
        # if parameters[0].value:  
            # with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    # parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
        # else:  
            # parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
    
        def unicode_csv_reader(utf8_data, dialect=csv.excel, **kwargs):
            csv_reader = csv.reader(utf8_data, dialect=dialect, **kwargs)
            for row in csv_reader:
                yield [unicode(cell, 'utf-8') for cell in row]
        import os 
        dir_path = os.path.dirname(os.path.realpath(__file__))
        reader = unicode_csv_reader(open(dir_path + r"\Oprindelse.txt"),delimiter = ";")
        journal = []
        journalkode = []
        for row in reader:
            journal.append(row[1])
            journalkode.append(row[2])

        
        
        with open(parameters[0].ValueAsText,"r") as infile:
            txt = infile.read()
            txt = re.sub(r' xmlns="[^"]+"',"",txt)
        
        tree = ET.fromstring(txt)
        
        nodes = tree.findall("Knude")
        
        outGDB = arcpy.env.scratchGDB
        env.workspace = outGDB
        outFC = arcpy.CreateFeatureclass_management(outGDB, "Knuder", "Point")
        arcpy.AddField_management(outFC,"Knudenavn","TEXT")
        arcpy.AddField_management(outFC,"Afloebstype","TEXT")
        arcpy.AddField_management(outFC,"Ledningstype","TEXT")
        arcpy.AddField_management(outFC,"Bundkote","FLOAT")
        arcpy.AddField_management(outFC,"Daekselkote","FLOAT")
        arcpy.AddField_management(outFC,"Terraenkote","FLOAT")
        arcpy.AddField_management(outFC,"BKOprind","TEXT")
        arcpy.AddField_management(outFC,"DKOprind","TEXT")
        arcpy.AddField_management(outFC,"TKOprind","TEXT")
        cur = arcpy.da.InsertCursor(outFC, ["SHAPE@XY","KnudeNavn","Afloebstype","Ledningstype","Bundkote","Daekselkote","Terraenkote","BKOprind","TKOprind","DKOprind"])
        ID = 0
        knudenavn = []
        
        oprindelseKoord = {
            "": "Skoen",
            "0": "Uoplyst",
            "1": "Skoennet",
            "2": "Projekt",
            "3": "Digitaliseret fra kort",
            "4": "Fotogrammetri",
            "5": "Landmaaling",
            "6": "TV-inspektion",
            "7": "Broendrapport",
            "8": "Nedstik",
            "20": "Broendfoto GP",
            "50": "Andet",
            "60": "Projektkontrol"
            }
        
        afloebkode = ["Waste Water", "Storm Drain", "Combined Sewer"]
        for nodei in range(len(nodes)):
            try:
                skip = False
                # arcpy.AddMessage(nodes[nodei].attrib['Knudenavn'])
                oprindelsTKJournal = nodes[nodei].find("OprTerraenkoteJournalnr")
                oprindelseBKJournal = nodes[nodei].find("OprindBundkoteJournalnr")
                try:
                    oprindelseDKJournal = nodes[nodei].find("DaekselItems").find("Daeksel").find("OprindKoteJournalnr")
                except:
                    oprindelseDKJournal = None
                    
                if oprindelseBKJournal != None:
                    oprindelseBKJournal = journalkode[[i for i,a in enumerate(journal) if a == oprindelseBKJournal.text][0]]
                    oprindelseBKJournal = oprindelseKoord[oprindelseBKJournal]
                if oprindelsTKJournal != None:
                    oprindelsTKJournal = journalkode[[i for i,a in enumerate(journal) if a == oprindelsTKJournal.text][0]]
                    oprindelsTKJournal = oprindelseKoord[oprindelsTKJournal]
                if oprindelseDKJournal != None:
                    oprindelseDKJournal = journalkode[[i for i,a in enumerate(journal) if a == oprindelseDKJournal.text][0]]
                    oprindelseDKJournal = oprindelseKoord[oprindelseDKJournal]
                        
                    
                try:
                    daekselkote = float(nodes[nodei].find("DaekselItems").find("Daeksel").find("Daekselkote").text)
                except:
                    daekselkote = None
                try:
                    bundkote = float(nodes[nodei].find("Bundkote").text)
                except:
                    bundkote = None
                try:
                    terraenkote = float(nodes[nodei].find("Terraenkote").text)
                except:
                    terraenkote = None
                if not nodes[nodei].find("TypeAfloebKode") == None and int(nodes[nodei].find("TypeAfloebKode").text) in range(1,4) and afloebkode[int(nodes[nodei].find("TypeAfloebKode").text)-1] in parameters[1].ValueAsText:
                    if (not nodes[nodei].find("KategoriAfloebKode") == None and int(nodes[nodei].find("KategoriAfloebKode").text) == 4) or (not nodes[nodei].find("Broend/BroendKode") == None and int(nodes[nodei].find("Broend/BroendKode").text) == 16):
                        if "Stikledning" in parameters[2].ValueAsText:
                            KategoriAfloebKode = "Stikledning"
                        else:
                            skip = True
                    elif not nodes[nodei].find("KategoriAfloebKode") == None and int(nodes[nodei].find("KategoriAfloebKode").text) == 1:
                        KategoriAfloebKode = "Hovedledning"
                    elif not nodes[nodei].find("KategoriAfloebKode") == None and int(nodes[nodei].find("KategoriAfloebKode").text) == 6:
                        skip = True
                    elif not "Ukendt" in  parameters[2].ValueAsText:
                        skip = True 
                    elif not nodes[nodei].find("KategoriAfloebKode") == None:
                        KategoriAfloebKode = nodes[nodei].find("KategoriAfloebKode").text
                    else:
                        KategoriAfloebKode = "Ukendt"
                    if not skip:
                        cur.insertRow([(float(nodes[nodei].find("XKoordinat").text),float(nodes[nodei].find("YKoordinat").text)),nodes[nodei].attrib['Knudenavn'],
                            afloebkode[int(nodes[nodei].find("TypeAfloebKode").text)-1],KategoriAfloebKode,bundkote,daekselkote,terraenkote,oprindelseBKJournal,oprindelseDKJournal,oprindelsTKJournal])
                    ID = ID+1
            except Exception as e:
                arcpy.AddMessage("Knude: %s\nError: %s (line %d)" % (nodes[nodei].attrib['Knudenavn'],str(e),sys.exc_info()[2].tb_lineno))
                
        del cur
        
        # cur = arcpy.da.UpdateCursor(outFC, ["Knudenavn"])
        # ID = 0
    
        # for row in cur:
            # row[0] = knudenavn[ID]
            # cur.updateRow(row)
            # ID = ID+1
        # del cur
        if parameters[3].ValueAsText == None:            
            outFC = "Knuder"
            addLayer = arcpy.mapping.Layer(outFC)
            mxd = arcpy.mapping.MapDocument("CURRENT")
            df = arcpy.mapping.ListDataFrames(mxd)[0] 
            arcpy.mapping.AddLayer(df, addLayer)
            arcpy.RefreshActiveView()
        else:
            arcpy.CopyFeatures_management ("Knuder", parameters[3].ValueAsText)
        
        return
        
class DDS2ArcGISRK(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "c) Find unconnected Dandas Catchments w/ Reduktionsknuder"
        self.description = "c) Find unconnected Dandas Catchments w/ Reduktionsknuder"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        
        knuderparameter = arcpy.Parameter(
            displayName="Dandas Knuder XML-file",
            name="dandas_knuder",
            datatype="File",
            parameterType="Required",
            direction="Input")
        knuderparameter.filter.list=["xml"]
        
        reduktionsknudeparameter = arcpy.Parameter(
            displayName="Reduktionsknuder",
            name="reduktionsknuder",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
                    
        omraaderparameter = arcpy.Parameter(
            displayName="Reduktionsomraade",
            name="reduktionsomraade",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        omraaderparameter.filter.type = "ValueList"  
        omraaderparameter.filter.list = []
                    
        params = [knuderparameter, reduktionsknudeparameter, omraaderparameter] 
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        return

    def updateMessages(self, parameters):
        if parameters[1].value:  
            try:
                with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
            except:
                pass
        else:  
            parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
        with open(parameters[0].ValueAsText,"r") as infile:
            txt = infile.read()
            txt = re.sub(r' xmlns="[^"]+"',"",txt)
        
        tree = ET.fromstring(txt)
        
        nodes = tree.findall("Knude")
        
        outGDB = arcpy.env.scratchGDB
        env.workspace = outGDB
        outFC = arcpy.CreateFeatureclass_management(outGDB, "Catchments", "Polygon")
        arcpy.AddField_management(outFC,"Knudenavn","TEXT")
        arcpy.AddField_management(outFC,"Forbundet","SHORT")
        cur = arcpy.da.InsertCursor(outFC, ["SHAPE@"])
        
        ID = 0
        knudenavn = []
        for nodei in range(len(nodes)):
            deloplandrows = nodes[nodei].find("DeloplandItems").findall("Delopland")
            for deloplandei in range(len(deloplandrows)):
                array = arcpy.Array()
                deloplandekoorrows = deloplandrows[deloplandei].find("DeloplandKoordItems").findall("DeloplandKoord")
                for deloplandekoori in range(len(deloplandekoorrows)):
                    array.add(arcpy.Point(float(deloplandekoorrows[deloplandekoori].find("Xkoordinat").text),float(deloplandekoorrows[deloplandekoori].find("Ykoordinat").text),ID))
                try:
                    cur.insertRow([arcpy.Polygon(array)])
                    ID = ID+1
                    knudenavn.append(nodes[nodei].attrib['Knudenavn'])
                except:
                    arcpy.AddMessage("Process failed to draw catchment")
        del cur
        
        cur = arcpy.da.UpdateCursor(outFC, ["Knudenavn"])
        ID = 0
    
        for row in cur:
            row[0] = knudenavn[ID]
            cur.updateRow(row)
            ID = ID+1
        del cur
        
        outFC = "Catchments"
        addLayer = arcpy.mapping.Layer(outFC)
        mxd = arcpy.mapping.MapDocument("CURRENT")  
        df = arcpy.mapping.ListDataFrames(mxd)[0]  
        arcpy.mapping.AddLayer(df, addLayer)
        arcpy.RefreshActiveView()
        
        #cur = arcpy.da.SearchCursor(parameters[1].ValueAsText, ["ReduktionsOmraade","NodeID"],"""[ReduktionsOmraade] = %s""" % (parameters[2].ValueAsText))
        where_clause = """[ReduktionOmraade] = '%s'""" % (parameters[2].ValueAsText)
        arcpy.AddMessage(where_clause)
        nodes = []
        with arcpy.da.SearchCursor(parameters[1].ValueAsText, ["ReduktionOmraade","NodeID"],where_clause=where_clause) as rows:
            for row in rows:
                nodes.append(row[1])
        
        where_clause="""Knudenavn IN (%s)""" % ("'" + "', '".join(nodes) + "'")
        arcpy.AddMessage(where_clause)
        cur = arcpy.da.UpdateCursor(outFC, ["Knudenavn","Forbundet"],where_clause)
        for row in cur:
            row[1] = 1
            cur.updateRow(row)
        del cur
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="Catchments", selection_type="NEW_SELECTION", where_clause = where_clause)

        return
        
class DDS2ArcGISMU(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "d) Find unconnected Dandas Catchments w/ Mike Urban Manholes"
        self.description = "d) Find unconnected Dandas Catchments w/ Mike Urban Manholes"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions
        
        knuderparameter = arcpy.Parameter(
            displayName="Dandas Knuder XML-file",
            name="dandas_knuder",
            datatype="File",
            parameterType="Required",
            direction="Input")
        knuderparameter.filter.list=["xml"]
        
        reduktionsknudeparameter = arcpy.Parameter(
            displayName="Manholes",
            name="reduktionsknuder",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
                    
        knudeparameter = arcpy.Parameter(
            displayName="Field for MUID",
            name="knudeparameter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        knudeparameter.filter.type = "ValueList"  
        knudeparameter.filter.list = []
                    
        params = [knuderparameter, reduktionsknudeparameter, knudeparameter] 
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):        
        # Set the default distance threshold to 1/100 of the larger of the width
        #  or height of the extent of the input features.  Do not set if there is no 
        #  input dataset yet, or the user has set a specific distance (Altered is true).
        #
        if parameters[1].altered:
            fields = arcpy.ListFields(parameters[1].value)
            fil = []
            for field in fields:
                fil.append(field.name)
                if field.name == "NodeID" or field.name == "MUID":
                     parameters[2].value = field.name

            parameters[2].filter.type = "ValueList"
            parameters[2].filter.list = fil
        return

    def updateMessages(self, parameters):
        if parameters[1].value:  
            try:
                with arcpy.da.SearchCursor(parameters[1].valueAsText, field_names="ReduktionOmraade") as rows:  
                    parameters[2].filter.list = sorted(list(set([row[0] for row in rows])))  
            except:
                pass
        else:  
            parameters[2].filter.list = []  
        return

    def execute(self, parameters, messages):
        with open(parameters[0].ValueAsText,"r") as infile:
            txt = infile.read()
            txt = re.sub(r' xmlns="[^"]+"',"",txt)
        
        tree = ET.fromstring(txt)
        
        nodes = tree.findall("Knude")
        
        outGDB = arcpy.env.scratchGDB
        env.workspace = outGDB
        outFC = arcpy.CreateFeatureclass_management(outGDB, "Catchments", "Polygon")
        arcpy.AddField_management(outFC,"Knudenavn","TEXT")
        arcpy.AddField_management(outFC,"Forbundet","SHORT")
        cur = arcpy.da.InsertCursor(outFC, ["SHAPE@"])
        
        ID = 0
        knudenavn = []
        for nodei in range(len(nodes)):
            deloplandrows = nodes[nodei].find("DeloplandItems").findall("Delopland")
            for deloplandei in range(len(deloplandrows)):
                array = arcpy.Array()
                deloplandekoorrows = deloplandrows[deloplandei].find("DeloplandKoordItems").findall("DeloplandKoord")
                for deloplandekoori in range(len(deloplandekoorrows)):
                    array.add(arcpy.Point(float(deloplandekoorrows[deloplandekoori].find("Xkoordinat").text),float(deloplandekoorrows[deloplandekoori].find("Ykoordinat").text),ID))
                try:
                    cur.insertRow([arcpy.Polygon(array)])
                    ID = ID+1
                    knudenavn.append(nodes[nodei].attrib['Knudenavn'])
                except:
                    arcpy.AddMessage("Process failed to draw catchment for node %s" % (nodes[nodei].attrib['Knudenavn']))
        del cur
        
        cur = arcpy.da.UpdateCursor(outFC, ["Knudenavn"])
        ID = 0
    
        for row in cur:
            row[0] = knudenavn[ID]
            cur.updateRow(row)
            ID = ID+1
        del cur
        
        
        outFC = "Catchments"
        addLayer = arcpy.mapping.Layer(outFC)
        mxd = arcpy.mapping.MapDocument("CURRENT")  
        df = arcpy.mapping.ListDataFrames(mxd)[0]  
        arcpy.mapping.AddLayer(df, addLayer)
        arcpy.RefreshActiveView()
        
        #cur = arcpy.da.SearchCursor(parameters[1].ValueAsText, ["ReduktionsOmraade","NodeID"],"""[ReduktionsOmraade] = %s""" % (parameters[2].ValueAsText))
        nodes = []
        with arcpy.da.SearchCursor(parameters[1].ValueAsText, [parameters[2].ValueAsText]) as rows:
            for row in rows:
                nodes.append(row[0])
        
        where_clause="""Knudenavn IN (%s)""" % ("'" + "', '".join(nodes) + "'")
        arcpy.AddMessage(where_clause)
        cur = arcpy.da.UpdateCursor(outFC, ["Knudenavn","Forbundet"],where_clause)
        for row in cur:
            row[1] = 1
            cur.updateRow(row)
        del cur
        arcpy.SelectLayerByAttribute_management(in_layer_or_view="Catchments", selection_type="NEW_SELECTION", where_clause = where_clause)

        return