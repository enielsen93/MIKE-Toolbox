"""
Created on Mon Jul 30 11:21:31 2018

@author: eni
"""

# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import arcpy
import numpy as np
import csv
import os
import xlwt
import re

from arcpy import env

def getAvailableFilename(filepath, parent = None):
    parent = os.path.basename(re.sub(r"\.[^\.\\]+$","", parent)).replace(".","_").replace("-","_").replace(" ","_").replace(",","_") if parent else None
    filepath = "%s\%s_%s" % (os.path.dirname(filepath), parent, os.path.basename(filepath)) if parent else filepath
    if arcpy.Exists(filepath):
        try:
            arcpy.Delete_management(filepath)
            return filepath
        except:
            i = 1
            while arcpy.Exists(filepath + "%d" % i):
                try:
                    arcpy.Delete_management(filepath + "%d" % i)
                    return filepath + "%d" % i
                except:
                    i += 1
            return filepath + "%d" % i
    else: 
        return filepath
        
class Toolbox(object):
    def __init__(self):
        self.label =  "(Obsolete) Compare Mike Urban Database with Spildevandsplan"
        self.alias  = "(Obsolete) Compare Mike Urban Database with Spildevandsplan"

        # List of tool classes associated with this toolbox
        self.tools = [GetImperviousnessSPV, SelectSPV, AnalyzeSPV] 

class SelectSPV(object):
    def __init__(self):
        self.label       = "Select SPV-catchments that overlap with MU catchments"
        self.description = "Select SPV-catchments that overlap with MU catchments"
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
        
        SPV_feature = arcpy.Parameter(
            displayName="Spildevandsplan (polygons)",
            name="SPV_feature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        SPV_feature.filter.list = ["Polygon"]
        
        percentageCoverage = arcpy.Parameter(
            displayName= "Percentage of SPV-catchment covered by MU catchments",  
            name="percentageCoverage",  
            datatype="GPDouble",  
            parameterType="Optional",  
            direction="Input") 
        percentageCoverage.value = 10
        
        parameters = [catchments, SPV_feature, percentageCoverage]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):        
        catchments = parameters[0].ValueAsText
        SPV_feature = parameters[1].ValueAsText
        percentageCoverage = float(parameters[2].ValueAsText)
        
        arcpy.Dissolve_management(in_features=SPV_feature, out_feature_class="SPV_feature_diss", dissolve_field="UserId", statistics_fields="", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")
        SPV_featureDiss = "SPV_feature_diss"
        
        SPV_intersect = arcpy.Intersect_analysis(in_features=[[catchments, 1],[SPV_featureDiss, 2]], out_feature_class="SPV_intersect", join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")
        
        SPV_catchments = set()
        with arcpy.da.SearchCursor(SPV_intersect,['UserId','SHAPE@AREA']) as intersectCursor:
            for intersectRow in intersectCursor:
                SPV_catchments.add(intersectRow[0])
        
        where_clause=r"""UserId IN (%s)""" % ("'" + "', '".join(SPV_catchments) + "'")
        
        SPV_catchments_select = set()
        with arcpy.da.SearchCursor(SPV_featureDiss,['UserId','SHAPE@AREA'],where_clause= where_clause) as spvCursor:
            for spvRow in spvCursor:
                area = spvRow[1]
                areaCatchments = 0
                with arcpy.da.SearchCursor(SPV_intersect,['UserId','SHAPE@AREA'], where_clause = "UserId = '%s'" % spvRow[0]) as intersectCursor:
                    for intersectRow in intersectCursor:
                        areaCatchments += intersectRow[1]
                if areaCatchments > area * percentageCoverage / 100:
                    SPV_catchments_select.add(spvRow[0])
        #arcpy.AddMessage(SPV_catchments_select)
        where_clause=r"""UserId IN (%s)""" % ("'" + "', '".join(SPV_catchments_select) + "'")
        arcpy.AddMessage(where_clause)
        arcpy.SelectLayerByAttribute_management(in_layer_or_view=SPV_feature, selection_type="NEW_SELECTION", where_clause = where_clause)
        where_clause=r"""UserId NOT IN (%s)""" % ("'" + "', '".join(SPV_catchments_select) + "'")
        arcpy.FeatureClassToFeatureClass_conversion (SPV_intersect, arcpy.env.scratchGDB , "Catchments_outside_SPV", where_clause = where_clause)
        addLayer = arcpy.mapping.Layer(arcpy.env.scratchGDB + "\Catchments_outside_SPV")
        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.mapping.AddLayer(df, addLayer)
        
        return
        
class GetImperviousnessSPV(object):
    def __init__(self):
        self.label       = "1) Compare Catchments with Spildevandsplan"
        self.description = "1) Compare Catchments with Spildevandsplan"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        CatchmentLayer = arcpy.Parameter(
            displayName="Catchment Layer (from Mike Urban Database)",
            name="database",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        CatchmentLayer.filter.list = ["Polygon"]
            
        SPV_feature = arcpy.Parameter(
            displayName="Spildevandsplan (polygons)",
            name="SPV_feature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        SPV_feature.filter.list = ["Polygon"]
            
        # includeWasteWater = arcpy.Parameter(
            # displayName="Include wastewater?",
            # name="includeWasteWater",
            # datatype="Boolean",
            # parameterType="Optional",
            # direction="Input")
        
        showWarnings = arcpy.Parameter(
            displayName="Show catchments that may not fit spildevandsplan",
            name="showWarnings",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
            
        OPL_name = arcpy.Parameter(
			displayName= "Field with SPV-opland name",  
			name="OPL_name",  
			datatype="GPString",  
			parameterType="Required",  
			direction="Input") 
                    
        parameters = [CatchmentLayer, SPV_feature, showWarnings, OPL_name]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[1].ValueAsText and not parameters[3].ValueAsText:
            parameters[3].filter.list = [f.name for f in arcpy.Describe(parameters[1].ValueAsText).fields]
            
        # if parameters[1].ValueAsText and not parameters[3].ValueAsText:
            # # parameters[3].filter.list = [f.name for f in arcpy.Describe(parameters[1].ValueAsText).fields]
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):
        CatchmentLayer = parameters[0].ValueAsText
        if ".sqlite" in arcpy.Describe(CatchmentLayer).catalogPath or ".mdb" in arcpy.Describe(CatchmentLayer).catalogPath:
            MU_database = re.findall(r"(.+(?:\.(mdb)|(sqlite)))", arcpy.Describe(CatchmentLayer).catalogPath)[0][0]
        else:
            MU_database = "not_available"
        SPV_feature = parameters[1].ValueAsText
        # includeWasteWater = parameters[2].ValueAsText
        showWarnings = parameters[2].Value
        OPL_name = parameters[3].ValueAsText

        mxd = arcpy.mapping.MapDocument("CURRENT")
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        arcpy.env.addOutputsToMap = False
        ms_Catchment = arcpy.CopyFeatures_management(CatchmentLayer, getAvailableFilename(arcpy.env.scratchGDB + "\ms_CatchmentImp", parent = MU_database))[0] if MU_database != "not_available" else CatchmentLayer
        # if not includeWasteWater:
        # with arcpy.da.UpdateCursor(ms_Catchment[0], ['NetTypeNo']) as cursor:
        # for row in cursor:
        # if row[0] == 3:
        # cursor.deleteRow()


        SPV_feature = arcpy.CopyFeatures_management(SPV_feature, getAvailableFilename(arcpy.env.scratchGDB + "\SPV_feature", parent = MU_database))[0]
        ms_CatchmentImpLayer = arcpy.MakeFeatureLayer_management(ms_Catchment, "ms_CatchmentImpLayer")

        # arcpy.JoinField_management(in_data=ms_Catchment, in_field="MUID", join_table=MU_database + r"\msm_HModA", join_field="CatchID", fields="ImpArea")

        layer_add = arcpy.mapping.Layer(ms_Catchment)
        if MU_database != "not_available":
            arcpy.mapping.AddLayer(df, layer_add,"TOP")
            updatelayer = arcpy.mapping.ListLayers(mxd, layer_add, df)[0]
            sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\ms_Catchment_Symbology.lyr")
            arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
            updatelayer.replaceDataSource(unicode(layer_add.workspacePath), 'FILEGDB_WORKSPACE', unicode(layer_add.datasetName))

        try:
            IntersectFeature = arcpy.Intersect_analysis(in_features=[[layer_add, 2],[SPV_feature, 1]], out_feature_class="GetImperviousnessIntersect", join_attributes="ONLY_FID", cluster_tolerance="-1 Unknown", output_type="INPUT")
        except:
            try:
                arcpy.AddWarning("Error upon running intersect analysis - attempting a repair on catchment and spildevandsplan layer")
                arcpy.RepairGeometry_management(SPV_feature)
                arcpy.RepairGeometry_management(layer_add)
                IntersectFeature = arcpy.Intersect_analysis(in_features=[[layer_add, 2],[SPV_feature, 1]], out_feature_class="GetImperviousnessIntersect", join_attributes="ONLY_FID", cluster_tolerance="-1 Unknown", output_type="INPUT")
            except Exception as e:
                arcpy.AddError("Exporting layers to Scratch Database")
                arcpy.CopyFeatures_management(layer_add, arcpy.env.scratchGDB + r"\addlayerDebug")
                arcpy.CopyFeatures_management(SPV_feature,arcpy.env.scratchGDB + r"\SPV_featureDebug")
                arcpy.AddError(arcpy.env.scratchGDB)
                arcpy.AddError([[layer_add, 2],[SPV_feature, 1]])
                raise e

        catchmentFeatureArea = {}
        for row in arcpy.da.SearchCursor(ms_Catchment,["OBJECTID","SHAPE@AREA"]):
            catchmentFeatureArea[row[0]-1] = row[1]
        # spvFeature = arcpy.env.scratchGDB + r"\SPV_feature"

        rows = int(arcpy.GetCount_management(IntersectFeature)[0])
        catchmentFID = np.empty(rows,dtype=np.int32)
        spvFID = np.empty(rows,dtype=np.int32)
        area = np.empty(rows)
        try:
            with arcpy.da.SearchCursor(IntersectFeature, ['FID_%s' % (arcpy.Describe(layer_add).name), 'FID_%s' % (os.path.basename(str(SPV_feature))),'SHAPE@AREA']) as cursor:
                for i,row in enumerate(cursor):
                    catchmentFID[i] = row[0]-1
                    spvFID[i] = row[1]
                    area[i] = row[2]
        except:
            arcpy.AddError(IntersectFeature)
            arcpy.AddError([f.name for f in arcpy.ListFields(IntersectFeature)])
            arcpy.AddError(['FID_%s' % (arcpy.Describe(layer_add).name),'FID_%s' % (os.path.basename(str(SPV_feature))),'SHAPE@AREA'])
            arcpy.AddMessage(arcpy.Describe(layer_add))
            raise

        catchmentUnique = np.unique(catchmentFID)
        catchmentSPVFID = np.empty(len(catchmentUnique),dtype=np.int32)
        catchmentPercInside = np.empty(len(catchmentUnique),dtype=np.float32)
        catchmentPercOther = np.empty(len(catchmentUnique),dtype=np.float32)

        for FIDi,FID in enumerate(catchmentUnique):
            idx = [i for i,v in enumerate(catchmentFID) if v==FID]
            catchmentSPVFID[FIDi] = spvFID[idx[np.argmax(area[idx])]]
            catchmentPercInside[FIDi] = np.max(area[idx])/catchmentFeatureArea[FID]*1e2
            catchmentPercOther[FIDi] = (np.sum(area[idx])-np.max(area[idx]))/catchmentFeatureArea[FID]*1e2


        arcpy.AddField_management(ms_Catchment, "SPVOpl", "TEXT", field_length=50)
        arcpy.AddField_management(ms_Catchment, "PercInside", "DOUBLE")
        arcpy.AddField_management(ms_Catchment, "PercInOther", "DOUBLE")
        SPVrows = int(arcpy.GetCount_management(SPV_feature)[0])
        SPVOplande = np.empty(SPVrows,dtype='object')
        SPVArea = np.empty(SPVrows)
        with arcpy.da.SearchCursor(SPV_feature, [OPL_name,'SHAPE@AREA']) as cursor:
            for i,row in enumerate(cursor):
                SPVOplande[i] = row[0]
                SPVArea[i] = row[1]

        SPVOplandeDict = {}
        SPVOplandeAreaDict = {}
        SPVOplandePEDict = {}
        catchmentFields = ['OBJECTID','SPVOpl','SHAPE@AREA','ImpArea','Area','Persons','PercInside', 'PercInOther']
        try:
            with arcpy.da.UpdateCursor(ms_Catchment, catchmentFields) as cursor:
                for row in cursor:
                    if row[0]-1 in catchmentUnique:
                        row[1] = SPVOplande[catchmentSPVFID[np.where(catchmentUnique==row[0]-1)[0][0]]-1]
                        row[6] = round(catchmentPercInside[np.where(catchmentUnique==row[0]-1)[0][0]],2)
                        row[7] = round(catchmentPercOther[np.where(catchmentUnique==row[0]-1)[0][0]],2)
                        cursor.updateRow(row)
                        # if row[5] == None:
                        # row[5] = 0
                        # if row[4] == 0:
                        # row[3] = -1
                        # if row[1] not in SPVOplandeDict:
                        # SPVOplandeDict[row[1]] = row[2]*row[3]/1e2
                        # SPVOplandeAreaDict[row[1]] = row[2]
                        # SPVOplandePEDict[row[1]] = row[5]
                        # else:
                        # SPVOplandeDict[row[1]] = SPVOplandeDict[row[1]] + row[2]*row[3]/1e2
                        # SPVOplandeAreaDict[row[1]] = SPVOplandeAreaDict[row[1]] + row[2]
                        # SPVOplandePEDict[row[1]] = SPVOplandePEDict[row[1]] + row[5]
        except Exception as e:
            if type(e) == TypeError:
                for i in [1,2,3]:
                    if row[i]==None:
                        arcpy.AddError("%s is type None in catchment features" % catchmentFields[i])
            # if "A column was specified that does not exist" in str(e):

            arcpy.AddMessage([f.name for f in arcpy.ListFields(ms_Catchment)])
            arcpy.AddMessage(ms_Catchment)
            raise e

        # addLayer = arcpy.mapping.Layer(ms_Catchment[0])
        # arcpy.mapping.AddLayer(df, addLayer, "TOP")
        # updatelayer = arcpy.mapping.ListLayers(mxd, os.path.basename(str(ms_Catchment[0])), df)[0]
        # sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\ms_Catchment_SPVOpland.lyr")
        # arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        # updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(addLayer.datasetName))

        # arcpy.mapping.AddLayer(df, arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\ms_Catchment_SPVOpland.lyr"))
        if showWarnings:
            arcpy.AddMessage(ms_Catchment)
            layer_add = arcpy.mapping.Layer(ms_Catchment)
            arcpy.mapping.AddLayer(df, layer_add, "TOP")
            updatelayer = arcpy.mapping.ListLayers(mxd, os.path.basename(str(ms_Catchment)), df)[0]
            sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Warnings.lyr")
            arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
            updatelayer.replaceDataSource(unicode(layer_add.workspacePath), 'FILEGDB_WORKSPACE', unicode(layer_add.datasetName))
        return

class AnalyzeSPV(object):
    def __init__(self):
        self.label       = "2) Export Excel Sheet of Catchments compared with Spildevandsplan"
        self.description = "2) Export Excel Sheet of Catchments compared with Spildevandsplan"
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        catchments = arcpy.Parameter(
            displayName="Catchment Shapefile with SPV-connections",
            name="database",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        SPV_feature = arcpy.Parameter(
            displayName="Spildevandsplan (polygons)",
            name="SPV_feature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        SPV_feature.filter.list = ["Polygon"]
        
        OPL_name = arcpy.Parameter(
			displayName= "Field with SPV-opland name",  
			name="OPL_name",  
			datatype="GPString",  
			parameterType="Required",  
			direction="Input") 
        
        imp_field = arcpy.Parameter(
			displayName= "Field with imperviousness in Spildevandsplan",  
			name="imp_field",  
			datatype="GPString",  
			parameterType="Optional",
			direction="Input") 
            
        sheet = arcpy.Parameter(
            displayName="Excel-sheet with SPV-summary",
            name="Excel_sheet",
            datatype="File",
            parameterType="Required",
            direction="Output")
        sheet.filter.list = ["xls"]
        
        includeWasteWater = arcpy.Parameter(
            displayName="Include wastewater?",
            name="includeWasteWater",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
            
        getBookmarkConnections = arcpy.Parameter(
            displayName="Data Driven Pages layer",
            name="getBookmarkConnections",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        getBookmarkConnections.filter.list = ["Polygon"]
        
        parameters = [catchments, SPV_feature, OPL_name, imp_field, sheet, includeWasteWater, getBookmarkConnections]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[1].ValueAsText and not parameters[2].ValueAsText:
            parameters[2].filter.list = [f.name for f in arcpy.Describe(parameters[1].ValueAsText).fields]
        if parameters[1].ValueAsText and not parameters[3].ValueAsText:
            parameters[3].filter.list = [f.name for f in arcpy.Describe(parameters[1].ValueAsText).fields]
        return

    def updateMessages(self, parameters): #optional
       
        return

    def execute(self, parameters, messages):
        catchments = parameters[0].ValueAsText
        SPV_feature = parameters[1].ValueAsText
        OPL_name = parameters[2].ValueAsText
        imp_field = parameters[3].ValueAsText
        sheet = parameters[4].ValueAsText
        includeWasteWater = parameters[5].ValueAsText
        getBookmarkConnections = parameters[6].ValueAsText
        
        class Opland:
            def __init__(self, area, imp_area, pe_total):
                self.area = area
                self.imp_area = imp_area
                self.pe_total = pe_total        
                
        oplande_model = {}
        with arcpy.da.SearchCursor(catchments, ['SPVOpl','SHAPE@AREA','ImpArea','Area','Persons']) as cursor:
            for row in cursor:
                SPVOpl = "None" if row[0] is None else row[0]
                PE = 0 if row[4] is None else row[4]
                area = row[3]*1e4 if row[3] else row[1]
                if SPVOpl == "None":
                    arcpy.AddMessage(area)
                    
                
                if SPVOpl not in oplande_model:
                    oplande_model[SPVOpl] = Opland(area, area*row[2]/1e2, PE)
                else:
                    oplande_model[SPVOpl].area += area
                    oplande_model[SPVOpl].imp_area += area*row[2]/1e2
                    oplande_model[SPVOpl].pe_total += PE
        
        oplande_SPV = {}
        with arcpy.da.SearchCursor(SPV_feature, [OPL_name, "SHAPE@AREA", imp_field]) as cursor:
            for row in cursor:
                try:
                    oplande_SPV[row[0]] = Opland(row[1]/1e4, (float(row[2].replace(",",".")) if isinstance(row[2], (str,unicode)) else row[2])*row[1]/1e4, 0)
                except Exception as e:
                    arcpy.AddError(row)
                    raise(e)

        oplande_SPV["None"] = Opland(0,0,0)
        
        if getBookmarkConnections:
            spvBookmarkIntersect = arcpy.Intersect_analysis(in_features=[[getBookmarkConnections, 1],[SPV_feature, 2]], out_feature_class="spvBookmarkIntersect", join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")
        
        if not sheet==None:
            wb = xlwt.Workbook()
            ws = wb.add_sheet('Oplandssammenligning')
            for col in np.arange(1,9):
                ws.col(col).width = int(2962*(8/9.44))
            ws.col(0).width = int(2962*(11/9.44))
            ws.col(9).width = int(2962*(14/9.44))
            xlwt.add_palette_colour("header", 0x21)
            wb.set_colour_RGB(0x21, 0, 157, 224)
            # xlwt.add_palette_colour("column1", 0x21)
            # wb.set_colour_RGB(0x22, 228, 223, 236)
            ws.write(0,0,"Oplandsnavn",xlwt.easyxf('font: bold on, color white; pattern: pattern solid, fore_colour header; '))
            ws.write(1,0,"",xlwt.easyxf('font: bold on, color white; pattern: pattern solid, fore_colour header'))
            ws.write_merge(0,0,1,2,u"Total areal [ha]",xlwt.easyxf('align: horiz center; font: bold on, color white; pattern: pattern solid, fore_colour header'))
            ws.write_merge(0,0,3,4,u"Befæstet areal [ha]",xlwt.easyxf('align: horiz center; font: bold on, color white; pattern: pattern solid, fore_colour header'))
            ws.write_merge(0,0,5,6,u"Afløbskoefficient [%]",xlwt.easyxf('align: horiz center; font: bold on, color white; pattern: pattern solid, fore_colour header'))
            ws.write_merge(0,0,7,8,u"PE",xlwt.easyxf('align: horiz center; font: bold on, color white; pattern: pattern solid, fore_colour header'))
            ws.write(0,10,"Hyperlink til kort",xlwt.easyxf('font: bold on'))
            row1Style = "align: horiz right; font: color white; pattern: pattern solid, fore_colour header; border: bottom thick"
            ws.write(1,1,"Model",xlwt.easyxf(row1Style))
            ws.write(1,2,"SPV",xlwt.easyxf(row1Style))
            ws.write(1,3,"Model",xlwt.easyxf(row1Style))
            ws.write(1,4,"SPV",xlwt.easyxf(row1Style))
            ws.write(1,5,"Model",xlwt.easyxf(row1Style))
            ws.write(1,6,"SPV",xlwt.easyxf(row1Style))
            ws.write(1,7,"Model",xlwt.easyxf(row1Style))
            ws.write(1,8,"SPV",xlwt.easyxf(row1Style))
            i = 2

            for key,value in sorted(oplande_model.iteritems()):
                key = str(key) if not type(key) is str else key
                if key not in oplande_SPV:
                    arcpy.AddMessage(str(key))
                    arcpy.AddWarning("Can't find %s in current SPV-plan" % (key))
                else:
                    default_style = xlwt.easyxf('', num_format_str = "0.0")
                    default_style_border = xlwt.easyxf('border: right thick', num_format_str = "0.0")
                    ws.write(i,0,key,xlwt.easyxf('font: bold on, color white; pattern: pattern solid, fore_colour header; border: right thick'))
                    ws.write(i,1, oplande_model[key].area/1e4, default_style)
                    ws.write(i,2, oplande_SPV[key].area, default_style_border)
                    ws.write(i,3,oplande_model[key].imp_area/1e4, default_style)
                    ws.write(i,4, oplande_SPV[key].imp_area, default_style_border)
                    ws.write(i,5, oplande_model[key].imp_area/oplande_model[key].area*1e2, xlwt.easyxf('', num_format_str = "0"))
                    if oplande_SPV[key].area > 0:
                        ws.write(i,6,oplande_SPV[key].imp_area/oplande_SPV[key].area*1e2, xlwt.easyxf('border: right thick', num_format_str = "0"))
                    else:
                        ws.write(i,6,0, xlwt.easyxf('border: right thick'))
                    ws.write(i,7,oplande_model[key].pe_total, default_style)
                    ws.write(i,8,oplande_SPV[key].pe_total, default_style_border)
                    if getBookmarkConnections:
                        areas = {}
                        areaSum = 0
                        pageNames = set()
                        with arcpy.da.SearchCursor(spvBookmarkIntersect,['PageName','UserId','SHAPE@AREA'], where_clause = "UserId = '%s'" % key) as cursor:
                            for row in cursor:
                                if row[0] in areas:
                                    areas[row[0]] += row[2]
                                else:
                                    areas[row[0]] = row[2]
                        areaSum = sum(areas.values())
                        
                        for pagei,pageName in enumerate(sorted(areas, key=areas.get, reverse=True)):
                            hyperlink = 'HYPERLINK("Kort\%s.pdf";"%s")' % (pageName,pageName)
                            ws.write(i,10+pagei,xlwt.Formula(hyperlink),xlwt.easyxf('font: underline on'))               
                        try:
                            maxPage = max(areas,key=areas.get)
                            if len(arcpy.ListFields(spvBookmarkIntersect,"Bynavn"))>0:
                                with arcpy.da.SearchCursor(spvBookmarkIntersect,['PageName','Bynavn'], where_clause = "PageName = '%s'" % maxPage) as cursor:
                                    for row in cursor:
                                        ws.write(i,9,row[1])
                                        break
                        except ValueError:
                            arcpy.AddWarning("Warning: Catchment %s is not found in any of the data driven pages" % key) 
                    else:
                        hyperlink = 'HYPERLINK("Kort\%s.pdf";"%s")' % (key,key)
                        ws.write(i,10,xlwt.Formula(hyperlink),xlwt.easyxf('font: underline on'))
                    i = i+1
            for j in range(1, 9):
                ws.write(i,j,"",xlwt.easyxf("border: top thick"))
            if ".xls" not in sheet:
                filename = sheet + ".xls"
            else:
                filename = sheet
            wb.save(filename)
        return
        