# -*- coding: utf-8 -*-
import os
import arcpy
import numpy as np
import re
from arcpy import env
arcpy.env.addOutputsToMap = False
import sys
import time
import mikeio
import bisect
from scipy.spatial import cKDTree
from datetime import timedelta
from scipy.interpolate import RegularGridInterpolator


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
        self.label =  "Tools for Mike Urban Flood"
        self.alias  = "Tools for Mike Urban Flood"

        # List of tool classes associated with this toolbox
        self.tools = [InterpolateToMesh, DFSUFloodStatisticsToRaster, DFSUToRaster, DFSUPlumeStatistics, DFSUToPolygons, DFS2ToRaster, CutFromMesh, CutFromDFSU, MeshToTIN]

    
        
class InterpolateToMesh(object):
    def __init__(self):
        self.label       = "Interpolate DHM to DFSU (ArcGIS Pro only)"
        self.description = "Interpolate DHM to DFSU (ArcGIS Pro only)"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DHMFile = arcpy.Parameter(
            displayName="Raster of terrain elevation",
            name="DHMFile",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        
        MeshFile = arcpy.Parameter(
            displayName="Input Mesh File",
            name="MeshFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MeshFile.filter.list = ["mesh"]
        
        MeshFileOutput = arcpy.Parameter(
            displayName="Output Mesh File with elevations interpolated",
            name="MeshFileOutput",
            datatype="File",
            parameterType="Required",
            direction="Output")
        MeshFileOutput.filter.list = ["mesh"]
        
        parameters = [DHMFile, MeshFile, MeshFileOutput]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[1].ValueAsText and not parameters[2].ValueAsText:
            parameters[2].value = parameters[1].ValueAsText.replace(".mesh","_interp.mesh")
        return

    def updateMessages(self, parameters): #optional
            
        return

    def execute(self, parameters, messages):
        tic = time.time()
        arcpy.AddMessage(parameters[2].Value)
        arcpy.AddMessage(type(parameters[2].Value))
        DHMFile = parameters[0].valueAsText
        MeshFile = parameters[1].valueAsText
        MeshFileOutput = parameters[2].valueAsText

        statusUpdate("Reading Mesh",tic)
        dfs = mikeio.dfsu.Mesh(MeshFile)

        statusUpdate("Reading Mesh: Getting Node Coordinates",tic)
        data = dfs.node_coordinates

        statusUpdate("Reading Raster",tic)
        DHMRaster = np.flip(arcpy.RasterToNumPyArray(DHMFile, nodata_to_value = 0),axis=0)

        DHMRasterInfo = np.array([float(arcpy.GetRasterProperties_management(DHMFile,"LEFT")[0].replace(",",".")),
                   float(arcpy.GetRasterProperties_management(DHMFile,"BOTTOM")[0].replace(",",".")),
                   float(arcpy.GetRasterProperties_management(DHMFile,"CELLSIZEX")[0].replace(",",".")),
                   float(arcpy.GetRasterProperties_management(DHMFile,"CELLSIZEY")[0].replace(",","."))])

        raster = arcpy.Raster(DHMFile)
        lower_left_corner = arcpy.Point(raster.extent.XMin, raster.extent.YMin)
        DHMRasterXs = lower_left_corner.X + np.arange(raster.width)*raster.meanCellWidth + raster.meanCellWidth/2
        DHMRasterYs = lower_left_corner.Y + np.arange(raster.height)*raster.meanCellHeight + raster.meanCellHeight/2

        # DHMRasterXs = np.arange(DHMRasterInfo[0],DHMRasterInfo[0]+DHMRaster.shape[1]*DHMRasterInfo[2],float(DHMRasterInfo[2])).astype(np.float32)
        # DHMRasterYs = np.arange(DHMRasterInfo[1],DHMRasterInfo[1]+DHMRaster.shape[0]*DHMRasterInfo[3],float(DHMRasterInfo[3])).astype(np.float32 )

        interp = RegularGridInterpolator((DHMRasterYs+raster.meanCellHeight, DHMRasterXs+raster.meanCellWidth), DHMRaster)

        # DHMRasterX, DHMRasterY = np.meshgrid(DHMRasterXs,DHMRasterYs)
        # DHMRasterXsFlat = DHMRasterX.flatten()
        # DHMRasterYsFlat = DHMRasterY.flatten()
        # DHMRasterFlat = DHMRaster.flatten()

        # idx = np.where(DHMRasterFlat!=0)
        # DHMRasterXsFlatSparse = DHMRasterXsFlat[idx]
        # DHMRasterYsFlatSparse = DHMRasterYsFlat[idx]
        # DHMRasterFlatSparse = DHMRasterFlat[idx]

        tic = time.time()
        bedLevelFlat = np.zeros((data.shape[0]),dtype = np.float64)

        statusUpdate("Interpolating raster to mesh (nearest neighbor)",tic)
        bedLevelFlat = interp(data[:,[1,0]])
        # for row in range(int(len(bedLevelFlat))):
        #     try:
        #         bedLevelFlat[row] = interp([bisect.bisect_left(DHMRasterYs, data[row,1]),bisect.bisect_left(DHMRasterXs, data[row,0])]
        #     except Exception as e:
        #         bedLevelFlat[row] = 0

        statusUpdate("Saving mesh",tic)
        dfs.zn = bedLevelFlat
        dfs.write(MeshFileOutput)
        return

class DFSUFloodStatisticsToRaster(object):
    def __init__(self):
        self.label       = "Convert DFSU Flood Statistics to Raster (ArcGIS Pro only)"
        self.description = "Convert DFSU Flood Statistics to Raster (ArcGIS Pro only)"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DFSUFile = arcpy.Parameter(
            displayName="Input DFSU File",
            name="DFSUFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        DFSUFile.filter.list = ["dfsu"]
        
        DFSUField = arcpy.Parameter(
            displayName="DFSU field to export",
            name="DFSUField",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        
        RasterFileOutput = arcpy.Parameter(
            displayName="Output Raster File",
            name="RasterFileOutput",
            datatype="File",
            parameterType="Required",
            direction="Output")
        RasterFileOutput.filter.list = ["tif"]
        
        clip_layers = arcpy.Parameter(
            displayName="Layers to clip dfsu by (e.g. set all raster cells inside building shapefile to NoData",
            name="clip_layers",
            multiValue="true",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            category="Additional Settings",
            direction="Input")
        clip_layers.filter.list = ["Polygon"]
        
        searchDistance = arcpy.Parameter(
            displayName="Search Distance to find nearest neighbor",
            name="searchDistance",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        searchDistance.value = 6.0
        
        raster_cell_size = arcpy.Parameter(
            displayName="Output Raster Cell Size",
            name="raster_cell_size",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        raster_cell_size.value = 1.0
        
        parameters = [DFSUFile, DFSUField, RasterFileOutput, clip_layers, searchDistance, raster_cell_size]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered:
            dfs = mikeio.dfsu.Dfsu2DH(parameters[0].valueAsText)
            statusUpdate(dfs.items,1)
            parameters[1].filter.list = [item.name for item in dfs.items]
        return

    def updateMessages(self, parameters): #optional
        if parameters[0].ValueAsText and not parameters[2].ValueAsText:
            parameters[2].Value = parameters[0].ValueAsText.replace(".dfsu",".tif")
        
        return

    def execute(self, parameters, messages):
        tic = time.time()
        DFSUFile = parameters[0].valueAsText
        DFSUField = parameters[1].valueAsText
        RasterFileOutput = parameters[2].valueAsText
        clip_layers = parameters[3].valueAsText
        searchDistance = parameters[4].Value
        raster_cell_size = parameters[5].Value

        clip_shapes = []
        if clip_layers:
            clip_layers_list = clip_layers.split(";")
            statusUpdate("Reading Clip Layers")
            for clip_layer in clip_layers:
                clip_layer_dissolved = arcpy.Dissolve_management(clip_layer, os.path.join("in_memory", os.path.splitext(
                    os.path.basename(clip_layer))[0]))[0]
                with arcpy.da.SearchCursor(clip_layer_dissolved, ["SHAPE@"]) as cursor:
                    for row in cursor:
                        clip_shapes.append(row[0])

        statusUpdate("Reading DFSU file", tic)
        dfs = mikeio.dfsu.Dfsu2DH(DFSUFile)

        statusUpdate("Retrieving element coordinates from DFSU file", tic)
        element_coordinates = dfs.geometry.element_coordinates
        dfs_read = dfs.read(items=[i for i, a in enumerate(dfs.items) if DFSUField == a.name])
        dfs_read_data = dfs_read.to_numpy()
        np.nan_to_num(dfs_read_data, copy=False)
        # arcpy.AddMessage((dfs_read_data, dfs_read_data.shape))

        statusUpdate("Isolating DFSU Elements with value", tic)

        elements_with_water = np.where(dfs_read_data[0, 0, :] > 0.003)[0] if DFSUField == "Maximum water depth" else \
        np.where(dfs_read_data[0, 0, :] != 0)[0]
        no_elements = False if len(elements_with_water) > 0 else True

        if no_elements:
            arcpy.AddError("Error: Could not find any DFSU Elements with a value different from zero")
        else:
            x_limit = [np.min(element_coordinates[elements_with_water, 0]) - searchDistance,
                       np.max(element_coordinates[elements_with_water, 0]) + searchDistance]
            y_limit = [np.min(element_coordinates[elements_with_water, 1]) - searchDistance,
                       np.max(element_coordinates[elements_with_water, 1]) + searchDistance]
            # arcpy.AddMessage((x_limit, y_limit))
            raster_xs_vector = np.arange(x_limit[0], x_limit[1], raster_cell_size)
            raster_ys_vector = np.arange(y_limit[0], y_limit[1], raster_cell_size)

            raster_x, raster_y = np.meshgrid(raster_xs_vector, raster_ys_vector)
            raster_depth = np.zeros(raster_x.shape)
            raster_x_flat = raster_x.flatten()
            raster_y_flat = raster_y.flatten()
            raster_depth_flat = np.zeros(raster_x.flatten().shape + tuple([1]))

            statusUpdate("Retrieving Raster Elements near water", tic)
            elements_searched = []
            idx = set()
            for element_i, element in enumerate(element_coordinates[elements_with_water]):
                ix = np.where(np.abs(element[0] - raster_x_flat) < searchDistance)[0]
                idx.update(ix[np.where(np.abs(element[1] - raster_y_flat[ix]) < searchDistance)[0]])

            idx_remove = []
            statusUpdate("Removing raster elements that overlap clip_layers", tic)
            for i in idx:
                point = arcpy.Point(raster_x_flat[i], raster_y_flat[i])
                for clip_shape in clip_shapes:
                    if clip_shape.geometry.contains(point):
                        idx_remove.append(i)
            for i in idx_remove:
                idx.remove(i)

            statusUpdate("Removing raster elements that are not contained inside DFSU-file", tic)
            idx_array = np.array(list(idx))
            raster_coord_in_mesh = idx_array[np.where(~dfs.geometry.contains(np.column_stack((raster_x_flat[idx_array],
                                                                                     raster_y_flat[
                                                                                         idx_array]))))]  # idx_list[np.where(dfs.contains(np.column_stack((raster_x_flat[idx_list],raster_y_flat[idx_list]))))]
            for i in raster_coord_in_mesh:
                idx.remove(i)

            statusUpdate("Creating KDTree", tic)
            elements_searchable = np.where((x_limit[0] < element_coordinates[:, 0]) &
                                           (element_coordinates[:, 0] < x_limit[1]) &
                                           (y_limit[0] < element_coordinates[:, 1]) &
                                           (element_coordinates[:, 1] < y_limit[1]))[0]
            dfsu_cKDTree = cKDTree(element_coordinates[elements_searchable, 0:2])

            statusUpdate("Interpolating DFSU to Raster (nearest neighbor)", tic)
            for i in idx:
                element_i = dfsu_cKDTree.query([raster_x_flat[i], raster_y_flat[i]])[1]
                raster_depth_flat[i] = dfs_read_data[0, 0, elements_searchable[element_i]]

            statusUpdate("Saving Raster", tic)
            raster_depth = raster_depth_flat.reshape(raster_depth.shape[0:2] + tuple([1]))

            raster_depth_compressed = raster_depth
            raster = arcpy.NumPyArrayToRaster(np.flip(raster_depth_compressed[:, :, 0], axis=0),
                                              lower_left_corner=arcpy.Point(x_limit[0], y_limit[0]),
                                              x_cell_size=raster_cell_size,
                                              y_cell_size=raster_cell_size,
                                              value_to_nodata=0)
            raster.save(RasterFileOutput)
        return
        
class DFSUToRaster(object):
    def __init__(self):
        self.description = "Convert DFSU to Raster (ArcGIS Pro only)"
        self.label       = "Convert DFSU to Raster (ArcGIS Pro only)"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DFSUFile = arcpy.Parameter(
            displayName="Input DFSU File",
            name="DFSUFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        DFSUFile.filter.list = ["dfsu"]
        
        DFSUFields = arcpy.Parameter(
            displayName="DFSU field to export",
            name="DFSUFields",
            datatype="String",
            multiValue="true",
            parameterType="Optional",
            direction="Input")
            
        DFSUFieldSummary = arcpy.Parameter(
            displayName="DFSU field statistics",
            name="DFSUFieldSummary",
            datatype="String",
            parameterType="Optional",
            direction="Input")
        DFSUFieldSummary.filter.list = ["Max","Last Timestep","Specific Timestep"]
        
        DFSUFieldTimeStep = arcpy.Parameter(
            displayName="Timestep from DFSU File to export",
            name="DFSUFieldTimeStep",
            datatype="String",
            multiValue="true",
            parameterType="Optional",
            direction="Input")
        
        RasterOutputFolder = arcpy.Parameter(
            displayName="Output Folder",
            name="RasterFileOutput",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")
        
        clip_layers = arcpy.Parameter(
            displayName="Layers to clip dfsu by (e.g. set all raster cells inside building shapefile to NoData",
            name="clip_layers",
            multiValue="true",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            category="Additional Settings",
            direction="Input")
        clip_layers.filter.list = ["Polygon"]
        
        searchDistance = arcpy.Parameter(
            displayName="Search Distance to find nearest neighbor",
            name="searchDistance",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        searchDistance.value = 6.0
        
        raster_cell_size = arcpy.Parameter(
            displayName="Output Raster Cell Size",
            name="raster_cell_size",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        raster_cell_size.value = 1.0
        
        parameters = [DFSUFile, DFSUFields, DFSUFieldSummary, DFSUFieldTimeStep, RasterOutputFolder, clip_layers, searchDistance, raster_cell_size]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered:
            dfs = mikeio.dfsu.Dfsu2DH(parameters[0].valueAsText)
            items = [item.name for item in dfs.items]
            if "Surface elevation" in items and not 'Total water depth' in items:
                items.append('Total water depth (derived)')
            parameters[1].filter.list = items
            if parameters[2].ValueAsText == "Specific Timestep":
                parameters[3].filter.list = [str(dfs.start_time+timedelta(seconds=dfs.timestep*i)) for i in range(dfs.n_timesteps+1)]

        return

    def updateMessages(self, parameters): #optional
        if parameters[0].ValueAsText and not parameters[4].ValueAsText:
            parameters[4].Value = os.path.dirname(parameters[0].ValueAsText)
        
        return

    def execute(self, parameters, messages):
        tic = time.time()
        DFSUFile = parameters[0].valueAsText
        DFSUFields = parameters[1].valueAsText
        DFSUFieldSummary = parameters[2].valueAsText
        DFSUFieldTimeStep = parameters[3].valueAsText
        RasterOutputFolder = parameters[4].valueAsText
        clip_layers = parameters[5].valueAsText
        searchDistance = parameters[6].Value
        raster_cell_size = parameters[7].Value
        DFSUFields = DFSUFields.split(";")
        derived_depth = False
        for field_i in range(len(DFSUFields)):
            DFSUFields[field_i] = DFSUFields[field_i].replace("'","")
            if DFSUFields[field_i] == "Total water depth (derived)":
                DFSUFields[field_i] = "Surface elevation"
                derived_depth = True
                derived_depth_field_i = field_i

        clip_shapes = []
        if clip_layers:
            clip_layers_list = clip_layers.split(";")
            statusUpdate("Reading Clip Layers")
            for clip_layer in clip_layers:
                clip_layer_dissolved = arcpy.Dissolve_management(clip_layer, os.path.join("in_memory",os.path.splitext(os.path.basename(clip_layer))[0]))[0]
                with arcpy.da.SearchCursor(clip_layer_dissolved, ["SHAPE@"]) as cursor:
                    for row in cursor:
                        clip_shapes.append(row[0])

        statusUpdate("Reading DFSU file", tic)
        dfs = mikeio.dfsu.Dfsu2DH(parameters[0].valueAsText)

        statusUpdate("Retrieving element coordinates from DFSU file", tic)
        element_coordinates = dfs.element_coordinates

        dfs_read = dfs.read(items=[i for i,item in enumerate(dfs.items) if item.name in DFSUFields])
        dfs_read_data = dfs_read.to_numpy()
        if derived_depth:
            dfs_read_data[derived_depth_field_i] = dfs_read_data[derived_depth_field_i] - element_coordinates[:, -1]

        for i in range(len(dfs_read_data)):
            np.nan_to_num(dfs_read_data[i], copy = False)
        #     dfs_read_data[derived_depth_field_i, :, :] = dfs_read_data[derived_depth_field_i, :, :] - element_coordinates[:, -1]

        statusUpdate("Isolating DFSU Elements with value",tic)
        elements_with_water = np.where(np.max(np.abs(dfs_read_data[0]),axis=0) > 0.01)[0] if DFSUFields[0] == "Total water depth" else np.where(np.max(np.abs(dfs_read_data[0]),axis=0))[0]
        no_elements = False if len(elements_with_water)>0 else True

        if no_elements:
            arcpy.AddError("Error: Could not find any DFSU Elements with a value different from zero")
        else:
            x_limit = [np.min(element_coordinates[elements_with_water,0]) - searchDistance,
                   np.max(element_coordinates[elements_with_water,0]) + searchDistance]
            y_limit = [np.min(element_coordinates[elements_with_water,1]) - searchDistance,
                   np.max(element_coordinates[elements_with_water,1]) + searchDistance]

            statusUpdate("Creating empty raster with extent: %d-%d and %d-%d" % (x_limit[0], x_limit[1], y_limit[0], y_limit[1]), tic)
            raster_xs_vector = np.arange(x_limit[0],x_limit[1],raster_cell_size)
            raster_ys_vector = np.arange(y_limit[0],y_limit[1],raster_cell_size)

            raster_x, raster_y = np.meshgrid(raster_xs_vector, raster_ys_vector)
            raster_depth = np.zeros(raster_x.shape)
            raster_x_flat = raster_x.flatten()
            raster_y_flat = raster_y.flatten()
            raster_depth_flat = np.zeros(raster_x.flatten().shape+tuple([len(dfs_read_data)]))

            statusUpdate("Retrieving Raster Elements near water", tic)
            arcpy.AddMessage("%d raster elements. %d elements with water" % (len(element_coordinates[elements_with_water]), len(raster_depth_flat)))
            raster_cKDTree = cKDTree(np.array((raster_x_flat, raster_y_flat)).transpose())

            idx = set()
            # for element_i, element in enumerate(element_coordinates[elements_with_water]):
            raster_points_near_water = raster_cKDTree.query_ball_point(element_coordinates[elements_with_water, :2],
                                                                       r=searchDistance)
            for points_list in raster_points_near_water:
                idx.update(points_list)

            idx_remove = []
            statusUpdate("Removing raster elements that overlap clip_layers",tic)
            if clip_layers:
                for i in idx:
                    point = arcpy.Point(raster_x_flat[i],raster_y_flat[i])
                    for clip_shape in clip_shapes:
                        if clip_shape.contains(point):
                            idx_remove.append(i)
                for i in idx_remove:
                    idx.remove(i)

            idx = list(idx)

            arcpy.CreateFeatureclass_management("in_memory", "ClipPolygon", "POLYGON")
            boundary_xy_table = dfs.geometry.boundary_polylines[1][0].xy
            polygons = arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table]))
            with arcpy.da.InsertCursor("in_memory\ClipPolygon", "SHAPE@") as cursor:
                cursor.insertRow([polygons])

            arcpy.CreateFeatureclass_management("in_memory", "CutPolygon", "POLYGON")
            cut_polygons = dfs.geometry.boundary_polylines[3]
            polygons = []
            for cut_polygon in cut_polygons:
                boundary_xy_table = cut_polygon.xy
                polygons.append(arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table])))
            with arcpy.da.InsertCursor("in_memory\CutPolygon", "SHAPE@") as cursor:
                for polygon in polygons:
                    cursor.insertRow([polygon])

            # arcpy.ddd.CreateTin(TIN, dfs.projection_string,
            #                     r"in_memory\nodes_Z Shape.Z Mass_Points; in_memory\ClipPolygon <None> Hard_Clip")

            statusUpdate("Creating KDTree",tic)
            elements_searchable = np.where((x_limit[0] < element_coordinates[:,0]) &
                                           (element_coordinates[:,0] < x_limit[1]) &
                                           (y_limit[0] < element_coordinates[:,1]) &
                                           (element_coordinates[:,1] < y_limit[1]))[0]
            dfsu_cKDTree = cKDTree(element_coordinates[:, 0:2])
            closest_elements = dfsu_cKDTree.query(np.array((raster_x_flat, raster_y_flat))[:, list(idx)].transpose())[1]
            statusUpdate("Interpolating DFSU to Raster (nearest neighbor)", tic)

            for raster_i, element_i in enumerate(closest_elements):
                for item_i in range(len(dfs_read_data)):
                    if DFSUFieldSummary == "Max":
                        raster_depth_flat[idx[raster_i], item_i] = np.max(dfs_read_data[item_i][:, element_i])
                    elif DFSUFieldSummary == "Last Timestep":
                        raster_depth_flat[idx[raster_i], item_i] = dfs_read_data[item_i][-1, element_i]

            raster_points_with_value = np.where(raster_depth_flat > 0)[0]
            arcpy.AddMessage(np.max(raster_depth_flat))
            idx_remove = ~dfs.geometry.contains(
                np.array((raster_x_flat, raster_y_flat))[:, raster_points_with_value].transpose())
            arcpy.AddMessage(idx_remove)
            if idx_remove.size>0:
                raster_depth_flat[
                np.array([raster_points_with_value[i] for i in [j for j, val in enumerate(idx_remove) if val]])] = 0

            for item_i in range(len(DFSUFields)):
                statusUpdate("Saving Raster", tic)
                raster_depth = raster_depth_flat[:,item_i].reshape(raster_depth.shape[0:2] + tuple([1]))

                raster_depth_compressed = raster_depth
                raster = arcpy.NumPyArrayToRaster(np.flip(raster_depth_compressed[:,:,0],axis=0), lower_left_corner = arcpy.Point(x_limit[0],y_limit[0]),
                                         x_cell_size = raster_cell_size,
                                         y_cell_size = raster_cell_size,
                                         value_to_nodata = 0)
                arcpy.AddMessage(DFSUFile.split("\\")[-1])

                raster_in_memory = os.path.join("in_memory", DFSUFile.split("\\")[-1].replace(".dfsu","_%s %s.tif" % (DFSUFields[item_i],DFSUFieldSummary)))
                raster.save(raster_in_memory)

                arcpy.Clip_management(
                    in_raster=in_memory,
                    rectangle="0 0 0 0",
                    out_raster=os.path.join(RasterOutputFolder, DFSUFile.split("\\")[-1].replace(".dfsu", "_%s %s.tif" % (
                                DFSUFields[item_i], DFSUFieldSummary))),
                    in_template_dataset=r"in_memory\ClipPolygon", nodata_value="-3.402823e+38",
                    clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT")


                raster.save(os.path.join(RasterOutputFolder, DFSUFile.split("\\")[-1].replace(".dfsu", "_%s %s.tif" % (
                DFSUFields[item_i], DFSUFieldSummary))))
        return

class DFS2ToRaster(object):
    def __init__(self):
        self.description = "Convert DFS2 to Raster (ArcGIS Pro only)"
        self.label       = "Convert DFS2 to Raster (ArcGIS Pro only)"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DFS2File = arcpy.Parameter(
            displayName="Input DFS2 File",
            name="DFS2File",
            datatype="File",
            parameterType="Required",
            direction="Input")
        DFS2File.filter.list = ["dfs2"]
        
        DFS2Fields = arcpy.Parameter(
            displayName="DFS2 field to export",
            name="DFS2Fields",
            datatype="String",
            multiValue="true",
            parameterType="Optional",
            direction="Input")
            
        RasterFileOutput = arcpy.Parameter(
            displayName="Output Raster File",
            name="RasterFileOutput",
            datatype="File",
            parameterType="Required",
            direction="Output")
        RasterFileOutput.filter.list = ["tif"]
            
        
        parameters = [DFS2File, DFS2Fields, RasterFileOutput]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].altered:
            dfs = mikeio.dfs2.Dfs2(parameters[0].valueAsText)
            parameters[1].filter.list = [item.name for item in dfs.items]
            
        return

    def updateMessages(self, parameters): #optional
        
        return

    def execute(self, parameters, messages):
        DFS2File = parameters[0].valueAsText
        DFS2Fields = parameters[1].valueAsText
        RasterFileOutput = parameters[2].valueAsText
        
        dfs2 = mikeio.dfs2.Dfs2(DFS2File)

        item = dfs2.items[0].name

        data = np.nan_to_num(dfs2.read(items = item).data[0],nan=-999)
        dfs_read = dfs2.read(items=[i for i,item in enumerate(dfs2.items) if item.name in DFS2Fields])


        lower_left_corner = arcpy.PointGeometry(arcpy.Point(dfs2.longitude, dfs2.latitude),"GCS_WGS_1984").projectAs("ETRS_1989_UTM_Zone_32N")[0]
        raster = arcpy.NumPyArrayToRaster(np.flip(data, axis=0), lower_left_corner = lower_left_corner, 
                                             x_cell_size = dfs2.dx,
                                             y_cell_size = dfs2.dy,
                                             value_to_nodata = -999)

        raster.save(RasterFileOutput)
        return
        
class DFSUPlumeStatistics(object):
    def __init__(self):
        self.label       = "Find and analyze volume of plumes from DFSU File"
        self.description = "Find and analyze volume of plumes from Flood Statistics File"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DFSUFile = arcpy.Parameter(
            displayName="Input DFSU File",
            name="DFSUFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        DFSUFile.filter.list = ["dfsu"]
        
        output_shape_file = arcpy.Parameter(
            displayName="Output shape file",
            name="output_shape_file",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")
        
        minimum_water_depth = arcpy.Parameter(
            displayName="Minimum water depth in order ot be included in a plume",
            name="minimum_water_depth",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        minimum_water_depth.value = 0.003
        
        parameters = [DFSUFile, output_shape_file, minimum_water_depth]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        
        return

    def updateMessages(self, parameters): #optional

        return

    def execute(self, parameters, messages):
        tic = time.time()
        DFSUFile = parameters[0].valueAsText
        output_shape_file = parameters[1].valueAsText
        minimum_water_depth = float(parameters[2].valueAsText)
        
        dfs = mikeio.dfsu.Dfsu2DH(DFSUFile)
        DFSUField = "Total water depth"

        # element_coordinates = dfs.element_coordinates
                
        dfs_read = dfs.read(items=[i for i,item in enumerate(dfs.items) if item.name in DFSUField])

        dfs_read_data = dfs_read.to_numpy()
        np.nan_to_num(dfs_read_data, copy = False)

        elements_with_water = np.where(dfs_read_data[0,-1,:]>minimum_water_depth)[0]

        element_table = dfs.element_table

        plumes_class = []
        plumes_nodes = []
        plumes_elements = []
        arcpy.SetProgressor("step", "Doing whatever...",
                            0, len(elements_with_water), 1)
        arcpy.AddMessage(len(elements_with_water))
        for element_i, element in enumerate(elements_with_water[:1000]):
            arcpy.SetProgressorPosition(element_i)
            matches_found = []
            for element_node in element_table[element]:
                matches = [i for i,nodes in enumerate(plumes_nodes) if element_node in nodes]
                for match in matches:
                    if match not in matches_found:
                        matches_found.append(match)
            if len(matches_found) > 0:
                plumes_elements[min(matches_found)].add(element)
                for node in element_table[element]:
                    plumes_nodes[min(matches_found)].add(node)
                for match in np.flip(np.sort(matches_found))[:-1]:
                    for node in plumes_nodes[match]:
                        plumes_nodes[min(matches_found)].add(node)
                    for element in plumes_elements[match]:
                        plumes_elements[min(matches_found)].add(element)
                    plumes_elements[min(matches_found)].add(element)
                    del plumes_class[match]
                    del plumes_nodes[match]
                    del plumes_elements[match]
            else:
                plume_class = np.max(plumes_class)+1 if len(plumes_class)>0 else 0
                plumes_class.append(plume_class)
                plumes_nodes.append(set(element_table[element]))
                plumes_elements.append(set([element]))
                
        output_file = arcpy.CreateFeatureclass_management(r"in_memory","boundary", "POLYGON")[0]
        arcpy.AddField_management(output_file, "Class", "LONG")
        arcpy.AddField_management(output_file, "Volume", "FLOAT", field_precision = 15, field_scale = 5)
        arcpy.AddField_management(output_file, "Area", "FLOAT", field_precision = 15, field_scale = 5)


        nodes_coordinates = dfs.node_coordinates

        with arcpy.da.InsertCursor(output_file,["SHAPE@","Class", "Volume"]) as cursor:
            for plume_i in range(len(plumes_class)):
                for element in plumes_elements[plume_i]:
                    node_coordinates = nodes_coordinates[element_table[element],:-1]
                    triangle = arcpy.Polygon(arcpy.Array([arcpy.Point(coords[0],coords[1]) for coords in node_coordinates]))
                    arcpy.AddMessage(dfs_read_data[0, -1, element])
                    cursor.insertRow([triangle, plume_i, dfs_read_data[-1,element]])

        with arcpy.da.UpdateCursor(output_file,["SHAPE@", "Volume", "Area"]) as cursor:
            for row in cursor:
                row[1] = row[1] * row[0].area
                row[2] = row[0].area
                cursor.updateRow(row)
        
        arcpy.Dissolve_management(output_file, output_shape_file , dissolve_field = ["Class"], statistics_fields = [["Volume", "SUM"],["Area","SUM"], ["Volume","MAX"]])
        
        def renameField(shape_file, old_field_name, new_field_name):
            arcpy.AddField_management(output_shape_file, new_field_name, "FLOAT", field_precision = 15, field_scale = 5)
            arcpy.CalculateField_management(output_shape_file, new_field_name, "!%s!" % (old_field_name), "PYTHON_9.3", "")
            arcpy.DeleteField_management(output_shape_file, old_field_name)
        
        field_name_rename_dictionairy = {"SUM_Volume": "Volume", "MAX_Volume": "Max_depth", "SUM_Area": "Area"}
        
        for old_field_name, new_field_name in field_name_rename_dictionairy.items():
            renameField(output_shape_file, old_field_name, new_field_name)
            
        # mxd = arcpy.mapping.MapDocument("CURRENT")
        # df = arcpy.mapping.ListDataFrames(mxd)[0]
        
        # group_name = os.path.splitext(os.path.basename(DFSUFile))[0]
        # empty_group_layer = [layer for layer in arcpy.mapping.ListLayers(mxd, group_name, df) if layer and layer.isGroupLayer is True]
        # if not empty_group_layer:
            # empty_group_mapped = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + r"\Data\EmptyGroup.lyr")
            # empty_group = arcpy.mapping.AddLayer(df, empty_group_mapped, "TOP")
            # empty_group_layer = arcpy.mapping.ListLayers(mxd, "Empty Group", df)[0]
            # empty_group_layer.name = group_name 

        # addLayer = arcpy.mapping.Layer(output_shape_file)
        # arcpy.mapping.AddLayerToGroup(df, empty_group_layer, addLayer, "TOP")
        # updatelayer = [layer for layer in arcpy.mapping.ListLayers(mxd, "", df) if layer.longName == "%s\Plume" % group_name]
        # sourcelayer = arcpy.mapping.Layer(os.path.dirname(os.path.realpath(__file__)) + "\Data\Template_Plume.lyr")
        # arcpy.mapping.UpdateLayer(df,updatelayer,sourcelayer,False)
        # updatelayer.replaceDataSource(unicode(addLayer.workspacePath), 'FILEGDB_WORKSPACE', unicode(output_shape_file.datasetName))
        return
        
class DFSUToPolygons(object):
    def __init__(self):
        self.label       = "Convert DFSU to Polygons"
        self.description = "Convert DFSU to Polygons"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        DFSUFile = arcpy.Parameter(
            displayName="Input DFSU File",
            name="DFSUFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        DFSUFile.filter.list = ["dfsu"]
        
        output_shape_file = arcpy.Parameter(
            displayName="Output shape file",
            name="output_shape_file",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")
        
        minimum_water_depth = arcpy.Parameter(
            displayName="Minimum water depth in order ot be included in a plume",
            name="minimum_water_depth",
            datatype="GPDouble",
            parameterType="Required",
            category="Additional Settings",
            direction="Input")
        minimum_water_depth.value = 0.003
        
        parameters = [DFSUFile, output_shape_file, minimum_water_depth]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        
        return

    def updateMessages(self, parameters): #optional
        
        return

    def execute(self, parameters, messages):
        tic = time.time()
        DFSUFile = parameters[0].valueAsText
        output_shape_file = parameters[1].valueAsText
        minimum_water_depth = float(parameters[2].valueAsText)
        
        statusUpdate("Reading DFSU-file",tic)
        dfs = mikeio.dfsu.Dfsu2DH(DFSUFile)
        DFSUField = "Total water depth"

        element_coordinates = dfs.geometry.element_coordinates
                
        dfs_read = dfs.read(items=[i for i,item in enumerate(dfs.items) if item.name in DFSUField])

        dfs_read_data = np.nan_to_num(dfs_read.to_numpy(), copy = False)
        
        statusUpdate("Finding wet elements",tic)
        elements_with_water = np.where(dfs_read_data[0,-1,:]>minimum_water_depth)[0]
        
        statusUpdate("Retrieving element table",tic)
        element_table = dfs.element_table
        
        statusUpdate("Creating polygon class",tic)
        output_file = arcpy.CreateFeatureclass_management(r"in_memory","boundary", "POLYGON")[0]
        arcpy.AddField_management(output_file, "Depth", "FLOAT", field_precision = 15, field_scale = 5)

        statusUpdate("Retrieving node coordinates",tic)
        nodes_coordinates = dfs.node_coordinates
        
        arcpy.SetProgressor("step", "Inserting wet polygons into feature class...",
                    0, len(elements_with_water), 1)
        with arcpy.da.InsertCursor(output_file,["SHAPE@", "Depth"]) as cursor:
            for element_i, element in enumerate(elements_with_water):
                    node_coordinates = nodes_coordinates[element_table[element],:-1]
                    triangle = arcpy.Polygon(arcpy.Array([arcpy.Point(coords[0],coords[1]) for coords in node_coordinates]))
                    cursor.insertRow([triangle, dfs_read_data[0,-1,element]])
                    arcpy.SetProgressorPosition(element_i)
        arcpy.CopyFeatures_management(output_file, output_shape_file)
        return
        
class CutFromMesh(object):
    def __init__(self):
        self.label       = "Remove Elements from Mesh based on shapefile"
        self.description = "Remove Elements from Mesh based on shapefile"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        MeshFile = arcpy.Parameter(
            displayName="Input Mesh File",
            name="MeshFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MeshFile.filter.list = ["mesh"]
        
        crop_shapefile = arcpy.Parameter(
            displayName="Input Shapefile (crop elements that overlap shapefile)",
            name="crop_shapefile",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        
        MeshFileOutput = arcpy.Parameter(
            displayName="Output Mesh File with elements removed",
            name="MeshFileOutput",
            datatype="File",
            parameterType="Required",
            direction="Output")
        MeshFileOutput.filter.list = ["mesh"]
        
        parameters = [MeshFile, crop_shapefile, MeshFileOutput]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].ValueAsText and not parameters[2].ValueAsText:
            parameters[2].value = parameters[0].ValueAsText.replace(".mesh","_cropped.mesh")
        return

    def updateMessages(self, parameters): #optional
        
        return

    def execute(self, parameters, messages):
        tic = time.time()
        MeshFile = parameters[0].valueAsText
        crop_shapefile = parameters[1].valueAsText
        MeshFileOutput = parameters[2].valueAsText
                
        try:
            arcpy.management.Delete("in_memory\Points")
            arcpy.management.Delete("in_memory\crop_dissolved")
        except Exception as e:
            pass

        mesh = mikeio.dfsu.Mesh(MeshFile)
        crop_dissolved = arcpy.Dissolve_management(crop_shapefile, "in_memory\crop_dissolved", multi_part = "MULTI_PART")[0]
            
        element_coordinates = mesh.element_coordinates[:,0:2]

        points_filepath = arcpy.CreateFeatureclass_management("in_memory", "Points", "POINT")[0]
        with arcpy.da.InsertCursor(points_filepath, "SHAPE@") as cursor:
            for point_shape in [arcpy.Point(element[0], element[1]) for element in element_coordinates]:
                cursor.insertRow([point_shape])
                
        points_within = arcpy.management.SelectLayerByLocation(points_filepath, overlap_type = "WITHIN", select_features = crop_dissolved, invert_spatial_relationship = "INVERT")[0]
        element_IDs = [ID-1 for ID in list(points_within.getSelectionSet())]

        mesh.write(MeshFileOutput, elements = element_IDs)
        
        return
        
class CutFromDFSU(object):
    def __init__(self):
        self.label       = "Remove Elements from DFSU based on shapefile"
        self.description = "Remove Elements from DFSU based on shapefile"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        MeshFile = arcpy.Parameter(
            displayName="Input DFSU File",
            name="MeshFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MeshFile.filter.list = ["DFSU"]
        
        crop_shapefile = arcpy.Parameter(
            displayName="Input Shapefile (crop elements that overlap shapefile)",
            name="crop_shapefile",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        
        MeshFileOutput = arcpy.Parameter(
            displayName="Output DFSU File with elements removed",
            name="MeshFileOutput",
            datatype="File",
            parameterType="Required",
            direction="Output")
        MeshFileOutput.filter.list = ["DFSU"]
        
        parameters = [MeshFile, crop_shapefile, MeshFileOutput]
        
        return parameters

    def isLicensed(self): #optional
        return True

    def updateParameters(self, parameters): #optional
        if parameters[0].ValueAsText and not parameters[2].ValueAsText:
            parameters[2].value = os.path.join(os.path.dirname(parameters[0].ValueAsText), "TIN")
        return

    def updateMessages(self, parameters): #optional
        
        return

    def execute(self, parameters, messages):
        tic = time.time()
        MeshFile = parameters[0].valueAsText
        crop_shapefile = parameters[1].valueAsText
        MeshFileOutput = parameters[2].valueAsText
                
        try:
            arcpy.management.Delete("in_memory\Points")
            arcpy.management.Delete("in_memory\crop_dissolved")
        except Exception as e:
            pass

        mesh = mikeio.dfsu.Dfsu(MeshFile)
        crop_dissolved = arcpy.Dissolve_management(crop_shapefile, "in_memory\crop_dissolved", multi_part = "MULTI_PART")[0]
            
        element_coordinates = mesh.element_coordinates[:,0:2]

        points_filepath = arcpy.CreateFeatureclass_management("in_memory", "Points", "POINT")[0]
        with arcpy.da.InsertCursor(points_filepath, "SHAPE@") as cursor:
            for point_shape in [arcpy.Point(element[0], element[1]) for element in element_coordinates]:
                cursor.insertRow([point_shape])
                
        points_within = arcpy.management.SelectLayerByLocation(points_filepath, overlap_type = "WITHIN", select_features = crop_dissolved)[0]
        element_IDs = [ID-1 for ID in list(points_within.getSelectionSet())]

        mesh.write(MeshFileOutput, data = mesh.read(elements = element_IDs), items = mesh.items, elements = element_IDs)
        
        return

# if __name__ == '__main__':
#     DHMFile = r"\\files\Projects\RWA2022N001XX\RWA2022N00174\Model\04_DTM\Viborg_Viborg (3).tif"
#     MeshFile = "C:\Offline\VOR_Status\Mesh\Mesh_v1.1.mesh"
#     MeshFileOutput = MeshFile.replace(".mesh","_test.mesh")
#
#     # statusUpdate("Reading Mesh", tic)
#     dfs = mikeio.dfsu.Mesh(MeshFile)
#
#     # statusUpdate("Reading Mesh: Getting Node Coordinates", tic)
#     data = dfs.node_coordinates
#     node_codes = dfs.codes
#
#     # statusUpdate("Reading Raster", tic)
#     DHMRaster = np.flip(arcpy.RasterToNumPyArray(DHMFile, nodata_to_value=0), axis=0)
#
#     DHMRasterInfo = np.array([float(arcpy.GetRasterProperties_management(DHMFile, "LEFT")[0].replace(",", ".")),
#                               float(arcpy.GetRasterProperties_management(DHMFile, "BOTTOM")[0].replace(",", ".")),
#                               float(arcpy.GetRasterProperties_management(DHMFile, "CELLSIZEX")[0].replace(",", ".")),
#                               float(arcpy.GetRasterProperties_management(DHMFile, "CELLSIZEY")[0].replace(",", "."))])
#     DHMRasterXs = np.arange(DHMRasterInfo[0], DHMRasterInfo[0] + DHMRaster.shape[1] * DHMRasterInfo[2],
#                             float(DHMRasterInfo[2])).astype(np.float32)
#     DHMRasterYs = np.arange(DHMRasterInfo[1], DHMRasterInfo[1] + DHMRaster.shape[0] * DHMRasterInfo[3],
#                             float(DHMRasterInfo[3])).astype(np.float32)
#
#     DHMRasterX, DHMRasterY = np.meshgrid(DHMRasterXs, DHMRasterYs)
#     DHMRasterXsFlat = DHMRasterX.flatten()
#     DHMRasterYsFlat = DHMRasterY.flatten()
#     DHMRasterFlat = DHMRaster.flatten()
#
#     idx = np.where(DHMRasterFlat != 0)
#     DHMRasterXsFlatSparse = DHMRasterXsFlat[idx]
#     DHMRasterYsFlatSparse = DHMRasterYsFlat[idx]
#     DHMRasterFlatSparse = DHMRasterFlat[idx]
#
#     tic = time.time()
#     bedLevelFlat = np.zeros((data.shape[0]), dtype=np.float64)
#
#     # statusUpdate("Interpolating raster to mesh (nearest neighbor)", tic)
#     for row in range(int(len(bedLevelFlat))):
#         bedLevelFlat[row] = DHMRaster[
#             bisect.bisect_left(DHMRasterYs, data[row, 1]), bisect.bisect_left(DHMRasterXs, data[row, 0])]
#
#     # statusUpdate("Saving mesh", tic)
#     dfs.zn = bedLevelFlat
#     dfs.write(MeshFileOutput)

class MeshToTIN(object):
    def __init__(self):
        self.label = "Convert Mesh to TIN (ArcGIS Pro only)"
        self.description = "Convert Mesh to TIN (ArcGIS Pro only)"
        self.canRunInBackground = True

    def getParameterInfo(self):
        # Define parameter definitions

        # Input Features parameter
        MeshFile = arcpy.Parameter(
            displayName="Input Mesh File",
            name="MeshFile",
            datatype="File",
            parameterType="Required",
            direction="Input")
        MeshFile.filter.list = ["mesh"]

        TIN = arcpy.Parameter(
            displayName="Output TIN",
            name="TIN",
            datatype="DETin",
            parameterType="Required",
            direction="Output")

        clip_TIN = arcpy.Parameter(
            displayName="Clip TIN to boundary polylines?",
            name="clip_TIN",
            datatype="Boolean",
            parameterType="optional",
            direction="Input")

        parameters = [MeshFile, TIN, clip_TIN]

        return parameters

    def isLicensed(self):  # optional
        return True

    def updateParameters(self, parameters):  # optional
        if parameters[0].ValueAsText and not parameters[1].ValueAsText:
            parameters[1].value = parameters[0].ValueAsText.replace(".mesh","")
        return

    def updateMessages(self, parameters):  # optional
        return

    def execute(self, parameters, messages):
        arcpy.SetProgressorLabel("Loading Mesh File")
        mesh_file = parameters[0].ValueAsText
        TIN = parameters[1].ValueAsText
        arcpy.AddMessage(TIN)
        clip_TIN = parameters[2].Value
        dfs = mikeio.dfsu.Mesh(mesh_file)
        arcpy.SetProgressorLabel("Reading Mesh File Element Coordinates")
        node_coordinates = dfs.node_coordinates

        arcpy.SetProgressorLabel("Creating Feature Class")
        arcpy.CreateFeatureclass_management("in_memory", "nodes_Z", "POINT", has_z="Enabled")
        # arcpy.AddField_management(os.path.join(nodesPath, nodesName), "Z", "DOUBLE")

        arcpy.SetProgressorLabel("Adding each element to the feature class")
        with arcpy.da.InsertCursor(os.path.join("in_memory", "nodes_Z"), ["SHAPE@"]) as cursor:
            for node in node_coordinates:
                cursor.insertRow([arcpy.Point(node[0], node[1], node[2])])

        if clip_TIN:
            arcpy.CreateFeatureclass_management("in_memory", "ClipPolygon", "POLYGON")
            boundary_xy_table = dfs.geometry.boundary_polylines[1][0].xy
            polygons = arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table]))
            with arcpy.da.InsertCursor("in_memory\ClipPolygon", "SHAPE@") as cursor:
                cursor.insertRow([polygons])

            arcpy.CreateFeatureclass_management("in_memory", "CutPolygon", "POLYGON")
            cut_polygons = dfs.geometry.boundary_polylines[3]
            polygons = []
            for cut_polygon in cut_polygons:
                boundary_xy_table = cut_polygon.xy
                polygons.append(arcpy.Polygon(arcpy.Array([arcpy.Point(xy[0], xy[1]) for xy in boundary_xy_table])))
            with arcpy.da.InsertCursor("in_memory\CutPolygon", "SHAPE@") as cursor:
                for polygon in polygons:
                    cursor.insertRow([polygon])

            arcpy.ddd.CreateTin(TIN, dfs.projection_string,
                                r"in_memory\nodes_Z Shape.Z Mass_Points; in_memory\ClipPolygon <None> Hard_Clip")
        else:
            arcpy.ddd.CreateTin(TIN, dfs.projection_string,
                                "%s Shape.Z Mass_Points")
        return

if __name__ == "__main__":
    # DHMFile = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\04_DTM\DTM_Status.tif"
    DHMFile = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\04_DTM\Løsningsforslag 2023-12-18.tif"
    mesh_files = [
                    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\VBT_Mesh_2023_11_29.mesh"
                ]
    mesh_files = np.flip(mesh_files)
    for MeshFile in mesh_files:
        print(MeshFile)

        MeshFileOutput = os.path.join(r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Nødoverløb med jordbassiner", os.path.basename(MeshFile).replace(".mesh", "_interp.mesh"))

        dfs = mikeio.dfsu.Mesh(MeshFile)

        data = dfs.node_coordinates

        DHMRaster = np.flip(arcpy.RasterToNumPyArray(DHMFile, nodata_to_value=0), axis=0)

        DHMRasterInfo = np.array([float(arcpy.GetRasterProperties_management(DHMFile, "LEFT")[0].replace(",", ".")),
                                  float(arcpy.GetRasterProperties_management(DHMFile, "BOTTOM")[0].replace(",", ".")),
                                  float(arcpy.GetRasterProperties_management(DHMFile, "CELLSIZEX")[0].replace(",", ".")),
                                  float(arcpy.GetRasterProperties_management(DHMFile, "CELLSIZEY")[0].replace(",", "."))])

        raster = arcpy.Raster(DHMFile)
        lower_left_corner = arcpy.Point(raster.extent.XMin, raster.extent.YMin)
        DHMRasterXs = lower_left_corner.X + np.arange(raster.width) * raster.meanCellWidth + raster.meanCellWidth/2
        DHMRasterYs = lower_left_corner.Y + np.arange(raster.height) * raster.meanCellHeight + raster.meanCellHeight/2

        # DHMRasterXs = np.arange(DHMRasterInfo[0],DHMRasterInfo[0]+DHMRaster.shape[1]*DHMRasterInfo[2],float(DHMRasterInfo[2])).astype(np.float32)
        # DHMRasterYs = np.arange(DHMRasterInfo[1],DHMRasterInfo[1]+DHMRaster.shape[0]*DHMRasterInfo[3],float(DHMRasterInfo[3])).astype(np.float32 )

        interp = RegularGridInterpolator((DHMRasterYs, DHMRasterXs),
                                         DHMRaster)

        # DHMRasterX, DHMRasterY = np.meshgrid(DHMRasterXs,DHMRasterYs)
        # DHMRasterXsFlat = DHMRasterX.flatten()
        # DHMRasterYsFlat = DHMRasterY.flatten()
        # DHMRasterFlat = DHMRaster.flatten()

        # idx = np.where(DHMRasterFlat!=0)
        # DHMRasterXsFlatSparse = DHMRasterXsFlat[idx]
        # DHMRasterYsFlatSparse = DHMRasterYsFlat[idx]
        # DHMRasterFlatSparse = DHMRasterFlat[idx]

        tic = time.time()
        bedLevelFlat = np.zeros((data.shape[0]), dtype=np.float64)

        bedLevelFlat = interp(data[:, [1, 0]])
        # for row in range(int(len(bedLevelFlat))):
        #     try:
        #         bedLevelFlat[row] = interp([bisect.bisect_left(DHMRasterYs, data[row,1]),bisect.bisect_left(DHMRasterXs, data[row,0])]
        #     except Exception as e:
        #         bedLevelFlat[row] = 0

        dfs.zn = bedLevelFlat
        dfs.write(MeshFileOutput)