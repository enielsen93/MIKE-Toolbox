# -*- coding: utf-8 -*-     
import arcpy
import pythonaddins
import numpy as np
import os
import warnings

class Node:
    def __init__(self, MUID, shape):
        self.MUID = MUID
        self.shape = shape

class ButtonClass21(object):
    """Implementation for TraceNetwork_addin.find_all_connected_catchments (Button)"""
    def __init__(self):
        self.enabled = True
        self.checked = False
    def onClick(self):
        nodes_MUID = [row[0] for row in arcpy.da.SearchCursor(group_layer.msm_Node, ["MUID"])]
        print(nodes_MUID)
        catchments = group_layer.graph.find_connected_catchments(nodes_MUID)
        catchments_MUID = [catchment.MUID for catchment in catchments]

        print("MUID IN ('%s')" % ("', '".join(catchments_MUID)))
        arcpy.SelectLayerByAttribute_management(group_layer.ms_Catchment.longName, "ADD_TO_SELECTION",
                                                "MUID IN ('%s')" % ("', '".join(catchments_MUID)))

        pass
        
class ButtonClass25(object):
    """Implementation for TraceNetwork_addin.find_all_unconnected_catchments (Button)"""
    def __init__(self):
        self.enabled = True
        self.checked = False
    def onClick(self):
        group_layer.graph._read_catchments()
        nodes_MUID = [row[0] for row in arcpy.da.SearchCursor(os.path.join(group_layer.msm_Node.workspacePath, "msm_Node"), ["MUID"])]
        catchments_MUID = [catchment.MUID for catchment in group_layer.graph.catchments_dict.values() if not catchment.nodeID or catchment.nodeID not in nodes_MUID]
        print("MUID IN ('%s')" % ("', '".join(catchments_MUID)))
        arcpy.SelectLayerByAttribute_management(group_layer.ms_Catchment.longName, "ADD_TO_SELECTION",
                                                "MUID IN ('%s')" % ("', '".join(catchments_MUID)))
        pass

class ButtonClass26(object):
        def __init__(self):
            self.enabled = True
            self.checked = False

        def onClick(self):
            nodes_MUID = [row[0] for row in arcpy.da.SearchCursor(group_layer.msm_Node, ["MUID"])]
            # print(nodes_MUID)
            nodes_in_path, links_in_path = group_layer.graph.trace_between(nodes_MUID)
            # print((nodes_in_path, links_in_path))
            if group_layer.msm_Node and group_layer.msm_Node.visible:
                arcpy.SelectLayerByAttribute_management(group_layer.msm_Node.longName, "NEW_SELECTION",
                                                        "MUID IN ('%s')" % ("', '".join(nodes_in_path)))

            if group_layer.msm_Link and group_layer.msm_Link.visible:
                arcpy.SelectLayerByAttribute_management(group_layer.msm_Link.longName, "ADD_TO_SELECTION",
                                                        "MUID IN ('%s')" % ("', '".join(links_in_path   )))

            pass

class GroupComboBoxClass1(object):
    """Implementation for TraceNetwork_addin.combobox (ComboBox)"""
    def __init__(self):
        self.items = ["item1", "item2","item1", "item2","item1", "item2"]
        self.editable = True
        self.enabled = True
        self.dropdownWidth = 'WWWWWWWWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWWWWWWWW'
        self.msm_Node = None
        self.msm_Link = None
        self.ms_Catchment = None
        self.layers_in_group = []
        self.map_only = ""
        self.points_muid = None
        self.points_xy = None

    def onSelChange(self, selection):
        self.selected_group = selection
        self.layers_in_group = [layer for layer in arcpy.mapping.ListLayers(self.mxd) if self.selected_group + "\\" in layer.longName]
        print([(layer.longName, layer.dataSource) for layer in self.layers_in_group])
        msm_Node_layer = [layer for layer in self.layers_in_group if
                    u"Brønd" in layer.name or "msm_Node" in layer.name]
        print(msm_Node_layer)
        self.msm_Node = msm_Node_layer[0] if msm_Node_layer else None

        msm_Link_layer = [layer for layer in self.layers_in_group if
         "Ledning" in layer.name or "msm_Link" in layer.name]
        print(msm_Link_layer)
        self.msm_Link = msm_Link_layer[0] if msm_Link_layer else None
        self.map_only = "link"
        # if msm_Link_layer:
            # self.map_only += ", link" if self.map_only else "link"

        ms_Catchment_layer = [layer for layer in self.layers_in_group if
         "Delopland" in layer.name or "Catchment" in layer.name]
        print(ms_Catchment_layer)
        self.ms_Catchment = ms_Catchment_layer[0] if ms_Catchment_layer else None

        link_layers = {"msm_Link": "link", "msm_Weir": "weir", "msm_Orifice":"orifice", "msm_Pump":"pump"}
        for link_layer in link_layers.keys():
            matching_layers = [layer for layer in self.layers_in_group if layer.isFeatureLayer and link_layer in layer.datasetName and layer.visible]
            print(matching_layers)
            if matching_layers:
                self.map_only += ", %s" % link_layers[link_layer] if self.map_only else "%s" % link_layers[link_layer]

        if any(["Passive Regulations" in [layer.name for layer in self.layers_in_group if layer.visible and layer.isFeatureLayer]]):
            ignore_regulations = True
        else:
            ignore_regulations = False

        links_MUID = [row[0] for row in arcpy.da.SearchCursor(self.msm_Link, ["MUID"])]
        duplicate_links = [MUID for MUID in links_MUID if links_MUID.count(MUID)>1]
        if duplicate_links:
            pythonaddins.MessageBox("Error: Links with identical MUIDs: ('%s')" % ("', '".join(duplicate_links)), "Error: Identical MUIDs", 0)
        
        import mikegraph
        pump_test = os.path.join(self.msm_Node.workspacePath, "msm_Pump")
        print("Graphing %s" % (self.msm_Node.workspacePath.replace("!delete!","")))
        if arcpy.Exists(pump_test):
            self.graph = mikegraph.MikeNetwork(self.msm_Node.workspacePath.replace("!delete!",""), map_only = self.map_only, ignore_regulations=ignore_regulations)
        else:
            print("Assuming it's not an MIKE Urban database as %s does not exist" % (pump_test))
            self.graph = mikegraph.MikeNetwork(nodes_and_links = [self.msm_Node.dataSource, self.msm_Link.dataSource])
        self.graph.map_network()

        self.points_muid = []
        self.points_xy = np.zeros((int(arcpy.management.GetCount(self.msm_Node.dataSource.replace("%",""))[0]),2))
        if int(arcpy.management.GetCount(self.msm_Node.dataSource.replace("%",""))[0]) == 0:
            print(self.msm_Node.dataSource.replace("%",""))
            print(self.msm_Node.dataSource.replace("%", ""))

        with arcpy.da.SearchCursor(self.msm_Node.dataSource.replace("%",""),
                                   ["MUID", "SHAPE@"]) as cursor:
            for i, row in enumerate(cursor):
                self.points_xy[i, :] = [row[1].firstPoint.X, row[1].firstPoint.Y]
                self.points_muid.append(row[0])


        pass
    def onEditChange(self, text):
        pass
    def onFocus(self, focused):
        if focused:
            self.mxd = arcpy.mapping.MapDocument("Current")
            self.df = arcpy.mapping.ListDataFrames(self.mxd)[0]

            group_layers = [layer for layer in arcpy.mapping.ListLayers(self.mxd) if layer.isGroupLayer and "Annotation" not in layer.name]
            self.items = []
            if len(group_layers) != 0:
                for layer in group_layers:
                    self.items.append(layer.longName)

        pass
    def onEnter(self):
        pass
    def refresh(self):
        pass

class TraceDownstreamToolClass18(object):
    """Implementation for TraceNetwork_addin.tool_1 (Tool)"""
    def __init__(self):
        self.enabled = True
        self.shape = "NONE" # Can set to "Line", "Circle" or "Rectangle" for interactive shape drawing and to activate the onLine/Polygon/Circle event sinks.
    def onMouseDown(self, x, y, button, shift):
        pass
    def onMouseDownMap(self, x, y, button, shift):
        pass
    def onMouseUp(self, x, y, button, shift):
        pass
    def onMouseUpMap(self, x, y, button, shift):
        pass
    def onMouseMove(self, x, y, button, shift):
        pass
    def onMouseMoveMap(self, x, y, button, shift):
        pass
    def onDblClick(self):
        pass
    def onKeyDown(self, keycode, shift):
        pass
    def onKeyUp(self, keycode, shift):
        pass
    def deactivate(self):
        pass
    def onCircle(self, circle_geometry):
        pass
    def onLine(self, line_geometry):
        pass
    def onRectangle(self, rectangle_geometry):
        pass

class TraceUpstreamToolClass16(object):
    """Implementation for TraceNetwork_addin.tool (Tool)"""
    def __init__(self):
        self.enabled = True
        self.shape = "NONE" # Can set to "Line", "Circle" or "Rectangle" for interactive shape drawing and to activate the onLine/Polygon/Circle event sinks.
    def onMouseDown(self, x, y, button, shift):
        pass
    def onMouseDownMap(self, x, y, button, shift):

        def findClosestNode(coord):
            distances = np.sum(np.abs(group_layer.points_xy-[coord[0], coord[1]]),axis=1)
            index_closest = np.argmin(distances)
            muid = group_layer.points_muid[index_closest]
            return muid

        try:
            print("Point Click [%1.2f, %1.2f]" % (x, y))
            target = findClosestNode([x,y])
        except Exception as e:
            warnings.warn("Could not find closest node.")
            print(group_layer.points_xy)
            raise(e)

        print("Tracing upstream from node %s" % (target))
        
        upstream_nodes = group_layer.graph.find_upstream_nodes([target])[0]
        
        print(upstream_nodes)
        print("MUID IN ('%s')" % ("', '".join(upstream_nodes)))
        if group_layer.msm_Node and group_layer.msm_Node.visible:
            arcpy.SelectLayerByAttribute_management(group_layer.msm_Node.longName, "ADD_TO_SELECTION", "MUID IN ('%s')" % ("', '".join(upstream_nodes)))
            print((group_layer.msm_Node.longName, "ADD_TO_SELECTION", "MUID IN ('%s')" % ("', '".join(upstream_nodes))))

        if group_layer.msm_Link and group_layer.msm_Link.visible:
            links_MUID = [link.MUID for link in group_layer.graph.network.links.values() if
                          link.fromnode in upstream_nodes and link.tonode in upstream_nodes]
            arcpy.SelectLayerByAttribute_management(group_layer.msm_Link.longName, "ADD_TO_SELECTION",
                                                    "MUID IN ('%s')" % ("', '".join(links_MUID)))

        if group_layer.ms_Catchment and group_layer.ms_Catchment.visible:
            total_catchments = []

            catchments = group_layer.graph.find_connected_catchments(upstream_nodes)
            for catchment in catchments:
                total_catchments.append(catchment.MUID)

            arcpy.SelectLayerByAttribute_management(group_layer.ms_Catchment.longName, "ADD_TO_SELECTION",
                                                    "MUID IN ('%s')" % ("', '".join(total_catchments)))

        pass
    def onMouseUp(self, x, y, button, shift):
        pass
    def onMouseUpMap(self, x, y, button, shift):
        pass
    def onMouseMove(self, x, y, button, shift):
        pass
    def onMouseMoveMap(self, x, y, button, shift):
        pass
    def onDblClick(self):
        pass
    def onKeyDown(self, keycode, shift):
        pass
    def onKeyUp(self, keycode, shift):
        pass
    def deactivate(self):
        pass
    def onCircle(self, circle_geometry):
        pass
    def onLine(self, line_geometry):
        pass
    def onRectangle(self, rectangle_geometry):
        pass