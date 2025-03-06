import arcpy
import pythonaddins
import numpy as np
import os
from scipy.spatial import cKDTree

def find_duplicates(strings):
    return list(set([x for x in strings if strings.count(x) > 1]))

class Node:
    def __init__(self, muid, x, y):
        self.muid = muid
        self.x = x
        self.y = y

class CatchmentLayer(object):
    """Implementation for CatchmentConnector.catchmentLayer (ComboBox)"""
    def __init__(self):
        self.items = ["item1", "item2", "item3", "item4", "item5"]
        self.editable = True
        self.enabled = True
        self.dropdownWidth = 'WWWWWWWWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWWWWWWWW'
        self.selectedLayer = None
    def onSelChange(self, selection):
        if selection != None:
            self.selectedLayer = selection
            catchLayer = \
            [layer for layer in arcpy.mapping.ListLayers(self.mxd) if layer.longName == self.selectedLayer][0]
            msm_CatchCon = os.path.join(catchLayer.workspacePath, "msm_CatchCon")

            duplicates = find_duplicates([row[0] for row in arcpy.da.SearchCursor(msm_CatchCon, ["MUID"])])
            if duplicates:
                # arcpy.AddMessage(["('%s')" % "', '".join(duplicates)])
                if pythonaddins.MessageBox(
                        ["Duplicates found in msm_CatchCon MUID: ('%s'). Continue?" % "', '".join(duplicates)],
                        "Continue", 1) != "OK":
                    exit()
            print("BOB11")
        else:
            self.selectedLayer = None
    def onEditChange(self, text):
        pass
    def onFocus(self, focused):
        if focused:
            self.mxd = arcpy.mapping.MapDocument("Current")
            layers = arcpy.mapping.ListLayers(self.mxd)
            self.items=[]
            if len(layers) != 0:
                for layer in layers:
                    if layer.isFeatureLayer and arcpy.Describe(layer).shapetype == "Polygon":
                        self.items.append(layer.longName)
    def onEnter(self):
        pass
    def refresh(self):
        pass

class ConnectCatchment(object):
    def __init__(self):
        self.enabled = True
        self.shape = "NONE" # Can set to "Line", "Circle" or "Rectangle" for interactive shape drawing and to activate the onLine/Polygon/Circle event sinks.
        self.mxd = arcpy.mapping.MapDocument("CURRENT")
    def onMouseDown(self, x, y, button, shift):
        pass
    def onMouseDownMap(self, x, y, button, shift):
        nodeLayer = [layer for layer in arcpy.mapping.ListLayers(self.mxd) if layer.longName == manholeLayer.selectedLayer][0]
        # if arcpy.Exists(os.path.join(nodeLayer.workspacePath, "msm_Node")):
        #     nodeLayer = os.path.join(nodeLayer.workspacePath, "msm_Node")

        # nodesCount = int(arcpy.GetCount_management(nodeLayer)[0])

        node_i = manholeLayer.find_closest_node([x,y])
        node = manholeLayer.nodes[node_i].muid

        # dist = np.sqrt(np.power(nodesX-x,2)+np.power(nodesY-y,2))
        # node = nodes[np.argmin(dist)]

        catchLayer = [layer for layer in arcpy.mapping.ListLayers(self.mxd) if layer.longName == catchmentLayer.selectedLayer][0]
        mike_urban_database = os.path.dirname(arcpy.Describe(catchLayer).catalogPath).replace("\mu_Geometry", "")
        is_sqlite = True if ".sqlite" in mike_urban_database else False

        catchments_selected = []
        ID_fields = ["MUID","FID","OBJECTID","OID"]
        fields = arcpy.ListFields(catchLayer)
        for ID_field in ID_fields:
            if ID_field.lower() in [field.name.lower() for field in fields]:
                break

        with arcpy.da.SearchCursor(catchLayer, [ID_field]) as cursor:
            for i,row in enumerate(cursor):
                catchments_selected.append(str(row[0]))

        if len(arcpy.ListFields(catchLayer,"NodeID"))>0:
            if pythonaddins.MessageBox("Assign %s to %s (SHP)?" % (", ".join(catchments_selected),node), "Confirm Assignment", 1) == "OK":
                arcpy.CalculateField_management(catchLayer, "NodeID", '"%s"' % node, "VB")
        else:
            if is_sqlite:
                pass
            else:
                msm_CatchConLink = os.path.join(catchLayer.workspacePath, "msm_CatchConLink")
                duplicates = find_duplicates([row[0] for row in arcpy.da.SearchCursor(msm_CatchConLink, ["MUID"])])
                if duplicates:
                    # arcpy.AddMessage(["('%s')" % "', '".join(duplicates)])
                    if pythonaddins.MessageBox(
                            ["Duplicates found in msm_CatchConLink MUID: ('%s'). Continue?" % "', '".join(duplicates)],
                            "Continue", 1) != "OK":
                        exit()
            if pythonaddins.MessageBox("Assign %s to %s?" % (", ".join(catchments_selected),node), "Confirm Assignment", 1) == "OK":
                msm_Node = os.path.join(catchLayer.workspacePath, "msm_Node")
                ms_Catchment = os.path.join(catchLayer.workspacePath, "ms_Catchment")
                msm_CatchCon = os.path.join(catchLayer.workspacePath, "msm_CatchCon")

                catchments_existing = []

                if is_sqlite:
                    with arcpy.da.SearchCursor(msm_CatchCon.replace("!delete!",""), ['CatchID', 'NodeID'], where_clause = "CatchID IN ('%s')" % ("', '".join(catchments_selected))) as cursor:
                        for row in cursor:
                            catchments_existing.append(row[0])

                    import sqlite3
                    conn = sqlite3.connect(mike_urban_database.replace("!delete!",""))
                    # conn.enable_load_extension(True)
                    # conn.execute('SELECT load_extension("mod_spatialite")')
                    cursor = conn.cursor()
                    if catchments_existing:
                        update_query = r"UPDATE msm_CatchCon SET nodeid = '%s' WHERE catchid IN ('%s')" % (node, "', '".join(catchments_existing))
                        print(update_query)
                        cursor.execute(update_query)
                    catchments_not_existing = list(set(catchments_selected) - set(catchments_existing))
                    for catchment in catchments_not_existing:
                        try:
                            insert_query = r"INSERT INTO msm_CatchCon (catchid, nodeid) VALUES ('%s', '%s')" % (node, catchment)
                            print(insert_query)
                            cursor.execute(insert_query)
                        except Exception as e:
                            pythonaddins.MessageBox(e.message + "\nCannot insert catchment connections in sqlite", "Error", 0)
                    conn.commit()
                    # for catchment in catchment_selected:
                        # insert_query = r"INSERT INTO msm_CatchCon (catchid, nodeid) VALUES (%s, %s)" % (catchment, node)
                    conn.close()
                else:
                    with arcpy.da.UpdateCursor(msm_CatchCon, ['CatchID', 'NodeID'], where_clause = "CatchID IN ('%s')" % ("', '".join(catchments_selected))) as cursor:
                        for row in cursor:
                            row[1] = node
                            catchments_existing.append(row[0])
                            cursor.updateRow(row)

                    for MUID in catchments_selected:
                        if MUID not in catchments_existing:
                            with arcpy.da.InsertCursor(msm_CatchCon, ['CatchID', 'NodeID']) as cursor:
                                cursor.insertRow([MUID,node])

                    where_clause = "MUID IN ('%s')" % "', '".join(catchments_selected)
                    catchments_coordinates = {row[0]: row[1] for row in
                                              arcpy.da.SearchCursor(ms_Catchment, ["MUID", "SHAPE@XY"],
                                                                    where_clause=where_clause)}
                    links = {row[0]: [row[1], row[2]] for row in
                             arcpy.da.SearchCursor(msm_CatchCon, ["MUID", "CatchID", "NodeID"],
                                                   where_clause=where_clause.replace("MUID", "CatchID")) if row[1] in catchments_selected}

                    nodes_MUIDs = [node for catchment, node in links.values()]

                    nodes_coordinates = {row[0]: row[1] for row in arcpy.da.SearchCursor(msm_Node, ["MUID", "SHAPE@XY"],
                                                                                         where_clause="MUID IN ('%s')" % ("', '".join(
                                                                                             nodes_MUIDs)))}
                    msm_CatchConLink_MUIDs = [link for link in links.values()]

                    with arcpy.da.UpdateCursor(msm_CatchConLink, ["CatchConID"], where_clause =
                                                "CatchConID IN (%s)" % ", ".join([str(key) for key in links.keys()])) as cursor:
                        for row in cursor:
                            cursor.deleteRow()

                    with arcpy.da.InsertCursor(msm_CatchConLink,
                                               ["SHAPE@", "CatchConID", "SHAPE_Length"]) as cursor:
                        for link_i, link in enumerate(links.keys()):
                            catchment, node = links[link]
                            catchment_coordinate = catchments_coordinates[catchment]
                            node_coordinate = nodes_coordinates[node]
                            shape = arcpy.Polyline(
                                arcpy.Array([arcpy.Point(catchment_coordinate[0], catchment_coordinate[1]),
                                             arcpy.Point(node_coordinate[0], node_coordinate[1])]))
                            row = [shape, link, shape.length]
                            cursor.insertRow(row)
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

class ManholeLayer(object):
    """Implementation for CatchmentConnector.manholeLayer (ComboBox)"""
    def __init__(self):
        self.items = ["item1", "item2", "item3", "item4", "item5"]
        self.editable = True
        self.enabled = True
        self.dropdownWidth = 'WWWWWWWWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWWWWWWWW'
        self.selectedLayer = None
        self.nodes = {}
        self.tree = None
    def onSelChange(self, selection):
        if selection != None:
            self.selectedLayer = selection
            self.nodes = []

            if len(arcpy.ListFields(self.selectedLayer, "MUID")) > 0:
                with arcpy.da.SearchCursor(self.selectedLayer, ['MUID', 'SHAPE@XY']) as cursor:
                    for i, row in enumerate(cursor):
                        self.nodes.append(Node(row[0], row[1][0], row[1][1]))
            else:
                with arcpy.da.SearchCursor(nodeLayer, ['NodeID', 'SHAPE@XY']) as cursor:
                    for i, row in enumerate(cursor):
                        self.nodes.append(Node(row[0], row[1][0], row[1][1]))

            self.nodes_kdtree = cKDTree([(node.x, node.y) for node in self.nodes])

        else: 
            self.selectedLayer = None
    def onEditChange(self, text):
        pass
    def onFocus(self, focused):
        if focused:
            self.mxd = arcpy.mapping.MapDocument("Current")
            layers = arcpy.mapping.ListLayers(self.mxd)
            self.items=[]
            if len(layers) != 0:
                for layer in layers:
                    if layer.isFeatureLayer and arcpy.Describe(layer).shapetype == "Point":
                        self.items.append(layer.longName)
    def onEnter(self):
        pass
    def refresh(self):
        pass

    def find_closest_node(self, coordinate):
        distance, index = self.nodes_kdtree.query([coordinate[0], coordinate[1]])
        return index