"""
Network module for mikegraph
"""

import os
import sys
import arcpy.da
import numpy as np
import re
from scipy.spatial import cKDTree
from .utils import calculate_full_flow

class PipeNetwork:
    """
        Create NetworkX graphs from MIKE+ pipe network databases.

        This class extracts pipe network data from MIKE+ databases and converts it into
        NetworkX graph structures for network analysis. Supports both .mdb and .sqlite
        database formats, and can work with pipes, weirs, orifices, and pumps.

        Parameters:
            mike_urban_database (str, optional): Path to Mike database file (.mdb or .sqlite)
            nodes_and_links (list, optional): List of [nodes_filepath, links_filepath] for
                direct file input instead of database
            map_only (str, optional): Specify what to map.
                Defaults to empty string.
            filter_sql_query (str, optional): SQL query to filter network elements.
                Defaults to None.

        Example:
            >>> # Using database file
            >>> pipes = PipeNetwork("stormwater_model.mdb")
        """
    def __init__(self, mike_urban_database = None, nodes_and_links = None, map_only = "", filter_sql_query = None):
        self.mike_urban_database = mike_urban_database
        if self.mike_urban_database:
            is_sqlite = True if ".sqlite" in self.mike_urban_database else False
            fromnode_fieldname = "FROMNODE" if ".mdb" in self.mike_urban_database else "fromnodeid"
            tonode_fieldname = "TONODE" if ".mdb" in self.mike_urban_database else "tonodeid"
            self.msm_Node = os.path.join(mike_urban_database,"msm_Node")
            self.msm_Link = os.path.join(mike_urban_database,"msm_Link")
            msm_Weir = os.path.join(mike_urban_database,"msm_Weir")
            msm_Orifice = os.path.join(mike_urban_database,"msm_Orifice")
            msm_Pump = os.path.join(mike_urban_database,"msm_Pump")
            map_only = map_only.lower()
        elif len(nodes_and_links)==2:
            self.msm_Node = nodes_and_links[0]
            self.msm_Link = nodes_and_links[1]
            fromnode_fieldname = "FROMNODE"
            tonode_fieldname = "TONODE"
            map_only = "links"
            is_sqlite = False
        else:
            raise(Exception("No MIKE Urban Database, or improper import nodes_and_links (must be list([nodes_filepath, links_filepath]))"))

        filter_sql_query = "" if not filter_sql_query or len(filter_sql_query)>2900 else filter_sql_query

        self.nodes = {}
#        print(arcpy.management.GetCount(msm_Node))
        self.points_xy = []

        self.points_muid = {}
        i = 0
        with arcpy.da.SearchCursor(self.msm_Node, ["MUID", "SHAPE@", "InvertLevel"]) as cursor:
            for row in cursor:
                if row[1] is not None:
                    self.nodes[row[0]] = self.Node(row[0], row[1], row[2] if row[2] else 0)
                    self.points_xy.append((row[1].firstPoint.X, row[1].firstPoint.Y))
                    self.points_muid[i] = row[0]
                    i += 1
        points_muid_set = set(self.points_muid.values())
        self._muid_to_index = {v: k for k, v in self.points_muid.items()}

        self.points_xy = np.array(self.points_xy)
        self.kdtree = cKDTree(self.points_xy)  # self.points_xy shape (N, 2)

        def validateNode(point, reference, search_radius = 0.1):
            distance = np.sum(reference-[point.X, point.Y])**2
            return distance < search_radius**2

        if map_only == "" or "link" in map_only:
            self.links = {}
#            getFromNodeRe = re.compile(r"(.+)l\d+")
            fields = ["MUID", "SHAPE@", 'Length', "SLOPE" if is_sqlite else "SLOPE_C", "Diameter", "uplevel", "dwlevel", fromnode_fieldname, tonode_fieldname] if fromnode_fieldname in [f.name for f in arcpy.ListFields(self.msm_Link)] else ["MUID", "SHAPE@", 'Length', "SLOPE" if is_sqlite else "SLOPE_C", "Diameter", "uplevel", "dwlevel"]
            with arcpy.da.SearchCursor(self.msm_Link, fields, where_clause = filter_sql_query) as cursor:
                fromnode_tonode_valid = True if fromnode_fieldname in fields else False
                for row in cursor:
                    if row[1] is not None:
                        self.links[row[0]] = self.Link(row[0])
                        if (fromnode_tonode_valid and row[7] and row[8] and
                                row[7] in points_muid_set and row[8] in points_muid_set and
                                validateNode(row[1].firstPoint, self.points_xy[self._muid_to_index[row[7]]]) and
                                validateNode(row[1].lastPoint, self.points_xy[self._muid_to_index[row[8]]])):
                            self.links[row[0]].fromnode = row[7]
                            self.links[row[0]].tonode = row[8]
                            self.links[row[0]].node_field_correct = True
                        else:
                            self.links[row[0]].fromnode = self.findClosestNode(row[1].firstPoint)
                            self.links[row[0]].tonode = self.findClosestNode(row[1].lastPoint)

                        self.links[row[0]].shape = row[1]
                        self.links[row[0]].length = row[2] if row[2] else row[1].length
                        self.links[row[0]].slope = row[3]
                        self.links[row[0]].diameter = row[4]

                        try:
                            self.links[row[0]].uplevel = row[5] if row[5] else self.nodes[self.links[row[0]].fromnode].invert_level
                            self.links[row[0]].dwlevel = row[6] if row[6] else self.nodes[self.links[row[0]].tonode].invert_level
                        except Exception as e:
                            print("Link %s does not have FromNode or ToNode (%s-%s)" % (row[0], self.links[row[0]].fromnode, self.links[row[0]].fromnode))
                            # raise(Exception("Link %s does not have FromNode or ToNode (%s-%s)" % (row[0], self.links[row[0]].fromnode, self.links[row[0]].fromnode)))

        if map_only == "" or "weir" in map_only:
            self.weirs = {}
            with arcpy.da.SearchCursor(msm_Weir, ["MUID", "SHAPE@"], where_clause = filter_sql_query) as cursor:
                for row in cursor:
                    self.weirs[row[0]] = self.Link(row[0])
                    self.weirs[row[0]].fromnode = self.findClosestNode(row[1].firstPoint)
                    self.weirs[row[0]].tonode = self.findClosestNode(row[1].lastPoint)
                    self.weirs[row[0]].length = row[1].length

        if map_only == "" or "pump" in map_only:
            self.pumps = {}
            with arcpy.da.SearchCursor(msm_Pump, ["MUID", "SHAPE@"], where_clause = filter_sql_query) as cursor:
                for row in cursor:
                    self.pumps[row[0]] = self.Link(row[0])
                    self.pumps[row[0]].fromnode = self.findClosestNode(row[1].firstPoint)
                    self.pumps[row[0]].tonode = self.findClosestNode(row[1].lastPoint)
                    self.pumps[row[0]].length = row[1].length

        if map_only == "" or "orifice" in map_only:
            self.orifices = {}
            with arcpy.da.SearchCursor(msm_Orifice, ["MUID", "SHAPE@"], where_clause = filter_sql_query) as cursor:
                for row in cursor:
                    self.orifices[row[0]] = self.Link(row[0])
                    self.orifices[row[0]].fromnode = self.findClosestNode(row[1].firstPoint)
                    self.orifices[row[0]].tonode = self.findClosestNode(row[1].lastPoint)
                    self.orifices[row[0]].length = row[1].length

    class Node:
        def __init__(self, MUID, shape, invertlevel):
            self.MUID = MUID
            self.shape = shape
            self.invert_level = invertlevel

    class Link:
        def __init__(self, MUID):
            self.MUID = MUID
            self._shape_3d = None

            self.fromnode = 1
            self.tonode = None
            self.length = None
            self.node_field_correct = False
            self.slope = None
            self.shape = None
            self.uplevel = None
            self.dwlevel = None

        @property
        def v_full(self):
            from .utils import calculate_full_flow
            try:
                v = calculate_full_flow(self.diameter, self.slope / 1e2, "PL") / ((self.diameter / 2) ** 2 * 3.1415)
            except Exception as e:
                v = 1
            return v

        @property
        def travel_time(self):
            return self.length / self.v_full

        def shape_3d(self, uplevel = None, dwlevel = None):
            if not uplevel:
                uplevel = self.uplevel
            if not dwlevel:
                dwlevel = self.dwlevel

            if self._shape_3d is None or (uplevel != self.uplevel and dwlevel != self.dwlevel):
                self._shape_3d = self._generate_shape_3d(uplevel, dwlevel)

            return self._shape_3d

        def _generate_shape_3d(self, uplevel, dwlevel):
            slope = (uplevel - dwlevel) / self.length

            linelist = []
            for part in self.shape:
                parts = []
                for part_i, point in enumerate(part):
                    if part_i == 0:
                        z = uplevel
                    elif part_i == len(part) - 1:
                        z = dwlevel
                    else:
                        total_distance = 0
                        point_geometries = [arcpy.PointGeometry(p) for p in part]
                        for i in range(1, part_i + 1):
                            total_distance += point_geometries[i - 1].distanceTo(point_geometries[i])
                        z = uplevel - total_distance * slope

                    parts.append(arcpy.Point(point.X, point.Y, z))
                linelist.append(parts)

            return arcpy.Polyline(arcpy.Array(linelist), None, True)

    def findClosestNode(self, point, search_radius=0.1):
        muid = None
        distance, index_closest = self.kdtree.query([point.X, point.Y], distance_upper_bound=search_radius)

        if distance < search_radius:
            muid = self.points_muid[index_closest]
        return muid

    def fixConnections(self, search_radius = 1):
        links_missing_fromnode = [link.MUID for link in self.links.values() if not link.fromnode]
        links_missing_tonode = [link.MUID for link in self.links.values() if not link.tonode]
        if "fromnodeid" in [field.name.lower() for field in arcpy.ListFields(self.msm_Link)]:
            links_missing_fromnode += [row[0] for row in arcpy.da.SearchCursor(self.msm_Link, ["MUID"], where_clause = "fromnodeid IS NULL")]
            links_missing_tonode += [row[0] for row in
                                     arcpy.da.SearchCursor(self.msm_Link, ["MUID"], where_clause="tonodeid IS NULL")]
        elif "fromnode" in [field.name.lower() for field in arcpy.ListFields(self.msm_Link)]:
            links_missing_fromnode += [row[0] for row in arcpy.da.SearchCursor(self.msm_Link, ["MUID"], where_clause = "fromnode = ''")]
            links_missing_tonode += [row[0] for row in
                                     arcpy.da.SearchCursor(self.msm_Link, ["MUID"], where_clause="tonode = ''")]


        print("Missing FromNode:")
        print(links_missing_fromnode)
        print(links_missing_tonode)

        class Point:
            def __init__(self, x, y):
                self.X = x
                self.Y = y

        if "sqlite" in arcpy.Describe(self.msm_Link).catalogPath:
            import sqlite3
            from shapely.geometry import LineString
            from shapely.wkt import loads, dumps
            conn_db1 = sqlite3.connect(os.path.dirname(arcpy.Describe(self.msm_Link).catalogPath))
            conn_db1.enable_load_extension(True)
            conn_db1.execute('SELECT load_extension("mod_spatialite")')
            # conn_db1.execute("SELECT sqlite_compileoption_used('ENABLE_RTREE')")
            # import rtree
            # conn_db1.load_extension('SQLITE_ENABLE_RTREE.dll')

            cursor = conn_db1.cursor()

            for link in links_missing_fromnode:
                cursor.execute("SELECT AsText(geometry) FROM msm_Link WHERE muid = '%s'" % (link))
                row = cursor.fetchone()
                line = loads(row[0])

                coords = list(line.coords)
                # coords[0] = self.findClosestNode(Point(coords[0][0],coords[0][1]), search_radius = search_radius)
                fromnode = self.findClosestNode(Point(coords[0][0],coords[0][1]), search_radius=search_radius)

                if not fromnode:
                    import warnings
                    warnings.warn("Could not find FromNode for Link %s with search radius of %d m" % (link, search_radius))
                else:
                    coords[0] = tuple(self.points_xy[
                        self.points_muid.index(fromnode)])

                    updated_line = LineString(coords)
                    updated_wkt = dumps(updated_line)

                    # cursor.execute(f"SELECT GeomFromText(?, 4326)", (updated_wkt,))
                    # print(row)
                    # print("GeomFromText('%s')" % (updated_wkt))
                    # cursor.execute("SELECT GeomFromText('%s')" % (updated_wkt))
                    print("UPDATE msm_Link SET geometry = GeomFromText('%s',-1) WHERE muid = '%s'" % (updated_wkt, link))
                    cursor.execute("UPDATE msm_Link SET geometry = GeomFromText('%s',-1) WHERE muid = '%s'" % (updated_wkt, link))
                    if "fromnodeid" in [field.name for field in arcpy.ListFields(self.msm_Link)]:
                        cursor.execute(
                            "UPDATE msm_Link SET fromnodeid = '%s' WHERE muid = '%s'" % (fromnode, link))

            for link in links_missing_tonode:
                cursor.execute("SELECT AsText(geometry) FROM msm_Link WHERE muid = '%s'" % (link))
                row = cursor.fetchone()
                line = loads(row[0])

                coords = list(line.coords)
                # coords[0] = self.findClosestNode(Point(coords[0][0],coords[0][1]), search_radius = search_radius)
                tonode = self.findClosestNode(Point(coords[-1][0],coords[-1][1]), search_radius=search_radius)
                if not tonode:
                    import warnings
                    warnings.warn("Could not find ToNode for Link %s with search radius of %d m" % (link, search_radius))
                else:
                    coords[-1] = tuple(self.points_xy[
                        self.points_muid.index(tonode)])
                    updated_line = LineString(coords)
                    updated_wkt = dumps(updated_line)

                    # cursor.execute(f"SELECT GeomFromText(?, 4326)", (updated_wkt,))
                    # print(row)
                    # print("GeomFromText('%s')" % (updated_wkt))
                    # cursor.execute("SELECT GeomFromText('%s')" % (updated_wkt))
                    print("UPDATE msm_Link SET geometry = GeomFromText('%s',-1) WHERE muid = '%s'" % (updated_wkt, link))
                    cursor.execute("UPDATE msm_Link SET geometry = GeomFromText('%s',-1) WHERE muid = '%s'" % (updated_wkt, link))
                    if "fromnodeid" in [field.name for field in arcpy.ListFields(self.msm_Link)]:
                        cursor.execute(
                            "UPDATE msm_Link SET tonodeid = '%s' WHERE muid = '%s'" % (tonode, link))
            conn_db1.commit()
            del cursor
        else:
            edit = arcpy.da.Editor(self.mike_urban_database)
            edit.startEditing(False, True)
            edit.startOperation()
            if links_missing_fromnode:
                with arcpy.da.UpdateCursor(self.msm_Link, ["MUID", "SHAPE@"], where_clause = "MUID IN ('%s')" % "', '".join(links_missing_fromnode)) as cursor:
                    for row in cursor:
                        points = [arcpy.Point(p.X, p.Y) for p in row[1].getPart(0)]
                        closest_node = self.findClosestNode(row[1].firstPoint, search_radius=search_radius)
                        if closest_node:
                            points[0] = arcpy.Point(*self.points_xy[
                                self.points_muid.index((closest_node))])
                            row[1] = arcpy.Polyline(arcpy.Array(points))
                            cursor.updateRow(row)
                        else:
                            print("Could not find a fromnode for link %s" % (row[0]))
            if links_missing_tonode:
                with arcpy.da.UpdateCursor(self.msm_Link, ["MUID", "SHAPE@"], where_clause = "MUID IN ('%s')" % "', '".join(links_missing_tonode)) as cursor:
                    for row in cursor:
                        points = [arcpy.Point(p.X, p.Y) for p in row[1].getPart(0)]
                        closest_node = self.findClosestNode(row[1].lastPoint, search_radius=search_radius)
                        if closest_node:
                            points[-1] = arcpy.Point(*self.points_xy[
                                self.points_muid.index((closest_node))])
                            row[1] = arcpy.Polyline(arcpy.Array(points))
                            cursor.updateRow(row)
                        else:
                            print("Could not find a tonode for link %s" % (row[0]))
            edit.stopOperation()
            edit.stopEditing(True)


if __name__ == "__main__":
    PipeNetwork("model.mdb")