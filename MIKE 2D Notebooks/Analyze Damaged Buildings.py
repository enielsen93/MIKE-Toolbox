import arcpy
import numpy as np
from arcpy import env
import os
import mikeio
from alive_progress import alive_bar
import math

# Set environment settings
env.overwriteOutput = True

# Define input and output raster paths
flood_dfsu_path = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\05_RESULTS\03_FLOOD\Status\2024-06-10\VBT_219_mGA_CDS_20_138_FMBaseDefault_2D_overland_reduced.dfsu"
buildings_path = r"C:\Papirkurv\Bygninger\dataset.shp"
cadastre_path = r"C:\Papirkurv\Jordstykker\dataset.shp"
building_field = "bbruuid"

MIKE_model = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\VBT_STATUS_011\VBT_STATUS_011.sqlite"
msm_Node = os.path.join(MIKE_model, "msm_Node")
nodes_output_path = "C:\Papirkurv\msm_Node"

arcpy.CreateFeatureclass_management(
    out_path=os.path.dirname(nodes_output_path),
    out_name=os.path.basename(nodes_output_path),
    geometry_type="POINT"
)

with arcpy.da.InsertCursor(nodes_output_path, ["SHAPE@"]) as insertcursor:
    with arcpy.da.SearchCursor(msm_Node, ["SHAPE@"]) as cursor:
        for row in cursor:
            insertcursor.insertRow(row)
arcpy.management.AddSpatialIndex(nodes_output_path)

erhvervsenhedspris = 2095 # DKK/m²
privatenhedspris = 1257 # DKK/m²
vejpris = 3 # DKK/m²

buidings_buffer_size = 0.5
buildings_damaged_length = 1.3
buildings_damaged_area = buidings_buffer_size*buildings_damaged_length

critical_output_path = r"C:\Papirkurv\critical_processed"
tracing_output_path = r"C:\Papirkurv\tracing_processed"
tracing_dissolved_output_path = r"C:\Papirkurv\tracing_dissolved_processed"

buildings_buffer_path = r"C:\Papirkurv\buildings_buffer.shp"
buildings_priced_path = r"C:\Papirkurv\priced_buildings.shp"
damage_area_path = r"C:\Papirkurv\damage_area.shp"
damage_area_filtered_path = r"C:\Papirkurv\damage_area_filtered.shp"
buildings_damaged_path = r"C:\Papirkurv\damaged_buildings.shp"
roads_path = r"in_memory\roads"
roads_damaged_path = r"C:\Papirkurv\damaged_roads.shp"

qflux_output_path = r"C:\Papirkurv\qflux.tif"
pflux_output_path = r"C:\Papirkurv\pflux.tif"
vector_field_output_path = r"C:\Papirkurv\vector_field.tif"



global cell_size

def process_raster(input_raster_path, output_raster_path, limit):
    global cell_size
    # Load raster into a NumPy array
    input_raster = arcpy.Raster(input_raster_path)
    raster_array = arcpy.RasterToNumPyArray(input_raster)

    # Set values above 0.1 to 1 and below or equal to 0.1 to 0
    raster_array = np.where(raster_array > limit, 1, 0)

    # Save the modified array to a new raster
    lower_left = arcpy.Point(input_raster.extent.XMin, input_raster.extent.YMin)
    cell_size = input_raster.meanCellWidth

    # Create output raster from the NumPy array
    output_raster = arcpy.NumPyArrayToRaster(raster_array, lower_left, cell_size, value_to_nodata=0)

    # Save the output raster
    output_raster.save(output_raster_path)

def process_dfsu(dfsu_path, output_polygon_path, limit_flood = None, limit_flux = None, limit_flux_volume = None):
    limits = [limit_flood, limit_flux]
    dfs = mikeio.dfsu.Dfsu2DH(dfsu_path)

    output_polygon = \
    arcpy.CreateFeatureclass_management(os.path.dirname(output_polygon_path), os.path.basename(output_polygon_path),
                                        "POLYGON", spatial_reference=dfs.projection_string)[0]
    # arcpy.AddField_management(output_polygon, "Critical", "STRING")
    arcpy.AddField_management(output_polygon, "Depth", "TEXT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "Flux", "TEXT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "Code", "SHORT")
    arcpy.AddField_management(output_polygon, "Public", "SHORT")

    if any([True for item in dfs.items if "Maximum" in item.name]):
        print("Assuming DFSU file is flood statistics")
        dfsu_items = ["Maximum water depth", "Maximum current speed"]
        dfs_read = dfs.read(items=[i for i, item in enumerate(dfs.items) if any([field.lower() in item.name.lower() for field in dfsu_items])])
        dfs_read_data = dfs_read.to_numpy()
        np.nan_to_num(dfs_read_data, copy=False)
        dfs_data_maximum = dfs_read_data
        maximum_flux = dfs_data_maximum[-1,0,:]
        maximum_depth = dfs_read_data[0,-1,:]
    else:
        print("Assuming DFSU file is dynamic")
        dfsu_items = ["Total water depth", "P Flux", "Q FLux"]
        dfs_read = dfs.read(
            items=[i for i, item in enumerate(dfs.items) if any([field.lower() in item.name.lower() for field in dfsu_items])])
        dfs_read_data = dfs_read.to_numpy()
        np.nan_to_num(dfs_read_data, copy=False)
        flux = np.sqrt(dfs_read_data[1,:,:]**2 + dfs_read_data[2,:,:]**2)
        maximum_flux = np.amax(flux, axis=0)
        maximum_depth = np.amax(dfs_read_data[0,:,:], 0)

    elements_above_limit = set(list(np.where(maximum_depth > limit_flood)[0]) + list(np.where(maximum_flux > limit_flux)[0]))
    element_table = dfs.element_table
    nodes_coordinates = dfs.node_coordinates
    codes = dfs.geometry.codes

    with alive_bar(len(elements_above_limit), force_tty=True) as bar:
        with arcpy.da.InsertCursor(output_polygon, ["SHAPE@", "Depth", "Flux", "Code"]) as cursor:
            for element in elements_above_limit:
                node_coordinates = nodes_coordinates[element_table[element], :-1]
                code = np.max(codes[element_table[element]])
                triangle = arcpy.Polygon(
                    arcpy.Array([arcpy.Point(coords[0], coords[1]) for coords in node_coordinates]))
                # arcpy.AddMessage(dfs_read_data.shape)
                # arcpy.AddMessage(dfs_read_data[0, -1, element])
                cursor.insertRow([triangle, maximum_depth[element], maximum_flux[element], code])
                bar()

    # with alive_bar(len(elements_with_water), force_tty=True) as bar:
    #     for element_i, element in enumerate(elements_above_limit):
    # for item_i, item in enumerate(dfsu_items):
    #    pass
    arcpy.management.AddSpatialIndex(output_polygon)

def flowArrows(dfsu_path, output_flux_path):
    dfs = mikeio.dfsu.Dfsu2DH(dfsu_path)

    output_polygon = \
    arcpy.CreateFeatureclass_management(os.path.dirname(output_flux_path), os.path.basename(output_flux_path),
                                        "POLYGON", spatial_reference=dfs.projection_string)[0]
    # arcpy.AddField_management(output_polygon, "Critical", "STRING")
    arcpy.AddField_management(output_polygon, "Depth", "FLOAT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "QFlux", "FLOAT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "PFlux", "FLOAT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "TotFlux", "FLOAT", field_precision=6, field_scale=3)
    arcpy.AddField_management(output_polygon, "Direction", "SHORT")

    print("Assuming DFSU file is dynamic")
    dfsu_items = ["Total water depth", "P Flux", "Q FLux"]
    dfs_read = dfs.read(
        items=[i for i, item in enumerate(dfs.items) if any([field.lower() in item.name.lower() for field in dfsu_items])])
    dfs_read_data = dfs_read.to_numpy()
    np.nan_to_num(dfs_read_data, copy=False)
    timesteps = np.diff([time.value / 1e9 for time in dfs.time])
    timesteps = timesteps.reshape(timesteps.shape[0], 1)
    total_flux = np.sum(timesteps * (np.sqrt(dfs_read_data[1,1:,:]**2 + dfs_read_data[2,1:,:]**2)),axis=0)
    mean_P_flux = np.mean(dfs_read_data[1, :, :], axis=0)
    mean_Q_flux = np.mean(dfs_read_data[2, :, :], axis=0)
    maximum_depth = np.amax(dfs_read_data[0, :, :], 0)

    elements_above_limit = set(list(np.where(total_flux > 0.1)[0]) + list(np.where(maximum_depth > 0.05)[0]))
    element_table = dfs.element_table
    nodes_coordinates = dfs.node_coordinates
    codes = dfs.geometry.codes

    with alive_bar(len(elements_above_limit), force_tty=True) as bar:
        with arcpy.da.InsertCursor(output_polygon, ["SHAPE@", "Depth", "QFlux", "PFlux", "TotFlux", "Direction"]) as cursor:
            for element in list(elements_above_limit):
                node_coordinates = nodes_coordinates[element_table[element], :-1]
                code = np.max(codes[element_table[element]])
                triangle = arcpy.Polygon(
                    arcpy.Array([arcpy.Point(coords[0], coords[1]) for coords in node_coordinates]))
                # arcpy.AddMessage(dfs_read_data.shape)
                # arcpy.AddMessage(dfs_read_data[0, -1, element])
                # print(maximum_depth[element])
                cursor.insertRow([triangle, float(maximum_depth[element]), float(mean_P_flux[element]), float(mean_Q_flux[element]), float(total_flux[element]), math.degrees(math.atan2(mean_Q_flux[element], mean_P_flux[element]))])
                bar()

    arcpy.management.AddSpatialIndex(output_polygon)


# analyze critical water depth and flux
process_dfsu(flood_dfsu_path, critical_output_path, 0.1, 25e-3)

# tracing water
flowArrows(flood_dfsu_path, tracing_output_path)

arcpy.Dissolve_management(
    in_features=tracing_output_path,
    out_feature_class=tracing_dissolved_output_path,
    multi_part="SINGLE_PART"  # Ensure no multipart features
)

spatial_join = arcpy.analysis.SpatialJoin(tracing_dissolved_output_path, nodes_output_path, "in_memory\spatial_join", join_type = "KEEP_COMMON", match_option = "CONTAINS")[0]

spatial_join2 = arcpy.analysis.SpatialJoin(critical_output_path, spatial_join, "in_memory\spatial_join2", join_type = "KEEP_COMMON", match_option = "WITHIN")[0]

FID = [row[0] for row in arcpy.da.SearchCursor(spatial_join2, ["TARGET_FID_1"])]

with arcpy.da.UpdateCursor(critical_output_path, ["FID", "Public"]) as updcursor:
    for row in updcursor:
        if row[0] in FID:
            row[1] = 1
            updcursor.updateRow(row)

arcpy.CheckOutExtension("Spatial")
arcpy.PolygonToRaster_conversion(tracing_output_path, "QFlux", qflux_output_path, "CELL_CENTER", cellsize = 0.4)
arcpy.PolygonToRaster_conversion(tracing_output_path, "PFlux", pflux_output_path, "CELL_CENTER", cellsize = 0.4)
arcpy.CheckInExtension("Spatial")
arcpy.CompositeBands_management(in_rasters=[qflux_output_path, pflux_output_path], out_raster=vector_field_output_path)

# Creating Buffer around buildings to search for nearby water depth and flux
arcpy.Buffer_analysis(in_features=buildings_path, out_feature_class=buildings_buffer_path, buffer_distance_or_field="%d Millimeters" % (buidings_buffer_size*1e3), line_side="OUTSIDE_ONLY", line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="PLANAR")

### BUILDINGS ###
# Clipping for buildings- used for finding areas where critical flood and flux overlap with the buffer of the buildings
# arcpy.Clip_analysis(in_features=buildings_buffer_path, clip_features=critical_output_path, out_feature_class=damage_area_path, cluster_tolerance="")
arcpy.Intersect_analysis(in_features=[buildings_buffer_path, critical_output_path], out_feature_class=damage_area_path, join_attributes="ALL", output_type="INPUT")

arcpy.Dissolve_management(in_features=damage_area_path, out_feature_class=damage_area_filtered_path, dissolve_field=["bbruuid", "public"], statistics_fields="FID_buildi COUNT", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")

# Repair Geometry of clip
arcpy.RepairGeometry_management(damage_area_filtered_path)

# Deleting all that are less than limit of area, saving all above to list
damaged_buildings = []
damaged_buildings_public = []
with arcpy.da.UpdateCursor(damage_area_filtered_path, ["COUNT_FID_", building_field, "public"]) as cursor:
    for row in cursor:
        if row[0] < 3:
            cursor.deleteRow()
        elif row[1] and not row[1] == ' ':
            if row[2] == 1:
                damaged_buildings_public.append(row[1])
            else:
                damaged_buildings.append(row[1])

arcpy.management.CopyFeatures(buildings_path, buildings_priced_path)

# Exporting buildings that are damaged
for field in ["enhpris", "pris"]:
    arcpy.AddField_management(in_table=buildings_priced_path, field_name=field, field_type="FLOAT")

with arcpy.da.UpdateCursor(buildings_priced_path, [building_field, "enhpris", "pris", 'byg040Bygn', 'byg039Bygn', 'byg041Beby']) as cursor:
    for row in cursor:
        erhvervsareal, boligareal, fodaftryk = [row[3], row[4], row[5]]
        if erhvervsareal < fodaftryk:
            row[1] = (erhvervsareal/fodaftryk) * erhvervsenhedspris + ((fodaftryk - erhvervsareal)/fodaftryk) * privatenhedspris
        else:
            row[1] = erhvervsenhedspris

        row[2] = row[1] * fodaftryk
        cursor.updateRow(row)

arcpy.management.CopyFeatures(in_features=buildings_priced_path, out_feature_class=buildings_damaged_path)#
arcpy.AddField_management(in_table=buildings_damaged_path, field_name="Public", field_type="SHORT")

with arcpy.da.UpdateCursor(buildings_damaged_path, [building_field, "Public"]) as cursor:
    for row in cursor:
        if row[0] in damaged_buildings_public:
            row[1] = 1
            cursor.updateRow(row)
        elif row[0] in damaged_buildings:
            row[1] = 0
            cursor.updateRow(row)
        else:
            cursor.deleteRow()

print("%d DKK" % (np.sum([row[0] for row in arcpy.da.SearchCursor(buildings_damaged_path, ["pris"], where_clause =
                                                                            "bbruuid IN ('744f7896-d9eb-4209-a3b6-14521a77db71', 'f9af46fb-eb48-4582-afee-fbeb153e0173', 'fbfa755f-940e-4f3e-9ffc-01a1d1759185', '2e37b4ce-c12f-4c8a-8b39-951f56ede430', '5a9dc789-2aca-4ede-bcca-332e73f4a7f9', '38995303-c4bc-4fdb-af90-983587e8e79d', '24f9bba7-2c51-4f64-b45b-28628732e77b', 'e3514705-d199-4687-bf78-5970028c6c6e', 'd5eb4921-c578-4dc9-a454-e3e8b8da9b26', '6e24b8e8-cf50-4cf8-96ae-b1d45b76401f', '71ccaafe-15eb-439a-8097-5f2982c50957', 'e569ec06-06b3-4ff9-b707-38f322d1b836', 'af27b999-4e23-4f70-aa3f-c71f9ab615e8', 'e42a0bd1-b9ae-4d9e-8abe-bcb66371c1e4', 'a81989b2-1290-42c5-8cde-6f718e3dc61a', '9d08aed6-c2b2-4186-8b91-62f321aae396', 'bc8c92ac-3c43-4a19-8b35-536d3fca92c7', 'f6c3e6d8-34ae-4ab3-a38f-d2738e3c3875', '96165878-279a-4aab-8180-1c4a9e85969c', '3c4cd372-f2bd-4b22-9b9a-d2d5db679612', '35e23204-5995-4532-aa74-0277e2dc4538', '27c7724a-d9ec-4c04-8a77-25e3c154e210', '8d60c348-a235-4a37-bead-c507ebf73ee4', 'c03d43c6-c758-4519-8738-4782b52f67e5', '8e8d830a-4644-43fc-85e7-972c7470c74c', 'f0d58c52-af6d-40a3-b396-b38792cce038', 'dfffc038-1a77-402b-bf61-1df23872e826', '8f1c3320-8ff5-44d0-9db5-67a7eea43fc0', 'd546a922-2f01-4c14-a8f0-d88fed838297', '08c942f7-0f98-4eb2-ba34-8e1dae5ae515', 'f3247cd3-2c6f-4e39-b2c3-9ac839610c9c', '94f3e2f7-40f8-4782-a202-ed51cf13a90e', '2644a9db-3df9-4f73-9bda-d133e16e12c6', '05ee267c-e2b1-4832-a442-db744cd37176', '95d01aef-3627-4a95-b08e-9689b609c05a', '704e5235-b24a-4f5f-a0bd-7cd5aee7d46a', '8df45842-79f3-4b0e-a0b9-12f390d291c6', '876cc666-668f-49ec-82d6-bdc78f385a78', '7c7cd0b1-2936-4c60-a8d3-409b4f8d6a6f', '055ce1bc-2bbc-4219-81bc-b64b2c4e6b08', '3b60a828-f50d-4c9b-b45f-6192bf097324', '9af435dc-a0a2-4dd8-baa1-7e35b52a637c', '891bf6ff-aa21-41e0-9cce-c8e068d6c059', 'b97585ba-b78f-4b68-af62-9fbaa356ae9e', 'e7ce62a6-e775-4bab-b5ac-2d1e64ce9b7b', '72d5abbe-283e-43e5-8bd2-3f4e07299c36', '4a4c5003-747b-42e6-b539-1dd31ae00012', '376ff359-3c7c-4743-8b48-4494fe8ddd02', '74998652-9d07-4a5d-9b2f-cdf4dc178f69', '9498ab35-3cda-46de-8c88-95a7e9891d2f', '3e1e39ad-fc6f-4bcd-a640-8caa7a64b2d2', '8fdb4b6c-7ae1-410e-8f6c-dbfb9cc262e8', '58ae140e-057b-4292-918a-b1ea36b1a231', 'e5fca670-9a19-4c51-89f7-78432dd3d05a', '40e65df9-4669-4b3f-91da-509369ba55f3', '6a97442c-89a8-4a26-9638-bca9f4c97260', '601fed9f-8125-44ea-b9ce-f4d84b75f1cd', '27822eb3-5526-49f9-9c86-e10fd256d13d', '711f53be-1d5f-4f98-9135-5cd26d78cfe3', '361c4223-3e2e-4e6b-9cbe-7fd32835ee07', 'b0a6c739-3f1c-4615-a154-fcc7f965327b', '4a20fabd-85ac-44c5-921f-5ffc2005aca3', '20bfd842-b1b3-4ad0-8d5a-6cbb41e214fe', 'dc85c306-9931-43dc-9fe0-6d61d248a33f', 'd85fe172-af4f-4703-a083-7f82db564534', '25569538-37bf-4197-93b2-2f85fdebe2b9', 'cd6cb193-8274-4c3b-af67-fba0e49604d7', '6aa5870e-81d8-491b-85fb-835194cedd26', '0ef97f7a-6051-4341-9ba1-1278360deb5a', '5e51d158-38ed-491e-b3c5-a63fb3eed01f', '53275866-b262-48dd-9204-7b930ed26d21', '2d600412-0e43-464e-9a2d-e688038e6650', '8895fd2a-ac73-4b08-87fe-fefdf5e5a5d3', '6401fe49-427f-4340-a358-90e35ce5e5ee', 'fcd6958e-c358-4066-aa64-efb6b57e78c9', '3cc6b39e-aabf-4fe6-ac19-bbeb97e443e1', 'd4522510-f74b-4a52-b24d-ed368c5bb50a', 'e1ed269d-fbed-4eac-88e5-1f463ad9b417', '74b66afe-110a-4c6b-9475-b6e2f26f4550', 'b7b0f817-9572-4798-9793-6ed3b15af815', '69a90216-8489-4f6b-ad3a-50873bbaa8c1', 'd8dfe2d6-3f8f-43c2-aae4-45561b912ed1')"
                                                                            )])))


### CADASTRE ###
# Clipping for buildings- used for finding areas where critical flood and flux overlap with the buffer of the buildings
arcpy.Select_analysis(in_features=cadastre_path, out_feature_class=roads_path, where_clause="vejareal>0")
arcpy.Clip_analysis(in_features=roads_path, clip_features=critical_output_path, out_feature_class=roads_damaged_path, cluster_tolerance="")

# Repair Geometry of clip
arcpy.RepairGeometry_management(roads_damaged_path)

for field in ["enhpris", "pris"]:
    arcpy.AddField_management(in_table=roads_damaged_path, field_name=field, field_type="FLOAT")

with arcpy.da.UpdateCursor(roads_damaged_path, ["vejareal", "registrere", "enhpris", "pris", "SHAPE@AREA"]) as cursor:
    for row in cursor:
        vejareal, total_areal = [row[0], row[1]]
        if vejareal > total_areal:
            vejareal = total_areal

        if total_areal > 0:
            row[2] = vejpris * (vejareal/total_areal)

        row[3] = row[2] * row[4]
        cursor.updateRow(row)