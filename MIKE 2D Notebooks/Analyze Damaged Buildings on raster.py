import arcpy
import numpy as np
from arcpy import env
import os

# Set environment settings
env.overwriteOutput = True

# Define input and output raster paths
water_depth_raster_path = r"C:\Papirkurv\VBT_244_Detailforslag\VBT_244_Detailforslag_12062024_Aarhus_C_depth.tif"
flux_raster_path = r"C:\Papirkurv\VBT_244_Detailforslag\VBT_244_Detailforslag_12062024_Aarhus_C_flux.tif"
buildings_path = r"C:\Papirkurv\Bygninger\dataset.shp"
cadastre_path = r"C:\Papirkurv\Jordstykker\dataset.shp"
building_field = "bbruuid"

erhvervsenhedspris = 2095 # DKK/m²
privatenhedspris = 1257 # DKK/m²
vejpris = 3 # DKK/m²

buidings_buffer_size = 0.5
buildings_damaged_length = 1.3
buildings_damaged_area = buidings_buffer_size*buildings_damaged_length

critical_water_depth_output_path = r"in_memory\critical_water_depth_processed.tif"
critical_flux_output_path = r"in_memory\critical_flux_processed.tif"

tracing_water_depth_output_path = r"in_memory\tracing_water_depth_processed.tif"
tracing_flux_output_path = r"in_memory\tracing_flux_processed.tif"

critical_mosaic_output_path = r"C:\Papirkurv\Critical_mosaic.tif"
tracing_mosaic_output_path = r"C:\Papirkurv\Tracing_mosaic.tif"
critical_mosaic_polygon_output_path = r"C:\Papirkurv\critical_mosaic_polygon.shp"
buildings_buffer_path = r"in_memory\buildings_buffer"
buildings_priced_path = r"C:\Papirkurv\priced_buildings.shp"
damage_area_path = r"C:\Papirkurv\damage_area.shp"
buildings_damaged_path = r"C:\Papirkurv\damaged_buildings.shp"
roads_path = r"in_memory\roads"
roads_damaged_path = r"C:\Papirkurv\damaged_roads.shp"

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


# analyze critical water depth and flux
process_raster(water_depth_raster_path, critical_water_depth_output_path, 0.1)
process_raster(flux_raster_path, critical_flux_output_path, 25e-3)

# analyze tracing water depth and flux
process_raster(water_depth_raster_path, tracing_water_depth_output_path, 5e-3)
process_raster(flux_raster_path, tracing_flux_output_path, 0.2)

arcpy.MosaicToNewRaster_management(
        input_rasters=[critical_water_depth_output_path, critical_flux_output_path],
        output_location=os.path.dirname(critical_mosaic_output_path),
        raster_dataset_name_with_extension=os.path.basename(critical_mosaic_output_path),
        pixel_type="8_BIT_UNSIGNED",
        cellsize=cell_size,
        number_of_bands=1,
        mosaic_method="MAXIMUM"
    )

arcpy.MosaicToNewRaster_management(
        input_rasters=[tracing_water_depth_output_path, tracing_flux_output_path],
        output_location=os.path.dirname(tracing_mosaic_output_path),
        raster_dataset_name_with_extension=os.path.basename(tracing_mosaic_output_path),
        pixel_type="8_BIT_UNSIGNED",
        cellsize=cell_size,
        number_of_bands=1,
        mosaic_method="MAXIMUM"
    )

# Creating Buffer around buildings to search for nearby water depth and flux
arcpy.Buffer_analysis(in_features=buildings_path, out_feature_class=buildings_buffer_path, buffer_distance_or_field="%d Millimeters" % (buidings_buffer_size*1e3), line_side="OUTSIDE_ONLY", line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="PLANAR")

# Converting the critical flood and flux to polygon, to use for clipping
arcpy.RasterToPolygon_conversion(in_raster=critical_mosaic_output_path, out_polygon_features=critical_mosaic_polygon_output_path, simplify="NO_SIMPLIFY", raster_field="Value", create_multipart_features="MULTIPLE_OUTER_PART", max_vertices_per_feature="")
arcpy.RepairGeometry_management(critical_mosaic_polygon_output_path)

# Deleting the area where there's no flood or flux above limit
with arcpy.da.UpdateCursor(critical_mosaic_polygon_output_path, ["gridcode"], where_clause = "gridcode = 0") as cursor:
    for row in cursor:
        cursor.deleteRow()


### BUILDINGS ###
# Clipping for buildings- used for finding areaas where critical flood and flux overlap with the buffer of the buildings
arcpy.Clip_analysis(in_features=buildings_buffer_path, clip_features=critical_mosaic_polygon_output_path, out_feature_class=damage_area_path, cluster_tolerance="")

# Repair Geometry of clip
arcpy.RepairGeometry_management(buildings_buffer_path)

# Deleting all that are less than limit of area, saving all above to list
damaged_buildings = []
with arcpy.da.UpdateCursor(damage_area_path, ["SHAPE@AREA", building_field]) as cursor:
    for row in cursor:
        if row[0] < buildings_damaged_area:
            cursor.deleteRow()
        elif row[1] and not row[1] == ' ':
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

arcpy.Select_analysis(in_features=buildings_priced_path, out_feature_class=buildings_damaged_path, where_clause="%s IN ('%s')" % (building_field, "', '".join(damaged_buildings)))

print("%d DKK" % (np.sum([row[0] for row in arcpy.da.SearchCursor(buildings_damaged_path, ["pris"], where_clause =
                                                                            "bbruuid IN ('744f7896-d9eb-4209-a3b6-14521a77db71', 'f9af46fb-eb48-4582-afee-fbeb153e0173', 'fbfa755f-940e-4f3e-9ffc-01a1d1759185', '2e37b4ce-c12f-4c8a-8b39-951f56ede430', '5a9dc789-2aca-4ede-bcca-332e73f4a7f9', '38995303-c4bc-4fdb-af90-983587e8e79d', '24f9bba7-2c51-4f64-b45b-28628732e77b', 'e3514705-d199-4687-bf78-5970028c6c6e', 'd5eb4921-c578-4dc9-a454-e3e8b8da9b26', '6e24b8e8-cf50-4cf8-96ae-b1d45b76401f', '71ccaafe-15eb-439a-8097-5f2982c50957', 'e569ec06-06b3-4ff9-b707-38f322d1b836', 'af27b999-4e23-4f70-aa3f-c71f9ab615e8', 'e42a0bd1-b9ae-4d9e-8abe-bcb66371c1e4', 'a81989b2-1290-42c5-8cde-6f718e3dc61a', '9d08aed6-c2b2-4186-8b91-62f321aae396', 'bc8c92ac-3c43-4a19-8b35-536d3fca92c7', 'f6c3e6d8-34ae-4ab3-a38f-d2738e3c3875', '96165878-279a-4aab-8180-1c4a9e85969c', '3c4cd372-f2bd-4b22-9b9a-d2d5db679612', '35e23204-5995-4532-aa74-0277e2dc4538', '27c7724a-d9ec-4c04-8a77-25e3c154e210', '8d60c348-a235-4a37-bead-c507ebf73ee4', 'c03d43c6-c758-4519-8738-4782b52f67e5', '8e8d830a-4644-43fc-85e7-972c7470c74c', 'f0d58c52-af6d-40a3-b396-b38792cce038', 'dfffc038-1a77-402b-bf61-1df23872e826', '8f1c3320-8ff5-44d0-9db5-67a7eea43fc0', 'd546a922-2f01-4c14-a8f0-d88fed838297', '08c942f7-0f98-4eb2-ba34-8e1dae5ae515', 'f3247cd3-2c6f-4e39-b2c3-9ac839610c9c', '94f3e2f7-40f8-4782-a202-ed51cf13a90e', '2644a9db-3df9-4f73-9bda-d133e16e12c6', '05ee267c-e2b1-4832-a442-db744cd37176', '95d01aef-3627-4a95-b08e-9689b609c05a', '704e5235-b24a-4f5f-a0bd-7cd5aee7d46a', '8df45842-79f3-4b0e-a0b9-12f390d291c6', '876cc666-668f-49ec-82d6-bdc78f385a78', '7c7cd0b1-2936-4c60-a8d3-409b4f8d6a6f', '055ce1bc-2bbc-4219-81bc-b64b2c4e6b08', '3b60a828-f50d-4c9b-b45f-6192bf097324', '9af435dc-a0a2-4dd8-baa1-7e35b52a637c', '891bf6ff-aa21-41e0-9cce-c8e068d6c059', 'b97585ba-b78f-4b68-af62-9fbaa356ae9e', 'e7ce62a6-e775-4bab-b5ac-2d1e64ce9b7b', '72d5abbe-283e-43e5-8bd2-3f4e07299c36', '4a4c5003-747b-42e6-b539-1dd31ae00012', '376ff359-3c7c-4743-8b48-4494fe8ddd02', '74998652-9d07-4a5d-9b2f-cdf4dc178f69', '9498ab35-3cda-46de-8c88-95a7e9891d2f', '3e1e39ad-fc6f-4bcd-a640-8caa7a64b2d2', '8fdb4b6c-7ae1-410e-8f6c-dbfb9cc262e8', '58ae140e-057b-4292-918a-b1ea36b1a231', 'e5fca670-9a19-4c51-89f7-78432dd3d05a', '40e65df9-4669-4b3f-91da-509369ba55f3', '6a97442c-89a8-4a26-9638-bca9f4c97260', '601fed9f-8125-44ea-b9ce-f4d84b75f1cd', '27822eb3-5526-49f9-9c86-e10fd256d13d', '711f53be-1d5f-4f98-9135-5cd26d78cfe3', '361c4223-3e2e-4e6b-9cbe-7fd32835ee07', 'b0a6c739-3f1c-4615-a154-fcc7f965327b', '4a20fabd-85ac-44c5-921f-5ffc2005aca3', '20bfd842-b1b3-4ad0-8d5a-6cbb41e214fe', 'dc85c306-9931-43dc-9fe0-6d61d248a33f', 'd85fe172-af4f-4703-a083-7f82db564534', '25569538-37bf-4197-93b2-2f85fdebe2b9', 'cd6cb193-8274-4c3b-af67-fba0e49604d7', '6aa5870e-81d8-491b-85fb-835194cedd26', '0ef97f7a-6051-4341-9ba1-1278360deb5a', '5e51d158-38ed-491e-b3c5-a63fb3eed01f', '53275866-b262-48dd-9204-7b930ed26d21', '2d600412-0e43-464e-9a2d-e688038e6650', '8895fd2a-ac73-4b08-87fe-fefdf5e5a5d3', '6401fe49-427f-4340-a358-90e35ce5e5ee', 'fcd6958e-c358-4066-aa64-efb6b57e78c9', '3cc6b39e-aabf-4fe6-ac19-bbeb97e443e1', 'd4522510-f74b-4a52-b24d-ed368c5bb50a', 'e1ed269d-fbed-4eac-88e5-1f463ad9b417', '74b66afe-110a-4c6b-9475-b6e2f26f4550', 'b7b0f817-9572-4798-9793-6ed3b15af815', '69a90216-8489-4f6b-ad3a-50873bbaa8c1', 'd8dfe2d6-3f8f-43c2-aae4-45561b912ed1')"
                                                                            )])))


### CADASTRE ###
# Clipping for buildings- used for finding areas where critical flood and flux overlap with the buffer of the buildings
arcpy.Select_analysis(in_features=cadastre_path, out_feature_class=roads_path, where_clause="vejareal>0")
arcpy.Clip_analysis(in_features=roads_path, clip_features=critical_mosaic_polygon_output_path, out_feature_class=roads_damaged_path, cluster_tolerance="")

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