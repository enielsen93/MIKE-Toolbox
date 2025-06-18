import arcpy
import numpy as np
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import os
import mikeio
import math

# Inputs
polyline_path = r"C:\Papirkurv\Trashcan.gdb\Cross_Section"
raster_path = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jyllands Alle\MIKE_REGNVAND\04_DTM\Terræn_Højbjerg.tif"
mesh_path = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jyllands Alle\MIKE_REGNVAND\07_2D\JYL_072_v090_Z.mesh"

# Load raster
raster = arcpy.Raster(raster_path)
arr = arcpy.RasterToNumPyArray(raster, nodata_to_value=np.nan)
left = raster.extent.XMin
top = raster.extent.YMax
cell_size = raster.meanCellWidth  # assumes square cells

dfs = mikeio.dfsu.Mesh(mesh_path)
node_coordinates = dfs.node_coordinates
element_coordinates = dfs.element_coordinates
def interp(coords):
    try:
        return element_coordinates[dfs.geometry.find_index(coords = coords), 2]
    except Exception as e:
        return np.nan

# Coordinate transform: world -> array row/col
def world_to_array(x, y):
    col = int((x - left) / cell_size)
    row = int((top - y) / cell_size)
    return row, col

# Create output folder for plots
output_folder = "C:\Papirkurv"
os.makedirs(output_folder, exist_ok=True)

# Read and process each polyline
with arcpy.da.SearchCursor(polyline_path, ["OID@", "SHAPE@", "curb_level", "slope"]) as cursor:
    for oid, geom, curb_level, slope in cursor:
        coords = [(p.X, p.Y) for part in geom for p in part if p]
        if len(coords) < 2:
            continue  # Skip degenerate

        # Interpolate at 0.2 m resolution
        line = LineString(coords)
        length = line.length
        resolution = 0.2
        num_points = max(2, int(length / resolution))
        distances = np.linspace(0, length, num_points)
        interp_points = [line.interpolate(d) for d in distances]
        xy_coords = [(pt.x, pt.y) for pt in interp_points]

        # Sample elevation
        elevations = []
        for x, y in xy_coords:
            row, col = world_to_array(x, y)
            try:
                val = arr[row, col]
            except IndexError:
                val = np.nan
            elevations.append(val)

        # If curb_level is missing, determine it from levees
        if True:# curb_level is None or (isinstance(curb_level, float) and math.isnan(curb_level)):
            N = 5  # Number of points at each end
            left_max = np.nanmax(elevations[:N])
            right_max = np.nanmax(elevations[-N:])
            curb_level = min(left_max, right_max)
            # print(f"OID {oid} | Calculated curb_level = {curb_level:.2f} m from levees")

        # Clip elevations to below curb_level
        clipped_z = np.minimum(elevations, curb_level)

        # Compute horizontal distances between points (assume constant spacing)
        dx = np.diff(distances)
        dz = np.diff(clipped_z)

        # Cross-sectional area using trapezoidal integration
        area = np.sum(0.5 * dx * (clipped_z[:-1] + clipped_z[1:] - 2 * curb_level) * -1)  # negative area becomes positive
        area = max(area, 0)

        # Wetted perimeter (2D length of submerged profile)
        submerged = np.logical_and(clipped_z[:-1] < curb_level, clipped_z[1:] < curb_level)
        perimeter = np.sum(np.sqrt(dx[submerged] ** 2 + dz[submerged] ** 2))

        # Manning calculation
        n = 1 / 70
        if perimeter > 0:
            R = area / perimeter
            q_full = (1 / n) * area * (R ** (2 / 3)) * (slope ** 0.5)
        else:
            q_full = 0.0

        import matplotlib.patches as patches

        # Start plot
        plt.figure(figsize=(8, 4))

        # Plot terrain profile
        plt.plot(distances, elevations, label="DTM Profile", color="teal", zorder=2)

        # Horizontal line at curb
        plt.axhline(y=curb_level, color="orange", linestyle="--", label="Curb Level", zorder=1)

        # Fill submerged area under curb level
        clipped_z = np.minimum(elevations, curb_level)
        plt.fill_between(distances, clipped_z, curb_level, where=clipped_z < curb_level,
                         color="skyblue", alpha=0.5, label="Flow Area", zorder=0)

        # Add title with hydraulics
        plt.title(f"Cross Section {oid} | A = {area:.2f} m² | Qₙ = {q_full:.2f} m³/s")

        plt.xlabel("Distance along cross section [m]")
        plt.ylabel("Elevation [m]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        # Save
        output_file = os.path.join(output_folder, f"profile_{oid}.png")
        plt.savefig(output_file, dpi=150)
        # plt.close()

        # Print or store
        print(f"OID {oid} | Area = {area:.3f} m² | Q_full = {q_full:.3f} m³/s")


# Read and process each polyline
with arcpy.da.SearchCursor(polyline_path, ["OID@", "SHAPE@", "curb_level", "slope"]) as cursor:
    for oid, geom, curb_level, slope in cursor:
        coords = [(p.X, p.Y) for part in geom for p in part if p]
        if len(coords) < 2:
            continue  # Skip degenerate

        # Interpolate at 0.2 m resolution
        line = LineString(coords)
        length = line.length
        resolution = 0.2
        num_points = max(2, int(length / resolution))
        distances = np.linspace(0, length, num_points)
        interp_points = [line.interpolate(d) for d in distances]
        xy_coords = [(pt.x, pt.y) for pt in interp_points]

        # Sample elevation using mesh
        elevations = []
        for x, y in xy_coords:
            z = interp((x, y))
            elevations.append(z[0])

        # If curb_level is missing, determine it from levees
        if True:#curb_level is None or (isinstance(curb_level, float) and math.isnan(curb_level)):
            N = 5  # Number of points at each end
            left_max = np.nanmax(elevations[:N])
            right_max = np.nanmax(elevations[-N:])
            curb_level = min(left_max, right_max)
            # print(f"OID {oid} | Calculated curb_level = {curb_level:.2f} m from levees")

        # Clip elevations to below curb_level
        clipped_z = np.minimum(elevations, curb_level)

        # Compute horizontal distances between points (assume constant spacing)
        dx = np.diff(distances)
        dz = np.diff(clipped_z)

        # Cross-sectional area using trapezoidal integration
        area = np.sum(0.5 * dx * (clipped_z[:-1] + clipped_z[1:] - 2 * curb_level) * -1)  # negative area becomes positive
        area = max(area, 0)

        # Wetted perimeter (2D length of submerged profile)
        submerged = np.logical_and(clipped_z[:-1] < curb_level, clipped_z[1:] < curb_level)
        perimeter = np.sum(np.sqrt(dx[submerged] ** 2 + dz[submerged] ** 2))

        # Manning calculation
        n = 1 / 70
        if perimeter > 0:
            R = area / perimeter
            q_full = (1 / n) * area * (R ** (2 / 3)) * (slope ** 0.5)
        else:
            q_full = 0.0

        import matplotlib.patches as patches

        # Start plot
        plt.figure(figsize=(8, 4))

        # Plot terrain profile
        plt.plot(distances, elevations, label="DTM Profile", color="teal", zorder=2)

        # Horizontal line at curb
        plt.axhline(y=curb_level, color="orange", linestyle="--", label="Curb Level", zorder=1)

        # Fill submerged area under curb level
        clipped_z = np.minimum(elevations, curb_level)
        plt.fill_between(distances, clipped_z, curb_level, where=clipped_z < curb_level,
                         color="skyblue", alpha=0.5, label="Flow Area", zorder=0)

        # Add title with hydraulics
        plt.title(f"Cross Section {oid} | A = {area:.2f} m² | Qₙ = {q_full:.2f} m³/s")

        plt.xlabel("Distance along cross section [m]")
        plt.ylabel("Elevation [m]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

        # Save
        output_file = os.path.join(output_folder, f"profile_{oid}_mesh04.png")
        plt.savefig(output_file, dpi=150)
        # plt.close()

        # Print or store
        print(f"OID {oid} | Area = {area:.3f} m² | Q_full = {q_full:.3f} m³/s")


plt.show()