print("Loading libraries")
import os
print("arcpy")
import arcpy
import numpy as np
import scipy.integrate
print("mikeio")
import mikeio
print("bisect")
import bisect
import matplotlib.pyplot as plt

print("Initializing")
figures = {}
# terrain_file = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\04_DTM\Vesterbro_Plan_16_T20_merged.tif"
terrain_files = [#r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\VBT_Mesh_5m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\04_DTM\Vesterbro_Plan_16__T20_Merged.tif",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_0_2m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_0_4m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_0_6m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_0_8m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_1_0m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_1_5m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_2_0m2_interp.mesh",
    r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Forslag 16\VBT_Mesh_5_0m2_interp.mesh"
    ]

polyline_filepath = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Profiler\Longitudinal Profile.shp"

interpolation = False

profiles = {}

class Profile:
    def __init__(self):
        self.name = None
        self.x = None
        self.z = None
        self.shape = None
        self.lysning_start = None
        self.lysning_slut = None
        self.fald = None
        self.start_area = None
        self.end_area = None
        self.skiller_i = None

def replace_none_with_closest(data):
    if any(np.isnan(data)):
        mask = np.isnan(data)
        data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
    return data

for terrain_file in terrain_files:
    print(terrain_file)
    profiles[os.path.basename(terrain_file)] = {}

    if ".tif" in terrain_file:
        is_rectangular = True

        raster = arcpy.Raster(terrain_file)
        raster_array = np.flip(arcpy.RasterToNumPyArray(terrain_file, nodata_to_value = 0),axis=0)

        lower_left_corner = arcpy.Point(raster.extent.XMin, raster.extent.YMin)
        raster_xs = lower_left_corner.X + np.arange(raster.width)*raster.meanCellWidth + raster.meanCellWidth/2
        raster_ys = lower_left_corner.Y + np.arange(raster.height)*raster.meanCellHeight + raster.meanCellHeight/2

        # if interpolation:
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator((raster_ys, raster_xs), raster_array, bounds_error = None)

        cell_size = raster.meanCellWidth
        arcpy.env.outputCoordinateSystem = arcpy.Describe(raster).spatialReference
        arcpy.env.overwriteOutput = True
    else:
        is_rectangular = False
        dfs = mikeio.dfsu.Mesh(terrain_file)
        node_coordinates = dfs.node_coordinates
        if interpolation:
            from scipy.interpolate import LinearNDInterpolator
            interp = LinearNDInterpolator(node_coordinates[:,:2,], node_coordinates[:,2])
        else:
            element_coordinates = dfs.element_coordinates
            def interp(coords):
                try:
                    return element_coordinates[dfs.geometry.find_index(coords = coords), 2]
                except Exception as e:
                    return np.nan

    polylines = {}
    i = 0
    for row in arcpy.da.SearchCursor(polyline_filepath, ["Navn", "SHAPE@", "Lys_start", "Lys_slut", "Fald"]):
        name = row[0]
        if name == ' ':
            name = "LP_%d" % (i)
            i += 1
        polylines[name] = row[1]

        profiles[os.path.basename(terrain_file)][name] = Profile()

        profile = profiles[os.path.basename(terrain_file)][name]
        profile.name = name
        profile.shape = row[1]
        profile.lysning_start = row[2]+0.1
        profile.lysning_slut = row[3]+0.1
        profile.fald = row[4]

    for name, shape in polylines.items():
        if name not in figures:
            figures[name] = plt.figure(figsize=(6.17*3,5*3))

        plt.figure(figures[name])

        points = []
        step = 0.1
        intervals = np.arange(0, shape.length + step, step)
        for x in intervals:
            point = shape.positionAlongLine(x)[0]
            points.append([point.X, point.Y]) if is_rectangular else points.append([point.X, point.Y])
        if interpolation:
            z = np.fliplr(points)
        else:
            if is_rectangular:
                z = interp(np.fliplr(points))#[raster_array[bisect.bisect_right(raster_ys, point[1])-1, bisect.bisect_right(raster_xs, point[0])-1] for point in points]
                z = replace_none_with_closest(np.array(z))
            else:
                z = []
                for point in points:
                    elev = interp([point[0], point[1]])
                    z.append(elev[0] if elev is not np.nan else np.nan)

                z = replace_none_with_closest(np.array(z))


                # z = np.array([float(value) for value in z])


        # plt.plot(np.sort(z), np.insert(scipy.integrate.cumtrapz(np.sort(z), intervals[np.argsort(z)]), 0, 0), label = i)
        profile = profiles[os.path.basename(terrain_file)][name]

        profile.x = intervals
        profile.z = np.array(z)
        # print((profile.x, profile.z))

        if is_rectangular:
            plt.step(intervals, z, 'k', linewidth = 2, label = os.path.basename(terrain_file))
        else:
            if len(plt.gcf().get_axes()[0].lines)>1:
                plt.gcf().get_axes()[0].lines[1].remove()
            plt.step(intervals, z, '--', label = os.path.basename(terrain_file))

        if profile.lysning_start:
            plt.hlines(profile.lysning_start, min(profile.x), max(profile.x), colors='r', linestyles='--',
                       alpha=0.5)
            plt.hlines(profile.lysning_slut, min(profile.x), max(profile.x), colors='b', linestyles='--', alpha=0.5)

            if not terrain_file == terrain_files[0]:
                profile.skiller_i = profiles[os.path.basename(terrain_files[0])][name].skiller_i
            else:
                distance = profile.x[-1] / 4
                skiller_idx = []
                for a in range(len(profile.x)):
                    if profile.x[a] > distance and profile.x[a] < profile.x[-1] - distance and profile.z[a] > max(
                            profile.lysning_start, profile.lysning_slut):
                        skiller_idx.append(a)
                if len(skiller_idx):
                    profile.skiller_i = skiller_idx[np.argmax(profile.z[skiller_idx])]

            if terrain_file == terrain_files[0]:
                if profile.skiller_i:
                    # plt.vlines(profile.x[profile.skiller_i], plt.ylim()[0], plt.ylim()[1], colors='k',
                    #            linestyles='--')
                    plt.fill_between(profile.x[0:profile.skiller_i], profile.z[0:profile.skiller_i],
                                     profile.lysning_start,
                                     where=profile.z[0:profile.skiller_i] < profile.lysning_start, step="pre",
                                     alpha=0.5)
                    plt.fill_between(profile.x[profile.skiller_i:], profile.z[profile.skiller_i:],
                                     profile.lysning_slut,
                                     where=profile.z[profile.skiller_i:] < profile.lysning_slut, step="pre",
                                     alpha=0.5)
                else:
                    try:
                        plt.fill_between(profile.x, profile.z, min(profile.lysning_start, profile.lysning_slut),
                                         where=profile.z < min(profile.lysning_start, profile.lysning_slut),
                                         step="pre", alpha=0.5)
                    except Exception as e:
                        print(e)

            if profile.skiller_i:
                area_curve = profile.lysning_start - profile.z[0:profile.skiller_i]
                area_curve[area_curve < 0] = 0
                profile.start_area = scipy.integrate.trapezoid(area_curve, profile.x[0:profile.skiller_i])

                area_curve = profile.lysning_slut - profile.z[profile.skiller_i:]
                area_curve[area_curve < 0] = 0
                profile.end_area = scipy.integrate.trapezoid(area_curve, profile.x[profile.skiller_i:])
            else:
                area_curve = min(profile.lysning_start, profile.lysning_slut) - profile.z
                area_curve[area_curve < 0] = 0
                profile.start_area = scipy.integrate.trapezoid(area_curve, profile.x)

        plt.ylim([np.min(profiles[os.path.basename(terrain_files[0])][profile.name].z), np.max(profiles[os.path.basename(terrain_files[0])][profile.name].z)])
        plt.grid(visible=True)
        plt.legend()
        plt.title(name)
        plt.savefig(
            r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Profiler\Analyse\%s-%s.png" % (profile.name, os.path.basename(terrain_file)), dpi = 200)
            # points.append(polyline.positionAlongLine(x)[0])
exit()

# for terrain_file in terrain_files:
#     for name, figure in figures.items():
#         plt.figure(figure)
#         profile = profiles[os.path.basename(terrain_file)][name]
#         # plt.step(profile.x, profile.z)
#
#         if profile.lysning_start:
#             plt.hlines(profile.lysning_start, min(profile.x), max(profile.x), colors='r', linestyles='--', alpha=0.5)
#             plt.hlines(profile.lysning_slut, min(profile.x), max(profile.x), colors='b', linestyles='--', alpha=0.5)
#
#             if not terrain_file == terrain_files[0]:
#                 profile.skiller_i = profiles[os.path.basename(terrain_files[0])][name].skiller_i
#             else:
#                 distance = profile.x[-1]/4
#                 skiller_idx = []
#                 for a in range(len(profile.x)):
#                     if profile.x[a] > distance and profile.x[a] < profile.x[-1] - distance and profile.z[a] > max(
#                             profile.lysning_start, profile.lysning_slut):
#                         skiller_idx.append(a)
#                 if len(skiller_idx):
#                     profile.skiller_i = skiller_idx[np.argmax(profile.z[skiller_idx])]
#
#             if terrain_file == terrain_files[0]:
#                 if profile.skiller_i:
#                     plt.vlines(profile.x[profile.skiller_i], plt.ylim()[0], plt.ylim()[1], colors = 'k', linestyles = '--')
#                     plt.fill_between(profile.x[0:profile.skiller_i], profile.z[0:profile.skiller_i], profile.lysning_start, where = profile.z[0:profile.skiller_i] < profile.lysning_start, step = "pre", alpha = 0.5)
#                     plt.fill_between(profile.x[profile.skiller_i:], profile.z[profile.skiller_i:], profile.lysning_slut,
#                                      where=profile.z[profile.skiller_i:] < profile.lysning_slut, step="pre", alpha = 0.5)
#                 else:
#                     try:
#                         plt.fill_between(profile.x, profile.z, min(profile.lysning_start, profile.lysning_slut),
#                                          where=profile.z < min(profile.lysning_start, profile.lysning_slut), step="pre", alpha = 0.5)
#                     except Exception as e:
#                         print(e)
#
#             if profile.skiller_i:
#                 area_curve = profile.lysning_start - profile.z[0:profile.skiller_i]
#                 area_curve[area_curve<0] = 0
#                 profile.start_area = scipy.integrate.trapezoid(area_curve, profile.x[0:profile.skiller_i])
#
#                 area_curve = profile.lysning_slut - profile.z[profile.skiller_i:]
#                 area_curve[area_curve<0] = 0
#                 profile.end_area = scipy.integrate.trapezoid(area_curve, profile.x[profile.skiller_i:])
#             else:
#                 area_curve = min(profile.lysning_start, profile.lysning_slut) - profile.z
#                 area_curve[area_curve < 0] = 0
#                 profile.start_area = scipy.integrate.trapezoid(area_curve, profile.x)
#
#         plt.grid()
#         plt.legend()
#         plt.title(name)
#         plt.savefig(
#             r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Analyse\%s.png" % name)

import matplotlib.colors as mcolors
cmap = plt.get_cmap('hsv')
norm = mcolors.Normalize(vmin=0, vmax=len(terrain_files))
def mean_squared_error(y_predicted, y_actual):
    return np.square(np.subtract(y_actual, y_predicted)).mean()
import math
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

plt.close("all")
for terrain_file_i in range(len(terrain_files)):
    plt.figure(terrain_files[terrain_file_i])
    differences = []
    for profile_i, profile in enumerate(profiles[os.path.basename(terrain_files[terrain_file_i])].values()):
        if not profile.name == "Janus_La_Cours_Gade_1":
            reference_start_area = profiles[os.path.basename(terrain_files[0])][profile.name].start_area
            if reference_start_area:
                reference_end_area = profiles[os.path.basename(terrain_files[0])][profile.name].end_area
                differences.append(profile.start_area/reference_start_area*1e2)

                plt.scatter(profile_i, profile.start_area / reference_start_area * 1e2, color="r", marker="_")

                if profile.end_area and profile.end_area>0 and reference_end_area and reference_end_area>0:
                    differences.append(profile.end_area/reference_end_area*1e2)
                    plt.scatter(profile_i, profile.end_area/reference_end_area*1e2, color = "b", marker = "_")

    plt.hlines(100, 0, plt.xlim()[-1], colors = 'k', linestyles ='--')
    plt.hlines(np.sqrt(mean_squared_error(np.array(differences), np.array(differences)*0+100))+100, 0, plt.xlim()[-1], colors='k', linestyles='--', alpha = 0.5)
    # plt.fill_between([0, plt.xlim()[-1]], 100 - np.array([1, 1]) * np.sqrt(mean_squared_error(np.array(differences), np.array(differences)*0+100)), np.array([1, 1]) * np.sqrt(mean_squared_error(np.array(differences), np.array(differences)*0+100))+100, facecolor = "k", alpha = 0.2)
    plt.title(os.path.basename(terrain_files[terrain_file_i]))
    plt.ylim([0, 250])
    plt.grid(which='both')
    plt.text(0.05, 0.95, "RMSE %d%%" % np.sqrt(mean_squared_error(np.array(differences), np.array(differences)*0+100)), transform=plt.gcf().get_axes()[0].transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    # print(differences)
    plt.savefig(
        r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\07_MESH\Profiler\Analyse\%s_area.png" % (os.path.basename(terrain_files[terrain_file_i])))

# for terrain_file in terrain_files:
#     for profile in profiles[os.path.basename(terrain_file)].values():
#         if profile.name == "Langelandsgade_Nord_1":
#             print((profile.name, profile.start_area, reference_start_area, profile.skiller_i))
    # print(profiles[terrain_file]["Langelandsgade_2"].start_area)


        # plt.show()

    # plt.close(figure)
#
# workspace = r'C:½\Papirkurv'
# output_fc = 'Slope3.shp'
# arcpy.env.overwriteOutput = True
# slope_fc = arcpy.CreateFeatureclass_management(workspace, output_fc, 'POINT')[0]
# arcpy.AddField_management(slope_fc, 'Slope', 'FLOAT')
#
# import pandas as pd
# for profile_name in profiles["Vesterbro_Plan_16_T20_merged.tif"].keys():
#     if "LP" in profile_name:
#         profile = profiles["Vesterbro_Plan_16_T20_merged.tif"][profile_name]
#         plt.figure(figsize=(6.17*3,5*3))
#         plt.step(profile.x, profile.z)
#         step = 500
#         # xs = np.arange(np.min(profile.x), np.max(profile.x)+step, step)
#         # series = pd.Series(profile.x, profile.z)
#         plt.xlim(np.min(profile.x), np.max(profile.x))
#         plt.ylim(np.min(profile.z), np.max(profile.z))
#         ax2 = plt.twinx()
#         slope = (np.array(profile.z[step:])-np.array(profile.z[:-step]))/(np.array(profile.x[step:])-np.array(profile.x[:-step]))*1e3
#         xs = np.array(profile.x[step:])
#         ax2.step(profile.x[:-step], (np.array(profile.z[step:])-np.array(profile.z[:-step]))/(np.array(profile.x[step:])-np.array(profile.x[:-step]))*1e3)
#         ax2.grid()
#
#         with arcpy.da.InsertCursor(slope_fc, ["SHAPE@","Slope"]) as cursor:
#             for x, slope in zip(xs[0::50],slope[0::50]):
#                 row = [profile.shape.positionAlongLine(x), -slope]
#                 cursor.insertRow(row)


# output_raster = arcpy.NumPyArrayToRaster(final_raster, lower_left_corner, cell_size, cell_size, 0)

# output_raster.save(r"C:\Users\ELNN\Downloads\Total.tif")



print("bob")