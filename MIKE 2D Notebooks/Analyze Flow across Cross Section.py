import mikeio
import arcpy
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from matplotlib.animation import PillowWriter
import os

dfsu_file = r"C:\Papirkurv\VBT_218_mGA_CDS_20_138_FMBaseDefault_2D_overland_reduced.dfsu"
# dfsu_file_GA = r"C:\Users\ELNN\OneDrive - Ramboll\Documents\Aarhus Vand\Vesterbro Torv\MIKE_URBAN\05_RESULTS\03_FLOOD\GA\2023-10-30\VBT_STATUS_002_GA_m21fm - Result Files\VBT_STATUS_CDS_20_138_FM_GABaseDefault_2D_overland.dfsu"
print(dfsu_file)
profiles_filepath = r"Profile_2.shp"


# # Customize the map
# m.drawcoastlines()
# m.drawcountries()
#
# # Show the map
# plt.show()

def replace_none_with_closest(data):
    if any(np.isnan(data)):
        mask = np.isnan(data)
        data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
    return data

dfs_NW = mikeio.dfsu.Dfsu2DH(dfsu_file)
# dfs_GA = mikeio.dfsu.Dfsu2DH(dfsu_file_GA)

class Profile:
    def __init__(self):
        self.name = None
        self.x = []
        self.z = []
        self.shape = None
        self.lysning_start = None
        self.lysning_slut = None
        self.fald = None
        self.start_area = None
        self.end_area = None
        self.skiller_i = None

profiles = []
dfs_read = dfs_NW.read(items=[i for i, item in enumerate(dfs_NW.items) if "Q flux" in item.name or "P flux" in item.name])
dfs_q_flux = dfs_read[0]
dfs_p_flux = dfs_read[1]

# dfs_GA_read = dfs_GA.read(items=[i for i, item in enumerate(dfs_GA.items) if "Q flux" in item.name or "P flux" in item.name])
# dfs_GA_q_flux = dfs_GA_read[0]
# dfs_GA_p_flux = dfs_GA_read[1]

# dfs_read = dfs.read(items=[i for i, item in enumerate(dfs.items) if "Total water depth" in item.name])
# dfs_total_water_depth = dfs_read[0]

plt.close("all")
with arcpy.da.SearchCursor(profiles_filepath, ["SHAPE@", "Address"]) as cursor:
    for row in cursor:
        print(row[1])
        fig, ax1 = plt.subplots(figsize=(6.17, 4))
        profile = Profile()
        profile.shape = row[0]

        step = 0.1

        points = []
        intervals = np.arange(0, profile.shape.length + step, step)
        # q = []
        # p = []
        # magn = []
        profile.total_flow = []
        profile.time = []

        for [dfs, q_flux, p_flux, color] in [[dfs_NW, dfs_q_flux, dfs_p_flux, 'b']]:
            q = np.empty((len(intervals), len(q_flux.time)))
            p = np.empty((len(intervals), len(q_flux.time)))
            magn = np.empty((len(intervals), len(q_flux.time)))
            depth = np.empty((len(intervals), len(q_flux.time)))
            surface_elevation = np.empty((len(intervals)))
            water_elevation = np.empty((len(intervals), len(q_flux.time)))
            for x_i, x in enumerate(intervals):
                point = profile.shape.positionAlongLine(x)[0]
                points.append([point.X, point.Y])
                # timestep = 30
                i = dfs.geometry.find_index(coords=[point.X, point.Y])[0]
                if i > 0:
                    q[x_i, :] = q_flux[:,i].values
                    p[x_i, :] = p_flux[:, i].values
                    p[x_i, :] = np.nan_to_num(p[x_i, :], copy = True)
                    q[x_i, :] = np.nan_to_num(q[x_i, :], copy = True)
                    for timestep_i in range(len(p[x_i, :])):
                        vec = np.array([p[x_i, timestep_i], q[x_i, timestep_i]])
                        magn[x_i, timestep_i] = np.sqrt(vec.dot(vec))

                    # magn.append(np.sqrt(vec.dot(vec)))
                    # print(np.sum(q))
                    # p = np.max(p_flux[timestep, i].values)
                    # print(np.sum(p))
            # for mi in range(magn.shape[1]):
            #     magn[:,mi] = replace_none_with_closest(magn[:,mi])

            magn = np.nan_to_num(magn)
            profile.total_flow.append([scipy.integrate.trapezoid(magn[:, mi], intervals) for mi in range(magn.shape[1])])
            profile.time.append(q_flux.time)
            ax1.fill_between(profile.time[-1], profile.total_flow[-1], step='pre', color = color, alpha = 0.4)

        ax1.grid()
        ax1.set_ylim([0, 2])
        ax1.set_ylabel(u"Vandføring [m³]")

        # ax2.step(intervals, surface_elevation, 'k--')
        # line = ax2.step(intervals, surface_elevation + depth[:, 0])
        # axvline = ax1.axvline(x = profile.time[0], color = 'k')
        #
        # ax2.set_ylim([ax2.get_ylim()[0], ax2.get_ylim()[1] + 0.1])

        # print("%s - %1.0f m³ (+ %1.0f m³)" % (row[1].split(",")[0], scipy.integrate.trapezoid(profile.total_flow[0][1:], [t.timestamp() for t in profile.time[0][1:]]),0))
        ax1.set_title("%s - %1.0f m³ (+ %1.0f m³)" % (row[1].split(",")[0], scipy.integrate.trapezoid(profile.total_flow[0][1:], [t.timestamp() for t in profile.time[0][1:]]),0))

        # def draw_frame(frame):
        #     line[0].set_ydata(surface_elevation + depth[:, frame])
        #     axvline.set_xdata(profile.time[frame])

        # ani = FuncAnimation(fig, draw_frame, frames=range(len(profile.time)), repeat=True)
        # ani.save(
        #     row[1].split(",")[0] + "_ani.gif", writer=PillowWriter(fps=5))
        # break
        # plt.show()
        # break

        # ax1.title("%s - %1.0f m³" % (row[1].split(",")[0], scipy.integrate.trapezoid(profile.total_flow[1:], [t.timestamp() for t in dfs_q_flux.time[1:]])))
        # if scipy.integrate.trapezoid(profile.total_flow[1:], [t.timestamp() for t in dfs_q_flux.time[1:]])>1e6:
        #     print("BOB")
        # if scipy.integrate.trapezoid(profile.total_flow[1:], [t.timestamp() for t in dfs_q_flux.time[1:]])>1e5:
        #     print("BOB")
        # profile.time = dfs_q_flux.time
        # plt.step(intervals, magn)
        # print(magn)
        # print("BOB")
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.basename(dfsu_file).split(".")[0])
        if not os.path.exists(folder):
            os.mkdir(folder)
        plt.savefig(os.path.join(folder, row[1].split(",")[0] + ".png"), dpi = 300)
        plt.close()
# plt.legend()
# plt.show()
# print("BOB")

# fc = arcpy.GenerateNearTable_analysis(profile.shape, r"C:\GIS\Adresser.shp", r"in_memory\NearTable")