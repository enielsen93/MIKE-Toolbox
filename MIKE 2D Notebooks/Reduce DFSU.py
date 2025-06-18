import mikeio
from datetime import datetime

dfsu_filepath = r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Jyllands Alle\MIKE_REGNVAND\JYL_072\JYL_072_FM_m21fm - Result Files\JYL_072_FM_CDS_5_144_240_2080BaseDefault_2D_overland.dfsu"

print("Reading %s" % dfsu_filepath)
dfs = mikeio.dfsu.Dfsu2DH(dfsu_filepath)

filter_start = datetime.strptime("%s 01:30" % dfs.time[0].strftime("%Y.%m.%d"), "%Y.%m.%d %H:%M")

step = 2#int(60/dfs.timestep) if 60>dfs.timestep else 1

# timesteps = dfs.time[-2:]
timesteps = dfs.time[dfs.time>filter_start]
timesteps = timesteps[0::step]
print("Reading %d timesteps" % len(timesteps))
data = dfs.read(time = timesteps)

dfsu_output_filepath = dfsu_filepath.replace(".dfsu", "_reduced.dfsu")
print("Writing %s" % dfsu_output_filepath)
mikeio.dfsu.write_dfsu(filename = dfsu_output_filepath, data = data)

# print("BOB")/