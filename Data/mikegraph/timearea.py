"""
Timearea module for mikegraph
"""
# import mikegraph
from io import StringIO
import pandas as pd
import numpy as np
import networkx as nx
import os
import importlib

def import_or_install(pkg_names):
    import importlib
    import site
    import tkinter as tk
    from tkinter import messagebox
    import subprocess
    imported = {}
    for pkg in pkg_names:
        try:
            imported[pkg] = importlib.import_module(pkg)
            continue
        except ImportError:
            pass

        # Check user site-packages
        user_site = site.getusersitepackages()
        candidate = os.path.join(user_site, pkg)
        if os.path.isdir(candidate):
            sys.path.insert(0, user_site)
            try:
                imported[pkg] = importlib.import_module(pkg)
                continue
            finally:
                sys.path.pop(0)

        # Not found: prompt user with tkinter
        import tkinter as tk
        root = tk.Tk()
        root.withdraw()  # hide main window
        msg = f"The library '{pkg}' is not installed.\nInstall now using ArcGIS Pro Python?"
        if messagebox.askokcancel("Missing Library", msg):
            propy_path = r"C:\Progra~1\ArcGIS\Pro\bin\Python\scripts\propy.bat"
            cmd = [propy_path, "-m", "pip", "install"] + pkg_names
            subprocess.check_call(cmd)
            # Try import again
            imported[pkg] = importlib.import_module(pkg)
        else:
            raise ImportError(f"{pkg} not installed and user declined installation.")

    return imported

class TimeAreaAnalyzer:
    """
       Load and process rainfall time series data from multiple file formats.

       This class handles rainfall data import from both MIKE DFS0 files and text-based
       formats (CSV/TSV). Automatically resamples data to 1-minute intervals and prepares
       it for hydrological analysis and modeling.

       Parameters:
           rain_filepath (str): Path to rainfall data file. Supports:
               - DFS0 files (.dfs0) - requires mikeio
               - Text files (CSV/TSV) with comma or tab delimiters
               - Files with European decimal notation (comma as decimal separator)

       Notes:
           - DFS0 files are read using mikeio and converted to pandas format
           - Text files expect intensity data starting from row 4 (skiprows=3)
           - European decimal notation (commas) is automatically converted to dots
           - All data is resampled to 60-second intervals using backward fill
           - 60 zero values are appended to the rainfall event array

       Example:
           >>> # Load DFS0 rainfall file
           >>> rain = TimeAreaAnalyzer("rainfall.dfs0")
       """
    def __init__(self, rain_filepath):
        if os.path.splitext(rain_filepath)[1].lower() == ".dfs0":
            try:
                libs = import_or_install(["mikeio"])
                mikeio = libs["mikeio"]
            except ImportError:
                raise ImportError(
                    "The 'mikeio' package is required but not installed.\n"
                    "You can install it by running:\n"
                    "   python -m pip install mikeio\n\n"
                    "where python is the path to your ArcGIS python.exe"
                )
            dfs0 = mikeio.open(rain_filepath)
            self.series = dfs0.read().to_pandas()
            self.series = self.series.resample("60s").bfill()
            self.rain_event = np.concatenate((self.series.values[:], np.zeros(60)))

        else:
            with open(rain_filepath, 'r') as f:
                txt = f.read()

            delimiter = r"  " if "  " in txt else r"\t"
            if "," in txt:
                txt = txt.replace(r",", r".")
                rain_filepath = StringIO(unicode(txt))

            self.series = pd.read_csv(rain_filepath, delimiter=delimiter, skiprows=3, names=["Intensity"], engine='python')
            self.series.index = pd.to_datetime(self.series.index)
            self.series = self.series.resample("60s").bfill()
            self.rain_event = np.concatenate((self.series.values[:, 0], np.zeros(60)))

        self.additional_discharge = {}
        self.scaling_factor = 1

    def rationelCurve(self, target, graph):
        sources = graph.find_upstream_nodes(target)[0]

        total_runoff = np.zeros(len(self.rain_event))
        for source in sources:
            if source in self.additional_discharge:
                total_runoff += self.additional_discharge[source]*1e3

            for catchment in graph.find_connected_catchments(source):
                total_runoff += self.rain_event/1e6*catchment.reduced_area*1e3*self.scaling_factor
        return total_runoff

    def timeareaCurve(self, target, graph):
        sources = graph.find_upstream_nodes(target)[0]

        total_runoff = np.zeros(len(self.rain_event))
        for source in sources:
            if source in self.additional_discharge:
                total_runoff += self.additional_discharge[source]*1e3

            travel_time = graph.travel_time(source, target)

            for catchment in graph.find_connected_catchments(source):
                runoff = np.zeros(len(self.rain_event))
                for time_i, rain_intensity in enumerate(self.rain_event):
                    time_i_adjusted = time_i - travel_time/60
                    rain = self.rain_event[
                           max(int(time_i_adjusted - catchment.concentration_time), 0):max(0, int(time_i_adjusted))]
                    runoff[time_i] = np.sum(rain) / catchment.concentration_time if rain.any() else 0
                total_runoff += runoff/1e6*catchment.reduced_area*1e3*self.scaling_factor

        return total_runoff

if __name__ == "__main__":
    # graph = mikegraph.MikeNetwork(r"C:\Users\elnn\OneDrive - Ramboll\Documents\Aarhus Vand\Soenderhoej\MIKE\MIKE_URBAN\_ORIGINAL\Viby_detailed_200101_40\Viby_detailed_200101_40.sqlite")
    # graph.map_network()

    rainseries = TimeAreaAnalyzer(r"C:\Papirkurv\Webinar presentation\02_RAIN\CDS_5 min_samlet.DFS0")
    # rainseries.additional_discharge = {"SEMI25":0.25}
    #
    # discharge_ta = rainseries.timeareaCurve(u'D43440R', graph)
    # # import timeit
    # # timeit.timeit(lambda: rainseries.timeareaCurve(u'Haslevangsvej_R719', graph), number = 5)
    # discharge_rat = rainseries.rationelCurve(u'D43440R', graph)
    # #
    # import matplotlib.pyplot as plt
    # #
    # plt.step(range(len(discharge_rat)), discharge_rat)
    # plt.step(range(len(discharge_ta)), discharge_ta)
    # plt.show()
    # print("PAUSE")
# [u'O23119R', u'O23117R-1', u'O23116R', u'O23117R', u'O23114R', u'O23134R', u'O23118R', u'O23115R']