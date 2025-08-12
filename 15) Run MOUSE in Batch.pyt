# Tool for reading DFS0 or KM2 files and creating LTS files from it
# Created by Emil Nielsen
# Contact: 
# E-mail: enielsen93@hotmail.com

import os
import sys
import numpy as np
thisFolder = os.path.dirname(__file__)
scriptFolder = os.path.join(thisFolder, r"scripts")
sys.path.append(scriptFolder)
import re
import copy
from shutil import copyfile
import subprocess
import time
import multiprocessing as mp

def run_mex(mouse_sim_launch, mex_file, parallel = False):
    run_cmd = r'"%s" "%s" "HD" "Run" "Close" "NoPrompt" "-wait"' % (mouse_sim_launch, mex_file)
    if parallel:
        subprocess.Popen(run_cmd)
    else: 
        subprocess.call(run_cmd)

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Batch Simulations"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [BatchSimulations]
    
    
class BatchSimulations(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Run Mex Files"
        self.description = "Run Mex Files. \n\nCreated by: Emil Nielsen \nContact: enielsen93@hotmail.com"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions
        LTSCount = arcpy.Parameter(
            displayName="Run this many simulations parallel?",
            name="LTSCount",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        LTSCount.value = 4
        
        RunoffFile = arcpy.Parameter(
            displayName="Runoff file (for making dupilcates)",
            name="RunoffFile",
            datatype="File",
            parameterType="Optional",
            direction="Input",)
        RunoffFile.filter.list=["crf", "res1d"]
        
        mex_file = arcpy.Parameter(
            displayName="Network Setup File",
            name="mex_file",
            datatype="File",
            parameterType="Required",
            multiValue="True",
            direction="Input")
        mex_file.filter.list=["mex", "m1dx"]
        
        parameter = arcpy.Parameter(
            displayName="Values for parameter (separate by commas)",
            name="parameter",
            datatype="String",
            parameterType="Required",
            direction="Input")

        mouse_sim_launch = arcpy.Parameter(
            displayName="MOUSE Sim Launch Executable path",
            name="mouse_sim_launch",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        mouse_sim_launch.filter.list=["exe"]
        
        mouse_sim_launch_paths = [r"C:\Program Files (x86)\DHI\2020\bin\x64\MOUSESimLaunch.exe"]
        for path in mouse_sim_launch_paths:
            for year in reversed(range(2010,2030)):
                if os.path.exists(path.replace("2020",str(year))):
                    mouse_sim_launch.value = path.replace("2020",str(year))
                    break

        mike1d_launch = arcpy.Parameter(
            displayName="MIKE1Ds Launch Executable path",
            name="mike1d_launch",
            datatype="File",
            parameterType="Optional",
            direction="Input")
        mike1d_launch.filter.list = ["exe"]

        mike1d_paths = [r"C:\Program Files (x86)\DHI\2020\bin\x64\DHI.Mike1D.Application.exe"]
        for path in mike1d_paths:
            for year in reversed(range(2010, 2030)):
                if os.path.exists(path.replace("2020", str(year))):
                    mike1d_launch.value = path.replace("2020", str(year))
                    break
        
        
        params = [LTSCount, RunoffFile, mex_file, parameter, mouse_sim_launch, mike1d_launch]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        if parameters[2].Values:
            for mex_file in parameters[2].ValueAsText.split(";"):
                if ".mex" in mex_file:
                    with open(mex_file,"r") as f:
                        mex_file_text = f.read()
                    if "-9999" not in mex_file_text:
                        parameters[2].setErrorMessage("No paramater '-9999' in mex-file")
        return

    def execute(self, parameters, messages):
        # Extract parameters from input list with clear variable names
        LTSCount = int(parameters[0].ValueAsText)                    # Max concurrent processes allowed
        RunoffFile = parameters[1].ValueAsText                       # Path to runoff file
        mex_files_original = [mex_file.replace("'", "") for mex_file in parameters[2].ValueAsText.split(";")]  # List of mex files, cleaned from quotes
        parameter_str = parameters[3].ValueAsText                     # Parameter string containing numeric values
        mouse_sim_launch = parameters[4].ValueAsText                  # Path to Mouse simulation executable
        mike1d_launch = parameters[5].ValueAsText                     # Path to Mike1D executable

        # Extract numeric parameters from the parameter string (e.g., ["10", "30", "50"])
        parameter_values = re.findall(r"[\d\.]+", parameter_str)

        # Prepare a list to store new mex file names created during the loop
        all_new_mex_files = []
        processes = []

        # Loop over each original mex file
        for mex_file in mex_files_original:
            # Read the mex file content into a list of lines
            with open(mex_file, "r") as f:
                mex_file_text = f.readlines()

            # Find line numbers containing the placeholder parameter "par1"
            par1_line_numbers = [lineno for lineno, line in enumerate(mex_file_text) if "-9999" in line]

            # Determine the line number to update the runoff or CRF file reference based on file content and RunoffFile presence
            if RunoffFile and "mex" in mex_file:
                # Use the line with CRF_file keyword if runoff file exists and mex file extension present
                crf_line_number = next(lineno for lineno, line in enumerate(mex_file_text) if "CRF_file" in line)
            else:
                # Otherwise, look for line containing "RR.res1d" as fallback
                crf_line_number = next(lineno for lineno, line in enumerate(mex_file_text) if "RR.res1d" in line)

            # Process each parameter value separately, creating new mex files customized for each parameter
            for idx, param_value in enumerate(parameter_values):
                if RunoffFile:
                    print("Copying %s" % RunoffFile)
                    if ".mex" in mex_file:
                        RunoffFileNew = "%s_Split%d.CRF" % (RunoffFile[:-4], idx + 1)
                    else:
                        RunoffFileNew = "%s_Split%d.res1d" % (RunoffFile[:-4], idx + 1)
                    copyfile(RunoffFile, RunoffFileNew)

                if ".mex" in mex_file:
                    mex_file_new = mex_file.replace(".mex", "_par1-%s.mex" % param_value)
                else:
                    mex_file_new = mex_file.replace(".m1dx", "_par1-%s.m1dx" % param_value)

                all_new_mex_files.append(mex_file_new)

                mex_file_new_text = copy.deepcopy(mex_file_text) # make a copy

                for lineno in par1_line_numbers:
                    mex_file_new_text[lineno] = re.sub("-9999", param_value, mex_file_new_text[lineno])

                if ".m1dx" in mex_file:
                    inside_rs = False
                    for i, line in enumerate(mex_file_new_text):
                        if "<ResultSpecification>" in line:
                            inside_rs = True
                        elif "</ResultSpecification>" in line:
                            inside_rs = False

                        if inside_rs and "<Path>" in line:
                            arcpy.AddMessage(line)
                            def repl(match):
                                content = match.group(1)
                                new_content = content.replace("Base.res1d", "_%sBase.res1d" % str(param_value).replace(".","_"))
                                return "<Path>%s</Path>" % new_content

                            mex_file_new_text[i] = re.sub(r"<Path>(.*?)</Path>", repl, line)
                            arcpy.AddMessage(line)
                            arcpy.AddMessage((r"<Path>(.*?)</Path>", repl, line))

                if RunoffFile:
                    mex_file_new_text[crf_line_number] = re.sub(r"'[^']*'", "'%s'" % RunoffFileNew,
                                                                mex_file_new_text[crf_line_number])
                else:
                    mex_file_new_text[crf_line_number] = re.sub(r"'[^']*'",
                                                                "'%s'" % mex_file_new.replace(".mex", ".CRF"),
                                                                mex_file_new_text[crf_line_number])

                with open(mex_file_new, 'w') as f:
                    f.writelines(mex_file_new_text)

                if mouse_sim_launch and ".mex" in mex_file:
                    while processes and sum(1 for proc in processes if proc.poll() is None) >= LTSCount:
                        time.sleep(5)

                    if not RunoffFile:
                        run_cmd = 'cmd.exe /k echo Running: "%s" "%s" "RO" "Run" "Close" "NoPrompt" "-wait" & "%s" "%s" "RO" "Run" "Close" "NoPrompt" "-wait"' % (
                            mouse_sim_launch, mex_file_new, mouse_sim_launch, mex_file_new
                        )
                        subprocess.check_output(run_cmd)

                    run_cmd = 'cmd.exe /k echo Running: "%s" "%s" "HD" "Run" "Close" "NoPrompt" "-wait" & "%s" "%s" "HD" "Run" "Close" "NoPrompt" "-wait"' % (
                        mouse_sim_launch, mex_file_new, mouse_sim_launch, mex_file_new
                    )
                    processes.append(subprocess.Popen(run_cmd))
                    time.sleep(1)

                elif mike1d_launch and "m1dx" in mex_file:
                    while processes and sum(1 for proc in processes if proc.poll() is None) >= LTSCount:
                        time.sleep(5)

                    run_cmd = '"%s" "%s"' % (mike1d_launch, mex_file_new)
                    processes.append(subprocess.Popen(run_cmd))
                    time.sleep(1)


