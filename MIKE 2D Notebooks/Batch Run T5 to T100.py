# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 18:29:13 2021

@author: mu
"""
import subprocess
import os
import mikeio
import numpy as np
import re
import copy
import arcpy

def run_mex(mouse_sim_launch, mex_file, simulation_type = "HD", timeout = None):
    run_cmd = r'"%s" "%s" "%s" "Run" "NoPrompt" "Close"' % (mouse_sim_launch, mex_file, simulation_type)
    subprocess.run(run_cmd, timeout = timeout)
    
mouse_sim_launch_paths = [r"C:\Program Files (x86)\DHI\2020\bin\x64\MOUSESimLaunch.exe"]
for path in mouse_sim_launch_paths:
    for year in reversed(range(2010,2030)):
        if os.path.exists(path.replace("2020",str(year))):
            mouse_sim_launch = path.replace("2020",str(year))
            break

CDS_folder = r"K:\Hydrauliske modeller\Regn\Valg af Regn 2018\CDS"
return_periods = [5, 10, 20, 50, 100]
safety_factors = {2050: [0.96, 0.97, 0.99, 1.00*1.25, 1.02*1.25], 2100: [1.08, 1.10, 1.12, 1.16*1.25, 1.19*1.25]}
CDS_files = [os.path.join(CDS_folder, file) for file in os.listdir(CDS_folder) if 
             any([True for return_period in return_periods if str(return_period) in file])]
CDS_files_sorted = []
for return_period in return_periods:
    i = [i for i,file in enumerate(CDS_files) if int(re.findall(r"CDS(\d+)aar", file)[0]) == return_period][0]
    CDS_files_sorted.append(CDS_files[i])

# mex_file = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test\Nebelvej.mex"
# mex_file = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test\Aarhusbakken.mex"
mex_file = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test\Funder.mex"
# mex_file = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test\Sejling Aadal.mex"
# mex_file = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test\Astrid Lindgrensvej.mex"
# mdb_file = r"K:\Hydrauliske modeller\Modeller\A0020-49 Astrid Lindgrensvej\Model\A0020-49 Astrid Lindgrensvej v1.5.mdb"

with open(mex_file,"r") as f:
    mex_file_text = f.readlines()
    
mex_file_text_CRF_lineno = [lineno for lineno,line in enumerate(mex_file_text) if "CRF_file" in line][0]
mex_file_text_model_A_lineno = [lineno for lineno,line in enumerate(mex_file_text) if "[Model_A]" in line][0]
mex_file_text_model_A_lineno_end = [lineno for lineno,line in enumerate(mex_file_text) if "EndSect  // Model_A" in line][0]
mex_file_text_TA_model_lineno = [lineno for lineno in np.arange(mex_file_text_model_A_lineno, mex_file_text_model_A_lineno_end) if "Line = " in mex_file_text[lineno] 
                                 and float(re.findall(r"', ([\d\.]+)", mex_file_text[lineno])[0])>0]
mex_file_text_dfs0_lineno = [lineno for lineno,line in enumerate(mex_file_text) if ".dfs0" in line][0]


mex_file_text_nodes_lineno = [lineno for lineno,line in enumerate(mex_file_text) if "Node =" in line]
mex_file_text[mex_file_text_nodes_lineno[0]-1]
replace_top_type_re = re.compile(r"(Node =[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,[^,]+,) [13](.+)")
for lineno in mex_file_text_nodes_lineno:
    if re.findall("Node = '[^']+', (\d)", mex_file_text[lineno])[0] == "1":
        # mex_file_text[lineno] = replace_top_type_re.sub(r"\1 3\2",mex_file_text[lineno])
        pass
    elif re.findall("Node = '[^']+', (\d)", mex_file_text[lineno])[0] == "3": 
        mex_file_text[lineno] = replace_top_type_re.sub(r"\1 1\2",mex_file_text[lineno])
        
re_get_row_values = re.compile(r">([^<]+)<\/TD>")
spilling_table = {}
total_spilled = {2050: [], 2100: []}
for year in safety_factors:
    for return_period_i in range(len(return_periods)):
        mex_name = "CDSY%03dT%03d.mex" % (year, return_periods[return_period_i])
        # print(mex_name)
        mex_file_new = os.path.join(os.path.dirname(mex_file), mex_name)
        mex_file_new_text = copy.deepcopy(mex_file_text)
        
        mex_file_new_text[mex_file_text_CRF_lineno] = re.sub("'.+'", 
                                                             "'%s.CRF'" % mex_file_new.replace(".mex",""), 
                                                             mex_file_new_text[mex_file_text_CRF_lineno])
        
        dfs0 = mikeio.Dfs0(CDS_files_sorted[return_period_i])
        mex_file_new_text[mex_file_text_dfs0_lineno] = re.sub(r"(.+)'[^']+.dfs0',[^,]+,[^,]+,(.+)", 
                                                              r"\1'%s', '%s', '%s',\2" % (CDS_files_sorted[return_period_i],
                                                              dfs0.items[0].type.display_name,
                                                              dfs0.items[0].name), 
                                                              mex_file_new_text[mex_file_text_dfs0_lineno])
        for lineno in mex_file_text_TA_model_lineno:
            mex_file_new_text[lineno] = re.sub(r"(.+',) [\d\.]+(.+)", r"\1 %1.2f\2" % (safety_factors[year][return_period_i]), 
                   mex_file_new_text[lineno])
        
        with open(mex_file_new,'w+') as f:
            f.writelines(mex_file_new_text)
        try:
            run_mex(mouse_sim_launch, mex_file_new,"RO")  
        except subprocess.TimeoutExpired:
            pass

        try:
            run_mex(mouse_sim_launch, mex_file_new,"HD")  
        except subprocess.TimeoutExpired:
            pass
        
        summary_filepath = os.path.join(os.path.dirname(mex_file), 
                                        "Summary_HD_%s" % (mex_name.replace(".mex", ".HTM")))
        html = np.array(open(summary_filepath).read().splitlines())
        spilling_table_lineno = [lineno for lineno, line in enumerate(html) if "Nodes - Volume spilled" in line][1]
        spilling_table_lineno_end = [lineno for lineno, line in enumerate(html) if r"/TABLE" in line and lineno>spilling_table_lineno][0]
        for line in html[range(spilling_table_lineno+4, spilling_table_lineno_end)]:
            match = re_get_row_values.findall(line)
            spilling_table[match[0]] = float(match[1])
        print(np.sum(list(spilling_table.values())))
        total_spilled[year].append(np.sum(list(spilling_table.values())))
        
        # output_folder = r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test"
        # scalgo_fc = arcpy.management.CreateFeatureclass(output_folder, 
        #                                     os.path.basename(mex_file_new.replace(".mex","Point.shp")), "POINT")[0]
        # # scalgo_fc = arcpy.management.CreateFeatureclass(r"K:\Hydrauliske modeller\Papirkurv\Scalgo-test", 
        #                                     # os.path.basename(mex_file_new.replace(".mex","Polygon.shp")), "POLYGON")[0]
        # arcpy.management.AddField(scalgo_fc, "Spill", "FLOAT", 15, 5)
        # nodes_dict = {}
        # with arcpy.da.InsertCursor(scalgo_fc, ["SHAPE@", "Spill"]) as insert_cursor:
        #     with arcpy.da.SearchCursor("K:\Hydrauliske modeller\Modeller\A0020-49 Astrid Lindgrensvej\Model\msm_Node.shp",
        #                                ["MUID","SHAPE@"]) as cursor:
        #         for row in cursor:
        #             if row[0] in spilling_table:
        #                 insert_cursor.insertRow([row[1], spilling_table[row[0]]])
                        
        # arcpy.Buffer_analysis(scalgo_fc, os.path.join(output_folder, os.path.basename(mex_file_new.replace(".mex","Polygon.shp"))), 
        #                                                   "Spill")
        