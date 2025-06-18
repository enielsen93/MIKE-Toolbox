# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:25:43 2021

@author: mu
"""

import mikeio
import numpy as np
filename = r"K:\Hydrauliske modeller\Modeller\A265_1 Kapacitetsanalyse Sydbyen\Mesh\A265_1_Mesh_GA_v2_1_cropped_interp_E12.mesh"

dfsu = mikeio.dfsu.Dfsu(filename)

# border_nodes = np.where(dfsu.codes == dfsu.boundary_codes)[0]
element_table = dfsu.element_table
unconnected_elements = []#np.ones((len(element_table)), dtype=bool)
for element1_i, element1 in enumerate(element_table):
    
    if len(np.intersect1d(element1, np.delete(element_table,element1_i)))==0:
        print(element1_i)
        unconnected_elements.append(element1_i)
    # for element2 in element_table:
        # if len(np.intersect1d(element1,element2))>1 and element1 is not element2:
            # unconnected_elements[element1_i] = True
            # continue

# for node in set(nodes_in_element_table):
#     if nodes_in_element_table
    

# delete_elements = []
# for element_i, element in enumerate(dfsu.element_table):
#     # print([True for node in element if node not in border_nodes])
#     if not True in [True for node in element if node not in border_nodes]:
#         delete_elements.append(element_i)
#         print(element_i)
#         print([True for node in element if node not in border_nodes])
        
els = dfsu.element_coordinates[unconnected_elements,0:2]