# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:58:40 2022

@author: janousu
"""

### RASTER INPUTS FOR SPAFHY DEMO
import os
#os.chdir(r'F:/SpaFHy_v1_Pallas_2D/PALLAS_RAW_DATA')
from spafhy_gis_inputs import create_catchment, write_AsciiGrid


# -- folder to data!
gispath = r'F:\SpaFHy_2D_2021\PALLAS_RAW_DATA\Lompolonjanka\16b'
# ID=1 is Lompolojanganoja
pgen={'catchment_id': 1}

gis = create_catchment(pgen, gispath, plotgrids=False, plotdistr=False)


for key in dict.keys(gis):
    fname = 'C:/SpaFHy_v1_Pallas_2D/testcase_input/new_parameters/' + key + '.dat'
    write_AsciiGrid(fname, gis[key], gis['info'])
    print(key + '= done')