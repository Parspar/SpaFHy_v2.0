# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:29:40 2022

@author: janousu
"""

from iotools import read_results
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from datetime import date
import os
from netCDF4 import Dataset #, date2num
import pandas as pd
#import pickle
#import seaborn as sns
from iotools import read_AsciiGrid, write_AsciiGrid
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from sklearn.metrics import mean_absolute_error as mae
from rasterio.plot import show
import rasterio
from rasterio.transform import from_origin
from raster_utils import read_pkrasteri_for_extent


#%%



# reading the stand results
outputfile_1d = r'D:\SpaFHy_2D_2021\testcase_input_1d_new.nc'
results_1d = read_results(outputfile_1d)

# reading the stand results
outputfile_2d = r'D:\SpaFHy_2D_2021\testcase_input_2d_new.nc'
results_2d = read_results(outputfile_2d)

# reading the catch results
outputfile_top = r'D:\SpaFHy_2D_2021\testcase_input_top_new_fixed.nc'
results_top = read_results(outputfile_top)

sar_temp = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment_ma3.nc'
sar_temp = Dataset(sar_temp, 'r')

sar_spat = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment.nc'
sar2_spat = Dataset(sar_spat, 'r')


