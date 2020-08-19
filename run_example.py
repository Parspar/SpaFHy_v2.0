# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:25:15 2019

@author: khaahti
"""

from model_driver import driver
from iotools import read_results
import matplotlib.pyplot as plt

# example of calling driver, reading results and plotting gwl

outputfile = driver(create_ncf=True, folder='testcase_input')

results = read_results(outputfile)

plt.figure()
ax = plt.subplot(4,1,1)
results['canopy_snow_water_equivalent'][:,0,:].plot.line(x='date')

plt.subplot(4,1,2, sharex=ax)
results['forcing_air_temperature'][:].plot.line(x='date')

plt.subplot(4,1,3, sharex=ax)
results['soil_ground_water_level'][:,5,:].plot.line(x='date')

plt.subplot(4,1,4, sharex=ax)
results['canopy_transpiration'][:,0,:].plot.line(x='date')

plt.figure()
results['canopy_leaf_area_index'][:,0,:].plot.line(x='date')

plt.figure()
results['soil_ground_water_level'][:,1:-1,1:-1].mean(dim='date').plot()

results.close()
