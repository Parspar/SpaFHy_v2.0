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
results['soil_rootzone_moisture'][:,5,:].plot.line(x='date')

plt.subplot(4,1,3, sharex=ax)
results['soil_ground_water_level'][:,5,:].plot.line(x='date')

plt.subplot(4,1,4, sharex=ax)
results['canopy_leaf_area_index'][:,5,:].plot.line(x='date')

plt.figure()
results['soil_drainage'].mean(dim='date').plot()

plt.figure()
results['soil_drainage'][:,1:-1,1:-1].mean(dim=['i','j']).plot()

plt.figure()
results['parameters_elevation'].plot()

results.close()