# -*- coding: utf-8 -*-
"""
Created on Mon May  3 13:08:47 2021

@author: janousu
"""

# figures

from iotools import read_results
import matplotlib.pyplot as plt

outputfile = 'results/testcase_input.nc'

results = read_results(outputfile)

# plt.figure()
# ax = plt.subplot(4,1,1)
# results['canopy_snow_water_equivalent'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,2, sharex=ax)
# results['soil_rootzone_moisture'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,3, sharex=ax)
# results['soil_ground_water_level'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,4, sharex=ax)
# results['canopy_leaf_area_index'][:,4,:].plot.line(x='date')

plt.figure()
results['soil_ground_water_level'][:,150,100].plot.line(x='date')

plt.figure()
results['soil_ground_water_level'][-1,:,:].plot()

plt.figure()
results['soil_drainage'][:,:,:].mean(dim='date').plot()

plt.figure()
results['parameters_elevation'][:,:].plot()

plt.figure()
results['soil_water_closure'][:,1:-1,1:-1].mean(dim='date').plot()

plt.figure()
results['soil_drainage'][:,1:-1,1:-1].mean(dim=['i','j']).plot()

plt.figure()
results['soil_water_closure'][:,1:-1,1:-1].mean(dim=['i','j']).plot()

plt.figure()
results['soil_rootzone_moisture'][:,150,100].plot.line(x='date')
results.close()