# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:25:15 2019

@author: khaahti
"""

from model_driver_stand import driver
from iotools import read_results
import matplotlib.pyplot as plt
import numpy as np

# example of calling driver, reading results and plotting gwl

outputfile = driver(create_ncf=True, folder='testcase_input')

results = read_results(outputfile)

results['bucket_water_storage'] = results['bucket_water_storage'] - results['bucket_water_storage'][0,:,:]

plt.figure(figsize=(25,15))
ax=plt.subplot(2,4,1)
results['parameters_elevation'][:,:].plot()
plt.subplot(2,4,3, sharex=ax, sharey=ax)
results['bucket_moisture_root'][-1,:,:].plot()
plt.subplot(2,4,5, sharex=ax, sharey=ax)
results['bucket_surface_runoff'][-1,:,:].plot()
plt.subplot(2,4,8, sharex=ax, sharey=ax)
results['bucket_water_closure'][-1,:,:].plot()

# plt.figure(figsize=(20,15))
# ax=plt.subplot(2,3,1)
# results['soil_ground_water_level'][30,:,:].plot(vmin=-3,vmax=3, cmap='RdBu')
# plt.subplot(2,3,2, sharex=ax, sharey=ax)
# results['soil_rootzone_moisture'][30,:,:].plot(vmin=0,vmax=1, cmap='RdBu')
# plt.subplot(2,3,3, sharex=ax, sharey=ax)
# results['soil_moisture_top'][30,:,:].plot(vmin=0,vmax=1, cmap='RdBu')
# plt.subplot(2,3,4, sharex=ax, sharey=ax)
# results['soil_ground_water_level'][-53,:,:].plot(vmin=-3,vmax=3, cmap='RdBu')
# plt.subplot(2,3,5, sharex=ax, sharey=ax)
# results['soil_rootzone_moisture'][-53,:,:].plot(vmin=0,vmax=1, cmap='RdBu')
# plt.subplot(2,3,6, sharex=ax, sharey=ax)
# results['soil_moisture_top'][-53,:,:].plot(vmin=0,vmax=1, cmap='RdBu')

plt.figure(figsize=(20,8))
ax=plt.subplot(2,1,1)
plt.plot(results['date'],results['bucket_drainage'].mean(['i','j']),label='bucket_drainage')
plt.legend()

plt.subplot(2,1,2,sharex=ax)
plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j']), label='bucket_water_storage')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])), label='+canopy_evaporation')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])), label='+canopy_transpiration')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])), label='+bucket_evaporation')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])), label='+bucket_surface_runoff')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])
                   + results['bucket_drainage'].mean(['i', 'j'])), label='+bucket_drainage')

plt.plot(results['date'],np.cumsum(results['forcing_precipitation']),'--k', label='forcing_precipitation')
plt.legend()





plt.figure(figsize=(20,5))
plt.plot(results['date'],results['soil_water_closure'].mean(['i','j']))