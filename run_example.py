# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 15:25:15 2019

@author: janousu
"""

from model_driver import driver
from iotools import read_results
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# example of calling driver, reading results and plotting gwl

outputfile = driver(create_ncf=True, folder='testcase_input')

#results = read_results(outputfile)
results = xr.open_dataset(outputfile)


'''
results['soil_ground_water_level_abs'] = results['soil_ground_water_level'][-1,:,:] + results['parameters_elevation']
results['soil_netflow_to_ditch'] = results['soil_netflow_to_ditch'] * results['parameters_cmask']
results['soil_lateral_netflow'] = results['soil_lateral_netflow'] * results['parameters_cmask']
results['soil_water_storage'] = results['soil_water_storage'] - results['soil_water_storage'][0,:,:]
results['bucket_water_storage'] = results['bucket_water_storage'] - results['bucket_water_storage'][0,:,:]
results['soil_water_storage'] = results['soil_water_storage'] + results['bucket_water_storage']

plt.figure(figsize=(25,12))
ax=plt.subplot(2,4,1)
results['soil_ground_water_level'][-1,:,:].plot(cmap='coolwarm_r', vmin=-4, vmax=4)
plt.subplot(2,4,2, sharex=ax, sharey=ax)
results['parameters_elevation'][:,:].plot()
plt.subplot(2,4,3, sharex=ax, sharey=ax)
results['soil_moisture_deep'][-1,:,:].plot(cmap='coolwarm_r', vmin=0, vmax=1)
plt.subplot(2,4,4, sharex=ax, sharey=ax)
results['bucket_moisture_root'][-1,:,:].plot(cmap='coolwarm_r', vmin=0, vmax=1)
plt.subplot(2,4,5, sharex=ax, sharey=ax)
results['soil_netflow_to_ditch'][-1,:,:].plot()
plt.subplot(2,4,6, sharex=ax, sharey=ax)
results['soil_lateral_netflow'][-1,:,:].plot()
plt.subplot(2,4,7, sharex=ax, sharey=ax)
results['bucket_surface_runoff'][-1,:,:].plot()
plt.subplot(2,4,8, sharex=ax, sharey=ax)
results['soil_water_closure'][-1,:,:].plot()
plt.savefig('TESTplots_part3.png')

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
plt.plot(results['date'],results['soil_netflow_to_ditch'].mean(['i','j']),label='soil_netflow_to_ditch')
plt.plot(results['date'],results['soil_netflow_to_ditch'].mean(['i','j'])
         +results['bucket_surface_runoff'].mean(['i','j']), label='+bucket_surface_runoff')
plt.legend()

plt.subplot(2,1,2,sharex=ax)
plt.plot(results['date'],results['soil_water_storage'].mean(['i','j']), label='soil+bucket_water_storage')

plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])), label='+canopy_evaporation')

plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])), label='+canopy_transpiration')

plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])), label='+bucket_evaporation')

plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])
                   + results['soil_netflow_to_ditch'].mean(['i','j'])), label='+soil_netflow_to_ditch')

plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i','j'])
                   + results['bucket_evaporation'].mean(['i','j'])
                   + results['soil_netflow_to_ditch'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])), label='+bucket_surface_runoff')

plt.plot(results['date'],np.cumsum(results['forcing_precipitation']),'--k', label='forcing_precipitation')
plt.legend()





plt.figure(figsize=(20,5))
plt.plot(results['date'],results['soil_water_closure'].mean(['i','j']))
'''
