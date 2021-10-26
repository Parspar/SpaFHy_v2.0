# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 13:53:07 2021

@author: janousu
"""

from model_driver_top import driver
from iotools import read_results
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
# example of calling driver, reading results and plotting gwl

results, spa, pcpy, psoil, ptopmodel, cmask = driver(create_ncf=False, folder='testcase_input')

#with open(r'run1.pk', 'wb') as f:
#    run = (results1,spa)
#    pk.dump(run, f)
#    

#%%


results['soil_ground_water_level_abs'] = results['soil_ground_water_level'][-1,:,:] + results['parameters_elevation']
results['soil_netflow_to_ditch'] = results['soil_netflow_to_ditch'] * results['parameters_cmask']
results['soil_lateral_netflow'] = results['soil_lateral_netflow'] * results['parameters_cmask']
results['soil_water_storage'] = results['soil_water_storage'] - results['soil_water_storage'][0,:,:]
results['bucket_water_storage'] = results['bucket_water_storage'] - results['bucket_water_storage'][0,:,:]
results['soil_water_storage'] = results['soil_water_storage'] + results['bucket_water_storage']

"""
twi-distribution
"""
xi = spa.top.xi # twi
xi = xi[xi>=0] # non-nan
X = spa.top.X # mean twi

plt.figure()
plt.hist(xi, 100)

"""
topmodel water storage change. Note: also topmodel returns now storage change
"""
def dSdt(s):
    """
    s - saturation deficit (m)
    returns storage change per unit area (mm d-1)
    """
    y = np.zeros(np.shape(s))
    y[1:] = -1e3 * np.diff(s)
    return y

S = results['top_saturation_deficit']
sat_areas = []
for k in range(0, len(S)):
    s = spa.top.local_s(S[k])
#    s[s >= 0] = 1.0
#    s[s < 0] = 0.0
    sat_areas.append(s)
    
plt.figure()
tres = {}
for var in ['top_drainage_in', 'top_returnflow', 'top_baseflow', 'top_storage_change']:
    y = np.cumsum(results[var])
    plt.plot(y, label=var)
    tres[var] = y
plt.legend(); plt.ylabel('mm')

# mbe check
plt.figure()
plt.plot(tres['top_drainage_in'], 'ro', label='top_drainage_in')
plt.plot(tres['top_returnflow'] + tres['top_baseflow'] + tres['top_storage_change'])



#plt.figure(figsize=(25,15))
#ax=plt.subplot(2,4,1)
#results['soil_ground_water_level'][-1,:,:].plot()
#plt.subplot(2,4,2, sharex=ax, sharey=ax)
#results['parameters_elevation'][:,:].plot()
#plt.subplot(2,4,3, sharex=ax, sharey=ax)
#results['soil_moisture_deep'][-1,:,:].plot()
#plt.subplot(2,4,4, sharex=ax, sharey=ax)
#results['bucket_moisture_root'][-1,:,:].plot()
#plt.subplot(2,4,5, sharex=ax, sharey=ax)
#results['soil_netflow_to_ditch'][-1,:,:].plot()
#plt.subplot(2,4,6, sharex=ax, sharey=ax)
#results['soil_lateral_netflow'][-1,:,:].plot()
#plt.subplot(2,4,7, sharex=ax, sharey=ax)
#results['bucket_surface_runoff'][-1,:,:].plot()
#plt.subplot(2,4,8, sharex=ax, sharey=ax)
#results['soil_water_closure'][-1,:,:].plot()

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
#
#plt.figure(figsize=(20,8))
#ax=plt.subplot(2,1,1)
#plt.plot(results['date'],results['soil_netflow_to_ditch'].mean(['i','j']),label='soil_netflow_to_ditch')
#plt.plot(results['date'],results['soil_netflow_to_ditch'].mean(['i','j'])
#         +results['bucket_surface_runoff'].mean(['i','j']), label='+bucket_surface_runoff')
#plt.legend()
#
#plt.subplot(2,1,2,sharex=ax)
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j']), label='soil+bucket_water_storage')
#
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
#         np.cumsum(results['canopy_evaporation'].mean(['i','j'])), label='+canopy_evaporation')
#
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
#         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
#                   + results['canopy_transpiration'].mean(['i','j'])), label='+canopy_transpiration')
#
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
#         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
#                   + results['canopy_transpiration'].mean(['i','j'])
#                   + results['bucket_evaporation'].mean(['i','j'])), label='+bucket_evaporation')
#
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
#         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
#                   + results['canopy_transpiration'].mean(['i','j'])
#                   + results['bucket_evaporation'].mean(['i','j'])
#                   + results['soil_netflow_to_ditch'].mean(['i','j'])), label='+soil_netflow_to_ditch')
#
#plt.plot(results['date'],results['soil_water_storage'].mean(['i','j'])+
#         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
#                   + results['canopy_transpiration'].mean(['i','j'])
#                   + results['bucket_evaporation'].mean(['i','j'])
#                   + results['soil_netflow_to_ditch'].mean(['i','j'])
#                   + results['bucket_surface_runoff'].mean(['i','j'])), label='+bucket_surface_runoff')
#
#plt.plot(results['date'],np.cumsum(results['forcing_precipitation']),'--k', label='forcing_precipitation')
#plt.legend()



