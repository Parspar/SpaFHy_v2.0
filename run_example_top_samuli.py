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
plt.plot(tres['top_returnflow'] + tres['top_baseflow'] + tres['top_storage_change'], label='out + dS')
plt.legend()