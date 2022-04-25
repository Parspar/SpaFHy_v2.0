# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:33:11 2022

@author: janousu
"""

# BUCKET WB

# IN
# bucket infiltration + bucket_returnflow
# OUT
# bucket_drainage + bucket_evaporation + bucket_surface_runoff

# dS
# bucket_water_storage

results['bucket_water_storage'] = results['bucket_water_storage'] - results['bucket_water_storage'][0,:,:]
results['canopy_water_storage'] = results['canopy_water_storage'] - results['canopy_water_storage'][0,:,:]
results['top_water_storage'] = results['top_storage_change']


#%%

# BUCKET MBE

plt.figure(figsize=(20,8))

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j']), label='bucket_water_storage')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['bucket_evaporation'].mean(['i','j'])), label='+bucket_evaporation')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['bucket_evaporation'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])), label='+bucket_surface_runoff')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['bucket_evaporation'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])
                   + results['bucket_drainage'].mean(['i','j'])), label='+bucket_drainage')

plt.plot(results['date'],results['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results['bucket_evaporation'].mean(['i','j'])
                   + results['bucket_surface_runoff'].mean(['i','j'])
                   + results['bucket_drainage'].mean(['i','j'])
                   + results['canopy_transpiration'].mean(['i', 'j'])), label='+canopy_transpiration')

plt.plot(results['date'],np.cumsum(results['bucket_potential_infiltration'].mean(['i', 'j'])) + np.cumsum(results['bucket_return_flow'].mean(['i', 'j'])),
         '--k',alpha=0.5, label='bucket_infiltration + bucket_return_flow')
plt.legend()

#%%

# CANOPY MBE

plt.figure(figsize=(20,8))

plt.plot(results['date'],results['canopy_water_storage'].mean(['i','j']), label='canopy_water_storage')

plt.plot(results['date'],results['canopy_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])), label='+canopy_evaporation')

plt.plot(results['date'],results['canopy_water_storage'].mean(['i','j'])+
         np.cumsum(results['canopy_evaporation'].mean(['i','j'])
                   + results['bucket_potential_infiltration'].mean(['i','j'])), label='+bucket_infiltration')

plt.plot(results['date'],np.cumsum(results['forcing_precipitation']),'--k', label='forcing_precipitation')

plt.legend()

#%%

# TOP MBE

plt.figure(figsize=(20,8))

plt.plot(results['date'],results['top_water_storage'], label='top_water_storage')

plt.plot(results['date'],results['top_water_storage'] +
         np.cumsum(results['top_baseflow']), label='+top_baseflow')

plt.plot(results['date'],results['top_water_storage'] +
         np.cumsum(results['top_baseflow']) +
                   + results['top_returnflow'], label='+top_returnflow')

plt.plot(results['date'],np.cumsum(results['bucket_drainage'].mean(['i', 'j'])),'--k', label='bucket_drainage')

plt.legend()