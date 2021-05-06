# -*- coding: utf-8 -*-
"""
Created on Mon May  3 13:08:47 2021

@author: janousu
"""

# SpaFHy_v1_Pallas_2D figures

from iotools import read_results
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from datetime import date
import os

# reading the results
outputfile = 'results/testcase_input.nc'
results = read_results(outputfile)

# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\Figs')

saveplots = True
today = date.today()


#%%
# plt.figure()
# ax = plt.subplot(4,1,1)
# results['canopy_snow_water_equivalent'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,2, sharex=ax)
# results['soil_rootzone_moisture'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,3, sharex=ax)
# results['soil_ground_water_level'][:,4,:].plot.line(x='date')
# plt.subplot(4,1,4, sharex=ax)
# results['canopy_leaf_area_index'][:,4,:].plot.line(x='date')
# figures type 1
'''
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
'''

#%%
# spatial figures

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(14,14));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[2][0]
ax5 = axs[2][1]
ax6 = axs[2][2]
ax7 = axs[1][0]
ax8 = axs[1][1]
ax9 = axs[1][2]

#fig.suptitle('Spatial results and parameters')


#d = np.random.normal(.4,2,(10,10))

#fig.suptitle('Volumetric water content', fontsize=15)
#norm = mcolors.TwoSlopeNorm(vmin=d.min(), vmax = d.max(), vcenter=0)
norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_ground_water_level'][-1,:,:]), 
                            vmax=4.0, vcenter=0)
im1 = ax1.imshow(results['soil_ground_water_level'][-1,:,:], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im1, ax=ax1)
ax1.title.set_text('GW level dates[-1]')


im2 = ax2.imshow(results['soil_water_closure'][-1,:,:], vmin=-0.0000000001, vmax=0.0000000001, cmap='bwr')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('soil water closure dates[-1]')


norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_drainage'][-1,:,:]), 
                            vmax=np.nanmax(results['soil_drainage'][-1,:,:]), vcenter=0)
im3 = ax3.imshow(results['soil_drainage'][-1,:,:], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im3, ax=ax3)
ax3.title.set_text('soil drainage dates[-1]')


im4 = ax4.imshow(results['parameters_elevation'], cmap='bwr')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('DEM')

im5 = ax5.imshow(results['parameters_soilclass'], cmap='bwr')
fig.colorbar(im5, ax=ax5)
ax5.title.set_text('soil class')


im6 = ax6.imshow(results['parameters_ditches'], cmap='binary')
fig.colorbar(im6, ax=ax6)
ax6.title.set_text('ditches and ponds')


im7 = ax7.imshow(results['soil_moisture_top'][-1,:,:], cmap='coolwarm_r')
fig.colorbar(im7, ax=ax7)
ax7.title.set_text('soil moisture top dates[-1]')

im8 = ax8.imshow(results['soil_rootzone_moisture'][-1,:,:]*results['parameters_cmask'], cmap='coolwarm_r')
fig.colorbar(im8, ax=ax8)
ax8.title.set_text('soil rootzone_moisture dates[-1]')

#im9 = ax9.imshow(results['parameters_lai_conif'] + results['parameters_lai_decid_max'], cmap='coolwarm_r')
#fig.colorbar(im9, ax=ax9)
#ax9.title.set_text('LAImax')

im9 = ax9.imshow(results['parameters_sitetype'], cmap='coolwarm_r')
fig.colorbar(im9, ax=ax9)
ax9.title.set_text('sitetype')


if saveplots == True:
    plt.savefig(f'spatial_plots_{today}.pdf')
    plt.savefig(f'spatial_plots_{today}.png')
        
#%%

# Closer to kenttärova
# Plotting
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7,3));
ax1 = axs[0]
ax2 = axs[1]

norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_ground_water_level'][-1,25:100,25:100]), 
                            vmax=np.nanmax(results['soil_ground_water_level'][-1,25:100,25:100]), vcenter=0)
im1 = ax1.imshow(results['soil_ground_water_level'][-1,25:100,25:100], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im1, ax=ax1)
ax1.title.set_text('GW level dates[-1]')

norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_ground_water_level'][-1,100:175,75:150]), 
                            vmax=np.nanmax(results['soil_ground_water_level'][-1,100:175,75:150]), vcenter=0)
im2 = ax2.imshow(results['soil_ground_water_level'][-1,100:175,75:150], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('GW level dates[-1]')

#%%
# temporal figures

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]

fig.suptitle('Temporal max, mean, min plots')


im1 = ax1.plot(np.nanmean(results['soil_ground_water_level'], axis=(1,2)))
ax1.plot(np.nanmax(results['soil_ground_water_level'], axis=(1,2)))
ax1.plot(np.nanmin(results['soil_ground_water_level'], axis=(1,2)))
ax1.title.set_text('GW level')

im2 = ax2.plot(np.nanmax(results['soil_water_closure'], axis=(1,2)))
ax2.plot(np.nanmean(results['soil_water_closure'], axis=(1,2)))
ax2.plot(np.nanmin(results['soil_water_closure'], axis=(1,2)))
ax2.title.set_text('soil water closure')

im3 = ax3.plot(np.nanmean(results['soil_drainage'], axis=(1,2)))
ax3.plot(np.nanmax(results['soil_drainage'], axis=(1,2)))
ax3.plot(np.nanmin(results['soil_drainage'], axis=(1,2)))
ax3.title.set_text('soil drainage')

if saveplots == True:
    plt.savefig(f'temporal_plots_{today}.pdf')
    plt.savefig(f'temporal_plots_{today}.png')
    
#%%
'''
# spatial figures 2

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(15,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[0][3]
ax5 = axs[1][0]
ax6 = axs[1][1]
ax7 = axs[1][2]

#fig.suptitle('Spatial results and parameters')


im1 = ax1.imshow(results['soil_pond_storage'][-1,:,:], cmap='bwr')
fig.colorbar(im1, ax=ax1)
ax1.title.set_text('Soil pond storage dates[-1]')



im2 = ax2.imshow(results['soil_moisture_top'][-1,:,:], cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('soil moisture top dates[-1]')


im3 = ax3.imshow(results['soil_rootzone_moisture'][-1,:,:], cmap='coolwarm_r')
fig.colorbar(im3, ax=ax3)
ax3.title.set_text('soil rootzone_moisture dates[-1]')


im4 = ax4.imshow(results['soil_surface_runoff'][-1,:,:], cmap='bwr')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('soil surface runoff')

im5 = ax5.imshow(results['parameters_elevation'], cmap='bwr')
fig.colorbar(im5, ax=ax5)
ax5.title.set_text('DEM')

im6 = ax6.imshow(results['parameters_soilclass'], cmap='bwr')
fig.colorbar(im6, ax=ax6)
ax6.title.set_text('soil class')

im7 = ax7.imshow(results['parameters_ditches'], cmap='binary')
fig.colorbar(im7, ax=ax7)
ax7.title.set_text('ditches and ponds')


if saveplots == True:
    plt.savefig(f'spatial_plots_2__{today}.pdf')
    plt.savefig(f'spatial_plots_2__{today}.png')

        
#%%
# temporal figures

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]

fig.suptitle('Temporal max, mean, min plots')

'''