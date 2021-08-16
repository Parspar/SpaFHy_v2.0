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
from netCDF4 import Dataset #, date2num
import pandas as pd
import pickle
import seaborn as sns
from iotools import read_AsciiGrid


# reading the results
outputfile = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_202108161216.nc'
results = read_results(outputfile)


sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_scd.nc'
sar = Dataset(sar_path, 'r')

dates_spa = []
for d in range(len(results['date'])):
    dates_spa.append(pd.to_datetime(str(results['date'][d])[36:46]))


# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

saveplots = True
today = date.today()

cmask = results['parameters_cmask']
sitetype = results['parameters_sitetype']
soilclass = results['parameters_soilclass']

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

plt_ind = 20

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
norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_ground_water_level'][plt_ind,:,:]), 
                            vmax=4.0, vcenter=0)
im1 = ax1.imshow(results['soil_ground_water_level'][plt_ind,:,:], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im1, ax=ax1)
ax1.title.set_text(f'GW level dates[{plt_ind}]')


im2 = ax2.imshow(results['soil_water_closure'][plt_ind,:,:], vmin=-0.0000000001, vmax=0.0000000001, cmap='bwr')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text(f'soil water closure dates[{plt_ind}]')


norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_lateral_netflow'][plt_ind,:,:]), 
                            vmax=np.nanmax(results['soil_lateral_netflow'][plt_ind,:,:]), vcenter=0)
im3 = ax3.imshow(results['soil_lateral_netflow'][-1,:,:], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im3, ax=ax3)
ax3.title.set_text(f'soil_lateral_netflow dates[{plt_ind}]')


im4 = ax4.imshow(results['parameters_elevation'], cmap='bwr')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('DEM')

im5 = ax5.imshow(results['parameters_soilclass'], cmap='bwr')
fig.colorbar(im5, ax=ax5)
ax5.title.set_text('soil class')


im6 = ax6.imshow(results['parameters_ditches'], cmap='binary')
fig.colorbar(im6, ax=ax6)
ax6.title.set_text('ditches and ponds')


im7 = ax7.imshow(results['bucket_moisture_top'][plt_ind,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im7, ax=ax7)
ax7.title.set_text(f'soil moisture top dates[{plt_ind}]')

im8 = ax8.imshow(results['bucket_moisture_root'][plt_ind,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im8, ax=ax8)
ax8.title.set_text(f'bucket_moisture_root dates[{plt_ind}]')

#im9 = ax9.imshow(results['parameters_lai_conif'] + results['parameters_lai_decid_max'], cmap='coolwarm_r')
#fig.colorbar(im9, ax=ax9)
#ax9.title.set_text('LAImax')

im9 = ax9.imshow(results['parameters_sitetype'], cmap='coolwarm_r')
fig.colorbar(im9, ax=ax9)
ax9.title.set_text('sitetype')

fig.suptitle('SpaFHy v2D')


if saveplots == True:
    plt.savefig(f'spatial_plots_{today}.pdf')
    plt.savefig(f'spatial_plots_{today}.png')
        

#%%
# temporal figures
'''
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
'''
    
#%%
'''
# temporal figures

# Plotting
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(12,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]
ax4 = axs[3]


im1 = ax1.plot(results['soil_rootzone_moisture'][:,50,50])
ax1.plot(results['soil_rootzone_moisture'][:,220,80])
ax1.title.set_text('soil moisture rootzone mire')
ax1.set_ylim([0,1])

im2 = ax2.plot(results['soil_moisture_top'][:,50,50])
ax2.plot(results['soil_moisture_top'][:,220,80])
ax2.title.set_text('soil moisture top mire')
ax2.set_ylim([0,1])

im3 = ax3.plot(results['soil_rootzone_moisture'][:,90,70])
ax3.plot(results['soil_rootzone_moisture'][:,125,125])
ax3.title.set_text('soil moisture rootzone mineral')
ax3.set_ylim([0,1])

im4 = ax4.plot(results['soil_moisture_top'][:,90,70])
ax4.plot(results['soil_moisture_top'][:,125,125])
ax4.title.set_text('soil moisture top mineral')
ax4.set_ylim([0,1])

fig.suptitle('SpaFHy v2D')
'''


#%%

# soil moist plots

# spafhy
wet_day = np.nansum(results['bucket_moisture_root'], axis=(1,2)).argmax()
dry_day = np.nansum(results['bucket_moisture_root'], axis=(1,2)).argmin()

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(18,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]


im1 = ax1.imshow(results['bucket_moisture_top'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im1, ax=ax1)
ax1.title.set_text('soil moisture top wet')

im2 = ax2.imshow(results['bucket_moisture_top'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('soil moisture top dry')

im3 = ax3.imshow(results['bucket_moisture_root'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im3, ax=ax3)
ax3.title.set_text('soil moisture root wet')

im4 = ax4.imshow(results['bucket_moisture_root'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('soil moisture root dry')

im5 = ax5.imshow(results['soil_ground_water_level'][wet_day,:,:], vmin=-4.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im5, ax=ax5)
ax5.title.set_text('soil ground water level root wet')

im6 = ax6.imshow(results['soil_ground_water_level'][dry_day,:,:], vmin=-4.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im6, ax=ax6)
ax6.title.set_text('soil ground water level root dry')


#%%

# sar soil moisture plots
import pandas as pd

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
kenttarova_loc = list([int(kenttarova_loc[0]), int(kenttarova_loc[1])])

# reading sar data
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_mean4.nc'
#sar = Dataset(sar_path, 'r')

sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = results['soil_rootzone_moisture']
spa_wliq_top = results['soil_moisture_top']

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 
dates_spa = pd.to_datetime(dates_spa, format='%Y%m%d') 

#spa dates to match sar dates
date_in_spa = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]

# driest and wettest days
spasum = np.nansum(spa_wliq, axis=(1,2))

# index in sar data
day_low = int(np.where(spasum == np.nanmin(spasum))[0])
day_hi = int(np.where(spasum == np.nanmax(spasum))[0])
#sar_low = 43
#sar_hi = 20
# day in sar data
low_date = dates_sar[day_low].strftime("%Y-%m-%d")
hi_date = dates_sar[day_hi].strftime("%Y-%m-%d")

# cropping for plots
xcrop = np.arange(20,170)
ycrop = np.arange(20,250)
sar_wliq = sar_wliq[:,ycrop,:]
sar_wliq = sar_wliq[:,:,xcrop]
spa_wliq = spa_wliq[:,ycrop,:]
spa_wliq = spa_wliq[:,:,xcrop]
spa_wliq_top = spa_wliq_top[:,ycrop,:]
spa_wliq_top = spa_wliq_top[:,:,xcrop]

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.imshow(sar_wliq[day_hi,:,:], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax1.title.set_text('SAR')
#ax1.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im2 = ax2.imshow(spa_wliq[day_hi, :,:], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax2.title.set_text('SPAFHY rootzone')
ax2.text(10, -15, f'Wet day : {hi_date}', fontsize=15)
#ax2.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im3 = ax3.imshow(sar_wliq[day_low, :,:], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax3.title.set_text('SAR')

im4 = ax4.imshow(spa_wliq[day_low, :,:], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax4.title.set_text('SPAFHY rootzone')
ax4.text(10, -15, f'Dry day : {low_date}', fontsize=15)


im5 = ax5.imshow(spa_wliq_top[day_hi, :, :], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax5.title.set_text('SPAFHY topsoil')

im6 = ax6.imshow(spa_wliq_top[day_low, :,:], cmap='coolwarm_r', vmin=0.0, vmax=1.0, aspect='equal')
ax6.title.set_text('SPAFHY topsoil')

ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.015, 0.7])
bar1 = fig.colorbar(im1, cax=cbar_ax)

fig.suptitle('SpaFHy v2D')

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_{today}.png')
    
#%%

# poikki ja pituusleikkaukset

x = 100
y = 90
sar_long = np.array(sar_wliq[:,:,x])
spa_long = np.array(spa_wliq[:,:,x])
spa_top_long = np.array(spa_wliq_top[:,:,x])

sar_cross = np.array(sar_wliq[:,y,:])
spa_cross = np.array(spa_wliq[:,y,:])
spa_top_cross = np.array(spa_wliq_top[:,y,:])


# Plotting
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(16,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[2][0]
ax6 = axs[2][1]

ax1.plot(sar_long[day_hi])
ax1.plot(spa_long[day_hi])
ax1.plot(spa_top_long[day_hi])
ax1.legend(['SAR', 'SPA', 'SPATOP'], ncol = 3)
ax1.title.set_text('wet day')


ax2.plot(sar_long[day_low])
ax2.plot(spa_long[day_low])
ax2.plot(spa_top_long[day_low])
ax2.title.set_text('dry day')


ax3.plot(sar_cross[day_hi])
ax3.plot(spa_cross[day_hi])
ax3.plot(spa_top_cross[day_hi])
#ax3.legend(['SAR', 'SPA', 'SPATOP'], ncol = 3)

ax4.plot(sar_cross[day_low])
ax4.plot(spa_cross[day_low])
ax4.plot(spa_top_cross[day_low])

ax1.text(200, 0.85, f'Cross x={x}', fontsize=12)
ax3.text(133, 0.85, f'Cross y={y}', fontsize=12)

ax5.imshow(sar_wliq[day_hi], cmap='coolwarm_r')
ax6.imshow(sar_wliq[day_low], cmap='coolwarm_r')
ax5.axvline(x=x,color='red')
ax5.axhline(y=y,color='red')
ax6.axvline(x=x,color='red')
ax6.axhline(y=y,color='red')
ax5.title.set_text('SAR')
ax6.title.set_text('SAR')

#ax5.axis("off")
#ax6.axis("off")

fig.suptitle('SpaFHy v2D')


if saveplots == True:
    plt.savefig(f'SAR_vs_SPAFHY_soilmoist_cross_{today}.pdf')
    plt.savefig(f'SAR_vs_SPAFHY_soilmoist_cross_{today}.png')

#%%
# check with 17 6 21 measurements

# sar soil moisture plots
import pandas as pd

# cell locations of kenttarova
measured, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\SOIL\GIS\soilmoist_170621.asc')
sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100

sar_wliq_mean = np.nanmean(sar_wliq[8:22], axis=(0))
sar_wliq_mean_flat = sar_wliq_mean.flatten()
measured_flat = measured.flatten()

compar = pd.DataFrame()
compar['meas'] = measured_flat[np.isfinite(measured_flat)]
compar['sar'] = sar_wliq_mean_flat[np.isfinite(measured_flat)]

# Plotting
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4));

ax.plot(compar['meas'])
ax.plot(compar['sar'])
ax.legend(['meas', 'SAR'])

#%%
'''
from matplotlib import gridspec

# plot it
fig = plt.figure(figsize=(8, 6))

gs = gridspec.GridSpec(6, 6)

ax0 = plt.subplot(gs[0,0])
ax0.plot(x, y)
ax1 = plt.subplot(gs[1,1])
ax1.plot(y, x)
ax2 = plt.subplot(gs[0,2])
ax2.plot(y, x)

ax3 = plt.subplot(gs[1,3])
ax3.plot(x, y)
ax4 = plt.subplot(gs[2:4,1])
ax4.plot(y, x)
ax5 = plt.subplot(gs[4:6,0])
ax5.plot(y, x)

ax6 = plt.subplot(gs[4:6,1])
#ax6.plot(x, y)
ax7 = plt.subplot(gs[4:6,2])
#ax7.plot(y, x)


plt.tight_layout()
plt.savefig('grid_figure.png')

plt.show()
'''

#%%

# GW levels low and high day
gwsum = np.nansum(results['soil_ground_water_level'], axis=(1,2))
hi = np.where(gwsum == gwsum.max())[0][0]
lo = np.where(gwsum == gwsum.min())[0][0]

#norm = mcolors.TwoSlopeNorm(vmin=results['soil_ground_water_level'][lo,:,:].min(), vmax = results['soil_ground_water_level'][hi,:,:].max(), vcenter=0)
norm = mcolors.TwoSlopeNorm(vmin=-4, vcenter=0, vmax = 4)

# Plotting
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,12));
ax1 = axs[0]
ax2 = axs[1]

im1 = ax1.imshow(results['soil_ground_water_level'][hi,:,:], cmap=plt.cm.RdBu, norm=norm)
ax2.imshow(results['soil_ground_water_level'][lo,:,:], cmap=plt.cm.RdBu, norm=norm)
ax1.title.set_text(f'GW level high {dates_spa[hi]}')
ax2.title.set_text(f'GW level low {dates_spa[lo]}')

ax1.axis("off")
ax2.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.35, 0.015, 0.3])
bar1 = fig.colorbar(im1, cax=cbar_ax)

if saveplots == True:
    plt.savefig(f'GW_levels_hi_low_{today}.pdf')
    plt.savefig(f'GW_levels_hi_low_{today}.png')

#%%

# Comparing GW level with two sets of results


norm = mcolors.TwoSlopeNorm(vmin=-4, vcenter=0, vmax = 2)
norm = mcolors.TwoSlopeNorm(vmin=results_raw['soil_ground_water_level'][lo,:,:].min(), vmax = results_raw['soil_ground_water_level'][hi,:,:].max(), vcenter=0)

raw_sum = np.nansum(results_raw['soil_ground_water_level'][hi,:,:]) + np.nansum(results_raw['soil_ground_water_level'][lo,:,:])
d4_sum = np.nansum(results_d4['soil_ground_water_level'][hi,:,:]) + np.nansum(results_d4['soil_ground_water_level'][lo,:,:])
d8_sum = np.nansum(results_d8['soil_ground_water_level'][hi,:,:]) + np.nansum(results_d8['soil_ground_water_level'][lo,:,:])

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12,16));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[2][0]
ax6 = axs[2][1]

im1 = ax1.imshow(results_raw['soil_ground_water_level'][hi,:,:], cmap=plt.cm.RdBu, norm=norm)
ax1.title.set_text(f'GW level high {dates_spa[hi]}')

ax2.imshow(results_raw['soil_ground_water_level'][lo,:,:], cmap=plt.cm.RdBu, norm=norm)
ax2.title.set_text(f'GW level low {dates_spa[lo]}')

ax3.imshow(results_d4['soil_ground_water_level'][hi,:,:], cmap=plt.cm.RdBu, norm=norm)
ax3.title.set_text(f'GW level high {dates_spa[hi]}')

ax4.imshow(results_d4['soil_ground_water_level'][lo,:,:], cmap=plt.cm.RdBu, norm=norm)
ax4.title.set_text(f'GW level low {dates_spa[lo]}')

ax5.imshow(results_d8['soil_ground_water_level'][hi,:,:], cmap=plt.cm.RdBu, norm=norm)
ax5.title.set_text(f'GW level high {dates_spa[hi]}')

ax6.imshow(results_d8['soil_ground_water_level'][lo,:,:], cmap=plt.cm.RdBu, norm=norm)
ax6.title.set_text(f'GW level low {dates_spa[lo]}')

ax1.text(200, -15, 'RAW DEM', fontsize=15)
ax3.text(200, -15, 'D4 FILL DEM', fontsize=15)
ax5.text(200, -15, 'D8 FILL DEM', fontsize=15)


ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.35, 0.015, 0.3])
bar1 = fig.colorbar(im1, cax=cbar_ax)

if saveplots == True:
    plt.savefig(f'GW_levels_hi_low_comp_dems_{today}.pdf')
    plt.savefig(f'GW_levels_hi_low_comp_dems_{today}.png')

#%%

# normalized by mean plots

# mean of each pixel
#spamean = np.nanmean(spa_wliq, axis=0)
#spatopmean = np.nanmean(spa_wliq_top, axis=0)
#sarmean = np.nanmean(sar_wliq, axis=0)

# mean of wet and dry days
spamean_wet = np.nanmean(spa_wliq[day_hi,:,:])
spatopmean_wet = np.nanmean(spa_wliq_top[day_hi,:,:])
sarmean_wet = np.nanmean(sar_wliq[day_hi,:,:])
spamean_dry = np.nanmean(spa_wliq[day_low,:,:])
spatopmean_dry = np.nanmean(spa_wliq_top[day_low,:,:])
sarmean_dry = np.nanmean(sar_wliq[day_low,:,:])

# median day of total sums
spasum = np.nansum(spa_wliq, axis=(1,2))
spamedian = spa_wliq[np.where(np.sort(np.nansum(spa_wliq, axis=(1,2)))[(int(len(spasum)/2))] == spasum)[0][0],:,:]
sarsum = np.nansum(sar_wliq, axis=(1,2))
sarmedian = sar_wliq[np.where(np.sort(np.nansum(sar_wliq, axis=(1,2)))[(int(len(sarsum)/2))] == sarsum)[0][0],:,:]
spatopsum = np.nansum(spa_wliq_top, axis=(1,2))
spatopmedian = spa_wliq[np.where(np.sort(np.nansum(spa_wliq_top, axis=(1,2)))[(int(len(spatopsum)/2))] == spatopsum)[0][0],:,:]

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.imshow(sar_wliq[day_hi,:,:]/sarmean_wet, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax1.title.set_text('SAR/SARwet_mean')
#ax1.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im2 = ax2.imshow(spa_wliq[day_hi, :,:]/spamean_wet, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax2.title.set_text('SPAFHYROOT/SPAFHYROOTwet_mean')
ax2.text(10, -15, f'Wet day : {hi_date}', fontsize=15)
#ax2.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im3 = ax3.imshow(sar_wliq[day_low, :,:]/sarmean_dry, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax3.title.set_text('SAR/SARdry_mean')

im4 = ax4.imshow(spa_wliq[day_low, :,:]/spamean_dry, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax4.title.set_text('SPAFHYROOT/SPAFHYROOTdry_mean')
ax4.text(10, -15, f'Dry day : {low_date}', fontsize=15)


im5 = ax5.imshow(spa_wliq_top[day_hi, :, :]/spatopmean_wet, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax5.title.set_text('SPAFHYTOP/SPAFHYTOPwet_mean')

im6 = ax6.imshow(spa_wliq_top[day_low, :,:]/spatopmean_dry, cmap='coolwarm_r', vmin=0.0, vmax=2.0, aspect='equal')
ax6.title.set_text('SPAFHYTOP/SPAFHYTOPdry_mean')

ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.015, 0.7])
bar1 = fig.colorbar(im1, cax=cbar_ax)

#ax1.text(10, -15, 'norm by mean of the day', fontsize=10)


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normmean_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normmean_{today}.png')

#%%

# normalized by median day plots

# mean of each pixel
#spamean = np.nanmean(spa_wliq, axis=0)
#spatopmean = np.nanmean(spa_wliq_top, axis=0)
#sarmean = np.nanmean(sar_wliq, axis=0)

# median day of total sums
spasum = np.nansum(spa_wliq, axis=(1,2))
spamedian = spa_wliq[np.where(np.sort(np.nansum(spa_wliq, axis=(1,2)))[(int(len(spasum)/2))] == spasum)[0][0],:,:]
sarsum = np.nansum(sar_wliq, axis=(1,2))
sarmedian = sar_wliq[np.where(np.sort(np.nansum(sar_wliq, axis=(1,2)))[(int(len(sarsum)/2))] == sarsum)[0][0],:,:]
spatopsum = np.nansum(spa_wliq_top, axis=(1,2))
spatopmedian = spa_wliq[np.where(np.sort(np.nansum(spa_wliq_top, axis=(1,2)))[(int(len(spatopsum)/2))] == spatopsum)[0][0],:,:]

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.imshow(sar_wliq[day_hi,:,:]/sarmedian, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax1.title.set_text('SAR/SAR_median')
#ax1.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im2 = ax2.imshow(spa_wliq[day_hi, :,:]/spamedian, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax2.title.set_text('SPAFHYROOT/SPAFHYROOT_median')
ax2.text(10, -15, f'Wet day : {hi_date}', fontsize=15)
#ax2.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im3 = ax3.imshow(sar_wliq[day_low, :,:]/sarmedian, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax3.title.set_text('SAR/SAR_median')

im4 = ax4.imshow(spa_wliq[day_low, :,:]/spamedian, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax4.title.set_text('SPAFHYROOT/SPAFHYROOT_median')
ax4.text(10, -15, f'Dry day : {low_date}', fontsize=15)


im5 = ax5.imshow(spa_wliq_top[day_hi, :, :]/spatopmedian, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax5.title.set_text('SPAFHYTOP/SPAFHYTOP_median')

im6 = ax6.imshow(spa_wliq_top[day_low, :,:]/spatopmedian, cmap='coolwarm_r', vmin=0.0, vmax=2.0, aspect='equal')
ax6.title.set_text('SPAFHYTOP/SPAFHYTOP_median')

ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.015, 0.7])
bar1 = fig.colorbar(im1, cax=cbar_ax)

#ax1.text(10, -15, 'norm by mean of the day', fontsize=10)


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normmedian_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normmedian_{today}.png')

#%%

# normalized by pixel mean

# mean of each pixel
spamean = np.nanmean(spa_wliq, axis=0)
spatopmean = np.nanmean(spa_wliq_top, axis=0)
sarmean = np.nanmean(sar_wliq, axis=0)

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.imshow(sar_wliq[day_hi,:,:]/sarmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax1.title.set_text('SAR/SARpixel_mean')
#ax1.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im2 = ax2.imshow(spa_wliq[day_hi, :,:]/spamean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax2.title.set_text('SPAFHYROOT/SPAFHYROOTpixel_mean')
ax2.text(10, -15, f'Wet day : {hi_date}', fontsize=15)
#ax2.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im3 = ax3.imshow(sar_wliq[day_low, :,:]/sarmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax3.title.set_text('SAR/SARpixel_mean')

im4 = ax4.imshow(spa_wliq[day_low, :,:]/spamean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax4.title.set_text('SPAFHYROOT/SPAFHYROOTpixel_mean')
ax4.text(10, -15, f'Dry day : {low_date}', fontsize=15)


im5 = ax5.imshow(spa_wliq_top[day_hi, :, :]/spatopmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax5.title.set_text('SPAFHYTOP/SPAFHYTOPpixel_mean')

im6 = ax6.imshow(spa_wliq_top[day_low, :,:]/spatopmean, cmap='coolwarm_r', vmin=0.0, vmax=2.0, aspect='equal')
ax6.title.set_text('SPAFHYTOP/SPAFHYTOPpixel_mean')

ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")

plt.tight_layout()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.015, 0.7])
bar1 = fig.colorbar(im1, cax=cbar_ax)

#ax1.text(10, -15, 'norm by mean of the day', fontsize=10)


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normpixel_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_normpixel_{today}.png')


#%%

# point examples from mineral and openmire
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5.nc'
#sar = Dataset(sar_path, 'r')

# soilscouts at Kenttarova
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
soil_file = 'soilscouts_s3_s5_s18.csv'
fp = os.path.join(folder, soil_file)
soilscout = pd.read_csv(fp, sep=';', date_parser=['time'])
soilscout['time'] = pd.to_datetime(soilscout['time'])

# ec observation data
ec_fp = r'C:\SpaFHy_v1_Pallas\data\obs\ec_soilmoist.csv'
ecmoist = pd.read_csv(ec_fp, sep=';', date_parser=['time'])
ecmoist['time'] = pd.to_datetime(ecmoist['time'])

soilm = soilscout.merge(ecmoist)

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[0]), int(kenttarova_loc[1])])
l_loc = [60, 60]
sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = results['soil_rootzone_moisture']
spa_wliq_top = results['soil_moisture_top']

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 


#spa dates to match sar dates
date_in_spa = []
date_in_soilm = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
    yx = np.where(soilm['time'] == dates_sar[i])[0][0]
    date_in_soilm.append(yx)
   
    
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]
soilm = soilm.loc[date_in_soilm]
soilm = soilm.reset_index()

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,8));
ax1 = axs[0]
ax2 = axs[1]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.plot(sar_wliq[:,k_loc[0],k_loc[1]])
#ax1.plot(spa_wliq[:,k_loc[0],k_loc[1]])
#ax1.plot(spa_wliq_top[:,k_loc[0],k_loc[1]])
ax1.plot(soilm['s3'])
ax1.plot(soilm['s5'])
ax1.plot(soilm['s18'])
ax1.title.set_text('Mineral')
ax1.legend(['SAR', 'SpaFHy rootzone', 'SpaFHy top', 's3 = -0.05m', 's5 = -0.60', 's18 = -0.3'], ncol = 6)


im2 = ax2.plot(sar_wliq[:,l_loc[0],l_loc[1]])
ax2.plot(spa_wliq[:,l_loc[0],l_loc[1]])
ax2.plot(spa_wliq_top[:,l_loc[0],l_loc[1]])
ax2.title.set_text('Open mire')
ax2.legend(['SAR', 'SpaFHy rootzone', 'SpaFHy top'], ncol = 3)

#%%

# point examples from mineral and openmire without SAR

# soilscouts at Kenttarova
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
soil_file = 'soilscouts_s3_s5_s18.csv'
fp = os.path.join(folder, soil_file)
soilscout = pd.read_csv(fp, sep=';', date_parser=['time'])
soilscout['time'] = pd.to_datetime(soilscout['time'])

# ec observation data
ec_fp = r'C:\SpaFHy_v1_Pallas\data\obs\ec_soilmoist.csv'
ecmoist = pd.read_csv(ec_fp, sep=';', date_parser=['time'])
ecmoist['time'] = pd.to_datetime(ecmoist['time'])

soilm = soilscout.merge(ecmoist)

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[0]), int(kenttarova_loc[1])])
l_loc = [60, 60]
sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = results['soil_rootzone_moisture']
spa_wliq_top = results['soil_moisture_top']

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 

'''
#spa dates to match sar dates
date_in_spa = []
date_in_soilm = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
    yx = np.where(soilm['time'] == dates_sar[i])[0][0]
    date_in_soilm.append(yx)
'''
spa_wliq_df = pd.DataFrame()
spa_wliq_df['spa_k'] = spa_wliq[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l'] = spa_wliq[:,l_loc[0],l_loc[1]]
spa_wliq_df['spatop_k'] = spa_wliq_top[:,k_loc[0],k_loc[1]]
spa_wliq_df['spatop_l'] = spa_wliq_top[:,l_loc[0],l_loc[1]]
spa_wliq_df['time'] = dates_spa

soilm = soilm.merge(spa_wliq_df)
soilm.index = soilm['time']
soilm = soilm[['s3', 's5', 's18', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B', 'spa_k', 'spa_l', 'spatop_k', 'spatop_l']]

soilm = soilm.loc[(soilm.index > '2018-04-01') & (soilm.index < '2019-12-01')]


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,8));
ax1 = axs[0]
ax2 = axs[1]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.plot(soilm['spa_k'])
ax1.plot(soilm['spatop_k'])
ax1.plot(soilm['s3'], alpha=0.4)
#ax1.plot(soilm['s5'], alpha=0.4)
ax1.plot(soilm['s18'], alpha=0.4)
#ax1.plot(soilm['SH-5A'], alpha=0.4)
ax1.plot(soilm['SH-5B'], alpha=0.4)
#ax1.plot(soilm['SH-20A'], alpha=0.4)
ax1.plot(soilm['SH-20B'], alpha=0.4)
ax1.title.set_text('Mineral')
ax1.legend(['spa root', 'spa top', 's3 = -0.05', 's18 = -0.3', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B'], ncol = 8)
ax1.set_ylim(0,0.6)

im2 = ax2.plot(soilm['spa_l'])
ax2.plot(soilm['spatop_l'])
ax2.title.set_text('Mire')
ax2.legend(['spa root', 'spa top'], ncol = 2)

fig.suptitle('SpaFHy v2D')


if saveplots == True:
    plt.savefig(f'pointplots_soilmoist_{today}.pdf')
    plt.savefig(f'pointplots_soilmoist_{today}.png')

#%%

# Q-Q plots of dry, wet and inbetween day

import numpy as np 
import pandas as pd

from iotools import read_AsciiGrid

norm = False

sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = np.array(results['soil_rootzone_moisture'])
spa_wliq_top = np.array(results['soil_moisture_top'])

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 


#spa dates to match sar dates
date_in_spa = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]

if norm == True:
    spa_wliq = spa_wliq/(np.nanmean(spa_wliq))
    spa_wliq_top = spa_wliq_top/(np.nanmean(spa_wliq_top))
    sar_wliq = sar_wliq/(np.nanmean(sar_wliq))

sar_flat_dry = sar_wliq[day_low,:,:].flatten()
#sar_flat[np.where(sar_flat <= 0)] = np.nan
spa_flat_dry = spa_wliq[day_low,:,:].flatten()
spa_top_flat_dry = spa_wliq_top[day_low,:,:].flatten()

flat_pd = pd.DataFrame()
flat_pd['sar_dry'] = sar_flat_dry
flat_pd['spa_dry'] = spa_flat_dry
flat_pd['spa_top_dry'] = spa_top_flat_dry

sar_flat_wet = sar_wliq[day_hi,:,:].flatten()
#sar_flat[np.where(sar_flat <= 0)] = np.nan
spa_flat_wet = spa_wliq[day_hi,:,:].flatten()
spa_top_flat_wet = spa_wliq_top[day_hi,:,:].flatten()

inb = int((day_hi+day_low)/2)
sar_flat_inb = sar_wliq[inb,:,:].flatten()
spa_flat_inb = spa_wliq[inb,:,:].flatten()
spa_top_flat_inb = spa_wliq_top[inb,:,:].flatten()

flat_pd['sar_wet'] = sar_flat_wet
flat_pd['spa_wet'] = spa_flat_wet
flat_pd['spa_top_wet'] = spa_top_flat_wet

flat_pd['sar_inb'] = sar_flat_inb
flat_pd['spa_inb'] = spa_flat_inb
flat_pd['spa_top_inb'] = spa_top_flat_inb

flat_pd = flat_pd.loc[np.isfinite(flat_pd['sar_dry']) & np.isfinite(flat_pd['spa_dry']) & np.isfinite(flat_pd['spa_top_dry'])]
#flat_pd = flat_pd.loc[(flat_pd['sar'] > 0) & (flat_pd['sar'] < 1)]

#g = sns.scatterplot(flat_pd['sar'], flat_pd['spa'], alpha=0.0001, s=2)
#g.set(ylim=(-0.1, 1.0))
#g.set(xlim=(-0.1, 1.0))


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[0][2]
ax6 = axs[1][2]

x1 = sns.regplot(ax=ax1, x=flat_pd['sar_dry'], y=flat_pd['spa_dry'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax1.set(ylim=(0, 1))
    ax1.set(xlim=(0, 1))
else:
    ax1.set(ylim=(0, 2.5))
    ax1.set(xlim=(0, 2.5))    
#ax1.set_title('Dry day')

x2 = sns.regplot(ax=ax2, x=flat_pd['sar_wet'], y=flat_pd['spa_wet'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax2.set(ylim=(0, 1))
    ax2.set(xlim=(0, 1))
else:
    ax2.set(ylim=(0, 2.5))
    ax2.set(xlim=(0, 2.5)) 
#ax2.set_title('Wet day')

x3 = sns.regplot(ax=ax3, x=flat_pd['sar_dry'], y=flat_pd['spa_top_dry'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax3.set(ylim=(0, 1))
    ax3.set(xlim=(0, 1))
else:
    ax3.set(ylim=(0, 2.5))
    ax3.set(xlim=(0, 2.5))     
#ax3.set_title('Dry day')

x4 = sns.regplot(ax=ax4, x=flat_pd['sar_wet'], y=flat_pd['spa_top_wet'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax4.set(ylim=(0, 1))
    ax4.set(xlim=(0, 1))
else:
    ax4.set(ylim=(0, 2.5))
    ax4.set(xlim=(0, 2.5))     
#ax4.set_title('Wet day')

x5 = sns.regplot(ax=ax5, x=flat_pd['sar_inb'], y=flat_pd['spa_inb'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax5.set(ylim=(0, 1))
    ax5.set(xlim=(0, 1))
else:
    ax5.set(ylim=(0, 2.5))
    ax5.set(xlim=(0, 2.5))     
#ax4.set_title('Wet day')

x6 = sns.regplot(ax=ax6, x=flat_pd['sar_inb'], y=flat_pd['spa_top_inb'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
if norm == False:
    ax6.set(ylim=(0, 1))
    ax6.set(xlim=(0, 1))
else:
    ax6.set(ylim=(0, 2.5))
    ax6.set(xlim=(0, 2.5))     
#ax4.set_title('Wet day')

if norm == True:
    fig.suptitle('SpaFHy v2D, norm by total means of each')
else:
    fig.suptitle('SpaFHy v2D')

if norm == False:
    if saveplots == True:
        plt.savefig(f'sar_spa_qq_wetdry_{today}.pdf')
        plt.savefig(f'sar_spa_qq_wetdry_{today}.png')
else: 
    if saveplots == True:
        plt.savefig(f'sar_spa_qq_wetdry_norm_{today}.pdf')
        plt.savefig(f'sar_spa_qq_wetdry_norm_{today}.png')
    

#%%

# QQ plots of the whole season

norm = False

sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = np.array(results['soil_rootzone_moisture'])
spa_wliq_top = np.array(results['soil_moisture_top'])

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 

#spa dates to match sar dates
date_in_spa = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]

sar_flat = sar_wliq[:,:,:].flatten()
spa_flat = spa_wliq[:,:,:].flatten()
spa_top_flat = spa_wliq_top[:,:,:].flatten()

flat_pd = pd.DataFrame()
flat_pd['sar'] = sar_flat
flat_pd['spa'] = spa_flat
flat_pd['spa_top'] = spa_top_flat
flat_pd = flat_pd.loc[np.isfinite(flat_pd['sar']) & np.isfinite(flat_pd['spa']) & np.isfinite(flat_pd['spa_top'])]


# Plotting
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4));
ax1 = axs[0]
ax2 = axs[1]

if norm == False:
    x1 = sns.regplot(ax=ax1, x=flat_pd['sar'], y=flat_pd['spa'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
    ax1.set(ylim=(0, 1))
    ax1.set(xlim=(0, 1))
    #ax1.set_title('Dry day')

    x2 = sns.regplot(ax=ax2, x=flat_pd['sar'], y=flat_pd['spa_top'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
    ax2.set(ylim=(0, 1))
    ax2.set(xlim=(0, 1))
    #ax2.set_title('Wet day')
    
    fig.suptitle('SpaFHy v2D')
    
    if saveplots == True:
        plt.savefig(f'sar_spa_qq_{today}.pdf')
        plt.savefig(f'sar_spa_qq_{today}.png')
    
else:
    x1 = sns.regplot(ax=ax1, x=flat_pd['sar']/(flat_pd['sar'].mean()), y=flat_pd['spa']/(flat_pd['spa'].mean()), scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
    ax1.set(ylim=(0, 2.5))
    ax1.set(xlim=(0, 2.5))
    #ax1.set_title('Dry day')

    x2 = sns.regplot(ax=ax2, x=flat_pd['sar']/(flat_pd['sar'].mean()), y=flat_pd['spa_top']/(flat_pd['spa_top'].mean()), scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
    ax2.set(ylim=(0, 2.5))
    ax2.set(xlim=(0, 2.5))
    #ax2.set_title('Wet day')  
    
    fig.suptitle('SpaFHy v2D, norm by total means of each')

    
    if saveplots == True:
        plt.savefig(f'sar_spa_qq_norm_{today}.pdf')
        plt.savefig(f'sar_spa_qq_norm_{today}.png')


#%%

# QQ plots of the whole season for sitetypes 1 and 4


sar_wliq = sar['soilmoisture']*np.array(results['parameters_cmask'])/100
spa_wliq = np.array(results['soil_rootzone_moisture'])
spa_wliq_top = np.array(results['soil_moisture_top'])

dates_sar = sar['time'][:]
dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 

#spa dates to match sar dates
date_in_spa = []
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]

sar_wliq_site1 = sar_wliq[day_hi,np.where(sitetype == 1)[0],np.where(sitetype == 1)[1]]
spa_wliq_site1 = spa_wliq[day_hi,np.where(sitetype == 1)[0],np.where(sitetype == 1)[1]]
spa_wliq_top_site1 = spa_wliq_top[day_hi,np.where(sitetype == 1)[0],np.where(sitetype == 1)[1]]

sar_flat_1 = sar_wliq_site1.flatten()
spa_flat_1 = spa_wliq_site1.flatten()
spa_top_flat_1 = spa_wliq_top_site1.flatten()

sar_wliq_site4 = sar_wliq[:,np.where(sitetype == 4)[0],np.where(sitetype == 4)[1]]
spa_wliq_site4 = spa_wliq[:,np.where(sitetype == 4)[0],np.where(sitetype == 4)[1]]
spa_wliq_top_site4 = spa_wliq_top[:,np.where(sitetype == 4)[0],np.where(sitetype == 4)[1]]

sar_flat_4 = sar_wliq_site4.flatten()
spa_flat_4 = spa_wliq_site4.flatten()
spa_top_flat_4 = spa_wliq_top_site4.flatten()

flat_pd1 = pd.DataFrame()
flat_pd4 = pd.DataFrame()

flat_pd1['sar'] = sar_flat_1
flat_pd1['spa'] = spa_flat_1
flat_pd1['spa_top'] = spa_top_flat_1
flat_pd4['sar'] = sar_flat_4
flat_pd4['spa'] = spa_flat_4
flat_pd4['spa_top'] = spa_top_flat_4
flat_pd1 = flat_pd1.loc[np.isfinite(flat_pd1['sar']) & np.isfinite(flat_pd1['spa']) & np.isfinite(flat_pd1['spa_top'])]
flat_pd4 = flat_pd4.loc[np.isfinite(flat_pd4['sar']) & np.isfinite(flat_pd4['spa']) & np.isfinite(flat_pd4['spa_top'])]


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]

x1 = sns.regplot(ax=ax1, x=flat_pd4['sar'], y=flat_pd4['spa'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
ax1.set(ylim=(0, 1))
ax1.set(xlim=(0, 1))
ax1.set_title('Mire')

x2 = sns.regplot(ax=ax2, x=flat_pd4['sar'], y=flat_pd4['spa_top'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
ax2.set(ylim=(0, 1))
ax2.set(xlim=(0, 1))
ax2.set_title('Mire')

x3 = sns.regplot(ax=ax3, x=flat_pd1['sar'], y=flat_pd1['spa'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
ax3.set(ylim=(0, 1))
ax3.set(xlim=(0, 1))
ax3.set_title('Mineral')

x4 = sns.regplot(ax=ax4, x=flat_pd1['sar'], y=flat_pd1['spa_top'], scatter_kws={'s':10, 'alpha':0.005}, line_kws={"color": "red"})
ax4.set(ylim=(0, 1))
ax4.set(xlim=(0, 1))
ax4.set_title('Mineral')

fig.suptitle('SpaFHy v2D')
    

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]

x1 = sns.scatterplot(ax=ax1, x=flat_pd4['sar'], y=flat_pd4['spa'], alpha=0.05)
ax1.set(ylim=(0, 1))
ax1.set(xlim=(0, 1))
ax1.set_title('Mire')

x2 = sns.scatterplot(ax=ax2, x=flat_pd4['sar'], y=flat_pd4['spa_top'], alpha=0.05)
ax2.set(ylim=(0, 1))
ax2.set(xlim=(0, 1))
ax2.set_title('Mire')

x3 = sns.scatterplot(ax=ax3, x=flat_pd1['sar'], y=flat_pd1['spa'], alpha=0.05)
ax3.set(ylim=(0, 1))
ax3.set(xlim=(0, 1))
ax3.set_title('Mineral')

x4 = sns.scatterplot(ax=ax4, x=flat_pd1['sar'], y=flat_pd1['spa_top'], alpha=0.05)
ax4.set(ylim=(0, 1))
ax4.set(xlim=(0, 1))
ax4.set_title('Mineral')

fig.suptitle('SpaFHy v2D')

   
#%%

# QQ plots of spa and sar

# Plotting
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12,4));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]

x1 = sns.scatterplot(ax=ax1, x=flat_pd['sar'], y=flat_pd['sar'], alpha=0.003)
ax1.set(ylim=(0, 1))
ax1.set(xlim=(0, 1))
#ax1.set_title('Dry day')

x2 = sns.scatterplot(ax=ax2, x=flat_pd['spa'], y=flat_pd['spa'], alpha=0.003)
ax2.set(ylim=(0, 1))
ax2.set(xlim=(0, 1))
#ax2.set_title('Wet day')

x3 = sns.scatterplot(ax=ax3, x=flat_pd['spa_top'], y=flat_pd['spa_top'], alpha=0.003)
ax3.set(ylim=(0, 1))
ax3.set(xlim=(0, 1))
#ax2.set_title('Wet day')
    
fig.suptitle('SpaFHy v2D')
    
if saveplots == True:
        #plt.savefig(f'qq_spa_spatop_sar_{today}.pdf')
        plt.savefig(f'qq_spa_spatop_sar_{today}.png')


#%%
# drainage plots

#catchment area in km2
area = len(np.array(results['parameters_cmask']).flatten()[np.isfinite(np.array(results['parameters_cmask'])).flatten()])*16*16*1e-6

# 

Qs = pd.DataFrame()
Qs['Qdrain'] = np.nanmean(results['soil_drainage']*cmask, axis=(1,2))
Qs['Qsurf'] = np.nanmean(results['soil_surface_runoff']*cmask, axis=(1,2))
Qs.index = dates_spa

Qm = pd.read_csv('C:\SpaFHy_v1_Pallas_2D\calibration\Runoffs1d_SVEcatchments_mmd.csv', 
                    sep=';', parse_dates = ['pvm'])
Qm.index = Qm['pvm']
Qm = Qm['1_Lompolojanganoja']
Qs['Qm'] = Qm


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,6));
ax1 = axs[0]
ax2 = axs[1]

ax1.plot(Qs['Qdrain'])
ax1.plot(Qs['Qsurf'])
ax1.plot(Qs['Qdrain']+Qs['Qsurf'])
ax1.set_ylim(-1,8)
ax1.title.set_text('SpaFHy2D')
ax1.legend(['Qdrain', 'Qsurf', 'Qdrain + Qsurf'], ncol = 3)
ax2.plot(Qs['Qm'])
ax2.set_ylim(-1,8)
ax2.title.set_text('Measured mm/d')


#plt.ylabel('Drainage mm/d')

if saveplots == True:
        #plt.savefig(f'spafhy2d_drainage_{today}.pdf')
        plt.savefig(f'spafhy2d_drainage_{today}.png')
 
    
#%%

# plot slope vs. sar
from matplotlib.colors import LogNorm

slope, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_slope.txt')
flowac, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_flow_accum_d8_area.txt')

#scmask = np.where(soilclass != 2)
sc2 = np.where(soilclass == 2)
st1 = np.where(sitetype == 1)

fdf = pd.DataFrame()
fdf['sar'] = sar_wliq[day_hi,sc2[0],sc2[1]].flatten()
fdf['slope'] = slope[sc2[0],sc2[1]].flatten()
fdf['flowac'] = flowac[sc2[0],sc2[1]].flatten()
fdf = fdf.loc[np.isfinite(fdf['sar']) & np.isfinite(fdf['slope'])]

# Plotting
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,6));
ax1 = axs[0]
ax2 = axs[1]

sns.regplot(ax=ax1, x=fdf['sar'], y=fdf['slope'], scatter_kws={'s':15, 'alpha':0.05}, line_kws={"color": "red"})
#ax1.set(ylim=(0, 10))
#ax1.set(xlim=(0.05,0.6))

sns.regplot(ax=ax2, x=fdf['sar'], y=fdf['flowac'], scatter_kws={'s':15, 'alpha':0.05}, line_kws={"color": "red"})
ax2.set_yscale('log') 
#ax1.set(xlim=(0.05,0.6))

