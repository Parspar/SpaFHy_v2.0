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


# reading the results
outputfile = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input.nc'
results = read_results(outputfile)

dates_spa = []
for d in range(len(results['date'])):
    dates_spa.append(pd.to_datetime(str(results['date'][d])[36:46]))


# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

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

plt_ind = 30

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


norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(results['soil_drainage'][plt_ind,:,:]), 
                            vmax=np.nanmax(results['soil_drainage'][plt_ind,:,:]), vcenter=0)
im3 = ax3.imshow(results['soil_drainage'][-1,:,:], cmap=plt.cm.RdBu, norm=norm)
fig.colorbar(im3, ax=ax3)
ax3.title.set_text(f'soil drainage dates[{plt_ind}]')


im4 = ax4.imshow(results['parameters_elevation'], cmap='bwr')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('DEM')

im5 = ax5.imshow(results['parameters_soilclass'], cmap='bwr')
fig.colorbar(im5, ax=ax5)
ax5.title.set_text('soil class')


im6 = ax6.imshow(results['parameters_ditches'], cmap='binary')
fig.colorbar(im6, ax=ax6)
ax6.title.set_text('ditches and ponds')


im7 = ax7.imshow(results['soil_moisture_top'][plt_ind,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im7, ax=ax7)
ax7.title.set_text(f'soil moisture top dates[{plt_ind}]')

im8 = ax8.imshow(results['soil_rootzone_moisture'][plt_ind,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im8, ax=ax8)
ax8.title.set_text(f'soil rootzone_moisture dates[{plt_ind}]')

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

im2 = ax2.plot(results['soil_moisture_top'][:,175,60])
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


#%%

# soil moist plots

# spafhy
wet_day = np.nansum(results['soil_rootzone_moisture'], axis=(1,2)).argmax()
dry_day = np.nansum(results['soil_rootzone_moisture'], axis=(1,2)).argmin()

sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
sar = Dataset(sar_path, 'r')


# Plotting
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]

im1 = ax1.imshow(results['soil_moisture_top'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im1, ax=ax1)
ax1.title.set_text('soil moisture top wet')

im2 = ax2.imshow(results['soil_moisture_top'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('soil moisture top dry')

im3 = ax3.imshow(results['soil_rootzone_moisture'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im3, ax=ax3)
ax3.title.set_text('soil moisture root wet')

im4 = ax4.imshow(results['soil_rootzone_moisture'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('soil moisture root dry')

#%%

# sar soil moisture plots
import pandas as pd

# define a big catchment mask
from iotools import read_AsciiGrid

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
kenttarova_loc = list([int(kenttarova_loc[0]), int(kenttarova_loc[1])])

# reading sar data
sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_mean4.nc'
sar = Dataset(sar_path, 'r')

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


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_{today}.png')
    
    
#%%

# normalized plots

# mean of each pixel
spamean = np.nanmean(spa_wliq, axis=0)
spatopmean = np.nanmean(spa_wliq_top, axis=0)
sarmean = np.nanmean(sar_wliq, axis=0)

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

im1 = ax1.imshow(sar_wliq[day_hi,:,:]/sarmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax1.title.set_text('SAR')
#ax1.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im2 = ax2.imshow(spa_wliq[day_hi, :,:]/spamean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax2.title.set_text('SPAFHY rootzone')
ax2.text(10, -15, f'Wet day : {hi_date}', fontsize=15)
#ax2.plot(kenttarova_loc[0], kenttarova_loc[1], marker='o', mec='b', mfc='k', alpha=0.8, ms=6.0)

im3 = ax3.imshow(sar_wliq[day_low, :,:]/sarmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax3.title.set_text('SAR')

im4 = ax4.imshow(spa_wliq[day_low, :,:]/spamean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax4.title.set_text('SPAFHY rootzone')
ax4.text(10, -15, f'Dry day : {low_date}', fontsize=15)


im5 = ax5.imshow(spa_wliq_top[day_hi, :, :]/spatopmean, vmin=0.0, vmax=2.0, cmap='coolwarm_r', aspect='equal')
ax5.title.set_text('SPAFHY topsoil')

im6 = ax6.imshow(spa_wliq_top[day_low, :,:]/spatopmean, cmap='coolwarm_r', vmin=0.0, vmax=2.0, aspect='equal')
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


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.8, 0.2, 0.02, 0.6])
#bar2 = fig.colorbar(im6, cax=cbar_ax)

#ax5.plot([1.1, 1.1], [1.1, -1.2], color='black', lw=1, transform=ax2.transAxes, clip_on=False)


if saveplots == True:
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_norm_{today}.pdf')
    plt.savefig(f'SAR_vs.SPAFHY_soilmoist_norm_{today}.png')


#%%

# point examples from mineral and openmire

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
for i in range(len(dates_sar)):
    ix = np.where(dates_spa == dates_sar[i])[0][0]
    date_in_spa.append(ix)
spa_wliq = spa_wliq[date_in_spa,:,:]
spa_wliq_top = spa_wliq_top[date_in_spa,:,:]

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,8));
ax1 = axs[0]
ax2 = axs[1]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.plot(sar_wliq[:,k_loc[0],k_loc[1]])
ax1.plot(spa_wliq[:,k_loc[0],k_loc[1]])
ax1.plot(spa_wliq_top[:,k_loc[0],k_loc[1]])
ax1.title.set_text('Mineral')
ax1.legend(['SAR', 'SpaFHy rootzone', 'SpaFHy top'], ncol = 3)


im2 = ax2.plot(sar_wliq[:,l_loc[0],l_loc[1]])
ax2.plot(spa_wliq[:,l_loc[0],l_loc[1]])
ax2.plot(spa_wliq_top[:,l_loc[0],l_loc[1]])
ax2.title.set_text('Open mire')
ax2.legend(['SAR', 'SpaFHy rootzone', 'SpaFHy top'], ncol = 3)


#%%

# Q-Q plots

import numpy as np 
import pylab 
import scipy.stats as stats
import pandas as pd

from iotools import read_AsciiGrid

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

inb = 35
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
ax1.set(ylim=(0, 1))
ax1.set(xlim=(0, 1))
#ax1.set_title('Dry day')

x2 = sns.regplot(ax=ax2, x=flat_pd['sar_wet'], y=flat_pd['spa_wet'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
ax2.set(ylim=(0, 1))
ax2.set(xlim=(0, 1))
#ax2.set_title('Wet day')

x3 = sns.regplot(ax=ax3, x=flat_pd['sar_dry'], y=flat_pd['spa_top_dry'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
ax3.set(ylim=(0, 1))
ax3.set(xlim=(0, 1))
#ax3.set_title('Dry day')

x4 = sns.regplot(ax=ax4, x=flat_pd['sar_wet'], y=flat_pd['spa_top_wet'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
ax4.set(ylim=(0, 1))
ax4.set(xlim=(0, 1))
#ax4.set_title('Wet day')

x5 = sns.regplot(ax=ax5, x=flat_pd['sar_inb'], y=flat_pd['spa_inb'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
ax5.set(ylim=(0, 1))
ax5.set(xlim=(0, 1))
#ax4.set_title('Wet day')

x6 = sns.regplot(ax=ax6, x=flat_pd['sar_inb'], y=flat_pd['spa_top_inb'], scatter_kws={'s':10, 'alpha':0.1}, line_kws={"color": "red"})
ax6.set(ylim=(0, 1))
ax6.set(xlim=(0, 1))
#ax4.set_title('Wet day')

fig.suptitle('SpaFHy v2D')

if saveplots == True:
    plt.savefig(f'sar_spa_qq_wetdry_{today}.pdf')
    plt.savefig(f'sar_spa_qq_wetdry_{today}.png')