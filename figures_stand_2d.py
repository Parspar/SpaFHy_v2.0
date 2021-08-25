# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 10:22:54 2021

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
#import pickle
#import seaborn as sns
from iotools import read_AsciiGrid


# reading the stand results
outputfile_stand = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_202108251256.nc'
results_stand = read_results(outputfile_stand)

# reading the stand results
outputfile_2d = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_2d.nc'
results_2d = read_results(outputfile_2d)

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_scd.nc'
#sar = Dataset(sar_path, 'r')

# water table at lompolonjänkä


# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[1])-2, int(kenttarova_loc[0])-2])
l_loc = [60, 60]

dates_spa = []
for d in range(len(results_stand['date'])):
    dates_spa.append(pd.to_datetime(str(results_stand['date'][d])[36:46]))

# forcing file
folder = r'C:\SpaFHy_v1_Pallas_2D\testcase_input\forcing'
ffile = 'Kenttarova_forcing.csv'
fp = os.path.join(folder, ffile)
forc = pd.read_csv(fp, sep=';', date_parser=['time'])
forc['time'] = pd.to_datetime(forc['time'])
forc.index = forc['time']
forc = forc[forc['time'].isin(dates_spa)]

ix_no_p = np.where(forc['rainfall'] == 0)[0]
ix_p = np.where(forc['rainfall'] > 0)[0]
#ix_no_p_d = forc.index[np.where(forc['rainfall'] == 0)[0]]

# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

saveplots = True
today = date.today()

cmask = results_stand['parameters_cmask']
sitetype = results_stand['parameters_sitetype']
soilclass = results_stand['parameters_soilclass']


# water balance check
# 2D
results_2d['soil_water_storage'] = results_2d['soil_water_storage'] - results_2d['soil_water_storage'][0,:,:] 
results_2d['bucket_water_storage'] =  results_2d['bucket_water_storage'] - results_2d['bucket_water_storage'][0,:,:]
results_2d['soil_water_storage'] = results_2d['soil_water_storage'] + results_2d['bucket_water_storage']
results_2d['soil_netflow_to_ditch'] = results_2d['soil_netflow_to_ditch'] * results_2d['parameters_cmask']

# stand
results_stand['bucket_water_storage'] = results_stand['bucket_water_storage'] - results_stand['bucket_water_storage'][0,:,:]


# results from spafhy topmodel to comparable form

ncf_file = r'C:\SpaFHy_v1_Pallas\results\C3.nc'
dat = Dataset(ncf_file, 'r')

bu_catch = dat['bu']
cpy_catch = dat['cpy']
top_catch = dat['top']

ix = np.where(np.array(bu_catch['Wliq'][:,100,100]) <= 1.0)[0][0]
Wliq = bu_catch['Wliq'][ix:,:,:]
Wliq_top = bu_catch['Wliq_top'][ix:,:,:]
Transpi = cpy_catch['Transpi'][ix:,:,:]
Efloor = cpy_catch['Efloor'][ix:,:,:]
Evap = cpy_catch['Evap'][ix:,:,:]



#%%

plt.figure(figsize=(20,8))
ax=plt.subplot(2,1,1)
plt.title('2D')
#plt.subplot(2,1,1,sharex=ax)
plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j']), 'b', label='water_storage')

plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results_2d['canopy_evaporation'].mean(['i','j'])), 'darkgreen', label='+canopy_evaporation')

plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results_2d['canopy_evaporation'].mean(['i','j'])
                   + results_2d['canopy_transpiration'].mean(['i','j'])),'lime', alpha=0.3, label='+canopy_transpiration')

plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results_2d['canopy_evaporation'].mean(['i','j'])
                   + results_2d['canopy_transpiration'].mean(['i','j'])
                   + results_2d['bucket_evaporation'].mean(['i','j'])),'turquoise', label='+bucket_evaporation')

plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results_2d['canopy_evaporation'].mean(['i','j'])
                   + results_2d['canopy_transpiration'].mean(['i','j'])
                   + results_2d['bucket_evaporation'].mean(['i','j'])
                   + results_2d['bucket_surface_runoff'].mean(['i','j'])),  'gold', label='+bucket_surface_runoff')

plt.plot(results_2d['date'],results_2d['soil_water_storage'].mean(['i','j'])+
         np.cumsum(results_2d['canopy_evaporation'].mean(['i','j'])
                   + results_2d['canopy_transpiration'].mean(['i','j'])
                   + results_2d['bucket_evaporation'].mean(['i','j'])
                   + results_2d['bucket_surface_runoff'].mean(['i','j'])
                   + results_2d['soil_netflow_to_ditch'].mean(['i','j'])), 'peru', label='+soil_netflow_to_ditch')

plt.plot(results_2d['date'],np.cumsum(results_2d['forcing_precipitation']),'k.', alpha=0.5, label='forcing_precipitation')
plt.legend()


plt.subplot(2,1,2,sharex=ax)
plt.title('stand')
plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j']),'b', label='water_storage')

plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results_stand['canopy_evaporation'].mean(['i','j'])),'darkgreen', label='+canopy_evaporation')

plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results_stand['canopy_evaporation'].mean(['i','j'])
                   + results_stand['canopy_transpiration'].mean(['i','j'])),'lime', alpha=0.3, label='+canopy_transpiration')

plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results_stand['canopy_evaporation'].mean(['i','j'])
                   + results_stand['canopy_transpiration'].mean(['i','j'])
                   + results_stand['bucket_evaporation'].mean(['i','j'])),'turquoise', label='+bucket_evaporation')

plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j'])+
         np.cumsum(results_stand['canopy_evaporation'].mean(['i','j'])
                   + results_stand['canopy_transpiration'].mean(['i','j'])
                   + results_stand['bucket_evaporation'].mean(['i','j'])
                   + results_stand['bucket_drainage'].mean(['i','j'])),  'peru',  label='+bucket_drainage')

#plt.plot(results_stand['date'],results_stand['bucket_water_storage'].mean(['i','j'])+
#         np.cumsum(results_stand['canopy_evaporation'].mean(['i','j'])
#                   + results_stand['canopy_transpiration'].mean(['i','j'])
#                   + results_stand['bucket_evaporation'].mean(['i','j'])
#                   + results_stand['bucket_drainage'].mean(['i','j'])
#                   + results_stand['bucket_surface_runoff'].mean(['i', 'j'])),  'y',alpha=0.3, label='+bucket_surface_runoff')

plt.plot(results_stand['date'],np.cumsum(results_stand['forcing_precipitation']),'k.', alpha=0.5, label='forcing_precipitation')
plt.legend()

#%%

# soil moist plots

# spafhy
wet_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmax()
dry_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmin()

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(24,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[0][3]
ax5 = axs[1][0]
ax6 = axs[1][1]
ax7 = axs[1][2]
ax8 = axs[1][3]

# 2D top moist
im1 = ax1.imshow(results_2d['bucket_moisture_top'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im1, ax=ax1)
ax1.title.set_text('2D top wet')

im2 = ax2.imshow(results_2d['bucket_moisture_top'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im2, ax=ax2)
ax2.title.set_text('2D top dry')

# 2D root moist
im3 = ax3.imshow(results_2d['bucket_moisture_root'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im3, ax=ax3)
ax3.title.set_text('2D root wet')

im4 = ax4.imshow(results_2d['bucket_moisture_root'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('2D root dry')

# stand top moist
im5 = ax5.imshow(results_stand['bucket_moisture_top'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im5, ax=ax5)
ax5.title.set_text('stand top wet')

im6 = ax6.imshow(results_stand['bucket_moisture_top'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im6, ax=ax6)
ax6.title.set_text('stand top dry')

# stand root moist
im7 = ax7.imshow(results_stand['bucket_moisture_root'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im7, ax=ax7)
ax7.title.set_text('stand root wet')

im8 = ax8.imshow(results_stand['bucket_moisture_root'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im8, ax=ax8)
ax8.title.set_text('stand root dry')

'''
ax1.axis("off")
ax2.axis("off")
ax3.axis("off")
ax4.axis("off")
ax5.axis("off")
ax6.axis("off")
ax7.axis("off")
ax8.axis("off")
'''

#%%

# soil moist plots without topsoil

# soil moist plots

# spafhy
wet_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmax()
dry_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmin()

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]

# 2D root moist
im1 = ax1.imshow(results_2d['bucket_moisture_root'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im3, ax=ax3)
ax1.title.set_text('2D root wet')

im2 = ax2.imshow(results_2d['bucket_moisture_root'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.title.set_text('2D root dry')

# stand root moist
im3 = ax3.imshow(results_stand['bucket_moisture_root'][wet_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im7, ax=ax7)
ax3.title.set_text('stand root wet')

im4 = ax4.imshow(results_stand['bucket_moisture_root'][dry_day,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im4, ax=ax4)
ax4.title.set_text('stand root dry')


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

#soilm = soilscout.merge(ecmoist)

soilm = pd.concat([ecmoist, soilscout]).sort_values('time').reset_index(drop=True)


l_loc = [60, 60]
spa_wliq_2d_root = results_2d['bucket_moisture_root']
spa_wliq_2d_top = results_2d['bucket_moisture_top']

spa_wliq_st_root = results_stand['bucket_moisture_root']
spa_wliq_st_top = results_stand['bucket_moisture_top']

#dates_sar = sar['time'][:]
#dates_sar = pd.to_datetime(dates_sar, format='%Y%m%d') 

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
# 2D
spa_wliq_df['spa_k_2d_root'] = spa_wliq_2d_root[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_2d_root'] = spa_wliq_2d_root[:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_2d_top'] = spa_wliq_2d_top[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_2d_top'] = spa_wliq_2d_top[:,l_loc[0],l_loc[1]]

# stand
spa_wliq_df['spa_k_st_root'] = spa_wliq_st_root[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_st_root'] = spa_wliq_st_root[:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_st_top'] = spa_wliq_st_top[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_st_top'] = spa_wliq_st_top[:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_ca_root'] = Wliq[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_root'] = Wliq[:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_ca_top'] = Wliq[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_top'] = Wliq[:,l_loc[0],l_loc[1]]

spa_wliq_df['time'] = dates_spa
#spa_wliq_df.index = dates_spa

#soilm = pd.concat([soilm, spa_wliq_df]).sort_values('time').reset_index(drop=True)
soilm = spa_wliq_df.merge(soilm)
soilm.index = soilm['time']
#soilm = soilm[['s3', 's5', 's18', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B', 'spa_k', 'spa_l', 'spatop_k', 'spatop_l']]

#soilm = soilm.loc[(soilm.index > '2018-04-01') & (soilm.index < '2019-12-01')]


poi = 1085
# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,4));
ax1 = axs[0]
ax2 = axs[1]

#fig.suptitle('Volumetric water content', fontsize=15)

im1 = ax1.plot(soilm['spa_k_2d_root'], alpha=0.8)
ax1.plot(soilm['spa_k_st_root'], alpha=0.8)
ax1.plot(soilm['spa_k_ca_root'], alpha=0.8)
#ax1.plot(soilm.index[poi], 0.05, marker='o', mec='k', mfc='r', alpha=0.5, ms=8.0)
#ax1.axvline(soilm.index[poi], 0, 0.6, label='pyplot vertical line')
ax1.axvline(soilm.index[poi], ymin=0, ymax=1, color='r', alpha=0.4)

#ax1.plot(soilm['s3'], 'r', alpha=0.4)
#ax1.plot(soilm['s5'], alpha=0.4)
#ax1.plot(soilm['s18'], 'r', alpha=0.9)
#ax1.plot(soilm['SH-5A'], 'g', alpha=0.6)
#ax1.plot(soilm['SH-5B'], 'g', alpha=0.2)
#ax1.plot(soilm['SH-20A'], 'g', alpha=0.8)
#ax1.plot(soilm['SH-20B'], 'g', alpha=0.95)
ax1.title.set_text('Mineral')
ax1.legend(['2D root', 'stand root', 'catch root', 'spatial plot'], ncol=2)# 's3 = -0.05', 's18 = -0.3', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B'], ncol = 8)
y = ax1.set_ylabel(r"${\Theta}$")
y.set_rotation(0)
ax1.set_ylim(0,0.6)

im2 = ax2.plot(soilm['spa_l_2d_root'], alpha=0.8)
ax2.plot(soilm['spa_l_st_root'], alpha=0.8)
ax2.plot(soilm['spa_l_ca_root'], alpha=0.8)
#ax2.plot(soilm.index[poi], 0.45, marker='o', mec='k', mfc='g', alpha=0.5, ms=8.0)
ax2.axvline(soilm.index[poi], ymin=0, ymax=1, color='k', alpha=0.4)
ax2.title.set_text('Mire')
ax2.set_ylim(0.35,1.0)
y = ax2.set_ylabel(r"${\Theta}$")
y.set_rotation(0)
ax2.legend(['2D root', 'stand root', 'catch root', 'spatial plot'], ncol=2)# 's3 = -0.05', 's18 = -0.3', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B'], ncol = 8)

#g.suptitle('SpaFHy stand vs. 2D')

#%% 

# where 2D is interesting? 
# kenttärovalla esim 2018 kesäkuussa

# soil moist plots

# spafhy
wet_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmax()
dry_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmin()

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')

# Plotting
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(21,6));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]


# 2D root moist
im1 = ax1.imshow(results_2d['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im1, ax=ax1)
ax1.scatter(k_loc[0],k_loc[1],color='r', alpha=0.4)
ax1.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax1.title.set_text(f'2D root {dates_spa[poi]}')

# stand moist
im2 = ax2.imshow(results_stand['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax2.scatter(k_loc[0],k_loc[1],color='r', alpha=0.4)
ax2.title.set_text(f'stand root {dates_spa[poi]}')

# catchment moist
im3 = ax3.imshow(Wliq[poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im3, ax=ax3)
ax3.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax3.scatter(k_loc[0],k_loc[1],color='r', alpha=0.4)
ax3.title.set_text(f'catch root {dates_spa[poi]}')



#%%

# et sim. obs
# ET obs
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
file = 'ec_et.csv'
fp = os.path.join(folder, file)
ec_full = pd.read_csv(fp, sep=';', date_parser=['time'])
ec_full['time'] = pd.to_datetime(ec_full['time'])
#ec_kr.index = ec_kr['time']
ec = ec_full[ec_full['time'].isin(dates_spa)]
ec.index = ec['time']
ec_full.index = ec_full['time']
ec['ET-1-KR'].iloc[ec['time'] < '2017-01-01'] = np.nan
ec['ET-1-LV'].iloc[ec['time'] < '2017-01-01'] = np.nan
ec['ET-1-KR'].iloc[ix_p] = np.nan
ec['ET-1-LV'].iloc[ix_p] = np.nan

#dry et
results_2d['dry_et'] = results_2d['canopy_evaporation'] + results_2d['bucket_evaporation'] + results_2d['canopy_transpiration']
results_2d['dry_et'][ix_p,:,:] = np.nan
results_stand['dry_et'] = results_stand['canopy_evaporation'] + results_stand['bucket_evaporation'] + results_stand['canopy_transpiration']
results_stand['dry_et'][ix_p,:,:] = np.nan
catch_dry_et = Transpi + Efloor + Evap
catch_dry_et[ix_p,:,:] = np.nan


# Plotting
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(20,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[2][0]
ax6 = axs[2][1]

ax1.set_title('Kenttärova')
ax1.plot(dates_spa, results_stand['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.7)
ax1.plot(ec['ET-1-KR'], 'k.', alpha=0.4)
ax1.legend(['stand dry et', 'ec et'])

ax2.set_title('Lompolo')
ax2.plot(dates_spa, results_stand['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.7)
ax2.plot(ec['ET-1-LV'], 'k.', alpha=0.4)
ax2.legend(['2d dry et', 'stand dry et', 'catch dry et', 'ec et'])

ax3.set_title('Kenttärova')
ax3.plot(dates_spa, results_2d['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.7)
ax3.plot(ec['ET-1-KR'], 'k.', alpha=0.4)
ax3.legend(['2d dry et', 'ec et'])

ax4.set_title('Lompolo')
ax4.plot(dates_spa, results_2d['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.7)
ax4.plot(ec['ET-1-LV'], 'k.', alpha=0.4)
ax4.legend(['2d dry et', 'ec et'])

ax5.set_title('Kenttärova')
ax5.plot(dates_spa, catch_dry_et[:,k_loc[0],k_loc[1]], alpha=0.7)
ax5.plot(ec['ET-1-KR'], 'k.', alpha=0.4)
ax5.legend(['catch dry et', 'ec et'])

ax6.set_title('Lompolo')
ax6.plot(dates_spa, catch_dry_et[:,l_loc[0],l_loc[1]], alpha=0.7)
ax6.plot(ec['ET-1-LV'], 'k.', alpha=0.4)
ax6.legend(['catch dry et', 'ec et'])




#%%
import seaborn as sns

# scatter et sim vs. obs
ET_flat = pd.DataFrame()
ET_flat['stand_kr_et'] = results_stand['dry_et'][:,k_loc[0],k_loc[0]]
ET_flat['stand_lv_et'] = results_stand['dry_et'][:,l_loc[0],l_loc[0]]
ET_flat['2D_kr_et'] = results_2d['dry_et'][:,k_loc[0],k_loc[0]]
ET_flat['2D_lv_et'] = results_2d['dry_et'][:,l_loc[0],l_loc[0]] 
ET_flat['catch_kr_et'] = catch_dry_et[:,k_loc[0],k_loc[0]] 
ET_flat['catch_lv_et'] = catch_dry_et[:,l_loc[0],l_loc[0]]
ET_flat['ET-1-KR'] = ec['ET-1-KR'].reset_index(drop=True)
ET_flat['ET-1-LV'] = ec['ET-1-LV'].reset_index(drop=True)
ET_flat = ET_flat.dropna()
ET_flat = ET_flat[(ET_flat > 0).all(1)]

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8,12));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[1][0]
ax4 = axs[1][1]
ax5 = axs[2][0]
ax6 = axs[2][1]

x1 = sns.regplot(ax=ax1, x=ET_flat['stand_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
x1.plot([1, 1], [4, 4])
ax1.set(ylim=(0, 5))
ax1.set(xlim=(0, 5))
x2 = sns.regplot(ax=ax2, x=ET_flat['stand_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
ax2.set(ylim=(0, 5))
ax2.set(xlim=(0, 5))

x3 = sns.regplot(ax=ax3, x=ET_flat['2D_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
x3.plot([1, 1], [4, 4])
ax3.set(ylim=(0, 5))
ax3.set(xlim=(0, 5))
x4 = sns.regplot(ax=ax4, x=ET_flat['2D_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
ax4.set(ylim=(0, 5))
ax4.set(xlim=(0, 5))

x5 = sns.regplot(ax=ax5, x=ET_flat['catch_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
x5.plot([1, 1], [4, 4])
ax5.set(ylim=(0, 5))
ax5.set(xlim=(0, 5))
x6 = sns.regplot(ax=ax6, x=ET_flat['catch_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
ax6.set(ylim=(0, 5))
ax6.set(xlim=(0, 5))

#%%

# self.LAI comparison