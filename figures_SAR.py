# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 12:34:49 2021

@author: janousu
"""


# SAR plots


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
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from sklearn.metrics import mean_absolute_error as mae


# reading the stand results
outputfile_stand = r'C:\SpaFHy_v1_Pallas_2D\results\testcase_input_2d_ditch03.nc'
results_stand = read_results(outputfile_stand)
cmask = results_stand['parameters_cmask']
soil = results_stand['parameters_soilclass']


sar_raw = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment.nc', 'r')

sar_ma3 = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma3.nc', 'r')

sar_ma5 = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5.nc', 'r')

sar_ma5_mean4 = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean4.nc', 'r')

sar_ma5_mean8 = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8.nc', 'r')

sar_ma5_mean8_sc = Dataset(r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_sc.nc', 'r')

sar_raw_np = np.zeros(shape=sar_raw['soilmoisture'].shape)
sar_ma3_np = np.zeros(shape=sar_raw['soilmoisture'].shape)
sar_ma5_np = np.zeros(shape=sar_raw['soilmoisture'].shape)
sar_ma5_mean4_np = np.zeros(shape=sar_raw['soilmoisture'].shape)
sar_ma5_mean8_np = np.zeros(shape=sar_raw['soilmoisture'].shape)
sar_ma5_mean8_sc_np = np.zeros(shape=sar_raw['soilmoisture'].shape)

for i in range(sar_raw['soilmoisture'].shape[0]):
    sar_raw_np[i] = sar_raw['soilmoisture'][i] * cmask
    sar_ma3_np[i] = sar_ma3['soilmoisture'][i] * cmask
    sar_ma5_np[i] = sar_ma5['soilmoisture'][i] * cmask
    sar_ma5_mean4_np[i] = sar_ma5_mean4['soilmoisture'][i] * cmask
    sar_ma5_mean8_np[i] = sar_ma5_mean8['soilmoisture'][i] * cmask
    sar_ma5_mean8_sc_np[i] = sar_ma5_mean8_sc['soilmoisture'][i] * cmask

sar_file = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_scd.nc'
sar = Dataset(sar_file, 'r')
#sar = np.array(sar['soilmoisture'])

sar_file2 = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment.nc'
sar2 = Dataset(sar_file2, 'r')
#sar2 = np.array(sar2['soilmoisture'])
# water table at lompolonjänkä

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[1]), int(kenttarova_loc[0])])
k_loc = [118, 136]
l_loc = [46, 54]

dates_sar = pd.to_datetime(sar['time'][:], format='%Y%m%d')


# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

saveplots = True
today = date.today()

# indexes for tighet plots
zx = np.arange(20, 171, 1)
zy = np.arange(20, 246, 1)

os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures\SAR_processing')


#%%

# ALKUTILANNE TEMPORAALINEN KOHINA

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,8));
ax1 = axs[0]
ax2 = axs[1]

ax1.plot(dates_sar, sar_raw['soilmoisture'][:,k_loc[0],k_loc[1]], label='raw')
ax1.set_title('Kenttärova')
ax2.plot(dates_sar, sar_raw['soilmoisture'][:,l_loc[0],l_loc[1]], label='raw')
ax2.set_title('Lompolonjänkkä')

ax1.grid(); ax2.grid()

plt.savefig('temporal_starting_point.pdf')

# %%

# ALKUTILANNE SPATIAALINEN KOHINA

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(20,8));
ax1 = axs[0][0]
ax2 = axs[1][0]
ax3 = axs[0][1]
ax4 = axs[1][1]

ax1.imshow(sar_raw['soilmoisture'][15], cmap='coolwarm_r', label='raw')
ax1.set_title(f'{dates_sar[15]}')
ax1.scatter(k_loc[1], k_loc[0], linewidth=3)
ax1.scatter(l_loc[1], l_loc[0], linewidth=3)
ax2.imshow(sar_raw['soilmoisture'][45], cmap='coolwarm_r', label='raw')
ax2.scatter(k_loc[1], k_loc[0], linewidth=3)
ax2.scatter(l_loc[1], l_loc[0], linewidth=3)
ax2.set_title(f'{dates_sar[45]}')

#ax1.axvline(y=k_loc[0],color='red')
ax1.axvline(x=l_loc[1],color='orange')
ax1.axvline(x=k_loc[1],color='blue')
ax2.axvline(x=l_loc[1],color='orange')
ax2.axvline(x=k_loc[1],color='blue')

ax3.plot(sar_raw['soilmoisture'][15, :, k_loc[1]], label='Kenttarova')
ax3.plot(sar_raw['soilmoisture'][15, :, l_loc[1]], label='Lompolo')
ax3.legend()
ax4.plot(sar_raw['soilmoisture'][45, :, k_loc[1]], label='Kenttarova')
ax4.plot(sar_raw['soilmoisture'][45, :, l_loc[1]], label='Lompolo')
ax4.legend()
plt.tight_layout()

ax3.grid(); ax4.grid()


plt.savefig('spatial_starting_point.pdf')


# %%
# timeseries of SAR at lompolonjankka and kenttarova

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,12));
ax1 = axs[0]
ax2 = axs[1]


ax1.plot(dates_sar, sar_raw['soilmoisture'][:,k_loc[0],k_loc[1]], label='raw')
ax1.plot(dates_sar, sar_ma3['soilmoisture'][:,k_loc[0],k_loc[1]], label='ma3')
ax1.plot(dates_sar, sar_ma5['soilmoisture'][:,k_loc[0],k_loc[1]], label='ma5')
#ax1.plot(dates_sar, sar_ma5_mean4['soilmoisture'][:,k_loc[0],k_loc[1]], label='ma5, mean4')
#ax1.plot(dates_sar, sar_ma5_mean8['soilmoisture'][:,k_loc[0],k_loc[1]],  label='ma5, mean8')
#ax1.plot(dates_sar, sar_ma5_mean8_scd['soilmoisture'][:,k_loc[0],k_loc[1]],  label='ma5, mean8, sc')
ax1.legend()
ax1.xaxis.set_tick_params(rotation=45)

ax2.plot(dates_sar, sar_raw['soilmoisture'][:,l_loc[0],l_loc[1]], label='raw')
ax2.plot(dates_sar, sar_ma3['soilmoisture'][:,l_loc[0],l_loc[1]], label='ma3')
ax2.plot(dates_sar, sar_ma5['soilmoisture'][:,l_loc[0],l_loc[1]], label='ma5')
#ax2.plot(dates_sar, sar_ma5_mean4['soilmoisture'][:,l_loc[0],l_loc[1]], label='ma5, mean4')
#ax2.plot(dates_sar, sar_ma5_mean8['soilmoisture'][:,l_loc[0],l_loc[1]],  label='ma5, mean8')
#ax2.plot(dates_sar, sar_ma5_mean8_scd['soilmoisture'][:,l_loc[0],l_loc[1]],  label='ma5, mean8, sc')
ax2.legend()
ax2.xaxis.set_tick_params(rotation=45)
ax1.set_title('Kenttärova')
ax2.set_title('Lompolonjänkkä')

ax1.grid(); ax2.grid()


plt.savefig('temporal_progress.pdf')


# %%

# SPATIAALINEN KOHINA LOPPU

fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(8,10));
ax1 = axs[0][0]
ax2 = axs[1][0]
ax3 = axs[2][0]
ax4 = axs[3][0]
ax5 = axs[0][1]
ax6 = axs[1][1]
ax7 = axs[2][1]
ax8 = axs[3][1]
ax9 = axs[4][0]
ax10 = axs[4][1]

ax1.imshow(sar_raw_np[15], cmap='coolwarm_r', vmin=0, vmax=100, label='raw')
ax1.set_title(f'{dates_sar[15]}')
ax1.scatter(k_loc[1], k_loc[0], linewidth=3)
ax1.scatter(l_loc[1], l_loc[0], linewidth=3)

ax2.imshow(sar_ma5_mean4_np[15], cmap='coolwarm_r', vmin=0, vmax=100, label='raw')
ax2.scatter(k_loc[1], k_loc[0], linewidth=3)
ax2.scatter(l_loc[1], l_loc[0], linewidth=3)

ax3.imshow(sar_ma5_mean8_np[15], cmap='coolwarm_r', vmin=0, vmax=100, label='raw')
ax3.scatter(k_loc[1], k_loc[0], linewidth=3)
ax3.scatter(l_loc[1], l_loc[0], linewidth=3)

ax4.imshow(sar_ma5_mean8_sc_np[15], cmap='coolwarm_r', vmin=0, vmax=100, label='raw')
ax4.scatter(k_loc[1], k_loc[0], linewidth=3)
ax4.scatter(l_loc[1], l_loc[0], linewidth=3)

#ax1.axvline(y=k_loc[0],color='red')
ax1.axvline(x=l_loc[1],color='orange')
ax1.axvline(x=k_loc[1],color='blue')
ax2.axvline(x=l_loc[1],color='orange')
ax2.axvline(x=k_loc[1],color='blue')

ax3.axvline(x=l_loc[1],color='orange')
ax3.axvline(x=k_loc[1],color='blue')
ax4.axvline(x=l_loc[1],color='orange')
ax4.axvline(x=k_loc[1],color='blue')

ax5.plot(sar_raw_np[15, :, k_loc[1]], label='Kenttarova')
ax5.plot(sar_raw_np[15, :, l_loc[1]], label='Lompolo')
ax5.legend()

ax6.plot(sar_ma5_mean4_np[15, :, k_loc[1]], label='Kenttarova')
ax6.plot(sar_ma5_mean4_np[15, :, l_loc[1]], label='Lompolo')

ax7.plot(sar_ma5_mean8_np[15, :, k_loc[1]], label='Kenttarova')
ax7.plot(sar_ma5_mean8_np[15, :, l_loc[1]], label='Lompolo')

ax8.plot(sar_ma5_mean8_sc_np[15, :, k_loc[1]], label='Kenttarova')
ax8.plot(sar_ma5_mean8_sc_np[15, :, l_loc[1]], label='Lompolo')

ax9.imshow(soil, label='raw')
ax9.scatter(k_loc[1], k_loc[0], linewidth=3)
ax9.scatter(l_loc[1], l_loc[0], linewidth=3)
ax9.axvline(x=l_loc[1],color='orange')
ax9.axvline(x=k_loc[1],color='blue')

ax10.plot(soil[:, k_loc[1]], label='Kenttarova')
ax10.plot(soil[:, l_loc[1]], label='Lompolo')

ax1.grid(); ax2.grid()
ax3.grid(); ax4.grid()
ax5.grid(); ax6.grid()
ax7.grid(); ax8.grid()
ax9.grid(); ax10.grid()

ax5.set_title('raw')
ax6.set_title('ma5_mean4')
ax7.set_title('ma5_mean8')
ax8.set_title('ma5_mean8_sc')
ax10.set_title('soilclass')

plt.tight_layout()

plt.savefig('spatial_progress.pdf')


