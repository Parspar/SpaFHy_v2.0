# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 17:14:12 2021

@author: janousu
"""

# debugging between 2D and original spafhy

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
outputfile_stand = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_202108251251.nc'
results_stand = read_results(outputfile_stand)

# reading the stand results
outputfile_2d = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_2d.nc'
results_2d = read_results(outputfile_2d)

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_scd.nc'
#sar = Dataset(sar_path, 'r')

# water table at lompolonjänkä

# forcing file
folder = r'C:\SpaFHy_v1_Pallas_2D\testcase_input\forcing'
ffile = 'Kenttarova_forcing.csv'
fp = os.path.join(folder, ffile)
forc = pd.read_csv(fp, sep=';', date_parser=['time'])
forc['time'] = pd.to_datetime(forc['time'])
forc.index = forc['time']

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[1])-2, int(kenttarova_loc[0])-2])
l_loc = [60, 60]


dates_spa = []
for d in range(len(results_stand['date'])):
    dates_spa.append(pd.to_datetime(str(results_stand['date'][d])[36:46]))


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

forc = forc[forc['time'].isin(dates_spa)]


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
# ET comparison kenttärova

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(18,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]
#ax4 = axs[3]

#fig.suptitle('Volumetric water content', fontsize=15)
ax1.set_title('Tr')
#ax1.plot(dates_spa, results_2d['canopy_transpiration'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax1.plot(dates_spa, results_stand['canopy_transpiration'][:,k_loc[0],k_loc[1]], 'r', alpha=0.6)
ax1.plot(dates_spa, Transpi[:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax1.legend(['2d', 'stand',' top'])

ax2.set_title('Ef')
#ax2.plot(dates_spa, results_2d['bucket_evaporation'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax2.plot(dates_spa, results_stand['bucket_evaporation'][:,k_loc[0],k_loc[1]], 'r', alpha=0.6)
ax2.plot(dates_spa, Efloor[:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
#ax2.legend(['2d', 'stand',' top'])

ax3.set_title('E')
#ax3.plot(dates_spa, results_2d['canopy_evaporation'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax3.plot(dates_spa, results_stand['canopy_evaporation'][:,k_loc[0],k_loc[1]], 'r', alpha=0.6)
ax3.plot(dates_spa, Evap[:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
#ax3.legend(['2d', 'stand',' top'])

#%%

# ET comparison suo
l_loc = [60, 60]

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(18,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]
#ax4 = axs[3]

#fig.suptitle('Volumetric water content', fontsize=15)
ax1.set_title('Tr')
#ax1.plot(dates_spa, results_2d['canopy_transpiration'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax1.plot(dates_spa, results_stand['canopy_transpiration'][:,l_loc[0],l_loc[1]], 'k.', alpha=0.6)
ax1.plot(dates_spa, Transpi[:,l_loc[0],l_loc[1]], 'b', alpha=0.6)
ax1.legend(['2d',' top'])

ax2.set_title('Ef')
#ax2.plot(dates_spa, results_2d['bucket_evaporation'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax2.plot(dates_spa, results_stand['bucket_evaporation'][:,l_loc[0],l_loc[1]], 'k.', alpha=0.6)
ax2.plot(dates_spa, Efloor[:,l_loc[0],l_loc[1]], 'b', alpha=0.6)
#ax2.legend(['2d', 'stand',' top'])

ax3.set_title('E')
#ax3.plot(dates_spa, results_2d['canopy_evaporation'][:,k_loc[0],k_loc[1]], 'b', alpha=0.6)
ax3.plot(dates_spa, results_stand['canopy_evaporation'][:,l_loc[0],l_loc[1]], 'k.', alpha=0.6)
ax3.plot(dates_spa, Evap[:,l_loc[0],l_loc[1]], 'b', alpha=0.6)
#ax3.legend(['2d', 'stand',' top'])

#%%

# self.LAI comparison

plt.plot(dates_spa, results_stand['canopy_fLAI'][:,l_loc[0],l_loc[1]], 'r', alpha=0.6)
plt.plot(dates_spa, cpy_catch['fLAI'][ix:,l_loc[0],l_loc[1]], 'k.')
plt.legend(['2d','top'])



#%%

plt.plot(dates_spa, results_2d['canopy_phenostate'][:,l_loc[0],l_loc[1]], 'k.')
plt.plot(dates_spa, cpy_catch['fPheno'][ix:,l_loc[0],l_loc[1]], 'r')
plt.plot(dates_spa, results_2d['canopy_transpiration'][:,l_loc[0],l_loc[1]], 'b', alpha=0.6)
plt.plot(dates_spa, Transpi[:,l_loc[0],l_loc[1]], 'k', alpha=0.3)

plt.legend(['2d_phenostate', 'catch_phenostate','2d_tr', 'catch_tr'])

#%%
# Plotting
fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(22,12));
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]
ax4 = axs[3]

tt = 600
tt1 = 800
ax1.plot(dates_spa[tt:tt1], cpy_catch['Evap'][ix+tt:ix+tt1,l_loc[0],l_loc[1]], 'k.')
ax1.plot(dates_spa[tt:tt1], results_stand['canopy_evaporation'][tt:tt1,l_loc[0],l_loc[1]], 'r')
#ax1.plot(dates_spa[tt:tt1], Wliq[tt:tt1,l_loc[0], l_loc[1]], 'g')
#ax1.plot(dates_spa[tt:tt1], results_stand['bucket_moisture_root'][tt:tt1,l_loc[0], l_loc[1]], 'r')
#ax2.plot(dates_spa[tt:tt1], forc['t_mean'][tt:tt1])
#plt.plot(dates_spa, top_catch['Sloc'][ix:,l_loc[0],l_loc[1]])
ax1.legend(['catch_evap', 'stand_evap'])

ax2.plot(dates_spa[tt:tt1], cpy_catch['SWE'][ix+tt:ix+tt1,l_loc[0],l_loc[1]]*1e-3, 'k.')
ax2.plot(dates_spa[tt:tt1], results_stand['canopy_snow_water_equivalent'][tt:tt1,l_loc[0],l_loc[1]]*1e-3, 'r')
#ax2.plot(dates_spa[tt:tt1], Wliq[tt:tt1,k_loc[0], k_loc[1]], 'g')
#ax2.plot(dates_spa[tt:tt1], results_stand['bucket_moisture_root'][tt:tt1,k_loc[0], k_loc[1]], 'r')
ax2.legend(['catch SWE', 'stand SWE'])

ax3.plot(dates_spa[tt:tt1], cpy_catch['WC'][ix+tt:ix+tt1,l_loc[0],l_loc[1]], 'k.')
ax3.plot(dates_spa[tt:tt1], results_stand['canopy_water_storage'][tt:tt1,l_loc[0],l_loc[1]], 'r')
ax3.legend(['catch WC', 'stand WC'])

ax4.plot(dates_spa[tt:tt1], cpy_catch['fLAI'][ix+tt:ix+tt1,l_loc[0],l_loc[1]], 'k.')
ax4.plot(dates_spa[tt:tt1], results_stand['canopy_fLAI'][tt:tt1,l_loc[0],l_loc[1]], 'r')
ax4.legend(['catch fLAI', 'stand fLAI'])

#%%

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,12));
ax1 = axs[0]
ax2 = axs[1]

tt = 600
ax1.plot(dates_spa[tt:], bu_catch['Drain'][ix+tt:,l_loc[0],l_loc[1]]*1e3, 'g')
ax1.plot(dates_spa[tt:], bu_catch['Retflow'][ix+tt:,l_loc[0],l_loc[1]]*1e3, 'k')
ax1.plot(dates_spa[tt:], results_stand['bucket_drainage'][tt:,l_loc[0],l_loc[1]], 'r')
ax1.plot(dates_spa[tt:], cpy_catch['SWE'][ix+tt:,l_loc[0],l_loc[1]]*1e-3, 'b')
#plt.plot(dates_spa, top_catch['Sloc'][ix:,l_loc[0],l_loc[1]])
#ax1.axvline(dates_spa[755], ymin=0, ymax=1, color='k.', alpha=0.4)
ax1.legend(['catch_drain', 'stand_drain'])

ax2.plot(dates_spa[tt:], Wliq[tt:,k_loc[0], k_loc[1]], 'g')
ax2.plot(dates_spa[tt:], results_stand['bucket_moisture_root'][tt:,k_loc[0], k_loc[1]], 'r')
ax2.legend(['catch_moist', 'stand_moist'])

#%%


