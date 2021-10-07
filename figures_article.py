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
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon



# reading the stand results
outputfile_stand = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_stand.nc'
results_stand = read_results(outputfile_stand)

# reading the stand results
outputfile_2d = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_2d.nc'
results_2d = read_results(outputfile_2d)

# reading the catch results
outputfile_catch = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_catch.nc'
results_catch = read_results(outputfile_catch)

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\16m_nc_spafhy_pallas\SAR_PALLAS_2019_mask2_16m_direct_catchment_ma5_mean8_scd.nc'
#sar = Dataset(sar_path, 'r')

# water table at lompolonjänkä


# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[1]), int(kenttarova_loc[0])])
k_loc = [118, 136]
l_loc = [46, 54]

dates_spa = []
for d in range(len(results_stand['date'])):
    dates_spa.append(pd.to_datetime(str(results_stand['date'][d])[36:46]))

# forcing file
folder = r'C:\SpaFHy_v1_Pallas_2D\testcase_input\forcing'
ffile = 'Kenttarova_forcing_new.csv'
fp = os.path.join(folder, ffile)
forc = pd.read_csv(fp, sep=';', date_parser=['time'])
forc['time'] = pd.to_datetime(forc['time'])
forc.index = forc['time']
forc = forc[forc['time'].isin(dates_spa)]
forc = forc[0:-105]
ix_no_p = np.where(results_stand['forcing_precipitation'] == 0)[0]
#ix_no_p = np.where(forc['rainfall'] == 0)[0]
ix_p = np.where(results_stand['forcing_precipitation'] > 0)[0]
#ix_p = np.where(forc['rainfall'] > 0)[0]
ix_pp = ix_p - 1
ix_pp = ix_pp[ix_pp >= 0]
#ix_no_p_d = forc.index[np.where(forc['rainfall'] == 0)[0]]

# specific discharge
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
ffile = 'Runoffs1d_SVEcatchments_mmd.csv'
fp = os.path.join(folder, ffile)
q = pd.read_csv(fp, sep=';', date_parser=['pvm'])
q = q.rename(columns={'pvm': 'time'})
q['time'] = pd.to_datetime(q['time'])
q.index = q['time']
q = q[q['time'].isin(dates_spa)]
q = q['1_Lompolojanganoja']

# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

saveplots = True
today = date.today()

cmask = results_stand['parameters_cmask']
sitetype = results_stand['parameters_sitetype']
soilclass = np.array(results_stand['parameters_soilclass'])
soilclass_copy = np.array(results_stand['parameters_soilclass'])

# indexes for tighet plots
zx = np.arange(20, 171, 1)
zy = np.arange(20, 246, 1)

# results from spafhy topmodel to comparable form
'''
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

'''

#%%

# defining first of july when there is no snowpack to be a starting date
dates_spa_str = []
for i in range(len(dates_spa)):
    dates_spa_str.append(str(dates_spa[i]))
    
first_july = np.where(pd.Series(dates_spa_str).str.contains('07-01') == True)[0][0]
#first_july_c = ix + first_july


# WB from first first of july
results_2d['WB_WS'] = ((results_2d['soil_water_storage'][first_july:] - results_2d['soil_water_storage'][first_july] 
                       + results_2d['bucket_water_storage'][first_july:] - results_2d['bucket_water_storage'][first_july]
                       + results_2d['canopy_snow_water_equivalent'][first_july:] - results_2d['canopy_snow_water_equivalent'][first_july]).mean(['i','j']))


results_2d['WB_ET'] = (np.cumsum(results_2d['canopy_evaporation'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_2d['canopy_transpiration'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_2d['bucket_evaporation'][first_july:].mean(['i','j'])))
                    
results_2d['WB_RO'] = (np.cumsum(results_2d['bucket_surface_runoff'][first_july:].mean(['i','j'])) + 
                       np.cumsum((results_2d['soil_netflow_to_ditch'][first_july:]*results_2d['parameters_cmask']).mean(['i','j'])))
results_2d['WB_P'] = np.cumsum(results_2d['forcing_precipitation'][first_july:])

results_2d['WB_date'] = results_2d['date'][first_july:]

# stand
results_stand['WB_WS'] = ((results_stand['bucket_water_storage'][first_july:] - results_stand['bucket_water_storage'][first_july]).mean(['i','j'])
                        + (results_stand['canopy_snow_water_equivalent'][first_july:] - results_stand['canopy_snow_water_equivalent'][first_july]).mean(['i','j']))

results_stand['WB_ET'] = (np.cumsum(results_stand['canopy_evaporation'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_stand['canopy_transpiration'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_stand['bucket_evaporation'][first_july:].mean(['i','j'])))
results_stand['WB_RO'] = (np.cumsum(results_stand['bucket_surface_runoff'][first_july:].mean(['i','j'])) + 
                       np.cumsum((results_stand['bucket_drainage'][first_july:]*results_2d['parameters_cmask']).mean(['i','j'])))
results_stand['WB_P'] = np.cumsum(results_stand['forcing_precipitation'][first_july:])

results_stand['WB_date'] = results_stand['date'][first_july:]


# catch

results_catch['WB_WS'] = ((results_catch['bucket_water_storage'][first_july:] - results_catch['bucket_water_storage'][first_july]).mean(['i','j'])
                        + (results_catch['canopy_snow_water_equivalent'][first_july:] - results_catch['canopy_snow_water_equivalent'][first_july]).mean(['i','j']))

results_catch['WB_ET'] = (np.cumsum(results_catch['canopy_evaporation'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_catch['canopy_transpiration'][first_july:].mean(['i','j'])) 
                       + np.cumsum(results_catch['bucket_evaporation'][first_july:].mean(['i','j'])))
results_catch['WB_RO'] = (np.cumsum(results_catch['bucket_surface_runoff'][first_july:].mean(['i','j'])) + 
                       np.cumsum((results_catch['top_baseflow'][first_july:]*results_2d['parameters_cmask']).mean(['i','j'])))
results_catch['WB_P'] = np.cumsum(results_catch['forcing_precipitation'][first_july:])

results_catch['WB_date'] = results_catch['date'][first_july:]

'''
#results_catch = results_stand.copy()
results_catch = pd.DataFrame()
results_catch['WB_WS'] = np.nanmean((bu_catch['WatSto'][first_july_c:]*1e3 - bu_catch['WatSto'][first_july_c]*1e3)
                          + cpy_catch['SWE'][first_july_c:] - cpy_catch['SWE'][first_july_c], axis=(1,2))

results_catch['WB_ET'] = (np.cumsum(np.nanmean(cpy_catch['Evap'][first_july_c:], axis=(1,2))) 
                       + np.cumsum(np.nanmean(cpy_catch['Transpi'][first_july_c:], axis=(1,2))) 
                       + np.cumsum(np.nanmean(cpy_catch['Efloor'][first_july_c:], axis=(1,2))))

results_catch['WB_RO'] = (np.cumsum(top_catch['Qt'][first_july_c:]*1e3))
                       
results_catch['WB_P'] = np.cumsum(results_stand['forcing_precipitation'][first_july:])

results_catch['WB_date'] = results_stand['date'][first_july:]
'''

nullref = [0]*len(results_stand['WB_date'])
nullrefc = [0]*len(results_catch['WB_date'])


#%%

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12,12));
ax1 = axs[2]
ax2 = axs[0]
ax3 = axs[1]

ax1.text(dates_spa[int(len(dates_spa)/2)], 2100, '2D')
ax1.plot(results_2d['WB_date'], 
         results_2d['WB_WS'], 'b', label='Water and snow storage')
ax1.plot(results_2d['WB_date'],
         results_2d['WB_WS']+ results_2d['WB_ET'], 
         'darkgreen', label='+ Evapotranspiration')
ax1.plot(results_2d['WB_date'],
         results_2d['WB_WS'] + results_2d['WB_ET'] + results_2d['WB_RO'],
         'brown', alpha=0.3, label='+ Discharge')
ax1.plot(results_2d['WB_date'],results_2d['WB_P'],'k.', alpha=0.5, markersize=3, label='Precipitation')
ax1.fill_between(results_2d['WB_date'], nullref, results_2d['WB_WS'], color='blue', alpha=0.2)
ax1.fill_between(results_2d['WB_date'], results_2d['WB_WS'], results_2d['WB_WS'] + results_2d['WB_ET'], color='green', alpha=0.2)
ax1.fill_between(results_2d['WB_date'], results_2d['WB_WS'] + results_2d['WB_ET'], results_2d['WB_WS'] + results_2d['WB_ET']
                 + results_2d['WB_RO'], color='brown', alpha=0.2)
ax1.set(ylim=(0, 2200))

#ax1.legend()

ax2.text(dates_spa[int(len(dates_spa)/2)], 2100, 'stand')
ax2.plot(results_stand['WB_date'], 
         results_stand['WB_WS'], 'b', label='Water and snow storage')
ax2.plot(results_stand['WB_date'],
         results_stand['WB_WS']+ results_stand['WB_ET'], 
         'darkgreen', label='+ Evapotranspiration')
ax2.plot(results_stand['WB_date'],
         results_stand['WB_WS'] + results_stand['WB_ET'] + results_stand['WB_RO'],
         'brown', alpha=0.3, label='+ Discharge')
ax2.plot(results_stand['WB_date'],results_stand['WB_P'],'k.', alpha=0.5, markersize=3, label='Precipitation')
ax2.fill_between(results_stand['WB_date'], nullref, results_stand['WB_WS'], color='blue', alpha=0.2)
ax2.fill_between(results_stand['WB_date'], results_stand['WB_WS'], results_stand['WB_WS'] + results_stand['WB_ET'], color='green', alpha=0.2)
ax2.fill_between(results_stand['WB_date'], results_stand['WB_WS'] + results_stand['WB_ET'], results_stand['WB_WS'] + results_stand['WB_ET']
                 + results_stand['WB_RO'], color='brown', alpha=0.2)
ax2.axes.get_xaxis().set_visible(False)
ax2.legend()
ax2.set(ylim=(0, 2200))


ax3.text(dates_spa[int(len(dates_spa)/2)], 2100, 'catch')
ax3.plot(results_catch['WB_date'], 
         results_catch['WB_WS'], 'b', label='Water and snow storage')
ax3.plot(results_catch['WB_date'],
         results_catch['WB_WS']+ results_catch['WB_ET'], 
         'darkgreen', label='+ Evapotranspiration')
ax3.plot(results_catch['WB_date'],
         results_catch['WB_WS'] + results_catch['WB_ET'] + results_catch['WB_RO'],
         'brown', alpha=0.3, label='+ Discharge')
ax3.plot(results_catch['WB_date'],results_catch['WB_P'],'k.', alpha=0.5, markersize=3, label='Precipitation')
ax3.fill_between(results_catch['WB_date'], nullrefc, results_catch['WB_WS'], color='blue', alpha=0.2)
ax3.fill_between(results_catch['WB_date'], results_catch['WB_WS'], results_catch['WB_WS'] + results_catch['WB_ET'], color='green', alpha=0.2)
ax3.fill_between(results_catch['WB_date'], results_catch['WB_WS'] + results_catch['WB_ET'], results_catch['WB_WS'] + results_catch['WB_ET']
                 + results_catch['WB_RO'], color='brown', alpha=0.2)
#ax3.legend()
ax3.axes.get_xaxis().set_visible(False)
ax3.set(ylim=(0, 2200))

plt.subplots_adjust(wspace=1, hspace=0)

if saveplots == True:
        plt.savefig(f'WB_part_{today}.pdf')
        plt.savefig(f'WB_part_{today}.png')


#%%

plt.plot(dates_spa, np.cumsum(q), label='obs')
plt.plot(dates_spa, results_catch['WB_RO'], label='catch')
plt.plot(dates_spa, results_stand['WB_RO'], label='stand')
plt.plot(dates_spa, results_2d['WB_RO'], label='2D')
plt.legend()





#%%

# soil moist plots without topsoil

# spafhy
wet_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmax()
dry_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmin()

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,12));
ax1 = axs[0][2]
ax2 = axs[1][2]

ax3 = axs[0][0]
ax4 = axs[1][0]

ax5 = axs[0][1]
ax6 = axs[1][1]

# 2D root moist
im1 = ax1.imshow(results_2d['bucket_moisture_root'][wet_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im3, ax=ax3)
ax1.title.set_text('2D root wet')

im2 = ax2.imshow(results_2d['bucket_moisture_root'][dry_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im2, cax=cbar_ax, label=r'$\theta$ m$^3$m$^{-3}$')

#fig.colorbar(im2, ax=ax2, orientation='horizontal')
ax2.title.set_text('2D root dry')

# stand root moist
im3 = ax3.imshow(results_stand['bucket_moisture_root'][wet_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im7, ax=ax7)
ax3.title.set_text('stand root wet')

im4 = ax4.imshow(results_stand['bucket_moisture_root'][dry_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
#fig.colorbar(im4, ax=ax4)
ax4.title.set_text('stand root dry')

im5 = ax5.imshow(results_catch['bucket_moisture_root'][wet_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
ax5.title.set_text('catch root wet')

im6 = ax6.imshow(results_catch['bucket_moisture_root'][dry_day,zy,zx], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
ax6.title.set_text('catch root dry')

ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
ax4.axis('off')
ax5.axis('off')
ax6.axis('off')


if saveplots == True:
        plt.savefig(f'theta_wet_dry_{today}.pdf')
        plt.savefig(f'theta_wet_dry_{today}.png')

#%%

# COMPARISON BETWEEN MODELS
# SUBSTRACTED IN JUNE AND SEPT

# spafhy
#wet_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmax()
#dry_day = np.nansum(results_2d['bucket_moisture_root'], axis=(1,2)).argmin()
june = np.where(pd.to_datetime(dates_spa) == '2021-06-17')[0][0]
sept = np.where(pd.to_datetime(dates_spa) == '2021-09-01')[0][0]

#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2_direct10_16.nc'
#sar_path = r'C:\PALLAS_RAW_DATA\SAR_maankosteus\processed\SAR_PALLAS_2019_mask2.nc'
#sar = Dataset(sar_path, 'r')

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14,12));
ax1 = axs[0][2]
ax2 = axs[1][2]

ax3 = axs[0][0]
ax4 = axs[1][0]

ax5 = axs[0][1]
ax6 = axs[1][1]

# 2D root moist
im1 = ax1.imshow(results_stand['bucket_moisture_root'][june,zy,zx]
                 - results_catch['bucket_moisture_root'][june,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
#fig.colorbar(im3, ax=ax3)
ax1.title.set_text('stand - catch (june)')

im2 = ax2.imshow(results_stand['bucket_moisture_root'][sept,zy,zx]
                 - results_catch['bucket_moisture_root'][sept,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im2, cax=cbar_ax, label=r'$\theta$ m$^3$m$^{-3}$')

#fig.colorbar(im2, ax=ax2, orientation='horizontal')
ax2.title.set_text('stand - catch (sept)')

# stand root moist
im3 = ax3.imshow(results_2d['bucket_moisture_root'][june,zy,zx] 
                 - results_stand['bucket_moisture_root'][june,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
#fig.colorbar(im7, ax=ax7)
ax3.title.set_text('2d - stand (june)')

im4 = ax4.imshow(results_2d['bucket_moisture_root'][sept,zy,zx] 
                 - results_stand['bucket_moisture_root'][sept,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
#fig.colorbar(im4, ax=ax4)
ax4.title.set_text('2d - stand (sept)')

im5 = ax5.imshow(results_2d['bucket_moisture_root'][june,zy,zx] 
                 - results_catch['bucket_moisture_root'][june,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
ax5.title.set_text('2d - catch (june)')

im6 = ax6.imshow(results_2d['bucket_moisture_root'][sept,zy,zx] 
                 - results_catch['bucket_moisture_root'][sept,zy,zx], vmin=-0.4, vmax=0.4, cmap='coolwarm_r')
ax6.title.set_text('2d - catch (sept)')

ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
ax4.axis('off')
ax5.axis('off')
ax6.axis('off')


if saveplots == True:
        plt.savefig(f'theta_diff_june_sept_{today}.pdf')
        plt.savefig(f'theta_diff_june_sept_{today}.png')
#%%

soilclass_2 = np.ravel(soilclass)
soilclass_4 = np.ravel(soilclass_copy)
r, c = np.shape(soilclass)
soilclass_2[soilclass_2 != 2] = np.nan
soilclass_2[soilclass_2 == 2] = 1
soilclass_4[soilclass_4 != 4] = np.nan
soilclass_4[soilclass_4 == 4] = 1

soilclass_2 = soilclass_2.reshape(r, c)
soilclass_4 = soilclass_4.reshape(r, c)

td_sum2 = np.nanmean(results_2d['bucket_moisture_root'] * soilclass_2, axis=(1,2))
stand_sum2 = np.nanmean(results_stand['bucket_moisture_root'] * soilclass_2, axis=(1,2))
catch_sum2 = np.nanmean(results_catch['bucket_moisture_root'] * soilclass_2, axis=(1,2))
td_sum4 = np.nanmean(results_2d['bucket_moisture_root'] * soilclass_4, axis=(1,2))
stand_sum4 = np.nanmean(results_stand['bucket_moisture_root'] * soilclass_4, axis=(1,2))
catch_sum4 = np.nanmean(results_catch['bucket_moisture_root'] * soilclass_4, axis=(1,2))

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12,8));
ax1 = axs[0]
ax2 = axs[1]
sns.set_style('whitegrid')

ax1.set_title('Mineral soil')
ax1.plot(dates_spa, td_sum2, label='2D')
ax1.plot(dates_spa, stand_sum2, label='stand')
ax1.plot(dates_spa, catch_sum2, label='catch')
ax1.legend(ncol=3)
ax1.set_ylabel(r'$\theta_{mod}$ (m$^3$m$^{-3}$)')

ax2.set_title('Peat soil')
ax2.plot(dates_spa, td_sum4, label='2D')
ax2.plot(dates_spa, stand_sum4, label='stand')
ax2.plot(dates_spa, catch_sum4, label='catch')
ax2.legend(ncol=3)
ax2.set_ylabel(r'$\theta_{mod}$ (m$^3$m$^{-3}$)')

if saveplots == True:
        plt.savefig(f'theta_mod_mean_ts_{today}.pdf')
        plt.savefig(f'theta_mod_mean_ts_{today}.png')


#%%

plt.hist(np.array(results_2d['bucket_moisture_root'][:,zy,zx]).flatten(), bins=np.arange(0.1,1,0.01)); plt.ylim(0, 1e7)
plt.hist(np.array(results_stand['bucket_moisture_root'][:,zy,zx]).flatten(), bins=np.arange(0.1,1,0.01)); plt.ylim(0, 1e7)
plt.hist(np.array(results_catch['bucket_moisture_root'][:,zy,zx]).flatten(), bins=np.arange(0.1,1,0.01)); plt.ylim(0, 1e7)


#%%
'''
# Fixing random state for reproducibility
np.random.seed(19680801)

# some random data
#x = np.random.randn(1000)
#y = np.random.randn(1000)

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.01
    xymax = max(np.nanmax(x), np.nanmax(y))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(0.1, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')
    

x = list(soilm['spa_k_st_root'][np.isfinite(soilm['spa_k_st_root'])]); y = list(soilm['mean'][np.isfinite(soilm['spa_k_st_root'])])
m, b = np.polyfit(x[0:500], y[0:500], 1)

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

# Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

# use the previously defined function
scatter_hist(x, y, ax, ax_histx, ax_histy)
ax.plot(x, m*x[0:500]+b)

plt.show()

#f3_ax2 = fig3.add_subplot(gs[0, 3])
#ax_histx = fig.add_subplot(gs[0, 3], sharex=ax)
#ax_histy = fig.add_subplot(gs[0, 3], sharey=ax)
#obtain m (slope) and b(intercept) of linear regression line

#add linear regression line to scatterplot 
#scatter_hist(x, y, ax, ax_histx, ax_histy)
'''

#%%

# day with average bucket water storage as reference

stand_av = np.nansum(results_stand['bucket_water_storage'], axis=(1,2))
stand_av_ix = np.where(stand_av == np.median(stand_av))[0][0]
catch_av = np.nansum(results_catch['bucket_water_storage'], axis=(1,2))
catch_av_ix = np.where(catch_av == np.median(catch_av))[0][0]
td_av = np.nansum(results_2d['bucket_water_storage'], axis=(1,2))
td_av_ix = np.where(td_av == np.median(td_av))[0][0]





#%%

# preparing soil moisture datas
# point examples from mineral and openmire
# soilscouts at Kenttarova
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
soil_file = 'soilm_kenttarova.csv'
fp = os.path.join(folder, soil_file)
soilm = pd.read_csv(fp, sep=';', date_parser=['time'])
soilm['time'] = pd.to_datetime(soilm['time'])
soilm['mean'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A']].mean(numeric_only=True, axis=1)
#soilm['sd_20'] = soilm[['SH-20A', 'SH-20B', 's18']].std(skipna=True, numeric_only=True, axis=1)
soilm['min'] = soilm[['SH-20A', 'SH-20B', 's18', 's3', 'SH-5A', 'SH-5B']].min(skipna=True, numeric_only=True, axis=1)
soilm['max'] = soilm[['SH-20A', 'SH-20B', 's18', 's3', 'SH-5A', 'SH-5B']].max(skipna=True, numeric_only=True, axis=1)
soilm['iq25'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A']].quantile(q=0.25, numeric_only=True, axis=1)
soilm['iq75'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A']].quantile(q=0.75, numeric_only=True, axis=1)

# when freezing temp
tneg = np.where(results_stand['forcing_air_temperature'] < 0)[0]
tpos =  np.where(results_stand['forcing_air_temperature'] >= 0)[0]

#l_loc = [60, 60]
spa_wliq_2d_root = results_2d['bucket_moisture_root']
spa_wliq_2d_top = results_2d['bucket_moisture_top']

spa_wliq_st_root = results_stand['bucket_moisture_root']
spa_wliq_st_top = results_stand['bucket_moisture_top']

spa_wliq_ca_root = results_stand['bucket_moisture_root']
spa_wliq_ca_top = results_stand['bucket_moisture_top']

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
spa_wliq_df['spa_k_ca_root'] = spa_wliq_ca_root[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_root'] = spa_wliq_ca_root[:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_ca_top'] = spa_wliq_ca_top[:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_top'] = spa_wliq_ca_top[:,l_loc[0],l_loc[1]]

spa_wliq_df['time'] = dates_spa
#spa_wliq_df.index = dates_spa

#soilm = pd.concat([soilm, spa_wliq_df]).sort_values('time').reset_index(drop=True)
soilm = spa_wliq_df.merge(soilm)
soilm.index = soilm['time']
soilm['month'] = pd.DatetimeIndex(soilm.index).month
winter_ix = np.where((soilm['month'] <= 4) | (soilm['month'] >= 11))[0]
soilm.iloc[winter_ix] = np.nan
#soilm = soilm[['s3', 's5', 's18', 'SH-5A', 'SH-5B', 'SH-20A', 'SH-20B', 'spa_k', 'spa_l', 'spatop_k', 'spatop_l']]

#soilm = soilm.loc[(soilm.index > '2018-04-01') & (soilm.index < '2019-12-01')]

poi = 1085
poi = np.where(results_stand['bucket_moisture_root'][:,k_loc[0], k_loc[1]] + 0.09 < results_2d['bucket_moisture_root'][:,k_loc[0], k_loc[1]])[0][0]

#%%

# preparing gw data




#%%
# Plotting soil moisture comparison

###
fig3 = plt.figure(constrained_layout=True, figsize=(14,7))
gs = fig3.add_gridspec(2, 4)
sns.set_style('whitegrid')

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Kenttärova')
im1 = f3_ax1.plot(soilm['spa_k_2d_root'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax1.plot(soilm['spa_k_ca_root'], 'red', alpha=0.6, label='stand/catch')
f3_ax1.plot(soilm['mean'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax1.fill_between(soilm.index, soilm['iq25'], soilm['iq75'], color='blue', alpha=0.2, label=r'IQR$_{obs}$')

f3_ax1.legend(ncol=5,bbox_to_anchor=(0.8, 1.3))
y = f3_ax1.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#y.set_rotation(0)
f3_ax1.set_ylim(0.1,0.5)
#f3_ax1.axes.get_xaxis().set_visible(False)

f3_ax2 = fig3.add_subplot(gs[0, 3])
x4 = sns.regplot(ax=f3_ax2, x=soilm['spa_k_st_root'], y=soilm['mean'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0.1, 0.45))
f3_ax2.set(xlim=(0.1, 0.45))
f3_ax2.yaxis.tick_right()
f3_ax2.set_ylabel(r'$\theta_{obs}$ (m$^3$m$^{-3}$)')
f3_ax2.set_xlabel(r'$\theta_{mod}$ (m$^3$m$^{-3}$)')



f3_ax3 = fig3.add_subplot(gs[1, :3])
f3_ax3.set_title('Lompolonjänkkä')
im2 = f3_ax3.plot(soilm['spa_l_2d_root'],  'g', alpha=0.7, label='2D')
f3_ax3.plot(soilm['spa_l_st_root'], 'red', alpha=0.6, label='stand/catch')
#f3_ax3.plot(soilm['spa_l_ca_root'],  'red', alpha=0.6, label='catch')
#ax2.plot(soilm.index[poi], 0.45, marker='o', mec='k', mfc='g', alpha=0.5, ms=8.0)
#ax2.axvline(soilm.index[poi], ymin=0, ymax=1, color='k', alpha=0.4)
f3_ax3.text(dates_spa[1], 1.03, 'Mire')
#ax2.title.set_text('Mire')
f3_ax3.set_ylim(0.4,1.0)
f3_ax3.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#f3_ax3.legend(ncol=5)


if saveplots == True:
        plt.savefig(f'theta_model_ts_{today}.pdf')
        plt.savefig(f'theta_model_ts_{today}.png')

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
ax1.scatter(k_loc[1],k_loc[0],color='r', alpha=0.4)
ax1.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax1.title.set_text(f'2D root {dates_spa[poi]}')

# stand moist
im2 = ax2.imshow(results_stand['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im2, ax=ax2)
ax2.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax2.scatter(k_loc[1],k_loc[0],color='r', alpha=0.4)
ax2.title.set_text(f'stand root {dates_spa[poi]}')

# catchment moist
im3 = ax3.imshow(results_catch['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
fig.colorbar(im3, ax=ax3)
ax3.scatter(l_loc[1],l_loc[0],color='k', alpha=0.4)
ax3.scatter(k_loc[1],k_loc[0],color='r', alpha=0.4)
ax3.title.set_text(f'catch root {dates_spa[poi]}')

#%%
# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16,16));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[1][0]
ax5 = axs[1][1]
ax6 = axs[1][2]

sns.regplot(ax=ax3, x=np.array(results_2d['bucket_moisture_root'][900:1000,:,:]).flatten(), 
            y=np.array(results_stand['bucket_moisture_root'][900:1000,:,:]).flatten(), 
            scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
sns.regplot(ax=ax2, x=np.array(results_2d['bucket_moisture_root'][poi,:,:]).flatten(), 
            y=np.array(results_catch['bucket_moisture_root'][poi,:,:]).flatten(), 
            scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
sns.regplot(ax=ax1, x=np.array(results_stand['bucket_moisture_root'][poi,:,:]).flatten(), 
            y=np.array(results_catch['bucket_moisture_root'][poi,:,:]).flatten(), 
            scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})

ax4.imshow(results_stand['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
ax5.imshow(results_catch['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')
ax6.imshow(results_2d['bucket_moisture_root'][poi,:,:], vmin=0.0, vmax=1.0, cmap='coolwarm_r')


#%%
# DRY ET DATA PREPARING

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
ec['ET-1-KR'].iloc[(ec['time'] < '2017-01-01') & (ec['time'] > '2016-01-01')] = np.nan
#ec['ET-1-LV'].iloc[(ec['time'] < '2017-01-01') & (ec['time'] > '2016-01-01')] = np.nan
ec['ET-1-KR'].iloc[ix_p] = np.nan; ec['ET-1-LV'].iloc[ix_p] = np.nan
ec['ET-1-KR'].iloc[ix_pp] = np.nan; ec['ET-1-LV'].iloc[ix_pp] = np.nan
ec['ET-1-KR'].loc[ec['ET-1-KR'] < 0] = np.nan
ec['ET-1-LV'].loc[ec['ET-1-LV'] < 0] = np.nan

#dry et
results_2d['dry_et'] = results_2d['bucket_evaporation'] + results_2d['canopy_transpiration'] #+ results_2d['canopy_evaporation']
#results_2d['dry_et'][ix_p,:,:] = np.nan
results_stand['dry_et'] = results_stand['bucket_evaporation'] + results_stand['canopy_transpiration'] #+ results_stand['canopy_evaporation'] 
#results_stand['dry_et'][ix_p,:,:] = np.nan
results_catch['dry_et'] = results_catch['bucket_evaporation'] + results_catch['canopy_transpiration'] #+ results_catch['canopy_evaporation']
#catch_dry_et[ix_p,:,:] = np.nan

# scatter et sim vs. obs
ET_flat = pd.DataFrame()
ET_flat['stand_kr_et'] = results_stand['dry_et'][:,k_loc[0],k_loc[0]]
ET_flat['stand_lv_et'] = results_stand['dry_et'][:,l_loc[0],l_loc[0]]
ET_flat['2D_kr_et'] = results_2d['dry_et'][:,k_loc[0],k_loc[0]]
ET_flat['2D_lv_et'] = results_2d['dry_et'][:,l_loc[0],l_loc[0]] 
ET_flat['catch_kr_et'] = results_catch['dry_et'][:,k_loc[0],k_loc[0]] 
ET_flat['catch_lv_et'] = results_catch['dry_et'][:,l_loc[0],l_loc[0]]
ET_flat['ET-1-KR'] = ec['ET-1-KR'].reset_index(drop=True)
ET_flat['ET-1-LV'] = ec['ET-1-LV'].reset_index(drop=True)
ET_flat = ET_flat.dropna()
ET_flat = ET_flat[(ET_flat > 0).all(1)]


#%%

# Plotting Kenttärova ET mod vs. obs

fig3 = plt.figure(constrained_layout=True, figsize=(16,10))
gs = fig3.add_gridspec(3, 4)

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Kenttärova - modelled vs. observed dry ET')
f3_ax1.plot(dates_spa, results_stand['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.9)
f3_ax1.plot(ec['ET-1-KR'], 'k.',  alpha=0.4,markersize=10)
f3_ax1.legend(['stand', 'ec'], loc='upper right')
f3_ax1.axes.get_xaxis().set_visible(False)
f3_ax1.set(ylim=(-0.5, 9))


f3_ax2 = fig3.add_subplot(gs[0, 3])
x2 = sns.regplot(ax=f3_ax2, x=ET_flat['stand_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0, 5))
f3_ax2.set(xlim=(0, 5))
#f3_ax2.axes.get_xaxis().set_visible(False)
f3_ax2.yaxis.tick_right()


f3_ax3 = fig3.add_subplot(gs[1, :3])
#f3_ax3.set_title('catch vs. obs')
f3_ax3.plot(dates_spa, results_catch['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.9)
f3_ax3.plot(ec['ET-1-KR'], 'k.',  alpha=0.4,markersize=10)
f3_ax3.legend(['catch', 'ec'], loc='upper right')
f3_ax3.axes.get_xaxis().set_visible(False)
f3_ax3.set(ylim=(-0.5, 9))

f3_ax4 = fig3.add_subplot(gs[1, 3])
x6 = sns.regplot(ax=f3_ax4, x=ET_flat['catch_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax4.set(ylim=(0, 5))
f3_ax4.set(xlim=(0, 5))
#f3_ax4.axes.get_xaxis().set_visible(False)
f3_ax4.yaxis.tick_right()


f3_ax5 = fig3.add_subplot(gs[2, :3])
#f3_ax5.set_title('2D vs. obs')
f3_ax5.plot(dates_spa, results_2d['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.9)
f3_ax5.plot(ec['ET-1-KR'], 'k.',  alpha=0.4,markersize=10)
f3_ax5.legend(['2d', 'ec'], loc='upper right')
f3_ax5.set(ylim=(-0.5, 9))


f3_ax6 = fig3.add_subplot(gs[2, 3])
x4 = sns.regplot(ax=f3_ax6, x=ET_flat['2D_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax6.set(ylim=(0, 5))
f3_ax6.set(xlim=(0, 5))
f3_ax6.yaxis.tick_right()

if saveplots == True:
        plt.savefig(f'qq_spa_spatop_sar_{today}.pdf')
        plt.savefig(f'ET_MOD_OBS_KR_{today}.png')

#%%
# Plotting Lompolonjänkkä ET mod vs. obs

fig3 = plt.figure(constrained_layout=True, figsize=(16,10))
gs = fig3.add_gridspec(3, 4)

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Lompolonjänkkä - modelled vs. observed dry ET')
f3_ax1.plot(dates_spa, results_stand['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.9)
f3_ax1.plot(ec['ET-1-LV'], 'k.', alpha=0.4,markersize=10)
f3_ax1.legend(['stand', 'ec'], loc='upper right')
f3_ax1.axes.get_xaxis().set_visible(False)
f3_ax1.axes.get_xaxis().set_visible(False)
f3_ax1.set(ylim=(-0.5, 9))

f3_ax2 = fig3.add_subplot(gs[0, 3])
x2 = sns.regplot(ax=f3_ax2, x=ET_flat['stand_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0, 5))
f3_ax2.set(xlim=(0, 5))
f3_ax2.yaxis.tick_right()


f3_ax3 = fig3.add_subplot(gs[1, :3])
#f3_ax3.set_title('catch vs. obs')
f3_ax3.plot(dates_spa, results_catch['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.9)
f3_ax3.plot(ec['ET-1-LV'], 'k.', alpha=0.4,markersize=10)
f3_ax3.legend(['catch', 'ec'], loc='upper right')
f3_ax3.axes.get_xaxis().set_visible(False)
f3_ax3.axes.get_xaxis().set_visible(False)
f3_ax3.set(ylim=(-0.5, 9))

f3_ax4 = fig3.add_subplot(gs[1, 3])
x6 = sns.regplot(ax=f3_ax4, x=ET_flat['catch_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax4.set(ylim=(0, 5))
f3_ax4.set(xlim=(0, 5))
f3_ax4.yaxis.tick_right()

f3_ax5 = fig3.add_subplot(gs[2, :3])
#f3_ax5.set_title('2D vs. obs')
f3_ax5.plot(dates_spa, results_2d['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.9)
f3_ax5.plot(ec['ET-1-LV'], 'k.',  alpha=0.4,markersize=10)
f3_ax5.legend(['2d', 'ec'], loc='upper right')
f3_ax5.set(ylim=(-0.5, 9))

f3_ax6 = fig3.add_subplot(gs[2, 3])
x4 = sns.regplot(ax=f3_ax6, x=ET_flat['2D_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax6.set(ylim=(0, 5))
f3_ax6.set(xlim=(0, 5))
f3_ax6.yaxis.tick_right()

if saveplots == True:
        plt.savefig(f'ET_MOD_OBS_LV_{today}.pdf')
        plt.savefig(f'ET_MOD_OBS_LV_{today}.png')
        
        
#%%

# KENTTÄROVA AND LOMPOLO VS STAND


fig3 = plt.figure(constrained_layout=True, figsize=(14,7))
gs = fig3.add_gridspec(2, 4)

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Kenttärova')
f3_ax1.plot(dates_spa, results_stand['dry_et'][:,k_loc[0],k_loc[1]], alpha=0.9)
#f3_ax1.plot(dates_spa, results_stand['forcing_vapor_pressure_deficit'])
f3_ax1.plot(ec['ET-1-KR'], 'k.',  alpha=0.4,markersize=10)
f3_ax1.legend(['stand', 'ec'], loc='upper right')
f3_ax1.set(ylim=(-0.5, 9))
f3_ax1.set_ylabel(r'ET mm d$^{-1}$')

f3_ax2 = fig3.add_subplot(gs[0, 3])
x2 = sns.regplot(ax=f3_ax2, x=ET_flat['stand_kr_et'], y=ET_flat['ET-1-KR'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0, 5))
f3_ax2.set(xlim=(0, 5))
#f3_ax2.axes.get_xaxis().set_visible(False)
f3_ax2.yaxis.tick_right()
f3_ax2.set_ylabel(r'ET$_{obs}$ (mm d$^{-1}$)')
f3_ax2.set_xlabel(r'ET$_{mod}$ (mm d$^{-1}$)')


# Plotting Lompolonjänkkä ET mod vs. obs

f3_ax3 = fig3.add_subplot(gs[1, :3])
f3_ax3.set_title('Lompolonjänkkä')
f3_ax3.plot(dates_spa, results_stand['dry_et'][:,l_loc[0],l_loc[1]], alpha=0.9)
f3_ax3.plot(ec['ET-1-LV'], 'k.', alpha=0.4,markersize=10)
f3_ax3.legend(['stand', 'ec'], loc='upper right')
f3_ax3.set(ylim=(-0.5, 9))
f3_ax3.set_ylabel(r'ET mm d$^{-1}$')

f3_ax4 = fig3.add_subplot(gs[1, 3])
x4 = sns.regplot(ax=f3_ax4, x=ET_flat['stand_lv_et'], y=ET_flat['ET-1-LV'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax4.set(ylim=(0, 5))
f3_ax4.set(xlim=(0, 5))
f3_ax4.yaxis.tick_right()
f3_ax4.set_ylabel(r'ET$_{obs}$ (mm d$^{-1}$)')
f3_ax4.set_xlabel(r'ET$_{mod}$ (mm d$^{-1}$)')

if saveplots == True:
        plt.savefig(f'ET_MOD_OBS_KR_LV_{today}.pdf')
        plt.savefig(f'ET_MOD_OBS_KR_LV_{today}.png')

#%%
# SWE file

fn = r'C:\SpaFHy_v1_Pallas\data\obs\SWE_survey_2018-02-22_2021-05-16.txt'
SWE_m = pd.read_csv(fn, skiprows=5, sep=';', parse_dates = ['date'], encoding='iso-8859-1')
SWE_m.index = SWE_m['date']
#SWE_m['SWE'].loc[SWE_m['date'] < '2018-07-07'] = np.nan
SWE_m = SWE_m[['SWE', 'SWE_sd', 'quality']]

SWE = pd.DataFrame()
SWE['mod_mean'] = np.nanmean(results_stand['canopy_snow_water_equivalent'], axis=(1,2))
SWE.index = dates_spa
SWE_all = SWE_m.join(SWE)
SWE_all.loc[np.isnan(SWE_all['SWE'])] = np.nan
SWE_all.loc[np.isnan(SWE_all['mod_mean'])] = np.nan
SWE_all = SWE_all.dropna()

#%%

fig3 = plt.figure(constrained_layout=True, figsize=(14,3.5))
gs = fig3.add_gridspec(1, 4)

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.plot(SWE_m['SWE'], 'k.', markersize=6)
#f3_ax1.plot(SWE_m['SWE'] + SWE_m['SWE_sd'], 'k', alpha=0.5)
#f3_ax1.plot(SWE_m['SWE'] - SWE_m['SWE_sd'], 'k', alpha=0.5)
f3_ax1.plot(dates_spa, np.nanmean(results_stand['canopy_snow_water_equivalent'], axis=(1,2)))
f3_ax1.plot(forc['snowdepth'], 'r.', markersize=2, alpha=0.2)
f3_ax1.legend(['obs SWE', 'mod SWE', 'obs HS'],ncol=3,  loc='upper left')
f3_ax1.set_ylabel('SWE (mm)')


f3_ax2 = fig3.add_subplot(gs[0, 3])
line1 = sns.regplot(ax=f3_ax2, x=SWE_all['SWE'], y=SWE_all['mod_mean'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
f3_ax2.set_ylabel('mod SWE (mm)')
f3_ax2.set_xlabel('obs SWE (mm)')
f3_ax2.yaxis.tick_right()
f3_ax2.yaxis.set_label_position("right")
sns.set_style('whitegrid')


#f3_ax2 = sns.lmplot(x="SWE", y="mod_mean", data=SWE_all)
#ax = lm.axes
f3_ax2.set(ylim=(0, 350))
f3_ax2.set(xlim=(0, 350))

if saveplots == True:
        plt.savefig(f'SWE_MOD_OBS_{today}.pdf')
        plt.savefig(f'SWE_MOD_OBS_{today}.png')
        
        
#%%

# specific discharge

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(14,6));
ax1 = axs[0]
ax2 = axs[1]

#ax1.plot(dates_spa, np.nanmean(results_stand['bucket_drainage'], axis=(1,2)) + np.nanmean(results_stand['bucket_surface_runoff'], axis=(1,2)))
#ax1.plot(q)
#ax1.legend(['stand SD', 'obs SD'], loc='upper right')
#ax1.set(ylim=(0, 10))
sns.set_style('whitegrid')

ax1.plot(dates_spa, results_catch['top_baseflow'] + np.nanmean(results_catch['bucket_surface_runoff'], axis=(1,2)))
ax1.plot(q)
ax1.legend(['catch', 'obs'], loc='upper right')
ax1.set(ylim=(-1, 25))
ax1.set_ylabel(r'Qf mm d$^{-1}$')
#ax1.axes.get_xaxis().set_visible(False)

ax2.plot(dates_spa, np.nanmean(results_2d['soil_netflow_to_ditch'], axis=(1,2)) + np.nanmean(results_2d['bucket_surface_runoff'], axis=(1,2)))
ax2.plot(q)
ax2.legend(['2D', 'obs'], loc='upper right')
ax2.set(ylim=(-1, 25))
ax2.set_ylabel(r'Qf mm d$^{-1}$')

if saveplots == True:
        plt.savefig(f'QF_MOD_OBS_{today}.pdf')
        plt.savefig(f'QF_MOD_OBS_{today}.png')
        
#%%
# area
area = len(np.where(np.isfinite(np.array(results_stand['parameters_cmask']).flatten()))[0]) * 16 * 16
# yearly water balance in meters for whole area
wbdf = pd.DataFrame()
wbdf['P'] = results_stand['forcing_precipitation']
wbdf.index = dates_spa
wbdf['Qmod'] = results_catch['top_baseflow'] + np.nanmean(results_catch['bucket_surface_runoff'], axis=(1,2))
wbdf['ETmod'] = np.nanmean(results_catch['canopy_evaporation'] + results_catch['canopy_transpiration'] + results_catch['bucket_evaporation'], axis=(1,2))
wbdf['ETdrymod'] = np.nanmean(results_catch['dry_et'], axis=(1,2))
wbdf['S'] = np.nanmean(results_catch['bucket_water_storage'], axis=(1,2)) + np.nanmean(results_catch['canopy_snow_water_equivalent'], axis=(1,2)) + np.nanmean(results_catch['canopy_water_storage'], axis=(1,2))
#wbdf['S'] = np.nanmean((results_catch['bucket_water_storage'] + results_catch['canopy_snow_water_equivalent'] + results_catch['canopy_water_storage'] + results_catch['bucket_pond_storage']), axis=(1,2))

wbdfy = pd.DataFrame()
wbdfy['P'] = wbdf['P'].resample('AS-SEP').sum()
wbdfy['Qmod'] = wbdf['Qmod'].resample('AS-SEP').sum()
wbdfy['ETdrymod'] = wbdf['ETdrymod'].resample('AS-SEP').sum()
wbdfy['ETmod'] = wbdf['ETmod'].resample('AS-SEP').sum()
wbdfy['S'] = np.nan
wbdfy['Qobs'] = q.resample('AS-SEP').sum()

wbdf2d = pd.DataFrame()
wbdf2d['P'] = results_2d['forcing_precipitation']
wbdf2d.index = dates_spa
wbdf2d['Qmod'] = np.nanmean(results_2d['soil_netflow_to_ditch'] + results_2d['bucket_surface_runoff'], axis=(1,2))
wbdf2d['ETdrymod'] = wbdf['ETdrymod'].resample('Y').sum()
wbdf2d['ETmod'] = np.nanmean(results_2d['canopy_evaporation'] + results_2d['canopy_transpiration'] + results_2d['bucket_evaporation'], axis=(1,2))
wbdf2d['S'] = np.nanmean(results_2d['bucket_water_storage'] + results_2d['canopy_snow_water_equivalent'] + results_2d['canopy_water_storage'] + results_2d['soil_water_storage'], axis=(1,2))

wbdf2dy = pd.DataFrame()
wbdf2dy['P'] = wbdf2d['P'].resample('AS-SEP').sum()
wbdf2dy['Qmod'] = wbdf2d['Qmod'].resample('AS-SEP').sum()
wbdf2dy['ETmod'] = wbdf2d['ETmod'].resample('AS-SEP').sum()
wbdf2dy['ETdrymod'] = wbdf2d['ETdrymod'].resample('AS-SEP').sum()
wbdf2dy['S'] = np.nan
wbdf2dy['Qobs'] = q.resample('AS-SEP').sum()

sday=1
eday=31
smonth=1
emonth=12
for i in range(len(wbdfy)):
    year = wbdfy.index[i].year
    #start = np.where((wbdf.index.year == year) & (wbdf.index.month == smonth) & (wbdf.index.day == sday))[0]
    #end = np.where((wbdf.index.year == year) & (wbdf.index.month == emonth) & (wbdf.index.day == eday))[0]
    start = np.where((wbdf.index.year == year) & (wbdf.index.month == 9) & (wbdf.index.day == 1))[0]
    end = np.where((wbdf.index.year == year + 1) & (wbdf.index.month == 9) & (wbdf.index.day == 1))[0]
    if len(start) + len(end) > 1:
        wbdfy['S'][i] = float(wbdf['S'][start]) - float(wbdf['S'][end])
        wbdf2dy['S'][i] = float(wbdf2d['S'][start]) - float(wbdf2d['S'][end])

wbdfy['closure'] = wbdfy['P'] - wbdfy['Qmod'] - wbdfy['ETmod'] + wbdfy['S']
wbdf2dy['closure'] = wbdf2dy['P'] - wbdf2dy['Qmod'] - wbdf2dy['ETmod'] + wbdf2dy['S']
wbdfy = wbdfy[1:-1] 
wbdf2dy = wbdf2dy[1:-1] 

wbplot = wbdfy.mean()

# Plotting
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,10));
ax1 = axs[0]
ax2 = axs[1]

ax1.bar('2013', wbdfy['P'][0], alpha=0.7, color='tab:blue', label=r'P$_{obs}$')
ax1.set_ylim(-850, 1000)
ax1.bar('2013', - wbdfy['ETmod'][0] - wbdfy['Qmod'][0], color='tab:green', alpha=0.7, label=r'ET$_{mod}$')
ax1.bar('2013', - wbdfy['Qmod'][0], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
ax1.scatter('2013', - wbdfy['Qobs'][0], s=200, marker='x', color='k', linewidth=4, zorder=2, label=r'$Q_{obs}$')
#ax1.scatter('2013', - wbdfy['ETmod'][0] - wbdfy['Qmod'][0] + wbdfy['S'][0] , s=200, marker='_', color='blue', linewidth=1, zorder=2, label=r'$+\Delta$S')
#ax1.legend(ncol=3)
ax1.legend(ncol=5,bbox_to_anchor=(1.6, 1.10))

ax1.bar('2014', wbdfy['P'][1], alpha=0.7, color='tab:blue', label=r'P$_{obs}$')
ax1.set_ylim(-850, 850)
ax1.bar('2014', - wbdfy['ETmod'][1] - wbdfy['Qmod'][1], color='tab:green', alpha=0.7, label=r'ET$_{mod}$')
ax1.bar('2014', - wbdfy['Qmod'][1], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
ax1.scatter('2014', - wbdfy['Qobs'][1], s=200, marker='x', color='k', linewidth=4, zorder=2, label=r'$Q_{obs}$')
#ax1.scatter('2014', - wbdfy['ETmod'][1] - wbdfy['Qmod'][1] + wbdfy['S'][1] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2015', wbdfy['P'][2], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2015', - wbdfy['ETmod'][2] - wbdfy['Qmod'][2], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2015', - wbdfy['Qmod'][2], alpha=0.7, color='tab:brown', label='Q')
ax1.scatter('2015', - wbdfy['Qobs'][2], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2015', - wbdfy['ETmod'][2] - wbdfy['Qmod'][2] + wbdfy['S'][2] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2016', wbdfy['P'][3], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2016', - wbdfy['ETmod'][3] - wbdfy['Qmod'][3], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2016', - wbdfy['Qmod'][3], alpha=0.7, color='tab:brown', label='Q')
ax1.scatter('2016', - wbdfy['Qobs'][3], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2016', - wbdfy['ETmod'][3] - wbdfy['Qmod'][3] + wbdfy['S'][3] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2017', wbdfy['P'][4], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2017', - wbdfy['ETmod'][4] - wbdfy['Qmod'][4], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2017', - wbdfy['Qmod'][4], alpha=0.7, color='tab:brown', label='Q')
#ax1.scatter('2017', - wbdfy['Qobs'][4], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2017', - wbdfy['ETmod'][4] - wbdfy['Qmod'][4] + wbdfy['S'][4] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2018', wbdfy['P'][5], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2018', - wbdfy['ETmod'][5] - wbdfy['Qmod'][5], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2018', - wbdfy['Qmod'][5], alpha=0.7, color='tab:brown', label='Q')
ax1.scatter('2018', - wbdfy['Qobs'][5], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2018', - wbdfy['ETmod'][5] - wbdfy['Qmod'][5] + wbdfy['S'][5] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2019', wbdfy['P'][6], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2019', - wbdfy['ETmod'][6] - wbdfy['Qmod'][6], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2019', - wbdfy['Qmod'][6], alpha=0.7, color='tab:brown', label='Q')
ax1.scatter('2019', - wbdfy['Qobs'][6], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2019', - wbdfy['ETmod'][6] - wbdfy['Qmod'][6] + wbdfy['S'][6] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax1.bar('2020', wbdfy['P'][7], alpha=0.7, color='tab:blue', label='P')
ax1.bar('2020', - wbdfy['ETmod'][7] - wbdfy['Qmod'][7], color='tab:green', alpha=0.7, label='ET')
ax1.bar('2020', - wbdfy['Qmod'][7], alpha=0.7, color='tab:brown', label='Q')
#ax1.scatter('2020', - wbdfy['Qobs'][7], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax1.scatter('2020', - wbdfy['ETmod'][7] - wbdfy['Qmod'][7] + wbdfy['S'][7] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

#
ax2.bar('2013', wbdf2dy['P'][0], alpha=0.7, color='tab:blue', label=r'P$_{obs}$')
ax2.set_ylim(-850, 1000)
#ax.bar('2014', - wbdfy['ETmod'][0] - wbplot['Qmod'], color='tab:brown', alpha=0.7, label='Q')
#ax.bar('2014', - wbdfy['ETmod'][0], alpha=0.7, color='tab:green', label='ET')
ax2.bar('2013', - wbdf2dy['ETmod'][0] - wbdf2dy['Qmod'][0], color='tab:green', alpha=0.7,label=r'ET$_{mod}$')
ax2.bar('2013', - wbdf2dy['Qmod'][0], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
ax2.scatter('2013', - wbdf2dy['Qobs'][0], s=200, marker='x', color='k', linewidth=4, zorder=2, label=r'$Q_{obs}$')
#ax2.scatter('2013', - wbdf2dy['ETmod'][0] - wbdf2dy['Qmod'][0] + wbdf2dy['S'][0] , s=200, marker='_', color='blue', linewidth=1, zorder=2, label=r'$+\Delta$S')
#ax2.legend(ncol=2)

ax2.bar('2014', wbdf2dy['P'][1], alpha=0.7, color='tab:blue', label=r'P$_{obs}$')
ax2.set_ylim(-850, 850)
#ax.bar('2014', - wbdfy['ETmod'][0] - wbplot['Qmod'], color='tab:brown', alpha=0.7, label='Q')
#ax.bar('2014', - wbdfy['ETmod'][0], alpha=0.7, color='tab:green', label='ET')
ax2.bar('2014', - wbdf2dy['ETmod'][1] - wbdf2dy['Qmod'][0], color='tab:green', alpha=0.7,label=r'ET$_{mod}$')
ax2.bar('2014', - wbdf2dy['Qmod'][1], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
ax2.scatter('2014', - wbdfy['Qobs'][1], s=200, marker='x', color='k', linewidth=4, zorder=2, label=r'$Q_{obs}$')
#ax2.scatter('2014', - wbdf2dy['ETmod'][1] - wbdf2dy['Qmod'][1] + wbdf2dy['S'][1] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2015', wbdf2dy['P'][2], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2015', - wbdf2dy['ETmod'][2] - wbdf2dy['Qmod'][1], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2015', - wbdf2dy['Qmod'][2], alpha=0.7, color='tab:brown', label='Q')
ax2.scatter('2015', - wbdf2dy['Qobs'][2], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2015', - wbdf2dy['ETmod'][2] - wbdf2dy['Qmod'][2] + wbdf2dy['S'][2] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2016', wbdf2dy['P'][3], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2016', - wbdf2dy['ETmod'][3] - wbdf2dy['Qmod'][2], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2016', - wbdf2dy['Qmod'][3], alpha=0.7, color='tab:brown', label='Q')
ax2.scatter('2016', - wbdf2dy['Qobs'][3], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2016', - wbdf2dy['ETmod'][3] - wbdf2dy['Qmod'][3] + wbdf2dy['S'][3] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2017', wbdf2dy['P'][4], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2017', - wbdf2dy['ETmod'][4] - wbdf2dy['Qmod'][3], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2017', - wbdf2dy['Qmod'][4], alpha=0.7, color='tab:brown', label='Q')
#ax2.scatter('2017', - wbdf2dy['Qobs'][4], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2017', - wbdf2dy['ETmod'][4] - wbdf2dy['Qmod'][4] + wbdf2dy['S'][4] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2018', wbdf2dy['P'][5], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2018', - wbdf2dy['ETmod'][5] - wbdf2dy['Qmod'][4], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2018', - wbdf2dy['Qmod'][5], alpha=0.7, color='tab:brown', label='Q')
ax2.scatter('2018', - wbdf2dy['Qobs'][5], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2018', - wbdf2dy['ETmod'][5] - wbdf2dy['Qmod'][5] + wbdf2dy['S'][5] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2019', wbdf2dy['P'][6], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2019', - wbdf2dy['ETmod'][6] - wbdf2dy['Qmod'][5], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2019', - wbdf2dy['Qmod'][6], alpha=0.7, color='tab:brown', label='Q')
ax2.scatter('2019', - wbdf2dy['Qobs'][6], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2019', - wbdf2dy['ETmod'][6] - wbdf2dy['Qmod'][6] + wbdf2dy['S'][6] , s=200, marker='_', color='blue', linewidth=1, zorder=2)

ax2.bar('2020', wbdf2dy['P'][7], alpha=0.7, color='tab:blue', label='P')
ax2.bar('2020', - wbdf2dy['ETmod'][7] - wbdf2dy['Qmod'][6], color='tab:green', alpha=0.7, label='ET')
ax2.bar('2020', - wbdf2dy['Qmod'][7], alpha=0.7, color='tab:brown', label='Q')
#ax2.scatter('2020', - wbdf2dy['Qobs'][7], s=200, marker='x', color='k', linewidth=4, zorder=2)
#ax2.scatter('2020', - wbdf2dy['ETmod'][7] - wbdf2dy['Qmod'][7] + wbdf2dy['S'][7] , s=200, marker='_', color='blue', linewidth=1, zorder=2)
ax1.grid(axis='y'); ax2.grid(axis='y')
ax1.set_ylabel('mm / year')
ax1.set_title('SpaFHy-catch')
ax2.set_title('SpaFHy-2D')

if saveplots == True:
        plt.savefig(f'WB_BARPLOT_{today}.pdf')
        plt.savefig(f'WB_BARPLOT_{today}.png')