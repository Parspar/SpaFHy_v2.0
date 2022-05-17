# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:38:00 2022

@author: janousu
"""

from iotools import read_results
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from datetime import date
import os
from netCDF4 import Dataset
import pandas as pd
from iotools import read_AsciiGrid, write_AsciiGrid
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from sklearn.metrics import mean_absolute_error as mae
from rasterio.plot import show
import rasterio
import xarray as xr
from raster_utils import read_pkrasteri_for_extent

#%%

## reading simulation results
# reading the 1D results
outputfile_1d = r'D:\SpaFHy_2D_2021\results_1d.nc'
res_1d = xr.open_dataset(outputfile_1d)

# reading the 2D results
outputfile_2d = r'D:\SpaFHy_2D_2021\results_2d.nc'
res_2d = xr.open_dataset(outputfile_2d)

# reading the top results
outputfile_top = r'D:\SpaFHy_2D_2021\results_top.nc'
res_top = xr.open_dataset(outputfile_top)

dates_spa = []
for d in range(len(res_1d['time'])):
    dates_spa.append(pd.to_datetime(str(res_1d['time'][d])[36:46]))
#dates_spa = pd.to_datetime(dates_spa)

# reading basic map
pkfp = 'C:\SpaFHy_v1_Pallas_2D/testcase_input/parameters/pkmosaic_clipped.tif'
bbox = [res_1d['lon'].min(), res_1d['lon'].max(), res_1d['lat'].min(), res_1d['lat'].max()]
pk, meta = read_pkrasteri_for_extent(pkfp, bbox, showfig=False)

# reading SAR files
sar_tempfile = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_TEMPINTERP.nc'
sar_temp = xr.open_dataset(sar_tempfile)

sar_spatfile = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_TEMPINTERP_SPATINTERP.nc'
sar_spat = xr.open_dataset(sar_spatfile)

#%%

# simple tests for simulation files
try:
    res_1d['time'] in res_2d['time'] in res_top['time']
    print('Simulation times are a match!')
except ValueError:
    print('Simulation times DO NOT match!')
try:
    res_1d['lat'] in res_2d['lat'] in res_top['lat'] and res_1d['lon'] in res_2d['lon'] in res_top['lon']
    print('Simulation extents are a match!')
except ValueError:
    print('Simulation extents DO NOT match!')

# simple tests for SAR files
try:
    sar_temp['time'] in sar_spat['time']
    print('SAR dates are a match!')
except ValueError:
    print('SAR dates do not match!')
try:
    sar_temp['lat'] in sar_spat['lat'] and sar_temp['lon'] in sar_spat['lon']
    print('SAR extents are a match!')
except ValueError:
    print('SAR extents DO NOT match!')

# simple tests for SAR + SpaFHy files
try:
    sar_temp['lat'] in sar_spat['lat'] in res_1d['lat'] and sar_temp['lon'] in sar_spat['lon'] in res_1d['lon']
    print('SAR and SpaFHy extents are a match!')
except ValueError:
    print('SAR and SpaFHy extents DO NOT match!')
try:
    sar_temp['time'] in sar_spat['time'] in res_1d['time']
    print('SAR and SpaFHy times are a match!')
except ValueError:
    print('SAR and SpaFHy times DO NOT match!')

#%%
import geopandas as gpd
# READING POINT OBSERVATIONS

# 1 discharge
fn1 = r'C:\SpaFHy_v1_Pallas_2D\obs\Runoffs1d_SVEcatchments_mmd.csv'
q = pd.read_csv(fn1, sep=';', index_col=0, parse_dates=True)
q = q.loc[q.index.isin(dates_spa)]

# 2 swe
fn2 = r'C:\SpaFHy_v1_Pallas_2D\obs\swe_mm.csv'
swe = pd.read_csv(fn2, sep=';', index_col=0, parse_dates=True)

# 3 kenttärova soil moisture
fn3 = r'C:\SpaFHy_v1_Pallas_2D\obs\theta_kenttarova.csv'
theta = pd.read_csv(fn3, sep=';', index_col=0, parse_dates=True)

# 4 groundwater levels
fn4 = r'C:\SpaFHy_v1_Pallas_2D\obs\gw_levels.csv'
gw = pd.read_csv(fn4, sep=';', index_col=0, parse_dates=True)

# 5 evapotranspiration
fn5 = r'C:\SpaFHy_v1_Pallas_2D\obs\ec_et.csv'
et = pd.read_csv(fn5, sep=';', index_col=0, parse_dates=True)

# 6 spatial soil moisture
fn6 = r'C:\SpaFHy_v1_Pallas_2D\obs\theta_spatial.nc'
theta_spat = xr.open_dataset(fn6)

# 7 spatial soil moisture geopandas
fn7 = r'C:\Users\janousu\OneDrive - Oulun yliopisto\Pallas data\DATABASE\SOIL MOISTURE AND TEMPERATURE\PROCESSING\Spatial\final_catchment_SM_mean_gpd_JP.csv'
theta_spat_gpd = pd.read_csv(fn7, sep=';', index_col=0, parse_dates=True)
theta_spat_gpd = gpd.GeoDataFrame(theta_spat_gpd)
theta_spat_gpd['geometry'] = gpd.GeoSeries.from_wkt(theta_spat_gpd['geometry'])

#%%

# parameters
today = date.today()
saveplots = True

# indexes for tighet plots
zx = np.arange(20, 171, 1)
zy = np.arange(20, 246, 1)

# defining important raster locations
ht = [118, 136] # hilltop
om = [46, 54]   # open mire


#%%

# creating pd dataframes for those sim vs. obs we want scatterplots of
# soil moisture at hilltop and open mire

temporaldf = pd.DataFrame()

temporaldf['2d_bucket_moisture_root_ht'] = res_2d['bucket_moisture_root'][:,ht[0],ht[1]]
temporaldf['top_bucket_moisture_root_ht'] = res_top['bucket_moisture_root'][:,ht[0],ht[1]]

temporaldf['2d_bucket_moisture_root_om'] = res_2d['bucket_moisture_root'][:,om[0],om[1]]
temporaldf['top_bucket_moisture_root_om'] = res_top['bucket_moisture_root'][:,om[0],om[1]]

temporaldf.index = dates_spa

temporaldf['obs_mean_moisture_root_ht'] = theta['mean_obs']
temporaldf['obs_min_moisture_root_ht'] = theta['min']
temporaldf['obs_max_moisture_root_ht'] = theta['max']

temporaldf['obs_swe_mean_ca'] = np.nanmean(res_1d['canopy_snow_water_equivalent'], axis=(1,2))


#%%
# because we do not want to assess winter soil moisture (no soil freezing modelled)
temporaldf.loc[temporaldf['obs_swe_mean_ca']  > 0, '2d_bucket_moisture_root_ht'] = np.nan
temporaldf.loc[temporaldf['obs_swe_mean_ca']  > 0, 'top_bucket_moisture_root_ht'] = np.nan
temporaldf.loc[temporaldf['obs_swe_mean_ca']  > 0, 'obs_mean_moisture_root_ht'] = np.nan
temporaldf.loc[temporaldf['obs_swe_mean_ca']  > 0, 'obs_min_moisture_root_ht'] = np.nan
temporaldf.loc[temporaldf['obs_swe_mean_ca']  > 0, 'obs_max_moisture_root_ht'] = np.nan

#%%

os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

# GIS PLOT

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(12,9));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[0][3]

ax5 = axs[1][0]
ax6 = axs[1][1]
ax7 = axs[1][2]
ax8 = axs[1][3]

ax9 = axs[2][0]
ax10 = axs[2][1]
ax11 = axs[2][2]
ax12 = axs[2][3]

res_2d['parameters_lai_conif'].plot(ax=ax1)


im1 = ax1.imshow(res_2d['parameters_lai_conif'][20:250,20:165])
ax1.set_title('LAI conif')
fig.colorbar(im1, ax=ax1)

im2 = ax2.imshow(res_2d['parameters_lai_decid_max'][20:250,20:165])
ax2.set_title('LAI decid max')
fig.colorbar(im2, ax=ax2)

im3 = ax3.imshow(res_2d['parameters_lai_shrub'][20:250,20:165])
ax3.set_title('LAI shrub')
fig.colorbar(im3, ax=ax3)

im4 = ax4.imshow(res_2d['parameters_lai_grass'][20:250,20:165])
ax4.set_title('LAI grass')
fig.colorbar(im4, ax=ax4)

im5 = ax5.imshow(res_2d['parameters_hc'][20:250,20:165])
ax5.set_title('canopy height')
fig.colorbar(im5, ax=ax5)

im6 = ax6.imshow(res_2d['parameters_cf'][20:250,20:165], label='canopy fraction')
ax6.set_title('canopy fraction')
fig.colorbar(im6, ax=ax6)

cmapsoil = plt.get_cmap('viridis', 4)
im7 = ax7.imshow(res_top['parameters_soilclass'][20:250,20:165], cmap=cmapsoil)
ax7.set_title('soilclass')
fig.colorbar(im7, ax=ax7)

cmapsite = plt.get_cmap('viridis', 4)
im8 = ax8.imshow(res_top['parameters_sitetype'][20:250,20:165], cmap=cmapsite)
ax8.set_title('sitetype')
fig.colorbar(im8, ax=ax8)

cmapditch = plt.get_cmap('viridis', 2)
im9 = ax9.imshow(res_top['parameters_ditches'][20:250,20:165], cmap=cmapditch)
ax9.set_title('streams/ditches')
cbar = fig.colorbar(im9, ax=ax9)
cbar.ax.locator_params(nbins=1)

im10 = ax10.imshow(res_top['parameters_elevation'][20:250,20:165])
ax10.set_title('elevation')
fig.colorbar(im10, ax=ax10)

im11 = ax11.imshow(res_top['parameters_twi'][20:250,20:165])
ax11.set_title('TWI')
fig.colorbar(im11, ax=ax11)

#im12 = ax12.imshow(results_catch['parameters_twi'][20:250,20:165])
#ax12.set_title('shading coefficient')
#fig.colorbar(im12, ax=ax12)

ax1.axis('off'); ax2.axis('off'); ax3.axis('off'); ax4.axis('off')
ax5.axis('off'); ax6.axis('off'); ax7.axis('off'); ax8.axis('off')
ax9.axis('off'); ax10.axis('off'); ax11.axis('off'); ax12.axis('off')

#plt.tight_layout()
#ax10.imshow(results_2d['parameters_twi'][20:250,20:165])

if saveplots == True:
        plt.savefig(f'GIS_rasters_{today}.pdf',bbox_inches='tight')
        plt.savefig(f'GIS_rasters_{today}.png',bbox_inches='tight')


#%%

# TEMPORAL SOIL MOISTURE ET HILLTOP (KENTTÄROVA)
# Plotting

period = pd.date_range('2014-05-01', '2021-10-01', freq='D')

ylims = [0.1,0.47]
###
fig = plt.figure(constrained_layout=True, figsize=(14,8))
gs = fig.add_gridspec(2, 4)
sns.set_style('whitegrid')

ax1 = fig.add_subplot(gs[0, :3])
ax1.set_title('Hilltop')
#im1 = temporaldf['2d_bucket_moisture_root_ht'].plot(ax=ax1, color='g', alpha=0.7, label='2D')
temporaldf['top_bucket_moisture_root_ht'].plot(ax=ax1, color='g', alpha=0.6, label='TOP')
temporaldf['obs_mean_moisture_root_ht'].plot(color='k', alpha=0.5, label=r'mean$_{obs}$')
ax1.fill_between(temporaldf.index, temporaldf['obs_min_moisture_root_ht'], temporaldf['obs_max_moisture_root_ht'],
                 color='blue', alpha=0.2, label=r'min/max$_{obs}$')
ax1.legend(ncol=5)

y = ax1.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
ax1.set_ylim(ylims)
#ax1.axes.get_xaxis().set_ticklabels([])

ax2 = fig.add_subplot(gs[0, 3])
im2 = sns.regplot(ax=ax2, x=temporaldf['top_bucket_moisture_root_ht'], y=temporaldf['obs_mean_moisture_root_ht'], scatter_kws={'s':50, 'alpha':0.2}, line_kws={"color": "red"})
ax2.set(ylim=(ylims))
ax2.set(xlim=(ylims))
ax2.yaxis.tick_right()
ax2.set_ylabel(r'$\theta_{obs}$ (m$^3$m$^{-3}$)')
ax2.set_xlabel(r'$\theta_{2D}$ (m$^3$m$^{-3}$)')
ax2.plot(ylims, ylims, 'k--')

'''
ax3 = fig.add_subplot(gs[1, :3])
ax3.set_title('Open mire')
im3 = temporaldf['2d_bucket_moisture_root_om'].plot(ax=ax3, color='g', alpha=0.7, label='2D')
temporaldf['top_bucket_moisture_root_om'].plot(ax=ax3, color='red', alpha=0.6, label='TOP')
ax1.legend(ncol=5)
ax3.set_ylim(0.35,0.9)
ax3.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#f3_ax3.legend(ncol=5)
'''

if saveplots == True:
        plt.savefig(f'theta_model_ts_{today}.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(f'theta_model_ts_{today}.png', bbox_inches='tight', dpi=300)


#%%

start = np.where(pd.to_datetime(dates_spa) == '2021-05-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2021-09-02')[0][0]

d = '2021-06-17'
doi = np.where(pd.to_datetime(dates_spa[start:end]) == d)[0][0]
doi_m = np.where(pd.to_datetime(theta_spat['time'][:].data) == d)[0][0]

alp=0.5
ylims = [0.2,0.88]
fig = plt.figure(figsize=(8,13))
gs = fig.add_gridspec(5, 2)

ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[1:4, 0:2])

ax1.plot(dates_spa[start:end], np.nanmean(res_top['bucket_moisture_root'][start:end], axis=(1,2)), linewidth=2, color='black', label='mean')
ax1.fill_between(dates_spa[start:end], np.nanquantile(res_top['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
                      np.nanquantile(res_top['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-quantile')
ax1.scatter(dates_spa[start:end][doi], np.nanmean(res_top['bucket_moisture_root'][start:end][doi]), s=70, color='blue')
ax1.legend()
ax1.set_title('Timeseries 2021')
ax1.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
ax1.grid()
ax1.grid()

ax1.xaxis.set_tick_params(rotation=20)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
res_top['bucket_moisture_root'][start:end][doi].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1])

plt.subplots_adjust(wspace=0.05, hspace=0.7)


plt.savefig(f'spatial_06_2021_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'spatial_06_2021_{today}.png', bbox_inches='tight', dpi=300)

#%%

d = '2021-09-01'
doi = np.where(pd.to_datetime(dates_spa[start:end]) == d)[0][0]
doi_m = np.where(pd.to_datetime(theta_spat['time'][:].data) == d)[0][0]

alp=0.5
ylims = [0.2,0.88]
fig = plt.figure(figsize=(14,8))
gs = fig.add_gridspec(1, 3)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

rasterio.plot.show(pk, transform=meta['transform'], ax=ax1);
im1 = res_1d['bucket_moisture_root'][start:end][doi].plot(ax=ax1, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax1, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = res_top['bucket_moisture_root'][start:end][doi].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1])

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3)

res_2d['bucket_moisture_root'][start:end][doi].plot(ax=ax3, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                    add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax3, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1])
                                               #theta_spat_gpd.loc[theta_spat_gpd.index == d][['SM_mean', 'geometry']].plot()



ax2.axes.get_yaxis().set_ticklabels([])
ax2.set_ylabel('')
ax3.axes.get_yaxis().set_ticklabels([])
ax3.set_ylabel('')

cax = plt.axes([0.1, 0.14, 0.8, 0.03])
cbar1 = plt.colorbar(im1, cax=cax, orientation='horizontal')

plt.subplots_adjust(wspace=0.1, hspace=0)


plt.savefig(f'spatial_mode_{d}_2021_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'spatial_mode_{d}_2021_{today}.png', bbox_inches='tight', dpi=300)

#%%

june_id = np.where(pd.to_datetime(dates_spa).month == 6)
july_id = np.where(pd.to_datetime(dates_spa).month == 7)
august_id = np.where(pd.to_datetime(dates_spa).month == 8)
sept_id = np.where(pd.to_datetime(dates_spa).month == 9)

# HISTOGRAMS

fig = plt.figure(figsize=(15,6))
gs = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 2])
ax8 = fig.add_subplot(gs[1, 3])

res_top['bucket_moisture_root'][june_id].plot.hist(ax=ax1, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][june_id].plot.hist(ax=ax1, bins=25, alpha=0.7, label='2D')
ax1.legend()

res_top['bucket_moisture_root'][july_id].plot.hist(ax=ax2, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][july_id].plot.hist(ax=ax2, bins=25, alpha=0.7, label='2D')

res_top['bucket_moisture_root'][august_id].plot.hist(ax=ax3, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][august_id].plot.hist(ax=ax3, bins=25, alpha=0.7, label='2D')

res_top['bucket_moisture_root'][august_id].plot.hist(ax=ax4, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][august_id].plot.hist(ax=ax4, bins=25, alpha=0.7, label='2D')

###

res_top['bucket_moisture_root'][june_id].plot.hist(ax=ax5, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][june_id].plot.hist(ax=ax5, bins=25, alpha=0.7, label='2D')
ax1.legend()

res_top['bucket_moisture_root'][july_id].plot.hist(ax=ax6, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][july_id].plot.hist(ax=ax6, bins=25, alpha=0.7, label='2D')

res_top['bucket_moisture_root'][august_id].plot.hist(ax=ax7, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][august_id].plot.hist(ax=ax7, bins=25, alpha=0.7, label='2D')

res_top['bucket_moisture_root'][august_id].plot.hist(ax=ax8, bins=25, alpha=0.7, label='TOP')
res_2d['bucket_moisture_root'][august_id].plot.hist(ax=ax8, bins=25, alpha=0.7, label='2D')