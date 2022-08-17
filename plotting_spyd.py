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

os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')

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
res_top['total_evapotranspiration'] = res_top['bucket_evaporation'] + res_top['canopy_transpiration'] + res_top['canopy_evaporation']

# raster parameters
soilclass = np.array(res_2d['parameters_soilclass'])
cmask = np.array(res_2d['parameters_cmask'])

dates_spa = []
for d in range(len(res_1d['time'])):
    dates_spa.append(pd.to_datetime(str(res_1d['time'][d])[36:46]))
dates_spa = pd.to_datetime(dates_spa)

# reading basic map
pkfp = 'C:\SpaFHy_v1_Pallas_2D/testcase_input/parameters/pkmosaic_clipped.tif'
bbox = [res_1d['lon'].min(), res_1d['lon'].max(), res_1d['lat'].min(), res_1d['lat'].max()]
pk, meta = read_pkrasteri_for_extent(pkfp, bbox, showfig=False)

# reading SAR files
sar_tempfile = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_TEMPINTERP.nc'
sar_temp = xr.open_dataset(sar_tempfile)
sar_temp = sar_temp * cmask

sar_spatfile = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_TEMPINTERP_SPATINTERP.nc'
sar_spat = xr.open_dataset(sar_spatfile)
sar_spat = sar_spat * cmask

dates_sar = []
for d in range(len(sar_temp['time'])):
    dates_sar.append(pd.to_datetime(str(sar_temp['time'][d])[36:46]))
dates_sar = pd.to_datetime(dates_sar)

sar_mask = np.array(sar_temp['theta'][0].copy())
sar_mask[sar_mask >= 0] = 1
sar_mask = sar_mask.reshape(sar_temp['theta'][0].shape) * cmask

# parameters
today = date.today()
saveplots = True

# distributed radiation
distradfile = r'C:\SpaFHy_v1_Pallas_2D/obs/rad_ds.nc'
distrad = xr.open_dataset(distradfile)

# indexes for tighet plots
zx = np.arange(20, 171, 1)
zy = np.arange(20, 246, 1)

# defining important raster locations
ht = [118,136] # hilltop
om = [46, 54]   # open mire

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

alp=0.75


# GIS PLOT

# Plotting
fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(12,10));
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



#im1 = ax1.imshow(res_2d['parameters_lai_conif'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax1);
im1 = res_2d['parameters_lai_conif'][20:250,20:165].plot(ax=ax1, alpha=alp, add_colorbar=False)
ax1.set_title('LAI conif')
fig.colorbar(im1, ax=ax1)

#im2 = ax2.imshow(res_2d['parameters_lai_decid_max'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
im2 = res_2d['parameters_lai_decid_max'][20:250,20:165].plot(ax=ax2, alpha=alp, add_colorbar=False)
ax2.set_title('LAI decid')
fig.colorbar(im2, ax=ax2)

#im3 = ax3.imshow(res_2d['parameters_lai_shrub'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
im3 = res_2d['parameters_lai_shrub'][20:250,20:165].plot(ax=ax3, alpha=alp, add_colorbar=False)
ax3.set_title('LAI shrub')
fig.colorbar(im3, ax=ax3)

#im4 = ax4.imshow(res_2d['parameters_lai_grass'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax4);
im4 = res_2d['parameters_lai_grass'][20:250,20:165].plot(ax=ax4, alpha=alp, add_colorbar=False)
ax4.set_title('LAI grass')
fig.colorbar(im4, ax=ax4)

#im5 = ax5.imshow(res_2d['parameters_hc'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax5);
im5 = res_2d['parameters_hc'][20:250,20:165].plot(ax=ax5, alpha=alp, add_colorbar=False)
ax5.set_title('canopy height')
fig.colorbar(im5, ax=ax5)

#im6 = ax6.imshow(res_2d['parameters_cf'][20:250,20:165], label='canopy fraction')
rasterio.plot.show(pk, transform=meta['transform'], ax=ax6);
im6 = res_2d['parameters_cf'][20:250,20:165].plot(ax=ax6, alpha=alp, add_colorbar=False)
ax6.set_title('canopy fraction')
fig.colorbar(im6, ax=ax6)

cmapsoil = plt.get_cmap('viridis', 4)
#im7 = ax7.imshow(res_2d['parameters_soilclass'][20:250,20:165], cmap=cmapsoil)
rasterio.plot.show(pk, transform=meta['transform'], ax=ax7);
im7 = res_2d['parameters_soilclass'][20:250,20:165].plot(ax=ax7, alpha=alp, cmap=cmapsoil, add_colorbar=False)
ax7.set_title('soilclass')
cbar = fig.colorbar(im7, ax=ax7)
cbar.ax.locator_params(nbins=4)

cmapsite = plt.get_cmap('viridis', 4)
#im8 = ax8.imshow(res_2d['parameters_sitetype'][20:250,20:165], cmap=cmapsite)
rasterio.plot.show(pk, transform=meta['transform'], ax=ax8);
im8 = res_2d['parameters_sitetype'][20:250,20:165].plot(ax=ax8, alpha=alp, cmap=cmapsite, add_colorbar=False)
ax8.set_title('sitetype')
cbar = fig.colorbar(im8, ax=ax8)
cbar.ax.locator_params(nbins=4)

cmapditch = plt.get_cmap('viridis', 2)
#im9 = ax9.imshow(res_2d['parameters_ditches'][20:250,20:165], cmap=cmapditch)
rasterio.plot.show(pk, transform=meta['transform'], ax=ax9);
im9 = res_2d['parameters_ditches'][20:250,20:165].plot(ax=ax9, alpha=alp, cmap=cmapditch, add_colorbar=False)
ax9.set_title('streams/ditches')
cbar = fig.colorbar(im9, ax=ax9)
cbar.ax.locator_params(nbins=1)

# = ax10.imshow(res_2d['parameters_elevation'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax10);
im10 = res_2d['parameters_elevation'][20:250,20:165].plot(ax=ax10, alpha=alp, add_colorbar=False)
ax10.set_title('elevation')
fig.colorbar(im10, ax=ax10)

#im11 = ax11.imshow(res_top['parameters_twi'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax11);
im11 = res_top['parameters_twi'][20:250,20:165].plot(ax=ax11, alpha=alp, add_colorbar=False)
ax11.set_title('TWI')
fig.colorbar(im11, ax=ax11)

#im12 = ax12.imshow(results_catch['parameters_twi'][20:250,20:165])
rasterio.plot.show(pk, transform=meta['transform'], ax=ax12);
im12 = (distrad['c_rad'].mean(dim='time')*cmask)[20:250,20:165].plot(ax=ax12, alpha=alp, add_colorbar=False)
ax12.set_title('shading coefficient')
fig.colorbar(im12, ax=ax12)

ax1.axis('off'); ax2.axis('off'); ax3.axis('off'); ax4.axis('off')
ax5.axis('off'); ax6.axis('off'); ax7.axis('off'); ax8.axis('off')
ax9.axis('off'); ax10.axis('off'); ax11.axis('off'); ax12.axis('off')

#plt.tight_layout()
#ax10.imshow(results_2d['parameters_twi'][20:250,20:165])

plt.subplots_adjust(wspace=0, hspace=0.15)

if saveplots == True:
        plt.savefig(f'GIS_rasters_{today}.pdf',bbox_inches='tight', dpi=300)
        plt.savefig(f'GIS_rasters_{today}.png',bbox_inches='tight', dpi=300)

#%%

alp=0.75

# Plotting
fig = plt.figure(constrained_layout=True, figsize=(16,6))
gs = fig.add_gridspec(2, 5)

ax0 = fig.add_subplot(gs[0:, 0:2])
ax1 = fig.add_subplot(gs[0, 2])
ax2 = fig.add_subplot(gs[0, 3])
ax3 = fig.add_subplot(gs[0, 4])
ax4 = fig.add_subplot(gs[1, 2])
ax5 = fig.add_subplot(gs[1, 3])
ax6 = fig.add_subplot(gs[1, 4])


rasterio.plot.show(pk, transform=meta['transform'], ax=ax0)
im0 = res_2d['parameters_cmask'].plot(ax=ax0, alpha=0.2, add_colorbar=False)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax1)
im1 = (res_2d['parameters_lai_conif']+res_2d['parameters_lai_decid_max']).plot(ax=ax1, alpha=alp, add_colorbar=False)
fig.colorbar(im1, ax=ax1)
# this creates colorbar
ax1.axes.get_xaxis().set_ticklabels([])
ax1.axes.get_yaxis().set_ticklabels([])
ax1.axis('off')
ax1.set_title(r'LAI conif [m$^2$/m$^2$]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
im2 = res_2d['parameters_hc'].plot(ax=ax2, alpha=alp, add_colorbar=False)
fig.colorbar(im2, ax=ax2)
ax2.set_title('Canopy height [m]')
ax2.axes.get_xaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticklabels([])
ax2.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
im3 = res_2d['parameters_cf'].plot(ax=ax3, alpha=alp, add_colorbar=False)
fig.colorbar(im3, ax=ax3)
ax3.set_title('Canopy fraction [m$^2$/m$^2$]')
ax3.axes.get_xaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])
ax3.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax4);
im4 = res_2d['parameters_soilclass'].plot(ax=ax4, alpha=alp, add_colorbar=False)
fig.colorbar(im4, ax=ax4)
ax4.set_title('Soilclass [-]')
ax4.axes.get_xaxis().set_ticklabels([])
ax4.axes.get_yaxis().set_ticklabels([])
ax4.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax5);
im5 = res_2d['parameters_elevation'].plot(ax=ax5, alpha=alp, add_colorbar=False)
fig.colorbar(im5, ax=ax5)
ax5.axes.get_yaxis().set_ticks([])
ax5.axes.get_xaxis().set_ticks([])
ax5.set_title('Elevation [m]')
ax5.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax6);
im6 = res_top['parameters_twi'].plot(ax=ax6, alpha=alp, add_colorbar=False)
fig.colorbar(im6, ax=ax6)
ax6.axes.get_yaxis().set_ticks([])
ax6.axes.get_xaxis().set_ticks([])
ax6.set_title('TWI [-]')
ax6.axis('off')

#plt.subplots_adjust(wspace=0.3, hspace=0.1)

plt.savefig(f'EGU_gis_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'EGU_gis_{today}.png', bbox_inches='tight', dpi=300)

#%%
import matplotlib as mpl

# change working dir
os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')
alp=0.75
shr=0.9
# show raster overlays
plt.close('all')

# Plotting
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(11,8));
ax1 = axs[0][0]
ax2 = axs[0][1]
ax3 = axs[0][2]
ax4 = axs[1][0]
ax5 = axs[1][1]
ax6 = axs[1][2]


rasterio.plot.show(pk, transform=meta['transform'], ax=ax1)
im1 = (res_2d['parameters_lai_conif']+res_2d['parameters_lai_decid_max']).plot(ax=ax1, alpha=alp, add_colorbar=False)
fig.colorbar(im1, ax=ax1)
# this creates colorbar
ax1.axes.get_xaxis().set_ticklabels([])
ax1.axes.get_yaxis().set_ticklabels([])
ax1.axis('off')
ax1.set_title(r'Lehtialaindeksi [m$^2$/m$^2$]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
im2 = res_2d['parameters_hc'].plot(ax=ax2, alpha=alp, add_colorbar=False)
fig.colorbar(im2, ax=ax2)
ax2.axes.get_xaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticklabels([])
ax2.axis('off')
ax2.set_title('Kasvuston korkeus [m]')


rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
im3 = res_2d['parameters_cf'].plot(ax=ax3, alpha=alp, add_colorbar=False)
fig.colorbar(im3, ax=ax3)
ax3.axes.get_xaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])
ax3.axis('off')
ax3.set_title('Kasvuston osuus [m$^2$/m$^2$]')




rasterio.plot.show(pk, transform=meta['transform'], ax=ax4);
im4 = res_2d['parameters_soilclass'].plot(ax=ax4, alpha=alp, add_colorbar=False)
cbar = fig.colorbar(im4, ax=ax4)
ax4.set_title('Pintamaaluokka [-]')
ax4.axes.get_xaxis().set_ticklabels([])
ax4.axes.get_yaxis().set_ticklabels([])
cbar.set_ticks([1, 2, 3, 4])
#cbar.ax.set_yticklabels(['karkea', 'keski', 'hieno', 'turve'])

ax4.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax5);
im5 = res_2d['parameters_elevation'].plot(ax=ax5, alpha=alp, add_colorbar=False)
fig.colorbar(im5, ax=ax5)
ax5.axes.get_yaxis().set_ticks([])
ax5.axes.get_xaxis().set_ticks([])
ax5.set_title('Korkeus [m]')
ax5.axis('off')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax6);
im6 = res_top['parameters_twi'].plot(ax=ax6, alpha=alp, add_colorbar=False)
fig.colorbar(im6, ax=ax6)
ax6.axes.get_yaxis().set_ticks([])
ax6.axes.get_xaxis().set_ticks([])
ax6.set_title('Kosteusindeksi [-]')
ax6.axis('off')

plt.subplots_adjust(wspace=0.05, hspace=0.05)

plt.savefig(f'vesitalous_gis_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'vesitalous_gis_{today}.png', bbox_inches='tight', dpi=300)


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

## ONE SPAFHY SOIL MOIST PLOT WITH OBSERVATIONS IN JUNE 2021

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

#%%

start = np.where(pd.to_datetime(dates_spa) == '2020-05-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2020-10-01')[0][0]
maxd = np.where(np.nansum(res_top['bucket_moisture_root'][start:end], axis=(1,2))
                == np.nansum(res_top['bucket_moisture_root'][start:end], axis=(1,2)).max())[0][0]
mind = np.where(np.nansum(res_top['bucket_moisture_root'][start:end], axis=(1,2))
                == np.nansum(res_top['bucket_moisture_root'][start:end], axis=(1,2)).min())[0][0]
mind = np.where(pd.to_datetime(dates_spa[start:end]) == '2020-09-02')[0][0]

alp=0.8
ylims = [0.1, 0.88]

fig = plt.figure(figsize=(8,8))
gs = fig.add_gridspec(5, 2)

ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[1:4, 0])
ax3 = fig.add_subplot(gs[1:4, 1])

ax1.plot(dates_spa[start:end], np.nanmean(res_top['bucket_moisture_root'][start:end], axis=(1,2)), linewidth=2, color='black', label='keskiarvo')
ax1.fill_between(dates_spa[start:end], np.nanquantile(res_top['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
                      np.nanquantile(res_top['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-kvantiili')
ax1.scatter(dates_spa[start:end][maxd], np.nanmean(res_top['bucket_moisture_root'][start:end][maxd]), s=70, color='blue')
ax1.scatter(dates_spa[start:end][mind], np.nanmean(res_top['bucket_moisture_root'][start:end][mind]), s=70, color='red')
ax1.legend()
ax1.set_title('Maankosteuden aikasarja')
ax1.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
ax1.grid()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = res_top['bucket_moisture_root'][start:end][maxd].plot(ax=ax2, cmap='coolwarm_r',
                                                      vmin=ylims[0], vmax=ylims[1], alpha=alp, add_colorbar=False)
#im1 = rasterio.plot.show(snap1, transform=snap1.transform, ax=ax2, cmap='coolwarm_r', vmin=0.1, vmax=0.88, alpha=alp)
#im1b = im1.get_images()[1]
ax2.text(384500, 7542000, f'{str(dates_spa[start:end][maxd])[0:10]}')
ax2.scatter(384400, 7542100, color='blue')
ax2.set_title('Esim. kostea päivä')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
#im2 = ax3.imshow(results_catch['bucket_moisture_root'][start:end][mind], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im3 = res_top['bucket_moisture_root'][start:end][mind].plot(ax=ax3, cmap='coolwarm_r',
                                                      vmin=ylims[0], vmax=ylims[1], alpha=alp, add_colorbar=False)
#im2 = rasterio.plot.show(snap2, transform=snap2.transform, ax=ax3, cmap='coolwarm_r', vmin=0.1, vmax=0.88, alpha=alp)
#im2b = im1.get_images()[1]
#fig.colorbar(im2b, ax=ax3, shrink=0.5)
ax3.text(384500, 7542000, f'{str(dates_spa[start:end][mind])[0:10]}')
ax3.scatter(384400, 7542100, color='red')
cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
fig.colorbar(im3, cax=cbar_ax, label=r'$\theta$ [m$^3$m$^{-3}]$')
ax3.set_title('Esim. kuiva päivä')

ax3.axes.get_yaxis().set_ticklabels([])
ax3.set_ylabel('')

plt.subplots_adjust(wspace=0.05, hspace=0.7)


plt.savefig(f'vesitalous_kosteus_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'vesitalous_kosteus_{today}.png', bbox_inches='tight', dpi=300)

#%%


## ONE SPAFHY SOIL MOIST PLOT WITH OBSERVATIONS IN JUNE 2021

start = np.where(pd.to_datetime(dates_spa) == '2020-01-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2021-01-01')[0][0]

d = '2020-12-31'
doi = np.where(pd.to_datetime(dates_spa[start:end]) == d)[0][0]
#doi_m = np.where(pd.to_datetime(theta_spat['time'][:].data) == d)[0][0]

alp=0.65
ylims = [200,320]
fig = plt.figure(figsize=(6,11))
gs = fig.add_gridspec(5, 2)

ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[1:4, 0:2])

ax1.plot(dates_spa[start:end], np.cumsum(np.nanmean(res_top['total_evapotranspiration'][start:end], axis=(1,2))),
         linewidth=2, color='black', label='keskiarvo')
ax1.fill_between(dates_spa[start:end], np.cumsum(np.nanquantile(res_top['total_evapotranspiration'][start:end], 0.8, axis=(1,2))),
                      np.cumsum(np.nanquantile(res_top['total_evapotranspiration'][start:end], 0.2, axis=(1,2))), alpha=0.3, label='20-kvantiili')

ax1.legend()
ax1.set_title('Kumulatiivinen ET')
ax1.set_ylabel(r'ET [mm]')
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
#ax1.grid()
ax1.grid()
ax1.xaxis.set_tick_params(rotation=20)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = (np.cumsum(res_top['total_evapotranspiration'][start:end], axis=0)[doi]*cmask).plot(
    ax=ax2, cmap='bwr_r', vmin=ylims[0], vmax=ylims[1], alpha=alp, add_colorbar=False)
ax2.set_title(f'Kumulatiivinen summa')
cbar_ax = fig.add_axes([0.85, 0.3, 0.03, 0.4])
fig.colorbar(im2, cax=cbar_ax, label=r'ET [mm]')
#theta_spat_gpd.loc[theta_spat_gpd.index == d,
#                   ['theta', 'geometry']].plot(column='theta',
#                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
#                                               alpha=alp, vmin=ylims[0], vmax=ylims[1])

plt.subplots_adjust(wspace=0.05, hspace=0.7)


plt.savefig(f'et_cum_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'et_cum_{today}.png', bbox_inches='tight', dpi=300)


#%%

# SPAFHY SPATIAL MODEL VERSION COMPARISON

start = np.where(pd.to_datetime(dates_spa) == '2021-05-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2021-09-02')[0][0]

d = '2021-06-17'
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
                                                          add_colorbar=False, label='MOD')
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax1, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax1.legend()
rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = res_top['bucket_moisture_root'][start:end][doi].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax2.legend()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3)

res_2d['bucket_moisture_root'][start:end][doi].plot(ax=ax3, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                    add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax3, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax3.legend()
                                              #theta_spat_gpd.loc[theta_spat_gpd.index == d][['SM_mean', 'geometry']].plot()



ax2.axes.get_yaxis().set_ticklabels([])
ax2.set_ylabel('')
ax3.axes.get_yaxis().set_ticklabels([])
ax3.set_ylabel('')

cax = plt.axes([0.1, 0.14, 0.8, 0.03])
cbar1 = plt.colorbar(im1, cax=cax, orientation='horizontal')
cbar1.set_label('Volumetric water content [m3/m3]')

plt.subplots_adjust(wspace=0.1, hspace=0)


plt.savefig(f'spatial_mode_{d}_2021_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'spatial_mode_{d}_2021_{today}.png', bbox_inches='tight', dpi=300)

#%%

# HISTOGRAMS

june_id = np.where(pd.to_datetime(dates_spa).month == 6)[0]
june_id_sar = np.where(pd.to_datetime(dates_sar).month == 6)[0]
june_id_sarmask = np.where((dates_spa.isin(dates_sar)) & (dates_spa.month == 6))[0]

july_id = np.where(pd.to_datetime(dates_spa).month == 7)[0]
july_id_sar = np.where(pd.to_datetime(dates_sar).month == 7)[0]
july_id_sarmask = np.where((dates_spa.isin(dates_sar)) & (dates_spa.month == 7))[0]

august_id = np.where(pd.to_datetime(dates_spa).month == 8)[0]
august_id_sar = np.where(pd.to_datetime(dates_sar).month == 8)[0]
august_id_sarmask = np.where((dates_spa.isin(dates_sar)) & (dates_spa.month == 8))[0]

sept_id = np.where(pd.to_datetime(dates_spa).month == 9)[0]
sept_id_sar = np.where(pd.to_datetime(dates_sar).month == 9)[0]
sept_id_sarmask = np.where((dates_spa.isin(dates_sar)) & (dates_spa.month == 9))[0]


peat_id = np.where(np.ravel(soilclass == 4))[0]
mineral_id = np.where(np.ravel(soilclass == 2))[0]

top_moist_peat = np.array(res_top['bucket_moisture_root']).reshape(res_top['bucket_moisture_root'].shape[0], -1)[:,peat_id]
d2_moist_peat = np.array(res_2d['bucket_moisture_root']).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,peat_id]
d1_moist_peat = np.array(res_1d['bucket_moisture_root']).reshape(res_1d['bucket_moisture_root'].shape[0], -1)[:,peat_id]
sar_moist_peat = np.array(sar_temp['theta']).reshape(sar_temp['theta'].shape[0], -1)[:,peat_id]
d2_sarmask_moist_peat = np.array(res_2d['bucket_moisture_root']*sar_mask).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,peat_id]


top_moist_mineral = np.array(res_top['bucket_moisture_root']).reshape(res_top['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
d2_moist_mineral = np.array(res_2d['bucket_moisture_root']).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
d1_moist_mineral = np.array(res_1d['bucket_moisture_root']).reshape(res_1d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
sar_moist_mineral = np.array(sar_temp['theta']).reshape(sar_temp['theta'].shape[0], -1)[:,mineral_id]
d2_sarmask_moist_mineral = np.array(res_2d['bucket_moisture_root']*sar_mask).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]

alp=0.6
peat_lims = [0.1,0.9]
min_lims = [0.1,0.9]

def histogram_match(data1, data2, lims,  bins=25):
    hobs,binobs = np.histogram(data1,bins=25, range=lims)
    hsim,binsim = np.histogram(data2,bins=25, range=lims)
    hobs=np.float64(hobs)
    hsim=np.float64(hsim)
    minima = np.minimum(hsim, hobs)
    gamma = round(np.sum(minima)/np.sum(hobs),2)
    return str(gamma)

'''
hobs,binobs = np.histogram(top_moist_mineral[peat_id],bins=25, range=min_lims)
hsim,binsim = np.histogram(d2_moist_mineral[peat_id],bins=25, range=min_lims)
#convert int to float, critical conversion for the result
hobs=np.float64(hobs)
hsim=np.float64(hsim)
#find the overlapping of two histogram
minima = np.minimum(hsim, hobs)
#compute the fraction of intersection area to the observed histogram area, hist intersection/overlap index
gamma = np.sum(minima)/np.sum(hobs)
'''
#%%

# HISTOGRAM PLOT

alp=0.6
peat_lims = [0.1,0.9]
min_lims = [0.1,0.9]

fig = plt.figure(figsize=(12,4))
gs = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 2])
ax8 = fig.add_subplot(gs[1, 3])

ax1.hist(d2_moist_mineral[june_id].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax1.hist(d1_moist_mineral[june_id].flatten(), bins=25, range=min_lims, alpha=alp, label='1D')
ax1.text(0.60, 20000,  'HM:')
ax1.text(0.7, 20000,  histogram_match(d2_moist_mineral[june_id].flatten(),
         d1_moist_mineral[june_id].flatten(), bins=25, lims=min_lims))
ax1.legend()
ax1.set_title('JUNE')
ax1.set_ylabel('occurences')

ax2.hist(d2_moist_mineral[july_id].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax2.hist(d1_moist_mineral[july_id].flatten(), bins=25, range=min_lims, alpha=alp, label='1D')
ax2.set_title('JULY')
ax2.text(0.60, 20000,  'HM:')
ax2.text(0.7, 20000,  histogram_match(d2_moist_mineral[july_id].flatten(),
         d1_moist_mineral[july_id].flatten(), bins=25, lims=min_lims))

ax3.hist(d2_moist_mineral[august_id].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax3.hist(d1_moist_mineral[august_id].flatten(), bins=25, range=min_lims, alpha=alp, label='1D')
ax3.set_title('AUGUST')
ax3.text(0.60, 20000,  'HM:')
ax3.text(0.7, 20000,  histogram_match(d2_moist_mineral[august_id].flatten(),
         d1_moist_mineral[august_id].flatten(), bins=25, lims=min_lims))

ax4.hist(d2_moist_mineral[sept_id].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax4.hist(d1_moist_mineral[sept_id].flatten(), bins=25, range=min_lims, alpha=alp, label='1D')
ax4.set_title('SEPTEMBER')
ax4.text(0.60, 20000,  'HM:')
ax4.text(0.7, 20000,  histogram_match(d2_moist_mineral[sept_id].flatten(),
         d1_moist_mineral[sept_id].flatten(), bins=25, lims=min_lims))

###
ax5.hist(d2_moist_peat[june_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax5.hist(d1_moist_peat[june_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='1D')
ax5.set_ylabel('occurences')
ax5.text(0.1, 20000,  'HM:')
ax5.text(0.2, 20000,  histogram_match(d2_moist_peat[june_id].flatten(),
         d1_moist_peat[june_id].flatten(), bins=25, lims=min_lims))

ax6.hist(d2_moist_peat[july_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax6.hist(d1_moist_peat[july_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='1D')
ax6.text(0.1, 20000,  'HM:')
ax6.text(0.2, 20000,  histogram_match(d2_moist_peat[july_id].flatten(),
         d1_moist_peat[july_id].flatten(), bins=25, lims=min_lims))

ax7.hist(d2_moist_peat[august_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax7.hist(d1_moist_peat[august_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='1D')
ax7.text(0.1, 20000,  'HM:')
ax7.text(0.2, 20000,  histogram_match(d2_moist_peat[august_id].flatten(),
         d1_moist_peat[august_id].flatten(), bins=25, lims=min_lims))

ax8.hist(d2_moist_peat[sept_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax8.hist(d1_moist_peat[sept_id].flatten(), bins=25, range=peat_lims, alpha=alp, label='1D')
ax8.text(0.1, 20000,  'HM:')
ax8.text(0.2, 20000,  histogram_match(d2_moist_peat[sept_id].flatten(),
         d1_moist_peat[sept_id].flatten(), bins=25, lims=min_lims))

ax1.axes.get_yaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])
ax4.axes.get_yaxis().set_ticklabels([])

#ax1.axes.get_xaxis().set_ticklabels([])
#ax2.axes.get_xaxis().set_ticklabels([])
#ax3.axes.get_xaxis().set_ticklabels([])
#ax4.axes.get_xaxis().set_ticklabels([])

ax5.axes.get_yaxis().set_ticklabels([])
ax6.axes.get_yaxis().set_ticklabels([])
ax7.axes.get_yaxis().set_ticklabels([])
ax8.axes.get_yaxis().set_ticklabels([])

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()
ax7.grid()
ax8.grid()

plt.subplots_adjust(wspace=0.05, hspace=0.2)
fig.text(0.5, 0.03, 'Rootzone volumetric water content [m3/m3]', ha='center')

fig.text(0.92, 0.6, 'MINERAL SOIL', rotation='vertical', ha='center')
fig.text(0.92, 0.27, 'PEAT SOIL', rotation='vertical', ha='center')



plt.savefig(f'histogram_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'histogram_{today}.png', bbox_inches='tight', dpi=300)

#%%
from pyproj import CRS
import geopandas as gpd
from shapely.geometry import Point
import cartopy.crs as ccrs

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
pallas = [67.995, 24.224]
hyytiala = [61.850, 24.283]
lons = np.array([24.224,24.283])
lats = np.array([67.995, 61.850])

# Define an orthographic projection, centered in Finland! from: http://www.statsmapsnpix.com/2019/09/globe-projections-and-insets-in-qgis.html
ortho = CRS.from_proj4("+proj=ortho +lat_0=60.00 +lon_0=23.0000 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
ortho = ccrs.Orthographic(central_longitude=23, central_latitude=60)
geo = ccrs.Geodetic()

points = ortho.transform_points(geo, lons, lats)

world.to_crs(ortho).plot(color='tab:blue', linewidth=0.5, edgecolor='white')
plt.scatter(points[0][0],points[0][1], color='red', s=4)
plt.scatter(points[1][0],points[1][1], color='red', s=4)

plt.text(-1.3e6,0.9e6, 'Pallas', fontsize='x-small')
plt.text(-1.8e6,0.1e6, 'Hyytiälä', fontsize='x-small')

# Remove x and y axis
plt.axis('off')

plt.savefig(f'world_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'world_{today}.png', bbox_inches='tight', dpi=300)

#%%


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


# SAR SPATIAL PLOTS
start_day = pd.to_datetime('2019-05-01')
end_day = pd.to_datetime('2019-09-15')

start = np.where(pd.to_datetime(dates_spa) == start_day)[0][0]
end = np.where(pd.to_datetime(dates_spa) == end_day)[0][0]

sar_wet_day = np.where(np.nansum(sar_temp['theta'], axis=(1,2)) == np.nansum(sar_temp['theta'], axis=(1,2)).max())
sar_dry_day = np.where(np.nansum(sar_temp['theta'], axis=(1,2)) == np.nansum(sar_temp['theta'], axis=(1,2)).min())

wet_day = pd.to_datetime(sar_temp['time'][sar_wet_day].data[0])
dry_day = pd.to_datetime(sar_temp['time'][sar_dry_day].data[0])

doi = np.where(pd.to_datetime(dates_spa[start:end]) == wet_day)[0][0]
d = pd.to_datetime(wet_day)

bbox = [7.5452e6, 7.5440e6, 383500, 384600]
bbox_id = [np.where(res_1d['lat'] == find_nearest(res_1d['lat'],bbox[0]))[0][0],
           np.where(res_1d['lat'] == find_nearest(res_1d['lat'],bbox[1]))[0][0],
           np.where(res_1d['lon'] == find_nearest(res_1d['lon'],bbox[2]))[0][0],
           np.where(res_1d['lon'] == find_nearest(res_1d['lon'],bbox[3]))[0][0]]

alp=0.5
ylims = [0.2,0.88]
ylimstemp = [0.08,0.5]
fig = plt.figure(figsize=(12,8))
gs = fig.add_gridspec(3, 2)

ax0 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1:, 0])
ax1 = fig.add_subplot(gs[1:, 1])
#ax3 = fig.add_subplot(gs[1:, 2])

#ax0.plot(dates_spa[start:end], np.nanmean(res_2d['bucket_moisture_root'][start:end], axis=(1,2)), linewidth=2, color='black', label='mean')
ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_root'][start:end, ht[0], ht[1]], linewidth=2, color='black', label='2D rootzone')
#ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_top'][start:end, ht[0], ht[1]], linewidth=2, color='green', label='2D topsoil')

#ax0.fill_between(dates_spa[start:end], np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
#                      np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-quantile')

#ax0.scatter(pd.to_datetime(sar_temp['time'][:].data), np.nanmean(sar_temp['theta'], axis=(1,2)), s=50, color='red', label='SAR')
ax0.scatter(pd.to_datetime(sar_temp['time'][:].data), sar_temp['theta'][:,ht[0], ht[1]], s=25, color='red', label='SAR')

ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'mean_obs'], alpha=0.8, label='OBS mean')
ax0.fill_between(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'min'],
                      theta.loc[theta[start_day:end_day].index, 'max'], alpha=0.4, label='OBS range')
ax0.legend(ncol=5)
ax0.set_title('Hilltop timeseries 2019')
ax0.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
ax0.grid()
ax0.xaxis.set_tick_params(rotation=15)
ax0.set_ylim(ylimstemp)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax1);
im1 = res_2d['bucket_moisture_root'][start:end][doi][bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax1, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='2D root')

theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax1, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=45,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax1.legend()
rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = sar_temp['theta'].sel(time=d)[bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=45,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax2.legend()
ax1.axes.get_yaxis().set_ticklabels([])
ax1.set_ylabel('')

'''
rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
im3 = res_2d['bucket_moisture_top'][start:end][doi][bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax3, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='2D')
theta_spat_gpd.loc[theta_spat_gpd.index == d,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax3, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=35,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax3.legend()
ax3.axes.get_yaxis().set_ticklabels([])
ax3.set_ylabel('')
'''

cax = plt.axes([0.1, 0.04, 0.8, 0.03]) # 4-tuple of floats rect = [left, bottom, width, height]. A new axes is added
cbar1 = plt.colorbar(im1, cax=cax, orientation='horizontal')
cbar1.set_label('Volumetric water content [m3/m3]')

plt.subplots_adjust(wspace=0.1, hspace=0.5)


plt.savefig(f'sar_spa_2_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'sar_spa_2_{today}.png', bbox_inches='tight', dpi=300)


#%%

alp=0.6
peat_lims = [0.1,0.9]
min_lims = [0.1,0.9]

fig = plt.figure(figsize=(12,4))
gs = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 2])
ax8 = fig.add_subplot(gs[1, 3])

ax1.hist(d2_sarmask_moist_mineral[june_id_sarmask].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax1.hist(sar_moist_mineral[june_id_sar].flatten(), bins=25, range=min_lims, alpha=alp, label='SAR')
ax1.text(0.60, 2000,  'HM:')
ax1.text(0.7, 2000,  histogram_match(d2_moist_mineral[june_id_sarmask].flatten(),
         sar_moist_mineral[june_id_sar].flatten(), bins=25, lims=min_lims))
ax1.legend()
ax1.set_title('JUNE')
ax1.set_ylabel('occurences')

ax2.hist(d2_sarmask_moist_mineral[july_id_sarmask].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax2.hist(sar_moist_mineral[july_id_sar].flatten(), bins=25, range=min_lims, alpha=alp, label='SAR')
ax2.set_title('JULY')
ax2.text(0.60, 2000,  'HM:')
ax2.text(0.7, 2000,  histogram_match(d2_moist_mineral[july_id_sarmask].flatten(),
         sar_moist_mineral[july_id_sar].flatten(), bins=25, lims=min_lims))

ax3.hist(d2_moist_mineral[august_id_sarmask].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax3.hist(sar_moist_mineral[august_id_sar].flatten(), bins=25, range=min_lims, alpha=alp, label='SAR')
ax3.set_title('AUGUST')
ax3.text(0.60, 2000,  'HM:')
ax3.text(0.7, 2000,  histogram_match(d2_moist_mineral[august_id_sarmask].flatten(),
         sar_moist_mineral[august_id_sar].flatten(), bins=25, lims=min_lims))

ax4.hist(d2_moist_mineral[sept_id_sarmask].flatten(), bins=25, range=min_lims, alpha=alp, label='2D')
ax4.hist(sar_moist_mineral[sept_id_sar].flatten(), bins=25, range=min_lims, alpha=alp, label='SAR')
ax4.set_title('SEPTEMBER')
ax4.text(0.60, 2000,  'HM:')
ax4.text(0.7, 2000,  histogram_match(d2_moist_mineral[sept_id_sarmask].flatten(),
         d1_moist_mineral[sept_id_sar].flatten(), bins=25, lims=min_lims))

###
ax5.hist(d2_moist_peat[june_id_sarmask].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax5.hist(sar_moist_peat[june_id_sar].flatten(), bins=25, range=peat_lims, alpha=alp, label='SAR')
ax5.set_ylabel('occurences')
ax5.text(0.1, 2000,  'HM:')
ax5.text(0.2, 2000,  histogram_match(d2_moist_peat[june_id_sarmask].flatten(),
         sar_moist_peat[june_id_sar].flatten(), bins=25, lims=min_lims))

ax6.hist(d2_moist_peat[july_id_sarmask].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax6.hist(sar_moist_peat[july_id_sar].flatten(), bins=25, range=peat_lims, alpha=alp, label='SAR')
ax6.text(0.1, 2000,  'HM:')
ax6.text(0.2, 2000,  histogram_match(d2_moist_peat[july_id_sarmask].flatten(),
         sar_moist_peat[july_id_sar].flatten(), bins=25, lims=min_lims))

ax7.hist(d2_moist_peat[august_id_sarmask].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax7.hist(sar_moist_peat[august_id_sar].flatten(), bins=25, range=peat_lims, alpha=alp, label='SAR')
ax7.text(0.1, 2000,  'HM:')
ax7.text(0.2, 2000,  histogram_match(d2_moist_peat[august_id_sarmask].flatten(),
         sar_moist_peat[august_id_sar].flatten(), bins=25, lims=min_lims))

ax8.hist(d2_moist_peat[sept_id_sarmask].flatten(), bins=25, range=peat_lims, alpha=alp, label='2D')
ax8.hist(sar_moist_peat[sept_id_sar].flatten(), bins=25, range=peat_lims, alpha=alp, label='SAR')
ax8.text(0.1, 2000,  'HM:')
ax8.text(0.2, 2000,  histogram_match(d2_moist_peat[sept_id_sarmask].flatten(),
         sar_moist_peat[sept_id_sar].flatten(), bins=25, lims=min_lims))

ax1.axes.get_yaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])
ax4.axes.get_yaxis().set_ticklabels([])

ax5.axes.get_yaxis().set_ticklabels([])
ax6.axes.get_yaxis().set_ticklabels([])
ax7.axes.get_yaxis().set_ticklabels([])
ax8.axes.get_yaxis().set_ticklabels([])

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()
ax7.grid()
ax8.grid()

plt.subplots_adjust(wspace=0.05, hspace=0.2)
fig.text(0.5, 0.03, 'Rootzone volumetric water content [m3/m3]', ha='center')

fig.text(0.92, 0.6, 'MINERAL SOIL', rotation='vertical', ha='center')
fig.text(0.92, 0.27, 'PEAT SOIL', rotation='vertical', ha='center')

plt.savefig(f'histogram_sarspa_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'histogram_sarspa_{today}.png', bbox_inches='tight', dpi=300)
