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
from metrics import R2_metrics, MBE_metrics, MAE_metrics, histogram_match

os.chdir(r'C:\SpaFHy_v1_Pallas_2D\figures')
plt.rcParams['axes.grid'] = True

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

'''
# reading SAR files
sar_tempfile = r'D:\SpaFHy_2D_2021\PALLAS_RAW_DATA\SAR_maankosteus\processed\m16_nc_spafhy_pallas\original\SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_all.nc'
sar_temp = xr.open_dataset(sar_tempfile)
sar_temp = sar_temp * cmask

sar_spatfile = r'D:\SpaFHy_2D_2021\PALLAS_RAW_DATA\SAR_maankosteus\processed\m16_nc_spafhy_pallas\interpolated_temporally\SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_TEMPINTERP_all.nc'
sar_spat = xr.open_dataset(sar_spatfile)
sar_spat = sar_spat * cmask
'''

# reading SAR files
sar_mfile = r'D:\SpaFHy_2D_2021\PALLAS_RAW_DATA\SAR_maankosteus\processed\m16_nc_spafhy_pallas\original\SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_morning.nc'
sar_m = xr.open_dataset(sar_mfile) * cmask

sar_efile = r'D:\SpaFHy_2D_2021\PALLAS_RAW_DATA\SAR_maankosteus\processed\m16_nc_spafhy_pallas\original\SAR_SM_PALLAS_2019_16M_TM35_CATCHMENT_evening.nc'
sar_e = xr.open_dataset(sar_efile) * cmask

dates_sar_m = []
dates_sar_e = []

for d in range(len(sar_m['time'])):
    dates_sar_m.append(pd.to_datetime(str(sar_m['time'][d])[36:46]))
    dates_sar_e.append(pd.to_datetime(str(sar_e['time'][d])[36:46]))

dates_sar_m = pd.to_datetime(dates_sar_m).sort_values()
dates_sar_e = pd.to_datetime(dates_sar_e).sort_values()

sar_mask = np.array(sar_m['theta'][0].copy())
sar_mask[sar_mask >= 0] = 1
sar_mask = sar_mask.reshape(sar_m['theta'][0].shape) * cmask

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
    print('Simulation times are a MATCH!')
except ValueError:
    print('Simulation times DO NOT match!')
try:
    res_1d['lat'] in res_2d['lat'] in res_top['lat'] and res_1d['lon'] in res_2d['lon'] in res_top['lon']
    print('Simulation extents are a MATCH!')
except ValueError:
    print('Simulation extents DO NOT match!')

# simple tests for SAR files # !!! if false does not print DO NOT
try:
    sar_m['time'] in sar_e['time']
    print('SAR dates are a MATCH!')
except ValueError:
    print('SAR dates DO NOT match!')
try:
    sar_m['lat'] in sar_e['lat'] and sar_m['lon'] in sar_e['lon']
    print('SAR extents are a MATCH!')
except ValueError:
    print('SAR extents DO NOT match!')

# simple tests for SAR + SpaFHy files
try:
    sar_m['lat'] in sar_e['lat'] in res_1d['lat'] and sar_m['lon'] in sar_e['lon'] in res_1d['lon']
    print('SAR and SpaFHy extents are a MATCH!')
except ValueError:
    print('SAR and SpaFHy extents DO NOT match!')
try:
    sar_m['time'] in sar_e['time'] in res_1d['time']
    print('SAR and SpaFHy times are a MATCH!')
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

##

# spatial manual measurements and closest grid-cell
spat_theta_pair = pd.DataFrame()
obs_list = []
mod_list_1d = []
mod_list_2d = []
mod_list_top = []

for i in theta_spat_gpd.index.unique():
    data = theta_spat_gpd.loc[i]
    #data = data.loc[data['z'] != 0]
    for j in data['ID_i'].unique():
        point = data.loc[data['ID_i'] == j]
        coords = point['geometry'][0]
        time = point.index[0]
        if len(point) > 1:
            theta_obs = float(point['theta'].mean())
        else:
            theta_obs = float(point['theta'])
        obs_list.append(theta_obs)
        mod_list_top.append(float(res_top['bucket_moisture_root'].sel(time=time,
                                                           lat=coords.y, lon=coords.x,
                                                           method='nearest')))
        mod_list_1d.append(float(res_1d['bucket_moisture_root'].sel(time=time,
                                                           lat=coords.y, lon=coords.x,
                                                           method='nearest')))
        mod_list_2d.append(float(res_2d['bucket_moisture_root'].sel(time=time,
                                                           lat=coords.y, lon=coords.x,
                                                           method='nearest')))
spat_theta_pair['obs'] = obs_list
spat_theta_pair['mod_1d'] = mod_list_1d
spat_theta_pair['mod_2d'] = mod_list_2d
spat_theta_pair['mod_top'] = mod_list_top

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


ax3 = fig.add_subplot(gs[1, :3])
ax3.set_title('Open mire')
im3 = temporaldf['2d_bucket_moisture_root_om'].plot(ax=ax3, color='g', alpha=0.7, label='2D')
temporaldf['top_bucket_moisture_root_om'].plot(ax=ax3, color='red', alpha=0.6, label='TOP')
ax1.legend(ncol=5)
ax3.set_ylim(0.35,0.9)
ax3.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#f3_ax3.legend(ncol=5)


if saveplots == True:
        plt.savefig(f'theta_model_ts_{today}.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(f'theta_model_ts_{today}.png', bbox_inches='tight', dpi=300)


#%%
'''
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
'''


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
'''
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


#plt.savefig(f'spatial_mode_{d}_2021_{today}.pdf', bbox_inches='tight', dpi=300)
#plt.savefig(f'spatial_mode_{d}_2021_{today}.png', bbox_inches='tight', dpi=300)
'''
#%%

# SPAFHY SPATIAL MODEL VERSION COMPARISON WET AND DRY DAY

start = np.where(pd.to_datetime(dates_spa) == '2021-05-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2021-09-02')[0][0]

d1 = '2021-06-17'
doi1 = np.where(pd.to_datetime(dates_spa[start:end]) == d1)[0][0]
doi_m1 = np.where(pd.to_datetime(theta_spat['time'][:].data) == d1)[0][0]
d2 = '2021-09-01'
doi2 = np.where(pd.to_datetime(dates_spa[start:end]) == d2)[0][0]
doi_m2 = np.where(pd.to_datetime(theta_spat['time'][:].data) == d2)[0][0]

alp=0.7
ylims = [0.15,0.88]
fig = plt.figure(figsize=(12,10))
gs = fig.add_gridspec(2, 3)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])

rasterio.plot.show(pk, transform=meta['transform'], ax=ax1);
im1 = res_1d['bucket_moisture_root'][start:end][doi1].plot(ax=ax1, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='MOD')
theta_spat_gpd.loc[theta_spat_gpd.index == d1,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax1, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax1.legend()
rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = res_top['bucket_moisture_root'][start:end][doi1].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d1,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax2.legend()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3)

res_2d['bucket_moisture_root'][start:end][doi1].plot(ax=ax3, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                    add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d1,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax3, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax3.legend()
                                              #theta_spat_gpd.loc[theta_spat_gpd.index == d][['SM_mean', 'geometry']].plot()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax4);
im4 = res_1d['bucket_moisture_root'][start:end][doi2].plot(ax=ax4, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='MOD')
theta_spat_gpd.loc[theta_spat_gpd.index == d2,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax4, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax4.legend()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax5);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im5 = res_top['bucket_moisture_root'][start:end][doi2].plot(ax=ax5, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d2,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax5, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax5.legend()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax6)

res_2d['bucket_moisture_root'][start:end][doi2].plot(ax=ax6, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                    add_colorbar=False)
theta_spat_gpd.loc[theta_spat_gpd.index == d2,
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax6, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
                                               markersize=25,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
# TEST to see all locations
#theta_spat_gpd.loc[(theta_spat_gpd.index != d1) & (theta_spat_gpd.index != d2)][['theta', 'geometry']].plot(column='theta',
#                                               ax=ax6, cmap='coolwarm_r', edgecolor='black', linewidth=0.5,
#                                              markersize=25,
#                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax6.legend()
                                              #theta_spat_gpd.loc[theta_spat_gpd.index == d][['SM_mean', 'geometry']].plot()

props = dict(facecolor='white', alpha=1)

ax1.text(0.07, 0.96,'1D', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes,
     bbox=props)
ax2.text(0.07, 0.96,'TOP', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes,
     bbox=props)
ax3.text(0.07, 0.96,'2D', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax3.transAxes,
     bbox=props)

ax4.text(0.07, 0.96,'1D', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax4.transAxes,
     bbox=props)
ax5.text(0.07, 0.96,'TOP', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax5.transAxes,
     bbox=props)
ax6.text(0.07, 0.96,'2D', size=10,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax6.transAxes,
     bbox=props)

ax2.axes.get_yaxis().set_ticklabels([])
ax2.set_ylabel('')
ax3.axes.get_yaxis().set_ticklabels([])
ax3.set_ylabel('')

ax1.axes.get_xaxis().set_ticklabels([])
ax1.set_xlabel('')
ax2.axes.get_xaxis().set_ticklabels([])
ax2.set_xlabel('')
ax3.axes.get_xaxis().set_ticklabels([])
ax3.set_xlabel('')

ax5.axes.get_yaxis().set_ticklabels([])
ax5.set_ylabel('')
ax6.axes.get_yaxis().set_ticklabels([])
ax6.set_ylabel('')

cax = plt.axes([0.12, 0.06, 0.78, 0.02])
cbar1 = plt.colorbar(im1, cax=cax, orientation='horizontal')
cbar1.set_label(r'Volumetric water content [m$^3$/m$^3$]')

plt.subplots_adjust(wspace=0.1, hspace=0.05)


plt.savefig(f'spatial_models_2021_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'spatial_models__2021_{today}.png', bbox_inches='tight', dpi=300)


#%%

r2_spa_theta_1d = R2_metrics(spat_theta_pair, columns=['obs', 'mod_1d'])
r2_spa_theta_top =  R2_metrics(spat_theta_pair, columns=['obs', 'mod_top'])
r2_spa_theta_2d =  R2_metrics(spat_theta_pair, columns=['obs', 'mod_2d'])

mbe_spa_theta_1d = MBE_metrics(spat_theta_pair, columns=['obs', 'mod_1d'])
mbe_spa_theta_top =  MBE_metrics(spat_theta_pair, columns=['obs', 'mod_top'])
mbe_spa_theta_2d =  MBE_metrics(spat_theta_pair, columns=['obs', 'mod_2d'])

mae_spa_theta_1d = MAE_metrics(spat_theta_pair, columns=['obs', 'mod_1d'])
mae_spa_theta_top =  MAE_metrics(spat_theta_pair, columns=['obs', 'mod_top'])
mae_spa_theta_2d =  MAE_metrics(spat_theta_pair, columns=['obs', 'mod_2d'])

# SCATTERPLOTS OF SPATIAL SOIL MOISTURE MEASUREMENTS AND SIMULATIONS
alp=0.5
ylims = [0.1,1.0]
text1 = [0.01, 0.8]

fig = plt.figure(figsize=(12,4))
gs = fig.add_gridspec(1, 3)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

ax1.scatter(spat_theta_pair['obs'], spat_theta_pair['mod_1d'], alpha=alp)
ax1.plot(ylims, ylims, 'r--', alpha=0.6)
ax2.scatter(spat_theta_pair['obs'], spat_theta_pair['mod_top'], alpha=alp)
ax2.plot(ylims, ylims, 'r--', alpha=0.6)
ax3.scatter(spat_theta_pair['obs'], spat_theta_pair['mod_2d'], alpha=alp)
ax3.plot(ylims, ylims, 'r--', alpha=0.6)

ax1.set_ylim(ylims)
ax1.set_xlim(ylims)
ax2.set_ylim(ylims)
ax2.set_xlim(ylims)
ax3.set_ylim(ylims)
ax3.set_xlim(ylims)

ax1.set_xlabel('OBS')
ax1.set_ylabel('MOD')
ax2.set_xlabel('OBS')
ax3.set_xlabel('OBS')

ax1.set_title('1D vs. OBS')
ax2.set_title('TOP vs. OBS')
ax3.set_title('2D vs. OBS')

ax1.text(text1[0], text1[1], f'R2:{r2_spa_theta_1d}\nMAE:{mae_spa_theta_1d}\nMBE:{mbe_spa_theta_1d}',
         fontsize=14, transform=ax1.transAxes)
ax2.text(text1[0], text1[1], f'R2:{r2_spa_theta_top}\nMAE:{mae_spa_theta_top}\nMBE:{mbe_spa_theta_top}',
         fontsize=14, transform=ax2.transAxes)
ax3.text(text1[0], text1[1], f'R2:{r2_spa_theta_2d}\nMAE:{mae_spa_theta_2d}\nMBE:{mbe_spa_theta_2d}',
         fontsize=14, transform=ax3.transAxes)

ax1.set_xticks(np.arange(ylims[0], ylims[1], 0.2))
ax1.set_yticks(np.arange(ylims[0], ylims[1], 0.2))
ax2.set_xticks(np.arange(ylims[0], ylims[1], 0.2))
ax2.set_yticks(np.arange(ylims[0], ylims[1], 0.2))
ax3.set_xticks(np.arange(ylims[0], ylims[1], 0.2))
ax3.set_yticks(np.arange(ylims[0], ylims[1], 0.2))

plt.subplots_adjust(wspace=0.12, hspace=0)

plt.savefig(f'spatial_theta_scat_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'spatial_theta_scat_{today}.png', bbox_inches='tight', dpi=300)

#%%

# HISTOGRAMS

june_id = np.where(pd.to_datetime(dates_spa).month == 6)[0]
june_id_sar = np.where(pd.to_datetime(dates_sar_m).month == 6)[0]
june_id_sarmask = np.where((dates_spa.isin(dates_sar_m)) & (dates_spa.month == 6))[0]

july_id = np.where(pd.to_datetime(dates_spa).month == 7)[0]
july_id_sar = np.where(pd.to_datetime(dates_sar_m).month == 7)[0]
july_id_sarmask = np.where((dates_spa.isin(dates_sar_m)) & (dates_spa.month == 7))[0]

august_id = np.where(pd.to_datetime(dates_spa).month == 8)[0]
august_id_sar = np.where(pd.to_datetime(dates_sar_m).month == 8)[0]
august_id_sarmask = np.where((dates_spa.isin(dates_sar_m)) & (dates_spa.month == 8))[0]

sept_id = np.where(pd.to_datetime(dates_spa).month == 9)[0]
sept_id_sar = np.where(pd.to_datetime(dates_sar_m).month == 9)[0]
sept_id_sarmask = np.where((dates_spa.isin(dates_sar_m)) & (dates_spa.month == 9))[0]


peat_id = np.where(np.ravel(soilclass == 4))[0]
mineral_id = np.where(np.ravel(soilclass == 2))[0]

top_moist_peat = np.array(res_top['bucket_moisture_root']).reshape(res_top['bucket_moisture_root'].shape[0], -1)[:,peat_id]
d2_moist_peat = np.array(res_2d['bucket_moisture_root']).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,peat_id]
d1_moist_peat = np.array(res_1d['bucket_moisture_root']).reshape(res_1d['bucket_moisture_root'].shape[0], -1)[:,peat_id]
sar_moist_peat = np.array(sar_m['theta']).reshape(sar_m['theta'].shape[0], -1)[:,peat_id]
d2_sarmask_moist_peat = np.array(res_2d['bucket_moisture_root']*sar_mask).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,peat_id]


top_moist_mineral = np.array(res_top['bucket_moisture_root']).reshape(res_top['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
d2_moist_mineral = np.array(res_2d['bucket_moisture_root']).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
d1_moist_mineral = np.array(res_1d['bucket_moisture_root']).reshape(res_1d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]
sar_moist_mineral = np.array(sar_m['theta']).reshape(sar_m['theta'].shape[0], -1)[:,mineral_id]
d2_sarmask_moist_mineral = np.array(res_2d['bucket_moisture_root']*sar_mask).reshape(res_2d['bucket_moisture_root'].shape[0], -1)[:,mineral_id]

alp=0.6
peat_lims = [0.1,0.9]
min_lims = [0.1,0.9]



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

#ax1.grid()
#ax2.grid()
#ax3.grid()
#ax4.grid()
#ax5.grid()
#ax6.grid()
#ax7.grid()
#ax8.grid()

plt.subplots_adjust(wspace=0.05, hspace=0.2)
fig.text(0.5, 0.03, 'Rootzone volumetric water content [m3/m3]', ha='center')

fig.text(0.92, 0.6, 'MINERAL SOIL', rotation='vertical', ha='center')
fig.text(0.92, 0.27, 'PEAT SOIL', rotation='vertical', ha='center')



plt.savefig(f'histogram_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'histogram_{today}.png', bbox_inches='tight', dpi=300)

#%%

# temporal sar model obs

# SAR SPATIAL PLOTS
start_day = pd.to_datetime('2019-05-01')
end_day = pd.to_datetime('2019-09-15')
window = pd.date_range(start_day, end_day, freq='D')

window_sar = window[window.isin(dates_sar_m)]


alp=0.5
meas_depth=10
ylims = [0.2,0.88]
ylimstemp = [0.08,0.5]
fig = plt.figure(figsize=(10,12))
gs = fig.add_gridspec(18, 2)

ax0 = fig.add_subplot(gs[0:3, :])
ax1 = fig.add_subplot(gs[4:7, 0])
ax2 = fig.add_subplot(gs[4:7, 1])
ax3 = fig.add_subplot(gs[7:10, 0])
ax4 = fig.add_subplot(gs[7:10, 1])
ax5 = fig.add_subplot(gs[10:13, 0])
ax6 = fig.add_subplot(gs[10:13, 1])
ax7 = fig.add_subplot(gs[13:16, 0])
ax8 = fig.add_subplot(gs[13:16, 1])

#ax0.plot(dates_spa[start:end], np.nanmean(res_2d['bucket_moisture_root'][start:end], axis=(1,2)), linewidth=2, color='black', label='mean')
ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_root'][start:end, ht[0], ht[1]], linewidth=1, color='black', label='2D rootzone')
#ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_top'][start:end, ht[0], ht[1]], linewidth=2, color='black', label='2D rootzone')

#ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_top'][start:end, ht[0], ht[1]], linewidth=2, color='green', label='2D topsoil')

#ax0.fill_between(dates_spa[start:end], np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
#                      np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-quantile')

#ax0.scatter(pd.to_datetime(sar_temp['time'][:].data), np.nanmean(sar_temp['theta'], axis=(1,2)), s=50, color='red', label='SAR')
ax0.scatter(pd.to_datetime(sar_m['time'][:].data), sar_m['theta'][:,ht[0], ht[1]], marker='o', s=15, alpha=0.5, color='tab:red', label='SAR (m)')
ax0.plot(pd.to_datetime(sar_m['time'][:].data[4:]),
         sar_m['theta'][:,ht[0], ht[1]].rolling(time=5, center=True).mean().dropna("time"), color='tab:red', alpha=0.5, label='SAR rolling mean')

ax0.scatter(pd.to_datetime(sar_e['time'][:].data), sar_e['theta'][:,ht[0], ht[1]], s=15, color='tab:red', marker='d', alpha=0.5, label='SAR (e)')
ax0.plot(pd.to_datetime(sar_e['time'][:].data[4:]),
         sar_e['theta'][:,ht[0], ht[1]].rolling(time=5, center=True).mean().dropna("time"), color='tab:red', alpha=0.5)
#ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'mean_obs'], alpha=0.8, label='OBS mean')
#ax0.fill_between(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'min'],
#                      theta.loc[theta[start_day:end_day].index, 'max'], alpha=0.4, label='OBS range')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 's3'], alpha=0.7, color='tab:blue', label='OBS')
#ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 's5'], alpha=0.8, color='tab:blue')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 's18'], alpha=0.7, color='tab:blue')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'SH-5A'], alpha=0.7, color='tab:blue')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'SH-5B'], alpha=0.7, color='tab:blue')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'SH-20A'], alpha=0.7, color='tab:blue')
ax0.plot(theta[start_day:end_day].index, theta.loc[theta[start_day:end_day].index, 'SH-20B'], alpha=0.7, color='tab:blue')

ax0.legend(ncol=5,bbox_to_anchor=(0.9,1.3))
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
#ax0.xaxis.set_tick_params(rotation=15)
ax0.set_ylim(ylimstemp)

# ax1
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].x[0], time=window, method='nearest').plot(
        ax=ax1, color='black')

ax1.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i13', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax1, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i13') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax1, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i13') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax1, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i13') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax1, alpha=0.7, color='tab:blue')
# ax2
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax2, color='black')

ax2.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i15', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax2, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i15') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax2, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i15') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax2, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i15') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax2, alpha=0.7, color='tab:blue')
# ax3
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax3, color='black')

ax3.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i20', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax3, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i20') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax3, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i20') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax3, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i20') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax3, alpha=0.7, color='tab:blue')

# ax4
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax4, color='black')

ax4.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i21', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax4, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i21') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax4, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i21') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax4, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i21') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax4, alpha=0.7, color='tab:blue')

# ax5
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax5, color='black')

ax5.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i22', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax5, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i22') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax5, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i22') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax5, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i22') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax5, alpha=0.7, color='tab:blue')

# ax6
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax6, color='black')

ax6.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i18', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax6, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i18') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax6, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i18') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax6,  alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i18') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax6,  alpha=0.7, color='tab:blue')

# ax7
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax7, color='black')

ax7.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i14', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax7, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i14') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax7, alpha=0.7,  color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i14') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax7, alpha=0.7,  color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i14') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax7, alpha=0.7,  color='tab:blue')


# ax8
res_2d['bucket_moisture_root'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].x[0],
    time=window, method='nearest').plot(ax=ax8, color='black')

ax8.scatter(window_sar, sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].x[0], time=window_sar, method='nearest'),
    color='tab:red', alpha=0.7, s=15)

sar_m['theta'].sel(
    lat=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].y[0],
    lon=theta_spat_gpd.loc[theta_spat_gpd['ID_i'] == 'i16', 'geometry'].x[0],
    time=window, method='nearest').rolling(time=5, center=True).mean().plot(ax=ax8, color='tab:red', alpha=0.7)

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i16') &
                   (theta_spat_gpd['z'] == meas_depth) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax8, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i16') &
                   (theta_spat_gpd['z'] == 20) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax8, alpha=0.7, color='tab:blue')

theta_spat_gpd.loc[(theta_spat_gpd['ID_i'] == 'i16') &
                   (theta_spat_gpd['z'] == 30) & (theta_spat_gpd.index >= start_day) &
                   (theta_spat_gpd.index <= end_day), 'theta'].plot(ax=ax8, alpha=0.7, color='tab:blue')

ax1.axes.xaxis.set_ticklabels([])
ax2.axes.xaxis.set_ticklabels([])
ax3.axes.xaxis.set_ticklabels([])
ax4.axes.xaxis.set_ticklabels([])
ax5.axes.xaxis.set_ticklabels([])
ax6.axes.xaxis.set_ticklabels([])

ax1.set_title('')
ax2.set_title('')
ax3.set_title('')
ax4.set_title('')
ax5.set_title('')
ax6.set_title('')
ax7.set_title('')
ax8.set_title('')

ax1.set_xlabel('')
ax2.set_xlabel('')
ax3.set_xlabel('')
ax4.set_xlabel('')
ax5.set_xlabel('')
ax6.set_xlabel('')
ax7.set_xlabel('')
ax8.set_xlabel('')

ax0.set_xticks(pd.date_range(start_day, end_day, freq='10D'))
ax1.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax2.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax3.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax4.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax5.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax6.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax7.set_xticks(pd.date_range(start_day, end_day, freq='MS'))
ax8.set_xticks(pd.date_range(start_day, end_day, freq='MS'))

ax0.xaxis.set_tick_params(rotation=15)
ax7.xaxis.set_tick_params(rotation=15)
ax8.xaxis.set_tick_params(rotation=15)

ax0.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
ax1.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
ax3.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
ax5.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
ax7.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')



plt.subplots_adjust(wspace=0.1, hspace=0.2)


plt.savefig(f'sar_spa_temp_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'sar_spa_temp_{today}.png', bbox_inches='tight', dpi=300)

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

sar_wet_day = np.where(np.nansum(sar_m['theta'], axis=(1,2)) == np.nansum(sar_m['theta'], axis=(1,2)).max())
sar_dry_day = np.where(np.nansum(sar_m['theta'], axis=(1,2)) == np.nansum(sar_m['theta'], axis=(1,2)).min())

wet_day = pd.to_datetime(sar_m['time'][sar_wet_day].data[0])
dry_day = pd.to_datetime(sar_m['time'][sar_dry_day].data[0])

doi = np.where(pd.to_datetime(dates_spa[start:end]) == wet_day)[0][0]
d = pd.to_datetime(wet_day)

d = pd.to_datetime(pd.to_datetime('2019-06-26'))
doi = np.where(pd.to_datetime(dates_spa[start:end]) == '2019-06-26')[0][0]

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
ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_root'][start:end, ht[0], ht[1]], linewidth=1, color='black', label='2D rootzone')
#ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_top'][start:end, ht[0], ht[1]], linewidth=2, color='black', label='2D rootzone')

#ax0.plot(dates_spa[start:end], res_2d['bucket_moisture_top'][start:end, ht[0], ht[1]], linewidth=2, color='green', label='2D topsoil')

#ax0.fill_between(dates_spa[start:end], np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
#                      np.nanquantile(res_2d['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-quantile')

#ax0.scatter(pd.to_datetime(sar_temp['time'][:].data), np.nanmean(sar_temp['theta'], axis=(1,2)), s=50, color='red', label='SAR')
ax0.scatter(pd.to_datetime(sar_m['time'][:].data), sar_m['theta'][:,ht[0], ht[1]], marker='d', s=15, alpha=0.7, color='red', label='SAR (m)')
ax0.plot(pd.to_datetime(sar_m['time'][:].data[4:]),
         sar_m['theta'][:,ht[0], ht[1]].rolling(time=5, center=True).mean().dropna("time"), color='red', alpha=0.7)
#ax0.plot(dates_sar_m, sar_m['theta'][:,ht[0], ht[1]], color='red', alpha=0.7, label='SAR')
ax0.scatter(pd.to_datetime(sar_e['time'][:].data), sar_e['theta'][:,ht[0], ht[1]], s=15, color='red', alpha=0.7, label='SAR (e)')
ax0.plot(pd.to_datetime(sar_e['time'][:].data[4:]),
         sar_e['theta'][:,ht[0], ht[1]].rolling(time=5, center=True).mean().dropna("time"), color='red', alpha=0.7)
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
im2 = sar_m['theta'].sel(time=d)[bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
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

import matplotlib.patches as patches


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# SAR SPATIAL PLOTS
start_day = pd.to_datetime('2019-05-01')
end_day = pd.to_datetime('2019-09-15')

start = np.where(pd.to_datetime(dates_spa) == start_day)[0][0]
end = np.where(pd.to_datetime(dates_spa) == end_day)[0][0]

sar_wet_day = np.where(np.nansum(sar_m['theta'], axis=(1,2)) == np.nansum(sar_m['theta'], axis=(1,2)).max())
sar_dry_day = np.where(np.nansum(sar_m['theta'], axis=(1,2)) == np.nansum(sar_m['theta'], axis=(1,2)).min())

wet_day = pd.to_datetime(sar_m['time'][sar_wet_day].data[0])
dry_day = pd.to_datetime(sar_m['time'][sar_dry_day].data[0])

doi = np.where(pd.to_datetime(dates_spa[start:end]) == wet_day)[0][0]
d = pd.to_datetime(wet_day)

date = '2019-07-08'
d = pd.to_datetime(pd.to_datetime(date))
doi = np.where(pd.to_datetime(dates_spa[start:end]) == date)[0][0]

bbox = [7.5452e6, 7.5443e6, 383300, 384800]
bbox_id = [np.where(res_1d['lat'] == find_nearest(res_1d['lat'],bbox[0]))[0][0],
           np.where(res_1d['lat'] == find_nearest(res_1d['lat'],bbox[1]))[0][0],
           np.where(res_1d['lon'] == find_nearest(res_1d['lon'],bbox[2]))[0][0],
           np.where(res_1d['lon'] == find_nearest(res_1d['lon'],bbox[3]))[0][0]]
xs = [bbox[2], bbox[3], bbox[3], bbox[2], bbox[2]]
ys = [bbox[0], bbox[0], bbox[1], bbox[1], bbox[0]]

props = dict(facecolor='white', alpha=1)


alp=0.7
ylims = [0.2,0.88]
ylimstemp = [0.08,0.5]
fig = plt.figure(figsize=(12,12))
gs = fig.add_gridspec(10, 10)

ax0 = fig.add_subplot(gs[0:7, 0:5])
ax1 = fig.add_subplot(gs[0:7, 5:10])
ax2 = fig.add_subplot(gs[7:10, 0:5])
ax3 = fig.add_subplot(gs[7:10, 5:10])

rasterio.plot.show(pk, transform=meta['transform'], ax=ax0);
im0 = res_2d['bucket_moisture_root'][start:end][doi].plot(ax=ax0, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='2D root')

theta_spat_gpd.loc[(theta_spat_gpd.index == d) & (theta_spat_gpd['z'] == 20),
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax0, cmap='coolwarm_r', edgecolor='black', linewidth=1,
                                               markersize=45,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')

ax0.plot(xs, ys, color="black")
ax0.text(0.06, 0.97,'2D', size=15,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax0.transAxes,
     bbox=props)

ax0.legend()


rasterio.plot.show(pk, transform=meta['transform'], ax=ax0)
rasterio.plot.show(pk, transform=meta['transform'], ax=ax1);
im1 = sar_m['theta'].sel(time=d).plot(ax=ax1, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)

theta_spat_gpd.loc[(theta_spat_gpd.index == d) & (theta_spat_gpd['z'] == 20),
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax1, cmap='coolwarm_r', edgecolor='black', linewidth=1,
                                               markersize=45,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax1.plot(xs, ys, color="black")
ax1.text(0.06, 0.97,'SAR', size=15,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes,
     bbox=props)
ax1.legend()

                                      # Create a Rectangle patch
#rect = patches.Rectangle((bbox[1], bbox[2]), bbox[0], bbox[3], linewidth=1, edgecolor='r', facecolor='none')
# Add the patch to the Axes
#ax0.add_patch(rect)
# Coordinates of rectangle vertices
# in clockwise order

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2)
im2 = res_2d['bucket_moisture_root'][start:end][doi][bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax2, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False, label='2D root')

theta_spat_gpd.loc[(theta_spat_gpd.index == d) & (theta_spat_gpd['z'] == 20),
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax2, cmap='coolwarm_r', edgecolor='black', linewidth=1,
                                               markersize=65,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax2.legend()
ax2.text(0.06, 0.93,'2D', size=15,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax2.transAxes,
     bbox=props)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im3 = sar_m['theta'].sel(time=d)[bbox_id[0]:bbox_id[1],bbox_id[2]:bbox_id[3]].plot(ax=ax3, cmap='coolwarm_r', vmin=ylims[0], vmax=ylims[1], alpha=alp,
                                                          add_colorbar=False)
theta_spat_gpd.loc[(theta_spat_gpd.index == d) & (theta_spat_gpd['z'] == 20),
                   ['theta', 'geometry']].plot(column='theta',
                                               ax=ax3, cmap='coolwarm_r', edgecolor='black', linewidth=1,
                                               markersize=65,
                                               alpha=alp, vmin=ylims[0], vmax=ylims[1], label='OBS')
ax3.legend()
ax3.text(0.06, 0.93,'SAR', size=15,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax3.transAxes,
     bbox=props)
ax1.axes.get_yaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])

ax0.set_xlabel('')
ax1.set_xlabel('')

ax1.set_ylabel('')
ax3.set_ylabel('')

cax = plt.axes([0.12, 0.06, 0.77, 0.015]) # 4-tuple of floats rect = [left, bottom, width, height]. A new axes is added
cbar1 = plt.colorbar(im1, cax=cax, orientation='horizontal')
cbar1.set_label('Volumetric water content [m3/m3]')

ax2.set_title('')
ax3.set_title('')

plt.subplots_adjust(wspace=0.05, hspace=0.05)


plt.savefig(f'sar_spa_spat_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'sar_spa_spat_{today}.png', bbox_inches='tight', dpi=300)


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

#ax1.grid()
#ax2.grid()
#ax3.grid()
#ax4.grid()
#ax5.grid()
#ax6.grid()
#ax7.grid()
#ax8.grid()

plt.subplots_adjust(wspace=0.05, hspace=0.2)
fig.text(0.5, 0.03, 'Rootzone volumetric water content [m3/m3]', ha='center')

fig.text(0.92, 0.6, 'MINERAL SOIL', rotation='vertical', ha='center')
fig.text(0.92, 0.27, 'PEAT SOIL', rotation='vertical', ha='center')

plt.savefig(f'histogram_sarspa_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'histogram_sarspa_{today}.png', bbox_inches='tight', dpi=300)
