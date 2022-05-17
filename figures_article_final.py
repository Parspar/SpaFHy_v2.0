# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 12:26:01 2021

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
from iotools import read_AsciiGrid, write_AsciiGrid
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from sklearn.metrics import mean_absolute_error as mae
from rasterio.plot import show
import rasterio
from rasterio.transform import from_origin
from raster_utils import read_pkrasteri_for_extent

# reading the stand results
outputfile_stand = r'D:\SpaFHy_2D_2021\testcase_input_1d_new.nc'
results_stand = read_results(outputfile_stand)

# reading the stand results
outputfile_2d = r'D:\SpaFHy_2D_2021\testcase_input_2d_new.nc'
results_2d = read_results(outputfile_2d)

# reading the catch results
outputfile_catch = r'D:\SpaFHy_2D_2021\testcase_input_top_new_fixed.nc'
results_catch = read_results(outputfile_catch)

sar_file = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment_ma3.nc'
sar = Dataset(sar_file, 'r')

sar_file2 = 'C:\SpaFHy_v1_Pallas_2D/obs/SAR_PALLAS_2019_mask2_16m_direct_catchment.nc'
sar2 = Dataset(sar_file2, 'r')

# water table at lompolonjänkä

# cell locations of kenttarova
kenttarova, _, _, _, _ = read_AsciiGrid(r'C:\PALLAS_RAW_DATA\Lompolonjanka\16b\sve_kenttarova_soilmoist.asc')
kenttarova_loc = np.where(kenttarova == 0)
k_loc = list([int(kenttarova_loc[1]), int(kenttarova_loc[0])])
k_loc = [118, 136]
l_loc = [46, 54]


# READING CMASK FROM INPUT FILES
cmaskfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\cmask.dat'
cmaskascii = read_AsciiGrid(cmaskfp)

cmask = rasterio.open(cmaskfp, 'r')
bbox = cmask.bounds
transform = from_origin(bbox[0], bbox[1], cmask.transform[0],  cmask.transform[0])
del cmask, cmaskfp

# reading basic map
pkfp = 'C:\SpaFHy_v1_Pallas_2D/testcase_input/parameters/pkmosaic_clipped.tif'
pk, meta = read_pkrasteri_for_extent(pkfp, bbox=bbox,showfig=True)


dates_spa = []
for d in range(len(results_stand['date'])):
    dates_spa.append(pd.to_datetime(str(results_stand['date'][d])[36:46]))

dates_sar = pd.to_datetime(sar['time'][:], format='%Y%m%d')

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
ix_pd = []
for p in range(len(ix_p)):
    ix_pd.append(dates_spa[p])
    ix_pp = ix_p - 1
ix_pp = ix_pp[ix_pp >= 0]
ix_ppd = []
for p in range(len(ix_pp)):
    ix_ppd.append(dates_spa[p])
#ix_no_p_d = forc.index[np.where(forc['rainfall'] == 0)[0]]
rainday = forc.index[forc['rainfall'] > 0]
rainday_1 = rainday.shift(1, freq='D')
wetdays = rainday.append(rainday_1).sort_values()

# specific discharge
folder = r'C:\SpaFHy_v1_Pallas\data\obs'
ffile = 'Runoffs1d_SVEcatchments_mmd_new.csv'
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

#%%

# GIS RASTERS INTO COORDINATES

# READING CMASK FROM INPUT FILES
laiconiffp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\LAI_conif.dat'
LAI_conif = rasterio.open(laiconiffp, 'r')

laidecidfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\LAI_decid.dat'
LAI_decid = rasterio.open(laidecidfp, 'r')

laishrubfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\LAI_shrub.dat'
LAI_shrub = rasterio.open(laishrubfp, 'r')

laigrassfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\LAI_grass.dat'
LAI_grass = rasterio.open(laigrassfp, 'r')

hcfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\hc.dat'
hc = rasterio.open(hcfp, 'r')

cffp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\cf.dat'
cf = rasterio.open(cffp, 'r')

soilclassfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\soil_id_peatsoils.dat'
soilclass = rasterio.open(soilclassfp, 'r')

ditchesfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\ditches.dat'
ditches = rasterio.open(ditchesfp, 'r')

demfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\dem_raw.dat'
dem = rasterio.open(demfp, 'r')

twifp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\twi.dat'
twi = rasterio.open(twifp, 'r')

#%%

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

laiplot = rasterio.plot.show(LAI_conif, transform=LAI_conif.transform, ax=ax1, alpha=alp)
# this creates colorbar
im1 = laiplot.get_images()[1]
#ax1[0].set_title('Site 1: Sat deficit doy 123-126')
fig.colorbar(im1, ax=ax1, shrink=shr)
ax1.axes.get_xaxis().set_ticks([])
ax1.axes.get_yaxis().set_ticks([])

ax1.set_title(r'Lehtialaindeksi [m$^2$/m$^2$]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
hcplot = rasterio.plot.show(hc, transform=hc.transform, ax=ax2, alpha=alp)
im2 = hcplot.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
fig.colorbar(im2, ax=ax2, shrink=shr)
ax2.axes.get_xaxis().set_ticks([])
ax2.axes.get_yaxis().set_ticks([])
ax2.set_title('Kasvuston korkeus [m]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
cfplot = rasterio.plot.show(cf, transform=cf.transform, ax=ax3, alpha=alp)
im3 = cfplot.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
fig.colorbar(im3, ax=ax3, shrink=shr)
ax3.axes.get_xaxis().set_ticks([])
ax3.axes.get_yaxis().set_ticks([])
ax3.set_title('Kasvuston osuus [m$^2$/m$^2$]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax4);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
scplot = rasterio.plot.show(soilclass, transform=soilclass.transform, ax=ax4, alpha=alp)
im4 = scplot.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
cbs = fig.colorbar(im4, ax=ax4, shrink=shr, ticks=[1,2,3,4])
ax4.set_title('Pintamaaluokka [-]')
ax4.axes.get_xaxis().set_ticks([])
ax4.axes.get_yaxis().set_ticks([])

rasterio.plot.show(pk, transform=meta['transform'], ax=ax5);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
demplot = rasterio.plot.show(dem, transform=dem.transform, ax=ax5, alpha=alp)
im5 = demplot.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
fig.colorbar(im5, ax=ax5, shrink=shr)
ax5.axes.get_yaxis().set_ticks([])
ax5.axes.get_xaxis().set_ticks([])
ax5.set_title('Korkeus [m]')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax6);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
twiplot = rasterio.plot.show(twi, transform=twi.transform, vmin=4, vmax=14, ax=ax6, alpha=alp)
im6 = twiplot.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
fig.colorbar(im6, ax=ax6, shrink=shr)
ax6.axes.get_yaxis().set_ticks([])
ax6.axes.get_xaxis().set_ticks([])
ax6.set_title('Kosteusindeksi [-]')

plt.subplots_adjust(wspace=0.05, hspace=0.05)

plt.savefig(f'vesitalous_gis_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'vesitalous_gis_{today}.png', bbox_inches='tight', dpi=300)

#%%

# vesitalous artikkeli tuloskuvat

results_catch['canopy_total_et'] = (results_catch['canopy_evaporation']
                                    + results_catch['canopy_transpiration']
                                    + results_catch['bucket_evaporation'])

start = np.where(pd.to_datetime(dates_spa) == '2020-05-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2020-10-01')[0][0]
maxd = np.where(np.nansum(results_catch['bucket_moisture_root'][start:end], axis=(1,2))
                == np.nansum(results_catch['bucket_moisture_root'][start:end], axis=(1,2)).max())[0][0]
mind = np.where(np.nansum(results_catch['bucket_moisture_root'][start:end], axis=(1,2))
                == np.nansum(results_catch['bucket_moisture_root'][start:end], axis=(1,2)).min())[0][0]
mind = np.where(pd.to_datetime(dates_spa[start:end]) == '2020-09-02')[0][0]

# writing into ascii and reading in with rasterio
snap1fp = r'C:\SpaFHy_v1_Pallas_2D\results\temp\snap1.dat'
snap2fp = r'C:\SpaFHy_v1_Pallas_2D\results\temp\snap2.dat'

snap1 = results_catch['bucket_moisture_root'][start:end][maxd]
snap2 = results_catch['bucket_moisture_root'][start:end][mind]

write_AsciiGrid(snap1fp, np.array(snap1), info=cmaskascii[1])
write_AsciiGrid(snap2fp, np.array(snap2), info=cmaskascii[1])

snap1 = rasterio.open(snap1fp, 'r')
snap2 = rasterio.open(snap2fp, 'r')

#%%

# Plotting

alp=0.8

fig = plt.figure(figsize=(8,8))
gs = fig.add_gridspec(5, 2)

ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[1:4, 0])
ax3 = fig.add_subplot(gs[1:4, 1])

ax1.plot(dates_spa[start:end], np.nanmean(results_catch['bucket_moisture_root'][start:end], axis=(1,2)), linewidth=2, color='black', label='keskiarvo')
ax1.fill_between(dates_spa[start:end], np.nanquantile(results_catch['bucket_moisture_root'][start:end], 0.8, axis=(1,2)),
                      np.nanquantile(results_catch['bucket_moisture_root'][start:end], 0.2, axis=(1,2)), alpha=0.3, label='20-kvantiili')
ax1.scatter(dates_spa[start:end][maxd], np.nanmean(results_catch['bucket_moisture_root'][start:end][maxd]), s=70, color='blue')
ax1.scatter(dates_spa[start:end][mind], np.nanmean(results_catch['bucket_moisture_root'][start:end][mind]), s=70, color='red')
ax1.legend()
ax1.set_title('Maankosteuden aikasarja 2020')
ax1.set_ylabel(r'$\theta$ [m$^3$/m$^3$]')
#ax1.text(dates_spa[start:end][0], 0.2, 'a')
ax1.grid()

rasterio.plot.show(pk, transform=meta['transform'], ax=ax2);
#im1 = ax2.imshow(results_catch['bucket_moisture_root'][start:end][maxd], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im1 = rasterio.plot.show(snap1, transform=snap1.transform, ax=ax2, cmap='coolwarm_r', vmin=0.1, vmax=0.88, alpha=alp)
im1b = im1.get_images()[1]
ax2.text(384500, 7542000, f'{str(dates_spa[start:end][maxd])[0:10]}')
ax2.scatter(384400, 7542100, color='blue')
ax2.set_title('Esim. kostea päivä')

rasterio.plot.show(pk, transform=meta['transform'], ax=ax3);
#im2 = ax3.imshow(results_catch['bucket_moisture_root'][start:end][mind], cmap='coolwarm_r', vmin=0.1, vmax=0.88)
im2 = rasterio.plot.show(snap2, transform=snap2.transform, ax=ax3, cmap='coolwarm_r', vmin=0.1, vmax=0.88, alpha=alp)
im2b = im1.get_images()[1]
#fig.colorbar(im2b, ax=ax3, shrink=0.5)
ax3.text(384500, 7542000, f'{str(dates_spa[start:end][mind])[0:10]}')
ax3.scatter(384400, 7542100, color='red')
cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
fig.colorbar(im2b, cax=cbar_ax, label=r'$\theta$ [m$^3$m$^{-3}]$')
ax3.set_title('Esim. kuiva päivä')

ax3.axes.get_yaxis().set_ticklabels([])

plt.subplots_adjust(wspace=0.05, hspace=0.7)


plt.savefig(f'vesitalous_kosteus_{today}.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'vesitalous_kosteus_{today}.png', bbox_inches='tight', dpi=300)

#%%

im3 = ax4.imshow(results_catch['canopy_transpiration'][maxd], vmin=0, vmax=2.5)
fig.colorbar(im3, ax=ax3, shrink=0.5)

im4 = ax4.imshow(results_catch['canopy_transpiration'][mind], vmin=0, vmax=2.5)
fig.colorbar(im4, ax=ax4, shrink=0.5)

#%%
####################################################
# GIS PLOT
####################################################

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

im1 = ax1.imshow(results_2d['parameters_lai_conif'][20:250,20:165])
ax1.set_title('LAI conif')
fig.colorbar(im1, ax=ax1)

im2 = ax2.imshow(results_2d['parameters_lai_decid_max'][20:250,20:165])
ax2.set_title('LAI decid max')
fig.colorbar(im2, ax=ax2)

im3 = ax3.imshow(results_2d['parameters_lai_shrub'][20:250,20:165])
ax3.set_title('LAI shrub')
fig.colorbar(im3, ax=ax3)

im4 = ax4.imshow(results_2d['parameters_lai_grass'][20:250,20:165])
ax4.set_title('LAI grass')
fig.colorbar(im4, ax=ax4)

im5 = ax5.imshow(results_2d['parameters_hc'][20:250,20:165])
ax5.set_title('canopy height')
fig.colorbar(im5, ax=ax5)

im6 = ax6.imshow(results_2d['parameters_cf'][20:250,20:165], label='canopy fraction')
ax6.set_title('canopy fraction')
fig.colorbar(im6, ax=ax6)

cmapsoil = plt.get_cmap('viridis', 4)
im7 = ax7.imshow(results_2d['parameters_soilclass'][20:250,20:165], cmap=cmapsoil)
ax7.set_title('soilclass')
fig.colorbar(im7, ax=ax7)

cmapsite = plt.get_cmap('viridis', 4)
im8 = ax8.imshow(results_2d['parameters_sitetype'][20:250,20:165], cmap=cmapsite)
ax8.set_title('sitetype')
fig.colorbar(im8, ax=ax8)

cmapditch = plt.get_cmap('viridis', 2)
im9 = ax9.imshow(results_2d['parameters_ditches'][20:250,20:165], cmap=cmapditch)
ax9.set_title('streams/ditches')
cbar = fig.colorbar(im9, ax=ax9)
cbar.ax.locator_params(nbins=1)

im10 = ax10.imshow(results_2d['parameters_elevation'][20:250,20:165])
ax10.set_title('elevation')
fig.colorbar(im10, ax=ax10)

im11 = ax11.imshow(results_catch['parameters_twi'][20:250,20:165])
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

####################################################
# WATER BALANCE ANNUAL PLOT
####################################################

# yearly water balance in meters for whole area
wbdf = pd.DataFrame()
wbdf['P'] = results_2d['forcing_precipitation']
wbdf['SWE'] = np.nanmean(results_stand['canopy_snowfall'], axis=(1,2))
wbdf.index = dates_spa
wbdf['Qmod'] = results_catch['top_baseflow'] + np.nanmean(results_catch['bucket_surface_runoff'], axis=(1,2))
wbdf['ETmod'] = np.nanmean(results_catch['canopy_evaporation'] + results_catch['canopy_transpiration'] + results_catch['bucket_evaporation'], axis=(1,2))
#wbdf['ETdrymod'] = np.nanmean(results_catch['dry_et'], axis=(1,2))
wbdf['S'] = np.nanmean(results_catch['bucket_water_storage'], axis=(1,2)) + np.nanmean(results_catch['canopy_snow_water_equivalent'], axis=(1,2)) + np.nanmean(results_catch['canopy_water_storage'], axis=(1,2))
#wbdf['S'] = np.nanmean((results_catch['bucket_water_storage'] + results_catch['canopy_snow_water_equivalent'] + results_catch['canopy_water_storage'] + results_catch['bucket_pond_storage']), axis=(1,2))

wbdfy = pd.DataFrame()
wbdfy['P'] = wbdf['P'].resample('AS-SEP').sum()
wbdfy['SWE'] = wbdf['SWE'].resample('AS-SEP').sum()
wbdfy['Qmod'] = wbdf['Qmod'].resample('AS-SEP').sum()
#wbdfy['ETdrymod'] = wbdf['ETdrymod'].resample('AS-SEP').sum()
wbdfy['ETmod'] = wbdf['ETmod'].resample('AS-SEP').sum()
wbdfy['S'] = np.nan
wbdfy['Qobs'] = q.resample('AS-SEP').sum()
wbdfy['Qobs'].loc[wbdfy['Qobs'] < 220] = np.nan

wbdf2d = pd.DataFrame()
wbdf2d['P'] = results_2d['forcing_precipitation']
wbdf2d['SWE'] = np.nanmean(results_2d['canopy_snowfall'], axis=(1,2))
wbdf2d.index = dates_spa
wbdf2d['Qmod'] = np.nanmean(results_2d['soil_netflow_to_ditch'] + results_2d['bucket_surface_runoff'], axis=(1,2))
#wbdf2d['ETdrymod'] = wbdf['ETdrymod'].resample('Y').sum()
wbdf2d['ETmod'] = np.nanmean(results_2d['canopy_evaporation'] + results_2d['canopy_transpiration'] + results_2d['bucket_evaporation'], axis=(1,2))
wbdf2d['S'] = np.nanmean(results_2d['bucket_water_storage'] + results_2d['canopy_snow_water_equivalent'] + results_2d['canopy_water_storage'] + results_2d['soil_water_storage'], axis=(1,2))

wbdf2dy = pd.DataFrame()
wbdf2dy['P'] = wbdf2d['P'].resample('AS-SEP').sum()
wbdf2dy['SWE'] = wbdf2d['SWE'].resample('AS-SEP').sum()
wbdf2dy['Qmod'] = wbdf2d['Qmod'].resample('AS-SEP').sum()
wbdf2dy['ETmod'] = wbdf2d['ETmod'].resample('AS-SEP').sum()
#wbdf2dy['ETdrymod'] = wbdf2d['ETdrymod'].resample('AS-SEP').sum()
wbdf2dy['S'] = np.nan
wbdf2dy['Qobs'] = q.resample('AS-SEP').sum()
wbdf2dy['Qobs'].loc[wbdf2dy['Qobs'] < 220] = np.nan

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
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(9,5));
ax1 = axs[0]
ax2 = axs[1]

for i in range(len(wbdfy)):
    ax1.bar(str(wbdfy.index[i].year), wbdfy['P'][i], alpha=0.8, color='tab:blue', label=r'P$_{obs}$')
    #ax1.bar(str(wbdfy.index[i].year), wbdfy['P'][i] - wbdfy['SWE'][i], alpha=0.8, color='tab:blue', label=r'SWE$_{mod}$')
    ax1.bar(str(wbdfy.index[i].year), - wbdfy['ETmod'][i] - wbdfy['Qmod'][i], color='tab:green', alpha=0.7, label=r'ET$_{mod}$')
    ax1.bar(str(wbdfy.index[i].year), - wbdfy['Qmod'][i], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
    ax1.scatter(str(wbdfy.index[i].year), wbdfy['SWE'][i], s=100, marker='_', color='b', linewidth=2, zorder=2, label=r'$SWE_{mod}$')
    ax1.scatter(str(wbdfy.index[i].year), - wbdfy['Qobs'][i], s=100, marker='x', color='k', linewidth=3, zorder=2, label=r'$Q_{obs}$')
    ax1.set_ylim(-1000, 1000)
    if i == 0:
        ax1.legend(ncol=5,bbox_to_anchor=(1.8, 1.2))

    ax2.bar(str(wbdf2dy.index[i].year), wbdf2dy['P'][i], alpha=0.8, color='tab:blue', label=r'P$_{obs}$')
    #ax2.bar(str(wbdf2dy.index[i].year), wbdf2dy['P'][i] - wbdf2dy['SWE'][i], alpha=0.8, color='tab:blue', label=r'SWE$_{mod}$')
    ax2.bar(str(wbdf2dy.index[i].year), - wbdf2dy['ETmod'][i] - wbdf2dy['Qmod'][i], color='tab:green', alpha=0.7, label=r'ET$_{mod}$')
    ax2.bar(str(wbdf2dy.index[i].year), - wbdf2dy['Qmod'][i], alpha=0.7, color='tab:brown', label=r'Q$_{mod}$')
    ax2.scatter(str(wbdf2dy.index[i].year), wbdf2dy['SWE'][i], s=100, marker='_', color='b', linewidth=2, zorder=2, label=r'$SWE_{mod}$')
    ax2.scatter(str(wbdf2dy.index[i].year), - wbdf2dy['Qobs'][i], s=100, marker='x', color='k', linewidth=3, zorder=2, label=r'$Q_{obs}$')
    ax2.set_ylim(-1000, 1000)

ax1.xaxis.set_tick_params(rotation=45)
ax2.xaxis.set_tick_params(rotation=45)

ax1.grid(); ax2.grid()
ax1.set_ylabel('mm / year')
ax1.set_title('SpaFHy-TOP')
ax2.set_title('SpaFHy-2D')

if saveplots == True:
        plt.savefig(f'WB_BARPLOT_{today}.pdf',bbox_inches='tight')
        plt.savefig(f'WB_BARPLOT_{today}.png',bbox_inches='tight')


#%%

####################################################
# SWE AND RUNOFF PLOT
####################################################

# SWE file
fn = r'C:\SpaFHy_v1_Pallas\data\obs\SWE_survey_2018-02-22_2021-05-16.txt'
SWE_m = pd.read_csv(fn, skiprows=5, sep=';', parse_dates = ['date'], encoding='iso-8859-1')
SWE_m = SWE_m.rename(columns={'date':'time'})
SWE_m.index = SWE_m['time']
#SWE_m['SWE'].loc[SWE_m['date'] < '2018-07-07'] = np.nan
SWE_m = SWE_m[['SWE', 'SWE_sd', 'quality']]
SWE_m = SWE_m[~SWE_m.index.duplicated()]

SWE = pd.DataFrame()
SWE['mod_mean'] = np.nanmean(results_stand['canopy_snow_water_equivalent'], axis=(1,2))
#SWE['mod_iq25'] = np.nanquantile(results_stand['canopy_snow_water_equivalent'], 0.1, axis=(1,2))
#SWE['mod_iq75'] = np.nanquantile(results_stand['canopy_snow_water_equivalent'], 0.9, axis=(1,2))
SWE['mod_iq25'] = np.nanmin(results_stand['canopy_snow_water_equivalent'], axis=(1,2))
SWE['mod_iq75'] = np.nanmax(results_stand['canopy_snow_water_equivalent'], axis=(1,2))
SWE['time'] = dates_spa
SWE.index = SWE['time']
SWE = SWE.drop(columns='time')

SWE['obs'] = np.nan
SWE['obs'] = SWE_m['SWE']
SWE['obs_sd'] = SWE_m['SWE_sd']

SWE = SWE[SWE[SWE['obs'].notna()].index[0]:]


# Q file
q_all = pd.DataFrame()
q_all['2D'] = np.nanmean(results_2d['soil_netflow_to_ditch'], axis=(1,2)) + np.nanmean(results_2d['bucket_surface_runoff'], axis=(1,2))
q_all['time'] = dates_spa
q_all.index = q_all['time']
q_all = q_all.drop(columns='time')
q_all['catch'] = results_catch['top_baseflow'] + np.nanmean(results_catch['bucket_surface_runoff'], axis=(1,2))
q_all['obs'] = q
q_all['stand'] = np.nanmean(results_stand['bucket_drainage'], axis=(1,2)) + np.nanmean(results_stand['bucket_surface_runoff'], axis=(1,2))

# these period lack manual confirming observation (tarkistusmittaus)
q_all['obs'].loc['2017-08-15':'2018-09-21'] = np.nan
q_all['obs'].loc['2018-12-22':'2019-04-09'] = np.nan

#q_all['mae_2d'] = abs(q_all['obs'] - q_all['2D'])
q_all['mae_catch'] = abs(q_all['obs'] - q_all['catch'])
q_all = q_all[q_all[q_all['obs'].notna()].index[0]:]
#mae_2d = round(np.nanmean(q_all['mae_2d']), 2)
#mae_catch = round(np.nanmean(q_all['mae_catch']), 2)

q_all_sc = q_all.dropna()
sampl2 = np.random.randint(low=0, high=len(q_all_sc), size=(1000,))
q_all_sc['time'] = q_all_sc.index
q_all_sc = q_all_sc.reset_index(drop=True)
q_all_sc = q_all_sc.iloc[sampl2]
q_all_sc.index = q_all_sc['time']
q_all_sc = q_all_sc.drop(columns='time')


### PLOT
fig4 = plt.figure(constrained_layout=False, figsize=(10,8))
gs = fig4.add_gridspec(3, 4)

ax3 = fig4.add_subplot(gs[0, :3])
ax1 = fig4.add_subplot(gs[1, :3])
ax2 = fig4.add_subplot(gs[2, :3])

ax6 = fig4.add_subplot(gs[0, 3])
ax4 = fig4.add_subplot(gs[1, 3])
ax5 = fig4.add_subplot(gs[2, 3])

ax4.yaxis.tick_right()
ax5.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax5.yaxis.set_label_position("right")

l1 = ax1.plot(q_all['catch'], label = r'Q$_{TOP}$')
l2 = ax1.plot(q_all['obs'], alpha=0.7, label=r'Q$_{obs}$')
ax1.set(ylim=(-1, 23))
ax1.set_ylabel(r'Qf (mm d$^{-1}$)')


l3 = ax2.plot(q_all['2D'], color='tab:green', label=r'Q$_{2D}$')
ax2.plot(q_all['obs'], color='tab:orange', label=r'Q$_{obs}$', alpha=0.7)
ax2.set(ylim=(-1, 29))
ax2.set_ylabel(r'Qf (mm d$^{-1}$)')

l4 = ax3.errorbar(SWE.index, SWE['obs'], SWE['obs_sd'],
                  color = 'black', alpha=0.5, linestyle='None', markersize=3, marker='o',
                  label=r'SWE$_{obs}$mean w/ std')
ax3.plot(SWE['mod_mean'], label=r'SWE$_{mod}$mean')
ax3.legend(ncol=5, bbox_to_anchor=(0.6, 1.22))
ax3.set_ylabel('SWE (mm)')
#ax3.set(ylim=(-1, 400))


l5 = sns.regplot(ax=ax6, x=SWE['obs'], y=SWE['mod_mean'], scatter_kws={'s':50, 'alpha':0.4}, line_kws={"color": "red"})
ax6.set_ylabel(r'SWE$_{mod}$')
ax6.set_xlabel(r'SWE$_{obs}$')
ax6.yaxis.tick_right()
ax6.yaxis.set_label_position("right")

ax1.legend(ncol=2,loc='upper left')

#ax2.legend(loc='upper left', framealpha=1)
ax2.legend(ncol=2, loc='upper left')

sns.regplot(ax=ax4, x=q_all_sc['obs'], y=q_all_sc['catch'], scatter_kws={'s':50, 'alpha':0.15}, line_kws={"color": "red"})
sns.regplot(ax=ax5, x=q_all_sc['obs'], y=q_all_sc['2D'], scatter_kws={'s':50, 'alpha':0.15}, line_kws={"color": "red"})
ax4.set_ylabel(r'Qf$_{mod}$ mm d$^{-1}$')
ax4.set_xlabel(r'Qf$_{obs}$ mm d$^{-1}$')
ax5.set_ylabel(r'Qf$_{mod}$ mm d$^{-1}$')
ax5.set_xlabel(r'Qf$_{obs}$ mm d$^{-1}$')
ax4.set(ylim=(0, 15))
ax4.set(xlim=(0, 15))
ax5.set(ylim=(0, 15))
ax5.set(xlim=(0, 15))

ax1.grid(); ax2.grid()
ax3.grid(); ax4.grid()
ax5.grid(); ax6.grid()

plt.tight_layout()
#plt.subplots_adjust(wspace=0, hspace=0)


if saveplots == True:
        plt.savefig(f'QF_MOD_OBS_{today}.pdf')
        plt.savefig(f'QF_MOD_OBS_{today}.png')

#%%

####################################################
# ET PLOT
####################################################

# DRY ET DATA PREPARING
# et sim. obs
# ET obs
file = r'C:\SpaFHy_v1_Pallas_2D\obs\ec_et.csv'
ec = pd.read_csv(file, sep=';', index_col=0, parse_dates=True)
ec = ec[ec.index.isin(dates_spa)]
ec['k_obs'].iloc[(ec.index < '2017-01-01') & (ec.index > '2016-01-01')] = np.nan
ec['k_obs'].loc[ec['k_obs'] < 0] = np.nan
ec['l_obs'].loc[ec['l_obs'] < 0] = np.nan

#dry et
results_2d['dry_et'] = results_2d['bucket_evaporation'] + results_2d['canopy_transpiration'] #+ results_2d['canopy_evaporation']
#results_2d['dry_et'][ix_p,:,:] = np.nan
results_stand['dry_et'] = results_stand['bucket_evaporation'] + results_stand['canopy_transpiration'] #+ results_stand['canopy_evaporation']
#results_stand['dry_et'][ix_p,:,:] = np.nan
results_catch['dry_et'] = results_catch['bucket_evaporation'] + results_catch['canopy_transpiration'] #+ results_catch['canopy_evaporation']
#catch_dry_et[ix_p,:,:] = np.nan

dry_et = pd.DataFrame()
dry_et['k_mod'] = results_2d['bucket_evaporation'][:,k_loc[0],k_loc[0]] + results_2d['canopy_transpiration'][:,k_loc[0],k_loc[0]]
dry_et['l_mod'] = results_2d['bucket_evaporation'][:,l_loc[0],l_loc[0]] + results_2d['canopy_transpiration'][:,l_loc[0],l_loc[0]]
dry_et.index = dates_spa

dry_et['k_obs'] = ec['k_obs']
dry_et['l_obs'] = ec['l_obs']
dry_et[dry_et.index.isin(wetdays)] = np.nan


fig3 = plt.figure(constrained_layout=True, figsize=(10,5))
gs = fig3.add_gridspec(2, 4)

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Kenttärova')
f3_ax1.plot(dry_et['k_mod'], alpha=1.0, label='mod')
#f3_ax1.plot(dates_spa, results_stand['forcing_vapor_pressure_deficit'])
f3_ax1.plot(dry_et['k_obs'], 'k.',  alpha=0.3,markersize=6, label='obs')
f3_ax1.legend(loc='upper left')
f3_ax1.set(ylim=(-0.5,  8.2))
f3_ax1.set_ylabel(r'ET mm d$^{-1}$')

f3_ax2 = fig3.add_subplot(gs[0, 3])
x2 = sns.regplot(ax=f3_ax2, x=dry_et['k_mod'], y=dry_et['k_obs'], scatter_kws={'s':30, 'alpha':0.15}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0,  8.2))
f3_ax2.set(xlim=(0,  8.2))
#f3_ax2.axes.get_xaxis().set_visible(False)
f3_ax2.yaxis.tick_right()
f3_ax2.set_ylabel(r'ET$_{obs}$ (mm d$^{-1}$)')
f3_ax2.set_xlabel(r'ET$_{mod}$ (mm d$^{-1}$)')


# Plotting Lompolonjänkkä ET mod vs. obs

f3_ax3 = fig3.add_subplot(gs[1, :3])
f3_ax3.set_title('Lompolonjänkkä', y=-0.3)
f3_ax3.plot(dry_et['l_mod'], alpha=1.0, label='mod')
f3_ax3.plot(dry_et['l_obs'], 'k.', alpha=0.3,markersize=6, label='obs')
f3_ax3.legend(loc='upper left')
f3_ax3.set(ylim=(-0.5, 8.2))
f3_ax3.set_ylabel(r'ET mm d$^{-1}$')

f3_ax4 = fig3.add_subplot(gs[1, 3])
x4 = sns.regplot(ax=f3_ax4, x=dry_et['l_mod'], y=dry_et['l_obs'], scatter_kws={'s':30, 'alpha':0.15}, line_kws={"color": "red"})
f3_ax4.set(ylim=(0,  8.2))
f3_ax4.set(xlim=(0,  8.2))
f3_ax4.yaxis.tick_right()
f3_ax4.set_ylabel(r'ET$_{obs}$ (mm d$^{-1}$)')
f3_ax4.set_xlabel(r'ET$_{mod}$ (mm d$^{-1}$)')

f3_ax1.grid(); f3_ax2.grid(); f3_ax3.grid(); f3_ax4.grid();

if saveplots == True:
        plt.savefig(f'ET_MOD_OBS_KR_LV_{today}.pdf', bbox_inches='tight')
        plt.savefig(f'ET_MOD_OBS_KR_LV_{today}.png', bbox_inches='tight')


#%%

####################################################
# TEMPORAL SOIL MOISTURE AND GROUNDWATER PLOT
####################################################

# preparing soil moisture datas
# point examples from mineral and openmire
# soilscouts at Kenttarova
file = r'C:\SpaFHy_v1_Pallas_2D\obs\soilm_kenttarova.csv'
soilm = pd.read_csv(file, sep=';', index_col=0, parse_dates=True)
lessthantwo = soilm.index[soilm.notna().sum(axis=1) < 2]
soilm['mean_obs'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A', 's3', 'SH-5B']].mean(numeric_only=True, axis=1)
soilm['mean_obs'][lessthantwo] = np.nan
#soilm['iq25_obs'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A', 's3', 'SH-5B']].quantile(q=0.25, numeric_only=True, axis=1)
#soilm['iq75_obs'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A', 's3', 'SH-5B']].quantile(q=0.75, numeric_only=True, axis=1)
soilm['min'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A', 's3', 'SH-5B']].min(numeric_only=True, axis=1)
soilm['max'] = soilm[['SH-20A', 'SH-20B', 's18', 'SH-5A', 's3', 'SH-5B']].max(numeric_only=True, axis=1)
#soilm['iq25_obs'][lessthantwo] = np.nan
#soilm['iq75_obs'][lessthantwo] = np.nan
soilm['min'][lessthantwo] = np.nan
soilm['max'][lessthantwo] = np.nan

#
gw_file = r'C:\SpaFHy_v1_Pallas_2D\obs\pallas_gw_levels.csv'
gw = pd.read_csv(gw_file, sep=';', index_col=0, parse_dates=True)
gwdf = pd.DataFrame()
gwdf['mod_kr'] = results_2d['soil_ground_water_level'][:,k_loc[0], k_loc[1]]
gwdf.index = dates_spa
gwdf['mod_lv'] = results_2d['soil_ground_water_level'][:,l_loc[0], l_loc[1]]
gwdf['obs_lv'] = gw['pzp01']
gwdf['obs_kr'] = gw['pvp11']

spa_wliq_df = pd.DataFrame()
# 2D
spa_wliq_df['spa_k_2d_root'] = results_2d['bucket_moisture_root'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_2d_root'] = results_2d['bucket_moisture_root'][:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_2d_top'] = results_2d['bucket_moisture_top'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_2d_top'] = results_2d['bucket_moisture_top'][:,l_loc[0],l_loc[1]]
# stand
spa_wliq_df['spa_k_st_root'] = results_stand['bucket_moisture_root'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_st_root'] = results_stand['bucket_moisture_root'][:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_st_top'] = results_stand['bucket_moisture_top'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_st_top'] = results_stand['bucket_moisture_top'][:,l_loc[0],l_loc[1]]
# catch
spa_wliq_df['spa_k_ca_root'] = results_catch['bucket_moisture_root'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_root'] = results_catch['bucket_moisture_root'][:,l_loc[0],l_loc[1]]
spa_wliq_df['spa_k_ca_top'] = results_catch['bucket_moisture_top'][:,k_loc[0],k_loc[1]]
spa_wliq_df['spa_l_ca_top'] = results_catch['bucket_moisture_top'][:,l_loc[0],l_loc[1]]
spa_wliq_df.index = dates_spa
# sar
sar_wliq = pd.DataFrame()
sar_wliq['sar_k'] = sar['soilmoisture'][:,k_loc[0], k_loc[1]]/100
sar_wliq['sar_l'] = sar['soilmoisture'][:, l_loc[0], l_loc[1]]/100
sar_wliq.index = dates_sar

for i in spa_wliq_df.columns:
    soilm[i] = spa_wliq_df[i]
for i in sar_wliq.columns:
    soilm[i] = sar_wliq[i]

soilm['month'] = pd.DatetimeIndex(soilm.index).month
winter_ix = np.where((soilm['month'] <= 4) | (soilm['month'] >= 11))[0]
soilm.iloc[winter_ix] = np.nan
gwdf.iloc[winter_ix] = np.nan

# sampling for the scatterplot
soilmsc = soilm[np.isfinite(soilm['mean_obs'])]
soilmsc['time'] = soilmsc.index
soilmsc = soilmsc.reset_index(drop=True)

sampl2 = np.random.randint(low=0, high=len(soilmsc), size=(1000,))
soilmsc = soilmsc.iloc[sampl2]
soilmsc.index = soilmsc['time']

soilmsummer = soilm.loc[(soilm['month'] == 5) | (soilm['month'] == 6) | (soilm['month'] == 7) | (soilm['month'] == 8) | (soilm['month'] == 9) | (soilm['month'] == 10)]
soilmsummer['time'] = soilmsummer.index
soilsummer = soilmsummer.reset_index(drop=True)



#%%


###
fig3 = plt.figure(constrained_layout=True, figsize=(12,3))
gs = fig3.add_gridspec(1, 4)
sns.set_style('whitegrid')

f3_ax1 = fig3.add_subplot(gs[0, :3])
im1 = f3_ax1.plot(soilm['spa_k_2d_root'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax1.plot(soilm['spa_k_ca_root'], 'red', alpha=0.6, label='stand/catch')
f3_ax1.plot(soilm['mean_obs'],  'k', alpha=0.5, label=r'mean$_{obs}$')
#f3_ax1.fill_between(soilm.index, soilm['iq25_obs'], soilm['iq75_obs'], color='blue', alpha=0.2, label=r'IQR$_{obs}$')
f3_ax1.fill_between(soilm.index, soilm['min'], soilm['max'], color='blue', alpha=0.2, label=r'range$_{obs}$')

f3_ax1.legend(ncol=5,bbox_to_anchor=(0.7, 1.2))
y = f3_ax1.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#y.set_rotation(0)
f3_ax1.set_ylim(0.1,0.5)
#f3_ax1.axes.get_xaxis().set_visible(False)


f3_ax2 = fig3.add_subplot(gs[0, 3])
x4 = sns.regplot(ax=f3_ax2, x=soilmsc['spa_k_st_root'], y=soilmsc['mean_obs'], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
f3_ax2.set(ylim=(0.1, 0.45))
f3_ax2.set(xlim=(0.1, 0.45))
f3_ax2.yaxis.tick_right()
f3_ax2.set_ylabel(r'$\theta_{obs}$ (m$^3$m$^{-3}$)')
f3_ax2.set_xlabel(r'$\theta_{mod}$ (m$^3$m$^{-3}$)')

if saveplots == True:
        plt.savefig(f'THETA_MOD_OBS_KR_{today}.pdf', bbox_inches='tight')
        plt.savefig(f'THETA_MOD_OBS_KR_{today}.png', bbox_inches='tight')

# %%

fig3 = plt.figure(constrained_layout=True, figsize=(12,12))
gs = fig3.add_gridspec(7, 4)
sns.set_style('whitegrid')

f3_ax1 = fig3.add_subplot(gs[1, :2])
im1 = f3_ax1.plot(soilm['spa_k_2d_root']['2011-05-01':'2011-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax1.plot(soilm['mean_obs']['2011-05-01':'2011-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
#f3_ax1.plot(soilm['sar_k']['2011-05-01':'2011-10-01'], 'red', alpha=0.7, label='SAR')
f3_ax1.fill_between(pd.date_range('2011-05-01','2011-10-01', freq='D'), soilm['min']['2011-05-01':'2011-10-01'],
                    soilm['max']['2011-05-01':'2011-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')

f3_ax1.legend(ncol=4,bbox_to_anchor=(1.0, 1.3))


f3_ax2 = fig3.add_subplot(gs[1, 2:4])
im2 = f3_ax2.plot(soilm['spa_k_2d_root']['2012-05-01':'2012-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax2.plot(soilm['mean_obs']['2012-05-01':'2012-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax2.fill_between(pd.date_range('2012-05-01','2012-10-01', freq='D'), soilm['min']['2012-05-01':'2012-10-01'],
                    soilm['max']['2012-05-01':'2012-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax3 = fig3.add_subplot(gs[2, :2])
im3 = f3_ax3.plot(soilm['spa_k_2d_root']['2013-05-01':'2013-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax3.plot(soilm['mean_obs']['2013-05-01':'2013-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax3.fill_between(pd.date_range('2013-05-01','2013-10-01', freq='D'), soilm['min']['2013-05-01':'2013-10-01'],
                    soilm['max']['2013-05-01':'2013-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax4 = fig3.add_subplot(gs[2, 2:4])
im4 = f3_ax4.plot(soilm['spa_k_2d_root']['2014-05-01':'2014-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax4.plot(soilm['mean_obs']['2014-05-01':'2014-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax4.fill_between(pd.date_range('2014-05-01','2014-10-01', freq='D'), soilm['min']['2014-05-01':'2014-10-01'],
                    soilm['max']['2014-05-01':'2014-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax5 = fig3.add_subplot(gs[3, :2])
im5 = f3_ax5.plot(soilm['spa_k_2d_root']['2015-05-01':'2015-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax5.plot(soilm['mean_obs']['2015-05-01':'2015-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax5.fill_between(pd.date_range('2015-05-01','2015-10-01', freq='D'), soilm['min']['2015-05-01':'2015-10-01'],
                    soilm['max']['2015-05-01':'2015-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax6 = fig3.add_subplot(gs[3, 2:4])
im6 = f3_ax6.plot(soilm['spa_k_2d_root']['2016-05-01':'2016-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax6.plot(soilm['mean_obs']['2016-05-01':'2016-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax6.fill_between(pd.date_range('2016-05-01','2016-10-01', freq='D'), soilm['min']['2016-05-01':'2016-10-01'],
                    soilm['max']['2016-05-01':'2016-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax7 = fig3.add_subplot(gs[4, :2])
im7 = f3_ax7.plot(soilm['spa_k_2d_root']['2017-05-01':'2017-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax7.plot(soilm['mean_obs']['2017-05-01':'2017-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax7.fill_between(pd.date_range('2017-05-01','2017-10-01', freq='D'), soilm['min']['2017-05-01':'2017-10-01'],
                    soilm['max']['2017-05-01':'2017-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax8 = fig3.add_subplot(gs[4, 2:4])
im8 = f3_ax8.plot(soilm['spa_k_2d_root']['2018-05-01':'2018-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax8.plot(soilm['mean_obs']['2018-05-01':'2018-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax8.fill_between(pd.date_range('2018-05-01','2018-10-01', freq='D'), soilm['min']['2018-05-01':'2018-10-01'],
                    soilm['max']['2018-05-01':'2018-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')


f3_ax9 = fig3.add_subplot(gs[5, :2])
im9 = f3_ax9.plot(soilm['spa_k_2d_root']['2019-05-01':'2019-10-01'], 'g', alpha=0.7, label='2D')
f3_ax9.plot(soilm['sar_k']['2019-05-01':'2019-10-01'], 'red', alpha=0.7, label='SAR')
f3_ax9.plot(soilm['mean_obs']['2019-05-01':'2019-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax9.fill_between(pd.date_range('2019-05-01','2019-10-01', freq='D'), soilm['min']['2019-05-01':'2019-10-01'],
                    soilm['max']['2019-05-01':'2019-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')
#f3_ax9.legend(ncol=4)


f3_ax10 = fig3.add_subplot(gs[5, 2:4])
im10 = f3_ax10.plot(soilm['spa_k_2d_root']['2020-05-01':'2020-10-01'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax10.plot(soilm['mean_obs']['2020-05-01':'2020-10-01'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax10.fill_between(pd.date_range('2020-05-01','2020-10-01', freq='D'), soilm['min']['2020-05-01':'2020-10-01'],
                    soilm['max']['2020-05-01':'2020-10-01'], color='blue', alpha=0.2, label=r'range$_{obs}$')

f3_ax11 = fig3.add_subplot(gs[6, :2])
im11 = f3_ax11.plot(soilm['spa_k_2d_root']['2021-05-01':'2021-09-06'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax11.plot(soilm['mean_obs']['2021-05-01':'2021-09-06'],  'k', alpha=0.5, label=r'mean$_{obs}$')
f3_ax11.fill_between(pd.date_range('2021-05-01','2021-09-06', freq='D'), soilm['min']['2021-05-01':'2021-09-06'],
                    soilm['max']['2021-05-01':'2021-09-06'], color='blue', alpha=0.2, label=r'range$_{obs}$')


if saveplots == True:
        plt.savefig(f'THETA_MOD_OBS_KR_ANNUAL_{today}.pdf', bbox_inches='tight')
        plt.savefig(f'THETA_MOD_OBS_KR_ANNUAL_{today}.png',bbox_inches='tight')

'''
f3_ax3 = fig3.add_subplot(gs[1, :3])
f3_ax3.set_title('Water level')
im2 = f3_ax3.plot(gwdf['mod_kr'],  'g', alpha=0.7, label='2D')
f3_ax3.plot(gwdf['obs_kr'], 'red', alpha=0.6, label='stand/catch')
#ax2.title.set_text('Mire')
f3_ax3.set_ylim(-2.5,0.5)
f3_ax3.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#f3_ax3.legend(ncol=5)
'''

#%%

###
fig3 = plt.figure(constrained_layout=True, figsize=(14,7))
gs = fig3.add_gridspec(2, 4)
sns.set_style('whitegrid')

f3_ax1 = fig3.add_subplot(gs[0, :3])
f3_ax1.set_title('Soil moisture')
im1 = f3_ax1.plot(soilm['spa_l_2d_root'], 'g', alpha=0.7, label='2D')
#f3_ax1.plot(soilm['spa_k_st_root'], 'blue', alpha=0.7, label='stand')
f3_ax1.plot(soilm['spa_l_ca_root'], 'red', alpha=0.6, label='stand/catch')


f3_ax1.legend(ncol=5,bbox_to_anchor=(0.8, 1.3))
y = f3_ax1.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#y.set_rotation(0)
f3_ax1.set_ylim(0.4,1.0)
#f3_ax1.axes.get_xaxis().set_visible(False)


f3_ax3 = fig3.add_subplot(gs[1, :3])
f3_ax3.set_title('Water level')
im2 = f3_ax3.plot(gwdf['mod_lv'],  'g', alpha=0.7, label='2D')
f3_ax3.plot(gwdf['obs_lv'], 'red', alpha=0.6, label='stand/catch')
#ax2.title.set_text('Mire')
f3_ax3.set_ylim(-0.2,0.2)
f3_ax3.set_ylabel(r'$\theta$ m$^3$m$^{-3}$')
#f3_ax3.legend(ncol=5)

