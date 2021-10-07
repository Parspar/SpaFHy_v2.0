# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:51:09 2021

@author: janousu
"""


from __future__ import print_function
print(__doc__)

import numpy as np
import sys,csv

#sys.path.append('./py_files')

#from SPAEF_metric import SPAEF
#from SPAEF_figures import plot_SPAEFstats, plot_maps

from iotools import read_results
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.gridspec as gridspec
from datetime import date
import seaborn as sns
from scipy.stats import variation,zscore
import math

today = date.today()


def spaef_metrics(data1, data2, bins, histrange):
    hobs,binobs = np.histogram(data1,bins,histrange)
    hsim,binsim = np.histogram(data2,bins,histrange)
    #convert int to float, critical conversion for the result
    hobs=np.float64(hobs)
    hsim=np.float64(hsim)
    #find the overlapping of two histogram      
    minima = np.minimum(hsim, hobs)
    #compute the fraction of intersection area to the observed histogram area, hist intersection/overlap index   
    hm = round(np.sum(minima)/np.sum(hobs),2)
    # pearson coef
    cc = round(np.corrcoef(data1[np.isfinite(data1)], data2[np.isfinite(data2)])[0,1],2)
    # variation
    cv = round(variation(data2[np.isfinite(data2)])/variation(data1[np.isfinite(data1)]),2)
    spaef = round(1 - (math.sqrt((cc - 1)**2 + (cv - 1)**2 + (hm - 1)**2)),2)
    return hm, cc, cv, spaef


# reading the stand results
outputfile_stand = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_catch.nc'
results_stand = read_results(outputfile_stand)

# reading the stand results
outputfile_2d = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_2d.nc'
results_2d = read_results(outputfile_2d)

# reading the catch results
outputfile_catch = 'C:\SpaFHy_v1_Pallas_2D/results/testcase_input_catch.nc'
results_catch = read_results(outputfile_catch)

ditches = results_2d['parameters_ditches']
ditch_mask = np.array(ditches)
ditch_mask[ditch_mask > -1] = 1
ditch_mask[ditches == -1] = np.nan

results_2d['bucket_moisture_root'] = results_2d['bucket_moisture_root'] * ditch_mask
results_catch['bucket_moisture_root'] = results_catch['bucket_moisture_root'] * ditch_mask
results_stand['bucket_moisture_root'] = results_stand['bucket_moisture_root'] * ditch_mask

dates_spa = []
for d in range(len(results_stand['date'])):
    dates_spa.append(pd.to_datetime(str(results_stand['date'][d])[36:46]))
    
start = np.where(pd.to_datetime(dates_spa) == '2021-08-01')[0][0]
end = np.where(pd.to_datetime(dates_spa) == '2021-09-01')[0][0]

soilclass = np.array(results_stand['parameters_soilclass'])
soilclass_copy = np.array(results_stand['parameters_soilclass'])
soilclass_copy2 = np.array(results_stand['parameters_soilclass'])
soilclass_2 = np.ravel(soilclass)
soilclass_4 = np.ravel(soilclass_copy)
soilclass_3 = np.ravel(soilclass_copy2)
r, c = np.shape(soilclass)
soilclass_2[soilclass_2 != 2] = np.nan
soilclass_2[soilclass_2 == 2] = 1
soilclass_4[soilclass_4 != 4] = np.nan
soilclass_4[soilclass_4 == 4] = 1
soilclass_3[soilclass_3 != 3] = np.nan
soilclass_3[soilclass_3 == 3] = 1
soilclass_2 = soilclass_2.reshape(r, c)
soilclass_4 = soilclass_4.reshape(r, c)
soilclass_3 = soilclass_3.reshape(r, c)


#%%

# plotting histogram and scatterplot on my own

may_ind = np.where(pd.to_datetime(dates_spa).month == 5)[0]
jun_ind = np.where(pd.to_datetime(dates_spa).month == 6)[0]
jul_ind = np.where(pd.to_datetime(dates_spa).month == 7)[0]
aug_ind = np.where(pd.to_datetime(dates_spa).month == 8)[0]
sep_ind = np.where(pd.to_datetime(dates_spa).month == 9)[0]

soilclasses = [soilclass_2, soilclass_4]

new_arr = np.zeros((2,4))


bins =100# np.linspace(100, 350, 100)
histrange = [0,1]
#hobs,binobs = np.histogram(obs,bins, histrange)
#hsim,binsim = np.histogram(sim,bins, histrange)

td_may2 = np.array(results_2d['bucket_moisture_root'][may_ind,:,:] * soilclass_2).flatten()
td_may2 = td_may2[~np.isnan(td_may2)]
td_may4 = np.array(results_2d['bucket_moisture_root'][may_ind,:,:] * soilclass_4).flatten()
td_may4 = td_may4[~np.isnan(td_may4)]

td_june2 = np.array(results_2d['bucket_moisture_root'][jun_ind,:,:] * soilclass_2).flatten()
td_june2 = td_june2[~np.isnan(td_june2)]
td_june4 = np.array(results_2d['bucket_moisture_root'][jun_ind,:,:] * soilclass_4).flatten()
td_june4 = td_june4[~np.isnan(td_june4)]

td_july2 = np.array(results_2d['bucket_moisture_root'][jul_ind,:,:] * soilclass_2).flatten()
td_july2 = td_july2[~np.isnan(td_july2)]
td_july4 = np.array(results_2d['bucket_moisture_root'][jul_ind,:,:] * soilclass_4).flatten()
td_july4 = td_july4[~np.isnan(td_july4)]

td_august2 = np.array(results_2d['bucket_moisture_root'][aug_ind,:,:] * soilclass_2).flatten()
td_august2 = td_august2[~np.isnan(td_august2)]
td_august4 = np.array(results_2d['bucket_moisture_root'][aug_ind,:,:] * soilclass_4).flatten()
td_august4 = td_august4[~np.isnan(td_august4)]

td_sept2 = np.array(results_2d['bucket_moisture_root'][sep_ind,:,:] * soilclass_2).flatten()
td_sept2 = td_sept2[~np.isnan(td_sept2)]
td_sept4 = np.array(results_2d['bucket_moisture_root'][sep_ind,:,:] * soilclass_4).flatten()
td_sept4 = td_sept4[~np.isnan(td_sept4)]

st_may2 = np.array(results_stand['bucket_moisture_root'][may_ind,:,:] * soilclass_2).flatten()
st_may2 = st_may2[~np.isnan(st_may2)]
st_may4 = np.array(results_stand['bucket_moisture_root'][may_ind,:,:] * soilclass_4).flatten()
st_may4 = st_may4[~np.isnan(st_may4)]

st_june2 = np.array(results_stand['bucket_moisture_root'][jun_ind,:,:] * soilclass_2).flatten()
st_june2 = st_june2[~np.isnan(st_june2)]
st_june4 = np.array(results_stand['bucket_moisture_root'][jun_ind,:,:] * soilclass_4).flatten()
st_june4 = st_june4[~np.isnan(st_june4)]

st_july2 = np.array(results_stand['bucket_moisture_root'][jul_ind,:,:] * soilclass_2).flatten()
st_july2 = st_july2[~np.isnan(st_july2)]
st_july4 = np.array(results_stand['bucket_moisture_root'][jul_ind,:,:] * soilclass_4).flatten()
st_july4 = st_july4[~np.isnan(st_july4)]

st_august2 = np.array(results_stand['bucket_moisture_root'][aug_ind,:,:] * soilclass_2).flatten()
st_august2 = st_august2[~np.isnan(st_august2)]
st_august4 = np.array(results_stand['bucket_moisture_root'][aug_ind,:,:] * soilclass_4).flatten()
st_august4 = st_august4[~np.isnan(st_august4)]

st_sept2 = np.array(results_stand['bucket_moisture_root'][sep_ind,:,:] * soilclass_2).flatten()
st_sept2 = st_sept2[~np.isnan(st_sept2)]
st_sept4 = np.array(results_stand['bucket_moisture_root'][sep_ind,:,:] * soilclass_4).flatten()
st_sept4 = st_sept4[~np.isnan(st_sept4)]

may2_metrics = spaef_metrics(td_may2, st_may2, bins=bins, histrange=histrange)
may4_metrics = spaef_metrics(td_may4, st_may4, bins=bins, histrange=histrange)
june2_metrics = spaef_metrics(td_june2, st_june2, bins=bins, histrange=histrange)
june4_metrics = spaef_metrics(td_june4, st_june4, bins=bins, histrange=histrange)
july2_metrics = spaef_metrics(td_july2, st_july2, bins=bins, histrange=histrange)
july4_metrics = spaef_metrics(td_july4, st_july4, bins=bins, histrange=histrange)
august2_metrics = spaef_metrics(td_august2, st_august2, bins=bins, histrange=histrange)
august4_metrics = spaef_metrics(td_august4, st_august4, bins=bins, histrange=histrange)
september2_metrics = spaef_metrics(td_sept2, st_sept2, bins=bins, histrange=histrange)
september4_metrics = spaef_metrics(td_sept4, st_sept4, bins=bins, histrange=histrange)

# downsampling
sz = 1000
sampl2 = np.random.randint(low=0, high=len(td_june2), size=(sz,))
sampl4 = np.random.randint(low=0, high=len(td_june4), size=(sz,))
sampls2 = np.random.randint(low=0, high=len(td_sept2), size=(sz,))
sampls4 = np.random.randint(low=0, high=len(td_sept4), size=(sz,))

# min max
minim2 = round(np.nanmin(np.concatenate((st_june2, st_july2, st_august2, st_sept2, td_june2, td_july2, td_august2, td_sept2))) - 0.05, 1)
maksim2 = round(np.nanmax(np.concatenate((st_june2, st_july2, st_august2, st_sept2, td_june2, td_july2, td_august2, td_sept2))) + 0.05, 1)
minim4 = round(np.nanmin(np.concatenate((st_june4, st_july4, st_august4, st_sept4, td_june4, td_july4, td_august4, td_sept4))) - 0.05, 1)
maksim4 = round(np.nanmax(np.concatenate((st_june4, st_july4, st_august4, st_sept4, td_june4, td_july4, td_august4, td_sept4))) + 0.05, 1)

histalpha = 0.5
scatalpha = 0.1

# Plotting
fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(16,16));
#sns.set_style('whitegrid')
#plt.grid()
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
ax13 = axs[3][0]
ax14 = axs[3][1]
ax15 = axs[3][2]
ax16 = axs[3][3]

#fig, (ax1, ax2) = plt.subplots(1, 2 ,figsize=(16,10))
#gs = gridspec.GridSpec(1, 2)#, width_ratios=[1, 1])#, height_ratios=[1, 1]) 
#ax1 = plt.subplot(gs[0])
ax1.hist(td_june2, bins, alpha=histalpha, label='2D')
ax1.hist(st_june2, bins, alpha=histalpha, label='stand')
ax1.set_xlim(minim2, maksim2)
ax1.legend(loc='upper left')
ax1.set_title('MINERAL SOIL')
ax1.set_ylabel('JUNE')
ax1.text(0.15, 250000, f'HM: {june2_metrics[0]}')


#ax2.scatter(td_june2,st_june2, alpha=0.015)
sns.regplot(ax=ax2, x=td_june2[sampl2], y=st_june2[sampl2], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax2.set_ylim(minim2, maksim2)
ax2.set_xlim(minim2, maksim2)
ax2.set_ylabel('stand')
ax2.set_xlabel('2D')
ax2.set_title('MINERAL SOIL')
#z = np.polyfit(td_june2, st_june2, 1)
#p = np.poly1d(z)
#ax2.plot(td_june2,p(td_june2),"r--")
ax2.text(0.15, 0.45, f'CC: {june2_metrics[1]}')
ax2.text(0.15, 0.40, f'CV: {str(june2_metrics[2])}')
ax2.text(0.15, 0.35, f'SPAEF: {june2_metrics[3]}')

ax3.hist(td_june4, bins, alpha=histalpha, label='obs')
ax3.hist(st_june4, bins, alpha=histalpha, label='sim')
ax3.set_xlim(minim4, maksim4)
#ax3.legend(loc='upper left')
ax3.set_title('PEAT SOIL')
ax3.text(0.4, 250000, f'HM: {june4_metrics[0]}')

#ax4.scatter(td_june4[sampl4],st_june4[sampl4], vmin=0, vmax=1, alpha=0.01)
sns.regplot(ax=ax4, x=td_june4[sampl4], y=st_june4[sampl4], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax4.set_ylim(minim4, maksim4)
ax4.set_xlim(minim4, maksim4)
ax4.set_ylabel('stand')
ax4.set_xlabel('2D')
ax4.set_title('PEAT SOIL')
ax4.text(0.4, 0.9, f'CC: {june4_metrics[1]}')
ax4.text(0.4, 0.85, f'CV: {str(june4_metrics[2])}')
ax4.text(0.4, 0.80, f'SPAEF: {str(june4_metrics[3])}')

ax5.hist(td_july2, bins, alpha=histalpha, label='obs')
ax5.hist(st_july2, bins, alpha=histalpha, label='sim')
ax5.set_xlim(minim2, maksim2)
#ax5.legend(loc='upper left')
ax5.set_ylabel('JULY')
ax5.text(0.15, 150000, f'HM: {july2_metrics[0]}')

#ax6.scatter(td_july2,st_july2, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax6, x=td_july2[sampl2], y=st_july2[sampl2], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax6.set_ylim(minim2, maksim2)
ax6.set_xlim(minim2, maksim2)
ax6.set_ylabel('stand')
ax6.set_xlabel('2D')
ax6.text(0.15, 0.45, f'CC: {july2_metrics[1]}')
ax6.text(0.15, 0.40, f'CV: {str(july2_metrics[2])}')
ax6.text(0.15, 0.35, f'SPAEF: {july2_metrics[3]}')

ax7.hist(td_july4, bins, alpha=histalpha, label='obs')
ax7.hist(st_july4, bins, alpha=histalpha, label='sim')
ax7.set_xlim(minim4, maksim4)
#ax7.legend(loc='upper left') 
ax7.text(0.4, 150000, f'HM: {july4_metrics[0]}')

#ax8.scatter(td_july4,st_july4, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax8, x=td_july4[sampl4], y=st_july4[sampl4], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax8.set_ylim(minim4, maksim4)
ax8.set_xlim(minim4, maksim4)
ax8.set_ylabel('stand')
ax8.set_xlabel('2D')
ax8.text(0.4, 0.9, f'CC: {july4_metrics[1]}')
ax8.text(0.4, 0.85, f'CV: {str(july4_metrics[2])}')
ax8.text(0.4, 0.80, f'SPAEF: {july4_metrics[3]}')

ax9.hist(td_august2, bins, alpha=histalpha, label='obs')
ax9.hist(st_august2, bins, alpha=histalpha, label='sim')
ax9.set_xlim(minim2, maksim2)
#ax9.legend(loc='upper left') 
ax9.set_ylabel('AUGUST')
ax9.text(0.15, 150000, f'HM: {august2_metrics[0]}')

#ax10.scatter(td_august2,st_august2, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax10, x=td_august2[sampl2], y=st_august2[sampl2], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax10.set_ylim(minim2, maksim2)
ax10.set_xlim(minim2, maksim2)
ax10.set_ylabel('stand')
ax10.set_xlabel('2D')
ax10.text(0.15, 0.45, f'CC: {august2_metrics[1]}')
ax10.text(0.15, 0.40, f'CV: {str(august2_metrics[2])}')
ax10.text(0.15, 0.35, f'SPAEF: {str(august2_metrics[3])}')

ax11.hist(td_august4, bins, alpha=histalpha, label='obs')
ax11.hist(st_august4, bins, alpha=histalpha, label='sim')
ax11.set_xlim(minim4, maksim4)
#ax11.legend(loc='upper left') 
ax11.text(0.4, 150000, f'HM: {august4_metrics[0]}')

#ax12.scatter(td_august4,st_august4, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax12, x=td_august4[sampl4], y=st_august4[sampl4], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax12.set_ylim(minim4, maksim4)
ax12.set_xlim(minim4, maksim4)
ax12.set_ylabel('stand')
ax12.set_xlabel('2D')
ax12.text(0.4, 0.9, f'CC: {august4_metrics[1]}')
ax12.text(0.4, 0.85, f'CV: {str(august4_metrics[2])}')
ax12.text(0.4, 0.80, f'SPAEF: {str(august4_metrics[3])}')


ax13.hist(td_sept2, bins, alpha=histalpha, label='obs')
ax13.hist(st_sept2, bins, alpha=histalpha, label='sim')
ax13.set_xlim(minim2, maksim2)
#ax9.legend(loc='upper left') 
ax13.set_ylabel('SEPTEMBER')
ax13.text(0.15, 250000, f'HM: {september2_metrics[0]}')

#ax14.scatter(td_sept2,st_sept2, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax14, x=td_sept2[sampls2], y=st_sept2[sampls2], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax14.set_ylim(minim2, maksim2)
ax14.set_xlim(minim2, maksim2)
ax14.set_ylabel('stand')
ax14.set_xlabel('2D')
ax14.text(0.15, 0.45, f'CC: {september2_metrics[1]}')
ax14.text(0.15, 0.40, f'CV: {str(september2_metrics[2])}')
ax14.text(0.15, 0.35, f'SPAEF: {september2_metrics[3]}')

ax15.hist(td_sept4, bins, alpha=histalpha, label='obs')
ax15.hist(st_sept4, bins, alpha=histalpha, label='sim')
ax15.set_xlim(minim4, maksim4)
#ax11.legend(loc='upper left') 
ax15.text(0.4, 200000, f'HM: {september4_metrics[0]}')

#ax16.scatter(td_sept4,st_sept4, vmin=0, vmax=1,alpha=0.01)
sns.regplot(ax=ax16, x=td_sept4[sampls4], y=st_sept4[sampls4], scatter_kws={'s':50, 'alpha':0.1}, line_kws={"color": "red"})
ax16.set_ylim(minim4, maksim4)
ax16.set_xlim(minim4, maksim4)
ax16.set_ylabel('stand')
ax16.set_xlabel('2D')
ax16.text(0.4, 0.9, f'CC: {september4_metrics[1]}')
ax16.text(0.4, 0.85, f'CV: {str(september4_metrics[2])}')
ax16.text(0.4, 0.80, f'SPAEF: {september4_metrics[3]}')

plt.tight_layout()
ax1.grid(); ax2.grid(); ax3.grid(); ax4.grid(); ax5.grid(); ax6.grid(); ax7.grid(); ax8.grid()
ax9.grid(); ax10.grid(); ax11.grid(); ax12.grid(); ax13.grid(); ax14.grid(); ax15.grid(); ax16.grid()

plt.savefig(f'SCATTERHIST_{today}.png')
