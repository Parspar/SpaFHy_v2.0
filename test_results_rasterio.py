# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:41:31 2021

@author: janousu
"""

import rasterio
from rasterio.transform import from_origin
from rasterio.plot import show
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt
from iotools import read_results
from raster_utils import read_pkrasteri_for_extent

# READING CMASK FROM INPUT FILES
cmaskfp = r'C:\Users\janousu\OneDrive - Oulun yliopisto\SpaFHy_v1_Pallas_2D\testcase_input\parameters\cmask.dat'
cmask = rasterio.open(cmaskfp, 'r')
bbox = cmask.bounds
mask = cmask.read_masks(1)
mask = mask / 255
mask[np.where(mask==0)] = np.nan
shape=cmask.shape

# reading basic map
fp2 = r'C:\PALLAS_RAW_DATA\MML\Peruskartta\pkmosaic.tif'
pk, meta = read_pkrasteri_for_extent(fp2, bbox=bbox,showfig=False)

# reading the stand results
outputfile_stand = r'D:\SpaFHy_2D_2021\testcase_input_1d_new.nc'
results_stand = read_results(outputfile_stand)
arr = results_stand['bucket_moisture_root'][1000]


transform = from_origin(bbox[0], bbox[1], cmask.transform[0],  cmask.transform[0])

new_dataset = rasterio.open('test5.tif', 'w', driver='GTiff',
                            height = shape[0], width = shape[1],
                            count=1, dtype=str(results_stand['parameters_cmask'].dtype),
                            crs='epsg:3067',
                            transform=transform)

new_dataset.write(arr, 1)
new_dataset.close()

# reading top results
outputfile_top = r'D:\SpaFHy_2D_2021\testcase_input_top_new.nc'
results_top = read_results(outputfile_top)
arr = results_top['bucket_moisture_root'][1000]

transform = from_origin(bbox[0], bbox[1], cmask.transform[0],  cmask.transform[0])

new_dataset = rasterio.open('test6.tif', 'w', driver='GTiff',
                            height = shape[0], width = shape[1],
                            count=1, dtype=str(results_stand['parameters_cmask'].dtype),
                            crs='epsg:3067',
                            transform=transform)

new_dataset.write(arr, 1)
new_dataset.close()


# Open the files:
raster = rasterio.open(r'C:\SpaFHy_v1_Pallas_2D\figures\test5.tif')
raster2 = rasterio.open(r'C:\SpaFHy_v1_Pallas_2D\figures\test6.tif')


#%%

# show raster overlays
plt.close('all')

fig1, ax1 = plt.subplots(2,1, figsize=(10,10))
rasterio.plot.show(pk, transform=meta['transform'], ax=ax1[0])

rr = rasterio.plot.show(raster.read(1) * mask, transform=cmask.transform, ax=ax1[0], alpha=0.5)
# this creates colorbar
im = rr.get_images()[1]
#ax1[0].set_title('Site 1: Sat deficit doy 123-126')
fig1.colorbar(im, ax=ax1[0], shrink=1)

rasterio.plot.show(pk, transform=meta['transform'], ax=ax1[1]);
#show(twi * mask, transform=r.transform, ax=ax1, alpha=0.5, vmin=5., vmax=10.0)
rr = rasterio.plot.show(raster2.read(1)  *mask, transform=cmask.transform, ax=ax1[1], alpha=0.5)
im = rr.get_images()[1]
#ax1[1].set_title('vol. moisture (m3m-3)')
fig1.colorbar(im, ax=ax1[1], shrink=1)
plt.show()

fig1.savefig('Site1_moisture.png', dpi=300)


