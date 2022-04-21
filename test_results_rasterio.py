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
# reading the stand results
outputfile_stand = r'D:\SpaFHy_2D_2021\testcase_input_1d_new.nc'
results_stand = read_results(outputfile_stand)


XLLCORNER = 382500.0
YLLCORNER = 7541702.0
CELLSIZE = 16.0
SHAPE = results_stand['bucket_moisture_root'][1000].shape
arr = results_stand['bucket_moisture_root'][1000]

transform = from_origin(XLLCORNER, YLLCORNER, CELLSIZE, CELLSIZE)

new_dataset = rasterio.open('test5.tif', 'w', driver='GTiff',
                            height = SHAPE[0], width = SHAPE[1],
                            count=1, dtype=str(results_stand['parameters_cmask'].dtype),
                            crs='epsg:3067',
                            transform=transform)

new_dataset.write(arr, 1)
new_dataset.close()


# Open the file:
raster = rasterio.open(r'C:\SpaFHy_v1_Pallas_2D\figures\test4.tif')

# Plot band 1
show((raster, 1))

point1 = Point(7539500, 384000)
point1 = Point(7539500, 384000)


fig, ax = plt.subplots(figsize=(12, 8))

image_hidden = ax.imshow(raster.read(1),
                         cmap='coolwarm_r',
                         vmin=0,
                         vmax=1)

ax.scatter(point1.y, point1.x)


image = rasterio.plot.show(raster,
                        transform=raster.transform,
                        ax=ax,
                        cmap='coolwarm_r',
                        vmin=0,
                        vmax=1)

fig.colorbar(image_hidden, ax=ax)
plt.show()
