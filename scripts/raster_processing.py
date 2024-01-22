# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:50:15 2022

@author: janousu
"""

import matplotlib.pyplot as plt
import rasterio
from rasterio.merge import merge
from rasterio.plot import show
import os
from pathlib import Path
import numpy as np
from PIL import Image
from raster_utils import convert_peruskarttarasteri, show_raster, window_from_extent, read_pkrasteri_for_extent


# Data dir
data_dir = r"C:\PALLAS_RAW_DATA\MML\Peruskartta"

fps = list(Path(data_dir).glob('*.png'))

# convert png to geotiff and return list of filenames
outfiles = []
for fp in fps:
    print(fp)
    f = convert_peruskarttarasteri(str(fp), epsg_code='3067')
    outfiles.append(f)

#%%


s= rasterio.open(outfiles[0], 'r')
fig1, ax1 = plt.subplots(1,1, figsize=(10,10))
show(s, transform=s.meta['transform'], ax=ax1)


#%%

pk_mosaic, out_trans = merge(outfiles)
out_meta = s.meta.copy()
out_meta.update({"driver": "GTiff",
                 "height": pk_mosaic.shape[1],
                 "width": pk_mosaic.shape[2],
                "transform": out_trans,
                })
#show(pk_mosaic)

with rasterio.open(r'C:\PALLAS_RAW_DATA\MML\Peruskartta\pkmosaic.tif', "w", **out_meta) as dest:
    dest.write(pk_mosaic)
print('writing ok!')