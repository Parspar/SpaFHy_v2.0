# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 19:39:45 2022
@author: 03081268
"""

import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
import os
import numpy as np
from PIL import Image

def convert_peruskarttarasteri(fp, epsg_code='3067'):
    """
    reads grayscale (0-255) png peruskarttarasteri and converts it to 3-band
    georeferenced tif.
    epsg_code 3067 # ETRS-TM35FIN
    """
    # read attributes
    crs = rasterio.crs.CRS({"init": "epsg:" + epsg_code})
    r = rasterio.open(fp, 'r+')
    r.crs = crs
    info = r.meta
    print(info)
    r.close()

    # read P-mode image, convert to rgb
    img = Image.open(fp)
    rgb = Image.new("RGB", img.size)
    rgb.paste(img)
    x = np.asarray(rgb, dtype=int)
    img.close()

    # create new rasteriodataset in tif-format
    outfile = fp.split('.')[0] + '.tif'
    rout = rasterio.open(outfile, 'w',
                         driver='GTiff',
                         height=x.shape[0],
                         width=x.shape[1],
                         count=3,
                         dtype='uint8',
                         crs=info['crs'],
                         transform=info['transform']
                         )

    print('writing...')
    for k in range(3):
        rout.write(x[:,:,k], k+1)

    rout.close()

    return outfile

def read_pkrasteri_for_extent(fp, bbox, showfig=False):
    """
    Read geotif 3-band raster for extent
    Args:
        fp - filepath
        bbox - rasterio bounding box or list [xmin, ymin, xmax, ymax]
        showfig - boolean
    Return
        pk - rasterio data
        meta - metadata
    """

    r = rasterio.open(fp, 'r')
    meta = r.meta.copy()
    aff = r.transform

    # window
    window = window_from_extent(bbox, aff)
    pk = r.read(window=window)

    # Update dataset metadata
    meta.update(height = window[0][1] - window[0][0],
                width = window[1][1] - window[1][0],
                transform = r.window_transform(window))

    if showfig:
        fig1, ax1 = plt.subplots(1,)
        show(pk, transform=meta['transform'], ax=ax1)

    return pk, meta


def window_from_extent(bbox, aff):
    """
    create window from extent and affine
    """
    xmin, xmax, ymin, ymax = bbox[0], bbox[2], bbox[1], bbox[3]
    col_start, row_start = ~aff * (xmin, ymax)
    col_stop, row_stop = ~aff * (xmax, ymin)

    return ((int(row_start), int(row_stop)), (int(col_start), int(col_stop)))

def show_raster(r):

    show(r.read(), transform=r.transform)