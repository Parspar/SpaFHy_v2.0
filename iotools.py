# -*- coding: utf-8 -*-
"""
Input/output utilities for SpaFHy.

Provides functions for:
    - Reading GIS raster inputs (ESRI ASCII grid) for canopy, soil, soilprofile2D and TOPMODEL
    - Preprocessing spatial parameter dictionaries for model initialization
    - Reading and formatting meteorological forcing data
    - Initializing and writing NetCDF4 output files
    - Stitching multiple sub-catchment NetCDF results into a single file
    - Utility ASCII grid read/write functions

@authors: slauniai, khaahti, jpnousu
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from soilprofile2D import gwl_Wsto, gwl_Wsto_vectorized
from koordinaattimuunnos import koordTG
from topmodel import twi as twicalc
import re
import importlib
import sys

eps = np.finfo(float).eps  # machine epsilon
workdir = os.getcwd()

# Global variables to store the imported parameters
pgen, pcpy, pbu, pspd, ptop, aux = None, None, None, None, None, None

def initialize_parameters(catchment, folder):
    """
    Dynamically import parameters module based on the catchment and initialize global variables.
    """
    global pgen, pcpy, pbu, pspd, ptop, aux

    # Dynamically import the correct parameters module
    parameters_module = importlib.import_module(f'parameters_{catchment}')

    # Initialize the parameters
    pgen, pcpy, pbu, pspd = parameters_module.parameters(folder=folder)
    ptop = parameters_module.ptopmodel()
    aux = parameters_module.auxiliary_grids()


def read_bu_gisdata(fpath, spatial_pbu, mask=None, plotgrids=False):
    """
    Reads bucket (soil) GIS raster grids and returns numpy 2D arrays.

    Args:
        fpath        (str):  Path to the GIS data folder.
        spatial_pbu  (dict): Flags indicating which fields to read spatially
                             (True = read from raster, False = use scalar).
        mask         :       Unused; reserved for future masking support.
        plotgrids   (bool):  If True, plots the loaded grids.

    Returns:
        gis (dict): Loaded raster arrays and grid metadata:
            'org_id'     [-]  organic layer soil class index
            'root_id'    [-]  root zone soil class index
            'dxy'        [m]  cell size
            'xllcorner'  [m]  lower-left x coordinate
            'yllcorner'  [m]  lower-left y coordinate
    """

    fpath = os.path.join(workdir, fpath)

    gis = {}

    # overlandflow
    if 'flowacc' in spatial_pbu:
        if spatial_pbu['flowacc']:
            flowacc, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['flowacc']))
            gis['flowacc'] = flowacc

    if 'fdir' in spatial_pbu:
        if spatial_pbu['fdir']:
            fdir, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['fdir']))
            gis['fdir'] = fdir

    if 'streams' in spatial_pbu:
        if spatial_pbu['streams']:
            streams, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['streams']))
            gis['streams'] = streams

    if 'lakes' in spatial_pbu:
        if spatial_pbu['lakes']:
            lakes, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['lakes']))
            gis['lakes'] = lakes

    # soil classification
    if 'org_id' in spatial_pbu:
        if spatial_pbu['org_id']:
            org_id, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['org_id']))
            gis['org_id'] = org_id

    if 'root_id' in spatial_pbu:
        if spatial_pbu['root_id']:
            root_id, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['root_id']))
            gis['root_id'] = root_id

    xllcorner = int(re.findall(r'\d+', info[2])[0])
    yllcorner = int(re.findall(r'\d+', info[3])[0])

    if plotgrids is True:
        plt.figure()
        plt.subplot(311); plt.imshow(root_id); plt.colorbar(); plt.title('root_id')
        plt.subplot(312); plt.imshow(org_id); plt.colorbar(); plt.title('org_id')

    gis.update({'dxy': cellsize})
    gis.update({'xllcorner': xllcorner,
                'yllcorner': yllcorner})
    
    return gis

def read_cpy_gisdata(fpath, spatial_pcpy, mask=None, plotgrids=False):
    """
    Reads canopy GIS raster grids and returns numpy 2D arrays.

    Args:
        fpath        (str):  Path to the GIS data folder.
        spatial_pcpy (dict): Flags indicating which fields to read spatially.
        mask         :       Unused; reserved for future masking support.
        plotgrids   (bool):  If True, plots the loaded grids.

    Returns:
        gis (dict): Loaded raster arrays:
            'canopy_height'   [m]       mean stand height
            'canopy_fraction' [-]       canopy closure fraction
            'LAI_conif'       [m2 m-2]  conifer annual max LAI
            'LAI_decid'       [m2 m-2]  deciduous annual max LAI
            'LAI_shrub'       [m2 m-2]  shrub LAI (derived as 0.1 * LAI_conif if not provided)
            'LAI_grass'       [m2 m-2]  grass LAI (derived as 0.5 * LAI_decid if not provided)
    """

    fpath = os.path.join(workdir, fpath)

    gis = {}
    
    # tree height [m]
    if 'canopy_height' in spatial_pcpy:
        if spatial_pcpy['canopy_height']:
            canopy_height, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['canopy_height']))
            gis['canopy_height'] = canopy_height

    # canopy closure [-]
    if 'canopy_fraction' in spatial_pcpy:
        if spatial_pcpy['canopy_fraction']:
            canopy_fraction, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['canopy_fraction']))
            gis['canopy_fraction'] = canopy_fraction

    if 'LAI_conif' in spatial_pcpy:
        if spatial_pcpy['LAI_conif']:
            LAI_conif, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_conif']))
            LAI_shrub = 0.1 * LAI_conif
            gis['LAI_conif'] = LAI_conif
            gis['LAI_shrub'] = LAI_shrub

    if 'LAI_decid' in spatial_pcpy:
        if spatial_pcpy['LAI_decid']:
            LAI_decid, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_decid']))
            LAI_grass = 0.5 * LAI_decid
            gis['LAI_decid'] = LAI_decid
            gis['LAI_grass'] = LAI_grass

    if 'LAI_shrub' in spatial_pcpy:
        if spatial_pcpy['LAI_shrub']:
            LAI_shrub, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_shrub']))
            gis['LAI_shrub'] = LAI_shrub

    if 'LAI_grass' in spatial_pcpy:
        if spatial_pcpy['LAI_grass']:
            LAI_grass, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_grass']))
            gis['LAI_grass'] = LAI_grass

    if plotgrids:

        plt.figure()
        #plt.subplot(221); plt.imshow(LAI_pine+LAI_spruce); plt.colorbar();
        plt.title('LAI conif (m2/m2)')
        plt.subplot(222); plt.imshow(LAI_decid); plt.colorbar();
        plt.title('LAI decid (m2/m2)')
        plt.subplot(223); plt.imshow(canopy_height); plt.colorbar(); plt.title('hc (m)')
        plt.subplot(224); plt.imshow(canopy_fraction); plt.colorbar(); plt.title('cf (-)')

    return gis

def read_ds_gisdata(fpath, spatial_pspd, mask=None, plotgrids=False):
    """
    Reads deep soil GIS raster grids and returns numpy 2D arrays.

    Args:
        fpath        (str):  Path to the GIS data folder.
        spatial_pspd (dict): Flags indicating which fields to read spatially.
        mask         :       Unused; reserved for future masking support.
        plotgrids   (bool):  If True, plots the loaded grids.

    Returns:
        gis (dict): Loaded raster arrays and grid metadata:
            'deep_id'      [-]  deep soil type index
            'deep_z'       [m]  depth to impermeable bottom
            'elevation'    [m]  surface DEM
            'streams'      [m]  stream/ditch water depth (negative where present)
            'stream_depth' [m]  spatially varying stream depth if available
            'lakes'        [m]  lake depth (negative where present)
            'dxy'          [m]  cell size
            'xllcorner'    [m]  lower-left x coordinate
            'yllcorner'    [m]  lower-left y coordinate
    """
    fpath = os.path.join(workdir, fpath)

    gis = {}

    # deep soil layer
    if 'deep_id' in spatial_pspd:
        if spatial_pspd['deep_id']:
            deep_id, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['deep_id']))
            gis['deep_id'] = deep_id

    # soil depth
    if 'deep_z' in spatial_pspd:
        if spatial_pspd['deep_z']:
            deep_z, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['deep_z']))
            gis['deep_z'] = deep_z

    # dem
    if 'elevation' in spatial_pspd:
        if spatial_pspd['elevation']:
            elevation, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['elevation']))
            gis['elevation'] = elevation

    # streams
    if 'streams' in spatial_pspd:
        if spatial_pspd['streams']:
            streams, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['streams']))
            if isinstance(pspd['stream_depth'], float):
                streams[(np.isfinite(streams)) & (streams != 0.0)] = pspd['stream_depth']
                streams[~np.isfinite(streams)] = 0.0
            else:
                streams[(np.isfinite(streams)) & (streams != 0.0)] = -1.0
                streams[streams != -1.0] = 0.0

    if 'stream_depth' in spatial_pspd:
        if spatial_pspd['stream_depth']:
            stream_depth, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['stream_depth']))
            streams = stream_depth.copy()
            streams[~np.isfinite(streams)] = 0.0
            #stream_depth[~np.isfinite(stream_depth)] = 0.0
            gis['stream_depth'] = stream_depth
        else:
            stream_depth = pspd['stream_depth']

    if 'stream_length' in spatial_pspd:
        if spatial_pspd['stream_length'] == True:
            stream_length, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['stream_length']))
            #stream_length[~np.isfinite(stream_length)] = 0.0
            gis['stream_length'] = stream_length
        else:
            stream_length = pspd['stream_length']

    if 'stream_distance' in spatial_pspd:
        if spatial_pspd['stream_distance'] == True:
            stream_distance, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['stream_distance']))
            #stream_distance[~np.isfinite(stream_distance)] = 0.0
            gis['stream_distance'] = stream_distance
        else:
            stream_distance = pspd['stream_distance']

    if 'stream_width' in spatial_pspd:
        if spatial_pspd['stream_width'] == True:
            stream_width, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['stream_width']))
            #stream_width[~np.isfinite(stream_width)] = 0.0
            gis['stream_width'] = stream_width
        else:
            stream_width = pspd['stream_width']

    if not (('streams' in spatial_pspd and spatial_pspd['streams']) or 
        ('stream_depth' in spatial_pspd and spatial_pspd['stream_depth'])):
        print('*** No stream file ***')
        streams = np.full_like(deep_id, 0.0)
        stream_depth = np.full_like(deep_id, 0.0)

    gis['streams'] = streams
    gis['stream_depth'] = stream_depth

    # lakes if available
    if 'lakes' in spatial_pspd:
        if spatial_pspd['lakes']:
            lakes, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['lakes']))
            lakes[np.isfinite(lakes)] = -1.0
    else:
        print('*** No lakes file ***')
        lakes = np.full_like(deep_id, 0.0)
    gis['lakes'] = lakes
    xllcorner = int(re.findall(r'\d+', info[2])[0])
    yllcorner = int(re.findall(r'\d+', info[3])[0])

    if plotgrids:
        plt.figure()
        plt.subplot(311); plt.imshow(deep_id); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(elevation); plt.colorbar(); plt.title('elevation')
        plt.subplot(313); plt.imshow(elevation); plt.colorbar();

    gis.update({'dxy': cellsize})
    gis.update({'xllcorner': xllcorner,
                'yllcorner': yllcorner})
    return gis

def read_top_gisdata(fpath, spatial_ptop, mask=None, plotgrids=False):
    """
    Reads TOPMODEL GIS raster grids and returns numpy 2D arrays.

    Args:
        fpath        (str):  Path to the GIS data folder.
        spatial_ptop (dict): Flags indicating which fields to read spatially.
        mask         :       Unused; reserved for future masking support.
        plotgrids   (bool):  If True, plots the loaded grids.

    Returns:
        gis (dict): Loaded raster arrays and grid metadata:
            'flowacc' [m]   flow accumulation per unit contour length
            'slope'   [deg] local slope
            'twi'     [-]   topographic wetness index
            'dxy'     [m]   cell size
    """
    fpath = os.path.join(workdir, fpath)

    gis = {}

    # flow accumulation
    if 'flow_accumulation' in spatial_ptop:
        if spatial_ptop['flow_accumulation']:
            flowacc, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['flow_accumulation']))
            gis['flowacc'] = flowacc

    # slope
    if 'slope' in spatial_ptop:
        if spatial_ptop['slope']:
            slope, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['slope']))
            gis['slope'] = slope

    if 'twi' in spatial_ptop:
        if spatial_ptop['twi']:
            twi, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['twi']))
            gis['twi'] = twi

    if plotgrids:
        plt.figure()
        plt.subplot(311); plt.imshow(slope); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(twi); plt.colorbar(); plt.title('elevation')
        plt.subplot(313); plt.imshow(flowacc); plt.colorbar();

    gis.update({'dxy': cellsize})

    return gis

def read_aux_gisdata(fpath, spatial_aux, mask=None):
    """
    Reads auxiliary GIS raster grids (catchment mask, streams, lakes).

    Args:
        fpath       (str):  Path to the GIS data folder.
        spatial_aux (dict): Keys indicating which auxiliary grids to load
                            ('cmask', 'streams', 'lakes').
        mask        :       Unused; reserved for future masking support.

    Returns:
        gis (dict): Loaded raster arrays and grid metadata:
            'cmask'   [-]  catchment mask (1 inside, NaN outside)
            'streams' [-]  stream locations
            'lakes'   [-]  lake locations
            'dxy'     [m]  cell size
    """
    fpath = os.path.join(workdir, fpath)

    gis = {}
    
    if 'cmask' in spatial_aux:
        cmask, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, aux['cmask']))
        gis['cmask'] = cmask

    if 'streams' in spatial_aux:
        streams, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, aux['streams']))
        gis['streams'] = streams

    if 'lakes' in spatial_aux:
        lakes, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, aux['lakes']))
        gis['lakes'] = lakes
    else:
        lakes = np.full_like(cmask, 0.0)
        gis['lakes'] = lakes

    gis.update({'dxy': cellsize})

    return gis


def read_spatial_forcing(fpath):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
    Returns:
        gis - dict of gis-data rasters
            cmask
            lat
            lon
            forcing_id

    """
    fpath = os.path.join(workdir, fpath)

    # latitude and longitude
    Ncoord, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'Ncoord.dat'))
    Ecoord, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'Ecoord.dat'))
    lat, lon = koordTG(Ncoord, Ecoord)

    forcing_id, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'forcing_id.dat'))

    # catchment mask cmask[i,j] == 1, np.NaN outside
    if os.path.isfile(os.path.join(fpath, 'cmask.dat')):
        cmask, _, _, _, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.dat'))
    else:
        cmask = np.ones(np.shape(Ncoord))

    # dict of all rasters
    gis = {'cmask': cmask, 'lat': lat, 'lon': lon, 'forcing_id': forcing_id}

    for key in gis.keys():
        gis[key] *= cmask

    return gis


def preprocess_budata(pbu, spatial_pbu, orgp, rootp, gisdata, spatial=True):
    """
    Builds the input parameter dictionary for BucketGrid initialization.

    Merges scalar parameter defaults with spatial raster data, then fills
    per-soil-type hydraulic properties from orgp and rootp lookup tables.

    Args:
        pbu        (dict): Bucket parameter defaults (scalar values).
        spatial_pbu (dict): Flags: True = read field from gisdata, False = use scalar.
        orgp       (dict): Organic layer hydraulic properties keyed by soil type.
        rootp      (dict): Root zone hydraulic properties keyed by soil type.
        gisdata    (dict): GIS raster arrays including 'org_id', 'root_id', 'dxy'.
        spatial    (bool): If True, fills per-soil-type properties from orgp/rootp.

    Returns:
        data (dict): Fully populated parameter arrays (same grid shape as gisdata)
                     ready for BucketGrid.__init__().
    """
    # create dict for initializing soil profile.
    # copy pbu into sdata and make each value np.array(np.shape(cmask))
    data = pbu.copy()
    spatial_data = spatial_pbu.copy()
    gridshape = np.ones(shape=gisdata['org_id'].shape) # replace org_id with something generic?

    #  Create a mask where gisdata['org_id'] is finite (not NaN)
    mask = np.isfinite(gisdata['org_id'])

    for key in data:
        if spatial_data[key]:
            data[key] = gisdata[key]
        else:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, np.nan)
            data[key][mask] = uni_value  # Assign uni_value only where the mask is True
            
    root_ids = []
    for key, value in rootp.items():
        if ~np.isnan(value['root_id']):
            root_ids.append(value['root_id'])

    org_ids = []
    for key, value in orgp.items():
        if ~np.isnan(value['org_id']):
            org_ids.append(value['org_id'])
    
    # --- Root soil IDs check ---
    used_root_ids = set(np.unique(data['root_id'][np.isfinite(gisdata['org_id'])]).tolist())
    defined_root_ids = set(root_ids)

    if not defined_root_ids >= used_root_ids:
        raise ValueError(
            f"Root soil IDs missing in parameters. Defined: {defined_root_ids}, Used: {used_root_ids}"
        )

    # --- Organic soil IDs check ---
    used_org_ids = set(np.unique(data['org_id'][np.isfinite(gisdata['org_id'])]).tolist())
    defined_org_ids = set(org_ids)

    if not defined_org_ids >= used_org_ids:
        raise ValueError(
            f"Organic soil IDs missing in parameters. Defined: {defined_org_ids}, Used: {used_org_ids}"
        )

    if spatial:
        for key, value in orgp.items():
            t = value['org_id']
            yx = np.where(data['org_id'] == t)
            data['org_depth'][yx] = value['org_depth']
            data['org_poros'][yx] = value['org_poros']
            data['org_fc'][yx] = value['org_fc']
            data['org_rw'][yx] = value['org_rw']
            data['org_ksat'][yx] = value['org_ksat']
            data['org_beta'][yx] = value['org_beta']

    if spatial:
        for key, value in rootp.items():
            t = value['root_id']
            yx = np.where(data['root_id'] == t)
            data['root_fc'][yx] = value['root_fc']
            data['root_ksat'][yx] = value['root_ksat']
            data['root_poros'][yx] = value['root_poros']
            data['root_wp'][yx] = value['root_wp']
            data['root_beta'][yx] = value['root_beta']
            data['root_alpha'][yx] = value['root_alpha']
            data['root_n'][yx] = value['root_n']            
            data['root_wr'][yx] = value['root_wr']

    data['dxy'] = gisdata['dxy']
    for key in ('flowacc', 'fdir', 'streams', 'lakes'):
        if key in gisdata:
            data[key] = gisdata[key]

    return data


def preprocess_dsdata(pspd, spatial_pspd, deepp, gisdata, spatial=True):
    """
    Builds the input parameter dictionary for SoilGrid_2Dflow initialization
    using soiltype-wise (uniform per soil class) interpolation functions.

    Merges scalar defaults with spatial rasters, assigns per-soil-type
    hydraulic lookup functions from gwl_Wsto.

    Args:
        pspd        (dict): Deep soil parameter defaults.
        spatial_pspd (dict): Flags: True = read from gisdata, False = use scalar.
        deepp       (dict): Per-soil-type hydraulic properties and pF curves.
        gisdata     (dict): GIS raster arrays including 'deep_id', 'elevation',
                            'streams', 'lakes', 'dxy'.
        spatial     (bool): If True, fills spatial arrays from gisdata.

    Returns:
        data (dict): Fully populated parameter arrays ready for
                     SoilGrid_2Dflow.__init__(), including scipy interp1d
                     lookup functions keyed by soil type.
    """
    data = pspd.copy()
    spatial_data = spatial_pspd.copy()
    gridshape = np.ones(shape=gisdata['deep_id'].shape)

    for key in data:
        if spatial_data[key]:
            data[key] = gisdata[key]
        else:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, uni_value)

    if not spatial:
        data['deep_id'] = pspd['deep_id']
    else:
        data['deep_id'] = gisdata['deep_id']
        data['elevation'] = gisdata['elevation']
        data['streams'] = gisdata['streams']
        data['lakes'] = gisdata['lakes']

    deep_ids = []
    for key, value in deepp.items():
        deep_ids.append(value['deep_id'])

    used_ids = set(np.unique(data['deep_id'][np.isfinite(gisdata['deep_id'])]).tolist())
    defined_ids = set(deep_ids)
    if not defined_ids >= used_ids:
        raise ValueError(
            f"Deep soil IDs missing. Defined: {defined_ids}, Used: {used_ids}"
        )

    data.update({'soiltype': np.empty(np.shape(gisdata['deep_id']),dtype=object)})
    data.update({'deep_z': np.empty(np.shape(gisdata['deep_id']),dtype=object)})

    for key, value in deepp.items():
        c = value['deep_id']
        ix = np.where(data['deep_id'] == c)
        data['soiltype'][ix] = key
        data['deep_z'][ix] = value['deep_z']
        # interpolation function between wsto and gwl
        value.update(gwl_Wsto(value['deep_z'], value['pF'], value['deep_ksat']))
        # interpolation function between root_wsto and gwl
        value.update(gwl_Wsto(value['deep_z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True))

    # stream depth corresponding to assigned parameter
    #data['streams'] = np.where((data['streams'] < -eps) | (data['lakes'] < -eps), pspd['stream_depth'], 0)
    #data['streams'] = np.where(data['streams'] < -eps, pspd['stream_depth'], 0)
    data['lakes'] = np.where(data['lakes'] < -eps, pspd['lake_depth'], 0)
    #data['streams'] = np.where(data['lakes'] < -eps, pspd['stream_depth'], 0)
    
    data['wtso_to_gwl'] = {soiltype: deepp[soiltype]['to_gwl'] for soiltype in deepp.keys()}
    data['gwl_to_wsto'] = {soiltype: deepp[soiltype]['to_wsto'] for soiltype in deepp.keys()}
    data['gwl_to_C'] = {soiltype: deepp[soiltype]['to_C'] for soiltype in deepp.keys()}
    data['gwl_to_Tr'] = {soiltype: deepp[soiltype]['to_Tr'] for soiltype in deepp.keys()}
    data['gwl_to_rootmoist'] = {soiltype: deepp[soiltype]['to_rootmoist'] for soiltype in deepp.keys()}
    data['dxy'] = gisdata['dxy']

    return data


def preprocess_dsdata_vec(pspd, spatial_pspd, deepp, gisdata, spatial=True):
    """
    Builds the input parameter dictionary for SoilGrid_2Dflow initialization
    using cell-wise (spatially varying) interpolation functions.

    When spatial_pspd['deep_z'] is True, soil depth varies per grid cell and
    gwl_Wsto_vectorized is used to build per-cell lookup functions. Otherwise
    falls back to soiltype-wise gwl_Wsto (same as preprocess_dsdata).

    Args:
        pspd         (dict): Deep soil parameter defaults.
        spatial_pspd (dict): Flags: True = read from gisdata, False = use scalar.
        deepp        (dict): Per-soil-type hydraulic properties and pF curves.
        gisdata      (dict): GIS raster arrays including 'deep_id', 'deep_z',
                             'elevation', 'streams', 'lakes', 'dxy'.
        spatial      (bool): If True, fills spatial arrays from gisdata.

    Returns:
        data (dict): Fully populated parameter arrays ready for
                     SoilGrid_2Dflow.__init__(), including scipy interp1d
                     lookup functions (per cell or per soil type depending on
                     spatial_pspd['deep_z']).
    """
    data = pspd.copy()
    spatial_data = spatial_pspd.copy()

    gridshape = np.ones(shape=gisdata['deep_id'].shape)

    for key in data:
        if spatial_data[key]:
            data[key] = gisdata[key]
        else:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, uni_value)

    if not spatial:
        data['deep_id'] = pspd['deep_id']
    else:
        data['deep_id'] = gisdata['deep_id']
        data['elevation'] = gisdata['elevation']
        data['streams'] = gisdata['streams']
        data['lakes'] = gisdata['lakes']

    deep_ids = []
    for key, value in deepp.items():
        deep_ids.append(value['deep_id'])

    used_ids = set(np.unique(data['deep_id'][np.isfinite(gisdata['deep_id'])]).tolist())
    defined_ids = set(deep_ids)
    if not defined_ids >= used_ids:
        raise ValueError(
            f"Deep soil IDs missing. Defined: {defined_ids}, Used: {used_ids}"
        )

    data.update({'soiltype': np.empty(np.shape(gisdata['deep_id']),dtype=object)})

    if not spatial_data['deep_z']:
        for key, value in deepp.items():
            c = value['deep_id']
            ix = np.where(data['deep_id'] == c)
            data['soiltype'][ix] = key
            data['deep_z'][ix] = value['deep_z']
            # interpolation function between wsto and gwl
            value.update(gwl_Wsto(value['deep_z'], value['pF'], -0.01, value['deep_ksat']))
            # interpolation function between root_wsto and gwl
            value.update(gwl_Wsto(value['deep_z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True))
        
        data['wtso_to_gwl'] = {soiltype: deepp[soiltype]['to_gwl'] for soiltype in deepp.keys()}
        data['gwl_to_wsto'] = {soiltype: deepp[soiltype]['to_wsto'] for soiltype in deepp.keys()}
        data['gwl_to_C'] = {soiltype: deepp[soiltype]['to_C'] for soiltype in deepp.keys()}
        data['gwl_to_Tr'] = {soiltype: deepp[soiltype]['to_Tr'] for soiltype in deepp.keys()}
        data['gwl_to_rootmoist'] = {soiltype: deepp[soiltype]['to_rootmoist'] for soiltype in deepp.keys()}

    elif spatial_data['deep_z']:
        # we have data['deep_id'] and data['z']
        max_nlyrs = 0    
        for key, value in deepp.items():
            nlyrs = len(value['deep_z'])
            max_nlyrs = max(max_nlyrs, nlyrs)
        # flattening
        deep_id_f = data['deep_id'].flatten()
        deep_z = data['deep_z']
        deep_z[deep_z < 5] = 5. # NOTE MINIMUM IS 5M DEPTH!
        deep_z_f = deep_z.flatten()
        # creating the arrays
        deep_zs = np.full((len(deep_id_f), max_nlyrs), np.nan)
        deep_ksats = np.full((len(deep_id_f), max_nlyrs), np.nan)

        # Initialize deep_pFs with dictionaries containing NaN values
        default_pF = {
            'ThetaS': [np.nan] * max_nlyrs,
            'ThetaR': [np.nan] * max_nlyrs,
            'alpha': [np.nan] * max_nlyrs,
            'n': [np.nan] * max_nlyrs
        }

        deep_pFs = np.array([default_pF.copy() for _ in range(len(deep_id_f))], dtype=object)
        
        # temporary arrays to store the interpolation functions
        temp_to_gwl = np.empty(len(deep_id_f), dtype=object)
        temp_to_wsto = np.empty(len(deep_id_f), dtype=object)
        temp_to_Tr = np.empty(len(deep_id_f), dtype=object)
        temp_to_C = np.empty(len(deep_id_f), dtype=object)
        temp_to_rootmoist = np.empty(len(deep_id_f), dtype=object)

        # making parametes into 2D arrays (total length and layered information)
        for key, value in deepp.items():
            mask = deep_id_f == value['deep_id']
            if np.any(mask):  # Only proceed if at least one match
                nlyrs = len(value['deep_z'])
                deep_zs[mask, :nlyrs] = value['deep_z']
                a = np.abs(deep_z_f[mask])*-1
                b = deep_zs[mask, nlyrs - 1]
                # Replace last layer. Cannot be smaller than smallest assigned 'z' in parameters
                deep_zs[mask, nlyrs - 1] = np.minimum(np.abs(deep_z_f[mask])*-1, deep_zs[mask, nlyrs - 1])
                # Cannot be smaller than -30.
                deep_ksats[mask, :nlyrs] = value['deep_ksat']
                deep_pFs[mask] = value['pF']

        mask = np.isfinite(deep_zs[:,0])
        ifs_v = gwl_Wsto_vectorized(deep_zs[mask], deep_pFs[mask], -0.2, deep_ksats[mask])
        ifs_r = gwl_Wsto_vectorized(value['deep_z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True)
        
        temp_to_gwl[mask] = ifs_v['to_gwl']
        temp_to_wsto[mask] = ifs_v['to_wsto']
        temp_to_C[mask] = ifs_v['to_C']
        temp_to_Tr[mask] = ifs_v['to_Tr']
        temp_to_rootmoist[mask] = ifs_r['to_rootmoist']

        data['wtso_to_gwl'] = temp_to_gwl.reshape(gridshape.shape)
        data['gwl_to_wsto'] = temp_to_wsto.reshape(gridshape.shape)
        data['gwl_to_C'] = temp_to_C.reshape(gridshape.shape)
        data['gwl_to_Tr'] = temp_to_Tr.reshape(gridshape.shape)
        data['gwl_to_rootmoist'] = temp_to_rootmoist.reshape(gridshape.shape)
        data['deep_z'] = deep_z_f.reshape(data['deep_z'].shape)

    data['lakes'] = np.where(data['lakes'] < -eps, pspd['lake_depth'], 0)

    data['dxy'] = gisdata['dxy']

    return data

def preprocess_cpydata(pcpy, spatial_pcpy, gisdata, spatial=True):
    """
    Builds the input parameter dictionary for CanopyGrid initialization.

    Replaces scalar values in pcpy['state'] with spatial raster arrays from
    gisdata wherever spatial_pcpy flags are True.

    Args:
        pcpy        (dict): Canopy parameter dict; 'state' sub-dict holds
                            initial state fields (LAI, canopy_height, etc.).
        spatial_pcpy (dict): Flags: True = read field from gisdata, False = use scalar.
        gisdata     (dict): GIS raster arrays (LAI_conif, LAI_decid, LAI_shrub,
                            LAI_grass, canopy_height, canopy_fraction).
        spatial     (bool): Unused; retained for interface consistency.

    Returns:
        pcpy (dict): Updated canopy parameter dict with spatial arrays in 'state'.
    """
    data = pcpy['state'].copy()
    spatial_data = spatial_pcpy.copy()
    gridshape = np.ones(shape=gisdata['LAI_conif'].shape)

    for key in data:
        if spatial_data[key]:
            data[key] = gisdata[key]
        else:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, uni_value)

    pcpy['state'] = data

    return pcpy

def preprocess_topdata(ptop, spatial_ptop, gisdata, spatial=True):
    """
    Builds the input parameter dictionary for Topmodel_Homogenous initialization.

    Inserts spatial raster arrays from gisdata into the ptop dict.

    Args:
        ptop        (dict): TOPMODEL parameter dict (from parameters_*.ptopmodel()).
        spatial_ptop (dict): Flags indicating which fields are spatial (unused here;
                             all GIS fields are always inserted).
        gisdata     (dict): GIS raster arrays: 'slope', 'flowacc', 'twi', 'dxy',
                            and optionally 'lat', 'lon'.
        spatial     (bool): Unused; retained for interface consistency.

    Returns:
        ptop (dict): Updated TOPMODEL parameter dict with spatial arrays inserted.
    """
    # inputs for topmodel initialization: update ptop using spatial data
    
    ptop['slope'] = gisdata['slope']
    ptop['flowacc'] = gisdata['flowacc']
    ptop['twi'] = gisdata['twi']
    ptop['dxy'] = gisdata['dxy']
    if {'lat','lon'}.issubset(gisdata.keys()):
        ptop['loc']['lat'] = gisdata['lat']
        ptop['loc']['lon'] = gisdata['lon']
    
    return ptop


def read_HESS2019_weather(start_date, end_date, sourcefile, CO2=380.0, U=2.0, ID=0):
    """
    reads FMI interpolated daily weather data from file
    IN:
        ID - sve catchment ID. set ID=0 if all data wanted
        start_date - 'yyyy-mm-dd'
        end_date - 'yyyy-mm-dd'
        sourcefile - optional
        CO2 - atm. CO2 concentration (float), optional
        U - wind speed, optional
    OUT:
        fmi - pd.dataframe with datetimeindex
    """

    # OmaTunniste;OmaItä;OmaPohjoinen;Kunta;siteid;vuosi;kk;paiva;longitude;latitude;t_mean;t_max;t_min;
    # rainfall;radiation;hpa;lamposumma_v;rainfall_v;lamposumma;lamposumma_cum
    # -site number
    # -date (yyyy mm dd)
    # -latitude (in KKJ coordinates, metres)
    # -longitude (in KKJ coordinates, metres)
    # -T_mean (degrees celcius)
    # -T_max (degrees celcius)
    # -T_min (degrees celcius)
    # -rainfall (mm)
    # -global radiation (per day in kJ/m2)
    # -H2O partial pressure (hPa)

    sourcefile = os.path.join(sourcefile)

    ID = int(ID)

    # import forcing data
    try:
        fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                          usecols=['OmaTunniste', 'Kunta', 'aika','vuosi','kk','paiva',
                          'longitude','latitude', 't_mean', 't_max', 't_min', 'rainfall',
                          'radiation', 'hpa', 'lamposumma_v', 'rainfall_v'],
                          parse_dates=['aika'],encoding="ISO-8859-1")

        fmi['aika'] = pd.to_datetime({'year': fmi['vuosi'],
                                    'month': fmi['kk'],
                                    'day': fmi['paiva']})

        fmi = fmi.rename(columns={'aika': 'date',
                                  'OmaTunniste': 'ID',
                                  't_mean': 'air_temperature',
                                  'rainfall': 'precipitation',
                                  'radiation': 'global_radiation',
                                  'hpa': 'h2o'})

        time = fmi['date']
    except:
        try:
            fmi = pd.read_csv(sourcefile, sep=';', header='infer',
                              usecols=['x','y','date','temp_avg','prec',
                              'wind_speed_avg','global_rad','vapour_press'],
                              parse_dates=['date'],encoding="ISO-8859-1")

            fmi = fmi.rename(columns={'temp_avg': 'air_temperature',
                                      'prec': 'precipitation',
                                      'global_rad': 'global_radiation',
                                      'vapour_press': 'h2o',
                                      'wind_speed_avg':'wind_speed'})
            time = pd.to_datetime(fmi['date'], format='%Y-%m-%d')
        except:
            raise ValueError('Problem reading forcing data')

    fmi.index = time
    # get desired period and catchment
    fmi = fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]

    if ID > 0:
        fmi = fmi[fmi['ID'] == ID]

    fmi['h2o'] = 1e-1*fmi['h2o']  # hPa-->kPa
    fmi['global_radiation'] = 1e3 / 86400.0*fmi['global_radiation']  # kJ/m2/d-1 to Wm-2
    fmi['par'] = 0.45*fmi['global_radiation']

    # saturated vapor pressure
    esa = 0.6112*np.exp(
            (17.67*fmi['air_temperature']) / (fmi['air_temperature'] + 273.16 - 29.66))  # kPa
    vpd = esa - fmi['h2o']  # kPa
    vpd[vpd < 0] = 0.0
    rh = 100.0*fmi['h2o'] / esa
    rh[rh < 0] = 0.0
    rh[rh > 100] = 100.0

    fmi['RH'] = rh
    fmi['esa'] = esa
    fmi['vapor_pressure_deficit'] = vpd

    fmi['doy'] = fmi.index.dayofyear
    # replace nan's in prec with 0.0
    fmi['precipitation'] = fmi['precipitation'].fillna(0.0)

    # add CO2 and wind speed concentration to dataframe
    if 'CO2' not in fmi:
        fmi['CO2'] = float(CO2)
    if 'wind_speed' not in fmi:
        fmi['wind_speed'] = float(U)

    fmi['wind_speed'] = fmi['wind_speed'].fillna(U)

#    print("NaN values in forcing data:")
#    print(fmi.isnull().any())

    dates = pd.date_range(start_date, end_date).tolist()
#    fmi = fmi.drop_duplicates(keep='first')
#    print(fmi[fmi.duplicated()])
    if len(dates) != len(fmi):
        print(str(len(dates) - len(fmi)) + ' days missing from forcing file, interpolated')
    forcing = pd.DataFrame(index=dates, columns=[])
    forcing = forcing.merge(fmi, how='outer', left_index=True, right_index=True)
    forcing = forcing.fillna(method='ffill')

    return forcing


def read_FMI_weather(start_date, end_date, sourcefile, U=2.0, ID=1, CO2=380.0):
    """
    Reads FMI observed daily weather data from file.

    Computes VPD from temperature and vapour pressure if not already present.

    Args:
        start_date (str):   Start date 'yyyy-mm-dd'.
        end_date   (str):   End date 'yyyy-mm-dd'.
        sourcefile (str):   Path to the CSV forcing file.
        U         (float):  Default wind speed [m s-1] if not in file. Default 2.0.
        ID          (int):  Unused; retained for interface consistency.
        CO2       (float):  Atmospheric CO2 concentration [ppm]. Default 380.

    Returns:
        fmi (pd.DataFrame): Daily forcing with DatetimeIndex and columns:
            'air_temperature'        [degC]
            'precipitation'          [mm d-1]
            'global_radiation'       [W m-2]
            'vapor_pressure_deficit' [kPa]
            'par'                    [W m-2]
            'wind_speed'             [m s-1]
            'CO2'                    [ppm]
            'doy'                    [-]
    """

    # OmaTunniste;OmaItä;OmaPohjoinen;Kunta;siteid;vuosi;kk;paiva;longitude;latitude;t_mean;t_max;t_min;
    # rainfall;radiation;hpa;lamposumma_v;rainfall_v;lamposumma;lamposumma_cum
    # -site number
    # -date (yyyy mm dd)
    # -latitude (in KKJ coordinates, metres)
    # -longitude (in KKJ coordinates, metres)
    # -T_mean (degrees celcius)
    # -T_max (degrees celcius)
    # -T_min (degrees celcius)
    # -rainfall (mm)
    # -global radiation (per day in kJ/m2)
    # -H2O partial pressure (hPa)

    sourcefile = os.path.join(sourcefile)
    #print('*** Simulation forced with:', sourcefile)
    ID = int(ID)

    # import forcing data
    fmi = pd.read_csv(sourcefile, sep=';', header='infer', index_col=0,
                      parse_dates=True ,encoding="ISO-8859-1")

    if 'PAR' not in fmi.columns:
        fmi['PAR'] = 0.5 * fmi['radiation']

    if 'hpa' in fmi.columns:
        fmi['h2o'] = 1e-3*fmi['hpa']  # -> kPa

    if not any(col in fmi.columns for col in ['vpd', 'VPD', 'vapor_pressure_deficit']):    
        # saturated vapor pressure
        esa = 0.6112*np.exp((17.67*fmi['t_mean']) / (fmi['t_mean'] + 273.16 - 29.66))  # kPa
        vpd = esa - fmi['h2o']  # kPa
        vpd[vpd < 0] = 0.0
        fmi['vpd'] = vpd


    fmi = fmi.rename(columns={'t_mean': 'air_temperature', 't_max': 'Tmax',
                              't_min': 'Tmin', 'rainfall': 'precipitation',
                              'radiation': 'global_radiation', 'lamposumma_v': 'dds', 
                              'rh': 'relative_humidity', 'vpd': 'vapor_pressure_deficit', 'PAR':'par'})
    fmi.index.names = ['date']

    # get desired period and catchment
    fmi = fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]

    fmi['doy'] = fmi.index.dayofyear
    # replace nan's in prec with 0.0
    #fmi.loc[fmi['precipitation'].isna(), 'Prec'] = 0.0
    if 'wind_speed' not in fmi:
        fmi['wind_speed'] = float(U)
    fmi['wind_speed'] = fmi['wind_speed'].fillna(float(U))
    # add CO2 concentration to dataframe
    fmi['CO2'] = float(CO2)
    
    return fmi

def initialize_netcdf(pgen, cmask, filepath, filename, description, gisinfo):
    """
    Initializes a NetCDF4 output file for the main simulation results.

    Creates dimensions (time, lat, lon), coordinate variables, and empty
    data variables for all entries listed in pgen['variables'].

    Args:
        pgen        (dict): General parameters; uses 'variables', 'end_date',
                            'spatial_forcing'.
        cmask      (array): Catchment mask defining grid shape.
        filepath     (str): Directory path for the output file.
        filename     (str): Output filename (e.g. '20240101_results.nc').
        description  (str): Description string written to file metadata.
        gisinfo     (dict): Grid info with 'xllcorner', 'yllcorner', 'dxy'.

    Returns:
        ncf (Dataset): Open netCDF4 Dataset handle.
        ff    (str):   Full path to the created file.
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    # dimensions
    date_dimension = None
    lat_shape, lon_shape = np.shape(cmask)

    xllcorner = gisinfo['xllcorner']
    yllcorner = gisinfo['yllcorner']
    cellsize = gisinfo['dxy']
    
    #xcoords = np.arange(xllcorner, (xllcorner + (lon_shape*cellsize)), cellsize)
    #xcoords = np.arange(xllcorner, (xllcorner + (lon_shape*cellsize)-cellsize), cellsize) # ?????
    #ycoords = np.arange(yllcorner, (yllcorner + (lat_shape*cellsize)), cellsize)
    #ycoords = np.arange(yllcorner, (yllcorner + (lat_shape*cellsize+cellsize)), cellsize) # ?????
    xcoords = np.linspace(xllcorner, xllcorner + (lon_shape - 1) * cellsize, lon_shape)
    ycoords = np.linspace(yllcorner, yllcorner + (lat_shape - 1) * cellsize, lat_shape)
    ycoords = np.flip(ycoords)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)
    
    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'SpaFHy results : ' + description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'modified SpaFHy'

    ncf.createDimension('time', date_dimension)
    ncf.createDimension('lat', lat_shape)
    ncf.createDimension('lon', lon_shape)

    date = ncf.createVariable('time', 'f8', ('time',))
    date.units = 'days since 0001-01-01 00:00:00.0'
    date.calendar = 'standard'
    tvec = pd.date_range(pgen['spinup_end'], pgen['end_date']).tolist()[1:]
    date[:] = date2num(tvec, units=date.units, calendar=date.calendar)

    ivar = ncf.createVariable('lat', 'f8', ('lat',))
    ivar.units = 'ETRS-TM35FIN'
    ivar[:] = ycoords

    jvar = ncf.createVariable('lon', 'f8', ('lon',))
    jvar.units = 'ETRS-TM35FIN'
    jvar[:] = xcoords

    for var in pgen['variables']:

        var_name = var[0]
        var_unit = var[1]

        if (var_name.split('_')[0] == 'forcing' and
            pgen['spatial_forcing'] == False):
            var_dim = ('time')
        elif (var_name.split('_')[0] == 'top' and var_name.split('_')[1] != 'local'):
            var_dim = ('time')
        elif var_name.split('_')[0] == 'parameters':
            var_dim = ('lat', 'lon')
        else:
            var_dim = ('time','lat', 'lon')

        variable = ncf.createVariable(
                var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff


def initialize_netcdf_spinup(pgen, cmask, filepath, filename, description, gisinfo):
    """
    Initializes a NetCDF4 output file for saving end-of-spinup model state.

    Saves a single timestep (the last day of the spinup period) for the
    state variables needed to restart a simulation. The filename '_spinup.nc'
    suffix is added automatically.

    Args:
        pgen        (dict): General parameters; uses 'simtype', 'end_date'.
        cmask      (array): Catchment mask defining grid shape.
        filepath     (str): Directory path for the output file.
        filename     (str): Base filename (suffix replaced with '_spinup.nc').
        description  (str): Description string written to file metadata.
        gisinfo     (dict): Grid info with 'xllcorner', 'yllcorner', 'dxy'.

    Returns:
        ncf (Dataset): Open netCDF4 Dataset handle.
        ff    (str):   Full path to the created file.
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    filename = filename[0:-3]+'_spinup.nc'
    # dimensions
    date_dimension = None
    lat_shape, lon_shape = np.shape(cmask)

    xllcorner = gisinfo['xllcorner']
    yllcorner = gisinfo['yllcorner']
    cellsize = gisinfo['dxy']

    xcoords = np.linspace(xllcorner, xllcorner + (lon_shape - 1) * cellsize, lon_shape)
    ycoords = np.linspace(yllcorner, yllcorner + (lat_shape - 1) * cellsize, lat_shape)
    ycoords = np.flip(ycoords)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'SpaFHy results : ' + description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'modified SpaFHy'

    ncf.createDimension('time', date_dimension)
    ncf.createDimension('lat', lat_shape)
    ncf.createDimension('lon', lon_shape)

    date = ncf.createVariable('time', 'f8', ('time',))
    date.units = 'days since 0001-01-01 00:00:00.0'
    date.calendar = 'standard'
    tvec = pd.to_datetime(pgen['end_date'])
    date[:] = date2num(tvec, units=date.units, calendar=date.calendar)

    ivar = ncf.createVariable('lat', 'f8', ('lat',))
    ivar.units = 'ETRS-TM35FIN'
    ivar[:] = ycoords

    jvar = ncf.createVariable('lon', 'f8', ('lon',))
    jvar.units = 'ETRS-TM35FIN'
    jvar[:] = xcoords

    # 1D run
    if pgen['simtype'] == '1D':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]']
                ]
            }
    # 2D run
    elif pgen['simtype'] == '2D':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]'],
                ['soil_ground_water_level', 'ground water level [m]']
                ]
            }
    # TOP run
    elif pgen['simtype'] == 'TOP':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]'],
                ['top_saturation_deficit', 'topmodel saturation deficit [m]]']
                ]
            }

    for var in state_variables['variables']:

        var_name = var[0]
        var_unit = var[1]

        if (var_name.split('_')[0] == 'top' and var_name.split('_')[1] != 'local'):
            var_dim = ('time')
        else:
            var_dim = ('time','lat', 'lon')
        variable = ncf.createVariable(
                var_name, 'f4', var_dim)

        variable.units = var_unit

    return ncf, ff


def write_ncf(results, ncf, steps=None):
    """
    Writes model simulation results to an open NetCDF4 file.

    Args:
        results (dict):          Result arrays keyed by variable name.
        ncf     (Dataset):       Open netCDF4 Dataset handle.
        steps   (tuple or None): (start, end) time indices for writing a slice.
                                 If None, writes the full array.
    """

    keys = results.keys()
    variables = ncf.variables.keys()

    for key in keys:

        if key in variables and key != 'time':
            if len(ncf[key].shape) > 2:
                if steps==None:
                    ncf[key][:,:,:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1],:,:] = results[key][0:steps[1]-steps[0],:,:]
            elif len(ncf[key].shape) > 1:
                ncf[key][:,:] = results[key]
            elif len(ncf[key].shape) == 1:
                ncf[key][steps[0]:steps[1]] = results[key][0:steps[1]-steps[0]]
            else:
                if steps==None:
                    ncf[key][:] = results[key]
                else:
                    ncf[key][steps[0]:steps[1]] = results[key][0:steps[1]-steps[0]]

def write_ncf_spinup(results, pgen, ncf_spinup, steps=None):
    """
    Writes the final spinup state to an open NetCDF4 spinup file.

    Saves only the last timestep of each state variable, based on simtype.

    Args:
        results    (dict):    Result arrays keyed by variable name.
        pgen       (dict):    General parameters; uses 'simtype'.
        ncf_spinup (Dataset): Open netCDF4 Dataset handle for the spinup file.
        steps      :          Unused; retained for interface consistency.
    """


    # 1D run
    if pgen['simtype'] == '1D':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]']
                ]
            }
    # 2D run
    elif pgen['simtype'] == '2D':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]'],
                ['soil_ground_water_level', 'ground water level [m]']
                ]
            }
    # TOP run
    elif pgen['simtype'] == 'TOP':
        state_variables = {'variables':[
                ['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                ['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                ['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]'],
                ['top_saturation_deficit', 'topmodel saturation deficit [m]]']
                ]
            }

    variables = ncf_spinup.variables.keys()
    for key in state_variables['variables']:
        var = key[0]
        if var in variables and var != 'time':
            if len(ncf_spinup[var].shape) > 2:
                ncf_spinup[var][:,:,:] = results[var][-1]
            elif len(ncf_spinup[var].shape) > 1:
                ncf_spinup[var][:,:] = results[var][-1]
            elif len(ncf_spinup[var].shape) == 1:
                ncf_spinup[var][:] = results[var][-1]
            else:
                ncf_spinup[var][:] = results[var][-1]


def read_AsciiGrid(fname, setnans=True):
    """
    Reads an ESRI ASCII raster grid file.

    Expected file format:
        ncols         750
        nrows         375
        xllcorner     350000
        yllcorner     6696000
        cellsize      16
        NODATA_value  -9999
        <data rows>

    Args:
        fname    (str):  Full path to the .asc or .dat file.
        setnans (bool):  If True, replaces NODATA values with np.nan. Default True.

    Returns:
        data     (array): 2D numpy array of raster values.
        info     (list):  First 6 header lines as strings.
        (xloc, yloc) (tuple): Lower-left corner coordinates [m].
        cellsize (float): Cell size [m].
        nodata   (float): NoData value (np.nan if setnans=True).
    """
    fid = open(fname, 'r')
    info = fid.readlines()[0:6]
    fid.close()

    # print info
    # conversion to float is needed for non-integers read from file...
    xloc = float(info[2].split(' ')[-1])
    yloc = float(info[3].split(' ')[-1])
    cellsize = float(info[4].split(' ')[-1])
    nodata = float(info[5].split(' ')[-1])

    # read rest to 2D numpy array
    data = np.loadtxt(fname, skiprows=6)

    if setnans:
        data[data == nodata] = np.nan
        nodata = np.nan

    data = np.array(data, ndmin=2)

    return data, info, (xloc, yloc), cellsize, nodata


def write_AsciiGrid(fname, data, info, fmt='%.18e'):
    """
    Writes a 2D numpy array to an ESRI ASCII raster grid file.

    Args:
        fname (str):   Full path for the output file.
        data (array):  2D numpy array to write (NaN values replaced by nodata).
        info  (list):  6-line header list from read_AsciiGrid.
        fmt    (str):  Format string for numpy.savetxt. Default '%.18e'.
    """

    # replace nans with nodatavalue according to info
    nodata = int(info[-1].split(' ')[-1])
    data[np.isnan(data)] = nodata
    # write info
    fid = open(fname, 'w')
    fid.writelines(info)
    fid.close()

    # write data
    fid = open(fname, 'a')
    np.savetxt(fid, data, fmt=fmt, delimiter=' ')
    fid.close()


def read_results(outputfile):
    """
    Opens simulation results netcdf4 dataset in xarray
    Args:
        outputfile (str): outputfilename
    Returns:
        results (xarray): simulation results from given outputfile
    """

    import xarray as xr

    result = xr.open_dataset(outputfile)
    try:
        result.coords['lat'] = -np.arange(0,result.dims['lat'])
    except KeyError:
        result.coords['i'] = -np.arange(0,result.dims['i'])
    try:
        result.coords['lon'] = np.arange(0,result.dims['lon'])
    except KeyError:
        result.coords['j'] = np.arange(0,result.dims['j'])

    return result

def create_input_GIS(fpath, plotgrids=False):
    """
    Reads raw Finnish forest inventory rasters, converts them to model inputs,
    and writes the results as ESRI ASCII grid files to fpath/inputs/.

    Converts needle/leaf biomass to LAI using species-specific SLA values,
    clips ditch spacing to [20, 200] m, and masks all grids to peat areas only.

    Args:
        fpath      (str):  Path to the folder containing raw input rasters.
        plotgrids (bool):  If True, plots all output grids. Default False.
    """
    fpath = os.path.join(workdir, fpath)

    # specific leaf area (m2/kg) for converting leaf mass to leaf area
    SLA = {'pine': 6.8, 'spruce': 4.7, 'decid': 14.0}  # Härkönen et al. 2015 BER 20, 181-195

    # mask, cmask == 1, np.nan outside
    cmask, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'cmask.asc'))

    # latitude, longitude arrays
    nrows, ncols = np.shape(cmask)
    pos = [round(poss) for poss in pos]
    lon0 = np.arange(pos[0], pos[0] + cellsize*ncols, cellsize)
    lat0 = np.arange(pos[1], pos[1] + cellsize*nrows, cellsize)
    lat0 = np.flipud(lat0)  # indexes count from top left, pos is bottom left

    cmask[np.isfinite(cmask)] = 1.0
    ix = np.where(cmask == 1.0)
    rows = [min(ix[0]), max(ix[0])+2]
    cols = [min(ix[1]), max(ix[1])+2]

    cmask = cmask[rows[0]:rows[-1],cols[0]:cols[-1]]

    lat = lat0[rows[0]:rows[-1]]
    lon = lon0[cols[0]:cols[-1]]

    # peat (only peat areas simulated)
    peat, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'suo_16m.asc'))
    r, c = np.shape(peat)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    peat = peat[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    peat[np.isfinite(peat)] = 1.0

    cmask = cmask * peat

    # needle/leaf biomass to LAI
    bmleaf_pine, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_manty_neulaset.asc'))
    r, c = np.shape(bmleaf_pine)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_pine = bmleaf_pine[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    bmleaf_spruce, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_kuusi_neulaset.asc'))
    r, c = np.shape(bmleaf_spruce)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_spruce = bmleaf_spruce[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    bmleaf_decid, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'bm_lehtip_neulaset.asc'))
    r, c = np.shape(bmleaf_decid)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    bmleaf_decid = bmleaf_decid[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]

    LAI_pine = 1e-3*bmleaf_pine*SLA['pine']  # 1e-3 converts 10kg/ha to kg/m2
    LAI_spruce = 1e-3*bmleaf_spruce*SLA['spruce']
    LAI_decid = 1e-3*bmleaf_decid*SLA['decid']


    # tree height
    hc, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'keskipituus.asc'))
    r, c = np.shape(hc)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    hc = hc[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    hc = 0.1*hc # m

    # canopy closure
    cf, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'latvuspeitto.asc'))
    r, c = np.shape(cf)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    cf = cf[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    cf = 1e-2*cf

    # ditch depth
    ditch_depth = cmask * 0.8

    # ditch spacing
    ditch_spacing, _, pos, cellsize, _ = read_AsciiGrid(os.path.join(fpath, 'ditch_spacing.asc'))
    r, c = np.shape(ditch_spacing)
    pos = [round(poss) for poss in pos]
    ix = r - (lat0[rows] - pos[1]) / cellsize
    iy = (lon0[cols] - pos[0]) / cellsize
    ditch_spacing = ditch_spacing[int(ix[0]):int(ix[1]),int(iy[0]):int(iy[-1])]
    ditch_spacing = np.minimum(200., ditch_spacing)
    ditch_spacing = np.maximum(20., ditch_spacing)

    # soil_id
    soil_id = cmask * 2.0

    nrows, ncols = np.shape(cmask)
    info = ['ncols         ' + str(nrows) + '\n',
            'nrows         ' + str(ncols) + '\n',
            'xllcorner     ' + str(lon[0]) + '\n',
            'yllcorner     ' + str(lat[0]) + '\n',
            'cellsize      ' + str(cellsize) + '\n',
            'NODATA_value  -9999\n']

    # dict of all rasters
    GisData = {'cmask': cmask, 'ditch_spacing': ditch_spacing * cmask, 'ditch_depth': ditch_depth * cmask,
               'LAI_pine': LAI_pine * cmask, 'LAI_spruce': LAI_spruce * cmask, 'LAI_decid': LAI_decid * cmask,
               'hc': hc * cmask, 'cf': cf * cmask, 'soil_id': soil_id}

    if plotgrids is True:
        xy = np.meshgrid(lon, lat)

        plt.figure(99,figsize=(12, 12))
        i=1
        for key, gdata in GisData.items():
            if key != 'cmask':
                plt.subplot(3,3,i)
                plt.pcolor(xy[0], xy[1], gdata)
                plt.colorbar()
                plt.title(key)
                i+=1

    fpath = os.path.join(fpath, 'inputs')
    if not os.path.isdir(fpath):
        os.mkdir(fpath)

    for key, gdata in GisData.items():
        write_AsciiGrid(os.path.join(fpath, key + '.dat'), gdata, info, fmt='%.6e')

def rw_FMI_files(sourcefiles, out_path, plot=False):
    """
    Reads multiple FMI interpolated daily weather CSV files, splits them by
    site, and writes one file per site to out_path.

    Args:
        sourcefiles (list): List of paths to FMI CSV source files.
        out_path     (str): Output directory path.
        plot        (bool): If True, plots each site's time series. Default False.

    Returns:
        fmi (pd.DataFrame): Combined dataframe of all sites and dates.
    """
    frames = []
    for sourcefile in sourcefiles:
        sourcefile = os.path.join(sourcefile)

        # import forcing data
        try:
            fmi = pd.read_csv(sourcefile, sep=',', header='infer',index_col=False,
                              usecols=['pvm','latitude','longitude','t_mean','t_max','t_min',
                                       'rainfall','radiation','hpa','site'],
                              parse_dates=['pvm'],encoding="ISO-8859-1")

            fmi = fmi.rename(columns={'pvm': 'date',
                                      't_mean': 'temp_avg',
                                      't_max': 'temp_max',
                                      't_min': 'temp_min',
                                      'rainfall': 'prec',
                                      'radiation': 'global_rad',
                                      'hpa': 'vapour_press',
                                      'longitude':'x',
                                      'latitude':'y'})
            fmi = fmi[fmi['date']<'2016-07-03']
        except:
            try:
                fmi = pd.read_csv(sourcefile, sep=',', header='infer',index_col=False,
                                  usecols=['x','y','date','temp_avg','temp_min','temp_max',
                                           'prec', 'wind_speed_avg','global_rad','vapour_press',
                                           'snow_depth','pot_evap','site'],
                                  parse_dates=['date'],encoding="ISO-8859-1")

                fmi = fmi.rename(columns={})
            except:
                raise ValueError('Problem reading forcing data')

        time = pd.to_datetime(fmi['date'], format='%Y-%m-%d')
        fmi.index=time

        frames.append(fmi.copy())

    fmi = pd.concat(frames, sort=False)

    sites = list(set(fmi['site']))
    sites.sort()
    index = 0
    readme = 'Indices of weather files'
    for site in sites:
        fmi[fmi['site']==site].to_csv(path_or_buf=out_path + 'weather_id_' + str(index) + '.csv', sep=';', na_rep='NaN', index=False)
        readme += '\n'+ str(index) + ':' + site
        index+=1
        if plot:
            fmi[fmi['site']==site].plot(subplots=True)
    outF = open(out_path + "weather_readme.txt", "w")
    print(readme, file=outF)
    outF.close()
    return fmi


def stitch_result_nc_files(root_directory, output_file, plot=False):
    """
    Merges multiple sub-catchment NetCDF result files into a single file.

    Walks root_directory recursively to find all .nc files, determines the
    union of lat/lon extents, and fills a new output NetCDF with data from
    each sub-catchment file. Cells not covered by any file remain NaN.

    Args:
        root_directory (str):  Root folder containing sub-catchment result .nc files.
        output_file    (str):  Path for the merged output .nc file.
        plot          (bool):  If True, plots bucket_moisture_root for each file
                               and the final merged result. Default False.
    """
    
    import xarray as xr

    def extract_lats_lons(nc_file):
        """Extracts latitudes and longitudes from a NetCDF file."""
        with xr.open_dataset(nc_file) as ds:
            lat = ds['lat'].values
            lon = ds['lon'].values
            cellsize = np.float32(np.abs(ds['lat'][1]-ds['lat'][0]))
        return lat, lon, cellsize

    def extract_time_and_variables(nc_file):
        """Extracts time dimensions and data variable names from the first NetCDF file."""
        with xr.open_dataset(nc_file) as ds:
            time = ds['time'].values
            variables = {var: ds[var].dims for var in ds.data_vars.keys()}
        return time, variables

    def find_global_lats_lons(root_directory):
        """Finds the latitude and longitude coordinates and extracts time and variable dimensions from the first NetCDF file."""
        all_latitudes = set()
        all_longitudes = set()
        time = None
        variables = None

        for dirpath, _, filenames in os.walk(root_directory):
            for file in filenames:
                if file.endswith(".nc"):
                    nc_file_path = os.path.join(dirpath, file)
                    lat, lon, cellsize = extract_lats_lons(nc_file_path)

                    # Collect all unique latitudes and longitudes
                    all_latitudes.update(lat)
                    all_longitudes.update(lon)

                    # Extract time and variable dimensions from the first file
                    if time is None:
                        time, variables = extract_time_and_variables(nc_file_path)

        # Convert sets to sorted arrays
        latitudes = np.sort(np.array(list(all_latitudes)))
        longitudes = np.sort(np.array(list(all_longitudes)))

        latitudes = np.arange(min(latitudes), max(latitudes) + cellsize, cellsize)
        longitudes = np.arange(min(longitudes), max(longitudes) + cellsize, cellsize)

        return latitudes, longitudes, time, variables

    def create_new_nc_file(temp_file, time, latitudes, longitudes, variables, plot):
        """Creates a new NetCDF file with time, latitude, and longitude dimensions, and retains original variable dimensions."""
        data_vars = {}

        # Initialize data variables based on dimensions from the first .nc file
        for dirpath, _, filenames in os.walk(root_directory):
            for file in filenames:
                if file.endswith(".nc"):
                    nc_file_path = os.path.join(dirpath, file)
                    with xr.open_dataset(nc_file_path) as ds:
                        for var in ds.data_vars.keys():
                            var_dims = ds[var].dims

                            # Create an empty array with NaNs for the variable based on dimensions
                            shape = [len(time) if dim == 'time' else len(latitudes) if dim == 'lat' else len(longitudes) if dim == 'lon' else 1 for dim in var_dims]
                            empty_array = np.full(shape, np.nan)

                            # Store the empty array with the variable's dimensions
                            if var not in data_vars:
                                data_vars[var] = (var_dims, empty_array)

        # Create the new NetCDF file
        with xr.Dataset(
            data_vars={var: (dims, data) for var, (dims, data) in data_vars.items()},
            coords={
                "time": time,
                "lat": latitudes,
                "lon": longitudes
            }
        ) as ds:
            ds.to_netcdf(temp_file)

    # Step 1: Find global latitudes, longitudes, time, and variables
    print('checking global lats and lons')
    latitudes, longitudes, time, variables = find_global_lats_lons(root_directory)
    print('found global lats and lons')

    # Step 2: Create the new NetCDF file
    print('creating the new nc file')
    temp_file = os.path.join(os.path.dirname(output_file), "temp.nc")
    create_new_nc_file(temp_file, time, latitudes, longitudes, variables, plot)
    print('created the new nc file')

    # Step 3: Open the newly created NetCDF file
    new_ds = xr.open_dataset(temp_file)

    i = 0
    # Step 4: Walk through the root_directory and process each .nc file
    for dirpath, _, filenames in os.walk(root_directory):
        for file in filenames:
            if file.endswith('.nc'):
                result_path = os.path.join(dirpath, file)

                # Read the result dataset
                result_ds = xr.open_dataset(result_path)
                if plot:
                    plt.figure(i)
                i += 1
                print('processing file', result_path)
                if (result_ds['lat'].min() < new_ds['lat'].min()) or (result_ds['lat'].max() > new_ds['lat'].max()):
                    print('RESULT LAT OUTSIDE')
                if (result_ds['lon'].min() < new_ds['lon'].min()) or (result_ds['lon'].max() > new_ds['lon'].max()):
                    print('RESULT LON OUTSIDE')

                # Ensure the result_ds has the same time dimension and data variables as new_ds
                assert np.all(new_ds.time == result_ds.time), "Time dimensions do not match"
                assert set(new_ds.data_vars) == set(result_ds.data_vars), "Data variables do not match"

                # Reindex lat and lon once # TEST
                aligned_result_ds = result_ds.reindex(lat=new_ds['lat'], lon=new_ds['lon'], method='nearest') # TEST
                # Apply lat/lon masks to handle values outside original bounds TEST
                lat_mask = (new_ds['lat'] >= result_ds['lat'].min()) & (new_ds['lat'] <= result_ds['lat'].max())
                lon_mask = (new_ds['lon'] >= result_ds['lon'].min()) & (new_ds['lon'] <= result_ds['lon'].max())

                '''
                # TEST WITHOUT LOOP
                # Step 5: Separate variables into those with lat/lon and those without lat/lon dimensions
                lat_lon_vars = [var for var in result_ds.data_vars if 'lat' in result_ds[var].dims and 'lon' in result_ds[var].dims]
                time_only_vars = [var for var in result_ds.data_vars if 'time' in result_ds[var].dims and 'lat' not in result_ds[var].dims and 'lon' not in result_ds[var].dims]
                # 1. Handle variables with lat/lon dimensions (Reindex, mask, and combine)
                if lat_lon_vars:
                    # Reindex lat/lon dimensions once for all variables with these dimensions
                    aligned_result_ds = result_ds[lat_lon_vars].reindex(lat=new_ds['lat'], lon=new_ds['lon'], method='nearest')

                    # Apply lat/lon masks
                    lat_mask = (new_ds['lat'] >= result_ds['lat'].min()) & (new_ds['lat'] <= result_ds['lat'].max())
                    lon_mask = (new_ds['lon'] >= result_ds['lon'].min()) & (new_ds['lon'] <= result_ds['lon'].max())

                    # Mask the reindexed data to handle values outside the original lat/lon bounds
                    masked_result_ds = aligned_result_ds.where(lat_mask & lon_mask, np.nan)

                    # Combine new_ds and masked_result_ds, filling missing values in new_ds
                    new_ds = new_ds.combine_first(masked_result_ds)

                # 2. Handle variables with only the time dimension
                if time_only_vars:
                    # These variables can be directly combined without reindexing
                    time_only_result_ds = result_ds[time_only_vars]
                    new_ds = new_ds.combine_first(time_only_result_ds)
                # TEST WITHOUT LOOP ENDS
                '''
                    
                # Step 5: Iterate over data variables in new_ds
                print('entering the fill loop')
            
                for var in new_ds.data_vars:
                    if var in result_ds:
                        try:
                            masked_result_ds = aligned_result_ds[var].where(lat_mask & lon_mask, np.nan)
                            # Fill values in new_ds with values from aligned_result_ds where not NaN
                            new_ds[var] = xr.where(np.isnan(new_ds[var]), masked_result_ds, new_ds[var])
                        except Exception as e:
                            print(f"Error processing variable '{var}' in file '{result_path}': {e}")
                        if (var == 'bucket_moisture_root') & (plot == True):
                            masked_result_ds[-1].plot(vmin=0, vmax=1)
                print('exited the fill loop')
                


    if plot:
        plt.figure(i)
        new_ds['bucket_moisture_root'][-1].plot(vmin=0, vmax=1)

    new_ds.to_netcdf(output_file)
    os.remove(temp_file)  # Remove the temporary file
            
    print('*** Finished ***')
