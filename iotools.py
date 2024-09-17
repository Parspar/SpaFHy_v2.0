# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:34:37 2016

@author: slauniai & khaahti

"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from soilprofile2D import gwl_Wsto
from koordinaattimuunnos import koordTG
from topmodel import twi as twicalc
import re
import importlib

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
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            org_id
            root_id
    """
    
    fpath = os.path.join(workdir, fpath)

    gis = {}

    # soil classification
    if 'org_id' in spatial_pbu:
        if spatial_pbu['org_id'] == True:
            org_id, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pbu['org_id']))
            gis['org_id'] = org_id
    
    if 'root_id' in spatial_pbu:
        if spatial_pbu['root_id'] == True:        
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
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            canopy_fraction - canopy closure (-)
            canopy_height - mean stand height (m)

    """

    fpath = os.path.join(workdir, fpath)

    gis = {}
    
    # tree height [m]
    if 'canopy_height' in spatial_pcpy:
        if spatial_pcpy['canopy_height'] == True:
            canopy_height, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['canopy_height']))
            gis['canopy_height'] = canopy_height
        
    # canopy closure [-]
    if 'canopy_fraction' in spatial_pcpy:
        if spatial_pcpy['canopy_fraction'] == True:
            canopy_fraction, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['canopy_fraction']))
            gis['canopy_fraction'] = canopy_fraction

    if 'LAI_conif' in spatial_pcpy:
        if spatial_pcpy['LAI_conif'] == True:
            LAI_conif, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_conif']))
            LAI_shrub = 0.1 * LAI_conif
            gis['LAI_conif'] = LAI_conif
            gis['LAI_shrub'] = LAI_shrub

    if 'LAI_decid' in spatial_pcpy:
        if spatial_pcpy['LAI_decid'] == True:
            LAI_decid, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_decid']))
            LAI_grass = 0.5 * LAI_decid
            gis['LAI_decid'] = LAI_decid
            gis['LAI_grass'] = LAI_grass        

    if 'LAI_shrub' in spatial_pcpy:
        if spatial_pcpy['LAI_shrub'] == True:
            LAI_shrub, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_shrub']))
            gis['LAI_shrub'] = LAI_shrub        

    if 'LAI_grass' in spatial_pcpy:
        if spatial_pcpy['LAI_grass'] == True:
            LAI_grass, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pcpy['state']['LAI_grass']))
            gis['LAI_grass'] = LAI_grass 
    
    if plotgrids is True:

        plt.figure()
        plt.subplot(221); plt.imshow(LAI_pine+LAI_spruce); plt.colorbar();
        plt.title('LAI conif (m2/m2)')
        plt.subplot(222); plt.imshow(LAI_decid); plt.colorbar();
        plt.title('LAI decid (m2/m2)')
        plt.subplot(223); plt.imshow(canopy_height); plt.colorbar(); plt.title('hc (m)')
        plt.subplot(224); plt.imshow(canopy_fraction); plt.colorbar(); plt.title('cf (-)')

    return gis

def read_ds_gisdata(fpath, spatial_pspd, mask=None, plotgrids=False):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
            deep_id
            elevation
            lakes
            streams
    """
    fpath = os.path.join(workdir, fpath)

    gis = {}
    
    # deep soil layer
    if 'deep_id' in spatial_pspd:
        if spatial_pspd['deep_id'] == True:
            deep_id, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['deep_id']))
            gis['deep_id'] = deep_id

    # dem
    if 'elevation' in spatial_pspd:
        if spatial_pspd['elevation'] == True:
            elevation, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['elevation']))
            bedrock = elevation - 5.0
            gis['elevation'] = elevation
            gis['bedrock'] = bedrock

    # streams
    if 'streams' in spatial_pspd:
        if spatial_pspd['streams'] == True:
            streams, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['streams']))
            streams[(np.isfinite(streams)) & (streams != 0.0)] = -1.0
            streams[streams != -1.0] = np.nan

    else:
        print('*** No stream file ***')
        streams = np.full_like(deepsoil, 0.0)
    gis['streams'] = streams
    
    # lakes if available
    if 'lakes' in spatial_pspd:
        if spatial_pspd['lakes'] == True:        
            lakes, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, pspd['lakes']))
            lakes[np.isfinite(lakes)] = -1.0
            #akes[lakes == np.nan] = 0.0
            #lakes = np.where(lakes == 0, np.nan, -1.0)
    else:
        print('*** No lakes file ***')
        lakes = np.full_like(deep_id, 0.0)
    gis['lakes'] = lakes

    xllcorner = int(re.findall(r'\d+', info[2])[0])
    yllcorner = int(re.findall(r'\d+', info[3])[0])

    if plotgrids is True:
        plt.figure()
        plt.subplot(311); plt.imshow(soilclass); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(elevation); plt.colorbar(); plt.title('elevation')
        plt.subplot(313); plt.imshow(elevation); plt.colorbar();

    gis.update({'dxy': cellsize})
    gis.update({'xllcorner': xllcorner,
                'yllcorner': yllcorner})
    return gis

def read_top_gisdata(fpath, spatial_ptop, mask=None, plotgrids=False):
    """
    reads gis-data grids and returns numpy 2d-arrays
    Args:
        fpath - relative path to data folder (str)
        plotgrids - True plots
    Returns:
        gis - dict of gis-data rasters
        flowacc - flow accumulation raster
        slope - slope raster
        twi - topographic wetness index
        cmask - catchment mask

    """
    fpath = os.path.join(workdir, fpath)

    gis = {}
    
    # flow accumulation
    if 'flow_accumulation' in spatial_ptop:
        if spatial_ptop['flow_accumulation'] == True:        
            flowacc, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['flow_accumulation']))
            gis['flowacc'] = flowacc
        
    # slope
    if 'slope' in spatial_ptop:
        if spatial_ptop['slope'] == True:                
            slope, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['slope']))
            gis['slope'] = slope

    if 'twi' in spatial_ptop:
        if spatial_ptop['twi'] == True:                        
            twi, info, _, cellsize, _ = read_AsciiGrid(os.path.join(fpath, ptop['twi']))
            gis['twi'] = twi
    
    if plotgrids is True:
        plt.figure()
        plt.subplot(311); plt.imshow(slope); plt.colorbar(); plt.title('soiltype')
        plt.subplot(312); plt.imshow(twi); plt.colorbar(); plt.title('elevation')
        plt.subplot(313); plt.imshow(flowacc); plt.colorbar();

    gis.update({'dxy': cellsize})

    return gis

def read_aux_gisdata(fpath, spatial_aux, mask=None):

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
    creates input dictionary for initializing SoilGrid
    Args:
        soil parameters
        soiltype parameters
        gisdata
            cmask
            soilclass
        spatial
    """
    # create dict for initializing soil profile.
    # copy pbu into sdata and make each value np.array(np.shape(cmask))
    data = pbu.copy()
    spatial_data = spatial_pbu.copy()
    gridshape = np.ones(shape=gisdata['org_id'].shape) # replace org_id with something generic?

    #  Create a mask where gisdata['org_id'] is finite (not NaN)
    mask = np.isfinite(gisdata['org_id'])

    for key in data:
        if spatial_data[key] == True:
            data[key] = gisdata[key]
        if spatial_data[key] == False:
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
    
    if set(root_ids) >= set(np.unique(data['root_id'][np.isfinite(gisdata['org_id'])]).tolist()):
        # no problems
        print('*** Defined root soil IDs:',set(root_ids), 'Used root soil IDs:',
              set(np.unique(data['root_id'][np.isfinite(gisdata['org_id'])]).tolist()), '***')
    else:
        print('*** Defined root soil IDs:', set(root_ids),  'Used root soil IDs:',
              set(np.unique(data['root_id'][np.isfinite(gisdata['org_id'])]).tolist()), '***')
    #    raise ValueError("Root soil id in inputs not specified in parameters.py")

    if set(org_ids) >= set(np.unique(data['org_id'][np.isfinite(gisdata['org_id'])]).tolist()):
        # no problems
        print('*** Defined organic soil IDs:',set(org_ids), 'Used organic soil IDs:',
              set(np.unique(data['org_id'][np.isfinite(gisdata['org_id'])]).tolist()), '***')
    else:
        print('*** Defined organic soil IDs:', set(org_ids), 'Used organic soil IDs:',
              set(np.unique(data['org_id'][np.isfinite(gisdata['org_id'])]).tolist()), '***')
    #    raise ValueError("Org soil id in inputs not specified in parameters.py")

    if spatial == True:
        for key, value in orgp.items():
            t = value['org_id']
            yx = np.where(data['org_id'] == t)
            data['org_depth'][yx] = value['org_depth']
            data['org_poros'][yx] = value['org_poros']
            data['org_fc'][yx] = value['org_fc']
            data['org_rw'][yx] = value['org_rw']
            data['org_ksat'][yx] = value['org_ksat']
            data['org_beta'][yx] = value['org_beta']

    if spatial == True:
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

    return data

def preprocess_dsdata(pspd, spatial_pspd, deepp, gisdata, spatial=True):
    """
    creates input dictionary for initializing SoilGrid
    Args:
        soil parameters
        soiltype parameters
        gisdata
            cmask
            soilclass
        spatial
    """
    # create dict for initializing soil profile.
    # copy pbu into sdata and make each value np.array(np.shape(cmask))
    data = pspd.copy()
    spatial_data = spatial_pspd.copy()
    gridshape = np.ones(shape=gisdata['deep_id'].shape)
   
    for key in data:
        if spatial_data[key] == True:
            data[key] = gisdata[key]
        if spatial_data[key] == False:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, uni_value)
    
    if spatial == False:
        data['deep_id'] = pspd['deep_id']
    else:
        data['deep_id'] = gisdata['deep_id']
        data['elevation'] = gisdata['elevation']
        data['bedrock'] = gisdata['bedrock']        
        data['streams'] = gisdata['streams']
        data['lakes'] = gisdata['lakes']

    deep_ids = []
    for key, value in deepp.items():
        deep_ids.append(value['deep_id'])
        
    if set(deep_ids) >= set(np.unique(data['deep_id'][np.isfinite(gisdata['deep_id'])]).tolist()):
        # no problems
        print('*** Defined deep soil IDs:',set(deep_ids), 'Used soil IDs:',
              set(np.unique(data['deep_id'][np.isfinite(gisdata['deep_id'])]).tolist()))
    else:
        print(set(deep_ids),set(np.unique(data['deep_id'][np.isfinite(gisdata['deep_id'])]).tolist()))
        #raise ValueError("Deep soil id in inputs not specified in parameters.py")

    data.update({'soiltype': np.empty(np.shape(gisdata['deep_id']),dtype=object)})

    for key, value in deepp.items():
        c = value['deep_id']
        ix = np.where(data['deep_id'] == c)
        data['soiltype'][ix] = key
        # interpolation function between wsto and gwl
        value.update(gwl_Wsto(value['deep_z'], value['pF'], value['deep_ksat']))
        # interpolation function between root_wsto and gwl
        value.update(gwl_Wsto(value['deep_z'][:2], {key: value['pF'][key][:2] for key in value['pF'].keys()}, root=True))

    # stream depth corresponding to assigned parameter
    #data['streams'] = np.where((data['streams'] < -eps) | (data['lakes'] < -eps), pspd['stream_depth'], 0)
    data['streams'] = np.where(data['streams'] < -eps, pspd['stream_depth'], 0)
    data['lakes'] = np.where(data['lakes'] < -eps, pspd['lake_depth'], 0)
    #data['streams'] = np.where(data['lakes'] < -eps, pspd['stream_depth'], 0)
    
    data['wtso_to_gwl'] = {soiltype: deepp[soiltype]['to_gwl'] for soiltype in deepp.keys()}
    data['gwl_to_wsto'] = {soiltype: deepp[soiltype]['to_wsto'] for soiltype in deepp.keys()}
    data['gwl_to_C'] = {soiltype: deepp[soiltype]['to_C'] for soiltype in deepp.keys()}
    data['gwl_to_Tr'] = {soiltype: deepp[soiltype]['to_Tr'] for soiltype in deepp.keys()}
    data['gwl_to_rootmoist'] = {soiltype: deepp[soiltype]['to_rootmoist'] for soiltype in deepp.keys()}
    data['dxy'] = gisdata['dxy']

    return data

def preprocess_cpydata(pcpy, spatial_pcpy, gisdata, spatial=True):
    """
    creates input dictionary for initializing CanopyGrid
    Args:
        canopy parameters
        gisdata
            cmask
            LAI_pine, LAI_spruce - pine and spruce LAI (m2m-2)
            LAI_conif - conifer total annual max LAI (m2m-2)
            LAI_dedid - deciduous annual max LAI (m2m-2)
            canopy_fraction - canopy closure (-)
            canopy_height - mean stand height (m)
            (lat, lon)
        spatial
    """
    # inputs for CanopyGrid initialization: update pcpy using spatial data
    data = pcpy['state'].copy()
    spatial_data = spatial_pcpy.copy()
    gridshape = np.ones(shape=gisdata['LAI_conif'].shape)

    for key in data:
        if spatial_data[key] == True:
            data[key] = gisdata[key]
        if spatial_data[key] == False:
            uni_value = data[key]
            data[key] = np.full_like(gridshape, uni_value)

    pcpy['state'] = data

    return pcpy

def preprocess_topdata(ptop, spatial_ptop, gisdata, spatial=True):
    """
    creates input dictionary for initializing CanopyGrid
    Args:
        topmodel parameters as in parameters.py 'pgen' and:
        gisdata
        flowacc - flow accumulation raster
        slope - slope raster
        twi - topographic wetness index
        cmask - catchment mask
            (lat, lon)
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
    reads FMI OBSERVED daily weather data from file
    IN:
        ID - sve catchment ID. set ID=0 if all data wanted
        start_date - 'yyyy-mm-dd'
        end_date - 'yyyy-mm-dd'
        sourcefile - optional
        CO2 - atm. CO2 concentration (float), optional
    OUT:
        fmi - pd.dataframe with datetimeindex
            fmi columns:['ID','Kunta','aika','lon','lat','T','Tmax','Tmin',
                         'Prec','Rg','h2o','dds','Prec_a','Par',
                         'RH','esa','VPD','doy']
            units: T, Tmin, Tmax, dds[degC], VPD, h2o,esa[kPa],
            Prec, Prec_a[mm], Rg,Par[Wm-2],lon,lat[deg]
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
    print('*** Simulation forced with:', sourcefile)
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
    netCDF4 format output file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        cmask
        filepath: path for saving results
        filename: filename
        description: description
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
    netCDF4 format output spinup file initialization

    Args:
        variables (list): list of variables to be saved in netCDF4
        cmask
        filepath: path for saving results
        filename: filename, spinup notion automatically added
        description: description
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

    xcoords = np.arange(xllcorner, (xllcorner + (lon_shape*cellsize)), cellsize)
    ycoords = np.arange(yllcorner, (yllcorner + (lat_shape*cellsize)), cellsize)
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
    Writes model simultaion results in netCDF4-file

    Args:
        index (int): model loop index
        results (dict): calculation results from group
        ncf (object): netCDF4-file handle
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
    Writes model simultaion results in netCDF4-file

    Args:
        index (int): model loop index
        results (dict): calculation results
        ncf (object): netCDF4-file handle
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
    reads AsciiGrid format in fixed format as below:
        ncols         750
        nrows         375
        xllcorner     350000
        yllcorner     6696000
        cellsize      16
        NODATA_value  -9999
        -9999 -9999 -9999 -9999 -9999
        -9999 4.694741 5.537514 4.551162
        -9999 4.759177 5.588773 4.767114
    IN:
        fname - filename (incl. path)
    OUT:
        data - 2D numpy array
        info - 6 first lines as list of strings
        (xloc,yloc) - lower left corner coordinates (tuple)
        cellsize - cellsize (in meters?)
        nodata - value of nodata in 'data'
    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np
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

    if setnans is True:
        data[data == nodata] = np.NaN
        nodata = np.NaN

    data = np.array(data, ndmin=2)

    return data, info, (xloc, yloc), cellsize, nodata


def write_AsciiGrid(fname, data, info, fmt='%.18e'):
    """ writes AsciiGrid format txt file
    IN:
        fname - filename
        data - data (numpy array)
        info - info-rows (list, 6rows)
        fmt - output formulation coding

    Samuli Launiainen Luke 7.9.2016
    """
    import numpy as np

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

    """
    fpath = os.path.join(workdir, fpath)

    # specific leaf area (m2/kg) for converting leaf mass to leaf area
    SLA = {'pine': 6.8, 'spruce': 4.7, 'decid': 14.0}  # Härkönen et al. 2015 BER 20, 181-195

    # mask, cmask == 1, np.NaN outside
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
    if os.path.isdir(fpath) == False:
        os.mkdir(fpath)

    for key, gdata in GisData.items():
        write_AsciiGrid(os.path.join(fpath, key + '.dat'), gdata, info, fmt='%.6e')

def rw_FMI_files(sourcefiles, out_path, plot=False):
    """
    reads and writes FMI interpolated daily weather data
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
