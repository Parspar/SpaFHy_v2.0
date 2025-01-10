# -*- coding: utf-8 -*-
"""
PARAMETERS
@author: slauniai & khaahti & jpnousu
"""

import pathlib
import time

def parameters(folder=''):

    pgen = {'description': 'final_run',  # description written in result file
            'simtype': '2D', # 1D, TOP, 2D,
            'start_date': '2013-01-01',  # '2011-01-01', for tests: '2020-01-01'
            'end_date': '2015-12-31', # 2017-12-31,
            #'spinup_file': r'F:\SpaFHy_2D_2021/testcase_input_202304051037_spinup.nc',
            'spinup_end': '2014-12-31',  # '2014-12-31', for tests: '2020-09-01' results after this are saved in result file
            'dt': 86400.0,
            'spatial_cpy': True,  # if False uses parameters from cpy['state']
            # else needs cf.dat, hc.dat, LAI_decid.dat, LAI_spruce.dat, LAI_pine.dat, (cmask.dat)
            'spatial_soil': True,  # if False uses soil_id, stream_depth from psp
            'spatial_deep': True,
            'org_drain': True, # organic layer drainage True/False            
            'topmodel': True,
            # else needs soil_id.dat, stream_depth.dat
            'spatial_forcing': False,  # if False uses forcing from forcing file with pgen['forcing_id'] and cpy['loc']
            'spatial_radiation_file': None, # if spatial radiation file, otherwise None
            # else needs Ncoord.dat, Ecoord.dat, forcing_id.dat
            'gis_folder': str(pathlib.Path(folder+r'/gis/25m_3case')),
            'forcing_file': str(pathlib.Path(folder+r'/forcing/FORCING.csv')),
            'forcing_id': 0,  # used if spatial_forcing == False
            'ncf_file': time.strftime('%Y%m%d%H%M') + r'.nc',  # timestamp to result file name to avoid saving problem when running repeatedly
            'cmask' : 'catchment_mask.asc',
            'mask': 18, # 'cmask/streams', 'cmask', 'streams', None
            #'results_folder': r'/scratch/project_2000908/nousu/SpaFHy_RESULTS',
            'results_folder': str(pathlib.Path(folder+r'/results')),
            'save_interval': 366, # interval for writing results to file (decreases need for memory during computation)
            'variables':[ # list of output variables (rows can be commented away if not all variables are of interest)
                    ['parameters_lai_conif', 'leaf area index of conifers [m2 m-2]'],
                    ['parameters_lai_decid_max', 'leaf area index of decidious trees [m2 m-2]'],
                    ['parameters_lai_shrub', 'leaf area index of shrubs [m2 m-2]'],
                    ['parameters_lai_grass', 'leaf area index of grass [m2 m-2]'],
                    ['parameters_canopy_height', 'canopy height [m]'],
                    ['parameters_canopy_fraction', 'canopy closure [-]'],
                    ['parameters_org_id', 'soil class index'],
                    ['parameters_root_id', 'soil class index'],
                    ['parameters_deep_id', 'soil class index'],
                    ['parameters_elevation', 'elevation from dem [m]'],
                    ['parameters_lat', 'latitude [deg]'],
                    ['parameters_lon', 'longitude [deg]'],
                    ['parameters_streams', 'streams'],
                    ['parameters_lakes', 'lakes'],
                    ['parameters_cmask', 'cmask'],
                    ['parameters_twi', 'twi'],
                    ['parameters_slope', 'slope'],
                    ['forcing_air_temperature', 'air temperature [degC]'],
                    ['forcing_relative_humidity', 'relative humidity [%]'],
                    ['forcing_precipitation', 'precipitation [mm d-1]'],
                    ['forcing_vapor_pressure_deficit', 'vapor pressure deficit [kPa]'],
                    ['forcing_global_radiation', 'global radiation [Wm-2]'],
                    ['forcing_wind_speed','wind speed [m s-1]'],
                    ['forcing_wind_direction','wind direction [degrees]'],                
                    #['bucket_pond_storage', 'pond storage [m]'],
                    ['bucket_moisture_top', 'volumetric water content of moss layer [m3 m-3]'],
                    ['bucket_moisture_root', 'volumetric water content of rootzone [m3 m-3]'],
                    #['bucket_psi_root', 'soil water potential of rootzone [MPa]'],                    
                    #['bucket_potential_infiltration', 'potential infiltration [mm d-1]'],
                    ['bucket_surface_runoff', 'surface runoff [mm d-1]'],
                    ['bucket_evaporation', 'evaporation from soil surface [mm d-1]'],
                    ['bucket_drainage', 'drainage from root layer [mm d-1]'],
                    ['bucket_water_storage', 'bucket water storage (top and root) [mm d-1]'],
                    #['bucket_water_storage_top', 'bucket water storage (top) [mm d-1]'],
                    #['bucket_water_storage_root', 'bucket water storage (root) [mm d-1]'],
                    #['bucket_storage_change', 'bucket water storage change (top and root) [mm d-1]'],
                    ['bucket_water_closure', 'bucket water balance error [mm d-1]'],
                    ['bucket_return_flow', 'return flow from deepzone to bucket [mm d-1]'],
                    #['deep_water_storage', 'soil water storage (deeplayer) [m]'],
                    ['deep_ground_water_level', 'ground water level [m]'],
                    ['deep_lateral_netflow', 'subsurface lateral netflow [mm d-1]'],
                    ['deep_netflow_to_ditch', 'netflow to stream [mm d-1]'],
                    ['deep_moisture_deep', 'volumetric water content of deepzone [m3 m-3]'],
                    ['deep_water_closure', 'soil water balance error [mm d-1]'],
                    #['canopy_interception', 'canopy interception [mm d-1]'],
                    ['canopy_evaporation', 'evaporation from interception storage [mm d-1]'],
                    ['canopy_transpiration','transpiration [mm d-1]'],
                    #['canopy_stomatal_conductance','stomatal conductance [m s-1]'],
                    #['canopy_throughfall', 'throughfall to moss or snow [mm d-1]'],
                    ['canopy_snow_water_equivalent', 'snow water equivalent [mm]'],
                    ['canopy_water_closure', 'canopy water balance error [mm d-1]'],
                    #['canopy_phenostate', 'canopy phenological state [-]'],
                    #['canopy_leaf_area_index', 'canopy leaf area index [m2 m-2]'],
                    #['canopy_degree_day_sum', 'sum of degree days [degC]'],
                    #['canopy_fLAI', 'state of LAI'],
                    #['canopy_water_storage', 'canopy intercepted water storage [mm d-1]'],
                    #['canopy_snowfall', 'canopy snowfall'],
                    #['top_baseflow', 'topmodel baseflow [mm d-1]'],
                    #['top_water_closure', 'topmodel water balance error [mm d-1]'],
                    #['top_returnflow', 'topmodel returnflow [mm d-1]'],
                    #['top_local_returnflow', 'topmodel local returnflow [mm d-1]'],
                    #['top_drainage_in', 'topmodel inflow from drainage [mm d-1]'],
                    #['top_saturation_deficit', 'topmodel saturation deficit [m]'],
                    #['top_local_saturation_deficit', 'topmodel local saturation deficit [mm]'],
                    #['top_saturated_area', 'topmodel saturated area [-]'],
                    #['top_storage_change', 'topmodel_water_storage_change [mm d-1]']
                    ]
             }

    f=1.0

    # canopygrid
    pcpy = {
            'flow' : {  # flow field
                     'zmeas': 10.0,
                     'zground': 0.5,
                     'zo_ground': 0.01
                     },
            'interc': {  # interception
                        'wmax': 1.5,  # storage capacity for rain (mm/LAI)
                        'wmaxsnow': 4.5,  # storage capacity for snow (mm/LAI)
                        },
            'snow': {  # degree-day snow model
                    'kmelt': 2.8934e-05,  # melt coefficient in open (mm/s)
                    'kfreeze': 5.79e-6,  # freezing coefficient (mm/s)
                    'r': 0.05  # maximum fraction of liquid in snow (-)
                    },
            'physpara': {  # canopy conductance
                        'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1) # MOVING THESE TO SPEC_PARA
                        'g1_conif': f * 2.1, # stomatal parameter, conifers
                        'g1_decid': f * 3.5, # stomatal parameter, deciduous
                        'g1_shrub': f * 3.0, # stomatal parameter, deciduous
                        'g1_grass': f * 5.0, # stomatal parameter, deciduous
                        'q50': 50.0, # light response parameter (Wm-2)
                        'kp': 0.6, # light attenuation parameter (-)
                        'rw': 0.20, # critical value for REW (-),
                        'rwmin': 0.02, # minimum relative conductance (-)
                        # soil evaporation
                        'gsoil': 1e-2 # soil surface conductance if soil is fully wet (m/s)
                        },
            'spec_para': {
                        'conif': {  'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                                    'g1': 2.1, # stomatal parameter
                                    'q50': 50.0, # light response parameter (Wm-2)
                                    'lai_cycle': False,
                                     },
                        'decid': {  'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                                    'g1': 3.5, # stomatal parameter
                                    'q50': 50.0, # light response parameter (Wm-2)
                                    'lai_cycle': True,
                                     },
                        'shrub': {  'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                                    'g1': 3.0, # stomatal parameter
                                    'q50': 50.0, # light response parameter (Wm-2)
                                    'lai_cycle': False,
                                     },
                        'grass': {  'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                                    'g1': 5.0, # stomatal parameter
                                    'q50': 50.0, # light response parameter (Wm-2)
                                    'lai_cycle': True,
                                     },
                        },
            'phenopara': {
                        # seasonal cycle of physiology: smax [degC], tau[d], xo[degC],fmin[-](residual photocapasity)
                        'smax': 18.5, # degC
                        'tau': 13.0, # days
                        'xo': -4.0, # degC
                        'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                        # deciduos phenology
                        'LAI_decid_min': 0.1, # minimum relative LAI (-)
                        'ddo': 45.0, # degree-days for bud-burst (5degC threshold)
                        'ddur': 23.0, # duration of leaf development (days)
                        'sdl': 9.0, # daylength for senescence start (h)
                        'sdur': 30.0, # duration of leaf senescence (days),
                         },
            'state': {  # spatial_cpy = False -> floats | spatial_cpy = True -> filenames in gispath
                       'LAI_conif': 'LAI_conif.asc', # conifer 1-sided LAI (m2 m-2)
                       'LAI_decid': 'LAI_decid.asc',  # maximum annual deciduous 1-sided LAI (m2 m-2)
                       'LAI_shrub': 0.1,
                       'LAI_grass': 0.2,
                       'canopy_height': 'canopy_height.asc', # canopy height (m)
                       'canopy_fraction': 'canopy_fraction.asc', # canopy closure fraction (-)
                       # initial state of canopy storage [mm] and snow water equivalent [mm]
                       'w': 0.0, # canopy storage mm
                       'swe': 0.0, # snow water equivalent mm
                       },
            'loc': {  # following coordinates used if spatial_forcing == False
                    'lat': 67.995,  # decimal degrees
                    'lon': 24.224
                    }
            }

    # soil profile (bucket)
    pbu = {
            # soil profile, following properties are used if spatial_soil = False
            # organic moss-humus layer
            'org_id': 'soil.asc', # uniform (float) OR path to grid in gispath (str)       
            'org_depth': 0.05, # depth of organic top layer (m)
            'org_poros': 0.448, # porosity (-)
            'org_fc': 0.33, # field capacity (-)
            'org_rw': 0.15, # critical vol. moisture content (-) for decreasing phase in Ef
            'org_ksat': 1E-04, # root zone hydraulic conductivity
            'org_beta': 6.0, # 
            'maxpond': 0.05, # max ponding depth (m)
            # rootzone layer
            'root_id': 'soil.asc', # uniform (float) OR path to grid in gispath (str)     
            'root_depth': 0.3, # depth of rootzone layer (m)
            'root_sat': 0.6, # saturation ratio (-)
            'root_fc': 0.33, # field capacity
            'root_poros': 0.448, #  porosity
            'root_wp': 0.13, # wilting point
            'root_ksat': 1e-05, # hydraulic conductivity
            'root_beta': 4.7, #
            'root_alpha': 4.48, #
            'root_n': 1.20, # 
            'root_wr': 0.0, #
            # initial states
            'org_sat': 1.0, # organic top layer saturation ratio (-)
            'pond_storage': 0.0,  # initial pond depth at surface [m]
            }

    # soil profile (2D, deep)
    pspd = {
            # deep soil profile, following properties are used if spatial_deep = False
            'deep_id': 'soil.asc', # uniform (float) OR path to grid in gispath (str)
            'elevation': 'inflated_dem.asc', # uniform (float) OR path to grid in gispath (str) 
            'streams': 'stream_mask_dig.asc',
            'lakes': 'lake_mask.asc',
            'deep_z': -15.0, #'soildepth.asc',
            'deep_poros': 0.41,
            'deep_wr': 0.05,
            'deep_alpha': 0.024,
            'deep_n': 1.2,
            'deep_ksat': 1E-05,
            # initial states
            'ground_water_level': -2.0,  # groundwater depth [m]
            'stream_depth': 'stream_depth.asc', # 'stream_depth.asc',   # initial stream water level relative to ground surface (currently not dynamic) [m]
            'lake_depth': -1.0  # initial lake water level relative to ground surface (currently not dynamic) [m]
            }

    return pgen, pcpy, pbu, pspd


def ptopmodel():
    """
    parameters of topmodel submodel
    """
    ptopmodel = {
            'dem': 'inflated_dem.asc',
            'flow_accumulation': 'flow_accumulation.asc',
            'slope': 'slope.asc',
            'twi': 'twi.asc',
            'dt': 86400.0, # timestep (s)
            'm': 0.025, # 0.025 calibrated by Samuli, scaling depth (m), testin 0.01
            'ko': 0.001, # transmissivity parameter (ms-1)
            'twi_cutoff': 97.5,  # cutoff of cumulative twi distribution (%)
            'so': 0.05 # initial saturation deficit (m)
           }
    return ptopmodel

def auxiliary_grids():
    """
    paths to auxiliary grids such as cmask, lakes, streams
    """
    grids = {
            'cmask':    'catchment_mask.asc',
            'streams':  'stream_mask_dig.asc',
            'lakes':    'lake_mask.asc'
            }
    return grids

def deep_properties():
    """
    Properties of soil profiles.
    Note z is elevation of lower boundary of layer (soil surface at 0.0),
    e.g. z = [-0.05, -0.15] means first layer tickness is 5 cm and second 10 cm.
    """
    deepp = {
        'Postglacial_sand':{ # Launiainen et al. 2021
            'deep_id': 1,
            'deep_z': [-1.2, -3.8, -4.0, -14.4, -16.0], # -16.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.41,0.41,0.41,0.41,0.41], # Launiainen et al. 2021
                    'ThetaR': [0.05,0.05,0.05,0.05,0.05], # Launiainen et al. 2021
                    'alpha': [0.024,0.024,0.024,0.024,0.024], # Launiainen et al. 2021
                    'n': [1.2,1.2,1.2,1.2,1.2]}, # Launiainen et al. 2021
            'deep_ksat': [1E-07,3E-5,1E-8,3E-5,1E-4],
                },
        'Glaciofluvial_sediment':{ # CoarseText
            'deep_id': 2,
            'deep_z': [-27.0], # -27.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.41], # Launiainen et al. 2021
                    'ThetaR': [0.05], # Launiainen et al. 2021
                    'alpha': [0.024], # Launiainen et al. 2021
                    'n': [1.2]}, # Launiainen et al. 2021
            'deep_ksat': [1E-04],
                },
        'Peat':{
            'deep_id': 3,
            'deep_z': [-5.0, -7.0, -17.0], # -17.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.89,0.89,0.89],
                    'ThetaR': [0.196,0.196,0.196],
                    'alpha': [0.072,0.072,0.072],
                    'n': [1.255,1.255,1.255]}, 
            'deep_ksat': [1E-06,1E-9,1E-6], 
                },
        'Postglacial_sand_gravel':{ # CoarseText
            'deep_id': 4, 
            'deep_z': [-0.8, -3.0, -14.0], # -14.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.41,0.41,0.41], # Launiainen et al. 2021
                    'ThetaR': [0.05,0.05,0.05], # Launiainen et al. 2021
                    'alpha': [0.024,0.024,0.024], # Launiainen et al. 2021
                    'n': [1.2,1.2,1.2]}, # Launiainen et al. 2021
            'deep_ksat': [3E-5, 1E-8, 1E-6],
                },
        'Clay_silt':{ # Fine
            'deep_id': 5,
            'deep_z': [-3, -23.0], # -23.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.6,0.6],
                    'ThetaR': [0.07,0.07],
                    'alpha': [0.018,0.018],
                    'n': [1.16,1.16]},
            'deep_ksat': [1E-08, 1E-6],
                },
        'Washed_sediment_gravel_boulders':{ # CoarseText
            'deep_id': 6,
            'deep_z': [-0.5, -2.5, -13.0], # -13.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.41,0.41,0.41], # Launiainen et al. 2021
                    'ThetaR': [0.05,0.05,0.05], # Launiainen et al. 2021
                    'alpha': [0.024,0.024,0.024], # Launiainen et al. 2021
                    'n': [1.2,1.2,1.2]}, # Launiainen et al. 2021
            'deep_ksat': [2E-06, 2E-6, 1E-6],
                },
        'Water':{ # Medium
            'deep_id': 7,
            'deep_z': [-17.0], # -17.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.43], # Launiainen et al. 2019
                    'ThetaR': [0.05], # Launiainen et al. 2019
                    'alpha': [0.024], # Launiainen et al. 2019
                    'n': [1.2]}, # # Launiainen et al. 2021
            'deep_ksat': [1E-05],
                },        
        'Moraine':{ # Medium
            'deep_id': 8,
            'deep_z': [-0.5, -2.5, -13.0], # -13.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.43,0.43,0.43], # Launiainen et al. 2019
                    'ThetaR': [0.05,0.05,0.05], # Launiainen et al. 2019
                    'alpha': [0.024,0.024,0.024], # Launiainen et al. 2019
                    'n': [1.2,1.2,1.2]}, # # Launiainen et al. 2021
            'deep_ksat': [2E-06, 2E-6, 1E-6],
                },        
        'Fill':{ # Medium
            'deep_id': 9,
            'deep_z': [-0.5, -2.5, -13.0], # -18.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.43,0.43,0.43], # Launiainen et al. 2019
                    'ThetaR': [0.05,0.05,0.05], # Launiainen et al. 2019
                    'alpha': [0.024,0.024,0.024], # Launiainen et al. 2019
                    'n': [1.2,1.2,1.2]}, # # Launiainen et al. 2021
            'deep_ksat': [2E-06, 2E-6, 1E-6],
                }, 
        'Bedrock':{
            'deep_id': 10,
            'deep_z': [-5.0], # -1.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.43], # Launiainen et al. 2019
                    'ThetaR': [0.05], # Launiainen et al. 2019
                    'alpha': [0.024], # Launiainen et al. 2019
                    'n': [1.2]}, # # Launiainen et al. 2021
            'deep_ksat': [1E-05],
                },         
        'Fluvial_sedimend_sand':{ # Fine
            'deep_id': 11,
            'deep_z': [-0.8, -3.0, -14.0], # -32.0
            'pF': {  # vanGenuchten water retention parameters
                    'ThetaS': [0.41,0.41,0.41],
                    'ThetaR': [0.05,0.05,0.05],
                    'alpha': [0.024,0.024,0.024],
                    'n': [1.2,1.2,1.2]},
            'deep_ksat': [3E-5, 1E-8, 1E-6],
                },           
        }
    return deepp


def root_properties():
    """
    swedish_soilmap_root
    """
    rootp = {
            'Postglacial_sand': # C4
                {
                 'root_id': 1,
                 'root_poros': 0.53,
                 'root_fc': 0.22, #0.24,
                 'root_wp': 0.06, #0.08,
                 'root_alpha': 0.37,
                 'root_beta': 4.0,                 
                 'root_n': 1.24,
                 'root_wr': 0.0,
                 'root_ksat': 5e-5,
                 },
            'Glaciofluvial_sediment': #  CoarseText
                {
                 'root_id': 2,
                 'root_poros': 0.41,
                 'root_fc': 0.21,
                 'root_wp': 0.10, 
                 'root_alpha': 0.024, # UNIT????
                 'root_beta': 3.1,                    
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 'root_ksat': 1e-4, 
                 },        
            'Peat':
                {
                 'root_id': 3,
                 'root_poros': 0.89,
                 'root_fc': 0.54, # Leppä et al. Spaghnum -10 kPa (-1m)
                 'root_wp': 0.22, # Leppä et al. Spaghnum -1500 kPa (-150m)
                 'root_alpha': 0.4, # per kPa
                 'root_beta': 4.0,                 
                 'root_n': 1.46,
                 'root_wr': 0.178,
                 'root_ksat': 1e-5,
                 },
            'Postglacial_sand_gravel': #  CoarseText
                {'root_id': 4,
                 'root_poros': 0.41,
                 'root_fc': 0.21,
                 'root_wp': 0.10, 
                 'root_alpha': 0.024, # UNIT????
                 'root_beta': 3.1,                    
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 'root_ksat': 1e-4, 
                 },
            'Clay_silt': #  C3
                {'root_id': 5,
                 'root_poros': 0.55,
                 'root_fc': 0.26,
                 'root_wp': 0.09,
                 'root_alpha': 0.448,
                 'root_beta': 4.0,                    
                 'root_n': 1.20,
                 'root_wr': 0.0,
                 'root_ksat': 1e-5,
                },
            'Washed_sediment_gravel_boulders': # CoarseText
                {
                 'root_id': 6,
                 'root_poros': 0.41,
                 'root_fc': 0.21,
                 'root_wp': 0.10, 
                 'root_alpha': 0.024, # UNIT????
                 'root_beta': 3.1,                    
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 'root_ksat': 1e-4, 
                 },
            'Water': # MediumText
                {
                 'root_id': 7,
                 'root_poros': 0.43,
                 'root_fc': 0.33,
                 'root_wp': 0.13,
                 'root_ksat': 1e-05,
                 'root_beta': 4.7,
                 'root_alpha': 0.024, # UNIT?
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 },        
            'Moraine': # C5
                {
                 'root_id': 8,
                 'root_poros': 0.48,
                 'root_fc': 0.14,
                 'root_wp': 0.04,
                 'root_alpha': 0.38,
                 'root_beta': 4.0,                 
                 'root_n': 1.42,
                 'root_wr': 0.03,
                 'root_ksat': 1e-4,
                 },
            'Fill':  # MediumText
                {
                 'root_id': 9,
                 'root_poros': 0.43,
                 'root_fc': 0.33,
                 'root_wp': 0.13,
                 'root_ksat': 1e-05,
                 'root_beta': 4.7,
                 'root_alpha': 0.024, # UNIT?
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 },
            'Bedrock': # 
                {
                 'root_id': 10,
                 'root_poros': 0.43,
                 'root_fc': 0.33,
                 'root_wp': 0.13,
                 'root_ksat': 1e-05,
                 'root_beta': 4.7,
                 'root_alpha': 0.024, # UNIT?
                 'root_n': 1.2,
                 'root_wr': 0.05,
                 #'root_depth': 0.10, # !! CHECK IF WORKS STRAIGHT
                 },
            'Fluvial_sedimend_sand': # C1
                {
                 'root_id': 11,
                 'root_poros': 0.58,
                 'root_fc': 0.30, # 0.34,
                 'root_wp': 0.13, # 0.11,
                 'root_alpha': 0.406,
                 'root_beta': 4.0,                    
                 'root_n': 1.17,
                 'root_wr': 0.0,
                 'root_ksat': 1e-6, #1e-06,
                 },
            }

    return rootp


def org_properties():
    """
    swedish_soilmap
    """
    orgp = {
        'Postglacial_sand':{
            'org_id': 1,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Glaciofluvial_sediment':{
            'org_id': 2,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Peat':{
            'org_id': 3,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.65,
            'org_rw': 0.3,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Postglacial_sand_gravel':{
            'org_id': 4,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Clay_silt':{
            'org_id': 5,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Washed_sediment_gravel_boulders':{
            'org_id': 6,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Water':{
            'org_id': 7,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },        
        'Moraine':{
            'org_id': 8,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Fill':{
            'org_id': 9,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Bedrock':{
            'org_id': 10,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        'Fluvial_sedimend_sand':{
            'org_id': 11,
            'org_depth': 0.05,
            'org_poros': 0.9,
            'org_fc': 0.3,
            'org_rw': 0.2,
            'org_ksat': 1E-03,
            'org_beta': 6.0
            },
        }
    return orgp


