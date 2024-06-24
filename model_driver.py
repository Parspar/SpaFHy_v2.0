# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:52:46 2019

@author: khaaahti
"""

import time
import numpy as np
import pandas as pd
from spafhy import SpaFHy
from iotools import read_FMI_weather, initialize_netcdf, initialize_netcdf_spinup, write_ncf, write_ncf_spinup
import matplotlib.pyplot as plt
import xarray as xr

eps = np.finfo(float).eps

def driver(create_ncf=False, create_spinup=False, output=True, folder=''):
    """
    Model driver: sets up model, runs it and saves results to file (create_ncf==True)
    or return dictionary of results.
    """

    """ set up model """
    running_time = time.time()

    # load and process parameters
    pgen, pcpy, pbu, pds, cmask, ptopmodel, gisinfo = preprocess_parameters(folder)
    
    # initialize SpaFHy
    spa = SpaFHy(pgen, pcpy, pbu, pds, ptopmodel)

    # read forcing data
    forcing = preprocess_forcing(pgen)

    Nsteps = len(forcing['date'])
    Nspin = (pd.to_datetime(pgen['spinup_end']) - pd.to_datetime(pgen['start_date'])).days + 1

    # results dictionary to accumulate simulation results
    # FOR ONE YEAR AT A TIME
    if create_ncf:
        save_interval = min(pgen['save_interval'], Nsteps - Nspin)
        results = _create_results(pgen, cmask, save_interval)
    else:
        save_interval = Nsteps - Nspin
        results = _create_results(pgen, cmask, Nsteps - Nspin)

    Nsaveresults = list(np.arange(Nspin, Nsteps, save_interval)[1:] - 1)

    # save parameters to results
    results = _append_results('parameters', pcpy['state'], results)
    results = _append_results('parameters', pbu, results)    
    results = _append_results('parameters', pcpy['loc'], results)
    if pgen['simtype'] == '2D':
        results = _append_results('parameters', pds, results)
    elif pgen['simtype'] == 'TOP':
        results = _append_results('parameters', ptopmodel, results)


    if create_ncf:
        ncf, outputfile = initialize_netcdf(
                pgen=pgen,
                cmask=cmask,
                filepath=pgen['results_folder'],
                filename=pgen['ncf_file'],
                description=pgen['description'],
                gisinfo=gisinfo)

    if create_spinup:
        ncf_spinup, outputfile_spinup = initialize_netcdf_spinup(
                pgen=pgen,
                cmask=cmask,
                filepath=pgen['results_folder'],
                filename=pgen['ncf_file'],
                description=pgen['description'],
                gisinfo=gisinfo)

    print('*** Running model ***')

    if pgen['simtype'] == '2D':
            print('*** 2D run')
    elif pgen['simtype'] == 'TOP':
            print('*** TOPMODEL run')
    elif pgen['simtype'] == '1D':
            print('*** 1D run')
            
    if pgen['org_drain'] == True:
            print('*** Bucket organic layer drains according to Campbell 1985')
    else:
            print('*** Bucket organic layer as in Launiainen et al., 2019')
        


    interval = 0
    Nsaved = Nspin - 1

    for k in range(0, Nsteps):
        if pgen['simtype'] == '2D':
            deep_results, canopy_results, bucket_results = spa.run_timestep(forcing.isel(date=k))
        elif pgen['simtype'] == 'TOP':
            top_results, canopy_results, bucket_results = spa.run_timestep(forcing.isel(date=k))
        elif pgen['simtype'] == '1D':
            canopy_results, bucket_results = spa.run_timestep(forcing.isel(date=k))

        if (k >= Nspin):  # save results after spinup done
            if pgen['simtype'] == '2D':
                results = _append_results('deep', deep_results, results, k - Nsaved - 1)
            elif pgen['simtype'] == 'TOP':
                results = _append_results('top', top_results, results, k - Nsaved - 1)
            results = _append_results('canopy', canopy_results, results, k - Nsaved - 1)
            results = _append_results('bucket', bucket_results, results, k - Nsaved - 1)

            if k in Nsaveresults and create_ncf:
                interval += 1
                print('*** Writing results to netCDF4-file, subset %.0f/%.0f ***' % (interval, len(Nsaveresults)+1))
                # save forcing to results
                results = _append_results('forcing', forcing[dict(
                        date=slice(Nsaved + 1, k + 1))], results)
                write_ncf(results=results, ncf=ncf, steps=[Nsaved + 1 - Nspin, k + 1 - Nspin])
                Nsaved = k
    if (create_spinup) and (k == Nsteps - 1):
        write_ncf_spinup(results=results, pgen=pgen, ncf_spinup=ncf_spinup)
        print('*** Writing spinup to netCDF4-file')


    if create_ncf:
        interval += 1
        print('*** Writing results to netCDF4-file, subset %.0f/%.0f ***' % (interval, len(Nsaveresults)+1))
        # save forcing to results
        results = _append_results('forcing', forcing[dict(
                date=slice(Nsaved + 1, k + 1))], results)
        write_ncf(results=results, ncf=ncf, steps=[Nsaved + 1 - Nspin, k + 1 - Nspin])
        ncf.close()
        print('--- Running time %.2f seconds ---' % (time.time() - running_time))
        print('--- Results are in file: ' + outputfile + ' ---')
        if output:
            return outputfile
    else:
        print('--- Running time %.2f seconds ---' % (time.time() - running_time))
        if output:
            return results, spa, pcpy, pbu, ptopmodel, cmask

def preprocess_parameters(folder=''):
    """
    Reading gisdata if applicable and preprocesses parameters
    """

    from iotools import read_bu_gisdata, read_ds_gisdata, read_cpy_gisdata, read_forcing_gisdata, read_top_gisdata
    from iotools import preprocess_budata, preprocess_dsdata, preprocess_cpydata, preprocess_topdata
    from parameters import root_properties, org_properties, deep_properties, parameters, root_properties_from_sitetype, ptopmodel

    pgen, pcpy, pbu, pspd = parameters(folder)
    
    spatial_pbu = {}
    if isinstance(pbu['org_id'], float): # org
        spatial_pbu['org'] = False
    if isinstance(pbu['org_id'], str):
        spatial_pbu['org'] = True
    if isinstance(pbu['root_id'], float): # root
        spatial_pbu['root'] = False
    if isinstance(pbu['root_id'], str):
        spatial_pbu['root'] = True

    spatial_pcpy = {}
    if isinstance(pcpy['state']['lai_conif'], float): # LAI_conif
        spatial_pcpy['lai_conif'] = False
    if isinstance(pcpy['state']['lai_conif'], str):
        spatial_pcpy['lai_conif'] = True
    if isinstance(pcpy['state']['lai_decid_max'], float): # LAI_decid
        spatial_pcpy['lai_decid'] = False
    if isinstance(pcpy['state']['lai_decid_max'], str):
        spatial_pcpy['lai_decid'] = True
    if isinstance(pcpy['state']['lai_grass'], float): # LAI_grass
        spatial_pcpy['lai_grass'] = False
    if isinstance(pcpy['state']['lai_grass'], str):
        spatial_pcpy['lai_grass'] = True
    if isinstance(pcpy['state']['lai_shrub'], float): # LAI_shrub
        spatial_pcpy['lai_shrub'] = False
    if isinstance(pcpy['state']['lai_shrub'], str):
        spatial_pcpy['lai_shrub'] = True              
    if isinstance(pcpy['state']['hc'], float): # canopy_height
        spatial_pcpy['hc'] = False
    if isinstance(pcpy['state']['hc'], str):
        spatial_pcpy['hc'] = True    
    if isinstance(pcpy['state']['cf'], float): # canopy_fraction
        spatial_pcpy['cf'] = False
    if isinstance(pcpy['state']['cf'], str):
        spatial_pcpy['cf'] = True

    spatial_pspd = {}
    if isinstance(pspd['deep_id'], float): # deep
        spatial_pspd['deep_id'] = False
    if isinstance(pspd['deep_id'], str):
        spatial_pspd['deep_id'] = True
    if isinstance(pspd['elevation'], float): # elevation
        spatial_pspd['elevation'] = False
    if isinstance(pspd['elevation'], str):
        spatial_pspd['elevation'] = True


    print('spatial_pcpy', spatial_pcpy)
    print('spatial_pbu', spatial_pbu)
    print('spatial_pspd', spatial_pspd)

    orgp = org_properties()
    rootp = root_properties_from_sitetype()
    deepp = deep_properties()
    gisdata = {}
    ptopmodel = ptopmodel()
   
    if pgen['simtype'] == '2D':
        if pgen['mask'] == 'stream': # to make sure streams are not masked with 2D
            pgen['mask'] = None
        elif pgen['mask'] == 'cmask/stream':
            pgen['mask'] = 'cmask'            
        
    if pgen['spatial_soil']:
        gisdata.update(read_bu_gisdata(pgen['gis_folder'], spatial_pbu=spatial_pbu, mask=pgen['mask']))
        
    if pgen['spatial_cpy']:
        gisdata.update(read_cpy_gisdata(pgen['gis_folder'], spatial_pcpy=spatial_pcpy, mask=pgen['mask']))

    if pgen['spatial_forcing']:
        gisdata.update(read_forcing_gisdata(pgen['gis_folder'], mask=pgen['mask']))
        pgen.update({'forcing_id': gisdata['forcing_id']})
        
    if (pgen['spatial_cpy'] == False and
        pgen['spatial_soil'] == False and
        pgen['spatial_forcing'] == False):
        gisdata = {'cmask': np.ones((1,1))}

    budata = preprocess_budata(pbu, orgp, rootp, gisdata, pgen['spatial_soil'])

    cpydata = preprocess_cpydata(pcpy, gisdata, pgen['spatial_cpy'])

    if pgen['simtype'] == 'TOP':
        gisdata.update(read_top_gisdata(pgen['gis_folder'], mask=pgen['mask']))
        ptopmodel = preprocess_topdata(ptopmodel, gisdata, spatial=True)

    if pgen['simtype'] == '2D':
        gisdata.update(read_ds_gisdata(pgen['gis_folder']))
        dsdata = preprocess_dsdata(pspd, deepp, gisdata, pgen['spatial_soil'])
    else:
        dsdata = pspd.copy() # dummy
        
    gisinfo = {}
    gisinfo['xllcorner'] = gisdata['xllcorner']
    gisinfo['yllcorner'] = gisdata['yllcorner']
    gisinfo['dxy'] = gisdata['dxy']

    # overwrites the state variables with the last timestep of spinup file (if given)
    try:
        spinup = xr.open_dataset(pgen['spinup_file'])
        cpydata['w'] = np.array(spinup['canopy_water_storage'][-1]) * 1e-3
        cpydata['swe'] = np.array(spinup['canopy_water_storage'][-1]) * 1e-3
        soildata['top_storage'] = np.array(spinup['bucket_water_storage_top'][-1]) * 1e-3
        soildata['root_storage'] = np.array(spinup['bucket_water_storage_root'][-1]) * 1e-3
        if pgen['simtype'] == '2D':
            soildata['ground_water_level'] = np.array(spinup['soil_ground_water_level'][-1])
        elif pgen['simtype'] == 'TOP':
            ptopmodel['so'] = np.array(spinup['top_saturation_deficit'][-1])

        print('*** State variables assigned from ', pgen['spinup_file'],  '***')
    except:
        print('*** State variables assigned from parameters.py ***')

    return pgen, cpydata, budata, dsdata, gisdata['cmask'], ptopmodel, gisinfo

def preprocess_forcing(pgen):
    """
    Reads forcing file(s) based on indices in pgen['forcing_id']
    Creates xarray dataset of forcing data
    """

    import xarray as xr

    indices = np.array(pgen['forcing_id'],ndmin=2)

    variables = ['doy',
                 'air_temperature',
                 'vapor_pressure_deficit',
                 'global_radiation',
                 'par',
                 'precipitation',
                 'CO2',
                 'relative_humidity',
                 'wind_direction',
                 'wind_speed']

    dims = ['date','i','j']
    dates = pd.date_range(pgen['start_date'], pgen['end_date']).tolist()
    empty_array = np.ones((len(dates),) + np.shape(indices)) * np.nan

    ddict = {var: (dims, empty_array.copy()) for var in variables}

    for index in np.unique(indices):
        if np.isnan(index):
            break
        fp = pgen['forcing_file'].replace('[forcing_id]',str(int(index)))
        df = read_FMI_weather(pgen['start_date'],
                              pgen['end_date'],
                              sourcefile=fp)
        ix = np.where(indices==index)
        for var in variables:
            ddict[var][1][:,ix[0],ix[1]] = np.matmul(
                    df[var].values.reshape(len(dates),1),np.ones((1,len(ix[0]))))
    ds = xr.Dataset(ddict, coords={'date': dates})

    return ds

def _create_results(pgen, cmask, Nsteps):
    """
    Creates results dictionary to accumulate simulation results
    """
    i, j = np.shape(cmask)

    results = {}

    for var in pgen['variables']:

        var_shape = []
        var_name = var[0]

        if var_name.split('_')[0] != 'parameters':
            var_shape.append(Nsteps)
        if (var_name.split('_')[0] != 'forcing' or
            pgen['spatial_forcing'] == True):
            if (var_name.split('_')[0] != 'top'):
                var_shape.append(i)
                var_shape.append(j)
            if (var_name.split('_')[0] == 'top' and var_name.split('_')[1] == 'local'):
                var_shape.append(i)
                var_shape.append(j)

        results[var_name] = np.full(var_shape, np.NAN)

    return results


def _append_results(group, step_results, results, step=None):
    """
    Adds results from each simulation steps to temporary results dictionary
    """

    results_keys = results.keys()

    for key in results_keys:
        var = key.split('_',1)[-1]

        if var in step_results and key.split('_',1)[0] == group:
            if group == 'forcing':
                res = step_results[var].values
            else:
                res = step_results[var]
            if step==None:
                results[key] = res
            else:
                results[key][step] = res
    return results


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', help='parameter folder', type=str)

    args = parser.parse_args()

    outputfile = driver(create_ncf=True, folder=args.folder)

    print(outputfile)