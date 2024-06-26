# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:18:57 2016

@author: slauniai & khaahti & jpnousu

"""

import numpy as np
import pandas as pd
from canopygrid import CanopyGrid
from bucketgrid import BucketGrid
from topmodel import Topmodel_Homogenous as Topmodel
from soilprofile2D import SoilGrid_2Dflow as SoilGrid

eps = np.finfo(float).eps  # machine epsilon

""" ************** SpaFHy v1.0 ************************************

Simple spatial hydrology and catchment water balance model.

CONSISTS OF FOUR CLASSES, defined in separate modules:
    CanopyGrid - vegetation and snowpack water storages and flows
    BucketGrid - topsoil bucket model (root zone / topsoil water storage)
    SoilGrid - 2D groundwater model (deepsoil water storage)
    Topmodel - integration to catchment scale using Topmodel -concept
HELPER FUNCTIONS:
    parameters - parameter definition file
    iotools - utility functions for data input & output

MAIN PROGRAM:
    spafhy is main program, call it as

NEEDS 2D gis rasters in ascii-grid format

CanopyGrid & BucketGrid can be initialized from gis-data or set to be spatially constant

ToDo:
    CanopyGrid:
        -include topographic shading to radiation received at canopy top
        -radiation-based snowmelt coefficient
        -add simple GPP-model; 2-step Farquhar or LUE-based approach
    BucketGrid:
        -make soil hydrologic properties more realistic e.g. using pedotransfer functions
        -kasvupaikkatyyppi (multi-NFI) --> soil properties
        -add soil frost model, simplest would be Stefan equation with coefficients modified based on snow insulation
          --> we need snow density algorithm: SWE <-----> depth
    Topmodel:
        -think of definging 'relative m & to grids' (soil-type & elevation-dependent?) and calibrate 'catchment averages'
        -topmodel gives 'saturated zone storage deficit in [m]'. This can be converted to gwl proxy (?) if:
        local water retention characteristics are known & hydrostatic equilibrium assumes.
        Look which water retention model was analytically integrable (Campbell, brooks-corey?)

    Graphics and analysis of results:
        -make ready functions


(C) Samuli Launiainen, Kersti Leppä, Jari-Pekka Nousu

VERSION 05.10.2018 / equations correspond to GMDD paper

"""



"""
******************************************************************************
            ----------- SpaFHy model class --------
******************************************************************************
"""


class SpaFHy():
    """
    SpaFHy model class
    """
    def __init__(self, pgen, pcpy, pbu, pds, ptop, flatten=True):

        self.dt = pgen['dt']  # s
        self.simtype = pgen['simtype']

        if self.simtype == '2D':
            flatten = False

        #cmask = pcpy['cmask']
        
        """
        flatten=True omits cells outside catchment
        
        if flatten:
            ix = np.where(np.isfinite(cmask))
    
            for key in pcpy['state']:
                pcpy['state'][key] = pcpy['state'][key][ix].copy()
                        
            for key in pbu['state']:
                pbu['state'][key] = pbu['state'][key][ix].copy()
                
            self.ix = ix  # indices to locate back to 2d grid
        """

        """--- initialize BucketGrid ---"""
        self.bu = BucketGrid(pbu, pgen['org_drain'])

        """--- initialize CanopyGrid ---"""
        self.cpy = CanopyGrid(pcpy, pcpy['state'], dist_rad_file=pgen['spatial_radiation_file'])

        if self.simtype == '2D':
            """--- initialize SoilGrid ---"""
            self.ds = SoilGrid(pds)
        elif self.simtype == 'TOP':
            """--- initialize Topmodel ---"""
            self.top = Topmodel(ptop)

        self.timestep = 1

    def run_timestep(self, forc):
        """
        Runs SpaFHy for one timestep starting from current state
        Args:
            forc - dictionary or pd.DataFrame containing forcing values for the timestep
            ncf - netCDF -file handle, for outputs
            flx - returns flux and state grids to caller as dict
            ave_flx - returns averaged fluxes and states to caller as dict
        Returns:
            dict of results from canopy and soil models
        """
        doy = forc['doy'].values
        ta = forc['air_temperature'].values
        vpd = forc['vapor_pressure_deficit'].values + eps
        rg = forc['global_radiation'].values
        par = forc['par'].values + eps
        prec = forc['precipitation'].values
        co2 = forc['CO2'].values
        u = forc['wind_speed'].values + eps

        print('Running timestep: ', self.timestep)
        self.timestep += 1

        if self.simtype == '2D':
            # run deep soil (2D) water balance
            RR = self.bu.drain
            deep_results = self.ds.run_timestep(
                dt=self.dt / 86400.,
                RR=RR)

            # run CanopyGrid
            canopy_results = self.cpy.run_timestep(
                    doy, self.dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                    beta=self.bu.Ree, Rew=self.bu.Rew, P=101300.0)

            # run BucketGrid
            bucket_results = self.bu.run_timestep(
                dt=self.dt,
                rr=1e-3*canopy_results['potential_infiltration'],
                tr=1e-3*canopy_results['transpiration'],
                evap=1e-3*canopy_results['forestfloor_evaporation'],
                retflow=self.ds.qr,
                airv_deep=self.ds.airv_deep) 

            return deep_results, canopy_results, bucket_results


        elif self.simtype == 'TOP':
            # run Topmodel water balance
            RR = self.bu._drainage_to_gw * self.top.CellArea / self.top.CatchmentArea
            #bu_airv = self.bu.Wair_root

            top_results = self.top.run_timestep(R=RR)

            # run CanopyGrid
            canopy_results = self.cpy.run_timestep(
                    doy, self.dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                    beta=self.bu.Ree, Rew=self.bu.Rew, P=101300.0)

            # run BucketGrid
            bucket_results = self.bu.run_timestep(
                dt=self.dt,
                rr=1e-3*canopy_results['potential_infiltration'],
                tr=1e-3*canopy_results['transpiration'],
                evap=1e-3*canopy_results['forestfloor_evaporation'],
                retflow=self.top.qr)

            return top_results, canopy_results, bucket_results

        elif self.simtype == '1D':

            # run CanopyGrid
            canopy_results = self.cpy.run_timestep(
                    doy, self.dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                    beta=self.bu.Ree, Rew=self.bu.Rew, P=101300.0)

            # run BucketGrid
            bucket_results = self.bu.run_timestep(
                dt=self.dt,
                rr=1e-3*canopy_results['potential_infiltration'],
                tr=1e-3*canopy_results['transpiration'],
                evap=1e-3*canopy_results['forestfloor_evaporation'],
                )

            return canopy_results, bucket_results