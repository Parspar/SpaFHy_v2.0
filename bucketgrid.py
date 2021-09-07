# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 10:38:59 2021

@author: janousu
"""

import numpy as np
eps = np.finfo(float).eps

class BucketGrid(object):
    """
    Two-layer soil water bucket model for gridded use in SpaFHy.
    """
    def __init__(self, spara, outputs=False):
        """
        Initializes BucketGrid:
        Args:
            REQUIRED:
            spara - dictionary of soil properties. keys - values np.arrays
                depth [m]
                poros [m3m-3]
                fc [m3m-3]
                wp [m3m-3]
                ksat [ms-1]
                beta [-]
                maxpond [-]
                org_depth [m]
                org_poros [m3m-3]
                org_fw [m3m-3]
                org_rw [m3m-3]
            
                pond_sto - initial pond storage [m]
                org_sat - initial saturation or organic layer [-]
                rootzone_sat - initial saturation of root zone [-]
            OPTIONAL:  
            outputs - True appends output grids to dict stored within object 

        CHANGES:
            05.05.2020 removed typo in watbal mbe computation and added outputs
        """
        
        """ set object properties. All will be 1d or 2d arrays of same shape """
        # above-ground pond storage [m]
        self.MaxPond = spara['maxpond']
        
        # top layer is interception storage, which capacity depends on its depth [m]
        # and field capacity
        self.D_top = spara['org_depth']     # depth, m3 m-3
        self.poros_top = spara['org_poros'] # porosity, m3 m-3
        self.Fc_top = spara['org_fc']       # field capacity m3 m-3
        self.rw_top = spara['org_rw']       # ree parameter m3 m-3
        self.MaxStoTop = self.Fc_top * self.D_top # maximum storage m 

        # root-zone layer properties
        self.D_root = spara['root_depth']             # depth, m
        self.poros_root = spara['root_poros']         # porosity, m3 m-3     
        self.Fc_root = spara['root_fc']               # field capacity, m3 m-3 
        self.Wp_root = spara['root_wp']               # wilting point, m3 m-3
        self.Ksat_root = spara['root_ksat']           # sat. hydr. cond., m s-1
        self.beta_root = spara['root_beta']           # hyd. cond. exponent, -
        # self.soilcode = spara['soilcode']   # soil type integer code        
        self.MaxWatStoRoot = self.D_root*self.poros_root  # maximum soil water storage, m

        """
        set buckets initial state: given as arrays
        """
        self.PondSto = np.minimum(spara['pond_storage'], self.MaxPond)
        
        # toplayer storage and relative conductance for evaporation
        self.WatStoTop = self.MaxStoTop * spara['org_sat']
        self.Wliq_top = self.poros_top *self.WatStoTop / (self.MaxStoTop + eps)
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)
        
        # root zone storage and relative extractable water
        self.WatStoRoot = np.minimum(spara['root_sat']*self.D_root*self.poros_root, self.D_root*self.poros_root)
        
        self.Wliq_root = self.poros_root*self.WatStoRoot / self.MaxWatStoRoot
        self.Wair_root = self.poros_root - self.Wliq_root
        self.Sat_root = self.Wliq_root/self.poros_root
        self.Rew = np.minimum((self.Wliq_root - self.Wp_root) / (self.Fc_root - self.Wp_root + eps), 1.0)
        
        # grid total drainage to ground water [m]
        self._drainage_to_gw = 0.0
        self.drain = np.full_like(self.Wliq_root, 0.0)
        self.drain[np.isnan(self.Wliq_root)] = np.nan
        self.retflow = np.full_like(self.Wliq_root, 0.0)

        '''# create dictionary of empty lists for saving results
        if outputs:
            self.results = {'Infil': [], 'Retflow': [], 'Drain': [], 'Roff': [], 'ET': [],
            'Mbe': [], 'Wliq': [], 'PondSto': [], 'Wliq_top': [], 'Ree': []}'''

    def run_timestep(self, dt=1.0, rr=0.0, tr=0.0, evap=0.0, airv_deep=1000.0, retflow=0.0):
        """
        Computes 2-layer bucket model water balance for one timestep dt
        Top layer is interception storage and contributes only to evap.
        Lower layer is rootzone and contributes only tr and creates drainage.
        Capillary interaction between layers is neglected and connection from bottom up
        is only in case of excess returnflow.
        Pond storage can exist above top layer.
        
        IN:
            dt [s]
            rr = potential infiltration [m]
            tr = transpiration from root zone [m]
            evap = evaporation from top layer [m]
            retflow = return flow from ground water [m]
        OUT: dict with 
            inflow [m] - total inflow to root zone
            roff [m] - surface runoff
            drain [m] - drainage from root zone
            tr [m] - transpiration from root zone
            mbe [m] - mass balance error

        """
        gridshape = np.shape(self.Wliq_root)  # rows, cols
    
        self.retflow = retflow
        if np.shape(self.retflow) != gridshape:
            self.retflow = self.retflow * np.ones(gridshape)
            
        if np.shape(rr) != gridshape:
            rr = rr * np.ones(gridshape)
        
        rr0 = rr.copy()
       
        # add current Pond storage to rr & update storage
        PondSto0 = self.PondSto.copy()
        rr += self.PondSto
        self.PondSto = np.zeros(gridshape)
        
        WatStoRoot0 = self.WatStoRoot.copy()
        WatStoTop0 = self.WatStoTop.copy()
        
        
        #top layer interception & water balance
        interc = np.maximum(0.0, (self.MaxStoTop - self.WatStoTop))\
                    * (1.0 - np.exp(-(rr / (self.MaxStoTop + eps))))
        
        self.WatStoTop = np.maximum(0.0, self.WatStoTop + interc)  
        evap = np.minimum(evap, self.WatStoTop)
        self.WatStoTop -= evap
      
        # infiltration to rootzone
        rr = rr - interc
                
        # ********* compute bottom layer (root zone) water balance ***********

        # transpiration removes water from rootzone
        tr = np.minimum(tr, self.WatStoRoot - eps)
        self.WatStoRoot -= tr
        
        # drainage: at gridcells where retflow > 0, set drain to zero.
        # This delays drying of cells which receive water from topmodel storage
        # ... and removes oscillation of water content at those cells.
        self.drain = np.minimum(self.hydrCond() * dt, np.maximum(0.0, (self.Wliq_root - self.Fc_root))*self.D_root)
        self.drain[self.retflow > 0.0] = 0.0
        #airv_deep = airv_deep * 1e3
        self.drain = np.minimum(self.drain, airv_deep)
        # inflow to root zone: restricted by potential inflow or available pore space
        Qin = self.retflow + rr  # m, pot. inflow
        
        inflow = np.minimum(Qin, self.MaxWatStoRoot - self.WatStoRoot + self.drain)
        
        dSto = (inflow - self.drain)
        self.WatStoRoot = np.minimum(self.MaxWatStoRoot, np.maximum(self.WatStoRoot + dSto, eps))
                
        # if inflow excess after filling rootzone, update first top layer storage
        exfil = Qin - inflow
        to_top_layer = np.minimum(exfil, self.MaxStoTop - self.WatStoTop - eps)
        # self.WatStoTop = self.WatStoTop + to_top_layer
        self.WatStoTop += to_top_layer
        
        # ... and then pond storage ...
        to_pond = np.minimum(exfil - to_top_layer, self.MaxPond - self.PondSto - eps)
        self.PondSto += to_pond
 
        # ... and route remaining to surface runoff
        roff = exfil - to_top_layer - to_pond
        
        # compute diagnostic state variables at root zone:
        self.setState()
        
        # update grid total drainage to ground water [m]
        self._drainage_to_gw = np.nansum(self.drain)
        
        # update self.drain into mm
        self.drain = self.drain * 1e3
        self.retflow = self.retflow * 1e3

        
        # mass balance error [m]
        mbe = (self.WatStoRoot - WatStoRoot0)  + (self.WatStoTop - WatStoTop0) + (self.PondSto - PondSto0) \
            - (rr0 + self.retflow - tr - evap - self.drain - roff)
            
        #print('uniq:', np.unique(self.Wliq_top), 'shape:', self.Wliq_top.shape)
        results = {
                'infiltration': inflow * 1e3,  # [mm d-1]
                'evaporation': evap * 1e3,  # [mm d-1]
                'transpiration': tr * 1e3,  # [mm d-1]
                'drainage': self.drain, #     !!!
                'surface_runoff': roff * 1e3, #  !!!
                'water_closure': mbe * 1e3,  # [mm d-1]
                'moisture_top': self.Wliq_top,  # [m3 m-3]
                'moisture_root': self.Wliq_root,  # [m3 m-3]
                'transpiration_limitation': self.Rew,  # [-] !!!
                'water_storage': (self.WatStoTop + self.WatStoRoot) * 1e3, # [mm]
                'return_flow': self.retflow # [mm d-1]
                }

        return results

        '''
        # append results to lists; use only for testing small grids!
        if hasattr(self, 'results'):
            self.results['Infil'].append(inflow - retflow)   # infiltration through top boundary
            self.results['Retflow'].append(retflow)         # return flow from below 
            self.results['Roff'].append(roff)       # surface runoff
            self.results['Drain'].append(drain)     # drainage
            self.results['ET'].append(tr + evap)            
            self.results['Mbe'].append(mbe)
            self.results['Wliq'].append(self.Wliq)
            self.results['PondSto'].append(self.PondSto)
            self.results['Wliq_top'].append(self.Wliq_top)
            self.results['Ree'].append(self.Ree)
        
        return inflow, roff, drain, tr, evap, mbe
        '''

    def setState(self):
        """ updates state variables"""
        # root zone
        self.Wliq_root = self.poros_root*self.WatStoRoot / self.MaxWatStoRoot
        self.Wair_root = self.poros_root - self.Wliq_root
        self.Sat_root = self.Wliq_root / self.poros_root
        self.Rew = np.maximum(0.0,
              np.minimum((self.Wliq_root - self.Wp_root) / (self.Fc_root - self.Wp_root + eps), 1.0))
        
        # organic top layer; maximum that can be hold is Fc
        self.Wliq_top = self.Fc_top * self.WatStoTop / (self.MaxStoTop + eps)
        self.Ree = self.relative_evaporation()
        
    def hydrCond(self):
        """
        returns hydraulic conductivity [ms-1] based on Campbell -formulation
        """
        k = self.Ksat_root*self.Sat_root**(2*self.beta_root + 3.0)
        return k

    def relative_evaporation(self):
        """
        returns relative evaporation rate from the organic top layer; loosely
        based on Launiainen et al. 2015 Ecol. Mod. Moss-module
        Returns:
            f - [-], array or grid of 
        """
        f = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))
        return f