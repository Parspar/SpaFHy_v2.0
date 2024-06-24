# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 10:38:59 2021

@author: janousu
"""

import numpy as np
eps = np.finfo(float).eps
import matplotlib.pyplot as plt
import sys

class BucketGrid(object):
    """
    Two-layer soil water bucket model for gridded use in SpaFHy.
    """
    def __init__(self, spara, org_drain):
        """
        Initializes BucketGrid:
        Args:
            REQUIRED:
            spara - dictionary of soil properties. keys - values are np.arrays
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

        """

        """ set object properties. All will be 1d or 2d arrays of same shape """
        # organic layer drainage True/False
        self.org_drain = org_drain
        
        # above-ground pond storage [m]
        self.MaxPond = spara['maxpond']

        # top layer is interception storage, which capacity depends on its depth [m]
        # and field capacity/porosity
        self.D_top = spara['org_depth']     # depth, m3 m-3
        self.poros_top = spara['org_poros'] # porosity, m3 m-3
        self.Fc_top = spara['org_fc']       # field capacity m3 m-3
        self.rw_top = spara['org_rw']       # ree parameter m3 m-3
        
        # top layer maximum storage depends on org_drain option:
        if self.org_drain == True: # different maximum for interception and total storages
            self.MaxStoTop = self.poros_top * self.D_top # maximum storage, m
            self.Ksat_top = spara['org_ksat']       # sat. hydr. cond., m s-1
            self.beta_top = spara['org_beta']       # hyd. cond. exponent, -
        elif self.org_drain == False: # same maximum for interception and total storages
            self.MaxStoTop = self.Fc_top * self.D_top # maximum storage, m
            
        # maximum interception storage
        self.MaxStoTopInt = self.Fc_top * self.D_top # maximum storage for interception, m
        
        # root-zone layer is a bucket, receives infiltration and returnflow, outflows
        # are transpiration and drainage
        self.D_root = spara['root_depth']             # depth, m
        self.poros_root = spara['root_poros']         # porosity, m3 m-3
        self.Fc_root = spara['root_fc']               # field capacity, m3 m-3
        self.Wp_root = spara['root_wp']               # wilting point, m3 m-3
        self.Ksat_root = spara['root_ksat']           # sat. hydr. cond., m s-1
        self.beta_root = spara['root_beta']           # hyd. cond. exponent, -
        self.alpha_root = spara['root_alpha']         # kPa-1
        self.n_root = spara['root_n']
        self.wr_root = spara['root_wr']               # m3m-3
        self.MaxStoRoot = self.D_root*self.poros_root # maximum soil water storage, mi
        
        """
        set buckets initial state
        """
        self.PondSto = np.minimum(spara['pond_storage'], self.MaxPond)

        # toplayer storage and relative conductance for evaporation
        self.WatStoTop = self.MaxStoTop * spara['org_sat']
        try:
            self.WatStoTop = spara['top_storage']#, self.MaxStoTop * spara['org_sat'])
        except:
            pass

        self.Wliq_top = (self.MaxStoTop / self.D_top) * self.WatStoTop / (self.MaxStoTop + eps)
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)

        # root zone storage and relative extractable water
        self.WatStoRoot = np.minimum(spara['root_sat']*self.D_root*self.poros_root, self.D_root*self.poros_root)
        try:
            self.WatStoRoot = spara['root_storage']#, spara['root_sat']*self.D_root*self.poros_root)
        except:
            pass
        self.Wliq_root = self.poros_root*self.WatStoRoot / self.MaxStoRoot
        self.Wair_root = self.poros_root - self.Wliq_root
        self.Sat_root = self.Wliq_root/self.poros_root
        self.Sat_top = self.Wliq_top/self.poros_top
        self.Rew = np.minimum((self.Wliq_root - self.Wp_root) / (self.Fc_root - self.Wp_root + eps), 1.0)
        
        # drainage to rootzone
        if self.org_drain == True:
            self.drain_top = np.full_like(self.Wliq_top, 0.0)
            self.drain_top[np.isnan(self.Wliq_top)] = np.nan
        
        # grid total drainage to ground water [m]
        self._drainage_to_gw = 0.0
        self.drain = np.full_like(self.Wliq_root, 0.0)
        self.drain[np.isnan(self.Wliq_root)] = np.nan
        self.retflow = np.full_like(self.Wliq_root, 0.0)

    def run_timestep(self, dt=1.0, rr=0.0, tr=0.0, evap=0.0, airv_deep=1000.0, retflow=0.0):
        """
        Computes 2-layer bucket model water balance for one timestep dt
        Top layer is interception storage: dW/dt = Interception - Evaporation + Recharge from returnflow
        Lower layer is rootzone: DW/dt = Infiltration - Transpiration - Drainae + Recharge from returnflow
        Capillary interaction between layers is neglected and connection from bottom up is only in case of excess returnflow.
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

        # add current Pond storage to rr & update storage, save intial conditions
        PondSto0 = self.PondSto.copy()
        rr += self.PondSto
        self.PondSto = np.zeros(gridshape)

        WatStoRoot0 = self.WatStoRoot.copy()
        WatStoTop0 = self.WatStoTop.copy()

        #top layer interception & water balance
        interc = np.maximum(0.0, (self.MaxStoTopInt - self.WatStoTop))\
                    * (1.0 - np.exp(-(rr / (self.MaxStoTopInt + eps))))
        
        self.WatStoTop = np.maximum(0.0, self.WatStoTop + interc)
        evap = np.minimum(evap, self.WatStoTop)
        self.WatStoTop -= evap
        
        if self.org_drain == True: # drainage according to Campbell 1985
            self.Wliq_top = (self.MaxStoTop / self.D_top) * self.WatStoTop / (self.MaxStoTop + eps)
            self.drain_top = np.minimum(self.hydrCond_top() * dt, np.maximum(0.0, (self.Wliq_top - self.Fc_top))*self.D_top)
            rr = rr - interc + self.drain_top # infiltration/drainage to rootzone
            self.WatStoTop -= self.drain_top
        elif self.org_drain == False: # organic layer as in Launiainen et al., 2019
            rr = rr - interc # infiltration to rootzone

        # ********* compute bottom layer (root zone) water balance ***********

        # transpiration removes water from rootzone
        tr = np.minimum(tr, self.WatStoRoot - eps)
        self.WatStoRoot -= tr

        # drainage: at gridcells where retflow > 0, set drain to zero.
        # This delays drying of cells which receive water from returnflow
        # ... and removes oscillation of water content at those cells.
        self.drain = np.minimum(self.hydrCond() * dt, np.maximum(0.0, (self.Wliq_root - self.Fc_root))*self.D_root)
        self.drain[self.retflow > 0.0] = 0.0
        self.drain = np.minimum(self.drain, airv_deep)

        # inflow to root zone: restricted by potential inflow or available pore space
        Qin = self.retflow + rr  # m, pot. inflow

        inflow = np.minimum(Qin, self.MaxStoRoot - self.WatStoRoot + self.drain)

        dSto = (inflow - self.drain)
        self.WatStoRoot = np.minimum(self.MaxStoRoot, np.maximum(self.WatStoRoot + dSto, eps))

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

        # storage change
        dStorage = (self.WatStoRoot - WatStoRoot0)  + (self.WatStoTop - WatStoTop0) + (self.PondSto - PondSto0)

        # mass balance error [m]
        mbe = (self.WatStoRoot - WatStoRoot0)  + (self.WatStoTop - WatStoTop0) + (self.PondSto - PondSto0) \
            - (rr0 + self.retflow - tr - evap - self.drain - roff)


        results = {
                'potential_infiltration': rr0 * 1e3,  # [mm d-1]
                'evaporation': evap * 1e3,  # [mm d-1]
                'transpiration': tr * 1e3,  # [mm d-1]
                'drainage': self.drain * 1e3, # [mm d-1]    !!!
                'surface_runoff': roff * 1e3, #  !!!
                'water_closure': mbe * 1e3,  # [mm d-1]
                'moisture_top': self.Wliq_top,  # [m3 m-3]
                'moisture_root': self.Wliq_root,  # [m3 m-3]
                'psi_root': self.Psi, # MPa
                'transpiration_limitation': self.Rew,  # [-] !!!
                'water_storage_root': self.WatStoRoot * 1e3, # [mm]
                'water_storage_top': self.WatStoTop * 1e3, # [mm]
                'pond_storage': self.PondSto * 1e3, # [mm]
                'water_storage': (self.WatStoTop + self.WatStoRoot) * 1e3, # [mm]
                'storage_change': dStorage * 1e3, # [mm]
                'return_flow': self.retflow * 1e3 # [mm d-1]
                }

        return results


    def setState(self):
        """ updates state variables"""
        # root zone
        self.Wliq_root = self.poros_root*self.WatStoRoot / self.MaxStoRoot
        self.Wair_root = self.poros_root - self.Wliq_root
        self.Sat_root = self.Wliq_root / self.poros_root
        self.Rew = np.maximum(0.0,
              np.minimum((self.Wliq_root - self.Wp_root) / (self.Fc_root - self.Wp_root + eps), 1.0))

        # organic top layer; maximum that can be hold is Fc or poros
        self.Wliq_top = (self.MaxStoTop / self.D_top) * self.WatStoTop / (self.MaxStoTop + eps) 
        self.Sat_top = self.Wliq_top / self.poros_top
        self.Ree = self.relative_evaporation()
        self.Wliq_top[self.D_top == 0] = np.NaN
        self.Ree[self.D_top == 0] = eps # vie canopyn Efloor'in nollaan (mass-balance consistency)
        self.Psi = self.theta_psi() # MPa

    def theta_psi(self):
        """
        converts vol. water content (m3m-3) to soil water potential (MPa)
        """
        n = self.n_root
        m = 1 - 1 / n
        # converts water content (m3m-3) to potential
        x = np.minimum(self.Wliq_root, self.poros_root)
        x = np.maximum(x, self.wr_root)  # checks limits

        #print('n', self.n_root[160,345])
        #print('poros', self.poros_root[160,345])
        #print('wr', self.wr_root[160,345])
        #print('alpha', self.alpha_root[160,345])
        
        s = (self.poros_root - self.wr_root) / ((x - self.wr_root) + eps)
        
        Psi = -1 / self.alpha_root*(s**(1.0 / m) - 1.0)**(1.0 / n)  # alpha defines the unit (kPa)
        Psi[Psi==np.NaN] = 0.0
        Psi = 1e-3*Psi # kPa to MPa
        #Psi[Psi<-3.0] = -3.0
        
        return Psi
        
    def hydrCond(self):
        """
        returns hydraulic conductivity [ms-1] based on Campbell -formulation
        """
        k = self.Ksat_root*self.Sat_root**(2*self.beta_root + 3.0)
        return k
    
    def hydrCond_top(self):
        """
        returns hydraulic conductivity [ms-1] based on Campbell -formulation
        """
        k = self.Ksat_top*self.Sat_top**(2*self.beta_top + 3.0)
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

