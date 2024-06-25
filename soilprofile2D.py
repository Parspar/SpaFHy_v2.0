# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:59:54 2019
@author: alauren

Modified by khaahti & jpnousu
"""

import numpy as np
from scipy.stats.mstats import gmean
from scipy.interpolate import interp1d
from scipy.sparse import diags, linalg
import matplotlib.pyplot as plt
eps = np.finfo(float).eps

apply_vectorized = np.vectorize(lambda f, x: f(x))

class SoilGrid_2Dflow(object):
    """
    2D soil water flow model based on Ari Lauren SUSI2D
    Simulates deep soil water storage, and drainage to ditches.
    """
    def __init__(self, spara):
        """
        Initializes SoilProfile2D:
        Args:
            spara (dict):
                'elevation': elevation [m]
                'streams': ditch water level [m], < 0 for ditches otherwise 0
                'dxy': cell horizontal length
                # scipy interpolation functions describing soil behavior
                'wtso_to_gwl'
                'gwl_to_wsto'
                'gwl_to_Tr'
                'gwl_to_C'
                'gwl_to_rootmoist'
                # initial states
                'ground_water_level': groundwater depth [m]
        """

        """ deep soil """

        # soil/peat type
        self.soiltype = spara['soiltype']

        # interpolated functions for soil column ground water depth vs. water storage, transmissivity etc.
        self.wsto_to_gwl = spara['wtso_to_gwl']
        self.gwl_to_wsto = spara['gwl_to_wsto']
        self.gwl_to_Tr = spara['gwl_to_Tr']
        self.gwl_to_C = spara['gwl_to_C']
        self.gwl_to_rootmoist = spara['gwl_to_rootmoist']

        # initial h (= gwl) and boundaries [m]
        self.ditch_h = spara['streams']
        self.h = spara['ground_water_level']
        # soil surface elevation and hydraulic head [m]
        self.ele = spara['elevation']
        self.H = self.ele + self.h

        # replace nans (values outside catchment area)
        self.H[np.isnan(self.H)] = -999
        #self.h[np.isnan(self.h)] = -999
        # water storage [m]
        self.Wsto_deep_max = np.full_like(self.h, 0.0)  # storage of fully saturated profile
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_deep_max[self.soiltype == key] = value(0.0)
        self.Wsto_deep = np.full_like(self.h, 0.0)  # storage corresponding to h
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_deep[self.soiltype == key] = value(self.h[self.soiltype == key])

        # air volume and returnflow
        self.airv_deep = np.maximum(0.0, self.Wsto_deep_max - self.Wsto_deep)
        self.qr = np.full_like(self.h, 0.0) # 
        
        # rootzone moisture [m3 m-3], parameters related to transpiration limit during dry conditions
        self.deepmoist = np.full_like(self.h, 0.0)
        self.deepmoist[np.isnan(self.h)] = np.nan
        self.deep_fc0 = np.full_like(self.h, 0.0)
        self.deep_fc1 = np.full_like(self.h, 0.0)
        self.deep_wp = np.full_like(self.h, 0.0)
        for key, value in self.gwl_to_rootmoist.items():
            self.deepmoist[self.soiltype == key] = value(self.h[self.soiltype == key])
            self.deep_fc0[self.soiltype == key] = value(-0.7 - 0.1)
            self.deep_fc1[self.soiltype == key] = value(-1.2 - 0.1)
            self.deep_wp[self.soiltype == key] = value(-150.0 - 0.1)

        """ parameters for 2D solution """
        # parameters for solving
        # 0.5 seems to work better when gwl is close to impermeable bottom
        # (probably because transmissivity does not switch between 0. and > 0 as much)
        self.implic = 0.5  # solving method: 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson

        # grid
        self.rows = np.shape(self.h)[0]
        self.cols = np.shape(self.h)[1]
        self.n = self.rows * self.cols  # length of flattened array
        self.dxy = spara['dxy']  # horizontal distance between nodes dx=dy [m]

        # create arrays needed in computation only once
        # previous time step neighboring hydraylic head H (West, East, North, South)
        self.HW = np.zeros((self.rows,self.cols))
        self.HE = np.zeros((self.rows,self.cols))
        self.HN = np.zeros((self.rows,self.cols))
        self.HS = np.zeros((self.rows,self.cols))
        # previous time step transmissivities (West, East, North, South)
        self.TrW0 = np.zeros((self.rows,self.cols))
        self.TrE0 = np.zeros((self.rows,self.cols))
        self.TrN0 = np.zeros((self.rows,self.cols))
        self.TrS0 = np.zeros((self.rows,self.cols))
        # current time step transmissivities (West, East, North, South)
        self.TrW1 = np.zeros((self.rows,self.cols))
        self.TrE1 = np.zeros((self.rows,self.cols))
        self.TrN1 = np.zeros((self.rows,self.cols))
        self.TrS1 = np.zeros((self.rows,self.cols))
        # computation matrix
        # self.A = np.zeros((self.n,self.n))

        self.CC = np.ones((self.rows,self.cols))
        self.Tr0 = np.zeros((self.rows,self.cols))
        self.Tr1 = np.zeros((self.rows,self.cols))
        self.Wtso1_deep = np.zeros((self.rows,self.cols))
        self.tmstep = 0
        self.conv99 = 99
        #self.totit = 0

    def rolling_window(self, a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

    def run_timestep(self, dt=1.0, RR=0.0):

        r""" Solves soil water storage in column assuming hydrostatic equilibrium.

        Args:
            dt (float): solution timestep [days!!!!]
            rr (float/array): potential infiltration [m]
            tr (float/array): transpiration from root zone [m]
            evap (float/array): evaporation from top layer [m]
        Returns:
            results (dict)

        """

        #***********REMIND: map of array*******************
        """
        #2D array: indices i row, j col
        #Flattened array: n from 0 to rows*cols: n=i*cols+j
        #West element: n=i*cols-1
        #East element: n=i*cols+1
        #North element: n=i*cols+j-cols
        #South element: n=i*cols-j+cols
        """
        self.tmstep += 1

        # for computing mass balance later, RR: drainage from bucketgrid
        S = RR
        S[np.isnan(S)] = 0.0

        state0 = self.Wsto_deep + S # Wsto_deep ????

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:-1]
        self.HE[:,:-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:-1,:]
        self.HS[:-1,:] = self.H[1:,:]
        
        # ravel 2D arrays
        HW = np.ravel(self.HW)
        HE = np.ravel(self.HE)
        HN = np.ravel(self.HN)
        HS = np.ravel(self.HS)
        H = np.ravel(self.H)
        Wsto_deep = np.ravel(self.Wsto_deep)
        ditch_h = np.ravel(self.ditch_h)
        ele = np.ravel(self.ele)

        # Ditches
        # calculate mean H of neighboring nodes to find out whether ditch is active (constant head)
        # from previous timestep to avoid switching boundary/no-boundary during iteration
        H_neighbours = ditch_h.copy()
        for k in np.where(ditch_h < -eps)[0]:
            H_ave = 0
            n_neigh = 0
            if k%self.cols != 0 and ditch_h[k-1] > -eps: # west non-ditch neighbor
                    H_ave += H[k-1]
                    n_neigh += 1
            if (k+1)%self.cols != 0 and ditch_h[k+1] > -eps: # east non-ditch neighbor
                    H_ave += H[k+1]
                    n_neigh += 1
            if k-self.cols >= 0 and  ditch_h[k-self.cols] > -eps: # north non-ditch neighbor
                    H_ave += H[k-self.cols]
                    n_neigh += 1
            if k+self.cols < self.n and ditch_h[k+self.cols] > -eps: # sounth non-ditch neighbor
                    H_ave += H[k+self.cols]
                    n_neigh += 1
            if n_neigh > 0:
                H_neighbours[k] = H_ave / n_neigh  # average of neighboring non-ditch nodes
            else:  # corners or nodes surrounded by ditches dont have neighbors, given its ditch depth
                H_neighbours[k] = ele[k] + ditch_h[k] + eps

        H_neighbours_2d = np.reshape(H_neighbours,(self.rows,self.cols))

        # Transmissivity of previous timestep [m2 d-1]
        # for ditch nodes that are active, transmissivity calculated based on mean H of
        # neighboring nodes, not ditch depth which would restrict tranmissivity too much
        # Whole profile depth still considered, but I think that is the usual way..
        H_for_Tr = np.where((self.ditch_h < -eps) & (H_neighbours_2d > self.ele + self.ditch_h),
                            H_neighbours_2d, self.H)
        for key, value in self.gwl_to_Tr.items():
            self.Tr0[self.soiltype == key] = value(H_for_Tr[self.soiltype == key] - self.ele[self.soiltype == key])
        # transmissivity at all four sides of the element is computed as geometric mean of surrounding element transimissivities
        TrTmpEW = gmean(self.rolling_window(self.Tr0, 2), -1)
        TrTmpNS = np.transpose(gmean(self.rolling_window(np.transpose(self.Tr0), 2), -1))
        self.TrW0[:,1:] = TrTmpEW
        self.TrE0[:,:-1] = TrTmpEW
        self.TrN0[1:,:] = TrTmpNS
        self.TrS0[:-1,:] = TrTmpNS
        del TrTmpEW, TrTmpNS

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:-1]
        self.HE[:,:-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:-1,:]
        self.HS[:-1,:] = self.H[1:,:]

        # ravel 2D arrays
        # to avoid reshaping, save in other variable
        TrW0 = np.ravel(self.TrW0)
        TrE0 = np.ravel(self.TrE0)
        TrN0 = np.ravel(self.TrN0)
        TrS0 = np.ravel(self.TrS0)

        # hydraulic heads, new iteration and old iteration
        Htmp = self.H.copy()
        Htmp1 = self.H.copy()

        # convergence criteria
        crit = 1e-3  # loosened this criteria from 1e-4, seems mass balance error remains resonable
        maxiter = 100

        for it in range(maxiter):

            # transmissivity [m2 d-1] to neighbouring cells with HTmp1
            # for ditch nodes that are active, transmissivity calculated based on mean H of
            # neighboring nodes, not ditch depth which would restrict tranmissivity too much
            H_for_Tr = np.where((self.ditch_h < -eps) & (H_neighbours_2d > self.ele + self.ditch_h),
                                H_neighbours_2d, Htmp)
            for key, value in self.gwl_to_Tr.items():
                self.Tr1[self.soiltype == key] = value(H_for_Tr[self.soiltype == key] - self.ele[self.soiltype == key])
            TrTmpEW = gmean(self.rolling_window(self.Tr1, 2),-1)
            TrTmpNS = np.transpose(gmean(self.rolling_window(np.transpose(self.Tr1), 2),-1))
            self.TrW1[:,1:] = TrTmpEW
            self.TrE1[:,:-1] = TrTmpEW
            self.TrN1[1:,:] = TrTmpNS
            self.TrS1[:-1,:]=TrTmpNS
            del TrTmpEW, TrTmpNS
            # ravel 2D arrays
            TrW1 = np.ravel(self.TrW1); TrE1= np.ravel(self.TrE1)
            TrN1 = np.ravel(self.TrN1); TrS1 = np.ravel(self.TrS1)

            # differential water capacity dSto/dh
            # CCtmp = self.CC.copy()
            for key, value in self.gwl_to_C.items():
                self.CC[self.soiltype == key] = value(Htmp[self.soiltype == key] - self.ele[self.soiltype == key])
            alfa = np.ravel(self.CC * self.dxy**2 / dt)
            # alfa = np.ravel((0.5*self.CC + 0.5*CCtmp) * self.dxy**2 / dt)

            # water storage
            for key, value in self.gwl_to_wsto.items():
                self.Wtso1_deep[self.soiltype == key] = value(Htmp[self.soiltype == key] - self.ele[self.soiltype == key])

            # Setup of diagonal sparse matrix
            a_d = self.implic * (TrW1 + TrE1 + TrN1 + TrS1) + alfa  # Diagonal
            a_w = -self.implic * TrW1[1:]  # West element
            a_e = -self.implic * TrE1[:-1]  # East element
            a_n = -self.implic * TrN1[self.cols:]  # North element
            a_s = -self.implic * TrS1[:self.n-self.cols]  # South element

            # Knowns: Right hand side of the eq
            Htmp = np.ravel(Htmp)
            hs = (np.ravel(S) * dt * self.dxy**2 + alfa * Htmp
                   - np.ravel(self.Wtso1_deep) * self.dxy**2 / dt + Wsto_deep * self.dxy**2 / dt
                  + (1.-self.implic) * (TrN0*HN) + (1.-self.implic) * (TrW0*HW)
                  - (1.-self.implic) * (TrN0 + TrW0 + TrE0 + TrS0) * H
                  + (1.-self.implic) * (TrE0*HE) + (1.-self.implic) * (TrS0*HS))

            # Ditches
            for k in np.where(ditch_h < -eps)[0]:
                # update a_x and hs when ditch as constant head
                if H_neighbours[k] > ele[k] + ditch_h[k]: # average of neighboring H above ditch depth
                    hs[k] = ele[k] + ditch_h[k]
                    a_d[k] = 1
                    if k%self.cols != 0:  # west node
                        a_w[k-1] = 0
                    if (k+1)%self.cols != 0:  # east node
                        a_e[k] = 0
                    if k-self.cols >= 0:  # north node
                        a_n[k-self.cols] = 0
                    if k+self.cols < self.n:  # south node
                        a_s[k] = 0

            A = diags([a_d, a_w, a_e, a_n, a_s], [0, -1, 1, -self.cols, self.cols],format='csc')

            # Solve: A*Htmp1 = hs
            Htmp1 = linalg.spsolve(A,hs)

            # testing, limit change
            Htmp1 = np.where(np.abs(Htmp1-Htmp)> 0.5, Htmp + 0.5*np.sign(Htmp1-Htmp), Htmp1)

            conv1 = np.max(np.abs(Htmp1 - Htmp))

            max_index = np.unravel_index(np.argmax(np.abs(Htmp1 - Htmp)),(self.rows,self.cols))

            # especially near profile bottom, solution oscillates so added these steps to avoid that
            if it > 40:
                Htmp = 0.25*Htmp1+0.75*Htmp
            elif it > 20:
                Htmp = 0.5*Htmp1+0.5*Htmp
            else:
                Htmp = Htmp1.copy()

            Htmp = np.reshape(Htmp,(self.rows,self.cols))

            # print to get sense what's happening when problems in convergence
            if it > 90:
                print('\t', it, conv1, max_index, self.ditch_h[max_index],
                      Htmp[max_index]-self.ele[max_index])

            if conv1 < crit:
                break
            # end of iteration loop
        if it == 99:
            self.conv99 +=1
        #self.totit += it
        print('Timestep:', self.tmstep, 'iterations', it, conv1, Htmp[max_index]-self.ele[max_index])
        
        # lateral flow is calculated in two parts: one depending on previous time step
        # and other on current time step (lateral flowsee 2/2). Their weighting depends
        # on self. implic
        # lateral flow 1/2
        lateral_flow = (self.implic*(self.TrW1*(self.H - self.HW)
                        + self.TrE1*(self.H - self.HE)
                        + self.TrN1*(self.H - self.HN)
                        + self.TrS1*(self.H - self.HS)))/ self.dxy**2

        """ update state """
        # soil profile
        self.H = Htmp.copy()      
        self.h = self.H - self.ele
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_deep[self.soiltype == key] = value(self.h[self.soiltype == key])
            
        
        # lateral flow 2/2
        lateral_flow += ((1-self.implic)*(self.TrW0*(self.H - self.HW)
                        + self.TrE0*(self.H - self.HE)
                        + self.TrN0*(self.H - self.HN)
                        + self.TrS0*(self.H - self.HS)))/ self.dxy**2


        """ new update state """
        # Let's limit head to 0 and assign rest as return flow to bucketgrid
        Wsto_before_qr = self.Wsto_deep.copy()
        self.h = np.minimum(0.0, self.h)
        self.H = self.h + self.ele
        self.H[np.isnan(self.H)] = -999

        # Updating the storage according to new head
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_deep[self.soiltype == key] = value(self.H[self.soiltype == key] - self.ele[self.soiltype == key])
        
        # The difference is the return flow to bucket grid
        self.qr = np.maximum(0.0, Wsto_before_qr - self.Wsto_deep)
        
        ##################################################
        for key, value in self.gwl_to_rootmoist.items():
            self.deepmoist[self.soiltype == key] = value(self.h[self.soiltype == key])
            
        # ditches are described as constant heads so the netflow to ditches can
        # be calculated from their mass balance
        netflow_to_ditch = (state0  - self.Wsto_deep - lateral_flow - self.qr)
        netflow_to_ditch = np.where(self.ditch_h < -eps, netflow_to_ditch, 0.0)

        # air volume
        self.airv_deep = np.maximum(0.0, self.Wsto_deep_max - self.Wsto_deep)
        
        # mass balance error [m]
        mbe = (state0  - self.Wsto_deep - lateral_flow - self.qr)

        mbe = np.where(self.ditch_h < -eps, 0.0, mbe)
        
        results = {
                'ground_water_level': self.h,  # [m]
                'lateral_netflow': -lateral_flow * 1e3,  # [mm d-1]
                'netflow_to_ditch': netflow_to_ditch * 1e3,  # [mm d-1]
                'water_closure': mbe * 1e3,  # [mm d-1]
                'moisture_deep': self.deepmoist,  # [m3 m-3]
                'water_storage': (self.Wsto_deep) * 1e3, # [mm]
                }

        return results

def gwl_Wsto(z, pF, Ksat=None, root=False):
    r""" Forms interpolated function for soil column ground water dpeth, < 0 [m], as a
    function of water storage [m] and vice versa + others

    Args:
        pF (dict of arrays):
            'ThetaS' saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'ThetaR' residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'alpha' air entry suction [cm\ :sup:`-1`]
            'n' pore size distribution [-]
        dz (np.arrays): soil conpartment thichness, node in center [m]
    Returns:
        (dict):
            'to_gwl': interpolated function for gwl(Wsto)
            'to_wsto': interpolated function for Wsto(gwl)
            'to_C'
            'to_Tr'
    """

    z = np.array(z) # profile depths
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:] # profile depths into profile thicknesses

    # finer grid for calculating wsto to avoid discontinuity in C (dWsto/dGWL)
    z_fine=np.arange(0,min(z),-0.01)-0.01
    dz_fine = z_fine*0.0 + 0.01
    z_mid_fine = dz_fine / 2 - np.cumsum(dz_fine)

    ix = np.zeros(len(z_fine))
    for depth in z:
        ix += np.where(z_fine < depth, 1, 0)

    pF_fine={}
    for key in pF.keys():
        pp = np.array([pF[key][int(ix[i])] for i in range(len(z_fine))])
        pF_fine.update({key: pp})

    # --------- connection between gwl and Wsto, Tr, C------------
    gwl = np.arange(1.0, -10., -1e-2)
    # solve water storage corresponding to gwls
    Wsto_deep = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine)
            + max(0.0,g) for g in gwl]  # water storage above ground surface == gwl
    # Wsto = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine) for g in gwl]  # old

    if root:
        Wsto_deep = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine) for g in gwl]
        Wsto_deep = Wsto_deep/sum(dz)
        GwlToWsto = interp1d(np.array(gwl), np.array(Wsto_deep), fill_value='extrapolate')
        return {'to_rootmoist': GwlToWsto}

    # solve transmissivity corresponding to gwls
    Tr = [transmissivity(dz, Ksat, g) * 86400. for g in gwl]  # [m2 d-1]

    # interpolate functions
    WstoToGwl = interp1d(np.array(Wsto_deep), np.array(gwl), fill_value='extrapolate')
    GwlToWsto = interp1d(np.array(gwl), np.array(Wsto_deep), fill_value='extrapolate')
    GwlToC = interp1d(np.array(gwl), np.array(np.gradient(Wsto_deep)/np.gradient(gwl)), fill_value='extrapolate')
    GwlToTr = interp1d(np.array(gwl), np.array(Tr), fill_value='extrapolate')

    #plt.figure(1)
    #plt.plot(np.array(gwl), np.array(np.gradient(Wsto_deep/np.gradient(gwl))))
    #plt.figure(2)
    #plt.plot(np.array(gwl), np.array(Tr))
    #plt.figure(3)
    #plt.plot(np.array(gwl), np.array(Wsto_deep))

    return {'to_gwl': WstoToGwl, 'to_wsto': GwlToWsto, 'to_C': GwlToC, 'to_Tr': GwlToTr}

def h_to_cellmoist(pF, h, dz):
    r""" Cell moisture based on vanGenuchten-Mualem soil water retention model.
    Partly saturated cells calculated as thickness weigthed average of
    saturated and unsaturated parts.

    Args:
        pF (dict):
            'ThetaS' (array): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'ThetaR' (array): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
            'alpha' (array): air entry suction [cm\ :sup:`-1`]
            'n' (array): pore size distribution [-]
        h (array): pressure head [m]
        dz (array): soil conpartment thichness, node in center [m]
    Returns:
        theta (array): volumetric water content of cell [m\ :sup:`3` m\ :sup:`-3`\ ]

    Kersti Haahti, Luke 8/1/2018
    """

    # water retention parameters
    Ts = np.array(pF['ThetaS'])
    Tr = np.array(pF['ThetaR'])
    alfa = np.array(pF['alpha'])
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)

    # moisture based on cell center head
    x = np.minimum(h, 0)
    theta = Tr + (Ts - Tr) / (1 + abs(alfa * 100 * x)**n)**m

    # correct moisture of partly saturated cells
    ix = np.where(abs(h) < dz/2)
    if len(Ts) == 1:
        ixx = 0
    else:
        ixx = ix
    # moisture of unsaturated part
    x[ix] = -(dz[ix]/2 - h[ix]) / 2
    theta[ix] = Tr[ixx] + (Ts[ixx] - Tr[ixx]) / (1 + abs(alfa[ixx] * 100 * x[ix])**n[ixx])**m[ixx]
    # total moisture as weighted average
    theta[ix] = (theta[ix] * (dz[ix]/2 - h[ix]) + Ts[ixx] * (dz[ix]/2 + h[ix])) / (dz[ix])

    return theta

def transmissivity(dz, Ksat, gwl):
    r""" Transmissivity of saturated layer.

    Args:
       dz (array):  soil conpartment thichness, node in center [m]
       Ksat (array): horizontal saturated hydr. cond. [ms-1]
       gwl (float): ground water level below surface, <0 [m]

    Returns:
       Qz_drain (array): drainage from each soil layer [m3 m-3 s-1]
    """
    z = dz / 2 - np.cumsum(dz)
    Tr = 0.0

    ib = sum(dz)

    # depth of saturated layer above impermeable bottom
    # Hdr = min(max(0, gwl + ib), ib)  # old
    Hdr = max(0, gwl + ib)  # not restricted to soil profile -> transmissivity increases when gwl above ground surface level

    """ drainage from saturated layers above ditch base """
    # layers above ditch bottom where drainage is possible
    ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -ib))

    if Hdr > 0:
        # saturated layer thickness [m]
        dz_sat = np.maximum(gwl - (z - dz / 2), 0)
        # transmissivity of layers  [m2 s-1]
        Trans = Ksat * dz_sat

        """ drainage from saturated layers above ditch base """
        # layers above ditch bottom where drainage is possible
        ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -ib))

        if ix.size > 0:
            dz_sat[ix[-1]] = dz_sat[ix[-1]] + (z[ix][-1] - dz[ix][-1] / 2 + ib)
            Trans[ix[-1]] = Ksat[ix[-1]] * dz_sat[ix[-1]]
            Tr = sum(Trans[ix])
    return Tr


def wrc(pF, theta=None, psi=None, draw_pF=False):
    """
    vanGenuchten-Mualem soil water retention model 

    References:
        Schaap and van Genuchten (2005). Vadose Zone 5:27-34
        van Genuchten, (1980). Soil Science Society of America Journal 44:892-898

    Args:
        pF (dict):
            ThetaS (float|array): saturated water content [m3 m-3]
            ThetaR (float|array): residual water content [m3 m-3]
            alpha (float|array): air entry suction [cm-1]
            n (float|array): pore size distribution [-]
        theta (float|array): vol. water content [m3 m-3]
        psi (float|array): water potential [m]
        draw_pF (bool): Draw pF-curve.
    Returns:
        y (float|array): water potential [m] or vol. water content [m3 m-3]. Returns None if only curve is drawn.

    """
    
    EPS = np.finfo(float).eps
    
    Ts = np.array(pF['ThetaS'])
    Tr = np.array(pF['ThetaR'])
    alfa = np.array(pF['alpha'])
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)

    def theta_psi(x):
        # converts water content [m3 m-3] to potential [m]]
        x = np.minimum(x, Ts)
        x = np.maximum(x, Tr)  # checks limits
        s = (Ts - Tr) / ((x - Tr) + EPS)
        Psi = -1e-2 / alfa*(s**(1.0 / m) - 1.0)**(1.0 / n)  # m
        Psi[np.isnan(Psi)] = 0.0
        return Psi

    def psi_theta(x):
        # converts water potential [m] to water content [m3 m-3]
        x = 100*np.minimum(x, 0)  # cm
        Th = Tr + (Ts - Tr) / (1 + abs(alfa*x)**n)**m
        return Th

    # --- convert between theta <-- --> psi
    if (theta != None).any():
        y = theta_psi(theta)  # 'Theta-->Psi'
    elif (psi != None).any():
        y = psi_theta(psi)  # 'Psi-->Theta'

    # draws pf-curve
    if draw_pF:
        Ts = Ts[0]; Tr = Tr[0]; alpha = alfa[0]; n = n[0]  
        xx = -np.logspace(-4, 5, 100)  # cm
        yy = psi_theta(xx)

        #  field capacity and wilting point
        fc = psi_theta(-1.0)
        wp = psi_theta(-150.0)

        fig = plt.figure(99)
        fig.suptitle('vanGenuchten-Mualem WRC', fontsize=16)
        ttext = r'$\theta_s=$' + str(Ts) + r', $\theta_r=$' + str(Tr) +\
                r', $\alpha=$' + str(alfa) + ',n=' + str(n)

        plt.title(ttext, fontsize=14)
        plt.semilogx(-xx, yy, 'g-')
        plt.semilogx(1, fc, 'ro', 150, wp, 'ro')  # fc, wp
        plt.text(1, 1.1*fc, 'FC'), plt.text(150, 1.2*wp, 'WP')
        plt.ylabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
        plt.xlabel('$\psi$ $(m)$', fontsize=14)
        plt.ylim(0.8*Tr, min(1, 1.1*Ts))

        del xx, yy
        y = None
    
    return y
