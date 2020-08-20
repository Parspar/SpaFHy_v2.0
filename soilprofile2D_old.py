# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:59:54 2019
@author: alauren

Modified by khaahti
"""

import numpy as np
from scipy.stats.mstats import gmean
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
eps = np.finfo(float).eps

apply_vectorized = np.vectorize(lambda f, x: f(x))

class SoilGrid_2Dflow(object):
    """
    2D soil water flow model based on Ari Lauren SUSI2D
    Simulates moss/organic layer with interception and evaporation,
    soil water storage, and drainage to ditches.
    """
    def __init__(self, spara):
        """
        Initializes SoilProfile2D:
        Args:
            spara (dict):
                'elevation': elevation [m]
                'ditch_depth': ditch depth [m] - ei vielä!
                'dxy': cell horizontal length
                # scipy interpolation functions describing soil behavior
                'wtso_to_gwl'
                'gwl_to_wsto'
                'gwl_to_Tr'
                'gwl_to_C'
                'gwl_to_rootmoist'
                # organic (moss) layer
                'org_depth': depth of organic top layer (m)
                'org_poros': porosity (-)
                'org_fc': field capacity (-)
                'org_rw': critical vol. moisture content (-) for decreasing phase in Ef
                # initial states
                'ground_water_level': groundwater depth [m]
                'org_sat': organic top layer saturation ratio (-)
        """

        """ moss/organic layer """
        # top layer is interception storage, which capacity depends on its depth [m]
        # and field capacity
        self.dz_top = spara['org_depth']  # depth, m3 m-3
        self.poros_top = spara['org_poros']  # porosity, m3 m-3
        self.fc_top = spara['org_fc']  # field capacity m3 m-3
        self.rw_top = spara['org_rw']  # ree parameter m3 m-3
        self.Wsto_top_max = self.fc_top * self.dz_top  # maximum storage m

        # initial state: toplayer storage and relative conductance for evaporation
        self.Wsto_top = self.Wsto_top_max * spara['org_sat']
        self.Wliq_top = self.poros_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(
                0.98*self.Wliq_top / self.rw_top, 1.0)) # relative evaporation rate (-)

        """ soil """
        # soil/peat type
        self.soiltype = spara['soiltype']

        # interpolated functions for soil column ground water dpeth vs. water storage, transmissivity etc.
        self.wsto_to_gwl = spara['wtso_to_gwl']
        self.gwl_to_wsto = spara['gwl_to_wsto']
        self.gwl_to_Tr = spara['gwl_to_Tr']
        self.gwl_to_C = spara['gwl_to_C']
        self.gwl_to_rootmoist = spara['gwl_to_rootmoist']

        # initial h (= gwl) and boundaries [m]
        self.h = spara['ground_water_level']
        self.h0 = -0.9  # depth of water level in ditch / canal, from input?!
        self.h[0,:] = self.h0
        self.h[-1,:] = self.h0
        self.h[:,0] = self.h0
        self.h[:,-1] = self.h0

        # soil surface elevation and hydraulic head [m]
        self.ele = spara['elevation']
        self.H = self.ele + self.h

        # water storage [m]
        self.Wsto_max = np.full_like(self.h, 0.0)  # storage of fully saturated profile
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_max[self.soiltype == key] = value(0.0)
        self.Wsto = np.full_like(self.h, 0.0)  # storage corresponding to h
        for key, value in self.gwl_to_wsto.items():
            self.Wsto[self.soiltype == key] = value(self.h[self.soiltype == key])

        # rootzone moisture [m3 m-3], parameters related to transpiration limit during dry conditions
        self.rootmoist = np.full_like(self.h, 0.0)
        self.root_fc0 = np.full_like(self.h, 0.0)
        self.root_fc1 = np.full_like(self.h, 0.0)
        self.root_wp = np.full_like(self.h, 0.0)
        for key, value in self.gwl_to_rootmoist.items():
            self.rootmoist[self.soiltype == key] = value(self.h[self.soiltype == key])
            self.root_fc0[self.soiltype == key] = value(-0.7 - 0.1)
            self.root_fc1[self.soiltype == key] = value(-1.2 - 0.1)
            self.root_wp[self.soiltype == key] = value(-150.0 - 0.1)

        self.Rew = 1.0

        """ parameters for 2D solution """
        # parameters for solving
        self.implic = 1  # solving method: 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
        self.DrIrr = False  # if true: constant head boundaries, if false changing boundary: gwl above canal, using constant head, if below below no flow

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
        self.A = np.zeros((self.n,self.n))

        self.CC = np.zeros((self.rows,self.cols))
        self.Tr0 = np.zeros((self.rows,self.cols))
        self.Tr1 = np.zeros((self.rows,self.cols))

    def rolling_window(self, a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

    def run_timestep(self, dt=1.0, rr=0.0, tr=0.0, evap=0.0):

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

        # for computing mass balance later
        state0 = self.Wsto + self.Wsto_top + rr

        # moss/organic layer interception & water balance
        interc = np.maximum(0.0, (self.Wsto_top_max - self.Wsto_top))\
                    * (1.0 - np.exp(-(rr / self.Wsto_top_max)))
        rr -= interc  # to soil profile
        self.Wsto_top += interc
        evap = np.minimum(evap, self.Wsto_top)
        self.Wsto_top -= evap

        # source/sink during dt [m]
        S = rr - tr
        # air volume [m]
        airv = self.Wsto_max - self.Wsto
        # remove water that cannot fit into soil
        S = np.minimum(S, airv)

        # inflow excess - either into moss or surface runoff
        exfil = rr - (S + tr)
        # water that can fit in top layer
        to_top_layer = np.minimum(exfil, self.Wsto_top_max - self.Wsto_top)
        self.Wsto_top += to_top_layer
        # route remaining to surface runoff
        surface_runoff = exfil - to_top_layer

        # Transmissivity of previous timestep [m2 d-1]
        for key, value in self.gwl_to_Tr.items():
            self.Tr0[self.soiltype == key] = value(self.h[self.soiltype == key])
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
        TrW0 = np.ravel(self.TrW0); TrE0 = np.ravel(self.TrE0)
        TrN0 = np.ravel(self.TrN0); TrS0 = np.ravel(self.TrS0)
        HW = np.ravel(self.HW); HE = np.ravel(self.HE)
        HN = np.ravel(self.HN); HS = np.ravel(self.HS)
        H = np.ravel(self.H)

        # hydrauli heads, new iteration and old iteration
        Htmp = self.H.copy()
        Htmp1 = self.H.copy()

        crit = 1e-6  # convergence criteria
        maxiter = 100

        for it in range(maxiter):

            # transmissivity [m2 d-1] to neighbouring cells with HTmp1
            for key, value in self.gwl_to_Tr.items():
                self.Tr1[self.soiltype == key] = value(Htmp1[self.soiltype == key] - self.ele[self.soiltype == key])
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
            for key, value in self.gwl_to_C.items():
                self.CC[self.soiltype == key] = value(Htmp1[self.soiltype == key] - self.ele[self.soiltype == key])
            alfa = np.ravel(self.CC * self.dxy**2 / dt)

            # Penta-diagonal matrix setup
            i,j = np.indices(self.A.shape)
            self.A[i==j] = self.implic * (TrW1 + TrE1 + TrN1 + TrS1) + alfa  # Diagonal
            self.A[i==j+1] = -self.implic * TrW1[1:]  # West element
            self.A[i==j-1] = -self.implic * TrE1[:-1]  # East element
            self.A[i==j+self.cols] = -self.implic * TrN1[self.cols:]  # North element
            self.A[i==j-self.cols] = -self.implic * TrS1[:self.n-self.cols]  # South element

            # Boundaries
            # North
            for k in range(0,self.cols):
                self.A[k,k] = self.implic * (2 * TrS1[k]) + alfa[k]
                self.A[k,k+1] = 0
                self.A[k,k-1] = 0
                self.A[k,k+self.cols] = -self.implic * 2 * TrS1[k]
            # South
            for k in range(self.n-self.cols,self.n-1):
                self.A[k,k] = self.implic * (2 * TrN1[k]) + alfa[k]
                self.A[k,k+1] = 0
                self.A[k, k-1] = 0
                self.A[k, k-self.cols] = -self.implic * 2 * TrN1[k]
            self.A[self.n-1,self.n-1] = self.implic * (2 * TrN1[self.n-1]) + alfa[self.n-1]
            self.A[self.n-1,self.n-2] = 0
            self.A[self.n-1,self.n-1-self.cols] = -self.implic * 2 * TrN1[self.n-1]
            # West
            for k in np.arange(self.cols, self.n-self.cols, self.cols, dtype='int'):
                self.A[k,k] = self.implic * (2 * TrE1[k]) + alfa[k]
                self.A[k, k-self.cols] = 0
                self.A[k, k+self.cols] = 0
                self.A[k,k+1] = -self.implic * 2 * TrE1[k]
            # East
            for k in np.arange(self.cols-1, self.n-self.cols, self.cols, dtype='int'):
                self.A[k,k]= self.implic*(2*TrW1[k])+alfa[k]
                if k-self.cols > 0:
                    self.A[k, k-self.cols]=0
                    self.A[k,k-1]=-self.implic*2*TrW1[k]
                    self.A[k, k+self.cols]=0

            # Knowns: Right hand side of the eq
            hs = (np.ravel(S) * dt * self.dxy**2 + alfa*H
                  + (1.-self.implic) * (TrN0*HN) + (1.-self.implic) * (TrW0*HW)
                  - (1.-self.implic) * (TrN0 + TrW0 + TrE0 + TrS0) * H
                  + (1.-self.implic) * (TrE0*HE) + (1.-self.implic) * (TrS0*HS))

            # Solve: A*Htmp1 = hs
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(self.A),hs])
            Htmp1 = np.reshape(Htmp1,(self.rows,self.cols))

            # Return boundaries --- !?
            Htmp1[0,:] = self.ele[0,:] + self.h0
            Htmp1[-1,:] = self.ele[-1,:] + self.h0
            Htmp1[:,0] = self.ele[:,0] + self.h0
            Htmp1[:,-1] = self.ele[:,-1] + self.h0

            Htmp1 = np.where(Htmp1 > self.ele, self.ele, Htmp1)
            conv1 = np.max(np.abs(Htmp1 - Htmp))
            Htmp = Htmp1.copy()
            if conv1 < crit:
                # print('iterations', it, np.average(self.ele[1:-1,1:-1]-Htmp[1:-1,1:-1]))
                break
            # end of iteration loop

        if self.DrIrr == True:
            Htmp1[0,:] = self.ele[0,:] + self.h0
            Htmp1[-1,:] = self.ele[-1,:] + self.h0
            Htmp1[:,0] = self.ele[:,0] + self.h0
            Htmp1[:,-1] = self.ele[:,-1] + self.h0
        else:
            # if neighboring node has lower H set ditch H equal to that (no-flow), otherwise set to water level in ditch (constant head)
            Htmp1[0,:] = np.where(Htmp1[1,:] < self.ele[0,:] + self.h0,
                                  Htmp1[1,:], self.ele[0,:] + self.h0)
            Htmp1[-1,:] = np.where(Htmp1[-2,:] < self.ele[-1,:] + self.h0,
                                  Htmp1[-2,:], self.ele[-1,:] + self.h0)
            Htmp1[:,0] = np.where(Htmp1[:,1] < self.ele[:,0] + self.h0,
                                  Htmp1[:,1], self.ele[:,0] + self.h0)
            Htmp1[:,-1] = np.where(Htmp1[:,-2] < self.ele[:,-1] + self.h0,
                                  Htmp1[:,-2], self.ele[:,-1] + self.h0)

        """ update state """
        # soil profile
        self.H = Htmp1.copy()
        self.h = self.H - self.ele
        for key, value in self.gwl_to_wsto.items():
            self.Wsto[self.soiltype == key] = value(self.h[self.soiltype == key])

        # organic top layer; maximum that can be hold is Fc
        self.Wliq_top = self.fc_top * self.Wsto_top / self.Wsto_top_max
        self.Ree = np.maximum(0.0, np.minimum(0.98*self.Wliq_top / self.rw_top, 1.0))

        for key, value in self.gwl_to_rootmoist.items():
            self.rootmoist[self.soiltype == key] = value(self.h[self.soiltype == key])

        # Koivusalo et al. 2008 HESS without wet side limit
        self.Rew = np.where(self.rootmoist > self.root_fc1,
                            np.minimum(1.0, 0.5*(1 + (self.rootmoist - self.root_fc1)/(self.root_fc0 - self.root_fc1))),
                            np.maximum(0.0, 0.5*(self.rootmoist - self.root_wp)/(self.root_fc1 - self.root_wp))
                            )

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:-1]
        self.HE[:,:-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:-1,:]
        self.HS[:-1,:] = self.H[1:,:]

        lateral_flow = (self.TrW1*(self.H - self.HW)
                        + self.TrE1*(self.H - self.HE)
                        + self.TrN1*(self.H - self.HN)
                        + self.TrS1*(self.H - self.HS)) / self.dxy**2

        # mass balance error [m]
        mbe = (state0  - self.Wsto_top - self.Wsto -
               tr - surface_runoff - evap - lateral_flow)

        results = {
                'ground_water_level': self.h,  # [m]
                'infiltration': (S - tr) * 1e3,  # [mm d-1]
                'surface_runoff': surface_runoff * 1e3,  # [mm d-1]
                'evaporation': evap * 1e3,  # [mm d-1]
                'drainage': lateral_flow * 1e3,  # [mm d-1]
                'moisture_top': self.Wliq_top,  # [m3 m-3]
                'water_closure': mbe * 1e3,  # [mm d-1]
                'transpiration_limitation': self.Rew,  # [-]
                'rootzone_moisture': self.rootmoist,  # [m3 m-3]
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

    z = np.array(z)
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:]
    z_mid = dz / 2 - np.cumsum(dz)

    # --------- connection between gwl and Wsto, Tr, C------------
    # gwl from ground surface gwl = 0 to gwl = -5
    gwl = np.arange(0.0, -5, -1e-3)
    gwl[-1] = -150
    # solve water storage corresponding to gwls
    Wsto = [sum(h_to_cellmoist(pF, g - z_mid, dz) * dz) for g in gwl]

    if root:
        Wsto = Wsto/sum(dz)
        GwlToWsto = interp1d(np.array(gwl), np.array(Wsto), fill_value='extrapolate')
        return {'to_rootmoist': GwlToWsto}

    # solve tranmissivity corresponding to gwls
    Tr = [transmissivity(dz, Ksat, g) * 86400. for g in gwl]  # [m2 d-1]

    # interpolate functions
    WstoToGwl = interp1d(np.array(Wsto), np.array(gwl), fill_value='extrapolate')
    GwlToWsto = interp1d(np.array(gwl), np.array(Wsto), fill_value='extrapolate')
    GwlToC = interp1d(np.array(gwl), np.array(np.gradient(Wsto)/np.gradient(gwl)), fill_value='extrapolate')
    GwlToTr = interp1d(np.array(gwl), np.array(Tr), fill_value='extrapolate')

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
    Ka = 0.0
    Tr = 0.0

    ib = max(abs(z))

    Hdr = min(max(0, gwl + ib), ib)  # depth of saturated layer above impermeable bottom

    """ drainage from saturated layers above ditch base """
    # layers above ditch bottom where drainage is possible
    ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -ib))

    if Hdr > 0:
        # saturated layer thickness [m]
        dz_sat = np.minimum(np.maximum(gwl - (z - dz / 2), 0), dz)
        # transmissivity of layers  [m2 s-1]
        Trans = Ksat * dz_sat

        """ drainage from saturated layers above ditch base """
        # layers above ditch bottom where drainage is possible
        ix = np.intersect1d(np.where((z - dz / 2)- gwl < 0), np.where(z + dz / 2 > -ib))

        if ix.size > 0:
            dz_sat[ix[-1]] = dz_sat[ix[-1]] + (z[ix][-1] - dz[ix][-1] / 2 + ib)
            if abs(sum(dz_sat[ix]) - Hdr) > eps:
                print(sum(dz_sat[ix]), Hdr, DitchDepth, gwl)
            Trans[ix[-1]] = Ksat[ix[-1]] * dz_sat[ix[-1]]
            Tr = sum(Trans[ix])
    return Tr

nan_function = interp1d(np.array([np.nan, np.nan]),
                        np.array([np.nan, np.nan]),
                        fill_value='extrapolate')