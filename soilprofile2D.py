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

class SoilGrid_2Dflow(object):
    """
    2D soil water flow model based on Ari Lauren SUSI2D
    Simulates deep soil water storage, and drainage to streams.
    """
    def __init__(self, spara):
        """
        Initializes SoilProfile2D:
        Args:
            spara (dict):
                'deep_id': deep soil id
                'elevation': elevation [m]
                'streams': stream water level [m], < 0 for streams otherwise 0
                'lakes': lake water level [m], < 0 for streams otherwise 0
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
        
        # catchment mask
        self.cmask = np.full_like(spara['deep_id'], np.nan)
        self.cmask[np.isfinite(spara['deep_id'])] = 1.0

        # stream parameters (currently np.nan for non-streams) 
        self.ditch_l = spara['stream_length'] # total stream length
        self.ditch_w = spara['stream_width'] # average stream width
        self.ditch_d = spara['stream_distance'] # average distance to stream (only for stream grid-cells)

        # interpolated functions for soil column groundwater depth vs. water storage, transmissivity etc.
        self.wsto_to_gwl = spara['wtso_to_gwl']
        self.gwl_to_wsto = spara['gwl_to_wsto']
        self.gwl_to_Tr = spara['gwl_to_Tr']
        self.gwl_to_C = spara['gwl_to_C']
        self.gwl_to_rootmoist = spara['gwl_to_rootmoist']

        # initial h (= gwl) and boundaries [m]
        self.ditch_h = spara['streams']
        self.lake_h = spara['lakes']
        self.gwl = spara['ground_water_level']
        # soil surface elevation and hydraulic head [m]
        self.ele = spara['elevation']
        self.H = self.ele + self.gwl
        # Identify lake interior and boundary cells
        lake_boundary = np.zeros_like(self.lake_h)
        self.lake_interior = np.zeros_like(self.lake_h)

        # grid
        self.rows = np.shape(self.gwl)[0]
        self.cols = np.shape(self.gwl)[1]
        self.n = self.rows * self.cols  # length of flattened array
        self.dxy = spara['dxy']  # horizontal distance between nodes dx=dy [m]
        
        # Iterate through each cell in lake_h to find if its lake boundary
        for i in range(self.rows):
            for j in range(self.cols):
                if self.lake_h[i, j] < -eps:
                    # Check if this cell is a boundary
                    is_boundary = False
                    # Check all 4 possible neighbors (west, east, north, south)
                    if (j > 0 and self.lake_h[i,j-1] == 0) or (j < self.cols - 1 and self.lake_h[i,j+1] == 0) or \
                       (i > 0 and self.lake_h[i-1,j] == 0) or (i < self.rows - 1 and self.lake_h[i+1,j] == 0):
                        is_boundary = True
                    if is_boundary:
                        lake_boundary[i,j] = 1       

        self.lake_interior[(lake_boundary != 1) & (self.lake_h < -eps)] = 1 # saving lake interior array

        # nan to lake interiors (lake interiors should not be solved)
        self.soiltype[self.lake_interior == 1] = np.nan
        self.cmask[self.lake_interior == 1] = np.nan
        self.gwl[self.lake_interior == 1] = np.nan
        self.H[self.lake_interior == 1] = -999

        # lower boundaries
        #print('spara[deep_z]', spara['deep_z'])
        self.deep_z = spara['deep_z']*-1
        self.bedrock_h = self.ele + self.deep_z
        self.bedrock_h[self.lake_interior == 1] = np.nan
        self.bedrock_h[self.lake_interior == 1] = -999
        self.bedrock_h = np.ravel(self.bedrock_h)

        # replace nans (values outside catchment area)
        self.H[np.isnan(self.H)] = -999

        # water storage [m]
        self.Wsto_deep_max = np.full_like(self.gwl, 0.0)  # storage of fully saturated profile
        self.Wsto_deep = np.full_like(self.gwl, 0.0)  

        # rootzone moisture [m3 m-3]
        #self.deepmoist = np.full_like(self.gwl, 0.0)

        # self.z_from_gis == True OR False
        # determines whether the deep_z and thus interpolation functions are made cell-wise (True) or soiltype-wise (False)
        # if-elif statements later in the code made accordingly
        self.z_from_gis = (
            isinstance(spara['wtso_to_gwl'], np.ndarray) and
            spara['deep_id'].shape == spara['wtso_to_gwl'].shape
            )
        
        # initial water storages according to gwl
        if self.z_from_gis == False: # soiltype-wise calculation
            for key, value in self.gwl_to_wsto.items():
                self.Wsto_deep_max[self.soiltype == key] = value(0.0)
                self.Wsto_deep[self.soiltype == key] = value(self.gwl[self.soiltype == key]) # storage corresponding to h
            #for key, value in self.gwl_to_rootmoist.items():
            #    self.deepmoist[self.soiltype == key] = value(self.gwl[self.soiltype == key])
        elif self.z_from_gis == True: # cell-wise calculation
            for i in range(self.gwl_to_wsto.shape[0]):
                for j in range(self.gwl_to_wsto.shape[1]):
                    if np.isfinite(self.cmask[i,j]): 
                        self.Wsto_deep_max[i,j] = self.gwl_to_wsto[i,j](0.0) # max storage with gwl = 0
                        self.Wsto_deep[i,j] = self.gwl_to_wsto[i,j](self.gwl[i,j]) # storage corresponding to h
                        #self.deepmoist[i,j] = self.gwl_to_rootmoist[i,j](self.gwl[i,j])
            
        #self.deepmoist[np.isnan(self.gwl)] = np.nan

        # air volume and returnflow
        self.airv_deep = np.maximum(0.0, self.Wsto_deep_max - self.Wsto_deep)
        #self.qr = np.full_like(self.gwl, 0.0)

        """ parameters for 2D solution """
        # parameters for solving
        # 0.5 seems to work better when gwl is close to impermeable bottom
        # (probably because transmissivity does not switch between 0. and > 0 as much)
        self.implic = 0.5  # solving method: 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson

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

        state0 = self.Wsto_deep + S  # [m]

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
        ditch_l = np.ravel(self.ditch_l)
        ditch_w = np.ravel(self.ditch_w)
        ditch_d = np.ravel(self.ditch_d)
        lake_h = np.ravel(self.lake_h)
        lake_interior = np.ravel(self.lake_interior)
        ele = np.ravel(self.ele)

        # ditch drainage [m] - outside iteration loop to avoid ditch switching on and off during iteration
        Ksat = 1E-04 * 86400.
        S_dd = np.where((ditch_h < -eps) & (ele + ditch_h < H), 
                        Ksat * ditch_l * ditch_w / ditch_d * (H - (ele + ditch_h)) * dt / self.dxy**2,
                        0.0)

        # Lakes
        # calculate mean H of neighboring nodes to find out whether lake is active (constant head)
        # from previous timestep to avoid switching boundary/no-boundary during iteration
        H_neighbours = lake_h.copy()
        for k in np.where(lake_h < -eps)[0]:
            H_ave = 0
            n_neigh = 0
            if k%self.cols != 0 and lake_h[k-1] > -eps: # west non-lake neighbor
                    H_ave += H[k-1]
                    n_neigh += 1
            if (k+1)%self.cols != 0 and lake_h[k+1] > -eps: # east non-lake neighbor
                    H_ave += H[k+1]
                    n_neigh += 1
            if k-self.cols >= 0 and  lake_h[k-self.cols] > -eps: # north non-lake neighbor
                    H_ave += H[k-self.cols]
                    n_neigh += 1
            if k+self.cols < self.n and lake_h[k+self.cols] > -eps: # south non-lake neighbor
                    H_ave += H[k+self.cols]
                    n_neigh += 1
            if n_neigh > 0:
                H_neighbours[k] = H_ave / n_neigh  # average of neighboring non-lake nodes
            else:  # corners or nodes surrounded by lakes dont have neighbors, given its ditch depth
                H_neighbours[k] = ele[k] + lake_h[k] + eps

        # lake interior do not have neighbours
        for k in np.where(lake_interior == 1)[0]:
            H_neighbours[k] = ele[k] + lake_h[k] + eps
        
        H_neighbours_2d = np.reshape(H_neighbours,(self.rows,self.cols)) # 2D array of lakes' neighbouring land cells' average hydraulic head

        # Transmissivity of previous timestep [m2 d-1]
        # for lake nodes that are active, transmissivity calculated based on mean H of
        # neighboring nodes, not lake depth which would restrict tranmissivity too much
        # Whole profile depth still considered, but I think that is the usual way..
        H_for_Tr = np.where((self.lake_h < -eps) & (H_neighbours_2d > self.ele + self.lake_h),
                            H_neighbours_2d, self.H)
        
        # transmissivities based on gwl
        if self.z_from_gis == False:
            for key, value in self.gwl_to_Tr.items():
                self.Tr0[self.soiltype == key] = value(H_for_Tr[self.soiltype == key] - self.ele[self.soiltype == key])
        elif self.z_from_gis == True:
            for i in range(self.gwl_to_Tr.shape[0]):
                for j in range(self.gwl_to_Tr.shape[1]):
                    if np.isfinite(self.cmask[i,j]): 
                        self.Tr0[i,j] = self.gwl_to_wsto[i,j](H_for_Tr[i,j] - self.ele[i,j])

        # transmissivity at all four sides of the element is computed as geometric mean of surrounding element transimissivities
        # is this actually at all four sides, or just along east-west and north-sound axes?
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

        # from previous timestep
        TrW1 = TrW0.copy()
        TrE1 = TrE0.copy()
        TrN1 = TrN0.copy()
        TrS1 = TrS0.copy()

        # hydraulic heads, new iteration and old iteration
        Htmp = self.H.copy()
        Htmp1 = self.H.copy()

        # convergence criteria
        crit = 1e-3  # seems mass balance error remains resonable
        maxiter = 100
        update_Tr_in_loop = True

        for it in range(maxiter):
            if update_Tr_in_loop:
                # transmissivity [m2 d-1] to neighbouring cells with HTmp1
                # for lake nodes that are active, transmissivity calculated based on mean H of
                # neighboring nodes, not lake depth which would restrict transmissivity too much
                H_for_Tr = np.where((self.lake_h < -eps) & (H_neighbours_2d > self.ele + self.lake_h),
                                    H_neighbours_2d, Htmp)
                # transmissivities based on gwl
                if self.z_from_gis == False:
                    for key, value in self.gwl_to_Tr.items():
                        self.Tr1[self.soiltype == key] = value(H_for_Tr[self.soiltype == key] - self.ele[self.soiltype == key])
                elif self.z_from_gis == True:
                    for i in range(self.gwl_to_Tr.shape[0]):
                        for j in range(self.gwl_to_Tr.shape[1]):
                            if np.isfinite(self.cmask[i,j]): 
                                self.Tr1[i,j] = self.gwl_to_wsto[i,j](H_for_Tr[i,j] - self.ele[i,j])            
                
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
            else:
                # from previous timestep
                TrW1 = TrW0.copy()
                TrE1 = TrE0.copy()
                TrN1 = TrN0.copy()
                TrS1 = TrS0.copy()

            # differential water capacity dSto/dh
            if self.z_from_gis == False:
                for key, value in self.gwl_to_C.items():
                    self.CC[self.soiltype == key] = value(Htmp[self.soiltype == key] - self.ele[self.soiltype == key])
                for key, value in self.gwl_to_wsto.items():
                    self.Wtso1_deep[self.soiltype == key] = value(Htmp[self.soiltype == key] - self.ele[self.soiltype == key])
            elif self.z_from_gis == True:
                for i in range(self.gwl_to_C.shape[0]):
                    for j in range(self.gwl_to_C.shape[1]):
                        if np.isfinite(self.cmask[i,j]): 
                            self.CC[i,j] = self.gwl_to_C[i,j](Htmp[i,j] - self.ele[i,j])
                            self.Wtso1_deep[i,j] = self.gwl_to_wsto[i,j](Htmp[i,j] - self.ele[i,j])

            alfa = np.ravel(self.CC * self.dxy**2 / dt)
            # alfa = np.ravel((0.5*self.CC + 0.5*CCtmp) * self.dxy**2 / dt)

            # Setup of diagonal sparse matrix
            a_d = self.implic * (TrW1 + TrE1 + TrN1 + TrS1) + alfa  # Diagonal
            a_w = -self.implic * TrW1[1:]  # West element
            a_e = -self.implic * TrE1[:-1]  # East element
            a_n = -self.implic * TrN1[self.cols:]  # North element
            a_s = -self.implic * TrS1[:self.n-self.cols]  # South element

            # Knowns: Right hand side of the eq
            Htmp = np.ravel(Htmp)  
            hs = (np.ravel(S) * self.dxy**2 / dt + alfa * Htmp - S_dd * self.dxy**2 / dt
                  - np.ravel(self.Wtso1_deep) * self.dxy**2 / dt + Wsto_deep * self.dxy**2 / dt
                  + (1.-self.implic) * (TrN0*HN) + (1.-self.implic) * (TrW0*HW)
                  - (1.-self.implic) * (TrN0 + TrW0 + TrE0 + TrS0) * H
                  + (1.-self.implic) * (TrE0*HE) + (1.-self.implic) * (TrS0*HS))

            # Lakes
            for k in np.where(lake_h < -eps)[0]:
                # update a_x and hs when lake as constant head
                if H_neighbours[k] > ele[k] + lake_h[k]: # average of neighboring H above lake depth
                    hs[k] = ele[k] + lake_h[k]
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

            #Htmp1 = np.fmax(self.bedrock_h, Htmp1) # limit to bedrock

            # testing, limit change
            #if np.any(np.abs(Htmp1-Htmp)) > 0.5:
                #print('Difference greater than 0.5')

            max_index_print = np.unravel_index(np.argmax(np.abs(Htmp1 - Htmp)),(self.rows,self.cols))
            Htmp_print = np.reshape(Htmp,(self.rows,self.cols))
            Htmp1_print = np.reshape(Htmp1,(self.rows,self.cols))
            
            #print('\t', 'iterations:', it,
            #        ' max_index:', max_index_print,
            #        ' H[max_index]', Htmp_print[max_index_print]-self.ele[max_index_print], 
            #        ' H1[max_index]', Htmp1_print[max_index_print]-self.ele[max_index_print])


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

            if np.any(np.abs(Htmp1-Htmp)) > 0.5:
                Htmp_print = np.reshape(Htmp,(self.rows,self.cols))
                Htmp1_print = np.reshape(Htmp1,(self.rows,self.cols))
                print('Difference greater than 0.5')
                print('\t', 'iterations:', it, ' con1:', conv1, 
                      ' max_index:', max_index, ' self.ditch_h[max_index]', self.ditch_h[max_index],
                      ' H[max_index]', Htmp_print[max_index]-self.ele[max_index], 
                      ' H1[max_index]', Htmp1_print[max_index]-self.ele[max_index])

            Htmp = np.reshape(Htmp,(self.rows,self.cols))

            # print to get sense what's happening when problems in convergence
            if it > 90:
                print('\t', 'iterations:', it, ' con1:', conv1, 
                      ' max_index:', max_index, ' self.ditch_h[max_index]', self.ditch_h[max_index],
                      ' H[max_index]', Htmp[max_index]-self.ele[max_index])
                
                Htmp_N_temp = Htmp[(max_index[0]-1, max_index[1])]
                Htmp_S_temp = Htmp[(max_index[0]+1, max_index[1])]
                Htmp_W_temp = Htmp[(max_index[0], max_index[1]-1)]
                Htmp_E_temp = Htmp[(max_index[0], max_index[1]+1)]

                H_N_temp = Htmp[(max_index[0]-1, max_index[1])] - self.ele[(max_index[0]-1, max_index[1])]
                H_S_temp = Htmp[(max_index[0]+1, max_index[1])] - self.ele[(max_index[0]+1, max_index[1])]
                H_W_temp = Htmp[(max_index[0], max_index[1]-1)] - self.ele[(max_index[0], max_index[1]-1)]
                H_E_temp = Htmp[(max_index[0], max_index[1]+1)] - self.ele[(max_index[0], max_index[1]+1)]

                TrN1_temp = np.reshape(TrN1,(self.rows,self.cols))[(max_index[0]-1, max_index[1])]
                TrS1_temp = np.reshape(TrS1,(self.rows,self.cols))[(max_index[0]+1, max_index[1])]
                TrW1_temp = np.reshape(TrW1,(self.rows,self.cols))[(max_index[0], max_index[1]-1)]
                TrE1_temp = np.reshape(TrE1,(self.rows,self.cols))[(max_index[0], max_index[1]+1)]

                #print('Max Htmp1-Htmp', np.max(np.abs(Htmp1-Htmp)))
                print('Htmp of max_index:', Htmp[max_index])
                print('Htmp of N neighbor:', Htmp_N_temp)
                print('Htmp of S neighbor:', Htmp_S_temp)
                print('Htmp of W neighbor:', Htmp_W_temp)
                print('Htmp of E neighbor:', Htmp_E_temp)

                print('H of max_index:', Htmp[max_index] - self.ele[max_index])
                print('H of N neighbor:', H_N_temp)
                print('H of S neighbor:', H_S_temp)
                print('H of W neighbor:', H_W_temp)
                print('H of E neighbor:', H_E_temp)                

                print('Tr of N neighbor:', TrN1_temp)
                print('Tr of S neighbor:', TrS1_temp)
                print('Tr of W neighbor:', TrW1_temp)
                print('Tr of E neighbor:', TrE1_temp)

            if conv1 < crit:
                break
            # end of iteration loop
        if it == 99:
            self.conv99 +=1
        #self.totit += it
        Htmp = np.reshape(Htmp,(self.rows,self.cols))

        print('Timestep:', self.tmstep, ', iterations:', it, ', conv1:', conv1, ', H[max_index]:', Htmp[max_index]-self.ele[max_index])
        
        # lateral flow [m d-1] is calculated in two parts: one depending on previous time step
        # and other on current time step (lateral flowsee 2/2). Their weighting depends
        # on self.implic
        # lateral flow 1/2 with old heads (and old transmissivities if updated inside the loop)
        # use 1-self.implic
        lateral_flow = ((1-self.implic)*(self.TrW0*(self.H - self.HW)
                        + self.TrE0*(self.H - self.HE)
                        + self.TrN0*(self.H - self.HN)
                        + self.TrS0*(self.H - self.HS)))/ self.dxy**2

        """ update state """
        # soil profile
        self.H = Htmp.copy()
        self.gwl = self.H - self.ele

        # water storages according to new gwl
        if self.z_from_gis == False:
            for key, value in self.gwl_to_wsto.items():
                self.Wsto_deep[self.soiltype == key] = value(self.gwl[self.soiltype == key])
        elif self.z_from_gis == True:
            for i in range(self.gwl_to_wsto.shape[0]):
                for j in range(self.gwl_to_wsto.shape[1]):
                    if np.isfinite(self.cmask[i,j]):
                        self.Wsto_deep[i,j] = self.gwl_to_wsto[i,j](self.gwl[i,j])

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:-1]
        self.HE[:,:-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:-1,:]
        self.HS[:-1,:] = self.H[1:,:]

        # lateral flow 2/2 [m d-1] with new heads (and new transmissivities if updated inside the loop)
        # use self.implic here
        lateral_flow += (self.implic*(self.TrW1*(self.H - self.HW)
                        + self.TrE1*(self.H - self.HE)
                        + self.TrN1*(self.H - self.HN)
                        + self.TrS1*(self.H - self.HS)))/ self.dxy**2
        
        Tr = np.nanmean([self.TrW1, self.TrE1, self.TrN1, self.TrS1], axis=0)
        Tr = np.reshape(Tr,(self.rows,self.cols))

        """ new update state """
        # state0 = water storage at timestep0 (including bucket drainage at timestep0)
        # Wsto_deep = water storage at timestep1
        # lateral_flow = lateral flow of each grid-cell as calculated in two parts earlier
        #   lateral_flow is positive when flow going out of the grid cells (saved as -lateral_flow)
        # mbe = state0 - Wsto_deep - lateral_flow * dt - S_dd (makes sense)

        # Let's limit head to 0 and assign rest as return flow to bucketgrid
        Wsto_before_qr = self.Wsto_deep.copy()

        # restricting gwl to 0 on land and, and to minimum of gwl and lake_h in lakes
        self.gwl = np.where(self.lake_h < -eps, np.minimum(self.gwl, self.lake_h), np.minimum(0.0, self.gwl))
        self.H = self.gwl + self.ele
        self.H[np.isnan(self.H)] = -999

        # Updating the storage according to new head
        if self.z_from_gis == False:
            for key, value in self.gwl_to_wsto.items():
                self.Wsto_deep[self.soiltype == key] = value(self.H[self.soiltype == key] - self.ele[self.soiltype == key])
            #for key, value in self.gwl_to_rootmoist.items():
            #    self.deepmoist[self.soiltype == key] = value(self.gwl[self.soiltype == key])
        elif self.z_from_gis == True:
            for i in range(self.gwl_to_wsto.shape[0]):
                for j in range(self.gwl_to_wsto.shape[1]):
                    if np.isfinite(self.cmask[i,j]): 
                        self.Wsto_deep[i,j] = self.gwl_to_wsto[i,j](self.H[i,j] - self.ele[i,j])  
                        #self.deepmoist[i,j] = self.gwl_to_rootmoist[i,j](self.gwl[i,j])

        # The difference is the return flow to bucketgrid [m]
        # qr = np.where(self.lake_h > -eps, Wsto_before_qr - self.Wsto_deep, 0.)  # should be zero for lake cells anyway?
        qr = Wsto_before_qr - self.Wsto_deep

        # Lakes are described as constant heads so the netflow to lakes can
        # be calculated from their mass balance [m]  !? does not work if head below lake ?!
        netflow_to_lake = np.where((self.lake_h < -eps) & (H_neighbours_2d > self.ele + self.lake_h),  # active lake cells
                                    state0 - np.reshape(S_dd,(self.rows,self.cols))  # ditches in lake cells?
                                    - self.Wsto_deep - qr - lateral_flow * dt, 0.0)

        # air volume
        self.airv_deep = np.maximum(0.0, self.Wsto_deep_max - self.Wsto_deep)
        
        # mass balance error [m] - now calculated also for lake cells which is relevant when H < lake_h
        mbe = (state0 - np.reshape(S_dd,(self.rows,self.cols)) - self.Wsto_deep - qr - lateral_flow * dt
               - netflow_to_lake * dt)
        mbe = np.where(self.lake_h < -eps, 0.0, mbe)

        # outputs multiplied by cmask
        h_out = self.gwl.copy() * self.cmask
        lateral_flow = lateral_flow * self.cmask
        netflow_to_lake = netflow_to_lake * self.cmask
        netflow_to_ditch = np.reshape(S_dd,(self.rows,self.cols)) * self.cmask
        mbe = mbe * self.cmask
        #deepmoist_out = self.deepmoist.copy() * self.cmask
        Wsto_deep_out = self.Wsto_deep.copy() * self.cmask

        results = {
                'ground_water_level': h_out,  # [m]
                'lateral_netflow': -lateral_flow * 1e3 / dt,  # [mm d-1]
                'netflow_to_ditch': netflow_to_ditch * 1e3 / dt,  # [mm d-1]
                'netflow_to_lake': netflow_to_lake * 1e3 / dt,  # [mm d-1]
                'water_closure': mbe * 1e3 / dt,  # [mm d-1]
                #'moisture_deep': deepmoist_out,  # [m3 m-3]
                'water_storage': Wsto_deep_out * 1e3, # [mm]
                'return_flow': qr * 1e3, # [mm],
                'transmissivity': Tr,  # [m2 d-1]
                }
        return results


def gwl_Wsto(z, pF, grid_step=-0.01, Ksat=None, root=False):
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
    z = np.array(z, dtype=np.float64) # profile depths
    dz = abs(z)
    dz[1:] = z[:-1] - z[1:] # profile depths into profile thicknesses

    # finer grid for calculating wsto to avoid discontinuity in C (dWsto/dGWL)
    z_fine= (np.arange(0, min(z), grid_step) - grid_step).astype(np.float64)
    dz_fine = z_fine*0.0 - grid_step
    z_mid_fine = dz_fine / 2 - np.cumsum(dz_fine)

    ix = np.zeros(len(z_fine), dtype=np.float64)

    for depth in z:
        # below makes sure floating point precision doesnt mess with the ix
        ix += np.where((z_fine < depth) & ~np.isclose(z_fine, depth, atol=1e-9), 1, 0)

    pF_fine={}
    for key in pF.keys():
        pp = []
        for i in range(len(z_fine)):
            pp.append(pF[key][int(ix[i])])
        pF_fine.update({key: np.array(pp)})

    # --------- connection between gwl and Wsto, Tr, C------------
    gwl = np.arange(1.0, min(z)-5, grid_step)
    # solve water storage corresponding to gwls
    Wsto_deep = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine)
            + max(0.0,g) for g in gwl]  # water storage above ground surface == gwl

    if root:
        Wsto_deep = [sum(h_to_cellmoist(pF_fine, g - z_mid_fine, dz_fine) * dz_fine) for g in gwl]
        Wsto_deep = Wsto_deep/sum(dz)
        GwlToWsto = interp1d(np.array(gwl), np.array(Wsto_deep), fill_value='extrapolate')
        return {'to_rootmoist': GwlToWsto}

    # solve transmissivity corresponding to gwls
    Tr = [transmissivity(dz, Ksat, g) * 86400. for g in gwl]  # [m2 d-1]

    #print('np.array(gwl).shape', np.array(gwl).shape)
    #print('np.array(Wsto_deep).shape', np.array(Wsto_deep).shape)

    # interpolate functions
    WstoToGwl = interp1d(np.array(Wsto_deep), np.array(gwl), fill_value='extrapolate')
    GwlToWsto = interp1d(np.array(gwl), np.array(Wsto_deep), fill_value='extrapolate')
    GwlToC = interp1d(np.array(gwl), np.array(np.gradient(Wsto_deep)/np.gradient(gwl)), fill_value='extrapolate')
    GwlToTr = interp1d(np.array(gwl), np.array(Tr), fill_value='extrapolate')
    
    plt.figure(1)
    plt.plot(np.array(gwl), np.array(np.gradient(Wsto_deep/np.gradient(gwl))))
    plt.figure(2)
    plt.plot(np.array(gwl), np.log10(np.array(Tr)))
    plt.plot(np.array(gwl), np.array(Tr))
    plt.figure(3)
    plt.plot(np.array(gwl), np.array(Wsto_deep))

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
       dz (array):  soil compartment thickness, node in center [m]
       Ksat (array): horizontal saturated hydr. cond. [ms-1]
       gwl (float): ground water level below surface, <0 [m]

    Returns:
       Tr (array): tranmissivity [m2 s-1]
    """
    
    # midpoint of cell, soil surface at 0
    z = dz / 2 - np.cumsum(dz)

    # saturated layer thickness [m], between [0, dz]
    dz_sat = np.minimum(np.maximum(gwl - (z - dz / 2), 0), dz)
    # in top cell, allow transmissivity to increases when gwl above ground surface level
    dz_sat[0] = np.maximum(gwl - (z[0] - dz[0] / 2), 0)  
    
    # transmissivity of layers  [m2 s-1]
    Trans = Ksat * dz_sat

    # sum over layers
    Tr = np.maximum(sum(Trans), 1e-4 / 86400)

    return Tr


def gwl_Wsto_vectorized(z, pF, grid_step=-0.01, Ksat=None, root=False):
    r""" Forms interpolated function for soil column ground water dpeth, < 0 [m], as a
    function of water storage [m] and vice versa + others

    Args:
        - pF (np.ndarray):
            - dict
                - 'ThetaS' (np.ndarray): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                - 'ThetaR' (np.ndarray): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                - 'alpha' (np.ndarray): air entry suction [cm\ :sup:`-1`]
                - 'n' (np.ndarray): pore size distribution [-]
        - z (np.ndarrays): soil compartment thichness, node in center [m]
    Returns:
        - (np.ndarray):
            - dict
                - 'to_gwl' (np.ndarray): interpolated function for gwl(Wsto)
                - 'to_wsto' (np.ndarray): interpolated function for Wsto(gwl)
                - 'to_C' (np.ndarray): interpolated function for C(Wsto)
                - 'to_Tr' (np.ndarray): interpolated function for Tr(gwl)
    """
    # Ensure z is a NumPy array
    z = np.array(z, dtype=np.float32)
    pF = np.array(pF)

    if z.ndim == 1:
        z = np.expand_dims(z, axis=0)
        pF = np.expand_dims(pF, axis=0)
        if Ksat is not None:
            Ksat = np.array(Ksat)
            Ksat = np.expand_dims(Ksat, axis=0)

    dz = np.abs(z)
    dz = np.hstack((dz[:, :1], np.diff(dz, axis=1)))
    
    if isinstance(grid_step, float):
        z_min = np.min(z, axis=1)
        max_len = int(np.abs(np.nanmin(z_min)) / np.abs(grid_step)) + 1
        z_fine = np.tile(np.arange(0., grid_step * max_len, grid_step), (z.shape[0], 1)) + grid_step
        z_fine = z_fine.astype(np.float32)
        z_fine[z_fine < z_min[:, None]] = np.nan
        dz_fine = z_fine*0.0 - grid_step
        z_mid_fine = dz_fine / 2 - np.cumsum(dz_fine, axis=1)
        ix = np.full((z_fine.shape), np.nan)
        # Expand z along the second axis to match z_fine's shape (broadcasting)
        z_expanded = np.expand_dims(z, axis=1)  # Shape: (rows, 1, cols)
        z_fine_expanded = np.expand_dims(z_fine, axis=2)  # Shape: (rows, fine_steps, 1)    

    elif grid_step == 'var':
        z_min = np.min(z, axis=1).astype(np.float32)
        z_min_min = np.nanmin(z_min)
        limits = [-1, -5, z_min_min]
        steps = [-0.01, -0.05, -0.3]
        z1 = np.arange(0+steps[0], limits[0], steps[0])
        z2 = np.arange(limits[0], limits[1], steps[1])
        z3 = np.arange(limits[1], limits[2] + steps[2], steps[2])  # Ensure we reach z_min
        # Combine all segments
        z_values = np.concatenate([z1, z2, z3])  # Ensure exact z_min
        z_fine = np.tile(z_values, (z.shape[0], 1))
        z_fine = z_fine.astype(np.float32)
        z_fine[z_fine < z_min[:, None]] = np.nan
        # Compute dz_fine
        dz_fine = np.abs(np.diff(z_fine, axis=1))  # Compute differences along the second axis
        # Insert the first element (z_fine[:, 0] - 0) at the beginning
        dz_fine = np.hstack([z_fine[:, [0]], dz_fine])
        z_mid_fine = dz_fine / 2 - np.cumsum(dz_fine, axis=1)
        ix = np.full((z_fine.shape), np.nan)
        # Expand z along the second axis to match z_fine's shape (broadcasting)
        z_expanded = np.expand_dims(z, axis=1)  # Shape: (rows, 1, cols)
        z_fine_expanded = np.expand_dims(z_fine, axis=2)  # Shape: (rows, fine_steps, 1)

    # Compute mask using broadcasting (row-wise comparison)
    mask = (z_fine_expanded < z_expanded) & ~np.isclose(z_fine_expanded, z_expanded, atol=1e-9)

    # Sum along the depth dimension to count how many times z_fine falls below z
    ix = np.sum(mask, axis=2).astype(np.float64)  # Convert to float to retain NaN compatibility

    pF_fine = {}

    for key in pF[0].keys():  # Iterate over each parameter in `pF`
        # Convert pF into an array ensuring consistent shapes
        try:
            pF_array = np.vstack([p[key] for p in pF])  # Ensures (rows, depths) shape
        except ValueError:  # If rows have different lengths, handle it gracefully
            max_depth = max(len(p[key]) for p in pF)  # Find the longest row
            pF_array = np.full((len(pF), max_depth), np.nan)  # Initialize padded array

            # Fill rows with actual values
            for i, p in enumerate(pF):
                pF_array[i, :len(p[key])] = p[key]

        # Ensure `ix` values are within valid range (clip to prevent indexing errors)
        ix_valid = np.clip(ix.astype(int), 0, pF_array.shape[1] - 1)

        # Assign values using vectorized indexing
        pF_fine[key] = np.take_along_axis(pF_array, ix_valid, axis=1)  # Shape: (rows, fine_steps)

    # --------- connection between gwl and Wsto, Tr, C------------
    if isinstance(grid_step, float):
        gwl = np.arange(1.0, min(z_min)-5, grid_step)
    elif grid_step == 'var':
        limits = [-1, -5, min(z_min)-5]
        steps = [-0.01, -0.05, -0.5]
        # First segment: 1m to 0m (step = 1.0m)
        z1 = np.arange(1, limits[0], steps[0])
        z2 = np.arange(limits[0], limits[1], steps[1])
        z3 = np.arange(limits[1], limits[2], steps[2])
        gwl = np.concatenate([z1, z2, z3])  # Ensure exact z_min

    Wsto_deep = np.stack([h_to_cellmoist_vectorized(pF_fine, g - z_mid_fine, dz_fine) + max(0.0, g) for g in gwl]).T

    if root:
        Wsto_deep = np.stack([h_to_cellmoist_vectorized(pF_fine, g - z_mid_fine, dz_fine) for g in gwl]).T
        Wsto_deep = Wsto_deep/np.nansum(dz, axis=1)
        #GwlToWsto = interp1d(np.array(gwl), np.array(Wsto_deep), fill_value='extrapolate')
        GwlToWsto = [interp1d(gwl, wsto_row, kind='linear', fill_value='extrapolate') for wsto_row in Wsto_deep]
        return {'to_rootmoist': GwlToWsto}

    Tr1 = np.stack([transmissivity_vectorized(dz, Ksat, g) * 86400. for g in gwl]).T

    # Generate interpolators for each row of Wsto_deep and Tr1 while keeping gwl the same
    WstoToGwl = [interp1d(wsto_row, gwl, kind='linear', fill_value='extrapolate') for wsto_row in Wsto_deep]
    GwlToWsto = [interp1d(gwl, wsto_row, kind='linear', fill_value='extrapolate') for wsto_row in Wsto_deep]
    GwlToC = [interp1d(gwl, np.gradient(wsto_row) / np.gradient(gwl), kind='linear', fill_value='extrapolate') for wsto_row in Wsto_deep]
    GwlToTr = [interp1d(gwl, tr_row, kind='linear', fill_value='extrapolate') for tr_row in Tr1]
    
    #plt.figure(1)
    #plt.plot(np.array(gwl), np.array(np.gradient(Wsto_deep[0])/np.gradient(gwl)), linestyle='--')
    #plt.figure(2)
    #plt.plot(np.array(gwl), np.log10(np.array(Tr1[0])), linestyle='--')
    #plt.plot(np.array(gwl), np.array(Tr1[0]), linestyle='--')
    #plt.figure(3)
    #plt.plot(np.array(gwl), np.array(Wsto_deep[0]), linestyle='--')

    return {'to_gwl': WstoToGwl, 'to_wsto': GwlToWsto, 'to_C': GwlToC, 'to_Tr': GwlToTr}

def h_to_cellmoist_vectorized(pF, h, dz):
    r""" Cell moisture based on vanGenuchten-Mualem soil water retention model.
    Partly saturated cells calculated as thickness weigthed average of
    saturated and unsaturated parts.

    Args:
        pF (np.ndarray):
            dict
                'ThetaS' (np.ndarray): saturated water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'ThetaR' (np.ndarray): residual water content [m\ :sup:`3` m\ :sup:`-3`\ ]
                'alpha' (np.ndarray): air entry suction [cm\ :sup:`-1`]
                'n' (np.ndarray): pore size distribution [-]
        h (float): pressure head [m]
        dz (np.ndarray): soil compartment thichness, node in center [m]
    Returns:
        theta (np.ndarray): Total volumetric water content of cell for given gwl
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
    ix = np.where(abs(h[0]) < dz/2)

    if Ts.shape[1] == 1:
        ixx = (np.array([0]), np.array([0]))  # Single index for 2D case
    else:
        ixx = ix
    # moisture of unsaturated part
    x[ix] = -(dz[ix]/2 - h[ix]) / 2
    theta[ix] = Tr[ixx] + (Ts[ixx] - Tr[ixx]) / (1 + abs(alfa[ixx] * 100 * x[ix])**n[ixx])**m[ixx]
    # total moisture as weighted average
    theta[ix] = (theta[ix] * (dz[ix]/2 - h[ix]) + Ts[ixx] * (dz[ix]/2 + h[ix])) / (dz[ix])
    # from vwc to total water content
    Wsto = theta * dz
    # 
    Wsto = np.nansum(Wsto, axis=1)
    
    return Wsto

def transmissivity_vectorized(dz, Ksat, gwl):
    r""" Vectorized transmissivity function for 2D inputs.

    Args:
       dz (np.ndarray):  Soil compartment thickness, node in center [m]
       Ksat (np.ndarray): Horizontal saturated hydraulic conductivity [m/s]
       gwl (float): Groundwater level below surface, <0 [m]

    Returns:
       Tr (np.ndarray): Transmissivity for each cell [m²/s]
    """

    # Compute midpoints of layers
    z = dz / 2 - np.cumsum(dz, axis=1)  # Shape: (n_cells, n_layers)

    # Compute saturated thickness for each layer, between [0, dz]
    dz_sat = np.minimum(np.maximum(gwl - (z - dz / 2), 0), dz)  # Shape: (n_cells, n_layers)
    # In top cell allow transmissivity to increases when gwl above ground surface level
    dz_sat[:, 0] = np.maximum(gwl - (z[:, 0] - dz[:, 0] / 2), 0) 

    # Compute transmissivity of each layer
    Trans = Ksat * dz_sat  # Shape: (n_cells, n_layers)

    # Sum transmissivity across layers
    Tr = np.maximum(np.sum(Trans, axis=1), 1e-4 / 86400)  # Shape: (n_cells,)

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
