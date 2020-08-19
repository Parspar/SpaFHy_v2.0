# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:59:54 2019
@author: alauren

Modified by khaahti
"""

import numpy as np
from scipy.stats.mstats import gmean

class SoilGrid_2Dflow(object):
    """
    2D soil water flow model based on Ari Lauren SUSI2D
    Simulates moss/organic layer with interception and evaporation,   !!! MOSS
    soil water storage, drainage to ditches and pond storage on top of soil.
    """
    def __init__(self, spara):
        """
        Initializes SoilProfile2D:
        Args:
            spara (dict):
                'elevation'
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
                'pond_storage_max': maximum pond depth [m]
                # initial states
                'ground_water_level': groundwater depth [m]
                'org_sat': organic top layer saturation ratio (-)
                'pond_storage': initial pond depth at surface [m]
        """
        # for testing
        self.maxEle = 0.0

        # parameters
        self.implic = 1  # solving method: 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
        self.DrIrr = False  # if true: constant head boundaries, if false changing boundary: gwl above canal, using constant head, if below below no flow

        # initial h (= gwl) and boundaries [m]
        self.h = spara['ground_water_level']
        h0 = -0.9  # depth of water level in ditch / canal, from input?!
        self.h[0,:] = h0
        self.h[-1,:] = h0
        self.h[:,0] = h0
        self.h[:,-1] = h0

        # soil surface elevation and hydraulic head [m]
        self.ele = np.full_like(self.h, 0.0)  # sould come from input?!
        self.H = self.ele + self.h

        # soil/peat type
        self.soiltype = spara['soiltype']

        # interpolated functions for soil column ground water dpeth vs. water storage, transmissivity etc.
        self.wsto_to_gwl = spara['wtso_to_gwl']
        self.gwl_to_wsto = spara['gwl_to_wsto']
        self.gwl_to_Tr = spara['gwl_to_Tr']
        self.gwl_to_C = spara['gwl_to_C']

        # water storage [m]
        self.Wsto_max = np.full_like(self.h, 0.0)  # storage of fully saturated profile
        for key, value in self.gwl_to_wsto.items():
            self.Wsto_max[self.soiltype == key] = value(0.0)
        self.Wsto = np.full_like(self.h, 0.0)  # storage corresponding to h
        for key, value in self.gwl_to_wsto.items():
            self.Wsto[self.soiltype == key] = value(self.h[self.soiltype == key])

        # transmissivity [m2 d-1]
        self.Tr0 = np.full_like(self.h, 0.0)
        for key, value in self.gwl_to_Tr.items():
            self.Tr0[self.soiltype == key] = value(self.h[self.soiltype == key])

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
        self.Tr1 = np.zeros((self.rows,self.cols))

        # needed for canopy
        self.Rew = 1.0
        self.Ree = 1.0

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
        # runoff components
        runo = 0.
        sruno = 0.

        # for testing
        h0ts = -0.9

        # source/sink [m s-1]
        S = rr - tr

        # hydrauli heads, new iteration and old iteration
        Htmp = self.H.copy()
        Htmp1 = self.H.copy()

        # air volume [m]
        airv = self.Wsto_max - self.Wsto
        # remove water that cannot fit into soil, add it to surface runoff
        S = np.where(S > airv, airv, S)
        sruno += np.where(S > airv, S - airv, 0.0)

        # needed here?
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

        # ravel 2D arrays
        self.TrW0 = np.ravel(self.TrW0)
        self.TrE0 = np.ravel(self.TrE0)
        self.TrN0 = np.ravel(self.TrN0)
        self.TrS0 = np.ravel(self.TrS0)

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:-1]
        self.HE[:,:-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:-1,:]
        self.HS[:-1,:] = self.H[1:,:]

        # ravel 2D arrays
        self.HW = np.ravel(self.HW)
        self.HE = np.ravel(self.HE)
        self.HN = np.ravel(self.HN)
        self.HS = np.ravel(self.HS)

        conv0 = 10000  # error
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
            self.TrW1 = np.ravel(self.TrW1)
            self.TrE1= np.ravel(self.TrE1)
            self.TrN1 = np.ravel(self.TrN1)
            self.TrS1 = np.ravel(self.TrS1)
            del TrTmpEW, TrTmpNS

            # differential water capacity dSto/dh
            for key, value in self.gwl_to_C.items():
                self.CC[self.soiltype == key] = value(Htmp1[self.soiltype == key] - self.ele[self.soiltype == key])
            alfa = np.ravel(self.CC * self.dxy**2 / dt)

            # Penta-diagonal matrix setup
            i,j = np.indices(self.A.shape)
            self.A[i==j] = self.implic * (self.TrW1 + self.TrE1 + self.TrN1 + self.TrS1) + alfa  # Diagonal
            self.A[i==j+1] = -self.implic * self.TrW1[1:]  # West element
            self.A[i==j-1] = -self.implic * self.TrE1[:-1]  # East element
            self.A[i==j+self.cols] = -self.implic * self.TrN1[self.cols:]  # North element
            self.A[i==j-self.cols] = -self.implic * self.TrS1[:self.n-self.cols]  # South element

            # Boundaries
            # North
            for k in range(0,self.cols):
                self.A[k,k] = self.implic * (2 * self.TrS1[k]) + alfa[k]
                self.A[k,k+1] = 0
                self.A[k,k-1] = 0
                self.A[k,k+self.cols] = -self.implic * 2 * self.TrS1[k]
            # South
            for k in range(self.n-self.cols,self.n-1):
                self.A[k,k] = self.implic * (2 * self.TrN1[k]) + alfa[k]
                self.A[k,k+1] = 0
                self.A[k, k-1] = 0
                self.A[k, k-self.cols] = -self.implic * 2 * self.TrN1[k]
            self.A[self.n-1,self.n-1] = self.implic * (2 * self.TrN1[self.n-1]) + alfa[self.n-1]
            self.A[self.n-1,self.n-2] = 0
            self.A[self.n-1,self.n-1-self.cols] = -self.implic * 2 * self.TrN1[self.n-1]
            # West
            for k in np.arange(self.cols, self.n-self.cols, self.cols, dtype='int'):
                self.A[k,k] = self.implic * (2 * self.TrE1[k]) + alfa[k]
                self.A[k, k-self.cols] = 0
                self.A[k, k+self.cols] = 0
                self.A[k,k+1] = -self.implic * 2 * self.TrE1[k]
            # East
            for k in np.arange(self.cols-1, self.n-self.cols, self.cols, dtype='int'):
                self.A[k,k]= self.implic*(2*self.TrW1[k])+alfa[k]
                if k-self.cols > 0:
                    self.A[k, k-self.cols]=0
                    self.A[k,k-1]=-self.implic*2*self.TrW1[k]
                    self.A[k, k+self.cols]=0

            # Knowns: Right hand side of the eq
            S = np.ravel(S)
            self.H = np.ravel(self.H)
            hs = (S * dt * self.dxy**2 + alfa*self.H
                  + (1.-self.implic) * (self.TrN0*self.HN) + (1.-self.implic)*(self.TrW0*self.HW)
                  - (1.-self.implic)*(self.TrN0 + self.TrW0 + self.TrE0 + self.TrS0)*self.H
                  + (1.-self.implic)*(self.TrE0*self.HE) + (1.-self.implic)*(self.TrS0*self.HS))

            # Solve: A*Htmp1 = hs
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(self.A),hs])
            Htmp1 = np.reshape(Htmp1,(self.rows,self.cols))

            # Return boundaries ?!
            Htmp1[0,:] = np.ones(self.cols)*self.maxEle+h0ts
            Htmp1[self.rows-1]=np.ones(self.cols)*self.maxEle+h0ts
            Htmp1[:,0]=np.ones(self.rows)*self.maxEle+h0ts
            Htmp1[:,self.cols-1]=np.ones(self.rows)*self.maxEle+h0ts

            self.TrW1 = np.reshape(self.TrW1,(self.rows,self.cols))
            self.TrE1 = np.reshape(self.TrE1,(self.rows,self.cols))
            self.TrN1 = np.reshape(self.TrN1,(self.rows,self.cols))
            self.TrS1 = np.reshape(self.TrS1,(self.rows,self.cols))

            self.H = np.reshape(self.H,(self.rows,self.cols))

            alfa = np.reshape(alfa,(self.rows,self.cols))
            Htmp1 = np.where(Htmp1>self.ele, self.ele,Htmp1)
            conv1 = np.max(np.abs(Htmp1 - Htmp))
            Htmp = Htmp1.copy()
            if conv1 < crit:
                print('iterations', it, np.average(self.ele[1:-1,1:-1]-Htmp[1:-1,1:-1]))
                break
            # end of iteration loop

        self.TrW0 = np.reshape(self.TrW0,(self.rows,self.cols))
        self.TrE0 = np.reshape(self.TrE0,(self.rows,self.cols))
        self.TrN0 = np.reshape(self.TrN0,(self.rows,self.cols))
        self.TrS0 = np.reshape(self.TrS0,(self.rows,self.cols))
        self.HW = np.reshape(self.HW,(self.rows,self.cols))
        self.HE = np.reshape(self.HE,(self.rows,self.cols))
        self.HN = np.reshape(self.HN,(self.rows,self.cols))
        self.HS = np.reshape(self.HS,(self.rows,self.cols))

        if self.DrIrr==True:
            Htmp1[0,:]=np.ones(self.cols)*self.maxEle+h0ts
            Htmp1[self.rows-1,:]=np.ones(self.cols)*self.maxEle+h0ts
            Htmp1[:,0]=np.ones(self.rows)*self.maxEle+h0ts
            Htmp1[:,self.cols-1]=np.ones(self.rows)*self.maxEle+h0ts
        else:
            Htmp1[0,:]=np.where(Htmp1[1,:]<np.ones(self.cols)*self.maxEle+h0ts,Htmp1[1,:],np.ones(self.cols)*self.maxEle+h0ts)
            Htmp1[self.rows-1,:]=np.where(Htmp1[self.rows-2,:]<np.ones(self.cols)*self.maxEle+h0ts,Htmp1[self.rows-2,:],np.ones(self.cols)*self.maxEle+h0ts)
            Htmp1[:,0]=np.where(Htmp1[:,1]<np.ones(self.rows)*self.maxEle+h0ts, Htmp1[:,1],np.ones(self.rows)*self.maxEle+h0ts)
            Htmp1[:,self.cols-1]=np.where(Htmp1[:,self.cols-2]<np.ones(self.rows)*self.maxEle+h0ts, Htmp1[:,self.cols-2],np.ones(self.rows)*self.maxEle+h0ts)

        # update state
        self.H = Htmp1.copy()
        for key, value in self.gwl_to_wsto.items():
            self.Wsto[self.soiltype == key] = value(Htmp1[self.soiltype == key] - self.ele[self.soiltype == key])

        results = {
                'pond_storage': 0.0 * self.H,  # [m]
                'ground_water_level': self.H - self.ele,  # [m]
                'infiltration': 0.0 * self.H,  # [mm d-1]
                'surface_runoff': 0.0 * self.H,  # [mm d-1]
                'evaporation': 0.0 * self.H,  # [mm d-1]
                'drainage': 0.0 * self.H,  # [mm d-1]
                'moisture_top': 0.0 * self.H,  # [m3 m-3]
                'water_closure': 0.0 * self.H,  # [mm d-1]
                'transpiration_limitation': 0.0 * self.H,  # [-]
                'rootzone_moisture': 0.0 * self.H,  # [m3 m-3]
                }

        return results