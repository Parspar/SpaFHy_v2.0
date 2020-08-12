# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:59:54 2019
@author: alauren

Modified by khaahti
"""

import numpy as np
from scipy.stats.mstats import gmean

from susi_utils import peat_hydrol_properties, CWTr, wrc

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

        self.h = spara['ground_water_level']
        self.ele = np.full_like(self.h, 0.0)
        self.rows = np.shape(self.h)[0]
        self.cols = np.shape(self.h)[1]
        self.dxy = spara['dxy']  # horizontal distance between nodes dx=dy [m]

        # print('Initializing 2D peat strip')
        # nLyrs = spara['nLyrs']                                                      # number of soil layers
        # dz = np.ones(nLyrs)*spara['dzLyr']                                          # thickness of layers, m
        # z = np.cumsum(dz)-dz/2.                                                     # depth of the layer center point, m
        # lenvp=len(spara['vonP top'])
        # vonP = np.ones(nLyrs)*spara['vonP bottom']; vonP[0:lenvp] = spara['vonP top']  # degree of  decomposition, von Post scale
        # ptype = spara['peat type bottom']*spara['nLyrs']
        # lenpt = len(spara['peat type']); ptype[0:lenpt] = spara['peat type']
        # self.pF, self.Ksat = peat_hydrol_properties(vonP, var='H', ptype=ptype)  # peat hydraulic properties after Päivänen 1973

        # for n in range(nLyrs):
        #     if z[n] < 0.31:
        #         self.Ksat[n]= self.Ksat[n]*spara['anisotropy']
        #     else:
        #         self.Ksat[n]= self.Ksat[n]*1.

        # self.hToSto, self.stoToGwl, self.hToTra, self.C, self.hToRat, self.hToAfp = CWTr(nLyrs, z, dz, self.pF,
        #                                                        self.Ksat, direction='negative') # interpolated storage, transmissivity, diff water capacity, and ratio between aifilled porosoty in rooting zone to total airf porosity  functions

              #**************Strip parameters****************
        self.implic = 1                                                              #Solving method: 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
        self.DrIrr=False                                                             #if true: constant head boundaries, if false changing boundary: gwl above canal, using constant head, if below below no flow
        self.dt= 1.                                                                   #time step, in days
        self.h0 =  -0.9                                                              #Initial conditions, gw at canal
        # self.hini = -0.2                                                            #initial h (gwl) in the compartment
        # self.dxy = spara['dxy']                                                               #square space element dx=dy in m
        # self.rows = spara['rows']
        # self.cols = spara['cols']
        self.n = self.rows*self.cols                                                 #extension of the strip
        print('    +Compartment length, width: ',  self.cols*self.dxy, ' m', self.rows*self.dxy,'m, ')

        #*********** Initial gwl********************
        # self.h=np.ones((self.rows,self.cols),float)*self.hini
        self.h[0]=np.ones(self.cols)*self.h0
        self.h[self.rows-1]=np.ones(self.cols)*self.h0
        self.h[:,0]=np.ones(self.rows)*self.h0
        self.h[:,self.cols-1]=np.ones(self.rows)*self.h0

        #************ Surface elevation with slopes x,y -direction*******************
        lev=1. ; slox=0.0; sloy=0.0                                        # ref point m, slope in x direction %, in y direction %
        y = np.linspace(0,self.dxy*self.rows*sloy/100., self.rows)                            # surface rise in y direction, m
        x = np.linspace(0,self.dxy*self.cols*slox/100., self.cols)                                 # surface rise in x direction, m
        xv,yv=np.meshgrid(x,y)                                                 # creating the surface change in x,y
        # self.ele --> todellinen korkeus
        self.ele= (xv+yv)+lev                                                       # creating the surface elevation
        #print 'Surface elevation in m, with respect to reference level'
        #print ele
        self.maxEle =np.max(self.ele)

        #************create arrays only once
        self.TrW0 = np.zeros((self.rows,self.cols)); self.TrE0 = np.zeros((self.rows,self.cols))          #Transmissivity West East, previous time step
        self.TrN0 = np.zeros((self.rows,self.cols)); self.TrS0 = np.zeros((self.rows,self.cols))          #Transmissivity North South, previous time step
        self.HW = np.zeros((self.rows,self.cols)); self.HE = np.zeros((self.rows,self.cols))              #Head West East, previous time step
        self.HN = np.zeros((self.rows,self.cols)); self.HS = np.zeros((self.rows,self.cols))              #Head North South, previous time step
        self.TrW1 = np.zeros((self.rows,self.cols)); self.TrE1 = np.zeros((self.rows,self.cols))          #Transmissivity West East, current time step
        self.TrN1 = np.zeros((self.rows,self.cols)); self.TrS1 = np.zeros((self.rows,self.cols))          #Transmissivity North South, current time step


        print('    +Max surface elevation', self.maxEle)
        print('    +Water level in canal ', self.h0, 'm -> In reference height H: ', self.maxEle+self.h0, ' m')
        print('2D peatland strip initialized')

    def reset_domain(self, length):
        self.A=np.zeros((self.n,self.n))                                       # computation matrix
        self.H=self.ele+self.hini
        self.h = np.ones(np.shape(self.h))*self.hini                                     # right hand side vector
        self.H=self.ele+self.h                                                 # head with respect to absolute reference level, m
        self.sruno=0.
        self.roff = 0.
        print('Resetting strip scenario')

    def rolling_window(self, a, window):
        shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
        strides = a.strides + (a.strides[-1],)
        return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

    def run_timestep(self,d,h0ts, p, moss):

        """
        IN:
            d day number
            h0ts boudary (ditch depth, m) in time series
            p rainfall-et m, arrayn n length
            moss as object
            all arrays flattened
        OUT:
            all arrays flattened
        """

        # update hydraulic head in ditch cells
        self.H[0]=np.ones(self.cols)*self.maxEle+h0ts
        self.H[self.rows-1]=np.ones(self.cols)*self.maxEle+h0ts
        self.H[:,0]=np.ones(self.rows)*self.maxEle+h0ts
        self.H[:,self.cols-1]=np.ones(self.rows)*self.maxEle+h0ts

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
        sruno=0.

        # source/sink [m]
        S = np.reshape(p, (self.rows, self.cols))

        # hydrauli heads, new iteration and old iteration
        Htmp = self.H.copy()
        Htmp1 = self.H.copy()

        # air volume [m]
        airv = self.hToSto(self.ele) - self.hToSto(Htmp-self.ele)
        # remove water that cannot fit into soil, add it to surface runoff
        S = np.where(S > airv, airv, S)
        sruno += np.where(S > airv, S - airv, 0.0)

        # transmissivity to West, East, North, South cell
        Tr0 = self.hToTra(self.H - self.ele)  # mid cell
        # transmissivity at all four sides of the element is computed as geometric mean of surrounding element transimissivities
        TrTmpEW =gmean(self.rolling_window(Tr0,2),-1)
        TrTmpNS = np.transpose(gmean(self.rolling_window(np.transpose(Tr0),2),-1))

        self.TrW0[:,1:] = TrTmpEW
        self.TrE0[:,:self.cols-1] = TrTmpEW
        self.TrN0[1:,:] = TrTmpNS
        self.TrS0[:self.rows-1,:] = TrTmpNS
        del TrTmpEW, TrTmpNS

        # ravel 2D arrays
        self.TrW0=np.ravel(self.TrW0)
        self.TrE0=np.ravel(self.TrE0)
        self.TrN0=np.ravel(self.TrN0)
        self.TrS0=np.ravel(self.TrS0)

        # Head in four neighbouring cells
        self.HW[:,1:] = self.H[:,:self.cols-1]
        self.HE[:,:self.cols-1] = self.H[:,1:]
        self.HN[1:,:] = self.H[:self.rows-1,:]
        self.HS[:self.rows-1,:] = self.H[1:,:]

        # ravel 2D arrays
        self.HW = np.ravel(self.HW)
        self.HE=np.ravel(self.HE)
        self.HN=np.ravel(self.HN)
        self.HS=np.ravel(self.HS)

        conv0 = 10000  # error
        crit = 1e-6  # convergence criteria
        maxiter = 100

        for it in range(maxiter):
            CC = self.C(Htmp1 - self.ele)  # differential water capacity dSto/dh

            # transmissivity to neighbouring cells with HTmp1
            Tr1 = self.hToTra(Htmp1 - self.ele)
            TrTmpEW =gmean(self.rolling_window(Tr1,2),-1)
            TrTmpNS = np.transpose(gmean(self.rolling_window(np.transpose(Tr1),2),-1))
            self.TrW1[:,1:] = TrTmpEW; self.TrE1[:,:self.cols-1] = TrTmpEW
            self.TrN1[1:,:] = TrTmpNS; self.TrS1[:self.rows-1,:]=TrTmpNS
            self.TrW1 = np.ravel(self.TrW1); self.TrE1= np.ravel(self.TrE1)
            self.TrN1 = np.ravel(self.TrN1); self.TrS1 = np.ravel(self.TrS1)
            del TrTmpEW, TrTmpNS

            alfa = np.ravel(CC * self.dxy**2 / self.dt)

            # Penta-diagonal matrix setup
            i,j = np.indices(self.A.shape)
            self.A[i==j]= self.implic*(self.TrW1+self.TrE1+self.TrN1+self.TrS1) +alfa                    # Diagonal
            self.A[i==j+1]= -self.implic*self.TrW1[1:]                                    # West element
            self.A[i==j-1]=-self.implic*self.TrE1[:self.n-1]                                   # East element
            self.A[i==j+self.cols]=-self.implic*self.TrN1[self.cols:]                               # North element
            self.A[i==j-self.cols]=-self.implic*self.TrS1[:self.n-self.cols]                             # South element
            #****Boundaries
            #*****Set north canal*************
            for k in range(0,self.cols):                                         #North canal
                self.A[k,k]=self.implic*(2*self.TrS1[k])+alfa[k]
                self.A[k,k+1]=0
                self.A[k, k-1]=0
                self.A[k,k+self.cols]=-self.implic*2*self.TrS1[k]
            #******Set south canal***************
            for k in range(self.n-self.cols,self.n-1):                                         #South canal
                self.A[k,k]=self.implic*(2*self.TrN1[k])+alfa[k]
                self.A[k,k+1]=0
                self.A[k, k-1]=0
                self.A[k, k-self.cols]=-self.implic*2*self.TrN1[k]

            self.A[self.n-1,self.n-1] = self.implic*(2*self.TrN1[self.n-1])+alfa[self.n-1]
            self.A[self.n-1,self.n-2]= 0
            self.A[self.n-1,self.n-1-self.cols]=-self.implic*2*self.TrN1[self.n-1]
            #********set west canal, first element
            for k in np.arange(self.cols,self.n-self.cols,self.cols,dtype='int'):
                self.A[k,k]= self.implic*(2*self.TrE1[k])+alfa[k]
                self.A[k, k-self.cols]=0
                self.A[k, k+self.cols]=0
                self.A[k,k+1]=-self.implic*2*self.TrE1[k]
            #********set east canal, last element*********************
            for k in np.arange(self.cols-1,self.n-self.cols,self.cols,dtype='int'):
                self.A[k,k]= self.implic*(2*self.TrW1[k])+alfa[k]
                if k-self.cols > 0:
                    self.A[k, k-self.cols]=0
                    self.A[k,k-1]=-self.implic*2*self.TrW1[k]
                    self.A[k, k+self.cols]=0

            #***********Knowns: Right hand side of the eq
            S=np.ravel(S); self.H =np.ravel(self.H)
            hs=S*self.dt*self.dxy**2 + alfa*self.H + (1.-self.implic)*(self.TrN0*self.HN) + (1.-self.implic)*(self.TrW0*self.HW) \
            - (1.-self.implic)*(self.TrN0 + self.TrW0 + self.TrE0 + self.TrS0)*self.H  \
            + (1.-self.implic)*(self.TrE0*self.HE) + (1.-self.implic)*(self.TrS0*self.HS)

            #************ Solve************
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(self.A),hs])
            Htmp1=np.reshape(Htmp1,(self.rows,self.cols));
            #******** Return boundaries********
            Htmp1[0]=np.ones(self.cols)*self.maxEle+h0ts; Htmp1[self.rows-1]=np.ones(self.cols)*self.maxEle+h0ts
            Htmp1[:,0]=np.ones(self.rows)*self.maxEle+h0ts; Htmp1[:,self.cols-1]=np.ones(self.rows)*self.maxEle+h0ts

            self.TrW1=np.reshape(self.TrW1,(self.rows,self.cols))
            self.TrE1=np.reshape(self.TrE1,(self.rows,self.cols))
            self.TrN1=np.reshape(self.TrN1,(self.rows,self.cols))
            self.TrS1=np.reshape(self.TrS1,(self.rows,self.cols))
            S=np.reshape(S,(self.rows,self.cols))
            self.H=np.reshape(self.H,(self.rows,self.cols))
            alfa = np.reshape(alfa,(self.rows,self.cols))
            Htmp1=np.where(Htmp1>self.ele, self.ele,Htmp1)
            conv1= np.max(np.abs(Htmp1-Htmp))
            Htmp=Htmp1.copy()
            if conv1<crit:
                if d%30==0: print('  - day #',d, 'iterations', it, np.average(self.ele[1:-1,1:-1]-Htmp[1:-1,1:-1]))
                break
            conv0=conv1.copy()
            # end of iteration loop

        self.TrW0=np.reshape(self.TrW0,(self.rows,self.cols))
        self.TrE0=np.reshape(self.TrE0,(self.rows,self.cols))
        self.TrN0=np.reshape(self.TrN0,(self.rows,self.cols))
        self.TrS0=np.reshape(self.TrS0,(self.rows,self.cols))
        self.HW=np.reshape(self.HW,(self.rows,self.cols))
        self.HE=np.reshape(self.HE,(self.rows,self.cols))
        self.HN=np.reshape(self.HN,(self.rows,self.cols))
        self.HS=np.reshape(self.HS,(self.rows,self.cols))
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

        self.H=Htmp1.copy()

        dwt = self.H-self.ele

        air_ratio = np.ones(np.shape(dwt))
        afp = np.ones(np.shape(dwt))

        return np.ravel(dwt), np.ravel(self.H), runo, np.ravel(air_ratio), np.ravel(afp)