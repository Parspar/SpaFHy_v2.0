# -*- coding: utf-8 -*-
"""
Gridded canopy and snow hydrology model for SpaFHy.

Computes water interception, evapotranspiration, and snowpack dynamics
at daily or sub-daily timesteps for gridded (2D array) applications.

Main components:
  - Phenology: seasonal LAI dynamics for deciduous species and grasses,
    driven by degree-day accumulation and daylength.
  - Canopy water interception and evaporation from canopy storage.
  - Snowpack: accumulation, melt, freezing, and liquid water retention.
  - Dry-canopy evapotranspiration: transpiration via Penman-Monteith with
    upscaled stomatal conductance, and forest floor evaporation.
  - Aerodynamic resistances: logarithmic above-canopy and exponential
    within-canopy wind profiles.

References:
    Launiainen et al. (2019). Hydrol. Earth Syst. Sci., 23, 3457-3480.
    Nousu et al. (2024). Hydrol. Earth Syst. Sci., 28, 4643-4666.
    Launiainen et al. (2016). Global Change Biol., doi:10.1111/gcb.13226.
    Leuning et al. (2008). Water Resour. Res., 44, W10419.
    Peltoniemi et al. (2015). Boreal Env. Res., 20, 1-20.
    Hedstrom & Pomeroy (1998). Hydrol. Process., 12, 1611-1625.

@authors: slauniainen, khaahti, jpnousu
"""

import numpy as np
import configparser
from netCDF4 import Dataset
eps = np.finfo(float).eps

class CanopyGrid():
    def __init__(self, cpara, state, dist_rad_file=None):
        """
        Initializes CanopyGrid object.

        Args:
            cpara (dict): Canopy parameter dictionary. Expected top-level keys:
                'loc'      - dict with 'lat' [deg] and 'lon' [deg]
                'physpara' - dict of physiological parameters (amax, g1_conif,
                            g1_decid, g1_grass, g1_shrub, kp, q50, rw, rwmin, gsoil)
                'phenopara'- dict of phenology parameters (LAI_decid_min, ddo, ddur,
                            sdl, sdur, tau, xo, smax, fmin)
                'interc'   - dict with 'wmax' [mm LAI-1] and 'wmaxsnow' [mm LAI-1]
                'snow'     - dict with 'kmelt' [mm s-1 degC-1], 'kfreeze', 'r' [-]
                'flow'     - dict with 'zmeas' [m], 'zground' [m], 'zo_ground' [m]

            state (dict): Initial state arrays (all 2D grids). Expected keys:
                'LAI_conif'      [m2 m-2]  conifer leaf area index
                'LAI_decid'      [m2 m-2]  deciduous leaf area index (maximum)
                'LAI_grass'      [m2 m-2]  grass leaf area index (maximum)
                'LAI_shrub'      [m2 m-2]  shrub leaf area index
                'canopy_height'  [m]       mean canopy height
                'canopy_fraction'[-]       canopy cover fraction
                'w'              [mm]      initial canopy water storage
                'swe'            [mm]      initial snow water equivalent

            dist_rad_file (str, optional): Path to NetCDF file containing 
                spatially distributed radiation coefficients ('c_rad'). 
                If None, spatially uniform radiation is used.

        NOTE:
            Currently the initialization assumes simulation start 1st Jan,
            and sets self._LAI_decid and self.X equal to minimum values.
            Also leaf-growth & senescence parameters are intialized to zero.
        """

        self.cmask = np.full_like(state['LAI_conif'], np.nan)
        self.cmask[np.isfinite(state['LAI_conif'])] = 1.0

        self.latitude = cpara['loc']['lat'] * self.cmask
        self.longitude = cpara['loc']['lon'] * self.cmask

        # physiology: transpi + floor evap
        self.physpara = cpara['physpara']

        # phenology
        self.phenopara = cpara['phenopara']

        # canopy parameters and state
        self.hc = state['canopy_height'] + eps
        self.cf = state['canopy_fraction'] + eps
        self._LAIconif = np.maximum(state['LAI_conif'], eps)  # m2m-2
        self._LAIdecid = state['LAI_decid'] * self.phenopara['LAI_decid_min']
        self._LAIgrass_max = state['LAI_grass']
        self._LAIgrass = state['LAI_grass'] * self.phenopara['LAI_decid_min']
        self._LAIshrub = np.maximum(state['LAI_shrub'], eps)
        self.LAI = self._LAIconif + self._LAIdecid + self._LAIshrub + self._LAIgrass
        self._LAIdecid_max = state['LAI_decid']  # m2m-2

        # senescence starts at first doy when daylength < self.phenopara['sdl']
        self.phenopara['sso'] = np.ones(np.shape(self.latitude))*np.nan
        doy = np.arange(1, 366)
        for lat in np.unique(self.latitude):
            if np.isnan(lat):
                break
            # senescence starts at first doy when daylength < self.phenopara['sdl']
            dl = daylength(lat, doy)
            ix = np.max(np.where(dl > self.phenopara['sdl']))
            self.phenopara['sso'][self.latitude == lat] = doy[ix]  # this is onset date for senescence
            del ix
            
        self.phenopara['sso'] = self.phenopara['sso'] * self.cmask

        self.wmax = cpara['interc']['wmax']
        self.wmaxsnow = cpara['interc']['wmaxsnow']
        self.Kmelt = cpara['snow']['kmelt']
        self.Kfreeze = cpara['snow']['kfreeze']
        self.R = cpara['snow']['r']  # max fraction of liquid water in snow

        # --- for computing aerodynamic resistances
        self.zmeas = cpara['flow']['zmeas']
        self.zground = cpara['flow']['zground'] # reference height above ground [m]
        self.zo_ground = cpara['flow']['zo_ground'] # ground roughness length [m]
        self.gsoil = self.physpara['gsoil']

        # --- state variables
        self.W = np.minimum(state['w'], self.wmax*self.LAI)
        self.SWE = state['swe']
        self.SWEi = self.SWE
        self.SWEl = np.zeros(np.shape(self.SWE))

        # deciduous leaf growth stage
        # NOTE: this assumes simulations start 1st Jan each year !!!
        self.DDsum = self.W * 0.0
        self.X = self.W * 0.0
        self._growth_stage = self.W * 0.0
        self._senesc_stage = self.W *0.0

        if dist_rad_file:
            self.rad_coeff = Dataset(dist_rad_file, 'r')
            self.distributed_radiation = True
            print('*** Distributed radiation used ***')
        else:
            self.distributed_radiation = False

    def run_timestep(self, doy, dt, Ta, Prec, Rg, Par, VPD, U=2.0, CO2=380.0, Rew=1.0, beta=1.0, P=101300.0):
        """
        Runs CanopyGrid instance for one timestep
        Args:
            doy  (int):              day of year [1-366]
            dt   (float):            timestep duration [s]
            Ta   (array or float):   air temperature [degC]
            Prec (array or float):   precipitation amount [mm per timestep]
            Rg   (array or float):   global radiation [W m-2]
            Par  (array or float):   photosynthetically active radiation [W m-2]
            VPD  (array or float):   vapor pressure deficit [kPa]
            U    (array or float):   wind speed at reference height above canopy [m s-1]
            CO2  (float):            atmospheric CO2 mixing ratio [ppm]
            Rew  (array or float):   relative extractable water in root zone [-]
            beta (array or float):   relative soil conductance for floor evaporation [-]
            P    (float):            ambient air pressure [Pa]

        Returns:
            dict with keys:
                'potential_infiltration'  [mm]: water reaching the soil surface
                'interception'            [mm]: water intercepted by canopy
                'evaporation'             [mm]: evaporation/sublimation from canopy storage
                'forestfloor_evaporation' [mm]: evaporation from forest floor
                'transpiration'           [mm]: canopy transpiration
                'throughfall'             [mm]: precipitation reaching field layer or snowpack
                'snow_water_equivalent'   [mm]: total SWE (ice + liquid)
                'water_closure'           [mm]: mass balance error (should be ~0)
                'phenostate'              [-]:  phenology modifier (0=dormant, 1=full leaf)
                'leaf_area_index'     [m2 m-2]: total one-sided LAI
                'stomatal_conductance'  [m s-1]: upscaled canopy conductance
                'degree_day_sum'        [degC]: cumulative degree-days since 1 Jan
                'water_storage'           [mm]: canopy water storage
                'snowfall'                [mm]: snowfall amount this timestep
                'rainfall'                [mm]: rainfall amount this timestep
        """

        flux_to_mm_d = 86400.0 / dt

        if self.distributed_radiation == True:
            Rg = update_distributed_radiation(rad_coeff=self.rad_coeff, doy=doy, Rg=Rg)
        # Rn = 0.7 * Rg #net radiation
        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
                            0.55) * Rg  # Launiainen et al. 2016 GCB, fit to Fig 2a
        
        """ --- update phenology: self.ddsum & self.X ---"""
        self._degreeDays(Ta, doy)
        fPheno = self._photoacclim(Ta)

        """ --- update deciduous leaf area index --- """
        laifract = self._lai_dynamics(doy)
        
        """ --- aerodynamic conductances --- """
        Ra, _, Ras, _, _, _ = aerodynamics(self.LAI, self.hc, U, w=0.01, zm=self.zmeas,
                                                  zg=self.zground, zos=self.zo_ground)

        """ --- interception, evaporation and snowpack --- """
        PotInf, Trfall, Evap, Interc, MBE, erate, unload, fact, Sfall, Rfall = self.canopy_water_snow(dt, Ta, Prec, Rn, VPD, Ra=Ra)

        """--- dry-canopy evapotranspiration [mm s-1] --- """
        Transpi, Efloor, Gc = self.dry_canopy_et(VPD, Par, Rn, Ta, Ra=Ra, Ras=Ras, CO2=CO2, Rew=Rew, beta=beta, fPheno=fPheno)

        Transpi = Transpi * dt
        Efloor = Efloor * dt
        #ET = Transpi + Efloor
        
        results = {
            # fluxes: scaled to mm d-1 regardless of timestep
            'potential_infiltration':  PotInf  * flux_to_mm_d,  # [mm d-1]
            'interception':            Interc  * flux_to_mm_d,  # [mm d-1]
            'evaporation':             Evap    * flux_to_mm_d,  # [mm d-1]
            'forestfloor_evaporation': Efloor  * flux_to_mm_d,  # [mm d-1]
            'transpiration':           Transpi * flux_to_mm_d,  # [mm d-1]
            'throughfall':             Trfall  * flux_to_mm_d,  # [mm d-1]
            'water_closure':           MBE     * flux_to_mm_d,  # [mm d-1]
            'snowfall':                Sfall   * flux_to_mm_d,  # [mm d-1]
            'rainfall':                Rfall   * flux_to_mm_d,  # [mm d-1]

            # states: instantaneous values, not scaled by dt
            'snow_water_equivalent': self.SWE,   # [mm]
            'water_storage':         self.W,     # [mm]
            'leaf_area_index':       self.LAI,   # [m2 m-2]
            'stomatal_conductance':  Gc,         # [m s-1]
            'phenostate':            fPheno,     # [-]
            'degree_day_sum':        self.DDsum, # [degC]
        }

        return results
    

    def _degreeDays(self, T, doy):
        """
        Updates the cumulative degree-day sum (DDsum) for the current timestep.
        Accumulation starts above a threshold temperature of 5 degC.
        DDsum is reset to zero at the start of each new year (doy == 1).

        Args:
            T   (array or float): daily mean air temperature [degC]
            doy (int):            day of year [1-366]
        """
        To = 5.0  # threshold temperature
        self.DDsum = self.DDsum + np.maximum(0.0, T - To)

        #reset at beginning of year
        self.DDsum[doy * self.cmask == 1] = 0.

    def _photoacclim(self, T):
        """
        Updates the temperature acclimation state and computes the phenology
        modifier for stomatal conductance. Follows Peltoniemi et al. (2015)
        Boreal Env. Res., 20, 1-20.

        The acclimation state X tracks a lagged temperature signal. The
        acclimation signal S is the excess of X above a threshold (xo).
        fPheno scales from fmin (winter minimum) to 1 (full acclimation).

        Args:
            T (array or float): daily mean air temperature [degC]

        Returns:
            fPheno (array): phenology modifier [-], range [fmin, 1.0]
        """

        self.X = self.X + 1.0 / self.phenopara['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.phenopara['xo'], 0.0)
        fPheno = np.maximum(self.phenopara['fmin'],
                            np.minimum(S / self.phenopara['smax'], 1.0))
        return fPheno

    def _lai_dynamics(self, doy):
        """
        Seasonal cycle of deciduous leaf area

        Args:
            doy (int): day of year [1-366]

        Returns:
            f (array): fractional LAI multiplier [-], ranging from lai_min to 1.0
        """

        lai_min = self.phenopara['LAI_decid_min']
        ddo = self.phenopara['ddo']
        ddur = self.phenopara['ddur']
        sso = self.phenopara['sso']
        sdur = self.phenopara['sdur']

        # growth phase
        self._growth_stage += 1.0 / ddur
        f = np.minimum(1.0, lai_min + (1.0 - lai_min) * self._growth_stage)

        # beginning of year
        ix = np.where(self.DDsum <= ddo)
        f[ix] = lai_min
        self._growth_stage[ix] = 0.
        self._senesc_stage[ix] = 0.

        # senescence phase
        ix = np.where(doy > sso)
        self._growth_stage[ix] = 0.
        self._senesc_stage[ix] += 1.0 / sdur
        f[ix] = 1.0 - (1.0 - lai_min) * np.minimum(1.0, self._senesc_stage[ix])

        # update self.LAIdecid and total LAI
        self._LAIdecid = self._LAIdecid_max * f
        self._LAIgrass = self._LAIgrass_max * f
        self.LAI = self._LAIconif + self._LAIdecid + self._LAIshrub + self._LAIgrass
        return f

    def dry_canopy_et(self, D, Qp, AE, Ta, Ra=25.0, Ras=250.0, CO2=380.0, Rew=1.0, beta=1.0, fPheno=1.0):
        """
        Computes transpiration and forest floor evaporation under dry-canopy
        conditions using a two-layer Penman-Monteith approach.

        Canopy conductance (Gc) is upscaled from leaf-level stomatal conductance
        using LAI, light attenuation, soil moisture, CO2, and phenology responses.

        Args:
            D      (array): vapor pressure deficit [kPa]
            Qp     (array): photosynthetically active radiation [W m-2]
            AE     (array): available energy / net radiation [W m-2]
            Ta     (array): air temperature [degC]
            Ra     (float or array): canopy aerodynamic resistance [s m-1]
            Ras    (float or array): soil aerodynamic resistance [s m-1]
            CO2    (float): atmospheric CO2 concentration [ppm]
            Rew    (array): relative extractable water [-]
            beta   (array): relative soil conductance for floor evaporation [-]
            fPheno (array): phenology modifier [-]

        Returns:
            Tr     (array): transpiration rate [mm s-1]
            Efloor (array): forest floor evaporation rate [mm s-1]
            Gc     (array): canopy conductance [m s-1]
        """

        # ---Amax and g1 as LAI -weighted average of conifers, decid, shrub and grass.
        # NOTE: currently all plant types share the same amax value.
        rhoa = 101300.0 / (8.31 * (Ta + 273.15)) # mol m-3
        Amax = 1./self.LAI * (self._LAIconif * self.physpara['amax']
                + self._LAIdecid *self.physpara['amax']
                             + self._LAIgrass * self.physpara['amax']
                             + self._LAIshrub * self.physpara['amax']) # umolm-2s-1

        g1 = 1./self.LAI * (self._LAIconif * self.physpara['g1_conif']
                + self._LAIdecid *self.physpara['g1_decid']
                           + self._LAIgrass * self.physpara['g1_grass']
                           + self._LAIshrub * self.physpara['g1_shrub'])

        kp = self.physpara['kp']  # (-) attenuation coefficient for PAR
        q50 = self.physpara['q50']  # Wm-2, half-sat. of leaf light response
        rw = self.physpara['rw']  # rew parameter
        rwmin = self.physpara['rwmin']  # rew parameter

        tau = np.exp(-kp * self.LAI)  # fraction of Qp at ground relative to canopy top

        """--- canopy conductance Gc (integrated stomatal conductance)----- """

        # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        fQ = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * self.LAI) + q50 + eps))

        # the next formulation is from Leuning et al., 2008 WRR for daily Gc; they refer to
        # Kelliher et al. 1995 AFM but the resulting equation is not exact integral of K95.
        # fQ = 1./ kp * np.log((Qp + q50) / (Qp*np.exp(-kp*self.LAI) + q50))

        # soil moisture response: Lagergren & Lindroth, xxxx"""
        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))

        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)

        # leaf level light-saturated gs (m/s)
        gs = np.minimum(1.6*(1.0 + g1 / np.sqrt(D))*Amax / 380. / rhoa, 0.1)  # large values if D -> 0

        # canopy conductance
        Gc = gs * fQ * fRew * fCO2 * fPheno
        Gc[np.isnan(Gc)] = eps
        
        """ --- transpiration rate --- """
        Tr = penman_monteith((1.-tau)*AE, 1e3*D, Ta, Gc, 1./Ra, units='mm')
        Tr[Tr < 0] = 0.0

        """--- forest floor evaporation rate--- """
        # soil conductance is function of relative water availability
        # gcs = 1. / self.soilrp * beta**2.0
        # beta = Wliq / FC; Best et al., 2011 Geosci. Model. Dev. JULES
        Gcs = self.gsoil

        Efloor = beta * penman_monteith(tau * AE, 1e3*D, Ta, Gcs, 1./Ras, units='mm')
        Efloor[self.SWE > 0] = 0.0  # no evaporation from floor if snow on ground or beta == 0
        
        return Tr, Efloor, Gc


    def canopy_water_snow(self, dt, T, Prec, AE, D, Ra=25.0, U=2.0):
        """
        Computes canopy water interception, evaporation/sublimation, and
        snowpack dynamics for one timestep.

        Precipitation phase is partitioned into rain and snow by temperature.
        Canopy interception follows the asymptotic approach of Hedstrom &
        Pomeroy (1998). Snowpack tracks separate ice (SWEi) and liquid (SWEl)
        fractions; melt and freeze follow degree-day formulations. Potential
        infiltration is liquid water in excess of the snowpack retention capacity.

        Args:
            dt   (float):          timestep duration [s]
            T    (array or float): air temperature [degC]
            Prec (array or float): precipitation amount [mm per timestep]
            AE   (array or float): available energy / net radiation [W m-2]
            D    (array or float): vapor pressure deficit [kPa]
            Ra   (array or float): canopy aerodynamic resistance [s m-1]
            U    (array or float): wind speed at canopy height [m s-1],
                                used for snow sublimation resistance

        Returns:
            PotInf (array): potential infiltration to soil [mm]
            Trfall (array): throughfall to field layer or snowpack [mm]
            Evap   (array): evaporation/sublimation from canopy storage [mm]
            Interc (array): canopy interception [mm]
            MBE    (array): mass balance error of canopy + snowpack [mm]
            erate  (array): potential evaporation/sublimation rate from canopy [mm]
            Unload (array): snow unloaded from canopy due to warming [mm]
            fS     (array): snowfall fraction of precipitation [-]
            Sfall  (array): snowfall amount [mm]
            Rfall  (array): rainfall amount [mm]
        """

        # quality of precipitation
        Tmin = 0.0  # 'C, below all is snow
        Tmax = 1.0  # 'C, above all is water
        Tmelt = 0.0  # 'C, T when melting starts

        # storage capacities mm
        Wmax = self.wmax * self.LAI
        Wmaxsnow = self.wmaxsnow * self.LAI

        # melting/freezing coefficients mm/s
        Kmelt = self.Kmelt - 1.64 * self.cf / dt  # Kuusisto E, 'Lumi Suomessa'
        Kfreeze = self.Kfreeze

        kp = self.physpara['kp']
        tau = np.exp(-kp*self.LAI)  # fraction of Rn at ground

        # inputs to arrays, needed for indexing later in the code
        gridshape = np.shape(self.LAI)  # rows, cols

        if np.shape(T) != gridshape:
            T = np.ones(gridshape) * T
            Prec = np.ones(gridshape) * Prec

        # latent heat of vaporization (Lv) and sublimation (Ls) J kg-1
        Lv = 1e3 * (3147.5 - 2.37 * (T + 273.15))
        Ls = Lv + 3.3e5

        # compute 'potential' evaporation / sublimation rates for each grid cell
        Ga = 1. / Ra  # aerodynamic conductance

        # resistance for snow sublimation adopted from:
        # Pomeroy et al. 1998 Hydrol proc; Essery et al. 2003 J. Climate;
        # Best et al. 2011 Geosci. Mod. Dev.
        # ri = (2/3*rhoi*r**2/Dw) / (Ce*Sh*W) == 7.68 / (Ce*Sh*W

        Ce = 0.01*((self.W + eps) / Wmaxsnow)**(-0.4)  # exposure coeff [-]
        Sh = (1.79 + 3.0*U**0.5)  # Sherwood number [-]

        gi = np.where(T <= Tmin, Sh*self.W*Ce / 7.68 + eps, 1e6) # m s-1
        Lambda = np.where(T <= Tmin, Ls, Lv)
        # evaporation of interception storage, mm
        erate = np.where(Prec==0,
                         dt / Lambda * penman_monteith((1.0 - tau)*AE, 1e3*D, T, gi, Ga, units='W'),
                         0.0)

        # ---state of precipitation [as water (fW) or as snow(fS)]
        fW = np.where(T >= Tmax, 1.0, 0.0)

        ix = ((T > Tmin) & (T < Tmax))
        fW[ix] = (T[ix] - Tmin) / (Tmax - Tmin)

        fS = 1.0 - fW

        sf = fS * Prec
        rf = fW * Prec

        """ --- initial conditions for calculating mass balance error --"""
        Wo = self.W  # canopy storage
        SWEo = self.SWE  # Snow water equivalent mm

        """ --------- Canopy water storage change -----"""
        # snow unloading from canopy, ensures also that seasonal LAI development does
        # not mess up computations
        Unload = np.where(T >= Tmax, np.maximum(self.W - Wmax, 0.0), 0.0)
        self.W = self.W - Unload

        # Interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom & Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.
        # above Tmin, interception capacity equals that of liquid precip
        Interc = np.where(T < Tmin,
                (Wmaxsnow - self.W)* (1.0 - np.exp(-(self.cf / Wmaxsnow) * Prec)),
                np.maximum(0.0, (Wmax - self.W))* (1.0 - np.exp(-(self.cf / Wmax) * Prec)))

        self.W = self.W + Interc  # new canopy storage, mm

        Trfall = Prec + Unload - Interc  # Throughfall to field layer or snowpack

        # evaporate from canopy and update storage
        Evap = np.minimum(erate, self.W)  # mm
        self.W = self.W - Evap

        """ Snowpack (in case no snow, all Trfall routed to floor) """
        # melting positive, freezing negative
        Melt_Freeze = np.where(T >= Tmelt,
                np.minimum(self.SWEi, Kmelt * dt * (T - Tmelt)),
                -np.minimum(self.SWEl, Kfreeze * dt * (Tmelt - T)))

        # amount of water as ice and liquid in snowpack
        Sice = np.maximum(0.0, self.SWEi + fS * Trfall - Melt_Freeze)
        Sliq = np.maximum(0.0, self.SWEl + fW * Trfall + Melt_Freeze)

        PotInf = np.maximum(0.0, Sliq - Sice * self.R)  # mm
        Sliq = np.maximum(0.0, Sliq - PotInf)  # mm, liquid water in snow

        # update Snowpack state variables
        self.SWEl = Sliq
        self.SWEi = Sice
        self.SWE = self.SWEl + self.SWEi

        # mass-balance error mm
        MBE = (self.W + self.SWE) - (Wo + SWEo) - (Prec - Evap - PotInf)

        return PotInf, Trfall, Evap, Interc, MBE, erate, Unload, fS + fW, sf, rf


""" *********** utility functions ******** """

# @staticmethod
def e_sat(T, P=101300, Lambda=2450e3):
    """
    Computes saturation vapor pressure, slope of the saturation vapor
    pressure curve, and the psychrometric constant.

    Args:
        T      (array or float): air temperature [degC]
        P      (float):          ambient pressure [Pa]
        Lambda (float):          latent heat of vaporization [J kg-1]

    Returns:
        esa (array or float): saturation vapor pressure [Pa]
        s   (array or float): slope of saturation vapor pressure curve [Pa K-1]
        g   (array or float): psychrometric constant [Pa K-1]
    """
    cp = 1004.67  # J/kg/K

    esa = 1e3 * 0.6112 * np.exp((17.67 * T) / (T + 273.16 - 29.66))  # Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * Lambda)
    return esa, s, g

# @staticmethod
def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W'):
    """
    Computes latent heat flux LE (Wm-2) i.e evapotranspiration rate ET (mm/s)
    from Penman-Monteith equation
    Args:
        AE    (array or float): available energy [W m-2]
        D     (array or float): vapor pressure deficit [Pa]
        T     (array or float): air temperature [degC]
        Gs    (array or float): surface (canopy or soil) conductance [m s-1]
        Ga    (array or float): aerodynamic conductance [m s-1]
        P     (float):          ambient pressure [Pa]
        units (str):            output unit: 'W' (W m-2), 'mm' (mm s-1),
                                or 'mol' (mol m-2 s-1)

    Returns:
        x (array or float): evapotranspiration rate in requested units,
                            clipped to >= 0
    """
    # --- constants
    cp = 1004.67  # J kg-1 K-1
    rho = 1.25  # kg m-3
    Mw = 18e-3  # kg mol-1
    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))
    _, s, g = e_sat(T, P, L)  # slope of sat. vapor pressure, psycrom const

    x = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / (Gs + eps)))  # Wm-2

    if units == 'mm':
        x = x / L  # kgm-2s-1 = mms-1
    if units == 'mol':
        x = x / L / Mw  # mol m-2 s-1

    x = np.maximum(x, 0.0)
    return x

def aerodynamics(LAI, hc, Uo, w=0.01, zm=2.0, zg=0.5, zos=0.01):
    """
    computes wind speed at ground and canopy + boundary layer conductances
    Computes wind speed at ground height assuming logarithmic profile above and
    exponential within canopy
    Args:
        LAI (array):  one-sided leaf/plant area index [m2 m-2]
        hc  (array):  canopy height [m]
        Uo  (array or float): mean wind speed at height zm above canopy [m s-1]
        w   (float):  leaf length scale [m]
        zm  (float):  wind speed measurement height above canopy top [m]
        zg  (float):  reference height above ground [m]
        zos (float):  forest floor roughness length [m]

    Returns:
        ra    (array): canopy aerodynamic resistance [s m-1]
        rb    (array): canopy boundary layer resistance [s m-1]
        ras   (array): forest floor aerodynamic resistance [s m-1]
        ustar (array): friction velocity [m s-1]
        Uh    (array): wind speed at canopy top [m s-1]
        Ug    (array): wind speed at reference height zg [m s-1]
    References:
        Cammalleri et al. (2010). Hydrol. Earth Syst. Sci.
        Massman (1987). Boundary-Layer Meteorol., 40, 179-197.
        Magnani et al. (1998). Plant Cell Environ.
        Yi (2008). [wind attenuation coefficient]
    """
    zm = hc + zm  # m
    kv = 0.4  # von Karman constant (-)
    beta = 285.0  # s/m, from Campbell & Norman eq. (7.33) x 42.0 molm-3
    alpha = LAI / 2.0  # wind attenuation coeff (Yi, 2008 eq. 23)
    d = 0.66*hc  # m
    zom = 0.123*hc  # m
    zov = 0.1*zom
    zosv = 0.1*zos

    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * kv / np.log((zm - d) / zom)
    Uh = ustar / kv * np.log((hc - d) / zom)

    # U(zg) from exponential wind profile
    zn = np.minimum(zg / hc, 1.0)  # zground can't be above canopy top
    Ug = Uh * np.exp(alpha*(zn - 1.0))

    # canopy aerodynamic & boundary-layer resistances (sm-1). Magnani et al. 1998 PCE eq. B1 & B5
    ra = 1./(kv**2.0 * Uo) * np.log((zm - d) / zom) * np.log((zm - d) / zov)
    rb = 1. / LAI * beta * ((w / Uh)*(alpha / (1.0 - np.exp(-alpha / 2.0))))**0.5

    # soil aerodynamic resistance (sm-1)
    ras = 1. / (kv**2.0*Ug) * (np.log(zg / zos))*np.log(zg / (zosv))

    ra = ra + rb
    return ra, rb, ras, ustar, Uh, Ug

def wind_profile(LAI, hc, Uo, z, zm=2.0, zg=0.2):
    """
    Computes wind speed at ground height assuming logarithmic profile above and
    hyperbolic cosine profile within canopy
    Args:
        LAI (float):          one-sided leaf/plant area index [m2 m-2]
        hc  (float):          canopy height [m]
        Uo  (float):          mean wind speed at height zm above canopy [m s-1]
        z   (array):          heights at which to compute wind speed [m]
        zm  (float):          wind speed measurement height above canopy top [m]
        zg  (float):          reference height above ground [m]

    Returns:
        U     (array): wind speed at each height in z [m s-1]
        ustar (float): friction velocity [m s-1]
        Uh    (float): wind speed at canopy top [m s-1]

    References:
        Cammalleri et al. (2010). Hydrol. Earth Syst. Sci.
        Massman (1987). Boundary-Layer Meteorol., 40, 179-197.
    """

    k = 0.4  # von Karman const
    Cd = 0.2  # drag coeff
    alpha = 1.5  # (-)

    zm = zm + hc
    d = 0.66*hc
    zom = 0.123*hc
    beta = 4.0 * Cd * LAI / (k**2.0*alpha**2.0)
    # solve ustar and U(hc) from log-profile above canopy
    ustar = Uo * k / np.log((zm - d) / zom)  # m/s

    U = np.ones(len(z))*np.nan

    # above canopy top wind profile is logarithmic
    U[z >= hc] = ustar / k * np.log((z[z >= hc] - d) / zom)

    # at canopy top, match log and exponential profiles
    Uh = ustar / k * np.log((hc - d) / zom)  # m/s

    # within canopy hyperbolic cosine profile
    U[z <= hc] = Uh * (np.cosh(beta * z[z <= hc] / hc) / np.cosh(beta))**0.5

    return U, ustar, Uh

def daylength(LAT, DOY):
    """
    Computes daylength from location and day of year.

    Args:
        LAT (float or array): latitude [degrees]
        DOY (int or array):   day of year [1-366]

    Returns:
        dl (float or array): daylength [hours]
    """
    CF = np.pi / 180.0  # conversion deg -->rad

    LAT = LAT*CF
    # ---> compute declination angle
    xx = 278.97 + 0.9856*DOY + 1.9165*np.sin((356.6 + 0.9856*DOY)*CF)
    DECL = np.arcsin(0.39785*np.sin(xx*CF))
    del xx

    # --- compute day length, the period when sun is above horizon
    # i.e. neglects civil twilight conditions
    cosZEN = 0.0
    value = cosZEN - np.sin(LAT)*np.sin(DECL) / (np.cos(LAT)*np.cos(DECL))
    value = np.clip(value, -1, 1)  # Clamp the value to the valid range (otherwise invalid value in np.arccos)  
    dl = 2.0*np.arccos(value) / CF / 15.0  # hours

    return dl

def update_distributed_radiation(rad_coeff, doy, Rg):
    """
    Scales global radiation by a spatially distributed radiation coefficient.

    For day 366 (leap years), the coefficient for day 365 is used as an
    approximation.

    Args:
        rad_coeff (netCDF4.Dataset): NetCDF dataset containing variable 'c_rad'
                                     with dimensions [doy, lat, lon]
        doy       (int):             day of year [1-366]
        Rg        (float):           spatially uniform global radiation [W m-2]

    Returns:
        Rg (array): spatially distributed global radiation [W m-2]
    """

    doy = int(doy)
    if doy <= 365:
        Rg = np.array(rad_coeff['c_rad'][doy-1,:,:]) * Rg
    elif doy == 366:
        Rg = np.array(rad_coeff['c_rad'][doy-2,:,:]) * Rg
    return Rg