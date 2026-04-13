# -*- coding: utf-8 -*-
"""
TOPMODEL implementation for catchment-scale groundwater and streamflow
simulation in SpaFHy.

Topmodel_Homogenous assumes spatially uniform effective soil depth (m) and
saturated hydraulic conductivity (ko), with hydrologic similarity determined
from the topographic wetness index (TWI = log(a / tan(b))).

The catchment average saturation deficit S [m] is the single state variable.
Baseflow follows an exponential recession; return flow (saturation excess)
occurs where the local deficit s becomes negative.

References:
    Beven, K.J. & Kirkby, M.J. (1979). A physically based, variable contributing
        area model of basin hydrology. Hydrol. Sci. Bull., 24(1), 43-69.
    Launiainen et al. (2019). Hydrol. Earth Syst. Sci., 23, 3457-3480.
    Nousu et al. (2024). Hydrol. Earth Syst. Sci., 28, 4643-4666.
    Tyystjärvi et al. (2022). For. Ecol. Manage., 522, 120447.

@authors: slauniai, jpnousu
"""

import numpy as np
eps = np.finfo(float).eps  # machine epsilon

class Topmodel_Homogenous():
    def __init__(self, pp, S_initial=None):
        """
        Initializes TOPMODEL assuming homogeneous effective soil depth and
        saturated hydraulic conductivity. Hydrologic similarity is determined
        from TWI = log(a / tan(b)).

        Args:
            pp (dict): Parameter dictionary with keys:
                'dt'          [s]     timestep duration
                'dxy'         [m]     grid cell size
                'ko'          [m s-1] saturated hydraulic conductivity
                'm'           [m]     effective soil depth (transmissivity decay factor)
                'twi_cutoff'  [%]     percentile used to clip the upper tail of the
                                      TWI distribution (removes stream network outliers)
                'so'          [m]     initial catchment average saturation deficit
                'twi'         [-]     2D array of topographic wetness index values
                'flowacc'     [m]     2D array of flow accumulation per unit contour length
                'slope'       [deg]   2D array of local slope
            S_initial (float, optional): Initial catchment average saturation deficit [m].
                Overrides pp['so'] if provided.
        """
        if S_initial is None:
            S_initial = pp['so']

        # importing grids from parameters
        cmask = pp['twi'].copy()
        cmask[np.isfinite(cmask)] = 1.0
        flowacc = pp['flowacc']
        slope = pp['slope']
        twi = pp['twi']
        
        # importing other parameters
        dxy  = pp['dxy'] # grid size
        self.M = pp['m'] # effective soil depth [m]
        self.dt = float(pp['dt']) # timestep 
        self.To = pp['ko']*self.dt # transmissivity at saturation

        # area of a given cell
        self.CellArea = dxy**2
        self.CatchmentArea = np.size(cmask[np.isfinite(cmask)])*self.CellArea
        self.qr = np.full_like(cmask, 0.0)
        

        # local and catchment average hydrologic similarity indices (xi=twi, X).
        # Set xi > twi_cutoff equal to cutoff value to remove the tail of the
        # TWI distribution. This mainly affects stream network cells, where
        # outlier TWI values degrade streamflow predictions.

        # local indices
        self.xi = twi

        # apply cutoff
        clim = np.percentile(self.xi[self.xi > 0], pp['twi_cutoff'])
        self.xi[self.xi > clim] = clim

        # catchment average indice
        self.X = 1.0 / self.CatchmentArea*np.nansum(self.xi*self.CellArea)

        # baseflow rate when catchment Smean=0.0
        self.Qo = self.To*np.exp(-self.X)

        # catchment average saturation deficit S [m] is the only state variable
        s = self.local_s(S_initial)
        s[s < 0] = 0.0
        self.S = np.nanmean(s)


    def local_s(self, Smean):
        """
        Computes local storage deficit from the catchment average deficit.

        Args:
            Smean (float): Catchment average saturation deficit [m].

        Returns:
            s (array): Local saturation deficit [m] at each grid cell.
        """
        s = Smean + self.M*(self.X - self.xi)
        return s

    def subsurfaceflow(self):
        """
        Computes subsurface (base)flow to the stream network based on the
        current catchment average saturation deficit.

        Returns:
            Qb (float): Baseflow per unit catchment area [m per timestep].
        """
        Qb = self.Qo*np.exp(-self.S / (self.M + eps))
        return Qb

    def run_timestep(self, R):
        """
        Runs one timestep: updates catchment average saturation deficit and
        returns catchment-scale water balance fluxes.

        Args:
            R (float): Recharge to the saturated zone [m per unit catchment area]
                per timestep. Typically the spatial mean of BucketGrid drainage.

        Returns:
            dict with keys:
                'baseflow'                [mm]: subsurface flow to stream network
                'returnflow'              [mm]: catchment average saturation-excess return flow
                'local_returnflow'        [mm]: gridded return flow
                'drainage_in'             [mm]: recharge input (R)
                'water_closure'           [mm]: mass balance error (should be ~0)
                'saturation_deficit'       [m]: updated catchment average saturation deficit
                'local_saturation_deficit'[mm]: gridded local saturation deficit
                'saturated_area'           [-]: fraction of catchment that is saturated
                'storage_change'          [mm]: change in saturated zone storage
        """
        # initial conditions
        So = self.S
        s = self.local_s(So)

        # subsurface flow, based on initial state
        Qb = self.subsurfaceflow()

        # update storage deficit and check where we have returnflow
        S = So + Qb - R
        s = self.local_s(S)

        # returnflow grid
        self.qr = -s
        self.qr[self.qr < 0] = 0.0  # returnflow grid, m

        # average returnflow per unit area
        Qr = np.nansum(self.qr)*self.CellArea / self.CatchmentArea

        # now all saturation excess is in Qr so update s and S.
        # Deficit increases when Qr is removed
        S = S + Qr
        self.S = S
        s = s + self.qr
        # saturated area fraction
        ix = np.where(s <= 0)

        fsat = len(ix[0])*self.CellArea / self.CatchmentArea

        # check mass balance
        dS = (So - self.S)
        dF = R - Qb - Qr
        mbe = dS - dF


        results = {
                'baseflow': Qb * 1e3,  # [mm d-1]
                'returnflow': Qr * 1e3, #[mm d-1]
                'local_returnflow': self.qr * 1e3, # [mm]
                'drainage_in': R * 1e3, #[mm d-1]
                'water_closure': mbe * 1e3,  # [mm]
                'saturation_deficit': self.S, # [m]
                'local_saturation_deficit': s * 1e3, # [mm]
                'saturated_area': fsat, #[-],
                'storage_change': dF *1e3 # [mm d-1]
                }

        return results

    
def twi(flowacc, dxy, slope_rad, twi_method):
    """
    Computes the topographic wetness index (TWI) grid.

    Supports two methods:
        'twi': Standard TWI (Launiainen et al. 2019):
               TWI = log(a / (dxy * tan(slope)))
        'swi': SAGA wetness index (Tyystjärvi et al. 2022), which modifies the
               specific catchment area using a slope-dependent correction to
               better represent wetness in flat terrain.

    Args:
        flowacc   (array): Flow accumulation per unit contour length [m].
        dxy       (float): Grid cell size [m].
        slope_rad (array): Local slope [radians].
        twi_method  (str): Method to use; either 'twi' or 'swi'.

    Returns:
        xi (array): Topographic wetness index [-].
    """
    
    from scipy.ndimage import maximum_filter
    eps = np.finfo(float).eps  # machine epsilon
    
    if twi_method == 'twi':
        xi = np.log(flowacc / dxy / (np.tan(slope_rad) + eps))
    elif twi_method == 'swi':
        footprint = np.array([[1, 1, 1], [1, 0, 1],[1, 1, 1]])
        scamax = maximum_filter(flowacc/dxy, footprint=footprint)
        scamax_mod = scamax * (1/15)**slope_rad*np.exp(slope_rad)**15
        scam = np.maximum(scamax_mod, flowacc/dxy)
        xi = np.log(scam / dxy / (np.tan(slope_rad) + eps))
    return xi