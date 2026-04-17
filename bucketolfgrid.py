# -*- coding: utf-8 -*-
"""
Combined overland flow routing and two-layer soil water bucket model for
gridded application in SpaFHy.

Cells are processed sequentially in flow-accumulation order (low to high).
Within each cell two steps execute in order:

  1. Overland flow: excess pond storage above MaxPondSto is routed to the D8
     downslope neighbour. Water body cells drain all pond storage directly to
     the ditch network.

  2. Bucket water balance: vertical processes (interception, evaporation,
     transpiration, Campbell gravitational drainage) are solved using the
     pond storage that remains after step 1 as an additional water input.

Processing cells from lowest to highest flow accumulation ensures that water
routed from an upstream cell is immediately available for infiltration when
that cell's downstream neighbour is processed in the same timestep.

The class interface is intentionally close to BucketGrid so that it can serve
as a drop-in replacement in SpaFHy simulations that require explicit overland
flow routing. The sequential Python loop is kept simple; Numba JIT compilation
can be added later if performance becomes a bottleneck.

References:
    Campbell, G.S. (1974). A simple method for determining unsaturated
        conductivity from moisture retention data. Soil Science, 117(6).
    Launiainen et al. (2019). Hydrol. Earth Syst. Sci., 23, 3457-3480.
    Nousu et al. (2024). Hydrol. Earth Syst. Sci., 28, 4643-4666.

@authors: jpnousu, slauniai
"""

import numpy as np
eps = np.finfo(float).eps


class BucketOLFGrid(object):
    """
    Combined overland flow routing and two-layer soil water bucket model for
    gridded use in SpaFHy.

    Overland flow and vertical soil water balance are solved together in a
    single cell-by-cell loop ordered by flow accumulation.
    """

    def __init__(self, spara, org_drain):
        """
        Args:
            spara (dict): Soil parameter dictionary. All values are np.arrays
                of the grid shape. Expected keys:

                Organic top layer:
                    'org_depth'    [m]       thickness of organic layer
                    'org_poros'    [m3 m-3]  porosity
                    'org_fc'       [m3 m-3]  field capacity
                    'org_rw'       [m3 m-3]  parameter for relative evaporation rate
                    'org_ksat'     [m s-1]   saturated hydraulic conductivity (org_drain=True only)
                    'org_beta'     [-]        Campbell exponent (org_drain=True only)

                Root zone layer:
                    'root_depth'   [m]       layer thickness
                    'root_poros'   [m3 m-3]  porosity
                    'root_fc'      [m3 m-3]  field capacity
                    'root_wp'      [m3 m-3]  wilting point
                    'root_ksat'    [m s-1]   saturated hydraulic conductivity
                    'root_beta'    [-]        Campbell exponent
                    'root_alpha'   [kPa-1]   van Genuchten alpha
                    'root_n'       [-]        van Genuchten n
                    'root_wr'      [m3 m-3]  residual water content

                Ponding:
                    'maxpond'      [m]       maximum above-ground pond storage

                Initial conditions:
                    'pond_storage' [m]       initial pond storage
                    'org_sat'      [-]       initial saturation of organic layer
                    'root_sat'     [-]       initial saturation of root zone
                    'top_storage'  [m]       (optional) overrides org_sat-based init
                    'root_storage' [m]       (optional) overrides root_sat-based init

                Overland flow routing (required when used as BucketOLFGrid):
                    'flowacc'  [m2]: flow accumulation raster
                    'fdir'     [-]:  D8 flow direction (ArcGIS encoding: N=64,
                                     NE=128, E=1, SE=2, S=4, SW=8, W=16, NW=32)
                    'streams'  [-]:  stream network raster (0 = no stream)
                    'lakes'    [-]:  lake raster (0 = no lake)

            org_drain (bool): If True, the organic top layer drains gravitationally
                to the root zone using Campbell hydraulic conductivity. If False,
                the organic layer acts as an interception store up to field capacity.
        """

        # --- overland flow / D8 routing setup ---
        flowacc   = spara['flowacc']
        self.fdir = spara['fdir']
        self.rows, self.cols = flowacc.shape

        streams = np.where(spara['streams'] == 0, np.nan, spara['streams'])
        lakes   = np.where(spara['lakes']   == 0, np.nan, spara['lakes'])
        self.water_bodies = np.where(
            np.isfinite(streams) | np.isfinite(lakes), 1.0, np.nan)

        self.valid_mask = ~np.isnan(flowacc) & ~np.isnan(self.fdir)
        flat_valid = np.flatnonzero(self.valid_mask)
        self.sorted_indices = flat_valid[np.argsort(flowacc.ravel()[flat_valid])]

        self.direction_offsets = {
            64:  (-1,  0),   # North
            128: (-1,  1),   # Northeast
            1:   ( 0,  1),   # East
            2:   ( 1,  1),   # Southeast
            4:   ( 1,  0),   # South
            8:   ( 1, -1),   # Southwest
            16:  ( 0, -1),   # West
            32:  (-1, -1),   # Northwest
        }

        # --- bucket model setup ---
        self.org_drain = org_drain
        self.MaxPond   = spara['maxpond']

        # organic top layer
        self.D_top     = spara['org_depth']
        self.poros_top = spara['org_poros']
        self.Fc_top    = spara['org_fc']
        self.rw_top    = spara['org_rw']

        if self.org_drain:
            self.MaxStoTop = self.poros_top * self.D_top
            self.Ksat_top  = spara['org_ksat']
            self.beta_top  = spara['org_beta']
        else:
            self.MaxStoTop = self.Fc_top * self.D_top

        self.MaxStoTopInt = self.Fc_top * self.D_top   # interception capacity

        # root zone layer
        self.D_root     = spara['root_depth']
        self.poros_root = spara['root_poros']
        self.Fc_root    = spara['root_fc']
        self.Wp_root    = spara['root_wp']
        self.Ksat_root  = spara['root_ksat']
        self.beta_root  = spara['root_beta']
        self.alpha_root = spara['root_alpha']
        self.n_root     = spara['root_n']
        self.wr_root    = spara['root_wr']
        self.MaxStoRoot = self.D_root * self.poros_root

        # initial states
        self.PondSto    = np.minimum(spara['pond_storage'], self.MaxPond)
        self.WatStoTop  = spara.get('top_storage',
                                    self.MaxStoTop * spara['org_sat'])
        self.WatStoRoot = spara.get('root_storage',
                                    np.minimum(spara['root_sat'] * self.MaxStoRoot,
                                               self.MaxStoRoot))

        # drainage state arrays
        if self.org_drain:
            self.drain_top = np.where(np.isfinite(self.WatStoTop), 0.0, np.nan)
        self.drain    = np.where(np.isfinite(self.WatStoRoot), 0.0, np.nan)
        self.retflow  = np.full_like(self.WatStoRoot, 0.0)
        self._drainage_to_gw = 0.0

        # initialise all diagnostic state variables
        self.setState()

    # ------------------------------------------------------------------
    # Diagnostic state
    # ------------------------------------------------------------------

    def setState(self):
        """
        Updates all diagnostic state variables from the primary storage arrays
        (WatStoRoot, WatStoTop). Called at the end of run_timestep using
        vectorised numpy operations over the full grid.

        Updates:
            Wliq_root, Wair_root, Sat_root, Rew (root zone)
            Wliq_top, Wair_top, Sat_top, Ree    (organic top layer)
            Psi                                  (matric potential, MPa)
        """
        # root zone
        self.Wliq_root = self.poros_root * self.WatStoRoot / self.MaxStoRoot
        self.Wair_root = np.maximum(0.0, self.MaxStoRoot - self.WatStoRoot)
        self.Sat_root  = self.Wliq_root / self.poros_root
        self.Rew = np.maximum(0.0,
            np.minimum((self.Wliq_root - self.Wp_root)
                       / (self.Fc_root - self.Wp_root + eps), 1.0))

        # organic top layer
        self.Wliq_top = ((self.MaxStoTop / self.D_top)
                         * self.WatStoTop / (self.MaxStoTop + eps))
        self.Sat_top  = self.Wliq_top / self.poros_top
        self.Ree = self.relative_evaporation()
        self.Wliq_top[self.D_top == 0] = np.nan
        self.Wair_top = np.maximum(0.0, self.MaxStoTopInt - self.WatStoTop)
        self.Ree[self.D_top == 0] = eps
        self.Psi = self.theta_psi()

    def theta_psi(self):
        """
        Computes soil water potential from volumetric water content using the
        van Genuchten (1980) retention curve.

        Returns:
            Psi (array): matric potential [MPa], <= 0
        """
        n = self.n_root
        m = 1.0 - 1.0 / n
        x = np.minimum(self.Wliq_root, self.poros_root)
        x = np.maximum(x, self.wr_root)
        s = (self.poros_root - self.wr_root) / ((x - self.wr_root) + eps)
        Psi = -1.0 / self.alpha_root * (s**(1.0 / m) - 1.0)**(1.0 / n)
        return 1e-3 * Psi   # kPa to MPa

    def relative_evaporation(self):
        """
        Returns relative evaporation rate from the organic top layer, loosely
        based on Launiainen et al. (2015) Ecol. Mod. moss module.

        Returns:
            f (array): relative evaporation rate [-], 0 to 1
        """
        return np.maximum(0.0,
                          np.minimum(0.98 * self.Wliq_top / self.rw_top, 1.0))

    # ------------------------------------------------------------------
    # Time step
    # ------------------------------------------------------------------

    def run_timestep(self, dt=86400.0, rr=0.0, tr=0.0, evap=0.0,
                     airv_deep=1000.0, retflow=0.0):
        """
        Runs the combined overland flow routing and two-layer bucket model for
        one timestep.

        Cells are processed in flow-accumulation order (low to high). For each
        cell:
          1. Bucket step: vertical water balance is solved; any excess beyond
             the top layer accumulates in pond storage (no MaxPond cap here).
          2. Overland flow step: excess pond storage (> MaxPond) is routed to
             the D8 downslope neighbour. Water body cells record pond storage
             as surface runoff and clear the pond.

        Args:
            dt        (float): timestep duration [s]
            rr        (array): potential infiltration input to top layer [m]
            tr        (array): transpiration demand from root zone [m]
            evap      (array): evaporative demand from top layer [m]
            airv_deep (array or float): available air volume in deep soil /
                groundwater layer [m]; limits gravitational drainage from root
                zone
            retflow   (array or float): return flow from groundwater to root
                zone [m]

        Returns:
            dict with keys (all grid arrays):
                'potential_infiltration' [mm d-1]: rainfall input to top layer
                'evaporation'           [mm d-1]: evaporation from top layer
                'transpiration'         [mm d-1]: transpiration from root zone
                'drainage'              [mm d-1]: gravitational drainage from root zone
                'surface_runoff'        [mm d-1]: overland flow reaching stream/lake cells
                                                  plus any water exiting the catchment boundary
                'lateral_netflow'       [mm d-1]: net lateral overland flow per cell
                                                  (positive = net inflow from upslope)
                'return_flow'           [mm d-1]: return flow from groundwater
                'water_closure'         [mm d-1]: mass balance error (should be ~0)
                'moisture_top'         [m3 m-3]: volumetric water content, top layer
                'moisture_root'        [m3 m-3]: volumetric water content, root zone
                'psi_root'                [MPa]: matric potential of root zone
                'transpiration_limitation'  [-]: relative extractable water (REW)
                'water_storage_root'       [mm]: root zone water storage
                'water_storage_top'        [mm]: top layer water storage
                'pond_storage'             [mm]: above-ground pond storage
                'water_storage'            [mm]: total soil water storage (top + root)
                'storage_change'           [mm]: change in total storage over dt
        """
        gridshape    = np.shape(self.WatStoRoot)
        flux_to_mm_d = 1e3 * (86400.0 / dt)

        # broadcast scalar inputs to full grids
        def _to_grid(x):
            return x * np.ones(gridshape) if np.shape(x) != gridshape else x

        rr        = _to_grid(rr)
        tr        = _to_grid(tr)
        evap      = _to_grid(evap)
        airv_deep = _to_grid(airv_deep)
        retflow   = _to_grid(retflow)
        self.retflow = retflow

        rr0 = rr.copy()   # original rainfall; saved for output and MBE

        # save initial storage states
        PondSto0    = self.PondSto.copy()
        WatStoTop0  = self.WatStoTop.copy()
        WatStoRoot0 = self.WatStoRoot.copy()

        # initialise per-cell output accumulators (NaN outside catchment)
        nan_grid = np.where(self.valid_mask, 0.0, np.nan)
        evap_out    = nan_grid.copy()
        tr_out      = nan_grid.copy()
        drain_out   = nan_grid.copy()
        roff_out    = nan_grid.copy()
        lateral_in  = nan_grid.copy()
        lateral_out = nan_grid.copy()

        # ----------------------------------------------------------------
        # Main cell loop: process cells from lowest to highest flowacc.
        # Order within each cell: Bucket first, OLF second.
        # This ensures water routed from an upstream cell (OLF step) arrives
        # at the downstream cell before its Bucket step runs, so routed water
        # can infiltrate in the same timestep.
        # ----------------------------------------------------------------
        for flat_idx in self.sorted_indices:
            r, c = divmod(flat_idx, self.cols)

            # ---- step 1: bucket vertical water balance ----

            # total water input: rainfall + pond storage received from upstream
            rr_cell = rr[r, c] + self.PondSto[r, c]
            self.PondSto[r, c] = 0.0

            # local parameter scalars (avoids repeated 2-D indexing inside expressions)
            D_top        = self.D_top[r, c]
            poros_top    = self.poros_top[r, c]
            Fc_top       = self.Fc_top[r, c]
            MaxStoTop    = self.MaxStoTop[r, c]
            MaxStoTopInt = self.MaxStoTopInt[r, c]
            MaxPond      = self.MaxPond[r, c]
            D_root       = self.D_root[r, c]
            poros_root   = self.poros_root[r, c]
            Fc_root      = self.Fc_root[r, c]
            MaxStoRoot   = self.MaxStoRoot[r, c]
            retflow_cell = retflow[r, c]
            airv_cell    = airv_deep[r, c]

            # top layer: interception
            interc = (max(0.0, MaxStoTopInt - self.WatStoTop[r, c])
                      * (1.0 - np.exp(-rr_cell / (MaxStoTopInt + eps))))
            self.WatStoTop[r, c] = max(0.0, self.WatStoTop[r, c] + interc)

            # top layer: evaporation
            evap_cell = min(evap[r, c], self.WatStoTop[r, c])
            self.WatStoTop[r, c] -= evap_cell

            # top layer: optional Campbell drainage to root zone
            if self.org_drain:
                Wliq_top  = ((MaxStoTop / D_top)
                             * self.WatStoTop[r, c] / (MaxStoTop + eps))
                Sat_top   = Wliq_top / poros_top
                k_top     = (self.Ksat_top[r, c]
                             * Sat_top**(2.0 * self.beta_top[r, c] + 3.0))
                drain_top = min(k_top * dt,
                                max(0.0, (Wliq_top - Fc_top)) * D_top)
                self.WatStoTop[r, c] -= drain_top
                rr_to_root = rr_cell - interc + drain_top
            else:
                rr_to_root = rr_cell - interc

            # root zone: transpiration
            tr_cell = min(tr[r, c], self.WatStoRoot[r, c] - eps)
            self.WatStoRoot[r, c] -= tr_cell

            # root zone: Campbell gravitational drainage
            Wliq_root  = poros_root * self.WatStoRoot[r, c] / MaxStoRoot
            Sat_root   = Wliq_root / poros_root
            k_root     = (self.Ksat_root[r, c]
                          * Sat_root**(2.0 * self.beta_root[r, c] + 3.0))
            drain_cell = min(k_root * dt,
                             max(0.0, (Wliq_root - Fc_root)) * D_root)

            # suppress drainage where return flow is active (avoids oscillation)
            if retflow_cell > 0.0:
                drain_cell = 0.0
            drain_cell = min(drain_cell, airv_cell)

            # root zone: inflow and storage update
            Qin    = retflow_cell + rr_to_root
            inflow = min(Qin, MaxStoRoot - self.WatStoRoot[r, c] + drain_cell)
            self.WatStoRoot[r, c] = min(MaxStoRoot,
                max(self.WatStoRoot[r, c] + inflow - drain_cell, eps))

            # excess water cascade: try top layer first, then pond.
            # No MaxPond cap here — excess above MaxPond stays in PondSto
            # so the OLF step below can route it downslope.
            exfil  = max(0.0, Qin - inflow)
            to_top = max(0.0, min(exfil, MaxStoTop - self.WatStoTop[r, c] - eps))
            self.WatStoTop[r, c] += to_top
            self.PondSto[r, c]   += max(0.0, exfil - to_top)

            # accumulate cell outputs
            evap_out[r, c]  = evap_cell
            tr_out[r, c]    = tr_cell
            drain_out[r, c] = drain_cell

            # ---- step 2: overland flow routing ----
            if np.isfinite(self.water_bodies[r, c]):
                # water body cell: all pond storage exits as surface runoff
                roff_out[r, c] += self.PondSto[r, c]
                self.PondSto[r, c] = 0.0

            elif self.PondSto[r, c] > MaxPond:
                excess   = self.PondSto[r, c] - MaxPond
                flow_dir = self.fdir[r, c]
                if flow_dir in self.direction_offsets:
                    dr, dc = self.direction_offsets[flow_dir]
                    nr, nc = r + dr, c + dc
                    if (0 <= nr < self.rows and 0 <= nc < self.cols
                            and self.valid_mask[nr, nc]):
                        self.PondSto[nr, nc] += excess
                        lateral_in[nr, nc]   += excess
                        lateral_out[r, c]    += excess
                        self.PondSto[r, c]    = MaxPond
                    else:
                        # no valid downslope neighbour: excess exits as surface runoff
                        roff_out[r, c]     += excess
                        self.PondSto[r, c]  = MaxPond

        # ----------------------------------------------------------------
        # Post-loop: update state and compute grid-wide diagnostics
        # ----------------------------------------------------------------

        self.drain = drain_out
        self._drainage_to_gw = np.nansum(self.drain)
        self.setState()

        lateral_netflow = lateral_in - lateral_out

        dStorage = ((self.WatStoRoot - WatStoRoot0)
                    + (self.WatStoTop  - WatStoTop0)
                    + (self.PondSto    - PondSto0))

        # per-cell mass balance error [m]; lateral terms cancel when summed
        mbe = (dStorage
               - (rr0 + retflow - evap_out - tr_out - drain_out - roff_out
                  + lateral_in - lateral_out))

        results = {
            'potential_infiltration': rr0              * flux_to_mm_d,  # [mm d-1]
            'evaporation':            evap_out          * flux_to_mm_d,  # [mm d-1]
            'transpiration':          tr_out            * flux_to_mm_d,  # [mm d-1]
            'drainage':               drain_out         * flux_to_mm_d,  # [mm d-1]
            'surface_runoff':         roff_out          * flux_to_mm_d,  # [mm d-1]
            'lateral_netflow':        lateral_netflow   * flux_to_mm_d,  # [mm d-1]
            'return_flow':            retflow           * flux_to_mm_d,  # [mm d-1]
            'water_closure':          mbe               * flux_to_mm_d,  # [mm d-1]
            'moisture_top':           self.Wliq_top,                     # [m3 m-3]
            'moisture_root':          self.Wliq_root,                    # [m3 m-3]
            'psi_root':               self.Psi,                          # [MPa]
            'transpiration_limitation': self.Rew,                        # [-]
            'water_storage_root':     self.WatStoRoot   * 1e3,           # [mm]
            'water_storage_top':      self.WatStoTop    * 1e3,           # [mm]
            'pond_storage':           self.PondSto      * 1e3,           # [mm]
            'water_storage':          (self.WatStoTop
                                       + self.WatStoRoot) * 1e3,         # [mm]
            'storage_change':         dStorage          * 1e3,           # [mm]
        }

        return results
