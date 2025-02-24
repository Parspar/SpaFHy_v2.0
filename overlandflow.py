# -*- coding: utf-8 -*-
"""
Created on Thu Jan 9 13:38:59 2025

@author: janousu
"""

import numpy as np
from typing import Dict, List, Tuple

class OverlandFlowModel(object):
    def __init__(self, gisdata: Dict, params: Dict):
        """
        Initialize the OverlandFlowModel with required constants.

        Parameters:
        - gidsdata (dict):
            - flowacc (np.ndarray): Flow accumulation raster
            - fdir (np.ndarray): Flow direction raster
            - streams (np.ndarray): Stream raster
            - lakes (np.ndarray): Lake raster
        - params (dict)
            - MaxPondSto (float): Maximum pond storage
        """

        # Max pond storage
        self.MaxPondSto = params['MaxPondSto']

        # initialize 2D rasters that are used in run
        self.flowacc = gisdata['flowacc']
        self.fdir = gisdata['fdir']
        self.rows, self.cols = self.flowacc.shape
        # only used in initialization
        streams = gisdata['streams']
        lakes = gisdata['lakes']
        # Convert 0 values in streams and lakes to np.nan
        streams = np.where(streams == 0, np.nan, streams)
        lakes = np.where(lakes == 0, np.nan, lakes)
        # Combine streams and lakes into a single water bodies array
        self.water_bodies = np.where(np.isfinite(
            streams) | np.isfinite(lakes), 1, np.nan)
        # Create a mask for valid cells (non-nan)
        self.valid_mask = ~np.isnan(self.flowacc) & ~np.isnan(self.fdir)
        # Get valid flat indices sorted by flow accumulation (low to high)
        flat_valid_indices = np.flatnonzero(self.valid_mask)
        self.sorted_indices = flat_valid_indices[np.argsort(
            self.flowacc.ravel()[flat_valid_indices])]

        # Define movement directions based on flow direction (fdir)
        self.direction_offsets = {
            64: (-1, 0),   # North
            128: (-1, 1),  # Northeast
            1: (0, 1),     # East
            2: (1, 1),     # Southeast
            4: (1, 0),     # South
            8: (1, -1),    # Southwest
            16: (0, -1),   # West
            32: (-1, -1)   # Northwest
        }

    def run_timestep(self, forcing: Dict):
        """
        Runs a single timestep of overland flow and updates the internal state.

        Parameters:
        - forcing (Dict):
            - pond_storage (np.ndarray): Current ponded water storage.
            - air_space (np.ndarray): Available airspace for infiltration.

        Returns:
        - results (Dict)
            - pond_storage (np.ndarray):
            - air_space (np.ndarray):
            - netflow_to_ditch (np.ndarray):
            - lateral_netflow (np.ndarray):
            - surface_runoff (np.ndarray):
            - mbe
        """

        # saving initial states for mbe calculation, no airspace nor pondsto where there's water
        airspace = np.where(self.water_bodies == 1, 0., forcing['air_space'])
        pondsto = np.where(self.water_bodies == 1, 0., forcing['pond_storage'])
        pondsto0 = pondsto.copy()
        airspace0 = airspace.copy()

        # Initialize flux arrays
        netflow_to_ditch = np.zeros_like(self.flowacc)*self.flowacc
        lateral_flow = np.zeros_like(self.flowacc)*self.flowacc
        mbe = np.zeros_like(self.flowacc)*self.flowacc

        # calculating cells in sorted order
        for flat_index in self.sorted_indices:
            # Convert flat index to 2D indices
            r, c = divmod(flat_index, self.cols)

            if pondsto[r, c] > self.MaxPondSto:  # calculate only cells with water
                # Compute ditch flow for water body cells before infiltration
                # Check if the cell is a stream or lake
                if np.isfinite(self.water_bodies[r, c]):
                    netflow_to_ditch[r, c] = pondsto[r, c]
                    # Remove water that runs off
                    pondsto[r, c] = 0.
                else:
                    # Calculate excess water beyond max pond storage
                    excess_water = max(pondsto[r, c] - self.MaxPondSto, 0.)
                    # Compute infiltration
                    infiltration = min(excess_water, airspace[r, c])
                    airspace[r, c] -= infiltration
                    pondsto[r, c] -= infiltration
                    # Check if excess water remains for overland flow
                    excess_water = max(pondsto[r, c] - self.MaxPondSto, 0.)
                    if excess_water > 0:
                        flow_direction = self.fdir[r, c]
                        if flow_direction in self.direction_offsets:
                            dr, dc = self.direction_offsets[flow_direction]
                            nr, nc = r + dr, c + dc  # Neighbor cell
                            # Check if neighbor is within bounds and valid
                            if 0 <= nr < self.rows and 0 <= nc < self.cols and self.valid_mask[nr, nc]:
                                pondsto[nr, nc] += excess_water
                                # Keep max pond storage at current cell
                                pondsto[r, c] = self.MaxPondSto

        # NOTE lateral flow should take into account 'remaining_ponds'. Check also mbe with those.
        # Computing lateral flow based on initial states
        lateral_flow = (pondsto - pondsto0) - (airspace - airspace0) #+ remaining_ponds

        # After processing all cells, remove any remaining excess water as surface runoff
        remaining_ponds = np.maximum(pondsto - self.MaxPondSto, 0.)
        pondsto -= remaining_ponds  # Remove excess water from pond storage

        # saving surface runoff
        surface_runoff = remaining_ponds + netflow_to_ditch

        # Mass balance error calculation using airspace and lateral flow
        mbe = (pondsto - pondsto0) - (airspace - airspace0) - lateral_flow + remaining_ponds

        results = {
                'pond_storage': pondsto * 1e3,  # [mm]
                'air_space': airspace * 1e3, # [mm]
                'netflow_to_ditch': netflow_to_ditch  * 1e3,  # [mm d-1]
                'lateral_netflow': lateral_flow  * 1e3, # [mm d-1]
                'surface_runoff': surface_runoff * 1e3, # [mm d-1]
                'mbe': mbe * 1e3, # [mm]
                }

        return results