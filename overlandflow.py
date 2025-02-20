# -*- coding: utf-8 -*-
"""
Created on Thu Jan 9 13:38:59 2025

@author: janousu
"""

import numpy as np

class OverlandFlowModel:
    def __init__(self, flowacc, fdir, streams, lakes, MaxPondSto=0.05):
        """
        Initialize the OverlandFlowModel with required constants.

        Parameters:
        - MaxPondSto (float): Maximum pond storage before water flows to the next cell.
        """

        # Max pond storage
        self.MaxPondSto = MaxPondSto

        # initialize flow rasters
        self.flowacc = flowacc
        self.fdir = fdir
        self.rows, self.cols = flowacc.shape

        # Convert 0 values in streams and lakes to np.nan
        streams = np.where(streams == 0, np.nan, streams)
        lakes = np.where(lakes == 0, np.nan, lakes)

        # Combine streams and lakes into a single water bodies array
        self.water_bodies = np.where(np.isfinite(
            streams) | np.isfinite(lakes), 1, np.nan)

        # Create a mask for valid cells (non-nan)
        self.valid_mask = ~np.isnan(flowacc) & ~np.isnan(fdir)

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

    def run_timestep(self, pondsto, airspace):
        """
        Runs a single timestep of overland flow and updates the internal state.

        Parameters:
        - pondsto (np.array): Current ponded water storage.
        - airspace (np.array): Available airspace for infiltration.

        Returns:
        - updated_pondsto (np.array): Updated pond storage.
        - updated_airspace (np.array): Updated airspace.
        - netflow_to_ditch (np.array): Netflow to ditch generated in this timestep.
        """

        # saving initial states for mbe calculation
        pondsto0 = pondsto.copy()
        airspace0 = airspace.copy()

        # Initialize arrays
        netflow_to_ditch = np.zeros_like(self.flowacc)
        mbe = np.zeros_like(self.flowacc)
        lateral_flow = np.zeros_like(self.flowacc)*self.flowacc

        # Get valid flat indices sorted by flow accumulation (low to high)
        flat_valid_indices = np.flatnonzero(self.valid_mask)
        sorted_indices = flat_valid_indices[np.argsort(
            self.flowacc.ravel()[flat_valid_indices])]

        # Process cells in sorted order
        for flat_index in sorted_indices:
            # Convert flat index to 2D indices
            r, c = divmod(flat_index, self.cols)

            if pondsto[r, c] > 0:  # Process cells with water

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
                                # Track lateral flow
                                #lateral_flow[r, c] = -excess_water  # Water leaving this cell
                                #lateral_flow[nr, nc] = excess_water  # Water entering neighbor cell

        lateral_flow = (pondsto0 - pondsto) - (airspace0 - airspace) - netflow_to_ditch

        # After processing all cells, remove any remaining excess water as surface runoff
        excess_water = np.maximum(pondsto - self.MaxPondSto, 0.)
        netflow_to_ditch += excess_water  # Add excess water to surface runoff
        pondsto -= excess_water  # Remove excess water from pond storage

        # Mass balance error calculation using airspace and lateral flow
        mbe = (pondsto0 - pondsto) - netflow_to_ditch - (airspace0 - airspace) - lateral_flow

        results = {
                'pond_storage': pondsto,  # [mm]
                'airspace': airspace, # [mm]
                'netflow_to_ditch': netflow_to_ditch,  # [mm d-1]
                'lateral_netflow': lateral_flow, # [mm d-1]
                'mbe': mbe, # [mm]
                }

        return results