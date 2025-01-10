# -*- coding: utf-8 -*-
"""
Created on Thu Jan 9 13:38:59 2025

@author: janousu
"""

import numpy as np
eps = np.finfo(float).eps
import matplotlib.pyplot as plt

class OverlandFlow(object):
    """
    Overland flow model for gridded use in SpaFHy.
    """
    def __init__(self, flowacc, flowdir):
        """
        Initializes OverlandFlow:
        Args:
            REQUIRED:
            spara - dictionary of ..... keys - values are np.arrays
    """
        
    self.flowacc = flowacc
    self.flowdir = flowdir

    # Define movement directions based on fdir
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
    
    self.rows, self.cols = self.flowacc.shape

    # Create a mask for valid cells (non-nan)
    self.valid_mask = ~np.isnan(self.flowacc) & ~np.isnan(self.fdir)

    # Flatten flowacc for sorting and maintain indices for valid cells only
    self.flat_indices = np.argsort(self.flowacc.ravel())  # Sort flowacc values (low to high)
    self.flat_indices = [i for i in self.flat_indices if self.valid_mask.ravel()[i]]  # Filter out invalid indices

    def run_timestep(excess_water, soil_airspace):
        # Iterate through cells in sorted order
        for flat_index in flat_indices:
            r, c = divmod(flat_index, cols)  # Convert flat index to 2D indices
            if pondsto[r, c] > 0:  # Only process cells with water
                # Determine infiltration capacity
                infiltration = min(pondsto[r, c], airspace[r, c])  # Amount of water that can infiltrate
                airspace[r, c] -= infiltration  # Reduce available air space
                pondsto[r, c] -= infiltration  # Reduce water in pondsto by infiltration

                # If there's still water left, flow it to the downstream cell
                if pondsto[r, c] > 0:
                    flow_direction = fdir[r, c]
                    if flow_direction in direction_offsets:
                        dr, dc = direction_offsets[flow_direction]
                        nr, nc = r + dr, c + dc  # Calculate the neighbor cell

                        # Check if the neighbor cell is within bounds and valid
                        if 0 <= nr < rows and 0 <= nc < cols and valid_mask[nr, nc]:
                            # Move water from current cell to the downstream cell
                            pondsto[nr, nc] += pondsto[r, c]
                            pondsto[r, c] = 0  # Set current cell water to 0 after flow

        return pondsto, airspace