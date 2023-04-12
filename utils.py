from scipy.ndimage import maximum_filter
import numpy as np
eps = np.finfo(float).eps  # machine epsilon

def twi(flowacc, dx, slope_rad, twi_method):
    """
    computes TWI according to 'standard' method (Launiainen et al. 2019) or
    as 'saga' wetness index (Tyystjärvi et al. 2022)
    """        
    if twi_method == 'standard':
        xi = np.log(flowacc / dx / (np.tan(slope_rad) + eps))
    elif twi_method == 'saga':
        footprint = np.array([[1, 1, 1], [1, 0, 1],[1, 1, 1]])
        scamax = maximum_filter(flowacc/dx, footprint=footprint)
        scamax_mod = scamax * (1/15)**slope_rad*np.exp(slope_rad)**15
        scam = np.maximum(scamax_mod, flowacc/dx)
        xi = np.log(scam / dx / (np.tan(slope_rad) + eps))
    return xi