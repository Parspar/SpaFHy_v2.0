from scipy.ndimage import maximum_filter
import numpy as np
eps = np.finfo(float).eps  # machine epsilon
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt

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

def wrc(pF: Dict, theta: np.ndarray=None, psi: np.ndarray=None, draw_pF: bool=False) -> np.ndarray:
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

    Ts = np.array(pF['ThetaS'])
    Tr = np.array(pF['ThetaR'])
    alfa = np.array(pF['alpha'])
    n = np.array(pF['n'])
    m = 1.0 - np.divide(1.0, n)

    def theta_psi(x):
        # converts water content [m3 m-3] to potential [m]]
        x = np.minimum(x, Ts)
        x = np.maximum(x, Tr)  # checks limits
        s = (Ts - Tr) / ((x - Tr) + eps)
        Psi = -1e-2 / alfa*(s**(1.0 / m) - 1.0)**(1.0 / n)  # m
        Psi[np.isnan(Psi)] = 0.0
        return Psi

    def psi_theta(x):
        # converts water potential [m] to water content [m3 m-3]
        x = 100*np.minimum(x, 0)  # cm
        Th = Tr + (Ts - Tr) / (1 + abs(alfa*x)**n)**m
        return Th

    # --- convert between theta <-- --> psi
    if theta:
        y = theta_psi(theta)  # 'Theta-->Psi'
    elif psi:
        y = psi_theta(psi)  # 'Psi-->Theta'

    # draws pf-curve
    if draw_pF:
        Ts = Ts[0]
        Tr = Tr[0]
        alfa = alfa[0]
        n = n[0]  
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

