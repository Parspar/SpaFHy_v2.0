# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:01:22 2022

@author: janousu
"""

import numpy as np
import pandas as pd
import sklearn.metrics
from scipy.stats import pearsonr

# SPAFHY METRICS


def histogram_match(data1, data2, lims,  bins=25):
    hobs,binobs = np.histogram(data1,bins=25, range=lims)
    hsim,binsim = np.histogram(data2,bins=25, range=lims)
    hobs=np.float64(hobs)
    hsim=np.float64(hsim)
    minima = np.minimum(hsim, hobs)
    gamma = round(np.sum(minima)/np.sum(hobs),2)
    return str(gamma)


def qq_plot_prep(data, columns=[]):
    data_prep = data[columns].reset_index(drop=True).dropna()
    data_prep[columns[0]] = np.sort(data_prep[columns[0]])
    data_prep[columns[1]] = np.sort(data_prep[columns[1]])
    return np.array(data_prep[columns[0]]), np.array(data_prep[columns[1]])


def R2_metrics(data, columns=[]):
    ind = ~data[columns].isnull().any(axis=1)
    R2 = sklearn.metrics.r2_score(data[ind][columns[0]], data[ind][columns[1]]).round(2)
    return R2

def MAE_metrics(data, columns=[]):
    ind = ~data[columns].isnull().any(axis=1)
    MAE = sklearn.metrics.mean_absolute_error(data[ind][columns[0]], data[ind][columns[1]]).round(2)
    return MAE

def MBE_metrics(data, columns=[]):
    ''' observations first column, simulations second column'''
    ind = ~data[columns].isnull().any(axis=1)
    MBE = round(np.mean(data[ind][columns[1]] - data[ind][columns [0]]),2)
    return MBE

def calc_kge(obs, sim):
    """Calculate the Kling-Gupta-Efficiency.
    
    Calculate the original KGE value following [1].

    Args:
        obs: Array of the observed values
        sim: Array of the simulated values

    Returns:
        The KGE value for the simulation, compared to the observation.

    Raises:
        ValueError: If the arrays are not of equal size or have non-numeric
            values.
        TypeError: If the arrays is not a supported datatype.
        RuntimeError: If the mean or the standard deviation of the observations
            equal 0.
    
    [1] Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). 
    Decomposition of the mean squared error and NSE performance criteria: 
    Implications for improving hydrological modelling. Journal of Hydrology, 
    377(1-2), 80-91.
    
    """
    # Validation check on the input arrays
    #obs = validate_array_input(obs, np.float64, 'obs')
    #sim = validate_array_input(sim, np.float64, 'sim')
    
    if len(obs) != len(sim):
        raise ValueError("Arrays must have the same size.")
     
    mean_obs = np.mean(obs)
    if mean_obs == 0:
        msg = "KGE not definied if the mean of the observations equals 0."
        raise RuntimeError(msg)
    
    std_obs = np.std(obs)
    if std_obs == 0:
        msg = ["KGE not definied if the standard deviation of the ",
               "observations equals 0."]
        raise RuntimeError("".join(msg))
    
    r = pearsonr(obs, sim)[0]
    alpha = np.std(sim) / std_obs
    beta = np.mean(sim) / mean_obs
    
    kge_val = 1 - np.sqrt((r-1)**2 + (alpha-1)**2 + (beta-1)**2)
    
    return kge_val 