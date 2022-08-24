# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:01:22 2022

@author: janousu
"""

import numpy as np
import pandas as pd
import sklearn.metrics

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