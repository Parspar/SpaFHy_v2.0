# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:01:22 2022

@author: janousu
"""

import numpy as np
import pandas as pd

# SPAFHY METRICS


def histogram_match(data1, data2, lims,  bins=25):
    hobs,binobs = np.histogram(data1,bins=25, range=lims)
    hsim,binsim = np.histogram(data2,bins=25, range=lims)
    hobs=np.float64(hobs)
    hsim=np.float64(hsim)
    minima = np.minimum(hsim, hobs)
    gamma = round(np.sum(minima)/np.sum(hobs),2)
    return str(gamma)