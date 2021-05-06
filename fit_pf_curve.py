# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:15:32 2021

@author: janousu
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
eps = np.finfo(float).eps  # machine epsilon


def fit_pF(head, watcont, fig=False, labels=None, percentage=False, kPa=False, neg_heads=None, bounds=None):
    """
    """

    if fig:
        plt.figure()
    head = np.array(head)
    if kPa:
        head = head * 10  # kPa -> cm
    if neg_heads:
        head = -1.0 * head

    vg_ini = (0.88, 0.09, 0.03, 1.3)
    if not bounds:

        bounds = ((eps, eps, eps, 1.0), (1.0, 1.0, 5.0, 5.0))

    van_g = lambda h, *p:   p[1] + (p[0] - p[1]) / (1. + (p[2] * h) **p[3]) **(1. - 1. / p[3])
    vgen_all = []

    for k in range(0, len(watcont)):
        Wcont = np.array(watcont[k])
        print(Wcont)
        ix = np.where(Wcont >= 0)
        print(ix)
        if percentage:
            Wcont[ix] = Wcont[ix] / 100  # % -> fraction
        if labels is None:
            label = ''
        else:
            label = labels[k] + ': '
        try:
            #print(head, Wcont)
            vgen, _ = curve_fit(van_g, head[ix], Wcont[ix], p0=vg_ini, bounds=bounds)
            label+='Ts=%5.3f, Tr=%5.3f, alfa=%5.3f, n=%5.3f' % tuple(vgen)
        except RuntimeError:
            vgen = [-1, -1, -1, -1]
            label+='No fit!'
        vgen_all.append(vgen)

        if fig:

            plt.semilogy(Wcont[ix], head[ix], '.')
            xx = np.logspace(-1, 5.0, 100)
            plt.semilogy(van_g(xx, *vgen), xx, '-',
                         label=label)
            #c += 1

    if fig:
        plt.xlabel(r'$\theta$  $(m^3m^{-3})$', fontsize=14)
        plt.ylabel('$-h$ $(cm)$', fontsize=14)
        plt.ylim(xx[0], xx[-1])
        plt.xlim(0.0, 1.0)
        #plt.legend(loc='upper right', fontsize=7)
        plt.legend(loc='best', bbox_to_anchor=(1.0, 0.5, 0.5, 0.5))

    return vgen_all