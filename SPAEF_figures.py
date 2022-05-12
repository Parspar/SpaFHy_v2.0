#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 11:33:33 2017
@ author:                  Mehmet Cüneyd Demirel 
@ author's webpage:        http://akademi.itu.edu.tr/demirelmc/
@ author's email id:       demirelmc@itu.edu.tr
@ author's website:        http://www.space.geus.dk/

A libray with Python functions for calculation of spatial efficiency (SPAEF) metric.

function:
    SPAEF : spatial efficiency   
"""

# import required modules
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
######################################################################################################################
def plot_SPAEFstats(sim,obs,SPAef, cc, alpha, histo,titlet, zs=True):
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 18,
        }
        
    newpath = r'.\temp_Maps' 
    if not os.path.exists(newpath):
        os.makedirs(newpath) 
        

    if zs==True:
        obs=(obs-np.nanmin(obs))/(np.nanmax(obs)-np.nanmin(obs))
        sim=(sim-np.nanmin(sim))/(np.nanmax(sim)-np.nanmin(sim))  
        
    bins =100# np.linspace(100, 350, 100)
    histrange = np.nanmin([sim, obs]), np.nanmax([sim,obs])
    hobs,binobs = np.histogram(obs,bins, histrange)
    hsim,binsim = np.histogram(sim,bins, histrange)

    fig, (ax1, ax2) = plt.subplots(1, 2 ,figsize=(16,10))
    gs = gridspec.GridSpec(1, 2)#, width_ratios=[1, 1])#, height_ratios=[1, 1]) 
    ax1 = plt.subplot(gs[0])
    ax1.hist(obs, bins, alpha=0.5, label='obs')
    ax1.hist(sim, bins, alpha=0.5, label='sim')
    ax1.legend(loc='upper left')
    ax1.text(min(obs), max(hsim)*2/3, r'$Histo Match$', fontdict=font)
    ax1.text(min(obs), max(hsim)*2/3-max(hsim)/13, '{:.2f}'.format(round(histo, 2)), fontdict=font)
    
    ax1.set_title("SPAEF is %.2f" %SPAef,fontdict=font)
    ax2 = plt.subplot(gs[1])
    ax2.scatter(obs,sim, alpha=0.01)
    plt.xlim(np.nanmin([sim, obs]) -0.05, np.nanmax([sim,obs])+0.05)
    plt.ylim(np.nanmin([sim, obs]) -0.05, np.nanmax([sim,obs])+0.05)
    
    #ax2.text(-0.011, -0.01, r'$CORR:  $', fontdict=font)
    #ax2.text(0.0195,-0.01, np.around(cc,2), fontdict=font)
    #ax2.text(-0.011, -0.05, r'$CV: $', fontdict=font)
    #ax2.text(0.0195,-0.05, np.around(alpha,2), fontdict=font)
       
    ax2.text(min(obs), (min(obs)*3/2)+0.05, r'$CORR: $', fontdict=font)
    ax2.text(min(obs), (min(obs)*3/2)-0.05, np.around(cc,2), fontdict=font)
    ax2.text(max(obs), max(obs)*2/3, r'$CV: $', fontdict=font)
    ax2.text(max(obs), max(obs)*2/3-max(obs)/13, np.around(alpha,2), fontdict=font)
    
    z = np.polyfit(obs, sim, 1)
    p = np.poly1d(z)
    ax2.plot(obs,p(obs),"r--")
    ax2.set_title("x = obs vs. y = sim", fontdict=font)

    plt.tight_layout()
    
    
    fig.savefig(newpath+'/'+titlet+'.png', dpi=60,bbox_inches='tight')

    plt.show()
    plt.close(fig)

print('figure saved')




def plot_maps(sim,obs,titlet):

    obs2=obs#(obs-np.nanmin(obs))/(np.nanmax(obs)-np.nanmin(obs))
    sim2=sim#(sim-np.nanmin(sim))/(np.nanmax(sim)-np.nanmin(sim))
    
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 18,
        }    
        
    font2 = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }  
    newpath = r'.\temp_Maps' 
    if not os.path.exists(newpath):
        os.makedirs(newpath) 

    fig = plt.figure(figsize=(16, 10)) 
 
    
    plt.subplot(1, 2, 1)
    x=plt.imshow(sim2, interpolation='none',cmap='Spectral_r')
    plt.axis('off') 
    plt.colorbar(x,orientation ='horizontal',fraction=0.04, pad = 0.001, shrink=0.80)
    plt.title('SIM - '+titlet, fontdict=font)

    plt.text(0,0,r'BIAS: %.0f ' %(100*np.around((np.nanmean(np.abs(sim-obs))/np.nanmean(obs)),2))+'%', fontdict=font2)
    
    plt.subplot(1, 2, 2)
    y=plt.imshow(obs2, interpolation='none',cmap='Spectral_r')
    plt.axis('off') 
    plt.colorbar(y,orientation ='horizontal',fraction=0.04, pad = 0.001, shrink=0.80)
    plt.title('OBS - '+titlet, fontdict=font)
    
    
    foldersim = './temp_Maps/'
    filenamesim=titlet+'.jpg'
    fig.savefig(foldersim+filenamesim, dpi=60,bbox_inches='tight')

    plt.close(fig)

print('figure saved')
