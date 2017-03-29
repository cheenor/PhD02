#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 08:53:23 2014

@author: Chenjh
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import collections
import math
import string
dirin='D:/MyPaper/PhD02/Data/'
dirout='D:/MyPaper/PhD02/PICSNEW/'
fnms=['GEO','GEO_LD','2B-CLDCLASS','2B-TAU','CPR5','ECMWF-AUX','CWC_RVOD',
        ]    #### MODIS 
v_in_file=[3,3,3,3,3,4,17]
nray=125 ### the number of vertial levels






lables_on_figure=['Radar Reflectivity (dBZe)','Cloud Fraction (%)',
                  'Cloud Scenario','Layer Optical Depth',
                  'Radar Reflectivity (dBZe)']
rgns=['ETP','WTP']
nr=2
pbnm='EventsProfile_cloudsat_'

###### begin open and reading
nf=len(fnms) 
for nm in range(0,nf-1):
    filenm=dirin+pbnm+fnms[nm]+'.txt'
    fopen=open(filenm)
    ff=fopen.readline()
    lineraw=[]
    onedim=[]
    vrnms=[]
    for line in ff:
        line=string.lstrip(line)
        lineraw.append(line.split(' '))
    if nm < 5: 
        for strs in lineraw[0]:
            if strs != '':
                vrnms.append(strs)
        for lnstr in lineraw[1:]:
            for strs in lnstr:
                if strs != '':
                    onedim.append(string.atof(strs))
    if nm > (5-1): 
        for strs in lineraw[1]:
            if strs != '':
                vrnms.append(strs)
        for lnstr in lineraw[2:2+nray-1]:
            for strs in lnstr:
                if strs != '':
                    onedim.append(string.atof(strs))
        for lnstr in lineraw[2+nray+1:]:
            for strs in lnstr:
                if strs != '':
                    onedim.append(string.atof(strs))   
    del lineraw
#### the target data [z,vars]     
    nv=v_in_file[nm]

    if nm <5:
        dataset=np.ndarray(shape=(nray,nv),dtype=float)
        for iz in range(0,nray-1):
            for iv in range(0,nv-1):
                idm=iv+iz*nv
                dataset[iz,iv]=onedim[idm]                    
    if nm >(5-1):
        dataset=np.ndarray(shape=(nray,nv,nr),dtype=float)        
        for iv in range(0,nv-1):
            for iz in range(0,nray-1):
                for ir in range(0,nr-1):
                    idm=ir+iz*nr+iv*(nv*nr)
                    dataset[iz,iv,ir]=onedim[idm]
    del onedim
###### plotting
    if nm <5:                    
        fig,ax = plt.figure(figsize=(8,10))
#    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
        colors=cm.gist_ncar
        drdata=[]
        curve = list(zip(dataset[iz,1:], dataset[iz,0]/1000.))
        drdata.append(curve)
        col = collections.LineCollection(drdata)
        ax.add_collection(col, autolim=True)
        col.set_color(colors)    
        ax.add_collection(col, autolim=True)
        ax.autoscale_view()
#    ax.set_title('Successive data offsets')
        ax.set_xlabel(lables_on_figure[nm])
        ax.set_ylabel('Height (km)')
        plt.legend((rgns[0],rgns[1]), loc='upper right')
        plt.show()
        plt.savefig(dirout+fnms[nm]+'.pdf')
        plt.close()
    if fnms[nm] == 'ECMWF-AUX': #### for ECMWF
        fig,ax = plt.figure(figsize=(8,10))
#    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
        colors=cm.gist_ncar        
        for iv in range(1,nv-1):
            drdata=[]
            curve = list(zip(dataset[iz,iv,:], dataset[iz,0,0]/1000.))
            drdata.append(curve)
            fig,ax = plt.figure(figsize=(8,10))
#    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
            colors=cm.gist_ncar
            col = collections.LineCollection(drdata)
            ax.add_collection(col, autolim=True)
            col.set_color(colors)    
            ax.add_collection(col, autolim=True)
            ax.autoscale_view()
#    ax.set_title('Successive data offsets')
            ax.set_xlabel(lables_on_figure[nm])
            ax.set_ylabel('Height (km)')
            plt.legend((rgns[0],rgns[1]), loc='upper right')
            plt.show()
            plt.savefig(dirout+fnms[nm]+vrnms[iv]+'.png',dpi=450)
            plt.close()
    if fnms[nm] == 'CWC_RVOD':
        fig,ax = plt.figure(figsize=(8,10))
#    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
        colors=cm.gist_ncar        
        for iv in range(1,nv-1):
            drdata=[]
            curve = list(zip(dataset[iz,iv,:], dataset[iz,0,0]/1000.))
            drdata.append(curve)
            fig,ax = plt.figure(figsize=(8,10))
#    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
            colors=cm.gist_ncar
            col = collections.LineCollection(drdata)
            ax.add_collection(col, autolim=True)
            col.set_color(colors)    
            ax.add_collection(col, autolim=True)
            ax.autoscale_view()
#    ax.set_title('Successive data offsets')
            ax.set_xlabel(lables_on_figure[nm])
            ax.set_ylabel('Height (km)')
            plt.legend((rgns[0],rgns[1]), loc='upper right')
            plt.show()
            plt.savefig(dirout+fnms[nm]+vrnms[iv]+'.png',dpi=450)
            plt.close()
        