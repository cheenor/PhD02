#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 13:58:30 2014

@author: Chenjh
"""
import numpy as np
import math
def uv2vordiv(uwnd,vwnd,lon,lat,div,vor):
    rr=6371393   # earth radius in m    
    pi=3.1415926535
    ny=len(lat)
    dy=(2*rr*pi/360.)*abs(lon[1]-lon[0])
    nx=len(lon)
    nt=len(uwnd[:,0,0,0])
    nz=len(uwnd[0,:,0,0])
    dx=np.ndarray(shape=(ny), dtype=float)
    for i in range(0,ny-1):
        dx[i]=(2*pi*rr*math.cos(pi*lat[i]/180)/360.)*abs(lat[1]-lat[0]) 
    for it in range(0,nt-1):
        for iz in range(0,nz-1):
            for ix in range(0,nx-1):
                for iy in range(0,ny-1):
                    if ix==0:
                        dux=uwnd[it,iz,iy,ix+1]-uwnd[it,iz,iy,ix]
                        dvx=vwnd[it,iz,iy,ix+1]-vwnd[it,iz,iy,ix]
                        dyy=dy                        
                    elif ix==nx-1:
                        dux=uwnd[it,iz,iy,ix]-uwnd[it,iz,iy,ix-1]
                        dvx=vwnd[it,iz,iy,ix]-vwnd[it,iz,iy,ix-1] 
                        dyy=-dy
                    else:
                        dux=uwnd[it,iz,iy,ix+1]-uwnd[it,iz,iy,ix-1]
                        dvx=vwnd[it,iz,iy,ix+1]-vwnd[it,iz,iy,ix-1] 
                        dyy=2*dy             
                    if iy==0:
                        dvy=vwnd[it,iz,iy+1,ix]-vwnd[it,iz,iy,ix]
                        duy=uwnd[it,iz,iy+1,ix]-uwnd[it,iz,iy,ix]
                        dxx=dx[iy]       
                    elif iy==ny-1:
                        dvy=vwnd[it,iz,iy,ix]-vwnd[it,iz,iy-1,ix]
                        duy=uwnd[it,iz,iy,ix]-uwnd[it,iz,iy-1,ix]
                        dxx=-dx[iy]
                    else:
                        dvy=vwnd[it,iz,iy+1,ix]-vwnd[it,iz,iy-1,ix]
                        duy=uwnd[it,iz,iy+1,ix]-uwnd[it,iz,iy-1,ix]
                        dxx=2*dx[iy]
                    div[it,iz,iy,ix]=dux/dxx+dvy/dyy
                    vor[it,iz,iy,ix]=dvx/dxx-duy/dyy
                    print div[it,iz,iy,ix]
            
    return div,vor
    