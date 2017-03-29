# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 18:17:29 2016

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
import matplotlib.cm as cm
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap 
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
###############################################################################
def readAscii2(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    linesplit=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit.append(line[:-1].split(' '))
    for lnstrs in linesplit:
        for strs in lnstrs:
            if strs!='':
                onedim.append(string.atof(strs))
    del linesplit
    f.close()
    return onedim
def readRain(fpath,iskp):
    #iskp  the total line skipped of the file
    # fpath   the full path of the file
    # usage: onedim=readAscii(fpaht,iskp)
    onedim=[]
    ystr=[]
    mstr=[]
    dstr=[]
    linesplit=[]
    datestr=[]
    f=open(fpath)
    ff=f.readlines()[iskp:]  ## first line in obs file is legend 
    for line in ff:
        line=string.lstrip(line)
        linesplit1=line[:-1].split(' ')
        datestr.append(linesplit1[0])
        linesplit.append(linesplit1[1:])
    for lnstrs in linesplit:
        #print lnstrs          
        for strs in lnstrs:
            if strs!='':
                #print strs
                onedim.append(string.atof(strs))
    for lnstrs in datestr:           
        ystr.append(string.atoi(lnstrs[0:4]))
        mstr.append(string.atoi(lnstrs[4:6]))
        dstr.append(string.atoi(lnstrs[6:8]))
    del linesplit1
    f.close()
    return ystr,mstr,dstr,onedim
def ave(tmp,iday):
    nz=len(tmp[0,:,0,0])
    ny=len(tmp[0,0,:,0])
    nx=len(tmp[0,0,0,:])
    print nz,ny,nx
    out=np.zeros(shape=(nz,ny,nx),dtype=float)
    for iz in range(0,nz):
        for ix in range(0,nx):
            for iy in range(0,ny):
                for it in range(iday,iday+1):
                    out[iz,iy,ix]=tmp[it,iz,iy,ix]               
    return out
def getdata(dirnc,iy,im,idy,nz,ny,nx):
    var=['uwnd','vwnd','sh']
    varin=['u','v','q']
    dt=datetime.datetime(iy,im,idy)
    daystr=dt.strftime('%j')
    iday=string.atoi(daystr)-1
    print iday
    out=np.zeros(shape=(len(var),nz,ny,nx),dtype=float)
    for i in range(0,3):
        for j in range(0,len(var)):
            filename='%d'%iy+'_'+var[j]+'_'+'%2.2i'%(i*6)+'.nc'
            fpath=dirnc+filename
            print fpath
            f=Dataset(fpath,'a')
            varname=varin[j]
            tmp=f.variables[varname][:]
            tmp1=ave(tmp,iday)
            out[j,:,:,:]=tmp1[:,:,:]+out[j,:,:,:]            
    return out[0,:,:,:]/4.0,out[1,:,:,:]/4.0,out[2,:,:,:]/4.0
dirin='D:/MyPaper/PhD02/Data/'
dirnc='K:/Data/ERA_interim/X2.5/'
dirout='D:/MyPaper/PhD02/Pics/'
regions=['ETP','WTP']
classes=['Dry','Rain','Normal','Norain','Bothrain']
ncls=len(classes)
###############get the lon and lat --------------------------------------------
fpath=dirnc+'1979_air_00.nc'
f=Dataset(fpath,'a')
lon=f.variables['longitude'][:]
lat=f.variables['latitude'][:]
air=f.variables['t'][:]
nx=len(lon)
ny=len(lat)
fpath=dirout+'MJJAS_'+"patterns_out_lon_lat.txt"
fout=open(fpath,"w")
for ix in range(0,nx):
    item="%f "%(lon[ix])
    fout.write(item)
for iy in range(0,ny):
    item="%f "%(lat[iy])
    fout.write(item)        
nz=len(air[0,:,0,0])
for ig in range(0,len(regions)):
    for icld in range(4,5):#ncls):
        if icld<2:
            fpath=dirin+'MJJAS_'+classes[icld]+'_Only_in_'+regions[ig]+'.txt'
            year,month,day,onedim=readRain(fpath,0)
            nt=int(len(onedim)/8)
            rain=np.zeros(shape=(nt),dtype=float)
            dry=np.zeros(shape=(3,nt),dtype=int)
            wet=np.zeros(shape=(3,nt),dtype=int)
            for i in range(0,nt):
                k=i*8
                rain[i]=onedim[k]
                for j in range(0,3):
                    dry[j,i]=onedim[k+j+2]
                    wet[j,i]=onedim[k+j+5]
        else:
            fpath=dirin+'MJJAS_'+classes[icld]+'_in_ETP+WTP'+'.txt'
            year,month,day,onedim=readRain(fpath,0)
            nt=len(year)
        udry=np.zeros(shape=(nz,ny,nx),dtype=float)
        vdry=np.zeros(shape=(nz,ny,nx),dtype=float)
        qdry=np.zeros(shape=(nz,ny,nx),dtype=float)
        contdry=0
        for i in range(0,nt):
            if year[i]>=1979 and month[i]>4 and month[i]<10:
                iyear=year[i]
                imonth=month[i]
                iday=day[i]
                contdry=contdry+1.
                u,v,q=getdata(dirnc,iyear,imonth,iday,nz,ny,nx)
                for iz in range(0,nz):
                    for ix in range(0,nx):
                        for iy in range(0,ny):            
                            udry[iz,iy,ix]=udry[iz,iy,ix]+u[iz,iy,ix]
                            vdry[iz,iy,ix]=vdry[iz,iy,ix]+v[iz,iy,ix]
                            qdry[iz,iy,ix]=qdry[iz,iy,ix]+q[iz,iy,ix]      
###############################################################################
        if icld<2:
            fpath=dirout+'MJJAS_'+classes[icld]+'_Only_in_'+regions[ig]+"_patterns_out_u_v_q.txt"
        else:
            fpath=dirout+'MJJAS_'+classes[icld]+"_ETP+WTP_patterns_out_u_v_q.txt"
        foutdry=open(fpath,"w")
        for iz in range(0,nz):
            for ix in range(0,nx):
                for iy in range(0,ny):                                  
                    udry[iz,iy,ix]=udry[iz,iy,ix]/contdry
                    item="%e "%(udry[iz,iy,ix])
                    foutdry.write(item)
                    vdry[iz,iy,ix]=vdry[iz,iy,ix]/contdry
                    item="%e "%(vdry[iz,iy,ix])
                    foutdry.write(item)
                    qdry[iz,iy,ix]=qdry[iz,iy,ix]/contdry
                    item="%e "%(qdry[iz,iy,ix]*1000)
                    foutdry.write(item)
                    foutdry.write("\n")   
        foutdry.close()