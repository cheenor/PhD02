#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 03 10:51:23 2016

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import math
import numpy as np
import matplotlib.cm as cm
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap 
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
###############################################################################
def readAscii(fpath,iskp):
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
def getmeanps():
    dirps='X:/Data/ERA_interim/SRFX2.5/'
    for ih in range(0,4):
        hourstr='%2.2d'%(ih*6)
        fpath=dirps+'1979-2014_ps_'+hourstr+'.nc'
        #print fpath
        f=Dataset(fpath,'a')
        tmp4=f.variables['sp'][:]
        #print tmp4.shape
        nx=len(tmp4[0,0,:])
        ny=len(tmp4[0,:,0])
        if ih==0:
            tmp=np.zeros(shape=(ny,nx),dtype=float)
        for ix in range(0,nx):
            for iy in range(0,ny):
                tmp2=tmp4[:,iy,ix]
                tmp[iy,ix]=tmp[iy,ix]+(tmp2.mean())/4.
    return tmp 
def myplot(dirout,lon,lat,rgname,typstr,datatype,u,v,data,iz,plv):
    nx=len(lon)
    ny=len(lat)
    nz=len(plv)
    ps=getmeanps() 
    for ix in range(0,nx):
        for iy in range(0,ny):
            #print ps[iy,ix]
            if plv[nz-1-iz]>=(ps[iy,ix]/100.):
                u[iz,iy,ix]=None
                v[iz,iy,ix]=None
                data[iz,iy,ix]=None                
    lon_grid=np.ndarray(shape=(ny,nx),dtype=float)
    lat_grid=np.ndarray(shape=(ny,nx),dtype=float)
    if plv[nz-1-iz]<10: 
        levstr='%1d'%plv[nz-1-iz]+'hPa'
    elif plv[nz-1-iz]<100: 
        levstr='%2.2d'%plv[nz-1-iz]+'hPa'
    elif plv[nz-1-iz]<1000: 
        levstr='%3.3d'%plv[nz-1-iz]+'hPa'
    else: 
        levstr='%4.4d'%plv[nz-1-iz]+'hPa'    
    if datatype=='q':        
        if typstr[0:4]=='diff':
            mycolors=cm.PiYG
            myscale=2
            mylevs=[-1.5,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,
                    0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.5]
        else: 
            mycolors=cm.Greens
            myscale=10
            if iz>27:
                mylevs=[2,4,6,8,10,12,14,16,18,20]
            elif iz>23:
                mylevs=[1,2,4,5,6,7,8,9,10]
            elif iz>19:
                mylevs=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
            else:
                mylevs=[0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3]
        #mylevs=np.linspace(0.5,20.5,39)   
    else:
        mycolors=cm.gist_rainbow
    for i in range(0,ny):
        for j in range(0,nx):
            lon_grid[i,j]=lon[j]
            lat_grid[i,j]=lat[i]   
    fig = plt.figure(figsize=(12,9)) ### figure 8"X8" inches，图片大小
    m=Basemap(projection='lcc', llcrnrlon=65.,llcrnrlat=5.,urcrnrlon=145.,urcrnrlat=50.,
          lat_1=15.,lat_2=45.,lat_0=35.,lon_0=105.,rsphere=6371200., resolution='l') 
    m.drawcoastlines() # 显示地图海岸线
    m.drawcountries() # 显示国界线
    m.drawstates() # 显示国界线
    parallels = np.arange(-90,90,10.) #地图上的纬度线，间隔是10，范围是-90到90
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14) #地图上画出
    # draw meridians
    meridians = np.arange(0.,360.,10.) #地图上的经度线，间隔是10，范围是-90到90
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14) # 地图上画出
    x, y = m(lon_grid, lat_grid) # compute map proj coordinates.
    cs = m.contourf(x,y,data[iz,:,:],levels=mylevs,cmap=mycolors,extend='both')
    #cbar = m.colorbar(cs,location='bottom',pad="5%") # 色标
    #cbar.set_label('($g$ $kg^{-3}$)',fontsize=16) # 色标下面标注单位  
    Q = m.quiver(x,y,u[iz,:,:],v[iz,:,:],units='inches',scale=myscale)
    qk = plt.quiverkey(Q, 0.98, -0.12, myscale,'%d'%myscale+r' $m$ $s^{-1}$',
                   fontproperties={'weight': 'bold'})
    if typstr=='diff':
        titstr='Differences of wind and moisture between the wet and dry conditions at '+levstr
    else:
        titstr='Averaged wind and moisture of the '+typstr+' condition '+'at '+levstr
    plt.title(titstr,fontsize=18)
    plt.subplots_adjust(left = 0.1,right=0.9, bottom = 0.12, top = 0.96)    
    cax = fig.add_axes([0.1, 0.085, 0.8, 0.035])
    cbar=fig.colorbar(cs, cax,extend='both',
             spacing='uniform', orientation='horizontal') 
    cbar.set_label('($g$ $kg^{-3}$)',fontsize=18) # 色标下面标注单位          
    plt.show() #显示图像
    plt.savefig(dirout+rgname+typstr+'_pattern_'+levstr+'.png',dpi=300)  #保存图片，这里保存为pdf，常用格式矢量格式 eps,ps等都是支持的
    # 非矢量格式，png，jpg，tiff等
    plt.close() # 关闭画图窗口
    return 
###############################################################################
dirnc='X:/Data/ERA_interim/X2.5/'
dirin='D:/MyPaper/PhD02/Pics/'
dirout='D:/MyPaper/PhD02/Pics/Patterns/'
regions=['ETP','WTP']
ng=len(regions)
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
nz=len(air[0,:,0,0])
plv=[1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.
     ,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.
     ,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.
     ,10.,7.,5.,3.,2.,1.]
uwd=np.zeros(shape=(2,2,nz,ny,nx),dtype=float)
vwd=np.zeros(shape=(2,2,nz,ny,nx),dtype=float)
qwd=np.zeros(shape=(2,2,nz,ny,nx),dtype=float)
ubth=np.zeros(shape=(3,nz,ny,nx),dtype=float)
vbth=np.zeros(shape=(3,nz,ny,nx),dtype=float)
qbth=np.zeros(shape=(3,nz,ny,nx),dtype=float)
for ip in range(0,ncls):
    if ip<2:
        for ig in range(0,ng):
            filename=classes[ip]+'_Only_in_'+regions[ig]+"_patterns_out_u_v_q.txt"
            fpath=dirin+filename
            onedim=readAscii(fpath,0)
            for iz in range(0,nz):
                for ix in range(0,nx):
                    for iy in range(0,ny):
                        kk=iz*nx*ny*3+ix*ny*3+iy*3                               
                        uwd[ig,ip,iz,iy,ix]=onedim[kk]
                        vwd[ig,ip,iz,iy,ix]=onedim[kk+1]
                        qwd[ig,ip,iz,iy,ix]=onedim[kk+2] # kg/kg to g/kg
            del onedim
            for iz in range(0,nz):
                myplot(dirout,lon,lat,regions[ig],' only '+classes[ip],'q', 
                    uwd[ig,ip,:,:,:],vwd[ig,ip,:,:,:],qwd[ig,ip,:,:,:],iz,plv)  
    else:
        filename=classes[ip]+"_ETP+WTP_patterns_out_u_v_q.txt"
        fpath=dirin+filename
        onedim=readAscii(fpath,0)
        for iz in range(0,nz):
            for ix in range(0,nx):
                for iy in range(0,ny):
                    kk=iz*nx*ny*3+ix*ny*3+iy*3                               
                    ubth[ip-2,iz,iy,ix]=onedim[kk]
                    vbth[ip-2,iz,iy,ix]=onedim[kk+1]
                    qbth[ip-2,iz,iy,ix]=onedim[kk+2] # kg/kg to g/kg
        del onedim
        for iz in range(0,nz):
            myplot(dirout,lon,lat,classes[ip],'both_ETP+WTP','q',
                ubth[ip-2,:,:,:],vbth[ip-2,:,:,:],qbth[ip-2,:,:,:],iz,plv)    
dfu=ubth[2,:,:,:]-ubth[1,:,:,:]
dfv=vbth[2,:,:,:]-vbth[1,:,:,:]
dfq=qbth[2,:,:,:]-qbth[1,:,:,:]    
for iz in range(0,nz):
    myplot(dirout,lon,lat,'ETP+WTP','diff_BOTH','q',dfu,dfv,dfq,iz,plv)
for ig in range(0,ng):
    dfu=uwd[ig,1,:,:,:]-uwd[ig,0,:,:,:]
    dfv=vwd[ig,1,:,:,:]-vwd[ig,0,:,:,:]
    dfq=qwd[ig,1,:,:,:]-qwd[ig,0,:,:,:]    
    for iz in range(0,nz):
        myplot(dirout,lon,lat,regions[ig],'diff_Only','q',dfu,dfv,dfq,iz,plv)        
    ###############################################################################
