#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 10:58:46 2014

@author: Chenjh
"""
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import uv2vrdvf as uv2dv #  this is a fuction written by Chen
import string
#####
whereisdata="Z:/DATA/LargeScale/NCEPR2/"
whereiseventdate='D:/MyPaper/PhD02/Data/'
outpic='D:/MyPaper/PhD02/PICSNEW/'
# the average period for normal
its=4 # May
ite=8 # Sep.
nzr=12  # there is 17 levels, but Do you want draw them all? Here, level 1 to 12

#open and read datasets
hgtfile=Dataset(whereisdata+"hgt.sfc.nc",'a') 
lon_hgt= hgtfile.variables['lon'][:]   # read 
tmp_hgt=hgtfile.variables['time'][:]
lat_hgt=hgtfile.variables['lat'][:]
height_all=hgtfile.variables['hgt'][:]  # read n
height=height_all[0,:,:]
del height_all
#
pss=Dataset(whereisdata+"pres.ps.sea.nc",'a') 
lon_ps= pss.variables['lon'][:]   # read 
tmp_ps=pss.variables['time'][:]
lat_ps=pss.variables['lat'][:]
ps_all=pss.variables['pres'][:]  # read 
ps=ps_all[2,:,:]
del ps_all
#
uwdf=Dataset(whereisdata+"uwnd.ymon.nc",'a') 
lon_uwd= uwdf.variables['lon'][:]   # read 
tmp_uwd=uwdf.variables['time'][:]
lat_uwd=uwdf.variables['lat'][:]
uwd=uwdf.variables['uwnd'][:]  # read 
#
vwdf=Dataset(whereisdata+"vwnd.ymon.nc",'a') 
lon_vwd= vwdf.variables['lon'][:]   # read 
tmp_vwd=vwdf.variables['time'][:]
lat_vwd=vwdf.variables['lat'][:]
vwd=vwdf.variables['vwnd'][:]  # read
#
airf=Dataset(whereisdata+"air.ymon.nc",'a') 
lon_air= airf.variables['lon'][:]   # read 
tmp_air=airf.variables['time'][:]
lat_air=airf.variables['lat'][:]
air=airf.variables['air'][:]  # read 
#
rhmf=Dataset(whereisdata+"rhum.ymon.nc",'a') 
lon_rh=rhmf.variables['lon'][:]   # read 
tmp_rh=rhmf.variables['time'][:]
lat_rh=rhmf.variables['lat'][:]
lev_rh=rhmf.variables['level'][:]
rh=rhmf.variables['rhum'][:]  # read 
# 
ny=len(lat_uwd)
nx=len(lon_uwd)
nt=len(uwd[:,0,0,0])
nz=len(uwd[0,:,0,0])
div=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)
vor=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)     
div,vor=uv2dv.uv2vordiv(uwd,vwd,lon_uwd,lat_uwd,div,vor) 
# the average period
divs=np.ndarray(shape=(nz,ny,nx), dtype=float)
vors=np.ndarray(shape=(nz,ny,nx), dtype=float) 
uwds=np.ndarray(shape=(nz,ny,nx), dtype=float) 
vwds=np.ndarray(shape=(nz,ny,nx), dtype=float) 
airs=np.ndarray(shape=(nz,ny,nx), dtype=float) 
rhms=np.ndarray(shape=(nz,ny,nx), dtype=float)  
for im in range(its,ite+1):
    for iz in range(0,nz-1):
        for iy in range(0,ny-1):
            for ix in range(0,nx-1):
                divs[iz,iy,ix]=divs[iz,iy,ix]+div[im,iz,iy,ix]
                vors[iz,iy,ix]=vors[iz,iy,ix]+vor[im,iz,iy,ix]
                uwds[iz,iy,ix]=uwds[iz,iy,ix]+uwd[im,iz,iy,ix]
                vwds[iz,iy,ix]=vwds[iz,iy,ix]+vwd[im,iz,iy,ix]
                airs[iz,iy,ix]=airs[iz,iy,ix]+air[im,iz,iy,ix]
                rhms[iz,iy,ix]=rhms[iz,iy,ix]+rh[im,iz,iy,ix]
divs=divs/(ite-its+1)
vors=vors/(ite-its+1)
uwds=uwds/(ite-its+1)
vwds=vwds/(ite-its+1)
airs=airs/(ite-its+1)
rhms=rhms/(ite-its+1)

for iz in range(0,nzr):
    lev="%04d"%(lev_rh[iz])
# start plotting
# what do you want
# fig1  wind temperature and rh
# fig2 div+vor    
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=10., urcrnrlon=140., urcrnrlat=55.,\
               lat_1=25.,lat_2=40., lon_0=105,  projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    lons,lats=np.meshgrid(lon_uwd,lat_uwd)
    x, y = m(lons,lats)
    clevs_air = np.arange(200,320,2)
    clevs_rh = np.arange(10,100,10)
    
    cs1 = m.contourf(x,y,airs[iz,:,:],clevs_air,cmap=cm.jet,animated=True)
    cs2 = m.contour(x,y,rhms[iz,:,:],clevs_rh,linewidths=1.5,\
                     colors='k')  #,animated=True)
    plt.clabel(cs2,inline=1,fmt='%1.1f',fontsize=12) 
# plot wind vectors on projection grid.
# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
    print lat_uwd[0]
    lat_rev=lat_uwd[::-1] 
    print lat_rev[72]
    print uwds[iz,0,0]
    uwdsx=uwds[iz,::-1,:]
    print uwdsx[72,0]
    vwdsx=vwds[iz,::-1,:]
    ugrid,newlons = shiftgrid(180.,uwdsx,lon_uwd,start=False)
    vgrid,newlons = shiftgrid(180.,vwdsx,lon_uwd,start=False)
# transform vectors to projection grid.       
    uproj,vproj,xx,yy = \
    m.transform_vector(ugrid,vgrid,newlons,lat_rev,40,30,returnxy=True,masked=True)
# now plot.
    Q = m.quiver(xx,yy,uproj,vproj,color='w')
# make quiver key.
    qk = plt.quiverkey(Q, 0.1, 0.1, 10, '10 m/s') #, labelpos='W')
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)
#    
    cb = m.colorbar(cs1) #,"bottom", size="5%", pad="2%")
    cb.set_label('Temperature (K)',fontsize=14)
    # set plot title
    ax.set_title('Temp., Rh and Wind Field at '+lev+'hPa')
    plt.show()
    plt.savefig(outpic+'Temp+Rh+Wind'+lev+'.pdf')
    plt.close()
# now figure 2 div+ vor   
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=10., urcrnrlon=140., urcrnrlat=55.,\
               lat_1=25.,lat_2=40., lon_0=105,  projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    lons,lats=np.meshgrid(lon_uwd,lat_uwd)
    x, y = m(lons,lats)
    clevs_div = np.arange(-3,3,0.1)
    clevs_vor = np.arange(-3,3,0.4)
#    clevs_ht=[3000.]
    cs1 = m.contourf(x,y,divs[iz,:,:]*100000.,clevs_div,cmap=cm.seismic,animated=True)
    cs2 = m.contour(x,y,vors[iz,:,:]*100000.,clevs_vor,linewidths=1.5,\
                     colors='k',animated=True)
    plt.clabel(cs2,inline=1,fmt='%1.1f',fontsize=12) 
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)
#    
    cb = m.colorbar(cs1) #,"bottom", size="5%", pad="2%")
    cb.set_label(r'$Disvergence 10^6 S^{-1}$', fontsize=14)
    # set plot title
    ax.set_title('Divergence and Vorticity at '+lev)
    plt.show()
    plt.savefig(outpic+'div+vor'+lev+'.pdf')
    plt.close()
del div,vor,uwd,vwd,air,rh
#  end plot of the comment condition
# now we start to hand the events    
### ETP first
filepath=whereiseventdate+'ETP_EventsDate_cloudsat.txt'
f=open(filepath)
evendate_onedim=[]
avg_rain_onedim=[]
linesplit=[]
ff=f.readlines()[1:]
nl=len(ff)
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            evendate_onedim.append(string.atof(strs))
            avg_rain_onedim.append(string.atof(strs))
eventdate=np.ndarray(shape=(nl,3), dtype=int) 
avg_rain =np.ndarray(shape=(nl), dtype=float) 
for i in range(0,nl-1):
    for j in range(0,3):
        ij=i*3+j
        eventdate[i,j]=int(evendate_onedim[ij])
    avg_rain[i]=avg_rain_onedim[i*3+2]
######################################################333
days=[31,28,31,30,31,30,40,30,30,31,30,31]
uwnd_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
vwnd_sum=np.ndarray(shape=(nz,ny,nx),dtype=float) 
air_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
rh_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
vor_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
div_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
isum=0
for i in range(0,nl-1):
    if eventdate[i,0]>1978:
        yearstr="%04d"%(eventdate[i,0])
        evn_uwd=Dataset(whereisdata+"uwnd."+yearstr+".nc","a")
        evn_vwd=Dataset(whereisdata+"vwnd."+yearstr+".nc","a")
        evn_air=Dataset(whereisdata+"air."+yearstr+".nc","a")
        evn_rh=Dataset(whereisdata+"rhum."+yearstr+".nc","a")
#------- read
        lon_evt= evn_uwd.variables['lon'][:]   # read 
        tmp_evt= evn_uwd.variables['time'][:]
        lat_evt= evn_uwd.variables['lat'][:]
        lev_evt= evn_uwd.variables['level'][:]
        uwnd_evt=evn_uwd.variables['uwnd'][:]  # read 
        vwnd_evt=evn_vwd.variables['vwnd'][:] 
        air_evt=evn_air.variables['air'][:] 
        rh_evt=evn_rh.variables['rhum'][:] 
        nt=len(uwnd_evt[:,0,0,0])
        nz=len(uwnd_evt[0,:,0,0])
        div_evt=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)
        vor_evt=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)     
        uv2dv.uv2vordiv(uwnd_evt,vwnd_evt,lon_evt,lat_evt,vor_evt,div_evt)
#######################################################3
        days[1]=28
        if eventdate[i,0] % 4 == 0 and eventdate[i,0] %100 != 0 \
            or eventdate[i,0] % 400 == 0:    
            days[1]=29
        idex=0
        for im in range(0,eventdate[i,1]-2):
            idex=idex+days[im]*4
        idex=idex+(eventdate[i,2]-1)*4
        for ix in range(idex+0,idex+4): 
            uwnd_sum=uwnd_sum+ uwnd_evt[ix,:,:,:] # read 
            vwnd_sum=vwnd_sum+ vwnd_evt[ix,:,:,:]
            air_sum=air_sum+ air_evt[ix,:,:,:]
            rh_sum=rh_sum+ rh_evt[ix,:,:,:]
            vor_sum=vor_sum+ vor_evt[ix,:,:,:]
            div_sum=div_sum+ div_evt[ix,:,:,:]
            isum=isum+1
        del div_evt,vor_evt,lon_evt,rh_evt # read 
        del tmp_evt,lat_evt,lev_evt,uwnd_evt,vwnd_evt,air_evt
uwnd_sum=uwnd_sum/isum 
vwnd_sum=vwnd_sum/isum 
air_sum=air_sum/isum 
rh_sum=rh_sum/isum 
vor_sum=vor_sum/isum 
div_sum=div_sum/isum  
#### plot evnets of ETP
for iz in range(0,nzr):
    lev="%04d"%(lev_evt[iz])
# start plotting
# what do you want
# fig1  wind temperature and rh
# fig2 div+vor    
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=10., urcrnrlon=140., urcrnrlat=55.,\
               lat_1=25.,lat_2=40., lon_0=105,  projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    lons,lats=np.meshgrid(lon_evt,lat_evt)
    x, y = m(lons,lats)
    clevs_air = np.arange(220,320,2.5)
    clevs_rh = np.arange(10,100,10)
#    clevs_ht=[3000.]
    cs1 = m.contourf(x,y,air_sum[iz,:,:],clevs_air,\
                    cmap=cm.jet,animated=True)
    cs2 = m.contour(x,y,rh_sum[iz,:,:],clevs_rh,linewidths=1.5,\
                     colors='k',animated=True)
# plot wind vectors on projection grid.
# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
    ugrid,newlons = shiftgrid(180.,uwnd_sum,lon_evt,start=False)
    vgrid,newlons = shiftgrid(180.,vwnd_sum,lon_evt,start=False)
# transform vectors to projection grid.
    uproj,vproj,xx,yy = \
    m.transform_vector(ugrid,vgrid,newlons,lat_evt,40,30,\
                        returnxy=True,masked=True)
# now plot.
    Q = m.quiver(xx,yy,uproj,vproj,scale=700)
# make quiver key.
    qk = plt.quiverkey(Q, 0.1, 0.1, 10, '10 m/s', labelpos='W')
#    
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)
    cb = m.colorbar(cs1,"bottom", size="5%", pad="2%")
    cb.set_label('K')
    # set plot title
    ax.set_title('Temp., Rh and Wind fild at '+lev)
    plt.show()
    plt.savefig(outpic+'ETP_Temp+Rh+Wind'+lev+'_Evnt.pdf')
    plt.close()
# now figure 2 div+ vor   
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=10., urcrnrlon=140., urcrnrlat=55.,\
               lat_1=25.,lat_2=40., lon_0=105,  projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    x, y = m(lon_uwd,lat_uwd)
    clevs_div = np.arange(-3,3,0.1)
    clevs_vor = np.arange(-3,3,0.4)
    clevs_ht=[3000.]
    cs1 = m.contourf(x,y,div_sum[iz,:,:]*100000.,clevs_div,\
                    cmap=cm.seismic,animated=True)
    cs2 = m.contour(x,y,vor_sum[iz,:,:]*100000.,clevs_vor,linewidths=1.5,\
                     colors='k',animated=True)
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)
#    
    cb = m.colorbar(cs1,"bottom", size="5%", pad="2%")
    cb.set_label(r'$10^6 S^-^1$')
    # set plot title
    ax.set_title('Divergence and Vorticity at '+lev)
    plt.show()
    plt.savefig(outpic+'ETP_div+vor'+lev+'Evnt.pdf')
    plt.close()  
### Then WTP
filepath=whereiseventdate+'WTP_EventsDate_cloudsat.txt'
f=open(filepath)
evendate_onedim=[]
avg_rain_onedim=[]
linesplit=[]
ff=f.readlines()[1:]
nl=len(ff)
for line in ff:
    line=string.lstrip(line)
    linesplit.append(line[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            evendate_onedim.append(string.atof(strs))
            avg_rain_onedim.append(string.atof(strs))
eventdate=np.ndarray(shape=(nl,3), dtype=int) 
avg_rain =np.ndarray(shape=(nl), dtype=float) 
for i in range(0,nl-1):
    for j in range(0,3):
        ij=i*3+j
        eventdate[i,j]=int(evendate_onedim[ij])
    avg_rain[i]=avg_rain_onedim[i*3+2]
######################################################333
days=[31,28,31,30,31,30,40,30,30,31,30,31]
uwnd_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
vwnd_sum=np.ndarray(shape=(nz,ny,nx),dtype=float) 
air_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
rh_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
vor_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
div_sum=np.ndarray(shape=(nz,ny,nx),dtype=float)
isum=0
for i in range(0,nl-1):
    if eventdate[i,0]>1978:
        yearstr="%04d"%(eventdate[i,0])
        evn_uwd=Dataset(whereisdata+"uwnd."+yearstr+".nc","a")
        evn_vwd=Dataset(whereisdata+"vwnd."+yearstr+".nc","a")
        evn_air=Dataset(whereisdata+"air."+yearstr+".nc","a")
        evn_rh=Dataset(whereisdata+"rhum."+yearstr+".nc","a")
#------- read
        lon_evt= evn_uwd.variables['lon'][:]   # read 
        tmp_evt= evn_uwd.variables['time'][:]
        lat_evt= evn_uwd.variables['lat'][:]
        lev_evt= evn_uwd.variables['level'][:]
        uwnd_evt=evn_uwd.variables['uwnd'][:]  # read 
        vwnd_evt=evn_vwd.variables['vwnd'][:] 
        air_evt=evn_air.variables['air'][:] 
        rh_evt=evn_rh.variables['rhum'][:] 
        nt=len(uwnd_evt[:,0,0,0])
        nz=len(uwnd_evt[0,:,0,0])
        div_evt=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)
        vor_evt=np.ndarray(shape=(nt,nz,ny,nx), dtype=float)     
        uv2dv.uv2vordiv(uwnd_evt,vwnd_evt,lon_evt,lat_evt,vor_evt,div_evt)
#######################################################3
        days[1]=28
        if eventdate[i,0] % 4 == 0 and eventdate[i,0] %100 != 0 \
            or eventdate[i,0] % 400 == 0:    
            days[1]=29
        idex=0
        for im in range(0,eventdate[i,1]-2):
            idex=idex+days[im]*4
        idex=idex+(eventdate[i,2]-1)*4
        for ix in range(idex+0,idex+4): 
            uwnd_sum=uwnd_sum+ uwnd_evt[ix,:,:,:] # read 
            vwnd_sum=vwnd_sum+ vwnd_evt[ix,:,:,:]
            air_sum=air_sum+ air_evt[ix,:,:,:]
            rh_sum=rh_sum+ rh_evt[ix,:,:,:]
            vor_sum=vor_sum+ vor_evt[ix,:,:,:]
            div_sum=div_sum+ div_evt[ix,:,:,:]
            isum=isum+1
        del div_evt,vor_evt,lon_evt,rh_evt # read 
        del tmp_evt,lat_evt,lev_evt,uwnd_evt,vwnd_evt,air_evt
uwnd_sum=uwnd_sum/isum 
vwnd_sum=vwnd_sum/isum 
air_sum=air_sum/isum 
rh_sum=rh_sum/isum 
vor_sum=vor_sum/isum 
div_sum=div_sum/isum  
#### plot evnets of ETP
for iz in range(0,nzr):
    lev="%04d"%(lev_evt[iz])
# start plotting
# what do you want
# fig1  wind temperature and rh
# fig2 div+vor    
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=10., urcrnrlon=140., urcrnrlat=55.,\
               lat_1=25.,lat_2=40., lon_0=105,  projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    lons,lats=np.meshgrid(lon_evt,lat_evt)
    x, y = m(lons,lats)
    clevs_air = np.arange(230,320,2)
    clevs_rh = np.arange(10,100,10)
#    clevs_ht=[3000.]
    cs1 = m.contourf(x,y,air_sum[iz,:,:],clevs_air,\
                    cmap=cm.jet,animated=True)
    cs2 = m.contour(x,y,rh_sum[iz,:,:],clevs_rh,linewidths=1.5,\
                     colors='k',animated=True)
# plot wind vectors on projection grid.
# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
    ugrid,newlons = shiftgrid(180.,uwnd_sum,lon_evt,start=False)
    vgrid,newlons = shiftgrid(180.,vwnd_sum,lon_evt,start=False)
# transform vectors to projection grid.
    uproj,vproj,xx,yy = \
    m.transform_vector(ugrid,vgrid,newlons,lat_evt,40,30,\
                        returnxy=True,masked=True)
# now plot.
    Q = m.quiver(xx,yy,uproj,vproj,scale=700)
# make quiver key.
    qk = plt.quiverkey(Q, 0.1, 0.1, 10, '10 m/s', labelpos='W')
#
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)    
    cb = m.colorbar(cs1) #,"bottom", size="5%", pad="2%")
    cb.set_label('K')
    # set plot title
    ax.set_title('Temp., Rh and Wind fild at '+lev)
    plt.show()
    plt.savefig(outpic+'WTP_Temp+Rh+Wind'+lev+'_Evnt.pdf')
    plt.close()
# now figure 2 div+ vor   
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    m=Basemap(llcrnrlon=70., llcrnrlat=15., urcrnrlon=135., urcrnrlat=45,\
    projection='lcc', resolution='l')
    m.drawcoastlines() #
    m.drawcountries() #
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14)
# draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14)
    x, y = m(lon_uwd,lat_uwd)
    clevs_div = np.arange(-3,3,0.1)
    clevs_vor = np.arange(-3,3,0.4)
    clevs_ht=[3000.]
    cs1 = m.contourf(x,y,div_sum[iz,:,:]*100000.,clevs_div,\
                    cmap=cm.seismic,animated=True)
    cs2 = m.contour(x,y,vor_sum[iz,:,:]*100000.,clevs_vor,linewidths=1.5,\
                     colors='k',animated=True)
    if iz<5:
        clevs_ht=[3000.,12000] 
        color_ht=['lightgrey','lightgrey','lightgrey']                
        cs3 = m.contourf(x,y,height,clevs_ht, alpha=1.0,\
                     colors=color_ht,antialiased=True,animated=True)
    else:
        clevs_ht=[3000.]
        color_ht=['lightgrey']
        cs3 = m.contour(x,y,height,clevs_ht,linewidths=3.5,\
                     colors= color_ht,animated=True)
#    
    cb = m.colorbar(cs1) #,"bottom", size="5%", pad="2%")
    cb.set_label(r'$10^6 S^-^1$')
    # set plot title
    ax.set_title('Divergence and Vorticity at '+lev)
    plt.show()
    plt.savefig(outpic+'WTP_div+vor'+lev+'Evnt.pdf')
    plt.close()
















































 