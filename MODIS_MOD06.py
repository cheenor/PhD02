#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 11:36:17 2015

@author: jhchen
"""
import gdal,ogr    ### osheo can't be imported with pyhdf
import osgeo
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
### this is a fuction to obtain the lon and lat in gdal  
def pixel2coord(x, y,geoform):
    """Returns global coordinates from pixel x, y coords"""
    xoff=geoform[0]    
    a=geoform[1]
    b=geoform[2]
    yoff=geoform[3]
    d=geoform[4]
    e=geoform[5]
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)
# xoff, a, b, yoff, d, e = ds.GetGeoTransform()
dirin="Z:/DATA/MODIS/Xiao/"
dirout="Z:/DATA/MODIS/"
fnm="MOD06_L2.A2011170.0015.006.2015054231551.hdf"
fpath=dirin+fnm
modis=gdal.Open(fpath)
#modisinfo=gdalinfo(fpath)
listvar=modis.GetSubDatasets()
fpath=dirout+"MOD06_L2.txt"
i=0
f=open(fpath,'w')
for nm in listvar:
    iband="%d "%i
    f.write(iband)
    itme="%s "%nm[1]
    f.write(itme)
    f.write('\n')
    i+=1
f.close()
lc_data=gdal.Open(listvar[34][0])
latrd1=gdal.Open(listvar[124][0])
lonrd1=gdal.Open(listvar[125][0])
latrd=latrd1.ReadAsArray()
lonrd=lonrd1.ReadAsArray()
st1=gdal.Open(listvar[0][0])
stm1=st1.ReadAsArray()
st2=gdal.Open(listvar[126][0])
stm2=st2.ReadAsArray()
lc=lc_data.ReadAsArray()
nx=lc_data.RasterXSize
ny=lc_data.RasterYSize
geoform=lc_data.GetGeoTransform()
lon=np.ndarray(shape=(nx,ny), dtype=float)
lat=np.ndarray(shape=(nx,ny), dtype=float)
for i in range(0,nx):
    for j in range(0,ny):
        lon[i,j],lat[i,j] =pixel2coord(i,j,geoform)
minlat=int(lat.min())
maxlat=int(lat.max())
minlon=int(lon.min())
maxlon=int(lon.max())
lonlabs=[]
latlabs=[]
for i in range(0,ny):
    itm="%d"%int(lonrd[i,0])
    lonlabs.append(itm)
for i in range(0,nx):
    itm="%d"%int(latrd[0,i])
    latlabs.append(itm)    
fig,ax = plt.subplots(nrows=1,ncols=1)
#ax.set_ylim(35,80)
#ax.set_xlim(250,350)
plt.imshow(lc[:,:])
ax.set_xticks(range(0,ny,100))
xticklabels = [lonlabs[nn] for nn in range(0,ny,100)] 
ax.set_xticklabels(xticklabels, size=14)
ax.set_yticks(range(0,nx,80))
yticklabels = [latlabs[nn] for nn in range(0,nx,80)] 
ax.set_yticklabels(yticklabels, size=14)
plt.colorbar()
plt.show()
 