#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 17:22:07 2015

@author: jhchen
"""
import h5py as hdf
from netCDF4 import Dataset
#from pyhdf.SD import SD,SDC
import gdal,ogr    ### osheo can't be imported with pyhdf
from osgeo.gdalnumeric import *  
from osgeo.gdalconst import *
import osgeo 
import numpy as np
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
fpath="E:/Data/2010_2/2010121055002_21318_CS_2B-FLXHR-LIDAR_GRANULE_P2_R04_E03.hdf"
toaobs=gdal.Open(fpath)
listvar=toaobs.GetSubDatasets()
#toaobs.GetSubDatasets()
#print toaobs.GetSubDatasets()
#print listvar[1]
lc_data=gdal.Open(listvar[29][0])
lc=lc_data.ReadAsArray()
fpath="D:/Leanring/gfsanl_4_20100622_0000_000.grb2"
gb2=gdal.Open(fpath)
geo=gb2.GetProjection()
lon1=gb2.GetRasterBand(4)
data1=BandReadAsArray(lon1)
lon2=gb2.GetRasterBand(10)
data2=BandReadAsArray(lon2)
geoform=gb2.GetGeoTransform()
nx=lon2.XSize
ny=lon2.YSize
lon=np.ndarray(shape=(nx,ny), dtype=float)
lat=np.ndarray(shape=(nx,ny), dtype=float)
for i in range(0,nx):
    for j in range(0,ny):
        print pixel2coord(i,j,geoform)
"""
f=open('D:/Leanring/halfdegree_grb2.txt','w')
iband="%s "%'BandNum'
f.write(iband)
itme="%s "%'GRIB_COMMENT'
f.write(itme)
itme="%s "%'GRIB_ELEMENT'
f.write(itme)
itme="%s "%'GRIB_UNIT'
f.write(itme)
f.write('\n')
cont=gb2.RasterCount
for i in range(1,300):
    gb2var=gb2.GetRasterBand(i).GetMetadata()
    iband="%d   "%i
    f.write(iband)
    itme="%s "%gb2var['GRIB_COMMENT']
    f.write(itme)
    itme="%s "%gb2var['GRIB_ELEMENT']
    f.write(itme)
    itme="%s "%gb2var['GRIB_UNIT']
    f.write(itme)
    f.write('\n')
    itme="%s "%gb2var['GRIB_SHORT_NAME']
    f.write(itme)
#    print gb2var['GRIB_COMMENT']
#print gb2var
#print gb2var.GetScale()
#print gb2var.GetUnitType()
f.close()
"""