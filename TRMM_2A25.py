#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 17:22:07 2015

@author: jhchen
"""
import h5py as hdf
from netCDF4 import Dataset
from pyhdf.SD import SD,SDC
from osgeo import gdal

fpath="E:/Data/TRMM/2A25/2006/2A25.20060101.46324.7.HDF"
toaobs=gdal.Open(fpath)
toaobs.GetSubDatasets()
