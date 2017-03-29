#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 09:03:46 2016

@author: chenjh
"""
import os
os.system("cls")

dirin='K:/DATA/CloudSat/201405-09/2014cldclass/'
dirout='K:/DATA/CloudSat/201405-09/'
for i in range(121,273):
    strs='%d'%i
    fpath=dirin+strs+'/'
    
    filenemes=os.listdir(fpath)
    l=len(filenemes)
    for a in range(0,l):
        newdir=dirout+filenemes[a]
        os.rename(fpath+filenemes[a],newdir)
        