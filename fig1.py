#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 20:13:48 2016

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
import string
import numpy as np
import datetime
import matplotlib.cm as cm
import calendar
from matplotlib.ticker import NullFormatter
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
def mosoon(avee,stde,avew,stdw,etp,wtp):
    esta=np.zeros(shape=(6),dtype=int)
    wsta=np.zeros(shape=(6),dtype=int)
    ewet1=avee+stde
    ewet2=avee+2*stde
    ewet3=avee+3*stde
    edry1=avee-stde
    edry2=avee-2*stde
    edry3=avee-3*stde
    wwet1=avee+stdw
    wwet2=avew+2*stdw
    wwet3=avew+3*stdw
    wdry1=avew-stdw
    wdry2=avew-2*stdw
    wdry3=avew-3*stdw
    if etp>ewet1:
        esta[3]=1
    if etp>ewet2:
        esta[4]=1    
    if etp>ewet3:
        esta[5]=1    
    if etp>edry2 and etp<edry1:
        esta[0]=1        
    if etp>edry3 and etp<edry2:
        esta[1]=1        
    if etp<edry3:
        esta[2]=1     
    if wtp>wwet1:
        wsta[3]=1
    if wtp>wwet2:
        wsta[4]=1    
    if wtp>wwet3:
        wsta[5]=1    
    if wtp>wdry2 and wtp<wdry1:
        wsta[0]=1        
    if wtp>wdry3 and wtp<wdry2:
        wsta[1]=1        
    if wtp<wdry3:
        wsta[2]=1
    #print esta,wsta
    return esta,wsta
def getim(idd,days,nday):
    nd=0
    for i in range(0,12):
        if idd>= nd and idd<=nd+days[i]:
            im=i
            return im
        else:
            nd=nd+days[i]
###############################################################################    
dirin='D:/Work/MyPaper/PhD02/Data/'
dirout='D:/Work/MyPaper/PhD02/Pics/'
fpath=dirin+'Fig1data.txt'
iskp=1
onedim=readAscii(fpath,iskp)
N=len(onedim)/3
etpx=[]
wtpx=[]
xstrs=[]
datastr=[]
etpmon=np.zeros(shape=(12),dtype=float)
wtpmon=np.zeros(shape=(12),dtype=float)
moncont=np.zeros(shape=(12),dtype=float)
ys=1961
nday1=0
nday2=365
nnn=0
im=0
for i in range(0,N):
    days=[31,28,31,30,31,30,31,31,30,31,30,31]
    tmpstr='%d'%onedim[i*3]
    if tmpstr[4:6] in ('05','06','07','08','09'):
        #print tmpstr[0:4]
        #print tmpstr[4:6]
        #print i
        #print tmpstr[4:6]
        xstrs.append(tmpstr[2:4])
        etpx.append(onedim[i*3+1])
        wtpx.append(onedim[i*3+2])
        datastr.append(tmpstr)
    if i>=nday2:
        ys=ys+1
        nday1=nday2
    if calendar.isleap(ys):
        days[1]=29
        nday=366
    else:
        days[1]=28
        nday=365
    nday2=nday1+nday
    if i>=nday1 and i<=nday2-1:
        im=getim(i-nday1,days,nday)
        #print im,i-nday1,nday1
        etpmon[im]=etpmon[im]+onedim[i*3+1]
        #print etpmon[im]+etp[i]
        wtpmon[im]=wtpmon[im]+onedim[i*3+2]
        moncont[im]=moncont[im]+1.0            
#    elif i>=nday2:
#        nday1=nday2
for im in range(0,12):
    etpmon[im]=etpmon[im]/moncont[im]
    wtpmon[im]=wtpmon[im]/moncont[im]
#print wtpmon
NN=len(wtpx)
etp=np.zeros(shape=(NN),dtype=float)
wtp=np.zeros(shape=(NN),dtype=float)
for i in range(0,NN):
    etp[i]=etpx[i]
    wtp[i]=wtpx[i]
avee=etp.mean()
avew=wtp.mean()
stde=etp.std()
stdw=wtp.std()
###############################################################################
fpath=dirin+'MJJAS_'+"Bothrain_in_ETP+WTP.txt"
f1=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Norain_in_ETP+WTP.txt"
f2=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Normal_in_ETP+WTP.txt"
f3=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Dry_Only_in_ETP.txt"
f4=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Dry_Only_in_WTP.txt"
f5=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Rain_Only_in_WTP.txt"
f6=open(fpath,'w')
fpath=dirin+'MJJAS_'+"Rain_Only_in_ETP.txt"
f7=open(fpath,'w')
for i in range(0,NN):
    esta,wsta=mosoon(avee,stde,avew,stdw,etp[i],wtp[i])
    #print esta,wsta
    es=esta.sum()
    ws=wsta.sum()
    if esta[3]> 0 and wsta[3]>0:
       f1.write(datastr[i]+' ')
       tmpout='%f '%etp[i]+'%d '%esta[3]+'%d '%esta[4]+'%d '%esta[5]
       f1.write(tmpout)
       tmpout='%f '%wtp[i]+'%d '%wsta[3]+'%d '%wsta[4]+'%d '%wsta[5]
       f1.write(tmpout)
       f1.write('\n')
    if esta[0]> 0 and wsta[0]>0:
       f2.write(datastr[i]+' ')
       tmpout='%f '%etp[i]+'%d '%esta[0]+'%d '%esta[1]+'%d '%esta[2]
       f2.write(tmpout)
       tmpout='%f '%wtp[i]+'%d '%wsta[0]+'%d '%wsta[1]+'%d '%wsta[2]
       f2.write(tmpout)
       f2.write('\n')
    if es==0 and ws==0:
       f3.write(datastr[i]+' ')
       tmpout='%f '%etp[i]+'%d '%esta[0]+'%d '%esta[1]+'%d '%esta[2]
       f3.write(tmpout)
       tmpout='%f '%wtp[i]+'%d '%wsta[0]+'%d '%wsta[1]+'%d '%wsta[2]
       f3.write(tmpout)
       f3.write('\n')
    if esta[0]> 0 and wsta[0]==0:
       f4.write(datastr[i]+' ')
       tmpout='%f '%etp[i]
       for ix in range(0,6):
           tmpout=tmpout+'%d '%esta[ix]
       f4.write(tmpout)
       f4.write('\n')
    if wsta[0]> 0 and esta[0]==0:
       f5.write(datastr[i]+' ')
       tmpout='%f '%wtp[i]
       for ix in range(0,6):
           tmpout=tmpout+'%d '%wsta[ix]
       f5.write(tmpout)
       f5.write('\n')
    if wsta[3]> 0 and esta[3]==0:
       f6.write(datastr[i]+' ')
       tmpout='%f '%wtp[i]
       for ix in range(0,6):
           tmpout=tmpout+'%d '%wsta[ix]
       f6.write(tmpout)
       f6.write('\n')
    if esta[3]> 0 and wsta[3]==0:
       f7.write(datastr[i]+' ')
       tmpout='%f '%etp[i]
       for ix in range(0,6):
           tmpout=tmpout+'%d '%esta[ix]
       f7.write(tmpout)
       f7.write('\n')
eline=avee+1*stde
wline=avew+1*stdw
nevs=0
for x in etp:
    if x>=eline:
        nevs=nevs+1
print nevs,'ETP'
nevs=0        
for x in wtp:
    if x>=wline:
        nevs=nevs+1
print nevs,'WTP'
conte=0
contw=0
print conte,contw
xx=range(0,NN)
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(15,6))
ax.plot(xx,etp,c='b',marker='o',ls='None',markeredgecolor='None',markersize=6,alpha=0.6)
ax.plot(xx,wtp,c='r',marker='*',ls='None',markeredgecolor='None',markersize=6,alpha=0.6)
x=[0,NN-1]
y=[eline,eline]
ax.plot(x,y,'b',lw=4.5)
x=[0,NN-1]
y=[wline,wline]
ax.plot(x,y,'r',lw=4.5)
rainstr='%0.2f'%eline
#ax.text(NN+1,eline,rainstr,fontsize=18)
rainstr='%0.2f'%wline
#ax.text(NN+1,wline,rainstr,fontsize=18)
plt.xlim([0, NN+1])
ax.set_xticks(range(0,NN+1,92*4))
xticklabels = [xstrs[nn] for nn in range(0,NN+1,92*4)] 
ax.set_xticklabels(xticklabels, size=16)
ymajorLocator   = MultipleLocator(3) 
ax.yaxis.set_major_locator(ymajorLocator)
#xmajorLocator   = MultipleLocator(4) 
#axx.xaxis.set_major_locator(xmajorLocator)
plt.ylabel('Precipitation '+r'($mm$ $d^{-1}$)',size=20)
plt.xlabel(r'Year $(1961-2012)$',size=20)
binwidth = 0.25
#ax.set_xlim((0, 10))
#ax.set_ylim((0, 10))
bins = np.arange(0, 10 + binwidth, binwidth)
left, width = 0.1, 0.70
bottom, height = 0.15, 0.6
left_h = left + width + 0.01
bottom_h =bottom+height
#rect_scatter = [left, bottom, width, height]
#rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.15, height]
#axHistx = plt.axes(rect_histx)
#axHisty = plt.axes(rect_histy)
nullfmt = NullFormatter()
#axHisty.yaxis.set_major_formatter(nullfmt)
#axHisty.hist(etp, bins=bins)
nbn=len(bins)
etpbins=np.zeros(shape=(nbn),dtype=float)
wtpbins=np.zeros(shape=(nbn),dtype=float)
tn=len(etp)
for i in range(0,tn):
    for j in range(0,nbn):
        if j<nbn-1:
            if etp[i]>=bins[j] and etp[i]<bins[j+1]:
                etpbins[j]=etpbins[j]+1./tn
                #print j
            if wtp[i]>=bins[j] and wtp[i]<bins[j+1]:
                wtpbins[j]=wtpbins[j]+1./tn
        elif j==nbn-1:
            if etp[i]>=bins[j-1]:
                etpbins[j]=etpbins[j]+1./tn
            if wtp[i]>bins[j-1]:
                wtpbins[j]=wtpbins[j]+1./tn            
a = plt.axes(rect_histy)
plt.barh(bins+0.125,etpbins*100,height=0.25,color='b',alpha=0.6,label='ETP')
plt.barh(bins+0.125,wtpbins*100,height=0.25,color='r',alpha=0.6,label='WTP')
plt.legend(frameon=False,loc=[0.2,1.1],fontsize='xx-large')
plt.yticks([])
ax.text(NN+2000,-1,'%',fontsize=18)
rect_histx=[left, bottom_h, width, 0.2]
aa = plt.axes(rect_histx)
months=np.arange(0, 12, 1)
plt.bar(months,etpmon,width=1,color='b',alpha=0.6)
plt.bar(months,wtpmon,width=1,color='r',alpha=0.6)
#plt.yticks([])
monstr=['Jan.','Feb.','March','April','May','June','July','Aug.',
        'Sept.','Oct.','Nov.','Dec.']
plt.xticks(months, monstr)
#axHisty.hist(etpbins, bins,color='b',orientation='horizontal')
fig.subplots_adjust(left=0.1,bottom=0.15,right=1-0.2,top=1-0.25)
plt.savefig(dirout+"fig1_rainfall_3.png",dpi=450)          
plt.show()
plt.close()
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()














