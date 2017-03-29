#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 08 16:54:53 2017

@author: chenjh
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
import datetime
import matplotlib.cm as cm
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
###############################################################################    
dirin='D:/MyPaper/PhD02/Pics/'
dirout='D:/MyPaper/PhD02/Pics/Patterns_MJJAS/'
fpath=dirin+'Fig1data.txt'
iskp=1
onedim=readAscii(fpath,iskp)
N=len(onedim)/3
etpx=[]
wtpx=[]
xstrs=[]
datastr=[]
for i in range(0,N):
    tmpstr='%d'%onedim[i*3]
    if tmpstr[4:6] in ('05','06','07','08','09'):
        #print tmpstr[0:4]
        #print tmpstr[4:6]
        xstrs.append(tmpstr[2:4])
        etpx.append(onedim[i*3+1])
        wtpx.append(onedim[i*3+2])
        datastr.append(tmpstr)
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
ax.text(NN+1,eline,rainstr,fontsize=18)
rainstr='%0.2f'%wline
ax.text(NN+1,wline,rainstr,fontsize=18)
plt.xlim([0, NN+1])
ax.set_xticks(range(0,NN+1,92*4))
xticklabels = [xstrs[nn] for nn in range(0,NN+1,92*4)] 
ax.set_xticklabels(xticklabels, size=16)
plt.ylabel('Precipitation '+r'($mm$ $d^{-1}$)',size=20)
plt.xlabel(r'Year $(1961-2012)$',size=20)
fig.subplots_adjust(left=0.08,bottom=0.15,right=1-0.05,top=1-0.05)
plt.savefig(dirout+"fig1_rainfall.png",dpi=300)          
plt.show()
plt.close()
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()
f7.close()













