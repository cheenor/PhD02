#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Probability Distribution Function
Created on Sat Feb 07 10:28:28 2015

@author: jhchen
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import numpy as np
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20

dirin ="D:/MyPaper/PhD02/Data/CloudSat/"
pic_out="D:/MyPaper/PhD02/Pics/"

nz=105
rgn=['ETP','WTP']
fnm=['-TotalCloudPDF_baseon_CLDCLASS_AllPeriod.txt',
     '-CloudPDF_baseon_CLDCLASS_AllPeriod.txt']
cldname=['Cirrus','Altostratus','Altocumulus','Stratus',
    'Stratuscumulus','Cumulus','Nimbostratus','Deep Convection']
discld=[0,3,5,7]  
raintype=['NoRain','Drizzle','Liquid+Soild']
#####################################################################
onedim1=[]
linesplit=[]
path=dirin+"Cloudsat_HeightBin.txt"
f=open(path)
ff=f.readlines()[1:1+nz] 
for line in ff:
    lines=string.lstrip(line)
    linesplit.append(lines[:-1].split(' '))
for lnstrs in linesplit:
    for strs in lnstrs:
        if strs!='':
            onedim1.append(string.atof(strs))
height=np.ndarray(shape=(2,nz),dtype=float)
for iz in range(0,nz):
    k=iz*2
    height[0,iz]=onedim1[k]/1000.    ### convert to km
    height[1,iz]=onedim1[k+1]/1000.  ### convert to km
for ig in range(0,2):
    for ifn in range(0,2):
        if ifn>0 :
            rainmoncontnz=np.ndarray(shape=(5,8,3),dtype=float)
            rainmonnz=np.ndarray(shape=(5,8,3,nz),dtype=float)
            totalmonnz=np.ndarray(shape=(5,8,nz),dtype=float)
            totalmoncontnz=np.ndarray(shape=(5,8),dtype=float)
            rainmoncont=np.ndarray(shape=(5,3),dtype=float)
            rainmon=np.ndarray(shape=(5,8,3),dtype=float)
            totalmoncont=np.ndarray(shape=(5),dtype=float)
            totalmon=np.ndarray(shape=(5,8),dtype=float)
            for imm in range(5,10):
                im=imm-5
                monstr="%02d"%(imm)
                path=dirin +rgn[ig]+"-"+monstr+fnm[1]
                f=open(path)
                ff=f.readlines()[1:]
                linesplit=[]
                lines=string.lstrip(ff[0])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            tmp.append(string.atof(strs))
                totalmoncont[im]=tmp[0]
                linesplit=[]
                lines=string.lstrip(ff[1])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            tmp.append(string.atof(strs))               
                for i in range(0,8):
                    totalmon[im,i]=tmp[i]            
                for irain in range(0,3):
                    linesplit=[]
                    lines=string.lstrip(ff[irain*2+2])
                    linesplit.append(lines[:-1].split(' '))
                    tmp=[]
                    for lnstrs in linesplit:
                        for strs in lnstrs:
                            if strs!='':
                                tmp.append(string.atof(strs))
                    rainmoncont[im,irain]=tmp[0]
                    linesplit=[]
                    lines=string.lstrip(ff[irain*2+3])
                    linesplit.append(lines[:-1].split(' '))
                    tmp=[]
                    for lnstrs in linesplit:
                        for strs in lnstrs:
                            if strs!='':
                                tmp.append(string.atof(strs))
                    for i in range(0,8):
                        rainmon[im,i,irain]=tmp[i]
###############################################################################  
                linesplit=[]        
                lines=string.lstrip(ff[8])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            tmp.append(string.atof(strs))
                for i in range(0,8):
                        totalmoncontnz[im,i]=tmp[i]

                for iz in range(0,nz):
                    linesplit=[]
                    lines=string.lstrip(ff[iz+9])
                    linesplit.append(lines[:-1].split(' '))
                    tmp=[]
                    for lnstrs in linesplit:
                        for strs in lnstrs:
                            if strs!='':
                                tmp.append(string.atof(strs))
                    for i in range(0,8):
                        totalmonnz[im,i,iz]=tmp[i]                      
                for j in range(0,3):
                    sk=j*(nz+1)+114
                    linesplit=[]
                    lines=string.lstrip(ff[sk])
                    linesplit.append(lines[:-1].split(' '))
                    tmp=[]
                    for lnstrs in linesplit:
                        for strs in lnstrs:
                            if strs!='':
                                tmp.append(string.atof(strs))
                    for i in range(0,8):
                        rainmoncontnz[im,i,j]=tmp[i]
                    for iz in range(0,nz):
                        linesplit=[]
                        lines=string.lstrip(ff[iz+sk+1])
                        linesplit.append(lines[:-1].split(' '))
                        tmp=[]
                        for lnstrs in linesplit:
                            for strs in lnstrs:
                                if strs!='':
                                    tmp.append(string.atof(strs))
                        for i in range(0,8):
                            rainmonnz[im,i,j,iz]=tmp[i]
                f.close()
        else:
            rainallcontnz=np.ndarray(shape=(8,3),dtype=float)
            rainallnz=np.ndarray(shape=(8,3,nz),dtype=float)
            totalallnz=np.ndarray(shape=(8,nz),dtype=float)
            totalallcontnz=np.ndarray(shape=(8),dtype=float)
            rainallcont=np.ndarray(shape=(3),dtype=float)
            rainall=np.ndarray(shape=(8,3),dtype=float)
            totalallcont=np.ndarray(shape=(1),dtype=float)
            totalall=np.ndarray(shape=(8),dtype=float)
            path=dirin +rgn[ig]+fnm[0]
            f=open(path)
            ff=f.readlines()[1:]
            linesplit=[]
            lines=string.lstrip(ff[0])
            linesplit.append(lines[:-1].split(' '))
            tmp=[]
            for lnstrs in linesplit:
                for strs in lnstrs:
                    if strs!='':
                        tmp.append(string.atof(strs))
            totalallcont[0]=tmp[0]
            linesplit=[]
            lines=string.lstrip(ff[1])
            linesplit.append(lines[:-1].split(' '))
            tmp=[]
            for lnstrs in linesplit:
                for strs in lnstrs:
                    if strs!='':
                        tmp.append(string.atof(strs))               
            for i in range(0,8):
                totalall[i]=tmp[i]            
            for irain in range(0,3):
                linesplit=[]
                lines=string.lstrip(ff[irain*2+2])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            tmp.append(string.atof(strs))
                rainallcont[irain]=tmp[0]
                linesplit=[]
                lines=string.lstrip(ff[irain*2+3])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                                tmp.append(string.atof(strs))
                    for i in range(0,8):
                        rainall[i,irain]=tmp[i]
###############################################################################  
                linesplit=[]        
                lines=string.lstrip(ff[8])
                linesplit.append(lines[:-1].split(' '))
                tmp=[]
                for lnstrs in linesplit:
                    for strs in lnstrs:
                        if strs!='':
                            tmp.append(string.atof(strs))
                for i in range(0,8):
                        totalallcontnz[i]=tmp[i]

                for iz in range(0,nz):
                    linesplit=[]
                    lines=string.lstrip(ff[iz+9])
                    linesplit.append(lines[:-1].split(' '))
                    tmp=[]
                    for lnstrs in linesplit:
                        for strs in lnstrs:
                            if strs!='':
                                tmp.append(string.atof(strs))
                    for i in range(0,8):
                        totalallnz[i,iz]=tmp[i]                      
                    for j in range(0,3):
                        sk=j*(nz+1)+114
                        linesplit=[]
                        lines=string.lstrip(ff[sk])
                        linesplit.append(lines[:-1].split(' '))
                        tmp=[]
                        for lnstrs in linesplit:
                            for strs in lnstrs:
                                if strs!='':
                                    tmp.append(string.atof(strs))
                        for i in range(0,8):
                            rainallcontnz[i,j]=tmp[i]
                        for iz in range(0,nz):
                            linesplit=[]
                            lines=string.lstrip(ff[iz+sk+1])
                            linesplit.append(lines[:-1].split(' '))
                            tmp=[]
                            for lnstrs in linesplit:
                                for strs in lnstrs:
                                    if strs!='':
                                        tmp.append(string.atof(strs))
                            for i in range(0,8):
                                rainallnz[i,j,iz]=tmp[i]
            f.close()
##############################################################################
#####            Plotting 
#    total   
    discld=[0,3,5,7]  
    colors = ['b', 'grey', 'orange','c','g','r','lime', 'm']
    labels=['Cr','As','Ac','St','Sc','Cu','Ns','DC']
    width=[2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5]    
    titlsize1=20
    titlsize2=24    
    fig,axes = plt.subplots(nrows=2,ncols=2,figsize=(13,10))
    mondtrs=['May','Jun.','Jul.','Aug.','Sept.']
    titlelbs=[r'($a$) ALL the Conditions',r'($b$) NoRain',r'($c$) Drizzle',r'($d$) Liquid+Soild']
    indx = np.arange(5)   ### may to Sept.
    for i in range(0,4):
        k=i+1
        ax=plt.subplot(2,2,k)
        ax.set_ylim(0,100)
        ax.set_xlim(0,5)
        paxx=[]
        if i==0 :
            for ic in range(0,8):
                if ic ==0 :
                    paxx.append(plt.bar(indx, totalmon[:,ic],width=1.0,color=colors[ic]))
                else:
                    botm=np.ndarray(shape=(5),dtype=float)
#                    botm=totalmon[:,ic-1]
                    for ii in range(0,5):
                        botm[ii]=0.
                        for iii in range(0,ic):                            
                            botm[ii]=totalmon[ii,iii]+botm[ii]
                    paxx.append(plt.bar(indx, totalmon[:,ic],width=1.0,color=colors[ic],bottom=botm))
#                    del botm
            plt.ylabel(r'Frequency (%)',fontsize=titlsize1)
            plt.title(titlelbs[i],fontsize=titlsize1)
            plt.xticks(indx+0.5, mondtrs)
            plt.yticks(np.arange(0,101,20))            
        else:
            j=i-1
            for ic in range(0,8):
                if ic ==0 :
                    paxx.append(plt.bar(indx,rainmon[:,ic,j],width=1.0,color=colors[ic]))
                else:
                    botmx=np.ndarray(shape=(5),dtype=float)
                    for ii in range(0,5):
                        botmx[ii]=0.
                        for iii in range(0,ic):                            
                            botmx[ii]=rainmon[ii,iii,j]+botmx[ii]
                    paxx.append(plt.bar(indx, rainmon[:,ic,j],width=1.0,color=colors[ic],bottom=botmx))
#                    del botm
            if i==2:
                plt.ylabel(r'Frequency (%)',fontsize=titlsize1)
            plt.title(titlelbs[i],fontsize=titlsize1)
            plt.xticks(indx+0.5, mondtrs)
            plt.yticks(np.arange(0,101,20))
            if i==1:
                plt.legend(paxx, labels,bbox_to_anchor=(1.015, 1), loc=2, borderaxespad=0.)
#    fig.suptitle(raintype[j],size=titlsize2)
    fig.subplots_adjust(hspace=0.4)
    plt.show()
    plt.savefig(pic_out+rgn[ig]+"_Bar_CloudPDF.pdf")
    plt.close()     
    """
    mondtrs=['May','June','July','August','September','All Period']
    fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(8,15))
    for i in range(0,6):
        k=i+1
        plt.subplot(3,2,k)
#    for ax in zip(axes):
        plt.axis('equal')
        if i<5 :
            tmpdis=[]
            tmplabel=[]
            tmpcolor=[]
            for ic in range(0,8):
                if totalmon[i,ic]!=0. :
                    tmpdis.append(totalmon[i,ic])
                    tmplabel.append(labels[ic])
                    tmpcolor.append(colors[ic])
            plt.pie(tmpdis,labels=tmplabel, colors=tmpcolor,autopct='%1.1f%%')
            plt.title(mondtrs[i],size=titlsize1)
        else:
            tmpdis=[]
            tmplabel=[]
            tmpcolor=[]
            for ic in range(0,8):
                if totalall[ic]!=0. :
                    tmpdis.append(totalall[ic])
                    tmplabel.append(labels[ic])
                    tmpcolor.append(colors[ic])
            plt.pie(tmpdis,labels=tmplabel, colors=tmpcolor,autopct='%1.1f%%')
#            plt.pie(totalall[:],labels=labels, colors=colors,autopct='%1.1f%%')
            plt.title(mondtrs[i],size=titlsize1)
#        i+=1    
    fig.suptitle('ALL the Conditions',size=titlsize2)
    fig.subplots_adjust(hspace=0.2)
    plt.show()
    plt.savefig(pic_out+rgn[ig]+'_Allcondition_CloudPDF.pdf')
    plt.close() 
    for j in range(0,3):
        fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(8,15))
        i=0
        for i in range(0,6):
#        for ax in zip(axes):
            plt.subplot(3,2,i+1)
            plt.axis('equal')
            if i<5 :
                tmpdis=[]
                tmplabel=[]
                tmpcolor=[]
                for ic in range(0,8):
                    if rainmon[i,ic,j]!=0 :
                        tmpdis.append(rainmon[i,ic,j])
                        tmplabel.append(labels[ic])
                        tmpcolor.append(colors[ic])
                plt.pie(tmpdis,labels=tmplabel, colors=tmpcolor,autopct='%1.1f%%')
#                plt.pie(rainmon[i,:,j],labels=labels, colors=colors,autopct='%1.1f%%')
                plt.title(mondtrs[i],size=titlsize1)
            else:
                tmpdis=[]
                tmplabel=[]
                tmpcolor=[]
                for ic in range(0,8):
                    if rainall[ic,j]!=0 :
                        tmpdis.append(rainall[ic,j])
                        tmplabel.append(labels[ic])
                        tmpcolor.append(colors[ic])
                plt.pie(tmpdis,labels=tmplabel, colors=tmpcolor,autopct='%1.1f%%')
#                plt.pie(rainall[:,j],labels=labels, colors=colors,autopct='%1.1f%%')
                plt.title(mondtrs[i],size=titlsize1)
#            i+=1    
        fig.suptitle(raintype[j],size=titlsize2)
        fig.subplots_adjust(hspace=0.2)
        plt.show()
        plt.savefig(pic_out+rgn[ig]+"_"+raintype[j]+'_CloudPDF.pdf')
        plt.close() 
    """
    mondtrs=['May','June','July','August','September','All Period']
    colors = ['b', 'k', 'orange','c','g','r','lime', 'm']
    labels=['Cr','As','Ac','St','Sc','Cu','Ns','DC']
    discld2=[0,1,3,5,7]  
    width=[3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5]
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20,8))
    for inr in range(0,2):
        for inc in range(0,3):
#    for ax in zip(axes):
#        axes[i]=plt.subplots(2,3,i+1)
            axes[inr,inc].set_ylim(0,20)
            axes[inr,inc].set_xlim(0,20)
            i=inr*3+inc
            if i <5 :
                for ii in range(0,8):    #discld: 
                    axes[inr,inc].plot(totalmonnz[i,ii,:],height[ig,:],label=labels[ii],
                        c=colors[ii],lw=width[ii])  #(5,8,nz)
                    axes[inr,inc].set_title(mondtrs[i],fontsize=titlsize1)
                    ylabs='Height'+r' ($km$)'
                    axes[inr,inc].set_ylabel(ylabs,size=titlsize1)
            else:
                for ii in range(0,8):    #discld: 
                    axes[inr,inc].plot(totalallnz[ii,:],height[ig,:],label=labels[ii],
                        c=colors[ii],lw=width[ii])  #(5,8,nz)
                    axes[inr,inc].set_title(mondtrs[i],fontsize=titlsize1)
                    axes[inr,inc].legend()
                    ylabs='Height'+r' ($km$)'
                    axes[inr,inc].set_ylabel(ylabs,size=titlsize1)
    fig.suptitle(raintype[j],size=titlsize2)
    fig.subplots_adjust(hspace=0.4)
    plt.show()
    plt.savefig(pic_out+rgn[ig]+'_Allcondition_CloudPDF_profiles.pdf')
    plt.close()
    for j in range(0,3):
        fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20,8))
        i=0
#        for i in range(0,6):
#        for ax in zip(axes):
        for inr in range(0,2):
            for inc in range(0,3):
#            axes[i]=plt.subplots(2,3,i+1)
                axes[inr,inc].set_ylim(0,20)
                axes[inr,inc].set_xlim(0,20)
                i=inr*3+inc
                if i <5 :
                    for ii in range(0,8):    #discld:  
                        axes[inr,inc].plot(rainmonnz[i,ii,j,:],height[ig,:],label=labels[ii],
                            c=colors[ii],lw=width[ii])  #(5,8,nz)
                        axes[inr,inc].set_title(mondtrs[i],fontsize=titlsize1)
                        ylabs='Height'+r' ($km$)'
                        axes[inr,inc].set_ylabel(ylabs,size=titlsize1)
                else:
                    for ii in range(0,8):    #discld: 
                        axes[inr,inc].plot(rainallnz[ii,j,:],height[ig,:],label=labels[ii],
                            c=colors[ii],lw=width[ii])  #(5,8,nz)
                        axes[inr,inc].set_title(mondtrs[i],fontsize=titlsize1)
                        axes[inr,inc].legend()
                        ylabs='Height'+r' ($km$)'
                        axes[inr,inc].set_ylabel(ylabs,fontsize=titlsize1)
        fig.suptitle(raintype[j],size=titlsize2)
        fig.subplots_adjust(hspace=0.4)
        plt.show()
        plt.savefig(pic_out+rgn[ig]+'_'+raintype[j]+'_CloudPDF_profiles.pdf')
        plt.close()        