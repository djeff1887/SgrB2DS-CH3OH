import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import glob
import math
import pdb
import matplotlib as mpl
from collections import OrderedDict
import sys

mpl.interactive(True)
plt.close('all')

def rolling_median(data,kernel):
    newlen=math.floor(len(data)/kernel)
    median=[]
    index=np.arange(newlen)
    temp=[]
    priorindex=0
    nantally=0
    for j in data:
        nantest=np.isnan(j)
        if nantest == True:
            if (priorindex+1) % kernel == 0 and priorindex != 0:
                med=np.nanmedian(temp)
                median.append(med)
                temp=[]
                nantally=0 
            else:
                nantally+=1
                pass
        elif nantally == 3:
            del radmed[priorindex+1]
            print(f'Deletion at index {priorindex+1}')
            nantally=0
            temp=[]
            temperr=[]
        else:
            dataindex=np.where(data==j)[0][0]
            if (dataindex+1) % kernel == 0 and dataindex != 0:
                temp.append(j)
                med=np.nanmedian(temp)
                median.append(med)
                temp=[]
                nantally=0
                #continueflag=True
            else:
                temp.append(j)
                #if continueflag:
                #   continue
        priorindex+=1#pdb.set_trace()
    if len(median) != newlen:
        pdb.set_trace()
    return median

def rolling_median2(data,kernel,error):
    newlen=math.floor(len(data)/kernel)
    median=[]
    mederr=[]
    index=np.arange(newlen)
    temp=[]
    temperr=[]
    priorindex=0
    nantally=0
    for j,e in zip(data,error):
        nantest=np.isnan(j)
        if nantest == True:
            if (priorindex+1) % kernel == 0 and priorindex != 0:
                med=np.nanmedian(temp)
                if np.any(temperr):
                    mee=np.nanmedian(temperr)
                    print('yay')
                else:
                    mee=np.nan
                    print('boo')
                median.append(med)
                mederr.append(mee)
                temp=[]
                temperr=[]
                nantally=0 
            else:
                nantally+=1
                pass
        elif nantally == 3:
            del radmed[priorindex+1]
            print(f'Deletion at index {priorindex+1}')
            nantally=0
            temp=[]
            temperr=[]
        else:
            dataindex=np.where(data==j)[0][0]
            if (dataindex+1) % kernel == 0 and dataindex != 0:
                temp.append(j)
                temperr.append(e)
                med=np.nanmedian(temp)
                mee=np.nanmedian(temperr)
                median.append(med)
                mederr.append(mee)
                temp=[]
                temperr=[]
                nantally=0
                #continueflag=True
            else:
                temp.append(j)
                temperr.append(e)
                #if continueflag:
                #   continue
        priorindex+=1#pdb.set_trace()
    if len(median) != newlen:
        pdb.set_trace()
    return median,mederr

pacmanfix=True
linestyles = {'solid':(0, ()),'densely dashdotdotted':(0, (3, 1, 1, 1, 1, 1)),'dotted':(0, (1, 5)),'densely dotted':(0, (1, 1)),'dashed':(0, (5, 5)),'densely dashed':(0, (5, 1)),'solid2':(0, ()),'dashdotted':(0, (3, 5, 1, 5)),'densely dashdotted':  (0, (3, 1, 1, 1)),'densely dashdotdotdot':(0, (3, 1, 1, 1, 1, 1))}

pacman='pacman_sep2023revolution'
dataversion='sep2023revolution'#homedict[source].replace('/','')
datadir=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'
pac_datadir='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/pacman_sep2023revolution/'

paths=glob.glob(datadir+'*radialavgabun*.txt')
errpaths=glob.glob(datadir+'*err_avgabun.txt')

if pacmanfix:
    print('Pacman Fix online')
    pacpath=glob.glob(pac_datadir+'*radialavgabun*.txt')
    errpacpath=glob.glob(pac_datadir+'*err_avgabun.txt')

    oldpath=paths[9]
    assert 'SgrB2S' in oldpath, 'SgrB2S not in old path'
    paths[9]=pacpath[0]
    comptable=Table.read(pac_datadir+'nov92023_allpropertytable.fits')
else:
    comptable=Table.read(datadir+'nov92023_allpropertytable.fits')

figpath='/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/'

sources=comptable['Source']
radii=comptable['Radius']
names={'DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9','SgrB2S':'SgrB2S'}
names2={'DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9','SgrB2S':'Sgr B2(S)'}

sortedpaths=[]
sortederrs=[]
for s in names.keys():
    tes=s+'_'
    for p in paths:
        if tes in p:
            sortedpaths.append(p)
for ss in names.keys():
    tess=ss+'_'
    for q in errpaths:
        if tess in q:
            sortederrs.append(q)

paths=sortedpaths
errpaths=sortederrs

plt.rcParams['figure.dpi']=300
plt.figure(figsize=(7,5))
for src,ls in zip(names.keys(),linestyles.keys()):
    test=src+'_'
    for i,e in zip(paths,errpaths):
        if test in i:
            sourceintable=np.where(sources==names[src])[0][0]
            profile=np.genfromtxt(i)
            proferr=np.genfromtxt(e)
            tempx=profile[0]
            listordered_centrtopix=tempx#[tempx<=radii[sourceintable]]
            radialabun=profile[1]
            radmed=rolling_median(listordered_centrtopix,3)
            abunmed,errorabunmed=rolling_median2(radialabun,3,proferr)
            upperabun=np.array(abunmed)+np.array(errorabunmed)
            lowerabun=np.array(abunmed)-np.array(errorabunmed)
            '''
            if len(np.where(lowerabun < 1e-9)[0]) != 0:
                pdb.set_trace()
            '''
            #print(f'{names2[src]}: {medtest[:5]}')
    plt.plot(radmed,abunmed[:len(radmed)],label=f'{names2[src]}',linestyle=linestyles[ls],)#yerr=errorabunmed[:len(radmed)])
    plt.fill_between(radmed,upperabun,lowerabun,alpha=0.2,)#color='blue')
plt.yscale('log')
plt.xlabel('$r$ (AU)',fontsize=14)
plt.ylabel('X(CH$_3$OH)',fontsize=14)
plt.ylim(ymin=5e-8)
plt.legend(bbox_to_anchor=(1.25, 1), loc="upper right")
#plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
figsavepath=figpath+f'radialavgabundiag_{pacman}_allsource_noradcut_fillbetween.pdf'
plt.savefig(figsavepath,bbox_inches='tight')
plt.show()
