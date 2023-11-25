import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import glob
import math
import pdb
import matplotlib as mpl
from collections import OrderedDict
from astropy.io import fits
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

fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
#for source in fielddict.keys():
source='DSi'
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':'/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/','DSi':'/nov2022continuumsanitycheck/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}
pixdict={'SgrB2S':(26,14),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}

home=base+homedict[source]

abunpath=home+'bootstrap_ch3ohabundance_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits'

data=np.squeeze(fits.getdata(abunpath))
'''
if source == 'SgrB2S':
    wcsobj=WCS(smooth_trotfits[0].header)

    regs = regions.Regions.read('/blue/adamginsburg/d.jeff/imaging_results/regfiles/roughsgrb2smassregion_ignoresHIIregion.reg')
    pixreg = regs[0].to_pixel(wcsobj)
    pixmask = pixreg.to_mask()
    
    abunds=pixmask.cutout(abunds,fill_value=np.nan)
'''
    
yy, xx = np.indices(data.shape)
center = pixdict[source]#np.unravel_index(np.argmax(data), data.shape)
rr = ((xx-center[1])**2 + (yy-center[0])**2)**0.5
rrs = np.unique(rr)

avgs, edges = np.histogram(rr, bins=rrs[:], weights=data)
norms, _ = np.histogram(rr, bins=rrs[:], weights=np.ones_like(data))

plt.plot((edges[:-1] + edges[1:])/2, avgs/norms)

#plt.plot(rr.flat, data.flat, ',')
#plt.ylim(ymin=5e-8)
plt.yscale('log')
plt.show()