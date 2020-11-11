import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import os
from astropy.modeling import models, fitting
import time
import pdb
import pickle

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
'''Collect constants and CH3OH-specific quantum parameters'''
print('Setting constants')
c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
mu_a=(0.896e-18*u.statC*u.cm).to('cm(3/2) g(1/2) s-1 cm')
R_i=1
f=1
Tbg=2.7355*u.K

pdb.set_trace()#Remember to update the below paths
home='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/correctedcubesandmaps/'
filepath='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/'
infile=open('NupperColDens/field1/testcore1/debug/allspwdict.obj','rb')
spwdict=pickle.load(infile)

fulltexmap=fits.getdata(home+'texmap_allspw_debug.fits')
testT=500*u.K
qrot_partfunc=Q_rot_asym(testT).to('')

print('Setting up and executing model fit')
testyshape=60
testxshape=60

texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))
nugsmap=fits.getdata(home+'alltransitions_nuppers.fits')
nugserrmap=fits.getdata(home+'alltransitions_nupper_error.fits')
mastereuks=np.loadtxt(home+'CH3OHmastereuks.txt')
testzshape=len(mastereuks)

for y in range(testyshape):
    print(f'Start Row {y} fitting')
    for x in range(testxshape):
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        print('Gathering nuppers and euks to fit')
        for z in range(testzshape):
            if nugsmap[z,y,x] <= 0:# or np.isnan(nugsmap[y,x,z]):
                continue
            else:
                nupperstofit.append(nugsmap[z,y,x])
                eukstofit.append(mastereuks[z])
                nuperrors.append(nugserrmap[z,y,x])
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
        else:
            log10nuerr=[]
            errstofit=[]
            for num in range(len(nupperstofit)):
                templog10=(1/nupperstofit[num])*nuperrors[num]
                temperrfit=1/templog10
                log10nuerr.append(templog10)
                errstofit.append(temperrfit)
                
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit),weights=errstofit)
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTex=-np.log10(np.e)/(fit_lin.slope)
            obsNtot=qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])
            dobsTex=(eukstofit[0]*u.K*np.log(10)*np.log(np.e))/(np.log(nupperstofit[0]/spwdict['spw2']['10_2--9_3-vt0']['degen'])-np.log(obsNtot/qrot_partfunc))**2
            
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTex.to('K').value
            
transmoment0=fits.open(spwdict['spw2']['10_2--9_3-vt0']['filename'])
transmom0header=transmoment0[0].header
            
primaryhdutex=fits.PrimaryHDU(texmap)
primaryhdutex.header=transmom0header
primaryhdutex.header['BTYPE']='Excitation temperature'
primaryhdutex.header['BUNIT']='K'
hdultex=fits.HDUList([primaryhdutex])
hdultex.writeto(filepath+'texmap_allspw_debug_weighted2.fits',overwrite=True)

primaryhduntot=fits.PrimaryHDU(ntotmap)
primaryhduntot.header=transmom0header
primaryhduntot.header['BTYPE']='Total column density'
primaryhduntot.header['BUNIT']='cm-2'
hdulntot=fits.HDUList([primaryhduntot])
hdulntot.writeto(filepath+'ntotmap_allspw_debug_weighted2.fits',overwrite=True)

primaryhdutexerr=fits.PrimaryHDU(texerrormap)
primaryhdutexerr.header=transmom0header
primaryhdutexerr.header['BTYPE']='Excitation temperature'
primaryhdutexerr.header['BUNIT']='K'
hdultexerror=fits.HDUList([primaryhdutexerr])
hdultexerror.writeto(filepath+'texmap_error_allspw_debug_weighted2.fits',overwrite=True)
            
plt.imshow(texmap,origin='lower')
plt.clim(10,500)
plt.show()