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
import matplotlib as mpl
from scipy.optimize import curve_fit as cf
from math import log10, floor
from astropy.stats import bootstrap
import astropy.stats
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters

mpl.interactive(True)

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
def line(x,slope,intercept):
    return slope*x+intercept

def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))
    
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

Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=150*u.GHz,
                                                         fmax=300*u.GHz)

source=os.getenv('envsource')#'DSIX'#os.getenv('envsource')#SgrB2S
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]

base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}'

sourcedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}

sourcepath=sourcedict[source]

savefigbase=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}'
savefighome=savefigbase+sourcepath

home=base+sourcedict[source]

print(f'Source: {source}')
print(f'Getting data from directory {home}')

'''
if os.path.isdir(rotdiagpath):
    print(f'Rotational diagram folder {rotdiagpath} already exists.')
    pass
else:
    print(f'Making rotational diagram folder {rotdiagpath}')
    os.mkdir(rotdiagpath)
    print('Directory created.\n')
'''
#filepath='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/pixelwiserotationaldiagrams/'
infile=open(home+'ch3ohlinesdict.obj','rb')#open(origsourcepath+'ch3ohlinesdict.obj','rb')
spwdict=pickle.load(infile)

fulltexmap=fits.getdata(home+'texmap_3sigma_allspw_withnans_weighted.fits')
#trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':150*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':175*u.K,'DSIX':150*u.K}
#testT=trotdict[source]
#qrot_partfunc=qrot(obsTrot)#Q_rot_asym(testT).to('')

print('Setting up and executing model fit')
testyshape=np.shape(fulltexmap)[0]
testxshape=np.shape(fulltexmap)[1]

texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))
ntoterrormap=np.empty((testyshape,testxshape))
nugsmap=fits.getdata(home+'alltransitions_nuppers.fits')
nugserrmap=fits.getdata(home+'alltransitions_nupper_error.fits')
allmaster=np.loadtxt(home+'mastereuksqnsfreqsdegens.txt',dtype=str)
mastereuks=[]
masterdegens=[]
masterqns=[]
eukshape=np.shape(mastereuks)
for master in range(len(allmaster[:,0])):
    mastereuks.append(float(allmaster[master,0]))
    masterdegens.append(float(allmaster[master,3]))
    masterqns.append(allmaster[master,1])
#masterdegens=np.loadtxt(home+'masterdegens.txt')
    
testzshape=len(mastereuks)

for y in np.arange(testyshape):
    print(f'Begin row: {y}')
    for x in np.arange(testxshape):
        #y=px[0]
        #x=px[1]
        #print(f'Fitting pixel y={y} x={x}')
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        qnsfitted=[]
        excludedqns=[]
        excludednuppers=[]
        excludedeuks=[]
        '''Assemble the list of N/g, Eupper, and errors for fitting'''
        for z in range(testzshape):
            if nugsmap[z,y,x] <= 0 or np.isnan(nugsmap[z,y,x]):
                continue
            elif nugsmap[z,y,x]/nugserrmap[z,y,x] == 1:
                #print(f'Excluded line detected: {masterqns[z]}')
                #print('Appending to exclusion lists')
                excludedqns.append(masterqns[z])
                excludednuppers.append(nugsmap[z,y,x])
                excludedeuks.append(mastereuks[z])
            else:
                nupperstofit.append(nugsmap[z,y,x])
                eukstofit.append(mastereuks[z])
                nuperrors.append(nugserrmap[z,y,x])
                qnsfitted.append(masterqns[z])
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            continue
            #print(f'Pixel ({y}, {x}) has no N_upper value > 0 or != nan')
        else:
            log10nuerr=[]
            log10variances=[]
            errstofit=[]
            for nup,nup_error in zip(nupperstofit,nuperrors):#,nuperrors,eukstofit,qnsfitted):
                #pdb.set_trace()
                templog10=nup_error/nup
                temperrfit=1/templog10
                #if np.isnan(templog10) or np.isnan(temperrfit):
                #    pdb.set_trace()
                log10nuerr.append(templog10)
                log10variances.append(templog10**2)
                errstofit.append(temperrfit)
                
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit), weights=errstofit)
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            qrot_partfunc=qrot(obsTrot)*u.dimensionless_unscaled
            
            bslist=[]
            for nug,euk,weight,var in zip(np.log10(nupperstofit),eukstofit,errstofit,log10variances):
                bslist.append((nug,euk,weight,var))
        
            numboots=1000
            bootresult=bootstrap(np.array(bslist),numboots)

            bootlines=[]
            bootTrots=[]
            bootInts=[]
            bootNums=[]
            bootSlopes=[]
            for boot in bootresult:
                tempfit=fit(linemod,boot[:,1],boot[:,0],weights=boot[:,2])
                bootlines.append(tempfit)
            for line in bootlines:
                tempbootTrot=-np.log10(np.e)/(line.slope)
                tempbootslope=line.slope.value
                tempbootNtot=qrot_partfunc*10**(line.intercept)
                tempbootintNtot=line.intercept.value
                if line.slope >= 0:#tempbootTrot >= 1000 or tempbootTrot <= 0:
                    continue
                else:
                    bootTrots.append(tempbootTrot)
                    bootInts.append(tempbootNtot)
                    bootNums.append(tempbootintNtot)
                    bootSlopes.append(tempbootslope)

            bootTstd=astropy.stats.mad_std(bootTrots)
            #slope_bootstd=astropy.stats.mad_std(bootSlopes)
            int_bootNstd=astropy.stats.mad_std(bootNums)
            interimfactor=fit_lin.intercept.value*(1-int_bootNstd)
            bootNstd=(qrot_partfunc*10**(interimfactor)).value
        

            #pdb.set_trace()
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            #obsTrot=-np.log10(np.e)/(fit_lin.slope)
            #obsTrotcf=-np.log10(np.e)/popt[0]
            #print(f'cf Trot: {obsTrotcf}')
            #obsNtot=qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])#eukstofit[0] is '5(1)-4(2)E1vt=0', eupper 55.87102 K
            altNtot=qrot_partfunc*10**(fit_lin.intercept)
            obsNtot=altNtot
        
            A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
            C=np.diagflat(log10variances)
            atc_1a=np.dot(np.dot(A.T, np.linalg.inv(C)), A)
            if np.linalg.det(atc_1a) == 0:
                print(f'Singular C matrix detected in pixel {y,x}')
                m_unc=np.nan
                b_unc=np.nan
            else:
                covmat = np.linalg.inv(atc_1a)
                m_unc = covmat[0,0]**0.5
                b_unc = covmat[1,1]**0.5
        
            dobsTrot=np.abs(m_unc/fit_lin.slope)*obsTrot*u.K
            dobsNtot=np.sqrt((qrot_partfunc*10**(fit_lin.intercept)*np.log(10)*b_unc)**2)*u.cm**-2
            snr=obsNtot/bootNstd#dobsNtot
        
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=bootTstd
            ntoterrormap[y,x]=bootNstd

texmaphdu=fits.open(home+'texmap_3sigma_allspw_withnans_weighted.fits')
templateheader=texmaphdu[0].header

bootntotphdu=fits.PrimaryHDU([ntoterrormap])
bootntotphdu.header=templateheader
bootntotphdu.header['BTYPE']='Total column density error (bootstrap)'
bootntotphdu.header['BUNIT']='cm-2'
bootntothdul=fits.HDUList([bootntotphdu])
bootntothdul.writeto(home+'error_ntot_intstd_boostrap1000_nonegativeslope.fits',overwrite=True) 
print(f'Saved Ntot error for {source} at {home}error_ntot_intstd_boostrap1000_nonegativeslope.fits')

bootntotphdu2=fits.PrimaryHDU([ntotmap])
bootntotphdu2.header=templateheader
bootntotphdu2.header['BTYPE']='Total column density (bootstrap)'
bootntotphdu2.header['BUNIT']='cm-2'
bootntothdul2=fits.HDUList([bootntotphdu2])
bootntothdul2.writeto(home+'bootstrap_ntot_intstd_boostrap1000_nonegativeslope.fits',overwrite=True) 
print(f'Saved Ntot for {source} at {home}ntot_intstd_boostrap1000_nonegativeslope.fits')

#Removed for ntot fix
boottexphdu=fits.PrimaryHDU([texerrormap])
boottexphdu.header=templateheader
boottexphdu.header['BTYPE']='Rotational temperature error (bootstrap)'
boottexphdu.header['BUNIT']='K'
boottexhdul=fits.HDUList([boottexphdu])
boottexhdul.writeto(home+'test_error_trot_boostrap10000_nonegativeslope.fits',overwrite=True)#'error_trot_boostrap1000_nonegativeslope.fits',overwrite=True)
print(f'Saved Trot error for {source} at {home}test_error_trot_boostrap10000_nonegativeslope.fits')


            
'''
plt.figure()
plt.imshow(texerrormap,origin='lower',vmin=0,vmax=500)
plt.show()

plt.figure()
plt.imshow(ntoterrormap,origin='lower',vmin=0,vmax=1e17)
plt.show()
'''
