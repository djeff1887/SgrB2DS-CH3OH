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

mpl.interactive(True)

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
def line(x,slope,intercept):
    return slope*x+intercept
    
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
#sourceid='SgrB2S'
#fnum=1
sourceid='DSii'
fnum=10

sourcepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/'
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial5_newexclusions3andexclusionsnotinfit/"#pixelwiserotationaldiagrams/"
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial3_3kmsslab_varyingslabwidth/"
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial2_5kmsslab/"
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial1/"
#sourcpath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/K_OctReimage_restfreqfix_newvelmask_newpeakamp'
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/checkthatutilitiesworks_K_OctReimage_restfreqfix_newvelmask_newpeakamp/"
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/"
#sourceid='DSv'
#fnum=10
#sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/field10originals_z0_000190713_5-6mhzwidth_stdfixes/"
'''
sourceid='SgrB2S'
fnum=1
sourcepath="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSii_iii/z_0_00017594380066803095_box1_5-6mhzwidth/"#'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/correctedcubesandmaps/'
'''
home=sourcepath
rotdiagpath=home+'pixelwiserotationaldiagrams/'

if os.path.isdir(rotdiagpath):
    print(f'Rotational diagram folder {rotdiagpath} already exists.')
    pass
else:
    print(f'Making rotational diagram folder {rotdiagpath}')
    os.mkdir(rotdiagpath)
    print('Directory created.\n')
#filepath='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/pixelwiserotationaldiagrams/'
infile=open(sourcepath+'ch3ohlinesdict.obj','rb')
spwdict=pickle.load(infile)

fulltexmap=fits.getdata(home+'texmap_3sigma_allspw_withnans_weighted.fits')
trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K}
testT=trotdict[sourceid]
qrot_partfunc=Q_rot_asym(testT).to('')

print('Setting up and executing model fit')
testyshape=np.shape(fulltexmap)[0]
testxshape=np.shape(fulltexmap)[1]

texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))
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

ypix=int(input('y coord:'))
xpix=int(input('x coord:'))
pixel=(ypix,xpix)
pixellist=list([pixel])

for px in pixellist:
    y=px[0]
    x=px[1]
    rotdiagfilename=rotdiagpath+f'rotdiag_pixel_{y}-{x}_weighted2.png'
    if os.path.isfile(rotdiagfilename):
        print(f'Pixel ({y},{x}) already has rotational diagram.')
        print(f'See {rotdiagfilename} for png image\n\n')
        print(f'Fitting pixel y={y} x={x}')
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        qnsfitted=[]
        excludedqns=[]
        excludednuppers=[]
        excludedeuks=[]
        for z in range(testzshape):
            if nugsmap[z,y,x] <= 0 or np.isnan(nugsmap[z,y,x]):
                continue
            elif nugsmap[z,y,x]/nugserrmap[z,y,x] == 1:
                print(f'Excluded line detected: {masterqns[z]}')
                print('Appending to exclusion lists')
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
            print(f'Pixel ({y}, {x}) has no N_upper value > 0 or != nan')
        else:
            log10nuerr=[]
            log10variances=[]
            errstofit=[]
            for nup,nup_error in zip(nupperstofit,nuperrors):#,nuperrors,eukstofit,qnsfitted):
                #pdb.set_trace()
                templog10=nup_error/nup
                '''
                if templog10 == 1:
                    print(f'Excluded line detected, log10(N/S)={templog10}')
                    print(f'Excluded line: {qnsfitted[num]}')
                    print('Removing from list of to-fit values')
                    nupperstofit.remove(originalnupperstofit[num])
                    qnsfitted.remove(qnsfitted[num])
                    eukstofit.remove(eukstofit[num])
                    #nuperrors.remove(nuperrors[num])
                else:
                '''
                temperrfit=1/templog10
                #if np.isnan(templog10) or np.isnan(temperrfit):
                #    pdb.set_trace()
                log10nuerr.append(templog10)
                log10variances.append(templog10**2)
                errstofit.append(temperrfit)
                
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit), weights=errstofit)
            #popt,pcov=cf(line,eukstofit,np.log10(nupperstofit),sigma=errstofit)
            #perr = np.sqrt(np.diag(pcov))
            
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            #obsTrotcf=-np.log10(np.e)/popt[0]
            #print(f'cf Trot: {obsTrotcf}')
            obsNtot=qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])
            
            A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
            C=np.diagflat(log10variances)
            covmat = np.linalg.inv(np.dot(np.dot(A.T, np.linalg.inv(C)), A))
            m_unc = covmat[0,0]**0.5
            b_unc = covmat[1,1]**0.5
            #m_unccf=pcov[0,0]**0.5
            #b_unccf=pcov[1,1]**0.5
            
            dobsTrot=np.abs(m_unc/fit_lin.slope)*obsTrot*u.K#(eukstofit[0]*u.K)/(np.log(nupperstofit[0]/masterdegens[0])-np.log(obsNtot/qrot_partfunc))**2#*nuperrors[0]
            
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTrot.to('K').value
            
            plt.clf()
            print('Begin plotting')
            plt.errorbar(eukstofit,np.log10(nupperstofit),yerr=log10nuerr,fmt='o')
            plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {round(obsTrot, 4)} $\pm$ {round(dobsTrot.value, 2)*u.K}\nobsNtot: {round(obsNtot.value,3)/u.cm**2}'))
            plt.scatter(excludedeuks,np.log10(excludednuppers),marker='v',color='red')
            plt.title(f'field{fnum} {sourceid} pixel ({y},{x}) CH$_3$OH Rotational Diagram')
            plt.xlabel(r'E$_u$ (K)')
            plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
            plt.legend()
            plt.show()
            print('Done.')
        continue
    else:
        print(f'Fitting pixel y={y} x={x}')
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        qnsfitted=[]
        excludedqns=[]
        excludednuppers=[]
        excludedeuks=[]
        for z in range(testzshape):
            if nugsmap[z,y,x] <= 0 or np.isnan(nugsmap[z,y,x]):
                continue
            elif nugsmap[z,y,x]/nugserrmap[z,y,x] == 1:
                print(f'Excluded line detected: {masterqns[z]}')
                print('Appending to exclusion lists')
                excludedqns.append(masterqns[z])
                excludednuppers.append(nugsmap[z,y,x])
                excludedeuks.append(mastereuks[z])
                continue
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
            print(f'Pixel ({y}, {x}) has no N_upper value > 0 or != nan')
        else:
            log10nuerr=[]
            errstofit=[]
            log10variances=[]
            for num in range(len(nupperstofit)):
                templog10=(1/nupperstofit[num])*nuperrors[num]
                temperrfit=1/templog10
                log10nuerr.append(templog10)
                errstofit.append(temperrfit)
                log10variances.append(templog10**2)
                
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit), weights=errstofit)
            popt,pcov=cf(line,eukstofit,np.log10(nupperstofit),sigma=errstofit)
            perr = np.sqrt(np.diag(pcov))
            
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            obsTrotcf=-np.log10(np.e)/popt[0]
            print(f'cf Trot: {obsTrotcf}')
            obsNtot=qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])
            
            A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
            C=np.diagflat(log10variances)
            covmat = np.linalg.inv(np.dot(np.dot(A.T, np.linalg.inv(C)), A))
            m_unc = covmat[0,0]**0.5
            b_unc = covmat[1,1]**0.5
            
            dobsTrot=np.abs(m_unc/fit_lin.slope)*obsTrot*u.K#(eukstofit[0]*u.K)/(np.log(nupperstofit[0]/masterdegens[0])-np.log(obsNtot/qrot_partfunc))**2#*nuperrors[0]
            
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTrot.to('K').value
            
            plt.clf()
            print('Begin plotting')
            plt.errorbar(eukstofit,np.log10(nupperstofit),yerr=log10nuerr,fmt='o')
            plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {round(obsTrot, 4)} $\pm$ {round(dobsTrot.value, 2)*u.K}\nobsNtot: {round(obsNtot.value,3)/u.cm**2}'))
            plt.scatter(excludedeuks,np.log10(excludednuppers),marker='v',color='red')
            plt.title(f'field{fnum} {sourceid} pixel ({y},{x}) CH$_3$OH Rotational Diagram')
            plt.xlabel(r'E$_u$ (K)')
            plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
            plt.legend()
            plt.savefig((rotdiagfilename),dpi=100,overwrite=True)
            plt.show()
            print('Done.')
            print(f'Diagram saved to {rotdiagfilename}') 
           
'''
log10nuerr=[]
for num in range(len(n_us)):
    templog10=(1/n_us[num])*n_uerr[num]
    log10nuerr.append(templog10)
plt.clf()
print('Begin plotting')
plt.errorbar(mastereuks,np.log10(n_us),yerr=log10nuerr,fmt='o')
plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {obsTex*u.K}\nobsNtot: {obsNtot/u.cm**2}'))
plt.title(f'field1 {spwdict.keys()} CH$_3$OH Rotational Diagram')
plt.xlabel(r'E$_u$ (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
plt.legend()
plt.savefig((fieldpath+'rotdiag.png'),dpi=100,overwrite=True)
plt.show()
print('cubes loopered.')
'''