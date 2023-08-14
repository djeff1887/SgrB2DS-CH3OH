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
#from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters

mpl.interactive(True)
plt.rcParams["figure.dpi"]=150
plt.close('all')

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

'''
Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=150*u.GHz,
                                                         fmax=300*u.GHz)
'''

source='DSi'
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]

origsourcedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/','DSiii':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/','DSiv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/','DSv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_trial1/','DSVI':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial2_16_6-16_7excluded/','DSVII':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_trial1_noexclusions/",'DSVIII':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVIII/Kfield3originals_175K_trial1_noexclusions/",'DSIX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/Kfield7originals_150K_trial1_noexclusions/'}
sourcedict={'SgrB2S':'/nov2022continuumsanitycheck/','DSi':'/nov2022continuumsanitycheck/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}#{'SgrB2S':'/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/Kfield10originals_noexclusions/','DSiii':'/Kfield10originals_noexclusions/','DSiv':'/Kfield10originals_noexclusions/','DSv':f'/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':'/Kfield2originals_trial3_8_6-8_7excluded/','DSVII':'/Kfield3originals_200K_trial1_noexclusions/','DSVIII':'/Kfield3originals_175K_trial1_noexclusions/','DSIX':f'/Kfield7originals_150K_trial1_noexclusions/'}
sourcepath=sourcedict[source]
origsourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/{sourcepath}'#origsourcedict[source]

savefigbase=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}'
savefighome=savefigbase+sourcepath

home=origsourcepath
rotdiagpath=savefighome+'pixelwiserotationaldiagrams/'

if os.path.isdir(rotdiagpath):
    print(f'Rotational diagram folder {rotdiagpath} already exists.')
    pass
else:
    print(f'Making rotational diagram folder {rotdiagpath}')
    os.mkdir(rotdiagpath)
    print('Directory created.\n')
#filepath='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/NupperColDens/field1/testcore1/debug/pixelwiserotationaldiagrams/'
infile=open(origsourcepath+'ch3ohlinesdict.obj','rb')
spwdict=pickle.load(infile)

fulltexmap=fits.getdata(home+'texmap_3sigma_allspw_withnans_weighted.fits')
trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':175*u.K,'DSIX':150*u.K}
testT=trotdict[source]
#qrot_partfunc=qrot(testT).to('')

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

ypix=36#37#69#36dsi#73#70#int(input('y coord:'))61 dsiring:36
xpix=42#41#58#40dsi#54#55#int(input('x coord:'))64 dsiring:35
pixel=(ypix,xpix)
pixellist=list([pixel])

for px in pixellist:
    y=px[0]
    x=px[1]
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
        obsTrot=-np.log10(np.e)/(fit_lin.slope)
        qrot_partfunc=9473.271845042498*u.dimensionless_unscaled#qrot(obsTrot*u.K).to('')
        
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
        slope_bootstd=astropy.stats.mad_std(bootSlopes)
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
        bogusobsNtot=qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])#eukstofit[0] is '5(1)-4(2)E1vt=0', eupper 55.87102 K
        altNtot=qrot_partfunc*10**(fit_lin.intercept)
        print(f'Alt Ntot: {altNtot}')
        obsNtot=altNtot
        
        upper=fit_lin.intercept.value+int_bootNstd
        lower=fit_lin.intercept.value-int_bootNstd
        upperNtot=qrot_partfunc*10**(upper)
        lowerNtot=qrot_partfunc*10**(lower)
        upperdiff=upperNtot-obsNtot
        lowerdiff=obsNtot-lowerNtot
        upperfrac=upperdiff/obsNtot
        lowerfrac=lowerdiff/obsNtot
        #pdb.set_trace()

        A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
        C=np.diagflat(log10variances)
        covmat = np.linalg.inv(np.dot(np.dot(A.T, np.linalg.inv(C)), A))
        m_unc = covmat[0,0]**0.5
        b_unc = covmat[1,1]**0.5
        #m_unccf=pcov[0,0]**0.5
        #b_unccf=pcov[1,1]**0.5
        
        dobsTrot=np.abs(m_unc/fit_lin.slope)*obsTrot*u.K#(eukstofit[0]*u.K)/(np.log(nupperstofit[0]/masterdegens[0])-np.log(obsNtot/qrot_partfunc))**2#*nuperrors[0]
        dobsNtot=np.sqrt((qrot_partfunc*10**(fit_lin.intercept)*np.log(10)*b_unc)**2)*u.cm**-2#np.sqrt((qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])*np.log(10)*eukstofit[0]*m_unc)**2+(qrot_partfunc*10**(np.log10(nugsmap[0,y,x])+fit_lin.slope*eukstofit[0])*(1/(nugsmap[0,y,x]*np.log(10)))*nuperrors[0])**2)*u.cm**-2
        print(f'dobsNtot: {dobsNtot.to("cm-2")}')
        print(f'Ntot S/N: {obsNtot/dobsNtot.value}')
        snr=obsNtot/bootNstd#dobsNtot
        
        texmap[y,x]=obsTrot
        ntotmap[y,x]=obsNtot
        texerrormap[y,x]=dobsTrot.to('K').value
        
        plt.figure(1)
        plt.clf()
        print('Begin plotting')
        tk='$T_{rot}$'
        ntot='log$_{10}(N_{tot})$'
        cm2='cm$^{-2}$'
        #strdobsntot=str(dobsNtot.value)[0]
        val_dntot=round((np.log10(bootNstd)/np.log10(obsNtot.value)),2)
        val_ntot=round(np.log10(obsNtot.value),2)#round((obsNtot.value/(1*10**int(np.log10(obsNtot.value)))),(len(str(round(snr.value)))-1))

        plt.errorbar(eukstofit,np.log10(nupperstofit),yerr=log10nuerr,fmt='o')
        plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'{tk}: {round(obsTrot, 2)} $\pm$ {round(bootTstd,2)} K\n{ntot}: {val_ntot} $\pm$ {round(int_bootNstd,2)}'))
        #plt.errorbar(eukstofit,np.log10(nupperstofit),yerr=log10nuerr,fmt='o')
        #plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {round(obsTrot, 4)} $\pm$ {round(dobsTrot.value, 2)*u.K}\nobsNtot: {round(obsNtot.value,3)/u.cm**2}'))
        for linmod in bootlines:
            plt.plot(linemod_euks,linmod(linemod_euks),color='black',alpha=0.03,zorder=0)
        #plt.scatter(excludedeuks,np.log10(excludednuppers),marker='v',color='red')
        #plt.title(f'field{fnum} {source} pixel ({y},{x}) CH$_3$OH Rotational Diagram')
        plt.xlabel(r'E$_u$ (K)')
        plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
        plt.legend()
        plt.savefig(f'debuggingrotationaldiagrams/JPLqrot_contfix_bootstrap_interr_{y}_{x}.png')
        #rotdiagpath+f'contfix_bootstrap_interr_{y}_{x}.png')
        plt.show()
        
        plt.figure(2)
        plt.clf()
        plt.hist(bootTrots,bins=np.linspace(0,np.max(bootTrots),100))
        plt.axvline(obsTrot,color='red',label=r'Observed $T_{rot}$')
        plt.axvline(np.median(bootTrots),color='cyan',label=r'Median $T_{rot}$')
        plt.axvline(np.mean(bootTrots),color='yellow',label=r'Mean $T_{rot}$')
        plt.axvline(np.percentile(bootTrots,16),color='orange',ls='--',label=r'1$\sigma$ interval')
        plt.axvline(np.percentile(bootTrots,84),color='orange',ls='--')
        plt.plot([np.mean(bootTrots)-bootTstd,np.mean(bootTrots)+bootTstd],[0,0])
        plt.xlabel(r'$T_{rot}$ (K)')
        plt.ylabel('Number of fits')
        plt.legend()
        plt.xlim(xmin=0,xmax=750)
        plt.savefig(f'debuggingrotationaldiagrams/JPLqrot_contfix_trothist_bootstrap_interr_{y}_{x}.png')
        #rotdiagpath+f'contfix_trothist_bootstrap_interr_{y}_{x}.png')
        plt.show()
        print('Done.')
    continue
