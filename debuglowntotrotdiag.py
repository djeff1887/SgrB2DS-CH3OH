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
from utilities import *
import sys
from pyspeckit.spectrum.models import lte_molecule

mpl.interactive(True)
plt.rcParams["figure.dpi"]=150
plt.close('all')

def line(x,slope,intercept):
    return slope*x+intercept

def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

def testntot(n_upper, g, q, euj, tex):
    return n_upper*(q/g)*np.exp(euj/(k*tex))
    
def nupper_of_ntot(Ntot, eupper, tex, Q_rot, degeneracy=1):
    """ Given an N_tot, E_upper, tex, Q_rot, and degeneracy for a single state, give N_upper

    Mangum & Shirley 2015 eqn 31

    .. math::

        N_{tot}/N_u = Q_{rot} / g_u \\exp(E_u/k T_{ex})

    Example:
        >>> import astropy.units as u
        >>> tex = 50*u.K
        >>> kkms = 100*u.K*u.km/u.s
        >>> from pyspeckit.spectrum.models import lte_molecule
        >>> freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(molecule_name='HNCO v=0', molecule_name_jpl='HNCO', fmin=87*u.GHz, fmax=88*u.GHz)
        >>> nupper = lte_molecule.nupper_of_kkms(kkms, freqs, 10**aij)
        >>> ntot = lte_molecule.ntot_of_nupper(nupper, EU*u.erg, tex, Q_rot=partfunc(tex), degeneracy=deg)
    """

    const = (Q_rot/degeneracy) * np.exp(eupper / (k*tex))
    nupper = Ntot / const

    return nupper
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

source='DSi'
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]

origsourcedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/','DSiii':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/','DSiv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/','DSv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_trial1/','DSVI':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial2_16_6-16_7excluded/','DSVII':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_trial1_noexclusions/",'DSVIII':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVIII/Kfield3originals_175K_trial1_noexclusions/",'DSIX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/Kfield7originals_150K_trial1_noexclusions/'}
sourcedict={'SgrB2S':'/nov2022continuumsanitycheck/','DSi':'/aug2023fulloverhaul/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}#
dopplershifts={'SgrB2S':0.000228,'DSi':0.0001865,'DSii':0.000163,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001662062062062062,'DSIX':0.00015453732389175085,'DS10':0.00015794099431186572}
z=dopplershifts[source]
z_vel=z*c

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

ypix=36
xpix=42
pixel=(ypix,xpix)

fulltexmap=fits.getdata(home+'texmap_3sigma_allspw_withnans_weighted.fits')
trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':175*u.K,'DSIX':150*u.K}
testT=fulltexmap[pixel[0],pixel[1]]*u.K#trotdict[source]
modelqrot_partfunc=qrot(testT)*u.dimensionless_unscaled

nch3oh=5e18*u.cm**-2
freqmin=216.8*u.GHz
freqmax=233.4*u.GHz
                     
nuppers=[]
radnuppers=[]
error_radnuppers=[]
trads=[]
pynuppers=[]
adamnuppers=[]

testntots=[]

Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=freqmin,
                                                         fmax=freqmax)
                     
methanol_table=Splatalogue.query_lines(freqmin, freqmax, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=['JPL'], show_upper_degeneracy=True)
minmethtable=utils.minimize_table(methanol_table)
mlines=((minmethtable['Freq']*10**9)/(1+z)*u.Hz).to('GHz')
mqns=minmethtable['QNs']
meuks=minmethtable['EU_K']*u.K
meujs=meuks*k
'''
for euk in meuks:
    meujs.append(KtoJ(euk))
'''
mdegs=methanol_table['Upper State Degeneracy']
mlog10aijs=minmethtable['log10_Aij']
maijs=10**mlog10aijs*u.s**-1

restlines=mlines*(1+z)
brightness=80*u.K
linewidth=5*u.km/u.s
trad=brightness*linewidth
error_trad=trad/3

for line,deg,euj,aij,qn in zip(restlines,mdegs,meujs,maijs,mqns):
    est_nupper=nupper_estimated(nch3oh,deg,modelqrot_partfunc,euj,testT).to('cm-2')
    adam_nupper=nupper_of_ntot(nch3oh,euj,testT,modelqrot_partfunc,deg).to('cm-2')
    temp_ntot=testntot(est_nupper,deg,modelqrot_partfunc,euj,testT).to('cm-2')
    
    linewidthghz=velocitytofreq(linewidth,line)
    phi_nu=lineprofile(linewidthghz,line,line)
    intertau=lte_molecule.line_tau(testT, nch3oh, modelqrot_partfunc, deg, line, euj, aij)
    est_tau=(intertau*phi_nu).to('')
    trad=t_rad(f,est_tau,line,testT).to('K')
    velintK=trad*linewidth
    error_velintK=velintK/3
    test_nupper,test_nuppererror=N_u(line,aij,velintK,error_velintK)
    py_nupper=lte_molecule.nupper_of_kkms(velintK,line,aij)
    
    nuppers.append(est_nupper.value)
    radnuppers.append(test_nupper.value)
    error_radnuppers.append(test_nuppererror.value)
    trads.append(trad.value)
    pynuppers.append(py_nupper.value)
    testntots.append(temp_ntot.value)
    adamnuppers.append(adam_nupper.value)
    
plt.figure()
plt.scatter(meuks.value,testntots)
plt.yscale('log')
plt.xlabel(r'$E_U$ [K]')
plt.ylabel(r'$N_{tot}$ [cm$^{-2}$]')
plt.show()
#sys.exit()


plt.figure()
plt.title(fr'Debugging rotational diagrams: T={int(testT.value)} K; Ntot={nch3oh}')
plt.scatter(meuks.value,(nuppers/mdegs),label='Desmond\'s Ntot-Nupper function (Mangum and Shirley)')
plt.scatter(meuks.value,(radnuppers/mdegs),color='orange',label='Ntot-Brightness-Desmond\'s nupper_of_kkms')
plt.scatter(meuks.value,(pynuppers/mdegs),color='green',label='pyspeckit nupper_of_kkms')
plt.scatter(meuks.value,(adamnuppers/mdegs),color='purple',marker='*',label='Adam\'s Mangum and Shirley')
plt.yscale('log')
plt.xlim(xmax=1000)
plt.ylabel(r'$N_{upper}/g$ [cm$^{-2}$]')
plt.xlabel(r'$E_U$ [K]')
plt.legend()
plt.show()

plt.figure()
plt.scatter(meuks.value,trads)
plt.xlim(xmax=1000)
plt.ylim(ymin=0.1)
plt.show()

print('Setting up and executing model fit')

y=pixel[0]
x=pixel[1]
print(f'Fitting pixel y={y} x={x}')
linemod=models.Linear1D(slope=1.0,intercept=14)
fit=fitting.LinearLSQFitter()
eukstofit=[]
nuperrors=[]

log10nuerr=[]
log10variances=[]
errstofit=[]

nugs=nuppers/mdegs

for nug in nugs:
    #pdb.set_trace()
    sigma3all=(nug/3)*u.cm**-2
    nuperrors.append(sigma3all.value)
    
    templog10=sigma3all/(nug*u.cm**-2)
    
    temperrfit=1/templog10
    #if np.isnan(templog10) or np.isnan(temperrfit):
    #    pdb.set_trace()
    log10nuerr.append(templog10)
    log10variances.append(templog10**2)
    errstofit.append(temperrfit.value)

fit_lin=fit(linemod,meuks.value,np.log10(nugs), weights=errstofit)
obsTrot=-np.log10(np.e)/(fit_lin.slope)
print(f'slope: {fit_lin.slope}, obsTrot: {obsTrot}')

qrot_partfunc=qrot(obsTrot*u.K)*u.dimensionless_unscaled#9473.271845042498*u.dimensionless_unscaled

bslist=[]
for nug,euk,weight,var in zip(np.log10(nugs),meuks,errstofit,log10variances):
    bslist.append((nug,euk.value,weight,var))

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
linemod_euks=np.linspace(min(meuks.value),max(meuks.value),100)
#print('Model fit complete')
#print('Compute obsTex and obsNtot')
#obsTrot=-np.log10(np.e)/(fit_lin.slope)
#obsTrotcf=-np.log10(np.e)/popt[0]
#print(f'cf Trot: {obsTrotcf}')
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

A=np.stack((meuks.value,np.ones_like(meuks.value)),axis=1)
C=np.diagflat(log10variances)
covmat = np.linalg.inv(np.dot(np.dot(A.T, np.linalg.inv(C)), A))
m_unc = covmat[0,0]**0.5
b_unc = covmat[1,1]**0.5
#m_unccf=pcov[0,0]**0.5
#b_unccf=pcov[1,1]**0.5

dobsTrot=np.abs(m_unc/fit_lin.slope)*obsTrot*u.K#(eukstofit[0]*u.K)/(np.log(nupperstofit[0]/masterdegens[0])-np.log(obsNtot/qrot_partfunc))**2#*nuperrors[0]
dobsNtot=np.sqrt((qrot_partfunc*10**(fit_lin.intercept)*np.log(10)*b_unc)**2)*u.cm**-2#np.sqrt((qrot_partfunc*10**
print(f'dobsNtot: {dobsNtot.to("cm-2")}')
print(f'Ntot S/N: {obsNtot/dobsNtot.value}')
snr=obsNtot/bootNstd#dobsNtot

plt.figure()
#plt.clf()
print('Begin plotting')
tk='$T_{rot}$'
ntot='log$_{10}(N_{tot})$'
cm2='cm$^{-2}$'
#strdobsntot=str(dobsNtot.value)[0]
val_dntot=round((np.log10(bootNstd)/np.log10(obsNtot.value)),2)
val_ntot=round(np.log10(obsNtot.value),2)#round((obsNtot.value/(1*10**int(np.log10(obsNtot.value)))),(len(str(round(snr.value)))-1))

plt.errorbar(meuks.value,np.log10(nugs),yerr=log10nuerr,fmt='o')
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
#plt.savefig(f'debuggingrotationaldiagrams/JPLqrot_5e18cm-2_bootstrap_interr_{y}_{x}.png')
#rotdiagpath+f'contfix_bootstrap_interr_{y}_{x}.png')
plt.show()

plt.figure()
#plt.clf()
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
#plt.savefig(f'debuggingrotationaldiagrams/JPLqrot_5e18cm-2_trothist_bootstrap_interr_{y}_{x}.png')
#rotdiagpath+f'contfix_trothist_bootstrap_interr_{y}_{x}.png')
plt.show()
print('Done.')