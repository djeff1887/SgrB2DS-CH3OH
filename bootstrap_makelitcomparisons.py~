import numpy as np
from astropy.table import QTable, Table, vstack
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib as mpl
import pdb

plt.close('all')
mpl.interactive(True)

def gieser_massprofile(nin,rout,rin,p):
    return (nin*(rout/rin)**p*((4/3)*np.pi*rout**3*mu)).to('solMass')

mu=2.8*u.Dalton

names=['WB 89789 SMM1','W51 e2e','IRAS 23O33 1', 'IRAS 23033 2', 'IRAS 23033 3', 'IRAS 23151 1','IRAS 23385 1','IRAS 23385 2','AFGL 2591 1','CepA HW2 1', 'CepA HW2 2','G084.95505 1','G094.6028 1','G100.38 1', 'G108.75 1', 'G108.75 2', 'G075.78 1','IRAS 21078 1', 'IRAS 21078 2', 'NGC 7538 IRS9 1', 'S87 IRS1 1', 'W3 H2O 3', 'W3 H2O 4', 'W3 IRS4 1']
wbsmm1trot=[245]*u.K
w51trot=[575,390,435]*u.K#e2e,north,e8e, Goddi+2020
coretrots=[114.9,167.2,160.8,129.2,239.5,226,159.9,238.2,170.5,169.0,183.8,68.2,111.2,82.9,176.8,135.5,121.9,200.5,112.0,176.6,166.4,173.0]*u.K
trotwb_err=[4]*u.K
trotcore_errs=[8.2,6.6,8.0,4.7,12.1,2.3,6.8,11.4,13.3,0.8,45.4,0.7,8.2,2.6,0.3,0.2,3.4,0.9,5.2,10.9,9.0,15.4]*u.K
xch3oh=[3.8e-7,2e-7]
err_xch3oh=[1.3e-7]
lbol=[8.4e3,2.3e4]*u.solLum

wbsmm1mass=[13]*u.solMass
w51mass=[225,290,310]*u.solMass
eta150_coremasses=[6.06,7.81,5.22,3.28,4.43,2.3,7.43,1.24,0.28,1.67,2.35,1.85,2.58,3.73,9.21,1.6,1.7,1.6,2.15,11.61,11.22,1.14]#*u.solMass
eta150_err_coremasses=[1.29,1.59,1.08,0.67,0.92,0.46,1.52,0.26,0.06,0.33,0.76,0.37,0.55,0.76,1.84,0.32,0.34,0.32,0.44,2.44,2.33,0.25]*u.solMass
eta150_snr=np.array(eta150_coremasses)/np.array(eta150_err_coremasses)

gieser2021kappa=0.9*u.cm**2/u.g
gieser2021eta=150
jeff2022kappa=0.00858941*u.cm**2/u.g#includes g/d=100, taken at 1.3 mm/230.60958308 GHz using Adam's kappa function
jeff2022eta=100

ratio_convert_to_eta100=((gieser2021kappa)/(jeff2022kappa*gieser2021eta)).to('')
innercoremasses=np.array(eta150_coremasses*ratio_convert_to_eta100)*u.solMass
err_innercoremasses=np.array(innercoremasses)/eta150_snr

innerradii=np.array([1837,1837,1837,1392,2299,2299,1394,289,289,2273,1584,1452,2044,2044,1642,612,612,1114,1005,770,770,780])*u.AU
nins=(innercoremasses/((4/3)*np.pi*innerradii**3*mu)).to('cm-3')

plindices=[np.inf,np.inf,0.34,0.48,0.96,0.27,0.14,0.21,0.25,0.31,0.56,1.38,1.48,0.41,0.56,0.23,0.29,0.43,0.59,0.38,0.45,0.59,0.45,0.54]
wbsmm1radius=[5363]*u.AU
w51radius=[2700,1200,2500]*u.AU
coreradii=[5720,4767,3813,2926,4345,4345,2926,931,776,6097,4434,3880,4767,5720,5055,1330,1330,2394,2439,1774,1774,1774]*u.AU

dindices=[2.22,1.79,1.39,2.25,2.23,2.4,2.12,2.25,1.99,0.83,0.77,1.83,2.14,1.9,2.07,1.67,1.31,2.25,2.03,1.76,1.66,1.94]
err_dindices=[0.11,0.21,0.09,0.011,0.07,0.07,0.08,0.03,0.11,0.06,0.27,0.23,0.09,0.04,0.06,0.07,0.05,0.06,0.08,0.11,0.09,0.09]

coremasses=[]
for n,outer,inner,dindex in zip(nins,coreradii,innerradii,dindices):
    coremass=gieser_massprofile(n,outer,inner,dindex)
    coremasses.append(coremass.value)
err_coremasses=coremasses/eta150_snr
#pdb.set_trace()

comptable=Table.read('t150_compositetransposedensitytable.fits')

dsmasses=np.array(list(comptable[4])[1:])
errormass=np.array(list(comptable[5])[1:])#list(comptable[5])[1:]
lums=list(comptable[6])[1:]
errorlum=np.array(list(comptable[7])[1:])#list(comptable[7])[1:]
temps=list(comptable[0])[1:]
errortemps=list(comptable[1])[1:]
abuns=np.array(list(comptable[16]))[1:]
abuns=list(map(float,abuns))
errorabun=np.array(list(comptable[17]))[1:]
errorabun=list(map(float,errorabun))
nh2s=list(comptable[2])[1:]
dsradii=list(comptable[8])[1:]
err_dsradii=np.array(list(comptable[9])[1:])/2
dsdensindex=list(comptable[18])[1:]
err_dsdens=list(comptable[19])[1:]

wbfmt='o'
w51fmt='x'
corefmt='^'
dsfmt='*'

plt.rcParams['figure.dpi']=150
fig=plt.figure()
plt.errorbar(wbsmm1mass.value,wbsmm1trot.value,yerr=trotwb_err.value,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51mass,w51trot,fmt=w51fmt,label='W51')
plt.errorbar(coremasses,coretrots.value,yerr=trotcore_errs.value,xerr=err_coremasses,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsmasses,temps,yerr=errortemps,xerr=errormass,fmt=dsfmt,label='DS Hot Cores')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('T$_{peak}$ (K)',fontsize=14)
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.legend()
plt.show()

plt.figure()
plt.errorbar(wbsmm1mass,wbsmm1radius,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51mass,w51radius,fmt=w51fmt,label='W51')
plt.errorbar(coremasses,coreradii,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsmasses,dsradii,xerr=errormass,fmt=dsfmt,label='DS Hot Cores')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('Radius (AU)',fontsize=14)
plt.legend()
plt.show()

plt.figure()
plt.errorbar(coremasses,dindices,yerr=err_dindices,fmt='^',color='green',label='CORE Catalogue')
plt.errorbar(dsmasses,dsdensindex,yerr=err_dsdens,fmt='*',color='red',label='DS Hot Cores')
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('$p$',fontsize=14)
plt.xscale('log')
plt.legend()
plt.show()

plt.figure()
plt.hist(dindices,label='CORE Catalogue')
plt.hist(dsdensindex,label='DS Hot Cores')
plt.xlabel('$p$',fontsize=14)
plt.ylabel('Counts',fontsize=14)
plt.legend()
plt.show()

plt.figure()
plt.scatter(coreradii,dindices,label='CORE Catalogue')
plt.scatter(dsradii,dsdensindex,label='DS Hot Cores')
plt.xlabel('radius')
plt.ylabel('$p$')
plt.show()

plt.figure()
plt.errorbar(wbsmm1radius.value,wbsmm1trot.value,yerr=trotwb_err.value,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51radius,w51trot,fmt=w51fmt,label='W51')
plt.errorbar(coreradii.value,coretrots.value,yerr=trotcore_errs.value,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsradii,temps,yerr=errortemps,fmt=dsfmt,label='DS Hot Cores')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Radius (AU)',fontsize=14)
plt.ylabel('T$_{peak}$ (K)',fontsize=14)
plt.legend()
plt.show()
