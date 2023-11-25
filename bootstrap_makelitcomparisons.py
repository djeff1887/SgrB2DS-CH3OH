import numpy as np
from astropy.table import QTable, Table, vstack
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib as mpl
import pdb
import glob
import sys
from utilities import *
from astropy.io import fits

plt.close('all')
mpl.interactive(True)

def gieser_massprofile(nin,rout,rin,p):
    return (nin*(rout/rin)**p*((4/3)*np.pi*rout**3*mu)).to('solMass')

mu=2.8*u.Dalton

dataversion='sep2023revolution'#homedict[source].replace('/','')
datadir=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'

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
eta150_nins=((np.array(eta150_coremasses)*u.solMass)/((4/3)*np.pi*innerradii**3*mu)).to('cm-3')

plindices=[np.inf,np.inf,0.34,0.48,0.96,0.27,0.14,0.21,0.25,0.31,0.56,1.38,1.48,0.41,0.56,0.23,0.29,0.43,0.59,0.38,0.45,0.59,0.45,0.54]
wbsmm1radius=[5363]*u.AU
w51radius=[2700,1200,2500]*u.AU
coreradii=[5720,4767,3813,2926,4345,4345,2926,931,776,6097,4434,3880,4767,5720,5055,1330,1330,2394,2439,1774,1774,1774]*u.AU

dindices=[2.22,1.79,1.39,2.25,2.23,2.4,2.12,2.25,1.99,0.83,0.77,1.83,2.14,1.9,2.07,1.67,1.31,2.25,2.03,1.76,1.66,1.94]
err_dindices=[0.11,0.21,0.09,0.011,0.07,0.07,0.08,0.03,0.11,0.06,0.27,0.23,0.09,0.04,0.06,0.07,0.05,0.06,0.08,0.11,0.09,0.09]

tindices=[0.34,0.48,0.96,0.27,0.14,0.21,0.25,0.31,0.56,1.35,1.48,0.41,0.56,0.23,0.29,0.43,0.59,0.38,0.45,0.59,0.45,0.54]
err_tindices=[0.11,0.2,0.07,0.11,0.05,0.04,0.08,0.02,0.1,0.05,0.27,0.22,0.09,0.02,0.05,0.05,0.01,0.04,0.07,0.06,0.03,0.08]

xch3oh_sgrb2n1=np.array([4.8e-5,1.1e-4])#N1S,1W Table 3 https://www.aanda.org/articles/aa/pdf/2022/09/aa43383-22.pdf
errorxch3oh_sgrb2n1=np.array([2.3e-5,0.5e-4])
nch3oh_sgrb2n=np.array([4e19,4e18,2.5e17,9.5e17])#N2-5
nh2_sgrb2n=np.array([1.8e24,0.5e24,2.9e24,1e24])
xch3oh_sgrb2n=(nch3oh_sgrb2n/nh2_sgrb2n)
xch3oh_wbsmm1=[3.8e-7]# Table 4 https://iopscience.iop.org/article/10.3847/1538-4357/ac289b/pdf
errorxch3oh_wbsmm1=[1.3e-7]
xch3oh_w51=[7.08e-6,2.338e-6,1.526e-6,3.2e-6,1.588e-6]#e2,e8,north-left,north-middle,north-right,mm14; data from Adam, published work is here https://iopscience.iop.org/article/10.3847/1538-4357/aa6bfa/pdf
xch3oh_orionkl=[2.2e-6]#Table 8 https://iopscience.iop.org/article/10.1088/0004-637X/787/2/112/pdf, comparable spatial resolution @ ~4000 AU res
trotforabun_w51=[318.276,212.7,189.821,228.801,216.67]
trotforabun_orionkl=[128]#Section 5.31 https://iopscience.iop.org/article/10.1088/0004-637X/787/2/112/pdf
sys.exit()
trot_sgrb2n1=[104,131]
errortrot_sgrb2n1=[1,2]
trot_sgrb2n=[160,170,190,180]
#pdb.set_trace()

coremasses=[]
eta150_fullmasses=[]
for n,e_n,outer,inner,dindex in zip(nins,eta150_nins,coreradii,innerradii,dindices):
    coremass=gieser_massprofile(n,outer,inner,dindex)
    eta150mass=gieser_massprofile(e_n,outer,inner,dindex)
    coremasses.append(coremass.value)
    eta150_fullmasses.append(eta150mass.value)
err_coremasses=coremasses/eta150_snr

percentdiff=np.abs((np.array(coremasses)-np.array(eta150_fullmasses))/(np.array(eta150_fullmasses)+np.array(coremasses)/2))#Doing percent difference with respect to the Gieser values, since we say "these values differ from the original values by XXX"
avg_ratio=np.mean(percentdiff)
print(f'Average percent difference between scaled and eta150 masses: {avg_ratio}')

comptable=Table.read(datadir+'pacman_lstarch3oh_densabunpowersummarytable.fits')
dsmasses=np.array(comptable['H_2 Mass'])#list(comptable[4])[1:])
errormass=np.array(comptable['H_2 Mass error'])#list(comptable[5])[1:])#list(comptable[5])[1:]
lums=np.array(comptable['Luminosity'])#list(comptable[6])[1:]
errorlum=np.array(comptable['Luminosity_error'])#list(comptable[7])[1:])#list(comptable[7])[1:]
temps=np.array(comptable['T_max'])#list(comptable[0])[1:]
errortemps=np.array(comptable['T_max_error'])#list(comptable[1])[1:]
abuns=np.array(comptable['X(CH3OH)'])#list(comptable[16]))[1:]
abuns=list(map(float,abuns))
errorabun=np.array(comptable['X(CH3OH) error'])#list(comptable[17]))[1:]
errorabun=list(map(float,errorabun))
nh2s=np.array(comptable['N(H_2) avg)'])#list(comptable[2])[1:]
dsradii=np.array(comptable['Radius'])#list(comptable[8])[1:]
err_dsradii=np.array(comptable['Radius_error'])#np.array(list(comptable[9])[1:])/2
dsdensindex=np.array(comptable['density_alpha'])#list(comptable[18])[1:]
err_dsdens=np.array(comptable['err_densityalpha'])#list(comptable[19])[1:]
dstempindex=np.array(comptable['alpha_2'])
dstempindex[3]=comptable['alpha_1'][3]
dstempindex[8]=comptable['alpha_1'][8]
dstempindex[9]=comptable['alpha_1'][9]
err_dstempindices=np.array(comptable['alpha_2 error'])
err_dstempindices[3]=comptable['alpha_1 error'][3]
err_dstempindices[8]=0.03#comptable['alpha_1 error'][8]

srealpeak=662.299
error_srealpeak=190.604


radialmassprofilepaths=glob.glob('*massinterior*')
names=['DSi_','DSii_','DSiii_','DSiv','DSv','DSVI_','DSVII_','DSVIII_','DSIX','SgrB2S']
orderedrmpps=[]
for src in names:
    for rmpp in radialmassprofilepaths:
        if src in rmpp:
            orderedrmpps.append(np.genfromtxt(rmpp))

src=0
edgeindices=[]#interiorrmpps=[]

for srcrmpp, rad in zip(orderedrmpps,dsradii):
    for i in srcrmpp[:,2]:
        if i > rad:
            index=np.where(srcrmpp[:,2]==rad)[0][0]
            edgeindices.append(index)#preint=np.delete(srcrmpp,np.arange(index,len(srcrmpp)))
            #interiorrmpps.append(preint)
            #print(f'edge:{index}')
            break
    src+=1
    #print(src)
    
#print(orderedrmpps)
#pdb.set_trace()

#pdb.set_trace()
meandsdensityindex=np.mean(dsdensindex)
meancoredensityindex=np.mean(dindices)
pdiff_dindex=np.abs(meandsdensityindex-meancoredensityindex)/np.mean([meandsdensityindex,meancoredensityindex])
err_meancoredindex=np.mean(err_dindices)
err_meandsdindex=np.mean(err_dsdens)

meandstempindex=np.mean(dstempindex)
meancoretempindex=np.mean(tindices)
pdiff_tindex=np.abs(meandstempindex-meancoretempindex)/np.mean([meandstempindex,meancoretempindex])
err_meandstempindex=np.mean(err_dstempindices)
err_meancoretempindex=np.mean(err_tindices)
print(f'Mean CORE density index: {meancoredensityindex} +/- {err_meancoredindex}')
print(f'Mean DS density index: {meandsdensityindex} +/- {err_meandsdindex}')
print(f'Difference between DS and CORE p indices: {pdiff_dindex*100}%\n')

print(f'Mean CORE temperature index: {meancoretempindex} +/- {err_meancoretempindex}')
print(f'Mean DS temperature index: {meandstempindex} +/- {err_meandstempindex}')
print(f'Difference between DS and CORE p indices: {pdiff_tindex*100}%')

trotsforabun_ds=np.array([])
errortrotsforabun_ds=np.array([])
abuntemp_s=0
errorabuntemp_s=0
for i, s in enumerate(sourcedict.keys()):
    if s == 'SgrB2S':
        fnum=fields[s]
        datafolder=sourcedict[s]
        base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{s}/{datafolder}'
        path_xch3oh=base+'bootstrap_ch3ohabundance_3sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
        path_trot=base+"bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
        path_errortrot=base+'test_error_trot_boostrap10000_nonegativeslope.fits'
        map_xch3oh=fits.getdata(path_xch3oh)
        map_trot=fits.getdata(path_trot)
        map_errortrot=np.squeeze(fits.getdata(path_errortrot))
        loc=np.where((map_xch3oh % abuns[9])==np.nanmin(map_xch3oh % abuns[9]))#Find the minimum mod to ID the peak abundance pixel 
        trot_peakabun=map_trot[loc[0][0],loc[1][0]]
        errortrot_peakabun=map_errortrot[loc[0][0],loc[1][0]]
        abuntemp_s=trot_peakabun
        errorabuntemp_s=errortrot_peakabun
    else:
        fnum=fields[s]
        datafolder=sourcedict[s]
        base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{s}/{datafolder}'
        path_xch3oh=base+'bootstrap_ch3ohabundance_3sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
        path_trot=base+"bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
        path_errortrot=base+'test_error_trot_boostrap10000_nonegativeslope.fits'
        map_xch3oh=fits.getdata(path_xch3oh)
        map_trot=fits.getdata(path_trot)
        map_errortrot=np.squeeze(fits.getdata(path_errortrot))
        loc=np.where((map_xch3oh % abuns[i-1])==np.nanmin(map_xch3oh % abuns[i-1]))
        trot_peakabun=map_trot[loc[0][0],loc[1][0]]
        errortrot_peakabun=map_errortrot[loc[0][0],loc[1][0]]
        trotsforabun_ds=np.append(trotsforabun_ds,trot_peakabun)
        errortrotsforabun_ds=np.append(errortrotsforabun_ds,errortrot_peakabun)
        
trotsforabun_ds=np.append(trotsforabun_ds,abuntemp_s)
errortrotsforabun_ds=np.append(errortrotsforabun_ds,errorabuntemp_s)

wbfmt='o'
w51fmt='x'
corefmt='^'
dsfmt='*'
n1fmt='p'
nfmt='s'
klfmt='P'
savefigpath='/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/LitComparisons/'

names=['Sgr B2(N)','Sgr B2(DS)','WB 89789 SMM1','W51','Orion KL']
colors={names[0]:'cyan',names[1]:'red',names[2]:'blue',names[3]:'orange',names[4]:'pink'}
numofregions=len(names)
regionxs=np.arange(numofregions)*2
all_xch3oh_sgrb2n=np.append(xch3oh_sgrb2n1,xch3oh_sgrb2n)
datasets=[all_xch3oh_sgrb2n,abuns,xch3oh_wbsmm1,xch3oh_w51,xch3oh_orionkl]
width = 0.25  # the width of the bars
multiplier = 0
'''
hist_range=[np.min(datasets[0]),np.max(datasets[0])]
binned_datasets = [np.histogram(d, range=hist_range, bins=numberofbins)[0] for d in datasets]
binned_maximums = np.max(binned_datasets, axis=1)
x_locations = np.arange(0, sum(binned_maximums), np.max(binned_maximums))

bin_edges = np.linspace(hist_range[0], hist_range[1], numberofbins + 1)
heights = np.diff(bin_edges)
centers = bin_edges[:-1]
'''

plt.rcParams['figure.dpi']=300

fig=plt.figure(figsize=(7,6))
plt.errorbar(trot_sgrb2n1,xch3oh_sgrb2n1,xerr=errortrot_sgrb2n1,yerr=errorxch3oh_sgrb2n1,marker=n1fmt,label='Sgr B2(N1)',color='cyan',ls='none')
plt.errorbar(trot_sgrb2n,xch3oh_sgrb2n,marker=nfmt,label='Sgr B2(N2-N5)',color='cyan',ls='none')
plt.errorbar(trotsforabun_ds,abuns,xerr=errortrotsforabun_ds,yerr=errorabun,fmt=dsfmt,label='Sgr B2(DS)',color='red')
plt.errorbar(wbsmm1trot.value,xch3oh_wbsmm1,xerr=trotwb_err.value,yerr=errorxch3oh_wbsmm1,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(trotforabun_w51,xch3oh_w51,marker=w51fmt,color='orange',label='W51',ls='none')
plt.errorbar(trotforabun_orionkl,xch3oh_orionkl,marker=klfmt,color='pink',label='Orion KL')
plt.yscale('log')
plt.xlabel('T$_{rot}$ (K)',fontsize=14)
plt.ylabel('X(CH$_3$OH)$_{peak}$',fontsize=14)
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_globalhotcore_xch3ohpeakvstrot.png')
plt.show()
#sys.exit()

widthofregionbars=0.5
fig, ax = plt.subplots(layout='constrained')
count=0

for dataset,regionx in zip(datasets,regionxs):
    assert len(datasets) == len(regionxs), 'Number of regions mismatched'
    firstbar=True
    numsources=len(dataset)
    sourceindex=np.arange(len(dataset))
    bindivision=numsources*width
    if numsources > 1:
        sourcexs=np.linspace((regionx-(widthofregionbars)),(regionx+(widthofregionbars)),numsources)
    else:
        sourcexs=[regionx]
    print(sourcexs)
    #pdb.set_trace()
    for obj,x in zip(dataset,sourcexs):
        offset = width * multiplier
        if firstbar:
            rects = ax.bar(x, obj, width,edgecolor='black',color=colors[names[count]],label=names[count])
            firstbar=False
        else:
            rects = ax.bar(x, obj, width,edgecolor='black',color=colors[names[count]])
        #ax.bar_label(rects, padding=3)
    count+=1
#ax.set_xticks(numofbins + width, names)
plt.yscale('log')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('X(CH$_3$OH)$_{peak}$')
ax.set_xticks(regionxs, names)
#ax.legend(loc='upper right')
plt.savefig(savefigpath+f'{dataversion}_globalhotcore_hist_xch3ohpeak.png')
plt.show()
sys.exit()
'''
fig, ax = plt.subplots()
for x_loc, binned_data in zip(x_locations, binned_datasets):
    lefts = x_loc - 0.5 * binned_data
    ax.barh(centers, binned_data, left=lefts,height=heights)

ax.set_xticks(x_locations, names)

ax.set_ylabel("XCH3OH")
#ax.set_xlabel("Data sets")

plt.show()
'''

fig=plt.figure(figsize=(7,6))
plt.errorbar(wbsmm1mass.value,wbsmm1trot.value,yerr=trotwb_err.value,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51mass,w51trot,fmt=w51fmt,label='W51')
plt.errorbar(coremasses,coretrots.value,yerr=trotcore_errs.value,xerr=err_coremasses,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsmasses[:8],temps[:8],yerr=errortemps[:8],xerr=errormass[:8],fmt=dsfmt,label='Sgr B2(DS1-DS9)')
plt.errorbar([dsmasses[9],dsmasses[9]],[temps[9],srealpeak],yerr=[errortemps[9],error_srealpeak],xerr=[errormass[9],errormass[9]],fmt=dsfmt,label='Sgr B2(S)')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('T$_{peak}$ (K)',fontsize=14)
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_tempvsmass.png')
plt.show()
sys.exit()

firstline=True
plt.figure(figsize=(7,6))
plt.errorbar(wbsmm1mass,wbsmm1radius,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51mass,w51radius,fmt=w51fmt,label='W51')
plt.errorbar(coremasses,coreradii,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsmasses,dsradii,xerr=errormass,fmt=dsfmt,label='DS Hot Cores')
for item,edge in zip(orderedrmpps,edgeindices):
    if firstline:
        plt.plot(item[:,0][:edge],item[:,2][:edge],ls='-',color='purple',linewidth=0.5,label='$r$-averaged Mass')
        firstline=False
    else:
        plt.plot(item[:,0][:edge],item[:,2][:edge],ls='-',color='purple',linewidth=0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('$R_{core}$ (AU)',fontsize=14)
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_radiusvsmass_lines.png')
plt.show()

plt.figure()
for item,edge in zip(orderedrmpps,edgeindices):
    mrsquared=(np.array(item[:,0])*u.solMass/(1*u.Da)).to('').value/((np.array(item[:,2])*u.AU).to('cm').value**2)
    #pdb.set_trace()
    plt.plot(item[:,2][:(edge+5)],mrsquared[:(edge+5)],ls='-',color='purple',linewidth=0.5,label='$r$-averaged Mass')
plt.yscale('log')
plt.xscale('log')
plt.show()
#sys.exit()

plt.figure()
plt.errorbar(coremasses,dindices,yerr=err_dindices,fmt='^',label='CORE Catalogue')#,color='green'
plt.errorbar(dsmasses,dsdensindex,yerr=err_dsdens,fmt='*',label='DS Hot Cores')#,color='red'
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('$p$',fontsize=14)
plt.xscale('log')
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_densityindexvsmass.png')
plt.show()

plt.figure()
plt.hist(dindices,label='CORE Catalogue')
plt.hist(dsdensindex,label='DS Hot Cores')
plt.xlabel('$p$',fontsize=14)
plt.ylabel('Counts',fontsize=14)
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_densityindexhistogram.png')
plt.show()

plt.figure()
coreradval=coreradii.value
plt.errorbar(coreradval,dindices,yerr=err_dindices,fmt=corefmt,label='CORE Catalogue')
plt.errorbar(dsradii,dsdensindex,yerr=err_dstempindices,fmt=dsfmt,label='DS Hot Cores')
plt.xlabel('$R_{core}$ (AU)')
plt.ylabel('$p$')
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_densityindexvsradius.png')
plt.show()

plt.figure()
plt.errorbar(wbsmm1radius.value,wbsmm1trot.value,yerr=trotwb_err.value,fmt=wbfmt,label='WB 89789 SMM1')
plt.errorbar(w51radius,w51trot,fmt=w51fmt,label='W51')
plt.errorbar(coreradii.value,coretrots.value,yerr=trotcore_errs.value,fmt=corefmt,label='CORE catalogue')
plt.errorbar(dsradii,temps,yerr=errortemps,fmt=dsfmt,label='DS Hot Cores')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('R$_{core}$ (AU)',fontsize=14)
plt.ylabel('T$_{peak}$ (K)',fontsize=14)
plt.savefig(savefigpath+f'{dataversion}_tempvsradius.png')
plt.legend()

plt.figure()
plt.errorbar(coretrots.value,tindices,yerr=err_tindices,xerr=trotcore_errs.value,fmt=corefmt,label='CORE catalogue')
plt.errorbar(temps,dstempindex,yerr=err_dstempindices,xerr=errortemps,fmt=dsfmt,label='DS Hot Cores')
plt.xlabel('T$_{peak}$ (K)',fontsize=14)
plt.ylabel('$q$')
plt.legend()
plt.savefig(savefigpath+f'{dataversion}_temperatureindexvstemp.png')
plt.show()

plt.figure()
plt.errorbar(temps,abuns,yerr=errorabun,xerr=errortemps,fmt=dsfmt,color='red')
plt.xlabel('T$_{peak}$ (K)',fontsize=14)
plt.ylabel('X(CH$_3$OH)$_{peak}$',fontsize=14)
#plt.xscale('log')
plt.yscale('log')
#plt.legend()
plt.savefig(savefigpath+f'{dataversion}_peakabundancevstemp.png')
plt.show()

