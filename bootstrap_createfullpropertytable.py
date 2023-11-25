from astropy.table import QTable, hstack
import astropy.units as u
import numpy as np
import math
import pdb
from astropy.io import fits
import sys

dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

source=['DSi','DSii','DSiii','DSiv','DSv','DSVI','DSVII','DSVIII','DSIX','SgrB2S']
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
homedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}

dataversion='sep2023revolution'
datadir=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'

pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}
sumtable=QTable.read(datadir+'hotcoresummarytable.fits')

radii=sumtable['Radius']

#order DS1,2,3,4,5,6,7,8,9,SgrB2S
ntot=[]
ntoterr=[]
abuns=[]#[2.7976e-7,1.30031e-7,1.91618e-8,3.80548e-8,1.97814e-8,1.4809e-7,4.26428e-8,5.21027e-8,7.227e-8,7.57974e-07]*u.dimensionless_unscaled
errabuns=[]#[1.17176e-8,5.47361e-9,6.53246e-10,1.80802e-9,1.50475e-9,6.03875e-9,1.97995e-9,2.54539e-9,3.88543e-9,3.10019e-8]*u.dimensionless_unscaled

for core in source:
    i=0
    fnum=fielddict[core]
    print(f'Source: {core}')
    base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{core}/'
    home=base+homedict[core]
    centralpix=pixdict[core]

    ntotmap=home+'bootstrap_smoothed_ntot_to_bolocamfeathercont_3sigma.fits'
    ntoterrmap=home+'bootstrap_smoothed_ntot_err.fits'
    abunmap=home+'bootstrap_ch3ohabundance_3sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
    abunerrmap=home+'bootstrap_ch3ohabundance_error_ntotintercept_bolocamfeather_smoothedtobolocam.fits'

    abundfits=fits.open(abunmap)
    abundances=abundfits[0].data
    err_abund=fits.getdata(abunerrmap)
    ntot=fits.getdata(abunmap)
    ntoterr=fits.getdata(abunerrmap)
    
    if core == 'SgrB2S':
        corerad=9000*u.AU#15000,radii[i]
    else:
        corerad=3000*u.AU
    cellsize=(np.abs(abundfits[0].header['CDELT1']*u.deg)).to('arcsec')
    #pixperbeam=(cntmbeam.sr/((cellsize**2).to('sr'))).value

    pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')
    r=math.ceil((corerad/pixtophysicalsize).to(''))

    yy2,xx2=np.indices(abundances.shape)
    rr2=((xx2-centralpix[1])**2+(yy2-centralpix[0])**2)**0.5
    mask2=rr2<r

    abunsinradius=abundances[rr2<r]
    abunerrinradius=err_abund[rr2<r]
    ntotinradius=ntot[rr2<r]
    ntoterrinradius=ntoterr[rr2<2]
    
    if core == 'SgrB2S' or core == 'DSVI':
        peakabundance=np.nanmax(abunsinradius)
        peakloc=np.where(abundances==peakabundance)
        erroronpeak=err_abund[int(peakloc[0]),int(peakloc[1])]
    else:
        peakabundance=abundances[centralpix[0],centralpix[1]]
        erroronpeak=err_abund[centralpix[0],centralpix[1]]

    peakntot=np.nanmax(ntotinradius)
    ntotpeakloc=np.where(ntot==peakntot)
    peakntoterr=ntoterr[int(ntotpeakloc[0]),int(ntotpeakloc[1])]
    print(peakntot)
    print(peakntoterr)
    pdb.set_trace()
    abuns.append(peakabundance)
    print(f'Abundance: {peakabundance}')
    errabuns.append(erroronpeak)
    print(f'Error: {erroronpeak}\n')
    #pdb.set_trace()
    i+=1

#sys.exit()#pdb.set_trace()
abuntable=QTable(data=[abuns],names=['X(CH3OH)'])
errabuntable=QTable(data=[errabuns],names=['X(CH3OH) error'])

powerlawtable=QTable.read('contsanitycheck_powerlawtable_bootmasked.fits')
densitytable=QTable.read('contsanitycheck_densityslopes_bootmasked.fits')

compositetable=hstack([sumtable,powerlawtable,abuntable,errabuntable,densitytable])

compositepath=datadir+'densabunpowersummarytable.fits'

print(f'Writing table to {compositepath}')
compositetable.write(compositepath,overwrite=True)

