from astropy.io import fits
import numpy as np
import astropy.units as u
import math 
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import visualization
import os
import pdb
from scipy.optimize import curve_fit as cf
from astropy.modeling import powerlaws, fitting
from operator import add,sub,truediv
from matplotlib import gridspec
from astropy.table import QTable, vstack, hstack
import radio_beam
import sys
import regions
from astropy.wcs import WCS

#Let's make the new ntot/dust criterion abundance images first
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
for source in fielddict.keys():
    fnum=fielddict[source]
    print(f'Source: {source}')
    base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
    homedict={'SgrB2S':'/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/','DSi':'/nov2022continuumsanitycheck/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}

    home=base+homedict[source]
    ntotmap=home+'bootstrap_smoothed_ntot_to_bolocamfeathercont_3sigma.fits'
    nh2map=home+'bootstrap_nh2map_3sigma_bolocamfeather_smoothedtobolocam.fits'
    abunmap=home+'bootstrap_ch3ohabundance_error_ntotintercept_bolocamfeather_smoothedtobolocam.fits'

    ntot=fits.getdata(ntotmap)
    nh2=fits.getdata(nh2map)
    abunfits=fits.open(abunmap)
    abun=abunfits[0].data

    ok=np.isfinite(ntot)*np.isfinite(nh2)

    maskabun=np.ma.masked_where(ok==False,abun)
    filledabun=maskabun.filled(fill_value=np.nan)
    maskabun2=np.ma.masked_where(filledabun>=1e-6,filledabun)
    filledabun=maskabun2.filled(fill_value=np.nan)
    #plt.imshow(filledabun,origin='lower',vmax=9e-7)
    templateheader=abunfits[0].header
    boottexphdu=fits.PrimaryHDU([filledabun])
    boottexphdu.header=templateheader
    boottexphdu.header['BTYPE']='CH3OH Abundance'
    boottexphdu.header['BUNIT']='cm-2/cm-2'
    boottexhdul=fits.HDUList([boottexphdu])
    boottexhdul.writeto(home+'bootstrap_ch3ohabundance_error_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits',overwrite=True)#'error_trot_boostrap1000_nonegativeslope.fits',overwrite=True)
    print(f'Saved Ntot/NH2 mask abundance for {source} at {home}bootstrap_ch3ohabundance_error_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits')



