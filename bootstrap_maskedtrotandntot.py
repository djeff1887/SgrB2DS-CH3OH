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

mpl.interactive(True)

plt.close('all')

fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
sources=list(fielddict.keys())
homedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}#{'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'Kfield3originals_200K_trial1_noexclusions/','DSVIII':'Kfield3originals_175K_trial1_noexclusions/','DSIX':'Kfield7originals_150K_trial1_noexclusions/'}

for source in sources:
    fnum=fielddict[source]
    print(f'Source: {source}')
    base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}'
    home=base+homedict[source]

    texmap=home+"texmap_3sigma_allspw_withnans_weighted.fits"
    texerrmap=home+'test_error_trot_boostrap10000_nonegativeslope.fits'#'error_trot_boostrap1000_nonegativeslope.fits'
    ntotmap=home+'bootstrap_smoothed_ntot_to_bolocamfeathercont.fits'#changed to "bootstrap" version after continuum sanity check, remove "bootstrap" if this causes a problem # we use the smoothed version since it needs to match the H2 map resolution
    ntoterrmap=home+'bootstrap_smoothed_ntot_err.fits'#'error_ntot_intstd_boostrap1000_nonegativeslope.fits'

    texmap=fits.open(texmap)
    texmapdata=texmap[0].data*u.K
    texerrdata=np.squeeze(fits.getdata(texerrmap))*u.K
    fitsntot=fits.open(ntotmap)
    ntots=fitsntot[0].data*u.cm**-2
    ntoterr=np.squeeze(fits.getdata(ntoterrmap))*u.cm**-2

    maskedtex=np.ma.masked_where(texmapdata < 3*texerrdata,texmapdata)
    filledmaskedtex=maskedtex.filled(np.nan)
    maskedntot=np.ma.masked_where(ntots < 3*ntoterr,ntots)
    filledmaskedntot=maskedntot.filled(np.nan)

    boottrotprimaryhdu=fits.PrimaryHDU(filledmaskedtex.value)
    boottrotprimaryhdu.header=texmap[0].header
    boottrothdul=fits.HDUList([boottrotprimaryhdu])
    boottrotpath=home+f'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
    print(f'Saving bootmasked temperature map at {boottrotpath}\n')
    boottrothdul.writeto(boottrotpath,overwrite=True)

    bootntotprimaryhdu=fits.PrimaryHDU(filledmaskedntot.value)
    bootntotprimaryhdu.header=fitsntot[0].header
    bootntothdul=fits.HDUList([bootntotprimaryhdu])
    bootntotpath=home+f'bootstrap_smoothed_ntot_to_bolocamfeathercont_3sigma.fits'
    print(f'Saving bootmasked column density map at {bootntotpath}\n')
    bootntothdul.writeto(bootntotpath,overwrite=True)
    '''
    plt.figure()
    plt.imshow(filledmaskedntot.value,origin='lower')
    plt.show()
    '''
