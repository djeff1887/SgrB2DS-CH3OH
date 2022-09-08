import numpy as np
import astropy.units as u
import astropy.constants as cnst
import pdb
import math
from astropy.table import Table
from astropy.io import fits

fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
homedict={'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'Kfield3originals_200K_trial1_noexclusions/','DSVIII':'Kfield3originals_175K_trial1_noexclusions/','DSIX':'Kfield7originals_150K_trial1_noexclusions/'}
sourcenamesfortable={'SgrB2S':'SgrB2S','DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9'}

h2upperdiffs={}
h2lowerdiffs={}
lumupperdiffs={}
lumlowerdiffs={}

for source in fielddict.keys():
#source='DSi'
    fnum=fielddict[source]
    base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
    home=base+homedict[source]
    sourcefortable=sourcenamesfortable[source]

    dGC=8.34*u.kpc

    pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}
    peakpix=pixdict[source]

    h2massmap=home+'h2massmap_3sigma_bolocamfeather_smoothedtobolocam.fits'
    h2masserrmap=home+'h2massmap_error_bolocamfeather_smoothedtobolocam.fits'
    lummap=home+'boltzmannlum_3sigma_bolocamfeather_smoothedtobolocam.fits'
    lumerrmap=home+'boltzmannlum_error_bolocamfeather_smoothedtobolocam.fits'

    h2masshdu=fits.open(h2massmap)

    compositetable=Table.read('compositetable_transpose_preservesnrsforsgrb2s_wip.fits')
    tablekeys=compositetable.keys()
    radius=compositetable[sourcefortable][8]*u.AU
    cellsize=(np.abs(h2masshdu[0].header['CDELT1']*u.deg)).to('arcsec')
    pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')
    radiusuncertainty=compositetable[sourcefortable][9]*u.AU

    h2mass=h2masshdu[0].data*u.solMass
    h2masserr=fits.getdata(h2masserrmap)*u.solMass
    lums=fits.getdata(lummap)*u.solLum
    lumserr=fits.getdata(lumerrmap)*u.solLum

    rpix=math.ceil((radius/pixtophysicalsize).to(''))
    rpix_inner=math.ceil(((radius-radiusuncertainty)/pixtophysicalsize).to(''))
    rpix_outer=math.ceil(((radius+radiusuncertainty)/pixtophysicalsize).to(''))

    yy,xx=np.indices(np.shape(h2mass))
    rr=((xx-peakpix[1])**2+(yy-peakpix[0])**2)**0.5
    #mask=rr<rpix

    h2massinrad=h2mass[rr<rpix]
    luminrad=lums[rr<rpix]

    h2mass_lowerbound=h2mass[rr<rpix_inner]
    lum_lowerbound=lums[rr<rpix_inner]

    h2mass_upperbound=h2mass[rr<rpix_outer]
    lum_upperbound=lums[rr<rpix_outer]

    lowerh2=np.nansum(h2mass_lowerbound)
    lowerlum=np.nansum(lum_lowerbound)

    upperh2=np.nansum(h2mass_upperbound)
    upperlum=np.nansum(lum_upperbound)

    trueh2=np.nansum(h2massinrad)
    truelum=np.nansum(luminrad)

    h2upperdiff=upperh2-trueh2
    h2lowerdiff=trueh2-lowerh2

    lumupperdiff=upperlum-truelum
    lumlowerdiff=truelum-lowerlum

    h2upperdiffs.update({source:h2upperdiff})
    h2lowerdiffs.update({source:h2lowerdiff})
    lumupperdiffs.update({source:lumupperdiff})
    lumlowerdiffs.update({source:lumlowerdiff})


