import numpy as np
from astropy.io import fits
import astropy.units as u
import scipy.constants as cnst
from astropy.wcs import WCS
import math
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import reproject
import radio_beam
import os
from astropy.convolution import convolve
import sys

mpl.interactive(True)
plt.close('all')

def kappa(nu, nu0=271.1*u.GHz, kappa0=0.0114*u.cm**2*u.g**-1, beta=1.75):
    """
    Compute the opacity $\kappa$ given a reference frequency (or wavelength)
    and a power law governing the opacity as a fuction of frequency:
    $$ \kappa = \kappa_0 \left(\\frac{\\nu}{\\nu_0}\\right)^{\\beta} $$
    The default kappa=0.0114 at 271.1 GHz comes from extrapolating the
    Ossenkopf & Henning 1994 opacities for the thin-ice-mantle, 10^6 year model
    anchored at 1.0 mm with an assumed beta of 1.75.
    Parameters
    ----------
    nu: astropy.Quantity [u.spectral() equivalent]
        The frequency at which to evaluate kappa
    nu0: astropy.Quantity [u.spectral() equivalent]
        The reference frequency at which $\kappa$ is defined
    kappa0: astropy.Quantity [cm^2/g]
        The dust opacity per gram of H2 along the line of sight.  Because of
        the H2 conversion, this factor implicitly includes a dust to gas ratio
        (usually assumed 100)
    beta: float
        The power-law index governing kappa as a function of nu
    """
    return (kappa0*(nu.to(u.GHz,u.spectral())/nu0.to(u.GHz,u.spectral()))**(beta)).to(u.cm**2/u.g)

def gasmass(snu,omega,nu,kappa,tex):
    return (snu*beamarea_sr*d**2*c**2)/(2*k*tex*nu**2*kappa)
    
def gassmasserr(snu,omega,nu,kappa,tex,texerr):
    return np.sqrt((((snu*beamarea_sr*d*c**2)/(k*tex*nu**2*kappa))*d_unc)**2+(((snu*beamarea_sr*d**2*c**2)/(2*k*tex**2*nu**2*kappa))*texerr)**2).to('solMass')

def luminosity(r,T):#Boltzmann luminosity
    return (4*np.pi*r**2)*sigmasb*T**4

def error_luminosity(r,T):
    return (16*np.pi*r**2)*sigmasb*T**3*rms_1mm_K

def soundspeed(m,T):
    return np.sqrt((k*T)/m)

def jeansmass(cs,rho):
    return (np.pi/6)*cs**3*G**(-3/2)*rho**(-1/2)

def jeanslength(cs,rho):
    return cs*G**(-1/2)*rho**(-1/2)

def tff(rho):#Free=fall time of gas
    return np.sqrt((3*np.pi)/(32*G*rho))
    
def opticaldepth(snu,nu,tex):#has sr-1 units
    return ((snu*c**2)/(2*k*tex*nu**2))
    
c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
d=8.34*u.kpc#Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
d_unc=0.16*u.kpc#Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
sigmasb=cnst.sigma*u.W/(u.m**2*u.K**4)
mu=2.8*u.Dalton
G=cnst.G*u.m**3/(u.kg*u.s**2)
molweight_ch3oh=(32.042*u.u).to('g')#https://pubchem.ncbi.nlm.nih.gov/compound/methanol
Tcmb=2.722*u.K#Plank Collaboration (2015) https://ui.adsabs.harvard.edu/abs/2016A%26A...594A..13P/abstract

cntminfile='/orange/adamginsburg/sgrb2/2017.1.00114.S/imaging_results/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter_feathered_with_bolocam.fits'
notreproj_cntmimage=fits.open(cntminfile)
print(f'Continuum image: {cntminfile}')
restfreq=notreproj_cntmimage[0].header['RESTFRQ']*u.Hz

source='DSi'#os.getenv('envsource')
assert source is not None; 'os.getenv didn\'t work'
print(f'Source: {source}')

fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]

sourcedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}
sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}'+sourcedict[source]

texmap=fits.open(sourcepath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits')
texerrfits=fits.open(sourcepath+'test_error_trot_boostrap10000_nonegativeslope.fits')

ntotfits=fits.open(sourcepath+'bootstrap_ntot_intstd_boostrap1000_nonegativeslope.fits')
ntotmap=np.squeeze(ntotfits[0].data)*u.cm**-2
ntoterrfits=fits.open(sourcepath+'error_ntot_intstd_boostrap1000_nonegativeslope.fits')
ntoterrmap=np.squeeze(ntoterrfits[0].data)*u.cm**-2
ch3ohSigmam_map=ntotmap*molweight_ch3oh

#snrfits=fits.open(sourcepath+'texmap_snr_allspw_weighted.fits')
snrmap=texmap[0].data/texerrfits[0].data#snrfits[0].data


texerrmap=np.squeeze(texerrfits[0].data)
#maskedtexerr=np.ma.masked_where(texerrmap < 1,texerrmap)
maskedtexerr=np.ma.masked_where(texerrmap >= 500,texerrmap)
temptexerrmap=maskedtexerr.filled(fill_value=np.nan)*u.K
cmbmasktexerr=np.ma.masked_where(temptexerrmap.value < Tcmb.value, temptexerrmap.value)
texerrmap=cmbmasktexerr.filled(fill_value=np.nan)*u.K

maskedntoterr=np.ma.masked_where(ntoterrmap.value >= 1e19, ntoterrmap)
tempmaskntoterr=maskedntoterr.filled(fill_value=np.nan)
maskedntoterr2=np.ma.masked_where(tempmaskntoterr.value <= 1e12, tempmaskntoterr.value)
ntoterrmap=maskedntoterr2.filled(fill_value=np.nan)*u.cm**-2

reprojtexmappath=sourcepath+'bootstrap_texmap_3sigma_allspw_withnans_weighted_reprojecttobolocamcont_smoothed.fits'
reprojtexerrpath=sourcepath+'boostrap_texerrmap_allspw_withnans_weighted_useintercept_3sigma_reprojecttobolocamcont_smoothed.fits'
reprojsnrpath=sourcepath+'bootstrap_texmap_snr_allspw_weighted_reprojecttobolocamcont_smoothed.fits'
reprojntotpath=sourcepath+'bootstrap_ntotmap_allspw_withnans_weighted_useintercept_3sigma_reprojecttobolocamcont_smoothed.fits'
reprojntoterrpath=sourcepath+'bootstrap_ntoterrmap_allspw_withnans_weighted_useintercept_reprojecttobolocamcont_smoothed.fits'
#reproj

newcontpath=sourcepath+'reprojectedcontinuum.fits'

cntmbeam=radio_beam.Beam.from_fits_header(notreproj_cntmimage[0].header)

if os.path.exists(reprojtexmappath):
    print(f'\nReprojected and smoothed images for source {source} already exist')
    print('Grabbing from respective paths')
    smoothed_trot_fits=fits.open(sourcepath+reprojtexmappath)
    smoothed_trot=smoothed_trot_fits[0].data
    smoothed_troterr=fits.open(sourcepath+reprojtexerrmappath)
    smoothed_snr=fits.open(sourcepath+reprojsnrpath)
    smoothed_ntot=fits.open(sourcepath+reprojntotpath)
    smoothed_ntoterr=fits.open(sourcepath+reprojntoterrpath)
    reprojandsmooth=False
else:
    print('Begin reprojection and smoothing')
    reprojandsmooth=True
    if os.path.exists(newcontpath):
        print(f'Reprojected continuum for {source} already exists.')
        print(f'Loading from {newcontpath}')
        cntmimage=fits.open(newcontpath)
    else:
        newcont,footprint=reproject.reproject_interp((notreproj_cntmimage[0].data.squeeze(),
                                                          WCS(notreproj_cntmimage[0].header).celestial),texmap[0].header)
        newcontprimaryHDU=fits.PrimaryHDU(newcont)
        newcontprimaryHDU.header=texmap[0].header
        newcontprimaryHDU.header['BUNIT']='Jy/beam'
        newcontprimaryHDU.header['BTYPE']='1mm continuum flux'
        newcontprimaryHDU.header['RESTFRQ']=restfreq.to('Hz').value
        newcontprimaryHDU.header['BMAJ']=cntmbeam.major.to('deg').value
        newcontprimaryHDU.header['BMIN']=cntmbeam.minor.to('deg').value
        #newcontprimaryHDU.header['BPA']=cntmbeam.
    
        newcontHDUList=fits.HDUList([newcontprimaryHDU])

        newcontHDUList.writeto(newcontpath,overwrite=True)
        cntmimage=newcontHDUList
    
    beam_tex=radio_beam.Beam.from_fits_header(texmap[0].header)
    beam_cntm=radio_beam.Beam.from_fits_header(notreproj_cntmimage[0].header)
    beam_deconv=beam_cntm.deconvolve(beam_tex)
    pixscale=WCS(texmap[0].header).proj_plane_pixel_area()**0.5
    okmask=np.isfinite(texmap[0].data)*np.isfinite(ntotmap.value)
    
    smoothed_trot = convolve(texmap[0].data, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    cmbmasktrot=np.ma.masked_less_equal(smoothed_trot,Tcmb.value)
    smoothed_trot=cmbmasktrot
    okmasktrot=np.ma.masked_where(okmask==False,smoothed_trot)
    smoothed_trot=okmasktrot.filled(fill_value=np.nan)

    smoothed_troterr = convolve(texerrmap, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    #pdb.set_trace()
    smoothed_snr = smoothed_trot/smoothed_troterr#convolve(snrmap, beam_deconv.as_kernel(pixscale))
    cmbmasksnr=np.ma.masked_where(smoothed_trot < Tcmb.value,smoothed_snr)
    smoothed_snr=cmbmasksnr.filled(fill_value=np.nan)

    smoothed_ntot = convolve(ntotmap, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    cmbmaskntot=np.ma.masked_where(smoothed_trot < Tcmb.value, smoothed_ntot)
    smoothed_ntot=cmbmaskntot.filled(fill_value=np.nan)
    okmaskntot=np.ma.masked_where(okmask == False, smoothed_ntot)
    smoothed_ntot=okmaskntot.filled(fill_value=np.nan)

    smoothed_ntoterr=convolve(ntoterrmap, beam_deconv.as_kernel(pixscale),fill_value=np.nan, nan_treatment='interpolate')
    cmbmaskntoterr=np.ma.masked_where(smoothed_trot < Tcmb.value, smoothed_ntoterr)
    smoothed_ntoterr=cmbmaskntoterr.filled(fill_value=np.nan)
    #pdb.set_trace()
'''
plt.figure(1)
plt.imshow(smoothed_ntoterr,origin='lower',vmin=0,vmax=2e17)
plt.show()

plt.figure(2)
plt.imshow(smoothed_troterr,origin='lower',vmin=0,vmax=300)
plt.show()

pdb.set_trace()
'''
texmapdata=smoothed_trot*u.K#texmap[0].data*u.K

cntmdata=cntmimage[0].data*u.Jy
cntmdatasqee=np.squeeze(cntmdata)
cntmstd=0.2*u.mJy#0.025*u.mJy#np.nanstd(cntmdatasqee)#Taken from 2016 ALMA Proposal here
cntmkappa=kappa(restfreq)
cntmcell=cntmimage[0].header['CDELT2']*u.deg
bmaj=cntmbeam.major#notreproj_cntmimage[0].header['BMAJ']*u.deg#Assert original continuum beam 2/2/2022
bmajpix=round((bmaj/cntmcell).value)
bmin=cntmbeam.minor#notreproj_cntmimage[0].header['BMIN']*u.deg#Assert original continuum beam 2/2/2022
bminpix=round((bmin/cntmcell).value)
beamarea_sr=cntmbeam.sr#np.pi*(bmaj/2)*(bmin/2)#Added pi and divide by 2 on 2/1/2022, they were missing prior to this date
beamarea_sqdeg=cntmbeam.sr.to('deg2')#beamarea_sqdeg.to('sr')
bmajtophyssize=(np.tan(bmaj)*d).to('AU')
bmintophyssize=(np.tan(bmin)*d).to('AU')
'''Can probably simplify beamarea_phys to d(np.tan(bmaj)*np.tan(bmin))'''
beamarea_phys=cntmbeam.beam_projected_area(d).to('AU2')#np.pi*(bmajtophyssize/2)*(bmintophyssize/2)#Added divide by 2 on 2/1/2022, it was missing prior to this date
#calcdomega=np.pi*(bmaj/2)*(bmin/2)

geomeanbeamarea=np.sqrt((bmajtophyssize*bmintophyssize)/2)
print(f'Geometric mean beam radius (continuum): {geomeanbeamarea}')

sys.exit()
    
cntmwcs=WCS(cntmimage[0].header)
texwcs=WCS(texmap[0].header)

equiv=u.brightness_temperature(cntmimage[0].header['RESTFRQ']*u.Hz)

rms_1mm_K=(cntmstd/beamarea_sr).to(u.K,equivalencies=equiv)

snr=3

if reprojandsmooth==True:
    smoothedtrotprimaryhdu=fits.PrimaryHDU(smoothed_trot)
    smoothedtrotprimaryhdu.header=texmap[0].header
    smoothedtrotprimaryhdu.header['BMAJ']=bmaj.to('deg').value
    smoothedtrotprimaryhdu.header['BMIN']=bmaj.to('deg').value
    smoothedtrothdul=fits.HDUList([smoothedtrotprimaryhdu])
    smoothedtrotpath=sourcepath+f'bootstrap_smoothed_trot_to_bolocamfeathercont.fits'
    print(f'Saving smoothed temperature map at {smoothedtrotpath}\n')
    smoothedtrothdul.writeto(smoothedtrotpath,overwrite=True)

    errsmoothedtrotprimaryhdu=fits.PrimaryHDU(smoothed_troterr.value)
    errsmoothedtrotprimaryhdu.header=texmap[0].header
    errsmoothedtrotprimaryhdu.header['BMAJ']=bmaj.to('deg').value
    errsmoothedtrotprimaryhdu.header['BMIN']=bmin.to('deg').value
    errsmoothedtrothdul=fits.HDUList([errsmoothedtrotprimaryhdu])
    errsmoothedtrotpath=sourcepath+f'bootstrap_smoothed_trot_err.fits'
    print(f'Saving smoothed temperature error map at {errsmoothedtrotpath}\n')
    errsmoothedtrothdul.writeto(errsmoothedtrotpath,overwrite=True)

    smoothedntotprimaryhdu=fits.PrimaryHDU(smoothed_ntot.value)
    smoothedntotprimaryhdu.header=texmap[0].header
    smoothedntotprimaryhdu.header['BMAJ']=bmaj.to('deg').value
    smoothedntotprimaryhdu.header['BMIN']=bmaj.to('deg').value
    smoothedntotprimaryhdu.header['BTYPE']='Column density (CH3OH)'
    smoothedntotprimaryhdu.header['BUNIT']='cm-2'
    smoothedntothdul=fits.HDUList([smoothedntotprimaryhdu])
    smoothedntotpath=sourcepath+f'bootstrap_smoothed_ntot_to_bolocamfeathercont.fits'
    print(f'Saving smoothed column density map at {smoothedntotpath}\n')
    smoothedntothdul.writeto(smoothedntotpath,overwrite=True)

    errsmoothedntotprimaryhdu=fits.PrimaryHDU(smoothed_ntoterr.value)
    errsmoothedntotprimaryhdu.header=texmap[0].header
    errsmoothedntotprimaryhdu.header['BMAJ']=bmaj.to('deg').value
    errsmoothedntotprimaryhdu.header['BMIN']=bmin.to('deg').value
    errsmoothedntotprimaryhdu.header['BTYPE']='Column density error (CH3OH)'
    errsmoothedntotprimaryhdu.header['BUNIT']='cm-2'
    errsmoothedntothdul=fits.HDUList([errsmoothedntotprimaryhdu])
    errsmoothedntotpath=sourcepath+f'bootstrap_smoothed_ntot_err.fits'
    print(f'Saving smoothed column density error map at {errsmoothedntotpath}\n')
    errsmoothedntothdul.writeto(errsmoothedntotpath,overwrite=True)


cntmsurfbright=cntmdata/beamarea_sr
err_cntmsurfbright=cntmstd/beamarea_sr

cntmsurfbrightK=cntmsurfbright.to(u.K,equivalencies=equiv)
err_Kcntmsurfbright=err_cntmsurfbright.to(u.K,equivalencies=equiv)

bettergasmass=(((beamarea_sr*d**2*c**2)/(2*k*restfreq**2*cntmkappa))*((cntmsurfbright)/texmapdata)).to('solMass')
err_bettergasmass=np.sqrt((((cntmsurfbright*beamarea_sr*d*c**2)/(k*texmapdata*restfreq**2*cntmkappa))*d_unc)**2+(((beamarea_sr*d**2*c**2)/(2*k*texmapdata*restfreq**2*cntmkappa))*err_cntmsurfbright)**2+(((-cntmsurfbright*beamarea_sr*d**2*c**2)/(2*k*texmapdata**2*restfreq**2*cntmkappa))*smoothed_troterr)**2).to('solMass')
snr_gasmass=(bettergasmass/err_bettergasmass)
masked_gasmass=np.ma.masked_where(snr_gasmass<3,bettergasmass)
sigma3_gasmass=masked_gasmass.filled(fill_value=np.nan)

betterSigmam=(bettergasmass/beamarea_phys).to('g cm-2')
err_betterSigmam=(err_bettergasmass/beamarea_phys).to('g cm-2')
snr_Sigmam=(betterSigmam/err_betterSigmam)
masked_Sigmam=np.ma.masked_where(snr_Sigmam<3,betterSigmam)
sigma3_Sigmam=masked_Sigmam.filled(fill_value=np.nan)

betterNh2=(betterSigmam/mu).to('cm-2')
err_betterNh2=(err_betterSigmam/mu).to('cm-2')
snr_Nh2=(betterNh2/err_betterNh2)
masked_Nh2=np.ma.masked_where(snr_Nh2<3,betterNh2)
sigma3_Nh2=masked_Nh2.filled(fill_value=np.nan)
nh2=sigma3_Nh2
nh2_unc=err_betterNh2

betterXch3oh=(smoothed_ntot/betterNh2).to('')
err_betterXch3oh=np.sqrt(((1/betterNh2)*smoothed_ntoterr)**2+((-smoothed_ntot/(betterNh2**2))*err_betterNh2)**2).to('')
snr_Xch3oh=(betterXch3oh/err_betterXch3oh)
masked_Xch3oh=np.ma.masked_where(snr_Xch3oh<3,betterXch3oh)
sigma3_Xch3oh=masked_Xch3oh.filled(fill_value=np.nan)
ch3ohabundance=betterXch3oh
ch3ohabundance_unc=err_betterXch3oh
ch3ohabundance_snr=snr_Xch3oh
ch3ohabundance_sigclip=sigma3_Xch3oh

bettertau=((c**2*cntmsurfbright)/(2*k*texmapdata*restfreq**2)).decompose()
err_bettertau=np.sqrt(((c**2/(2*k*texmapdata*restfreq**2))*err_cntmsurfbright)**2+(((-cntmsurfbright*c**2)/(2*k*texmapdata**2*restfreq**2))*smoothed_troterr)**2).decompose()
snr_tau=(bettertau/err_bettertau)
masked_tau=np.ma.masked_where(snr_tau<3,bettertau)
sigma3_tau=masked_tau.filled(fill_value=np.nan)
tau=sigma3_tau

betterlum=(4*np.pi*(geomeanbeamarea)**2*sigmasb*cntmsurfbrightK**4).to('solLum')
err_betterlum=np.sqrt(((16*np.pi*geomeanbeamarea**2)*sigmasb*cntmsurfbrightK**3*err_Kcntmsurfbright)**2).to('solLum')
snr_lum=(betterlum/err_betterlum)
masked_lum=np.ma.masked_where(snr_lum<3,betterlum)
sigma3_lum=masked_lum.filled(fill_value=np.nan)
lums=sigma3_lum
lumserr=err_betterlum

'''
ch3ohlum=(4*np.pi*(geomeanbeamarea)**2*sigmasb*(smoothed_trot*u.K)**4).to('solLum')
err_ch3ohlum=np.sqrt(((16*np.pi*geomeanbeamarea**2)*sigmasb*(smoothed_trot*u.K)**3*(smoothed_troterr))**2).to('solLum')
snr_ch3ohlum=(ch3ohlum/err_ch3ohlum)
masked_ch3ohlum=np.ma.masked_where(snr_ch3ohlum<3,ch3ohlum)
sigma3_ch3ohlum=masked_ch3ohlum.filled(fill_value=np.nan)
ch3ohlums=sigma3_ch3ohlum
ch3ohlumserr=err_ch3ohlum
'''
#pdb.set_trace()
  
gasmassprimaryhdu=fits.PrimaryHDU(sigma3_gasmass.value)
gasmassprimaryhdu.header=texmap[0].header
gasmassprimaryhdu.header['BTYPE']='H2 Mass'
gasmassprimaryhdu.header['BUNIT']='solMass'
gasmasshdul=fits.HDUList([gasmassprimaryhdu])
gasmassfitspath=sourcepath+f'bootstrap_h2massmap_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving H2 mass map at {gasmassfitspath}')
gasmasshdul.writeto(gasmassfitspath,overwrite=True)

errgasmassprimaryhdu=fits.PrimaryHDU(err_bettergasmass.value)
errgasmassprimaryhdu.header=texmap[0].header
errgasmassprimaryhdu.header['BTYPE']='H2 Mass error'
errgasmassprimaryhdu.header['BUNIT']='solMass'
errgasmasshdul=fits.HDUList([errgasmassprimaryhdu])
errgasmassfitspath=sourcepath+f'bootstrap_h2massmap_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving H2 mass error map at {errgasmassfitspath}')
errgasmasshdul.writeto(errgasmassfitspath,overwrite=True)

nh2primaryhdu=fits.PrimaryHDU(nh2.value)
nh2primaryhdu.header=texmap[0].header
nh2primaryhdu.header['BTYPE']='N(H2)'
nh2primaryhdu.header['BUNIT']='cm-2'
nh2hdul=fits.HDUList([nh2primaryhdu])
nh2fitspath=sourcepath+f'bootstrap_nh2map_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving N(H2) map at {nh2fitspath}\n')
nh2hdul.writeto(nh2fitspath,overwrite=True)

nh2uncprimaryhdu=fits.PrimaryHDU(nh2_unc.value)
nh2uncprimaryhdu.header=texmap[0].header
nh2uncprimaryhdu.header['BTYPE']='N(H2) error'
nh2uncprimaryhdu.header['BUNIT']='cm-2'
nh2unchdul=fits.HDUList([nh2uncprimaryhdu])
nh2uncfitspath=sourcepath+f'bootstrap_nh2map_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving N(H2) map at {nh2uncfitspath}\n')
nh2unchdul.writeto(nh2uncfitspath,overwrite=True)

ch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance.value)
ch3ohabundprimaryhdu.header=texmap[0].header
ch3ohabundprimaryhdu.header['BTYPE']='Abundance (CH3OH)'
ch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
ch3ohabundhdul=fits.HDUList([ch3ohabundprimaryhdu])
ch3ohabundfitspath=sourcepath+f'bootstrap_ch3ohabundance_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance map at {ch3ohabundfitspath}\n')
ch3ohabundhdul.writeto(ch3ohabundfitspath,overwrite=True)

errch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_unc.value)
errch3ohabundprimaryhdu.header=texmap[0].header
errch3ohabundprimaryhdu.header['BTYPE']='Abundance error (CH3OH)'
errch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
errch3ohabundhdul=fits.HDUList([errch3ohabundprimaryhdu])
errch3ohabundfitspath=sourcepath+f'bootstrap_ch3ohabundance_error_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance error map at {errch3ohabundfitspath}\n')
errch3ohabundhdul.writeto(errch3ohabundfitspath,overwrite=True)

snrch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_snr.value)
snrch3ohabundprimaryhdu.header=texmap[0].header
snrch3ohabundprimaryhdu.header['BTYPE']='Abundance SNR (CH3OH)'
snrch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
snrch3ohabundhdul=fits.HDUList([snrch3ohabundprimaryhdu])
snrch3ohabundfitspath=sourcepath+f'bootstrap_ch3ohabundance_snr_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance SNR map at {snrch3ohabundfitspath}\n')
snrch3ohabundhdul.writeto(snrch3ohabundfitspath,overwrite=True)

sigclipch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_sigclip.value)
sigclipch3ohabundprimaryhdu.header=texmap[0].header
sigclipch3ohabundprimaryhdu.header['BTYPE']='Abundance (CH3OH)'
sigclipch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
sigclipch3ohabundhdul=fits.HDUList([sigclipch3ohabundprimaryhdu])
sigclipch3ohabundfitspath=sourcepath+f'bootstrap_ch3ohabundance_{snr}sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance map at {sigclipch3ohabundfitspath}\n')
sigclipch3ohabundhdul.writeto(sigclipch3ohabundfitspath,overwrite=True)

lumprimaryhdu=fits.PrimaryHDU(lums.value)
lumprimaryhdu.header=texmap[0].header
lumprimaryhdu.header['BTYPE']='Luminosity'
lumprimaryhdu.header['BUNIT']='Lsun'
lumhdul=fits.HDUList([lumprimaryhdu])
lumfitspath=sourcepath+f'bootstrap_boltzmannlum_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving luminosity map at {lumfitspath}\n')
lumhdul.writeto(lumfitspath,overwrite=True)

lumerrprimaryhdu=fits.PrimaryHDU(lumserr.value)
lumerrprimaryhdu.header=texmap[0].header
lumerrprimaryhdu.header['BTYPE']='Luminosity error'
lumerrprimaryhdu.header['BUNIT']='Lsun'
lumerrhdul=fits.HDUList([lumerrprimaryhdu])
lumerrfitspath=sourcepath+f'bootstrap_boltzmannlum_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving luminosity error map at {lumerrfitspath}\n')
lumerrhdul.writeto(lumerrfitspath,overwrite=True)

tauprimaryhdu=fits.PrimaryHDU(tau.value)
tauprimaryhdu.header=texmap[0].header
tauprimaryhdu.header['BTYPE']='Optical depth (tau)'
tauprimaryhdu.header['BUNIT']='unitless'
tauhdul=fits.HDUList([tauprimaryhdu])
taufitspath=sourcepath+f'bootstrap_opticaldepth_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving optical depth map at {taufitspath}\n')
tauhdul.writeto(taufitspath,overwrite=True)

errtauprimaryhdu=fits.PrimaryHDU(err_bettertau.value)
errtauprimaryhdu.header=texmap[0].header
errtauprimaryhdu.header['BTYPE']='Optical depth (tau) error'
errtauprimaryhdu.header['BUNIT']='unitless'
errtauhdul=fits.HDUList([errtauprimaryhdu])
errtaufitspath=sourcepath+f'bootstrap_opticaldepth_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving optical depth error map at {errtaufitspath}\n')
errtauhdul.writeto(errtaufitspath,overwrite=True)
        
plt.imshow(ch3ohabundance_sigclip.value,origin='lower',norm=mpl.colors.LogNorm(),cmap='cividis')
plt.colorbar(pad=0)
plt.show()
