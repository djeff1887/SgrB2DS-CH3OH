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

source='SgrB2S'
print(f'Source: {source}')

sourcedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/Kfield10originals_noexclusions/",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_200K_trial1_noexclusions/','DSVIII':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVIII/Kfield3originals_175K_trial1_noexclusions/','DSIX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/Kfield7originals_150K_trial1_noexclusions/','DSX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSX/Kfield7originals_100K_trial1_noexclusions/'}
sourcepath=sourcedict[source]

texmap=fits.open(sourcepath+'texmap_3sigma_allspw_withnans_weighted.fits')

ntotfits=fits.open(sourcepath+'ntotmap_allspw_withnans_weighted_useintercept_3sigma.fits')
ntotmap=ntotfits[0].data*u.cm**-2
ntoterrfits=fits.open(sourcepath+'ntoterrmap_allspw_withnans_weighted_useintercept.fits')
ntoterrmap=ntoterrfits[0].data*u.cm**-2
ch3ohSigmam_map=ntotmap*molweight_ch3oh

snrfits=fits.open(sourcepath+'texmap_snr_allspw_weighted.fits')
snrmap=snrfits[0].data

texerrfits=fits.open(sourcepath+'texmap_error_allspw_withnans_weighted.fits')
texerrmap=texerrfits[0].data
maskedtexerr=np.ma.masked_where(texerrmap > 500,texerrmap)
texerrmap=maskedtexerr.filled(fill_value=np.nan)*u.K
cmbmasktexerr=np.ma.masked_where(texmap[0].data < Tcmb.value, texerrmap)
texerrmap=cmbmasktexerr.filled(fill_value=np.nan)

reprojtexmappath=sourcepath+'texmap_3sigma_allspw_withnans_weighted_reprojecttobolocamcont_smoothed.fits'
reprojtexerrpath=sourcepath+'ntotmap_allspw_withnans_weighted_useintercept_3sigma_reprojecttobolocamcont_smoothed.fits'
reprojsnrpath=sourcepath+'texmap_snr_allspw_weighted_reprojecttobolocamcont_smoothed.fits'
reprojntotpath=sourcepath+'ntotmap_allspw_withnans_weighted_useintercept_3sigma_reprojecttobolocamcont_smoothed.fits'
reprojntoterrpath=sourcepath+'ntoterrmap_allspw_withnans_weighted_useintercept_reprojecttobolocamcont_smoothed.fits'
#reproj

newcontpath=sourcepath+'reprojectedcontinuum.fits'

if os.path.exists(reprojtexmappath):
    print(f'\nReprojected and smoothed images for source {source} already exist')
    print('Grabbing from respective paths')
    smoothed_trot_fits=fits.open(sourcepath+reprojtexmappath)
    smoothed_trot=smoothed_trot_fits[0].data
    smoothed_troterr=fits.open(sourcepath+reprojtexerrmappath)
    smoothed_snr=fits.open(sourcepath+reprojsnrpath)
    smoothed_ntot=fits.open(sourcepath+reprojntotpath)
    smoothed_ntoterr=fits.open(sourcepath+reprojntoterrpath)
else:
    print('Begin reprojection and smoothing')
    if os.path.exists(newcontpath):
        print(f'Reprojected continuum for {source} already exists.')
        print(f'Loading from {newcontpath}')
        cntmimage=fits.open(newcontpath)
    else:
        newcont,footprint=reproject.reproject_interp((notreproj_cntmimage[0].data.squeeze(),WCS(notreproj_cntmimage[0].header).celestial),texmap[0].header)
        newcontprimaryHDU=fits.PrimaryHDU(newcont)
        newcontprimaryHDU.header=texmap[0].header
        newcontprimaryHDU.header['BUNIT']='Jy/beam'
        newcontprimaryHDU.header['BTYPE']='1mm continuum flux'
        newcontprimaryHDU.header['RESTFRQ']=restfreq.to('Hz').value
    
        newcontHDUList=fits.HDUList([newcontprimaryHDU])

        newcontHDUList.writeto(newcontpath,overwrite=True)
        cntmimage=newcontHDUList

    #pdb.set_trace()
    
    beam_tex=radio_beam.Beam.from_fits_header(texmap[0].header)
    beam_cntm=radio_beam.Beam.from_fits_header(notreproj_cntmimage[0].header)
    beam_deconv=beam_cntm.deconvolve(beam_tex)
    pixscale=WCS(texmap[0].header).proj_plane_pixel_area()**0.5
    smoothed_trot = convolve(texmap[0].data, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    cmbmasktrot=np.ma.masked_less_equal(smoothed_trot,Tcmb.value)
    smoothed_trot=cmbmasktrot#smoothed_trot.filled(fill_value=np.nan)

    smoothed_troterr = convolve(texerrmap, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    #pdb.set_trace()
    smoothed_snr = smoothed_trot/smoothed_troterr#convolve(snrmap, beam_deconv.as_kernel(pixscale))
    cmbmasksnr=np.ma.masked_where(smoothed_trot < Tcmb.value,smoothed_snr)
    smoothed_snr=cmbmasksnr.filled(fill_value=np.nan)

    smoothed_ntot = convolve(ntotmap, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    cmbmaskntot=np.ma.masked_where(smoothed_trot < Tcmb.value,smoothed_ntot)
    smoothed_ntot=cmbmaskntot.filled(fill_value=np.nan)

    smoothed_ntoterr=convolve(ntoterrmap, beam_deconv.as_kernel(pixscale),nan_treatment='interpolate')
    cmbmaskntoterr=np.ma.masked_where(smoothed_trot < Tcmb.value, smoothed_ntoterr)
    smoothed_ntoterr=cmbmaskntoterr.filled(fill_value=np.nan)

#pdb.set_trace()

texmapdata=smoothed_trot*u.K#texmap[0].data*u.K

cntmdata=cntmimage[0].data*u.Jy
cntmdatasqee=np.squeeze(cntmdata)
cntmstd=0.2*u.mJy#0.025*u.mJy#np.nanstd(cntmdatasqee)#Taken from 2016 ALMA Proposal here
cntmkappa=kappa(restfreq)
cntmcell=cntmimage[0].header['CDELT2']*u.deg
bmaj=notreproj_cntmimage[0].header['BMAJ']*u.deg#Assert original continuum beam 2/2/2022
bmajpix=round((bmaj/cntmcell).value)
bmin=notreproj_cntmimage[0].header['BMIN']*u.deg#Assert original continuum beam 2/2/2022
bminpix=round((bmin/cntmcell).value)
beamarea_sqdeg=np.pi*(bmaj/2)*(bmin/2)#Added pi and divide by 2 on 2/1/2022, they were missing prior to this date
beamarea_sr=beamarea_sqdeg.to('sr')
bmajtophyssize=(np.tan(bmaj)*d).to('AU')
bmintophyssize=(np.tan(bmin)*d).to('AU')
'''Can probably simplify beamarea_phys to d(np.tan(bmaj)*np.tan(bmin))'''
beamarea_phys=np.pi*(bmajtophyssize/2)*(bmintophyssize/2)#Added divide by 2 on 2/1/2022, it was missing prior to this date
calcdomega=np.pi*(bmaj/2)*(bmin/2)
    
cntmwcs=WCS(cntmimage[0].header)
texwcs=WCS(texmap[0].header)

equiv=u.brightness_temperature(cntmimage[0].header['RESTFRQ']*u.Hz)

rms_1mm_K=(cntmstd/beamarea_sr).to(u.K,equivalencies=equiv)

snr=3

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

betterlum=(4*np.pi*(bmajtophyssize/2)**2*sigmasb*cntmsurfbrightK**4).to('solLum')
err_betterlum=np.sqrt(((16*np.pi*bmajtophyssize**2)*sigmasb*cntmsurfbrightK**3*err_Kcntmsurfbright)**2).to('solLum')
snr_lum=(betterlum/err_betterlum)
masked_lum=np.ma.masked_where(snr_lum<3,betterlum)
sigma3_lum=masked_lum.filled(fill_value=np.nan)
lums=sigma3_lum
lumserr=err_betterlum

#pdb.set_trace()

'''
for y in range(np.shape(texmapdata)[0]):
    print(f'Start Row {y}')
    for x in range(np.shape(texmapdata)[1]):
        hotspottex=texmapdata[y,x]
        targetpixtoworld=texwcs.all_pix2world([(x,y)],0)#wcs_pix2world(([0,0],texmappix),ra_dec_order=True)
        #pix2world2=texwcs.pixel_to_world([x,y],0)
        #hotspottex=texmapdata[texmappix[0],texmappix[1]]
        #hmm=[[targetpixtoworld.ra[0].value,targetpixtoworld.dec[0].value,0,0],[pix2world2.ra[1].value,pix2world2.dec[1].value,0,0]]
        empty=np.empty((1,4))
        empty[:,0]=targetpixtoworld[:,0]
        empty[:,1]=targetpixtoworld[:,1]
        empty[:,2]=0
        empty[:,3]=0
        cntmworldtopix=cntmwcs.all_world2pix(empty,0,ra_dec_order=True)#world_to_array_index([266.83628632*u.deg,-28.39626318*u.deg,225265023886.0*u.Hz,1*u.m],[0,0,0,0],[1,1,1,1])
        cntmxpix=round(float(cntmworldtopix[:,0]))#math.floor(cntmworldtopix[1][0]-1)
        cntmypix=round(float(cntmworldtopix[:,1]))-1#math.floor(cntmworldtopix[0][1])
        
        hotspotjy=cntmdatasqee[cntmypix,cntmxpix]
        hotspotjysr=hotspotjy/beamarea_sr#The equations we're using start from the Plank function, so this gets us to spectral radiance units
        #ntoterr=(1/snrmap[y,x])*ntotmap[y,x]
        
        if (hotspotjysr*beamarea_sr) < snr*cntmstd or snrmap[y,x] < 3:
            #print('Negative continuum flux detected')
            #print(f'Continuum pixel: ({cntmypix},{cntmxpix})')
            #print(f'Corresponding Texmap pixel: ({y},{x})')
            nh2[y,x]=np.nan
            nh2_unc[y,x]=np.nan
            nh2_snr[y,x]=np.nan
            ch3ohabundance[y,x]=np.nan
            lums[y,x]=np.nan
            lumserr[y,x]=np.nan
            tau[y,x]=np.nan
            ch3ohabundance_unc[y,x]=np.nan
            ch3ohabundance_snr[y,x]=np.nan
            ch3ohabundance_sigclip[y,x]=np.nan
            #print(f'Pixel ({y},{x}) flux < 3*sigma in continuum')
            #pdb.set_trace()
            continue

        #print(f'Results for region in {source}\nTex pixel {texmappix}):\nCntm pixel: {cntmypix,cntmxpix}')
        #print(f'target flux density: {cntmdatasqee[cntmypix,cntmxpix].to("mJy")}, target temperature {hotspottex}')
        jysrtoK=hotspotjysr.to(u.K,equivalencies=equiv)
        targettau=opticaldepth(hotspotjysr,restfreq,hotspottex).decompose()#.to('')# units are disagreeing at the moment, but pretty sure this is the optical depth in spite of the dangling sr-1 term
        targetlum=luminosity(bmajtophyssize,jysrtoK).to("solLum")#per Ginsburg+2017, this uses the peak of the continuum, not methanol
        targetlumerr=error_luminosity(bmajtophyssize,jysrtoK).to('solLum')
        targetgasmass=gasmass(hotspotjysr,beamarea_sr,restfreq,cntmkappa,hotspottex)
        gassmass_unc=gassmasserr(hotspotjysr,beamarea_sr,restfreq,cntmkappa,hotspottex,texerrmap[y,x])
        targetSigmam=((targetgasmass/(np.pi*beamarea_phys)).to('g cm-2'))
        Sigmam_unc=np.abs(gassmass_unc/targetgasmass)*targetSigmam
        targetSigmam_reduced=(targetSigmam/mu).to('cm-2')
        Sigmam_reduced_unc=(np.abs(Sigmam_unc/targetSigmam)*targetSigmam_reduced).to('cm-2')
        Sigmam_reduced_snr=(targetSigmam_reduced/Sigmam_reduced_unc).to('')
        #pdb.set_trace()
        if targetSigmam_reduced < 0:
            print('Negative surface density detected')
            print(f'Continuum flux density: {hotspotjysr}')
            print(f'Brightness: {jysrtoK}')
            print(f'Gass mass: {targetgasmass.to("solMass")}')
            print(f'Sigma_m: {targetSigmam}')
            print(f'Reduced Sigma_m: {targetSigmam_reduced}')
            pdb.set_trace()
        if Sigmam_reduced_snr < snr:
            print(f'Pixel ({y},{x}) Sigma_mass < 3sigma')
            nh2[y,x]=np.nan
            nh2_unc[y,x]=Sigmam_reduced_unc.value
            nh2_snr[y,x]=Sigmam_reduced_snr.value
            ch3ohabundance[y,x]=np.nan
            lums[y,x]=np.nan
            lumserr[y,x]=np.nan
            tau[y,x]=np.nan
            ch3ohabundance_unc[y,x]=np.nan
            ch3ohabundance_snr[y,x]=np.nan
            ch3ohabundance_sigclip[y,x]=np.nan
            continue
        targetch3ohabund=(ntotmap[y,x]/targetSigmam_reduced).to('')
        ch3ohabun_unc=np.sqrt(((1/targetSigmam_reduced)*ntoterrmap[y,x])**2+((-ntotmap[y,x]/(targetSigmam_reduced**2))*Sigmam_reduced_unc)**2).to('')#((1/Sigmam_reduced_snr)*targetch3ohabund).to('')
        ch3ohabun_snr=targetch3ohabund/ch3ohabun_unc
        
        nh2[y,x]=targetSigmam_reduced.value
        nh2_unc[y,x]=Sigmam_reduced_unc.value
        nh2_snr[y,x]=Sigmam_reduced_snr.value
        ch3ohabundance[y,x]=targetch3ohabund.value
        ch3ohabundance_unc[y,x]=ch3ohabun_unc.value
        ch3ohabundance_snr[y,x]=ch3ohabun_snr.value
        lums[y,x]=targetlum.value
        lumserr[y,x]=targetlumerr.value
        tau[y,x]=targettau.value
        
        if ch3ohabun_snr >= 3:
            ch3ohabundance_sigclip[y,x]=targetch3ohabund.value
        else:
            ch3ohabundance_sigclip[y,x]=np.nan
   
'''     
gasmassprimaryhdu=fits.PrimaryHDU(sigma3_gasmass.value)
gasmassprimaryhdu.header=texmap[0].header
gasmassprimaryhdu.header['BTYPE']='H2 Mass'
gasmassprimaryhdu.header['BUNIT']='solMass'
gasmasshdul=fits.HDUList([gasmassprimaryhdu])
gasmassfitspath=sourcepath+f'h2massmap_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving H2 mass map at {gasmassfitspath}')
gasmasshdul.writeto(gasmassfitspath,overwrite=True)

errgasmassprimaryhdu=fits.PrimaryHDU(err_bettergasmass.value)
errgasmassprimaryhdu.header=texmap[0].header
errgasmassprimaryhdu.header['BTYPE']='H2 Mass error'
errgasmassprimaryhdu.header['BUNIT']='solMass'
errgasmasshdul=fits.HDUList([errgasmassprimaryhdu])
errgasmassfitspath=sourcepath+f'h2massmap_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving H2 mass error map at {errgasmassfitspath}')
errgasmasshdul.writeto(errgasmassfitspath,overwrite=True)

nh2primaryhdu=fits.PrimaryHDU(nh2.value)
nh2primaryhdu.header=texmap[0].header
nh2primaryhdu.header['BTYPE']='N(H2)'
nh2primaryhdu.header['BUNIT']='cm-2'
nh2hdul=fits.HDUList([nh2primaryhdu])
nh2fitspath=sourcepath+f'nh2map_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving N(H2) map at {nh2fitspath}\n')
nh2hdul.writeto(nh2fitspath,overwrite=True)

nh2uncprimaryhdu=fits.PrimaryHDU(nh2_unc.value)
nh2uncprimaryhdu.header=texmap[0].header
nh2uncprimaryhdu.header['BTYPE']='N(H2) error'
nh2uncprimaryhdu.header['BUNIT']='cm-2'
nh2unchdul=fits.HDUList([nh2uncprimaryhdu])
nh2uncfitspath=sourcepath+f'nh2map_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving N(H2) map at {nh2uncfitspath}\n')
nh2unchdul.writeto(nh2uncfitspath,overwrite=True)

ch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance.value)
ch3ohabundprimaryhdu.header=texmap[0].header
ch3ohabundprimaryhdu.header['BTYPE']='Abundance (CH3OH)'
ch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
ch3ohabundhdul=fits.HDUList([ch3ohabundprimaryhdu])
ch3ohabundfitspath=sourcepath+f'ch3ohabundance_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance map at {ch3ohabundfitspath}\n')
ch3ohabundhdul.writeto(ch3ohabundfitspath,overwrite=True)

errch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_unc.value)
errch3ohabundprimaryhdu.header=texmap[0].header
errch3ohabundprimaryhdu.header['BTYPE']='Abundance error (CH3OH)'
errch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
errch3ohabundhdul=fits.HDUList([errch3ohabundprimaryhdu])
errch3ohabundfitspath=sourcepath+f'ch3ohabundance_error_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance error map at {errch3ohabundfitspath}\n')
errch3ohabundhdul.writeto(errch3ohabundfitspath,overwrite=True)

snrch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_snr.value)
snrch3ohabundprimaryhdu.header=texmap[0].header
snrch3ohabundprimaryhdu.header['BTYPE']='Abundance SNR (CH3OH)'
snrch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
snrch3ohabundhdul=fits.HDUList([snrch3ohabundprimaryhdu])
snrch3ohabundfitspath=sourcepath+f'ch3ohabundance_snr_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance SNR map at {snrch3ohabundfitspath}\n')
snrch3ohabundhdul.writeto(snrch3ohabundfitspath,overwrite=True)

sigclipch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance_sigclip.value)
sigclipch3ohabundprimaryhdu.header=texmap[0].header
sigclipch3ohabundprimaryhdu.header['BTYPE']='Abundance (CH3OH)'
sigclipch3ohabundprimaryhdu.header['BUNIT']='cm-2'
sigclipch3ohabundhdul=fits.HDUList([sigclipch3ohabundprimaryhdu])
sigclipch3ohabundfitspath=sourcepath+f'ch3ohabundance_{snr}sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving CH3OH abundance map at {sigclipch3ohabundfitspath}\n')
sigclipch3ohabundhdul.writeto(sigclipch3ohabundfitspath,overwrite=True)

lumprimaryhdu=fits.PrimaryHDU(lums.value)
lumprimaryhdu.header=texmap[0].header
lumprimaryhdu.header['BTYPE']='Luminosity'
lumprimaryhdu.header['BUNIT']='Lsun'
lumhdul=fits.HDUList([lumprimaryhdu])
lumfitspath=sourcepath+f'boltzmannlum_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving luminosity map at {lumfitspath}\n')
lumhdul.writeto(lumfitspath,overwrite=True)

lumerrprimaryhdu=fits.PrimaryHDU(lumserr.value)
lumerrprimaryhdu.header=texmap[0].header
lumerrprimaryhdu.header['BTYPE']='Luminosity error'
lumerrprimaryhdu.header['BUNIT']='Lsun'
lumerrhdul=fits.HDUList([lumerrprimaryhdu])
lumerrfitspath=sourcepath+f'boltzmannlum_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving luminosity error map at {lumerrfitspath}\n')
lumerrhdul.writeto(lumerrfitspath,overwrite=True)

tauprimaryhdu=fits.PrimaryHDU(tau.value)
tauprimaryhdu.header=texmap[0].header
tauprimaryhdu.header['BTYPE']='Optical depth (tau)'
tauprimaryhdu.header['BUNIT']='unitless'
tauhdul=fits.HDUList([tauprimaryhdu])
taufitspath=sourcepath+f'opticaldepth_{snr}sigma_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving optical depth map at {taufitspath}\n')
tauhdul.writeto(taufitspath,overwrite=True)

errtauprimaryhdu=fits.PrimaryHDU(err_bettertau.value)
errtauprimaryhdu.header=texmap[0].header
errtauprimaryhdu.header['BTYPE']='Optical depth (tau) error'
errtauprimaryhdu.header['BUNIT']='unitless'
errtauhdul=fits.HDUList([errtauprimaryhdu])
errtaufitspath=sourcepath+f'opticaldepth_error_bolocamfeather_smoothedtobolocam.fits'
print(f'Saving optical depth error map at {errtaufitspath}\n')
errtauhdul.writeto(errtaufitspath,overwrite=True)
        
plt.imshow(ch3ohabundance_sigclip.value,origin='lower',norm=mpl.colors.LogNorm(),cmap='cividis')
plt.colorbar(pad=0)
plt.show()
