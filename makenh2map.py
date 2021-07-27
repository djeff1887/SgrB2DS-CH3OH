import numpy as np
from astropy.io import fits
import astropy.units as u
import scipy.constants as cnst
from astropy.wcs import WCS
import math
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.interactive(True)

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

def soundspeed(m,T):
    return np.sqrt((k*T)/m)

def jeansmass(cs,rho):
    return (np.pi/6)*cs**3*G**(-3/2)*rho**(-1/2)

def jeanslength(cs,rho):
    return cs*G**(-1/2)*rho**(-1/2)

def tff(rho):#Free=fall time of gas
    return np.sqrt((3*np.pi)/(32*G*rho))
    
def opticaldepth(snu,nu,tex):
    return (snu*c**2)/(2*k*tex*nu**2)
    
c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
d=8.34*u.kpc#Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
d_unc=0.16*u.kpc#Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
sigmasb=cnst.sigma*u.W/(u.m**2*u.K**4)
mu=2.8*u.Dalton
G=cnst.G*u.m**3/(u.kg*u.s**2)
molweight_ch3oh=(32.042*u.u).to('g')#https://pubchem.ncbi.nlm.nih.gov/compound/methanol

cntminfile='/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter.image.tt0.pbcor.fits'
cntmimage=fits.open(cntminfile)
print(f'Continuum image: {cntminfile}')
restfreq=cntmimage[0].header['RESTFRQ']*u.Hz
source='SgrB2S'
print(f'Source: {source}')

#texmapdict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/texmap_5transmask_3sigma_allspw_withnans_weighted.fits",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial5_newexclusions3andexclusionsnotinfit/texmap_5transmask_3sigma_allspw_withnans_weighted.fits",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"}
sourcedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/Kfield10originals_noexclusions/",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_150K_trial2/"}
sourcepath=sourcedict[source]

texmap=fits.open(sourcepath+'texmap_5transmask_3sigma_allspw_withnans_weighted.fits')

ntotmap=fits.getdata(sourcepath+'ntotmap_allspw_withnans_weighted.fits')*u.cm**-2
ch3ohSigmam_map=ntotmap*molweight_ch3oh

snrmapdict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/"}
snrmap=fits.getdata(sourcepath+'texmap_snr_allspw_weighted.fits')

texerrmapdict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/"}
texerrmap=fits.getdata(sourcepath+'texmap_error_allspw_withnans_weighted.fits')*u.K

texmapdata=texmap[0].data*u.K

cntmdata=cntmimage[0].data*u.Jy
cntmdatasqee=np.squeeze(cntmdata)
cntmstd=0.025*u.mJy#np.nanstd(cntmdatasqee)#Taken from 2016 ALMA Proposal here
cntmkappa=kappa(restfreq)
cntmcell=cntmimage[0].header['CDELT2']*u.deg
bmaj=cntmimage[0].header['BMAJ']*u.deg
bmajpix=round((bmaj/cntmcell).value)
bmin=cntmimage[0].header['BMIN']*u.deg
bminpix=round((bmin/cntmcell).value)
beamarea_sqdeg=bmaj*bmin
beamarea_sr=beamarea_sqdeg.to('sr')
bmajtophyssize=(np.tan(bmaj)*d).to('AU')
bmintophyssize=(np.tan(bmin)*d).to('AU')
'''Can probably simplify beamarea_phys to d(np.tan(bmaj)*np.tan(bmin))'''
beamarea_phys=np.pi*bmajtophyssize*bmintophyssize
calcdomega=np.pi*bmaj*bmin
    
cntmwcs=WCS(cntmimage[0].header)
texwcs=WCS(texmap[0].header)

equiv=u.brightness_temperature(cntmimage[0].header['RESTFRQ']*u.Hz)

nh2=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
nh2_unc=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
nh2_snr=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
ch3ohabundance=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
ch3ohabundance_unc=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
ch3ohabundance_snr=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
lums=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))
tau=np.empty((np.shape(texmapdata)[0],np.shape(texmapdata)[1]))

snr=3

for y in range(np.shape(texmapdata)[0]):
    print(f'Start Row {y}')
    for x in range(np.shape(texmapdata)[1]):
        hotspottex=texmapdata[y,x]
        targetpixtoworld=texwcs.pixel_to_world([0,0],[y,x])#wcs_pix2world(([0,0],texmappix),ra_dec_order=True)
        pix2world2=texwcs.pixel_to_world([y,x],[0,0])
        #hotspottex=texmapdata[texmappix[0],texmappix[1]]
        hmm=[[targetpixtoworld.ra[0].value,targetpixtoworld.dec[0].value,0,0],
             [pix2world2.ra[1].value,pix2world2.dec[1].value,0,0]]
        cntmworldtopix=cntmwcs.wcs_world2pix(hmm,1,ra_dec_order=True)#world_to_array_index([266.83628632*u.deg,-28.39626318*u.deg,225265023886.0*u.Hz,1*u.m],[0,0,0,0],[1,1,1,1])
        cntmxpix=math.floor(cntmworldtopix[1][0]-1)
        cntmypix=math.floor(cntmworldtopix[0][1])
        
        hotspotjy=cntmdatasqee[cntmypix,cntmxpix]
        hotspotjysr=hotspotjy/beamarea_sr#The equations we're using start from the Plank function, so this gets us to spectral radiance units
        #ntoterr=(1/snrmap[y,x])*ntotmap[y,x]
        
        
        if (hotspotjysr*beamarea_sr) < snr*cntmstd or snrmap[y,x] < 3:
            #print('Negative continuum flux detected')
            #print(f'Continuum pixel: ({cntmypix},{cntmxpix})')
            #print(f'Corresponding Texmap pixel: ({y},{x})')
            nh2[y,x]=np.nan
            nh2_unc[y,x]=np.nan
            ch3ohabundance[y,x]=np.nan
            lums[y,x]=np.nan
            tau[y,x]=np.nan
            continue

        #print(f'Results for region in {source}\nTex pixel {texmappix}):\nCntm pixel: {cntmypix,cntmxpix}')
        #print(f'target flux density: {cntmdatasqee[cntmypix,cntmxpix].to("mJy")}, target temperature {hotspottex}')
        jysrtoK=hotspotjysr.to(u.K,equivalencies=equiv)
        targettau=opticaldepth(hotspotjy,restfreq,hotspottex).to('')
        targetlum=luminosity(bmajtophyssize,jysrtoK).to("solLum")#per Ginsburg+2017, this uses the peak of the continuum, not methanol
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
            print(f'Pixel ({y},{x}) Sigma_m < 3sigma')
            nh2[y,x]=np.nan
            nh2_unc[y,x]=Sigmam_reduced_unc.value
            nh2_snr[y,x]=Sigmam_reduced_snr.value
            ch3ohabundance[y,x]=np.nan
            lums[y,x]=np.nan
            tau[y,x]=np.nan
            
            continue
        targetch3ohabund=(ntotmap[y,x]/targetSigmam_reduced).to('')
        ch3ohabun_unc=((1/Sigmam_reduced_snr)*targetch3ohabund).to('')
        ch3ohabun_snr=targetch3ohabund/ch3ohabun_unc
        
        nh2[y,x]=targetSigmam_reduced.value
        nh2_unc[y,x]=Sigmam_reduced_unc.value
        nh2_snr[y,x]=Sigmam_reduced_snr.value
        ch3ohabundance[y,x]=targetch3ohabund.value
        ch3ohabundance_unc[y,x]=ch3ohabun_unc.value
        ch3ohabundance_snr[y,x]=ch3ohabun_snr.value
        lums[y,x]=targetlum.value
        tau[y,x]=targettau.value
        
nh2primaryhdu=fits.PrimaryHDU(nh2)
nh2primaryhdu.header=texmap[0].header
nh2primaryhdu.header['BTYPE']='N(H2)'
nh2primaryhdu.header['BUNIT']='cm-2'
nh2hdul=fits.HDUList([nh2primaryhdu])
nh2fitspath=sourcepath+f'nh2map_{snr}sigmacontandsurfacedensity.fits'
print(f'Saving N(H2) map at {nh2fitspath}\n')
nh2hdul.writeto(nh2fitspath,overwrite=True)

nh2uncprimaryhdu=fits.PrimaryHDU(nh2_unc)
nh2uncprimaryhdu.header=texmap[0].header
nh2uncprimaryhdu.header['BTYPE']='N(H2) error'
nh2uncprimaryhdu.header['BUNIT']='cm-2'
nh2unchdul=fits.HDUList([nh2uncprimaryhdu])
nh2uncfitspath=sourcepath+f'nh2map_error.fits'
print(f'Saving N(H2) map at {nh2uncfitspath}\n')
nh2unchdul.writeto(nh2uncfitspath,overwrite=True)

ch3ohabundprimaryhdu=fits.PrimaryHDU(ch3ohabundance)
ch3ohabundprimaryhdu.header=texmap[0].header
ch3ohabundprimaryhdu.header['BTYPE']='Abundance (CH3OH)'
ch3ohabundprimaryhdu.header['BUNIT']='cm-2/cm-2'
ch3ohabundhdul=fits.HDUList([ch3ohabundprimaryhdu])
ch3ohabundfitspath=sourcepath+f'ch3ohabundance.fits'
print(f'Saving CH3OH abundance map at {ch3ohabundfitspath}\n')
ch3ohabundhdul.writeto(ch3ohabundfitspath,overwrite=True)

lumprimaryhdu=fits.PrimaryHDU(lums)
lumprimaryhdu.header=texmap[0].header
lumprimaryhdu.header['BTYPE']='Luminosity'
lumprimaryhdu.header['BUNIT']='Lsun'
lumhdul=fits.HDUList([lumprimaryhdu])
lumfitspath=sourcepath+f'boltzmannlum.fits'
print(f'Saving luminosity map at {lumfitspath}\n')
lumhdul.writeto(lumfitspath,overwrite=True)

tauprimaryhdu=fits.PrimaryHDU(tau)
tauprimaryhdu.header=texmap[0].header
tauprimaryhdu.header['BTYPE']='Optical depth (tau)'
tauprimaryhdu.header['BUNIT']='unitless'
tauhdul=fits.HDUList([tauprimaryhdu])
taufitspath=sourcepath+f'opticaldepth.fits'
print(f'Saving optical depth map at {taufitspath}\n')
tauhdul.writeto(taufitspath,overwrite=True)
        
plt.imshow(ch3ohabundance,origin='lower',norm=mpl.colors.LogNorm(),cmap='cividis')
plt.colorbar(pad=0)
plt.show()