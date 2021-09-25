import numpy as np
from astropy.io import fits
import astropy.units as u
import scipy.constants as cnst
from astropy.wcs import WCS
import math
import pdb

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
    return(snu*beamarea_sr*d**2*c**2)/(2*k*tex*nu**2*kappa)

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


c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
d=8.34*u.kpc##Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
sigmasb=cnst.sigma*u.W/(u.m**2*u.K**4)
mu=2.8*u.Dalton
G=cnst.G*u.m**3/(u.kg*u.s**2)

cntminfile='/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter.image.tt0.pbcor.fits'
cntmimage=fits.open(cntminfile)
print(f'Continuum image: {cntminfile}')
source='DSii'

texmapdict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/texmap_5transmask_3sigma_allspw_withnans_weighted.fits",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial5_newexclusions3andexclusionsnotinfit/texmap_5transmask_3sigma_allspw_withnans_weighted.fits",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial3_8_6-8_7excluded/texmap_3sigma_allspw_withnans_weighted.fits"}
texmap=fits.open(texmapdict[source])

texmapdata=texmap[0].data*u.K

cntmdata=cntmimage[0].data*u.Jy
cntmdatasqee=np.squeeze(cntmdata)
restfreq=cntmimage[0].header['RESTFRQ']*u.Hz
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

'''Round down for both x and y when converting from DS9 (I think)'''
#texmappix=[62,62]#DSVI cont peak
texmappix=[22,24]#DSii cont peak
#texmappix=[36,42]#DSi cont peak
#texmappix=[36,43]#DSi field10 hotspot
#texmappix=[35,44]#DSi hotspot
#texmappix=[69,58]#South-of-hotspot SgrB2S cont peak
#texmappix=[59,62]#Southern molecular ridge SgrB2S cont peak
#texmappix=[73,54]#post statcont SgrB2S hotspot
#texmappix=[73,56]#SgrB2S hotspot
#texmappix=[69,58]#SgrB2S hotspot-adjacent
#texmappix=[79,63]#SgrB2S HII hotpost

hotspottex=texmapdata[texmappix[0],texmappix[1]]
targetpixtoworld=texwcs.pixel_to_world([0,0],texmappix)#wcs_pix2world(([0,0],texmappix),ra_dec_order=True)
pix2world2=texwcs.pixel_to_world(texmappix,[0,0])
hotspottex=texmapdata[texmappix[0],texmappix[1]]
#hmm=[[0,0,0,0],[266.8315372,-28.3972150,0,0]]#DSi
#hmm=[[0,0,0,0],[266.8354046,-28.3960255,0,0]]#SgrB2S
#hmm=[[0,0,0,0],[266.8353844,-28.3960720,0,0]]#SgrB2S hotspot-adjacent
hmm=[[targetpixtoworld.ra[0].value,targetpixtoworld.dec[0].value,0,0],
     [pix2world2.ra[1].value,pix2world2.dec[1].value,0,0]]
#pdb.set_trace()
cntmworldtopix=cntmwcs.wcs_world2pix(hmm,1,ra_dec_order=True)#world_to_array_index([266.83628632*u.deg,-28.39626318*u.deg,225265023886.0*u.Hz,1*u.m],[0,0,0,0],[1,1,1,1])
cntmxpix=math.floor(cntmworldtopix[1][0]-1)
cntmypix=math.floor(cntmworldtopix[0][1])

def circle(data,ycenter,xcenter):
    edgex=[]
    edgey=[]
    cntminbeam=[]
    print('Begin scan for continuum values inside beam')
    for y in range(np.shape(cntmdatasqee)[0]):
        for x in range(np.shape(cntmdatasqee)[1]):
            if ((y-ycenter)**2/(bminpix**2))+((x-xcenter)**2/(bmajpix**2)) <= 1**2:
                cntminbeam.append(cntmdatasqee[y,x].value)
                edgex.append(x)
                edgey.append(y)
                #centrtopix.append((np.sqrt((y-texpeakpix[0])**2+(x-texpeakpix[1])**2)*pixtophysicalsize).value)
            else:
                pass
    
    return np.vstack((edgex,edgey)),cntminbeam

#cntmpixinbeam,cntmbeam=circle(cntmdatasqee,cntmypix,cntmxpix)
hotspotjy=cntmdatasqee[cntmypix,cntmxpix]/beamarea_sr#The equations we're using start from the Plank function, so this gets us to spectral radiance units
#hotspotjy=np.sum(cntmbeam)*u.Jy/beamarea_sr
print(f'Results for region in {source}\nTex pixel {texmappix}):\nCntm pixel: {cntmypix,cntmxpix}')
print(f'target flux density: {cntmdatasqee[cntmypix,cntmxpix].to("mJy")}, target temperature {hotspottex}')
equiv=u.brightness_temperature(cntmimage[0].header['RESTFRQ']*u.Hz)
jysrtoK=hotspotjy.to(u.K,equivalencies=equiv)
targetlum=luminosity(bmajtophyssize,jysrtoK).to("solLum")#per Ginsburg+2017, this uses the peak of the continuum, not methanol
targetgasmass=gasmass(hotspotjy,beamarea_sr,restfreq,cntmkappa,hotspottex)
targetcs=soundspeed(mu, hotspottex)
targetdustmass=targetgasmass/100
targetsigmam=((targetgasmass/(beamarea_phys)).to('g cm-2'))
targetrho=((targetgasmass/((4/3)*np.pi*bmajtophyssize**3)).to('g cm-3'))
tarMj=jeansmass(targetcs,targetrho)
tarlambdaj=jeanslength(targetcs,targetrho)
tartff=tff(targetrho)
print("Dust mass: {:.4}".format(targetdustmass.to('solMass')))
print('Gas mass: {:.6}'.format((targetgasmass).to("solMass")))
print('Sound speed: {:.4}'.format(targetcs.to("km s-1")))
print('Gas surface density: {:.4}'.format((targetsigmam/mu).to("cm-2")))
print('Jeans mass: {:.4}'.format(tarMj.to("solMass")))
print('Jeans radius: {:.6}'.format((tarlambdaj/2).to("AU")))
print('Free-fall time: {:.4}'.format(tartff.to("kyr")))
print('Boltzmann luminosity: {:.4}'.format(targetlum))