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

plt.rcParams['figure.dpi']=300

mu=2.8*u.Dalton

def circle(data,ycenter,xcenter,rad):
    edgex=[]
    edgey=[]
    for i in range(np.shape(data)[0]):
        for j in range(np.shape(data)[1]):
            if (j-xcenter)**2+(i-ycenter)**2==rad**2:
                edgex.append(j)
                edgey.append(i)
    return np.vstack((edgex,edgey))

def submid_profile(r):
    n=0.75
    c=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**n
    return c*((r/pixtophysicalsize)**-n)
    
def sublinear_profile(r):
    n=0.5
    c=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**n
    return c*((r/pixtophysicalsize)**-n)

def linear_profile(r):
    n=1
    c=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**n
    return c*((r/pixtophysicalsize)**-n)
    
def quadratic_profile(r):
    n=1.25
    c=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**n
    return c*((r/pixtophysicalsize)**-n)

def powerlaw_profile(x,a,n):
    return a*((x/pixtophysicalsize.value)**-n)

def round_to_1(x):
    return round(x, -int(math.floor(math.log10(abs(x)))))
    
source="SgrB2S"#os.getenv('SOURCE')#
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}
home=base+homedict[source]
fighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'
figpath=fighome+homedict[source]

dataversion='pacman_sep2023revolution'#homedict[source].replace('/','')
datadir=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'

if not os.path.exists(figpath):
    os.makedirs(figpath)
    print(f'Creating figpath {figpath}')
else:
    print(f'Figpath {figpath} already exists.')
    
if not os.path.exists(datadir):
    print(f'Creating data directory {datadir}')
    os.mkdir(datadir)
    print('Done')
else:
    print(f'Data directory {datadir} already exists.')

#notransmask=['DSv','DSVI','DSVII','DSVIII','DSIX']
#if source == 'SgrB2S':
texmap=home+"bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
#else:
#texmap=home+"texmap_5transmask_3sigma_allspw_withnans_weighted.fits"

texerrmap=home+'test_error_trot_boostrap10000_nonegativeslope.fits'#'error_trot_boostrap1000_nonegativeslope.fits'#"texmap_error_allspw_withnans_weighted.fits"
#snrmap=home+"texmap_snr_allspw_weighted.fits"
abunmap=home+'bootstrap_ch3ohabundance_3sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'#'bootstrap_ch3ohabundance_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
#ntotnh2_abunerrpath=home+'bootstrap_ch3ohabundance_error_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
abunsnrmap=home+'bootstrap_ch3ohabundance_snr_ntotintercept_bolocamfeather_smoothedtobolocam.fits'#'bootstrap_ch3ohabundance_snr_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
nh2map=home+'bootstrap_nh2map_3sigma_bolocamfeather_smoothedtobolocam.fits'
nh2errormap=home+'bootstrap_nh2map_error_bolocamfeather_smoothedtobolocam.fits'
lummap=home+'bootstrap_boltzmannlum_3sigma_bolocamfeather_smoothedtobolocam.fits'
lumerrmap=home+'bootstrap_boltzmannlum_error_bolocamfeather_smoothedtobolocam.fits'
ntotmap=home+'bootstrap_smoothed_ntot_to_bolocamfeathercont_3sigma.fits'
ntoterrmap=home+'bootstrap_smoothed_ntot_err.fits'
h2massmap=home+'bootstrap_h2massmap_3sigma_bolocamfeather_smoothedtobolocam.fits'
h2masserrmap=home+'bootstrap_h2massmap_error_bolocamfeather_smoothedtobolocam.fits'
smoothedtroterrmap=home+'bootstrap_smoothed_trot_err.fits'
smoothedtrotmap=home+'bootstrap_smoothed_trot_to_bolocamfeathercont.fits'

'Removed stopgap during initial qrot fix'
texmap=fits.open(texmap)
texmapdata=texmap[0].data*u.K
texerrdata=np.squeeze(fits.getdata(texerrmap))*u.K
snrs=(texmapdata/texerrdata).value#fits.getdata(snrmap)
abunds=np.squeeze(fits.getdata(abunmap))
#stopgap_errabun=np.squeeze(fits.getdata(ntotnh2_abunerrpath))
snr_abund=fits.getdata(abunsnrmap)#abunds/stopgap_errabun#fits.getdata(abunsnrmap)
nh2s=fits.getdata(nh2map)*u.cm**-2
nh2s_error=fits.getdata(nh2errormap)*u.cm**-2
lums=fits.getdata(lummap)*u.solLum
lumserr=fits.getdata(lumerrmap)*u.solLum
ntots=fits.getdata(ntotmap)*u.cm**-2
ntoterr=np.squeeze(fits.getdata(ntoterrmap))*u.cm**-2
h2mass=fits.getdata(h2massmap)*u.solMass
h2masserr=fits.getdata(h2masserrmap)*u.solMass
smooth_trotfits=fits.open(smoothedtrotmap)
smooth_trot=smooth_trotfits[0].data*u.K
smooth_trot_err=fits.getdata(smoothedtroterrmap)*u.K

'''These used to say pixmask.cutout, but that was a mistake that didn't appropriately mask things outside of the pacman mask'''
if source == 'SgrB2S':
    wcsobj=WCS(smooth_trotfits[0].header)

    regs = regions.Regions.read('/blue/adamginsburg/d.jeff/imaging_results/regfiles/roughsgrb2smassregion_ignoresHIIregion.reg')
    pixreg = regs[0].to_pixel(wcsobj)
    pixmask = pixreg.to_mask()
    
    texmapdata=pixmask.multiply(texmapdata,fill_value=np.nan)
    '''
    texmask=np.ma.masked_where(texmapdata > 1000, texmapdata)
    texmapdata=texmask.filled()
    '''
    texerrdata=pixmask.multiply(texerrdata,fill_value=np.nan)
    snrs=np.squeeze(texmapdata/texerrdata)
    abunds=pixmask.multiply(abunds,fill_value=np.nan)
    snr_abund=pixmask.multiply(snr_abund,fill_value=np.nan)
    nh2s=pixmask.multiply(nh2s,fill_value=np.nan)
    nh2s_error=pixmask.multiply(nh2s_error,fill_value=np.nan)
    lums=pixmask.multiply(lums,fill_value=np.nan)
    lumserr=pixmask.multiply(lumserr,fill_value=np.nan)
    ntots=pixmask.multiply(ntots,fill_value=np.nan)
    ntoterr=pixmask.multiply(ntoterr,fill_value=np.nan)
    h2mass=pixmask.multiply(h2mass,fill_value=np.nan)
    h2masserr=pixmask.multiply(h2masserr,fill_value=np.nan)
    smooth_trot=pixmask.multiply(smooth_trot,fill_value=np.nan)
    smooth_trot_err=pixmask.multiply(smooth_trot_err,fill_value=np.nan)

dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

cntmbeam=radio_beam.Beam.from_fits_header(smooth_trotfits[0].header)
trotbeam=radio_beam.Beam.from_fits_header(texmap[0].header)

cellsize=(np.abs(texmap[0].header['CDELT1']*u.deg)).to('arcsec')
pixperbeam=(cntmbeam.sr/((cellsize**2).to('sr'))).value

pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')

print(pixtophysicalsize)

bmaj=trotbeam.major#texmap[0].header['BMAJ']*u.deg
bmajpix=round((bmaj/cellsize).value)
bmin=trotbeam.minor#texmap[0].header['BMIN']*u.deg
bminpix=round((bmin/cellsize).value)
beamarea_sr=trotbeam.sr#bmaj*bmin
beamarea_sqdeg=beamarea_sr.to('deg2')
bmajtophyssize=(np.tan(bmaj)*dGC).to('AU')
bmintophyssize=(np.tan(bmin)*dGC).to('AU')
'''Can probably simplify beamarea_phys to d(np.tan(bmaj)*np.tan(bmin))'''
beamarea_phys=trotbeam.beam_projected_area(dGC)#np.pi*bmajtophyssize*bmintophyssize

#pdb.set_trace()

cntmbmaj=cntmbeam.major#3.629587176773e-05*u.deg
cntmbmajtoAU=(np.tan(cntmbmaj)*dGC).to('AU')

pixdict={'SgrB2S':(26,14),'DSi':(36,40),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}#y,x; DSiii was 24,24;S was 73,54 - contfix S was 69,58/'DSi':(36,42)
sgrb2scentralpix=(25,25)#contfix pix, prefinal was (66,70)/
#SgrB2S tpeak is 73,54
#nh2dict={'DSiii':(27,27)}
#ntotdict={'SgrB2S':(63,71)'}

texpeakpix=pixdict[source]
#nh2peakpix=nh2dict[source]
#(36,43)#DSi hotspot
#texpeakpix=(73,56)#SgrB2S hotspot
#x-1, y-1 from DS9

print(f'Center p: {texmapdata[texpeakpix[0],texpeakpix[1]]}')

#r=35 #for 15,000 AU
#pixradius=math.ceil((0.08*u.pc/pixtophysicalsize).to(''))
rdict={'SgrB2S':9500*u.AU,'DSi':8000*u.AU,'DSii':5500*u.AU,'DSiii':6000*u.AU,'DSiv':6000*u.AU,'DSv':7000*u.AU,'DSVII':6000*u.AU,'DSVIII':5400*u.AU,'DSIX':8000*u.AU}#1-6400,4-8500,5-4000,7-6600,S-12000
rdictkeys=rdict.keys()
if source not in rdictkeys:
    r_phys=10000*u.AU
else:
    r_phys=rdict[source]#r*pixtophysicalsize.to('pc')
r=math.ceil((r_phys/pixtophysicalsize).to(''))
print(f'physical radius: {r_phys}')

xpixs=np.arange(texpeakpix[0],(texpeakpix[0]+r))

texinradius=[]
xinradius=[]
yinradius=[]
centrtopix=[]
snrsinradius=[]
abundinradius=[]
abundsnrinradius=[]
nh2inradius=[]
lumsinradius=[]
massinradius=[]
ntotsinradius=[]
lookformax=[]

yy,xx=np.indices(texmapdata.shape)
rr=((xx-texpeakpix[1])**2+(yy-texpeakpix[0])**2)**0.5
mask=rr<r

yy2,xx2=np.indices(smooth_trot.shape)# These are all for the abundances, since the abundance peak in SgrB2S isn't at the temperature peak
rr2=((xx2-texpeakpix[1])**2+(yy2-texpeakpix[0])**2)**0.5
S_rr=((xx2-sgrb2scentralpix[1])**2+(yy2-sgrb2scentralpix[0])**2)**0.5# I may need to try using this on the abundance vs radius plot instead of the rr2
mask2=rr2<r
S_mask=S_rr<r

centrtopix=(rr[mask]*pixtophysicalsize).value
yinradius=rr[0]
xinradius=rr[1]
texinradius=texmapdata[rr<r].value
texerrinradius=texerrdata[rr<r].value
snrsinradius=snrs[rr<r]/5
abundinradius=abunds[rr2<r]
abundsnrinradius=snr_abund[rr2<r]
nh2inradius=nh2s[rr2<r].value
nh2errorinradius=nh2s_error[rr2<r].value
nh2snrinradius=(nh2inradius/nh2errorinradius)
lumsinradius=lums[rr2<r].value
lumerrinradius=lumserr[rr2<r].value
ntotsinradius=ntots[rr2<r].value
massinradius=h2mass[rr2<r].value#(nh2s[rr<r]*mu*beamarea_phys).to('solMass').value
masserrinradius=h2masserr[rr2<r].value#(nh2s_error[rr<r]*mu*beamarea_phys).to('solMass').value

S_massinradius=h2mass[S_rr<r].value
S_masserrinradius=h2masserr[S_rr<r].value

#pdb.set_trace()

lookformax=texmapdata[rr<10**0.5].value
lookformax_err=texerrdata[rr<10**0.5].value

abunderrinradius=(1/np.array(abundsnrinradius))*np.array(abundinradius)

teststack=np.stack((centrtopix,massinradius,snrsinradius,lumsinradius,nh2inradius,nh2errorinradius,texinradius,lumerrinradius,masserrinradius,abundinradius,abunderrinradius),axis=1)
teststack=teststack[teststack[:,0].argsort()]
set_centrtopix=set(teststack[:,0])
listordered_centrtopix=list(set_centrtopix)
listordered_centrtopix.sort()

avglist=[]
radialdensitylist=[]
err_radialdens=[]
mrcubedlist=[]
avgtexlist=[]
avgtexerrlist=[]
avginversesigma=[]
radialmaxtex=[]
radialmintex=[]
radialmasserr=[]
radialabun=[]
radialabunerr=[]

onlymassinteriorandradius=False
massinterior_radius=[]

totalmass=np.nansum(teststack[:,1])
sanitycheckmass=0
edge=None
localminmass=False


for bin in listordered_centrtopix:
    tempmass=[]
    tempsnr=[]
    tempmasssnr=[]
    temptex=[]
    temptexerr=[]
    tempinvsig=[]
    massesinterior=[]
    masserrinterior=[]
    tempmasserr=[]
    tempabun=[]
    tempabunerr=[]
    for data in teststack:
        #print(f'Bin: {bin} AU')
        if bin == data[0]:
            if np.isnan(data[1]):
                #print(f'Nan present in bin {data[0]}')
                continue
            else:
                tempmass.append(data[1])
                tempmasssnr.append(data[1]/data[8])
                tempmasserr.append(data[8])
                tempabun.append(data[9])
                tempabunerr.append(data[10])
                if np.isnan(data[6])==False:
                    if np.isinf(data[2])==False:
                        truesnr=data[2]*5
                        temptex.append(data[6])
                        temptexerr.append(data[6]/truesnr)
                        tempinvsig.append(truesnr/data[6])
                        tempsnr.append(truesnr)
                    
        if bin >= data[0]:
            massesinterior.append(data[1])
            masserrinterior.append(data[8])
        else:
            pass
    if len(tempmass)==0:
        avglist.append(0)
    else:
        avg=np.average(tempmass,weights=tempmasssnr)
        massinbin=np.nansum(tempmass)*u.solMass
        masserrinbin=np.sqrt(np.sum(np.square(tempmasserr)))*np.sqrt(pixperbeam)
        binvolume=(4/3)*np.pi*(bin*u.AU)**3

        massinteriorsum=np.nansum(massesinterior)*u.solMass
        masserrinteriorsum=np.sqrt(np.nansum(np.square(masserrinterior)))*u.solMass*np.sqrt(pixperbeam)
        rawdensityinbin=massinteriorsum/binvolume
        err_binrawdensity=np.sqrt(masserrinteriorsum/binvolume)**2#np.sqrt(binvolume**-1*masserrinbin*u.solMass)**2
        numberdensityinbin=(rawdensityinbin/mu).to('cm-3')
        err_binnumberdensity=(np.sqrt(mu**-1*err_binrawdensity)**2).to('cm-3')

        mrcubed=((massinbin/mu).to(''))/((bin*u.AU).to('cm'))**3
        
        if len(temptex)==0:
            print(f'Problem radius: {bin}')
            pdb.set_trace()
        else:
            pass
        avgtex=np.average(temptex,weights=tempsnr)
        #if np.isnan(avgtex):
        #    pdb.set_trace()
        #if len(avgtexterr)==0:
        #    avgtexerr=np.nan
        #else:
        avgtexerr=np.sqrt(np.sum(np.square(temptexerr)))#np.average(temptexerr)
        #pdb.set_trace()
        avginvsig=np.average(tempinvsig)
        
        if bin == 0:
            tempmaxtex=temptex[0]+temptexerr[0]
            tempmintex=temptex[0]-temptexerr[0]
            #pdb.set_trace()
        else:
            tempmaxtex=np.nanmax(temptex)
            tempmintex=np.nanmin(temptex)
        
        if avgtex<=150:#massinteriorsum >= 0.5*totalmass:
            if edge == None:
                #if source == 'DSii' or source == 'DSVI':
                #    sanitycheckmass=massinteriorsum
                #    localminmass=True
                #else:
                edge=bin
                sanitycheckmass=massinteriorsum
            else:
                pass
        else:
            pass
        
        ok=np.isfinite(tempabun)*np.isfinite(tempabunerr)
        if False in ok:
            tempabun=np.array(tempabun)[ok]
            tempabunerr=np.array(tempabunerr)[ok]
        if len(tempabun) == 0:
            pdb.set_trace()
            avgabuninrad=0
            err_avgabuninrad=np.nan
        avgabuninrad=np.average(tempabun, weights=(np.array(tempabun)/np.array(tempabunerr)))
        err_avgabuninrad=np.max(tempabun)-np.min(tempabun)#np.nanmean(tempabunerr)#This is the range in values, like the radial temperature profile
        
        if np.isfinite(avgabuninrad) == False:
            pdb.set_trace()

        avglist.append(avg)
        avgtexlist.append(avgtex)
        avgtexerrlist.append(avgtexerr)
        avginversesigma.append(avginvsig)
        radialmaxtex.append(tempmaxtex)
        radialmintex.append(tempmintex)
        radialdensitylist.append(numberdensityinbin.value)
        mrcubedlist.append(mrcubed.value)
        err_radialdens.append(err_binnumberdensity.value)
        radialabun.append(avgabuninrad)
        radialabunerr.append(err_avgabuninrad)
        
        massinterior_radius.append((massinteriorsum.value,masserrinteriorsum.value,bin))
    #pdb.set_trace()

avgabunerrpath=datadir+f'{source}_err_avgabun.txt'
np.savetxt(avgabunerrpath,np.array(radialabunerr))

if onlymassinteriorandradius:
    outlist=np.array(massinterior_radius)
    np.savetxt(datadir+f'{source}_massinterior_radius.txt',outlist)
    sys.exit()

if edge == None:#for hot sources where T never falls below 150 K
    mintemp=np.nanmin(avgtexlist)
    pix_mintemp=np.where(avgtexlist==mintemp)[0]
    edge=listordered_centrtopix[pix_mintemp[0]]

#pdb.set_trace()
fillwidth=np.copy(avgtexerrlist)#[x/2 for x in avgtexerrlist]
upperfill=radialmaxtex#list( map(add,avgtexlist,fillwidth))#radialmaxtex
lowerfill=radialmintex#list( map(sub,avgtexlist,fillwidth))#radialmintex#

massestosum=[]
masserrtoprop=[]
lumstosum=[]
lumerrtoprop=[]
nh2tomean=[]
nh2errortomean=[]
#pdb.set_trace()
if source == 'SgrB2S':#Selects masses,luminosities, nh2s, distances, and temperatures for use in radial plot
    
    rr2_sgrb2s=centrtopix#(rr2*pixtophysicalsize).value#list(pixmask.get_values(rr2*pixtophysicalsize).value)#rr2 and not rr here?
    '''
    xx_sgrb2s=[]
    yy_sgrb2s=[]
    for why in yinradius:
        for ex in xinradius:
            if regions.PixCoord(why,ex) in pixreg:
                xx_sgrb2s.append(ex)
                yy_sgrb2s.append(why)
            else:
                pass
    '''
    trotsforabunds=texinradius#list(pixmask.get_values(texmapdata.value))

onlytexabund=False

if onlytexabund:
    print('only data for X vs T will be outputted')
    if source == 'SgrB2S':
        np.savetxt(datadir+f'{source}_abuns.txt',abundinradius)
        np.savetxt(datadir+f'{source}_errabuns.txt',(np.array(abundinradius)/np.array(abundsnrinradius)))
        np.savetxt(datadir+f'{source}_tex.txt',trotsforabunds)
        np.savetxt(datadir+f'{source}_errtex.txt',(texerrinradius))
        print('Done')
        sys.exit()
    else:
        np.savetxt(datadir+f'{source}_abuns.txt',abundinradius)
        np.savetxt(datadir+f'{source}_errabuns.txt',(np.array(abundinradius)/np.array(abundsnrinradius)))
        np.savetxt(datadir+f'{source}_tex.txt',texinradius)
        np.savetxt(datadir+f'{source}_errtex.txt',(texerrinradius))
        print('Done')
        sys.exit()
        
for data2 in teststack:
    if data2[0] <= edge:
        massestosum.append(data2[1])
        masserrtoprop.append(data2[8])
        lumstosum.append(data2[3])
        lumerrtoprop.append(data2[7])
        nh2tomean.append(data2[4])
        nh2errortomean.append(data2[5])
    else:
        break
masssum=np.nansum(massestosum)
propmasserr=np.sqrt(np.nansum(np.square(masserrtoprop)))*np.sqrt(pixperbeam)#np.sum(masserrtoprop)
lumsum=np.nansum(lumstosum)
proplumerr=np.sqrt(np.nansum(np.square(lumerrtoprop)))*np.sqrt(pixperbeam)#np.sum(lumerrtoprop)
nh2mean=np.nanmean(nh2tomean)
nh2errormean=np.nanmean(nh2errortomean)

'''ds1 was 210'''
powerlaw_normpairs={'SgrB2S':(275,3500),'DSi':(210,3500),'DSii':(130,5000),'DSiii':(150,3000),'DSiv':(150,4500),'DSv':(160,2000),'DSVI':(130,4000),'DSVII':(200,2000),'DSVIII':(75,5000),'DSIX':(250,1600)}#7-75,4000
powerlawpair=powerlaw_normpairs[source]
fiducial_index=0.75
fiducial_norm=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**fiducial_index
lowerbound_initialguess={'SgrB2S':0.25,'DSi':0.25,'DSii':0.1,'DSiii':0.25,'DSiv':0.25,'DSv':0.25,'DSVI':0.25,'DSVII':0.25,'DSVIII':0.25,'DSIX':0.25}
bpl_alpha2_initialguess={'SgrB2S':0.25,'DSi':0.5,'DSii':0.5,'DSiii':0.25,'DSiv':0.5,'DSv':0.25,'DSVI':0.25,'DSVII':0.5,'DSVIII':0.25,'DSIX':0.6}
#ds9 was 0.25
if source == 'DSVII':
    inputamp=250#150
elif source == 'DSIX':
    inputamp=100
elif source == 'SgrB2S':
    inputamp=300
#elif source == 'DSv':
#    inputamp=250
elif source == 'DSIX':
    inputamp=250
else:
    inputamp=100
    
if source == 'DSIX':
    a1guess=0.01
    xbreakguess=1000
elif source == 'DSiv':
    a1guess=0.25
    xbreakguess=1000
elif source == 'DSVII':
    a1guess=-0.25
    xbreakguess=1000
else:
    a1guess=1
    xbreakguess=1000

fitter=fitting.LevMarLSQFitter()
fitter2=fitting.LevMarLSQFitter()
fitter3=fitting.LevMarLSQFitter()
base_bpl=powerlaws.BrokenPowerLaw1D(amplitude=inputamp,x_break=xbreakguess,alpha_1=a1guess,alpha_2=bpl_alpha2_initialguess[source])#original was -1 for all
base_spl=powerlaws.PowerLaw1D(amplitude=np.mean(radialdensitylist[1:]))
tex_spl=powerlaws.PowerLaw1D(amplitude=np.mean(avgtexlist))

if source=='DSIX':
    innerradius=9
    outerradius=len(listordered_centrtopix)-1#None

elif source=='DSiv':
    innerradius=1
    outerradius=None

elif source=='DSVII':
    innerradius=1
    outerradius=None

elif source == 'DSi':
    innerradius=1
    outerradius=119

elif source=='DSii':
    innerradius=1
    outerradius=len(radialdensitylist)-1
elif source == 'SgrB2S':
    innerradius=1
    outerradius=len(radialdensitylist)-4
else:
    innerradius=1
    outerradius=None

radiustofit=listordered_centrtopix[innerradius:outerradius]
textofit=avgtexlist[innerradius:outerradius]
texerrtofit=avgtexerrlist[innerradius:outerradius]
weightstofit=avginversesigma[innerradius:outerradius]
denstofit=radialdensitylist[innerradius:outerradius]#[innerradius:outerradius]
densweightstofit=1/(np.array(err_radialdens[innerradius:outerradius]))

#pdb.set_trace()
fit_spl=fitter3(tex_spl, radiustofit, textofit,weights=weightstofit)
#pdb.set_trace()
splperr=np.sqrt(np.diag(fitter3.fit_info['param_cov']))
#popt,pcov=cf(powerlaw_profile,radiustofit,textofit,sigma=texerrtofit,bounds=([0,lowerbound_initialguess[source]],[1000,1.25]))#[1:197] for ds4 (some nans may be present further out),[5:] for ds7, [9:] for ds9
fit_pl=fitter(base_bpl,radiustofit,textofit,weights=weightstofit)
perr=np.sqrt(np.diag(fitter.fit_info['param_cov']))

if source == 'DSiv'or 'DSv':
    fit_dens=fitter2(base_spl,radiustofit[:(len(radiustofit)-1)],denstofit[:(len(radiustofit)-1)],weights=densweightstofit[:(len(radiustofit)-1)])
#elif source == 'DSiv':

else:
    fit_dens=fitter2(base_spl,radiustofit,denstofit,weights=densweightstofit)
derr=np.sqrt(np.diag(fitter2.fit_info['param_cov']))

print(f'Sum: {masssum} +/- {propmasserr} Msun')
print(f'Core radius: {edge} +/- {cntmbmajtoAU/2}')
print(f'Core luminosity: {lumsum} +/- {proplumerr} Lsun')
print(f'Average H2 column in {edge} AU radius: {nh2mean} +/- {nh2errormean}')
plottexmax=np.nanmax(lookformax)+15
maxtexindex=np.where(lookformax==np.nanmax(lookformax))
maxtexerror=lookformax_err[maxtexindex]
print(f'Max Tex: {np.max(lookformax)} +/- {float(maxtexerror)} K')


copy_centrtopix=np.copy(listordered_centrtopix)
fittedtex=[]
    
plt.rcParams["figure.dpi"]=300

sourcenamesfortable={'SgrB2S':'SgrB2S','DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9'}

densalpha=fit_dens.alpha.value
errdensalpha=round_to_1(derr[2])
onlydens=False

plt.figure()
if source == 'DSiv':
    plt.errorbar(listordered_centrtopix,radialdensitylist,yerr=err_radialdens,label='Data',fmt='o')
    plt.plot(listordered_centrtopix,fit_dens(listordered_centrtopix),label=f'$p$={round(densalpha,(len(str(errdensalpha))-2))} \u00B1 {errdensalpha}',zorder=5)
    plt.axvline(edge,ls='--',label='$R_{core}$',color='black')
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('$n$ (cm$^{-3}$)',fontsize=14)
    plt.legend()
    plt.show()

else:
    plt.errorbar(listordered_centrtopix,radialdensitylist,yerr=err_radialdens,label='Data',fmt='o')
    plt.plot(listordered_centrtopix,fit_dens(listordered_centrtopix),label=f'$p$={round(densalpha,(len(str(errdensalpha))-2))} \u00B1 {errdensalpha}',zorder=5)
    plt.axvline(edge,ls='--',label='$R_{core}$',color='black')
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('$n$ (cm$^{-3}$)',fontsize=14)
    plt.legend()

    densplotpath=figpath+f'{dataversion}_densityprofile_trotntotbootmasked.png'
    print(f'\nSaving to {densplotpath}')
    plt.savefig(densplotpath)#,overwrite=True)
    plt.show()

densityslopepath=datadir+f'densityslopes_bootmasked.fits'

if os.path.exists(densityslopepath):
    dtable=QTable.read(densityslopepath)
    if densalpha in dtable['density_alpha']:
        print(f'alpha value {densalpha} already exists in table.')
    else:
        print(f'Appending alpha {densalpha} +/- {errdensalpha} to table.')
        densalpha=[densalpha,errdensalpha]
        preinsert=QTable(rows=[densalpha],names=['density_alpha','err_densityalpha'])
        densstack=vstack([dtable,preinsert])
        print('Appending complete.')
        #pdb.set_trace()
        print(f'Saving to {densityslopepath}')
        densstack.write(densityslopepath, overwrite=True)
        print('Done')
else:
    print('No density table found')
    print('Creating density table')
    densalpha=[densalpha,errdensalpha]
    dtable=QTable(rows=[densalpha],names=['density_alpha','err_densityalpha'])
    print('Table created')
    print(f'Writing to {densityslopepath}')
    dtable.write(densityslopepath)
    print('Done')

if onlydens:
    sys.exit()
else:
    pass

onlyradabun=False
if source == 'SgrB2S':
    plt.figure()
    plt.scatter(trotsforabunds,abundinradius,s=5,c=nh2inradius,norm=mpl.colors.LogNorm())
    plt.yscale('log')
    plt.xlabel('$T_{rot}$ (K)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.xlim(xmax=400)#,xmin=80)
    #plt.colorbar(pad=0,label='Luminosity (Lsun)')
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'N(CH$_3$OH) (cm$^{-2}$)')##
    figsavepath=figpath+f'{dataversion}_texabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    #pdb.set_trace()
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()

    plt.figure()
    plt.scatter(listordered_centrtopix,radialabun,s=5,)#c=nh2inradius,norm=mpl.colors.LogNorm())#abundinradius
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
    figsavepath=figpath+f'{dataversion}_real_radialavgabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    #pdb.set_trace()
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()
    
    savetxt=np.array([listordered_centrtopix,radialabun])
    np.savetxt(datadir+f'{source}_radialavgabun.txt',savetxt)

    if onlyradabun:
        sys.exit()

    plt.figure()
    plt.scatter(rr2_sgrb2s,nh2inradius,s=5,c=texinradius,vmax=plottexmax,cmap='inferno')
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('N(H$_2$) (cm$^{-2}$)',fontsize=14)
    plt.colorbar(pad=0,label='T$_{rot}$ (K)')#'T$_K$ (K)')
    figsavepath=figpath+f'{dataversion}_radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    #pdb.set_trace()
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()
else:
    plt.figure()
    plt.scatter(texinradius,abundinradius,s=5,c=nh2inradius,norm=mpl.colors.LogNorm())
    plt.yscale('log')
    plt.xlabel('$T_{rot}$ (K)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.xlim(xmax=plottexmax)
    #plt.colorbar(pad=0,label='Luminosity (Lsun)')
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'N(CH$_3$OH) (cm$^{-2}$)')##
    figsavepath=figpath+f'{dataversion}_texabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()

    plt.figure()
    plt.scatter(listordered_centrtopix,radialabun,s=5)#,c=nh2inradius,norm=mpl.colors.LogNorm())#abundinradius
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    #plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
    figsavepath=figpath+f'{dataversion}_real_radialavgabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()

    savetxt=np.array([listordered_centrtopix,radialabun])
    np.savetxt(datadir+f'{source}_radialavgabun.txt',savetxt)

    if onlyradabun:
        sys.exit()
    
    plt.figure()
    plt.scatter(centrtopix,nh2inradius,s=5,c=texinradius,vmax=plottexmax,cmap='inferno')
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('N(H$_2$) (cm$^{-2}$)',fontsize=14)
    plt.colorbar(pad=0,label='T$_{rot}$ (K)')#'T$_K$ (K)')
    figsavepath=figpath+f'{dataversion}_radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()





fig=plt.figure()#figsize=(14,7))

spec = gridspec.GridSpec(ncols=1, nrows=2, wspace=0.5,
                         hspace=.0, height_ratios=[2, 0.5])

ax0=fig.add_subplot(spec[0])

ax1=fig.add_subplot(spec[1],sharex=ax0)

residual_bpl=list( map(sub,avgtexlist,fit_pl(listordered_centrtopix)))
residual_spl=list( map(sub,avgtexlist,fit_spl(listordered_centrtopix)))#fittedtex))
avgsnr=np.array(list( map(truediv,avgtexlist,avgtexerrlist)))
avgsnrforplot=avgsnr/2

ax1.axhline(y=0,ls='--',color='black')
if source == 'DSiv':
    ax1.plot(copy_centrtopix[:197],residual_spl[:197],color='orange',label='single')#[:197] for ds4,[1:] for ds7
    ax1.plot(copy_centrtopix[:197],residual_bpl[:197],color='red',label='broken')#[:197] for ds4,
    ax0.scatter(listordered_centrtopix,avgtexlist)#[:197] for ds4
    ax0.fill_between(listordered_centrtopix[:197],upperfill[:197],lowerfill[:197],alpha=0.2,color='blue')#[:197] for ds

    ax0.plot(copy_centrtopix,fit_spl(copy_centrtopix),color='orange',label=f'$\\alpha$={round(fit_spl.alpha.value,2)} \u00B1 {round(splperr[2],2)}',zorder=1)
    ax0.plot(copy_centrtopix,fit_pl(copy_centrtopix),color='red',ls='-',zorder=2,label=f'$\\alpha_1$={round(fit_pl.alpha_1.value,2)} \u00B1 {round(perr[2],2)}\n$\\alpha_2$={round(fit_pl.alpha_2.value,2)} \u00B1 {round(perr[3],2)}\n$r_{{break}}$={round(fit_pl.x_break.value)} \u00B1 {round(perr[1])}')#[1:] for ds7
    ax1.set_xlabel('$r$ (AU)',fontsize=14)
    ax0.set_ylabel('T$_{rot}$ (K)',fontsize=14)
    ax1.set_ylabel('Residuals',fontsize=10)
    ax0.set_ylim(ymax=plottexmax+50)#(max(upperfill)+30))
    ax0.legend()
    ax0.tick_params(direction='in')
    ax0.tick_params(axis='x',labelcolor='w')
    ax1.tick_params(axis='x',top=True,direction='in')
    figsavepath=figpath+f'{dataversion}_radialavgtex_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'

    plt.tight_layout()

    plt.savefig(figsavepath)#,overwrite=True)

    plt.show()

elif source == 'DSVII':
    ax1.plot(copy_centrtopix[1:],residual_spl[1:],color='orange',label='single')#[:197] for ds4,[1:] for ds7
    ax1.plot(copy_centrtopix[1:],residual_bpl[1:],color='red',label='broken')#[:197] for ds4

    ax0.scatter(listordered_centrtopix,avgtexlist)#[:197] for ds4
    ax0.fill_between(listordered_centrtopix,upperfill,lowerfill,alpha=0.2,color='blue')#[:197] for ds4
    ax0.plot(copy_centrtopix,fit_spl(copy_centrtopix),color='orange',label=f'$\\alpha$={round(fit_spl.alpha.value,2)} \u00B1 {round(splperr[2],2)}',zorder=1)
    ax0.plot(copy_centrtopix[1:],fit_pl(copy_centrtopix)[1:],color='red',ls='-',zorder=2,label=f'$\\alpha_1$={round(fit_pl.alpha_1.value,2)} \u00B1 {round(perr[2],2)}\n$\\alpha_2$={round(fit_pl.alpha_2.value,2)} \u00B1 {round(perr[3],2)}\n$r_{{break}}$={round(fit_pl.x_break.value)} \u00B1 {round(perr[1])}')#[1:] for ds7
    ax1.set_xlabel('$r$ (AU)',fontsize=14)
    ax0.set_ylabel('T$_{rot}$ (K)',fontsize=14)
    ax1.set_ylabel('Residuals',fontsize=10)
    ax0.set_ylim(ymax=(max(upperfill)+30))
    ax0.legend()
    ax0.tick_params(direction='in')
    ax0.tick_params(axis='x',labelcolor='w')
    ax1.tick_params(axis='x',top=True,direction='in')
    figsavepath=figpath+f'{dataversion}_radialavgtex_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'

    plt.tight_layout()

    plt.savefig(figsavepath)#,overwrite=True)

    plt.show()
else:
    ax1.plot(copy_centrtopix,residual_spl,color='orange',label='single')#[:197] for ds4,[1:] for ds7
    ax1.plot(copy_centrtopix,residual_bpl,color='red',label='broken')#[:197] for ds4,

    ax0.scatter(listordered_centrtopix,avgtexlist)#[:197] for ds4
    ax0.fill_between(listordered_centrtopix,upperfill,lowerfill,alpha=0.2,color='blue')#[:197] for ds4
    ax0.plot(copy_centrtopix,fit_spl(copy_centrtopix),color='orange',label=f'$\\alpha$={round(fit_spl.alpha.value,2)} \u00B1 {round(splperr[2],2)}',zorder=1)
    ax0.plot(copy_centrtopix,fit_pl(copy_centrtopix),color='red',ls='-',zorder=2,label=f'$\\alpha_1$={round(fit_pl.alpha_1.value,2)} \u00B1 {round(perr[2],2)}\n$\\alpha_2$={round(fit_pl.alpha_2.value,2)} \u00B1 {round(perr[3],2)}\n$r_{{break}}$={round(fit_pl.x_break.value)} \u00B1 {round(perr[1])}')#[1:] for ds7
    ax1.set_xlabel('$r$ (AU)',fontsize=14)
    ax0.set_ylabel('T$_{rot}$ (K)',fontsize=14)
    ax1.set_ylabel('Residuals',fontsize=10)
    if source == 'SgrB2S':
        ax0.set_ylim(ymax=(575),ymin=100)#max(upperfill)+30),ymin=100)
    else:
        ax0.set_ylim(ymax=(max(upperfill)+30))
    ax0.legend()
    ax0.tick_params(direction='in')
    ax0.tick_params(axis='x',labelcolor='w')
    ax1.tick_params(axis='x',top=True,direction='in')
    figsavepath=figpath+f'{dataversion}_radialavgtex_LMLSQ_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'

    plt.tight_layout()

    plt.savefig(figsavepath)

    plt.show()
    
#sys.exit()    
if source=='SgrB2S':
    plt.figure()
    plt.scatter(nh2inradius,ntotsinradius,s=5,c=texinradius,vmax=plottexmax,cmap='inferno')
    '''
    x=[min(nh2inradius),max(nh2inradius)]
    y1=(9.5e-8*np.array(x))
    y2=(1.5e-7*np.array(x))
    y3=(5e-8*np.array(x))
    plt.plot(x,np.transpose([y1,y2,y3]))
    '''
    plt.yscale('log')
    plt.xscale('log')
    plt.colorbar(pad=0,label='T$_{rot}$ (K)')
    plt.xlabel('N(H$_2$) (cm$^{-2}$)')
    plt.ylabel('N(CH$_3$OH) (cm$^{-2}$)')
    figsavepath=figpath+f'{dataversion}_nch3ohvsnh2_smoothed.png'
    plt.savefig(figsavepath)#,overwrite=True)
    plt.show()

else:
    plt.figure()
    plt.scatter(nh2inradius,ntotsinradius,s=5,c=texinradius,cmap='inferno')
    '''
    x=[min(nh2inradius),max(nh2inradius)]
    y1=(9.5e-8*np.array(x))
    y2=(1.5e-7*np.array(x))
    y3=(5e-8*np.array(x))
    plt.plot(x,np.transpose([y1,y2,y3]))
    '''
    plt.yscale('log')
    plt.xscale('log')
    plt.colorbar(pad=0,label='T$_{rot}$ (K)')
    plt.xlabel('N(H$_2$) (cm$^{-2}$)')
    plt.ylabel('N(CH$_3$OH) (cm$^{-2}$)')
    figsavepath=figpath+f'{dataversion}_nch3ohvsnh2_smoothed.png'
    plt.savefig(figsavepath)
    plt.show()


onlypowerlaw=True
powerlawpath=datadir+'powerlawtable_bootmasked.fits'
pwrlwprams=[(round(fit_pl.alpha_1.value,2)*u.dimensionless_unscaled),(round(perr[2],2)*u.dimensionless_unscaled),(round(fit_pl.alpha_2.value,2)*u.dimensionless_unscaled),(round(perr[3],2)*u.dimensionless_unscaled),(round(fit_pl.x_break.value)*u.AU),(round(perr[1])*u.AU)]
powerlawparams=QTable(rows=[pwrlwprams],names=['alpha_1','alpha_1 error','alpha_2','alpha_2 error','x_break','x_break error'])

if os.path.exists(powerlawpath):
    print('\nStart power law table procedure')
    print(f'Power law table already exists at {powerlawpath}')
    pwrlwtable=QTable.read(powerlawpath)
    if pwrlwprams[4] in pwrlwtable['x_break']:
        print('Current parameter set already exists in table.')
        #print('Exiting...')
        #sys.exit()
    else:
        print('Appending new values')
        pwrlwstack=vstack([pwrlwtable,powerlawparams])
        print(f'Saving composite table to {powerlawpath}')
        pwrlwstack.write(powerlawpath,overwrite=True)
        print('Power law table update complete\n')
else:
    print('No power law table in current directory')
    if source != 'DSi':
        print('Please change source to DS1 to ensure proper stack')
        sys.exit()
    if source == 'SgrB2S':
        print('Nov 2023 pacman format in effect')
        print('Please make sure to change this, if it\'s not Nov 2023')
        powerlawparams.write(powerlawpath)
        print('Pacman powerlaw table created')
    else:
        print('Creating new power law table')
        powerlawparams.write(powerlawpath)
        print('New power law table created\n')

if onlypowerlaw:
    print('Only power law parameters requested.')
    print('Exiting...')
    sys.exit()
else:
    pass

sumtablepath=datadir+f'hotcoresummarytable.fits'
if os.path.exists(sumtablepath):
    print('\nUpdating summary table with new values')
    sumtable=QTable.read(sumtablepath)
    if sourcenamesfortable[source] in sumtable['Source']:
        print(f'Data for source {source} already exists.')
        sourceindex=int(np.where(sumtable['Source']==sourcenamesfortable[source])[0])
        print(sumtable[sourceindex])
        update=input('Overwrite? [y/n]\n')
        if update=='y':
              props=[sourcenamesfortable[source],(float(lookformax[maxtexindex])*u.K),(float(maxtexerror)*u.K),(nh2mean*u.cm**-2),(nh2errormean*u.cm**-2),
                     (masssum*u.solMass),(propmasserr*u.solMass),(lumsum*u.solLum),(proplumerr*u.solLum),(edge*u.AU),(cntmbmajtoAU/2)]
              sumtable[sourceindex]=tuple(props)
              print(f'{source} source properties updated.')
              print(f'\nSaving at {sumtablepath}')
              sumtable.write(sumtablepath,overwrite=True)
              print('Update complete')

        else:
              print(f'{source} source properties unchanged.')
    else:
        print(f'Adding new source {source} to summary table.')
        props=[sourcenamesfortable[source],(float(lookformax[maxtexindex])*u.K),(float(maxtexerror)*u.K),(nh2mean*u.cm**-2),(nh2errormean*u.cm**-2),
               (masssum*u.solMass),(propmasserr*u.solMass),(lumsum*u.solLum),(proplumerr*u.solLum),(edge*u.AU),(cntmbmajtoAU/2)]
        tempqtable=QTable(rows=[props],names=['Source','T_max','T_max_error','N(H_2) avg)','N(H_2) error','H_2 Mass','H_2 Mass error','Luminosity','Luminosity_error','Radius','Radius_error'])
        outtable=vstack([sumtable,tempqtable])
        print(f'{source} source properties added to table.')
        print(f'\nSaving at {sumtablepath}')
        outtable.write(sumtablepath,overwrite=True)
        print('Update complete')

else:
    newsumtable=input('Create new summary table? (y/n): ')
    if 'y' == newsumtable:
        print(f'Creating new summary table with {source} data')
        props=[sourcenamesfortable[source],(float(lookformax[maxtexindex])*u.K),(float(maxtexerror)*u.K),(nh2mean*u.cm**-2),(nh2errormean*u.cm**-2),
               (masssum*u.solMass),(propmasserr*u.solMass),(lumsum*u.solLum),(proplumerr*u.solLum),(edge*u.AU),(cntmbmajtoAU/2)]
        tempqtable=QTable(rows=[props],names=['Source','T_max','T_max_error','N(H_2) avg)','N(H_2) error','H_2 Mass','H_2 Mass error','Luminosity','Luminosity_error','Radius','Radius_error'])
        print(f'\nSaving at {sumtablepath}')
        outtable=tempqtable
        outtable.write(sumtablepath,overwrite=True)
        print('New summary table created')
    else:
        print(f'Source {source} results not saved')


'''
edgepoints=circle(texmapdata,texpeakpix[0],texpeakpix[1],r)

plt.imshow(texmapdata.value,origin='lower',vmax=550,vmin=10)
plt.scatter(edgepoints[0],edgepoints[1],color='orange')
#plt.scatter(xinradius,yinradius,color='orange')
plt.show()
'''
