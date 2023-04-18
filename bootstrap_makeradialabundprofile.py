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
import astropy.stats

mpl.interactive(True)

plt.close('all')

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
    
source='SgrB2S'#os.getenv('SOURCE')#
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':'/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/','DSi':'/nov2022continuumsanitycheck/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}#{'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'Kfield3originals_200K_trial1_noexclusions/','DSVIII':'Kfield3originals_175K_trial1_noexclusions/','DSIX':'Kfield7originals_150K_trial1_noexclusions/'}#base+'field10originals_z0_000186431_5-6mhzwidth_stdfixes/'
home=base+homedict[source]
fighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'
figpath=fighome+homedict[source]
if not os.path.exists(figpath):
    os.makedirs(figpath)
    print(f'Creating figpath {figpath}')
else:
    print(f'Figpath {figpath} already exists.')

#notransmask=['DSv','DSVI','DSVII','DSVIII','DSIX']
#if source == 'SgrB2S':
texmap=home+"bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
#else:
#texmap=home+"texmap_5transmask_3sigma_allspw_withnans_weighted.fits"

texerrmap=home+'error_trot_boostrap1000_nonegativeslope.fits'#"texmap_error_allspw_withnans_weighted.fits"
#snrmap=home+"texmap_snr_allspw_weighted.fits"
abunmap=home+'bootstrap_ch3ohabundance_3sigma_ntotintercept_intstd_bolocamfeather_smoothedtobolocam.fits'#'bootstrap_ch3ohabundance_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
ntotnh2_abunerrpath=home+'bootstrap_ch3ohabundance_error_ntotnh2mask_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
abunsnrmap=home+'bootstrap_ch3ohabundance_snr_intstd_ntotintercept_bolocamfeather_smoothedtobolocam.fits'#'bootstrap_ch3ohabundance_snr_ntotintercept_bolocamfeather_smoothedtobolocam.fits'
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


texmap=fits.open(texmap)
texmapdata=texmap[0].data*u.K
texerrdata=np.squeeze(fits.getdata(texerrmap))*u.K
snrs=(texmapdata/texerrdata).value#fits.getdata(snrmap)
abunds=np.squeeze(fits.getdata(abunmap))
stopgap_errabun=np.squeeze(fits.getdata(ntotnh2_abunerrpath))
snr_abund=abunds/stopgap_errabun#fits.getdata(abunsnrmap)
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

if source == 'SgrB2S':
    wcsobj=WCS(smooth_trotfits[0].header)

    regs = regions.Regions.read('/blue/adamginsburg/d.jeff/imaging_results/regfiles/roughsgrb2smassregion_ignoresHIIregion.reg')
    pixreg = regs[0].to_pixel(wcsobj)
    pixmask = pixreg.to_mask()
    
    texmapdata=pixmask.cutout(texmapdata,fill_value=np.nan)
    '''
    texmask=np.ma.masked_where(texmapdata > 1000, texmapdata)
    texmapdata=texmask.filled()
    '''
    texerrdata=pixmask.cutout(texerrdata,fill_value=np.nan)
    snrs=np.squeeze(texmapdata/texerrdata)
    abunds=pixmask.cutout(abunds,fill_value=np.nan)
    snr_abund=pixmask.cutout(snr_abund,fill_value=np.nan)
    nh2s=pixmask.cutout(nh2s,fill_value=np.nan)
    nh2s_error=pixmask.cutout(nh2s_error,fill_value=np.nan)
    lums=pixmask.cutout(lums,fill_value=np.nan)
    lumserr=pixmask.cutout(lumserr,fill_value=np.nan)
    ntots=pixmask.cutout(ntots,fill_value=np.nan)
    ntoterr=pixmask.cutout(ntoterr,fill_value=np.nan)
    h2mass=pixmask.cutout(h2mass,fill_value=np.nan)
    h2masserr=pixmask.cutout(h2masserr,fill_value=np.nan)
    smooth_trot=pixmask.cutout(smooth_trot,fill_value=np.nan)
    smooth_trot_err=pixmask.cutout(smooth_trot_err,fill_value=np.nan)
    

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

pixdict={'SgrB2S':(26,14),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}#y,x; DSiii was 24,24;S was 73,54 - contfix S was 69,58
sgrb2scentralpix=(25,25)#contfix pix, prefinal was (66,70)
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
rdict={'SgrB2S':10000*u.AU,'DSi':7500*u.AU,'DSii':8700*u.AU,'DSiii':6000*u.AU,'DSiv':8000*u.AU,'DSv':3500*u.AU,'DSVII':6000*u.AU,'DSVIII':5700*u.AU,'DSIX':5000*u.AU}#1-6400,4-8500,5-4000,7-6600,S-12000
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
trad=180

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
        avgabuninrad=np.average(tempabun, weights=(np.array(tempabun)/np.array(tempabunerr)))
        if len(tempabun) < 2:
            err_avgabuninrad=tempabunerr[0]
        else:
            err_avgabuninrad=astropy.stats.mad_std(tempabun)#np.max(tempabun)-np.min(tempabun)#np.nanmean(tempabunerr)#This is the range in values, like the radial temperature profile
        
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

avgabunerrpath=f'{source}_err_intstd_avgabun.txt'
np.savetxt(avgabunerrpath,np.array(radialabunerr))

if source == 'SgrB2S':
    rr2_sgrb2s=centrtopix
    trotsforabunds=texinradius

plottexmax=np.nanmax(lookformax)+10
maxtexindex=np.where(lookformax==np.nanmax(lookformax))
maxtexerror=lookformax_err[maxtexindex]

plt.rcParams['figure.dpi']=150

if source == 'SgrB2S':
    plt.figure(figsize=(7,5))
    plt.scatter(trotsforabunds,abundinradius,s=5,c=nh2inradius,norm=mpl.colors.LogNorm())
    plt.yscale('log')
    plt.xlabel('$T_{rot}$ (K)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.xlim(xmax=plottexmax,xmin=80)
    #plt.colorbar(pad=0,label='Luminosity (Lsun)')
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'N(CH$_3$OH) (cm$^{-2}$)')##
    figsavepath=figpath+f'texabundiag_contsanitycheck_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    #pdb.set_trace()
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()

    plt.figure()
    plt.scatter(listordered_centrtopix,radialabun,s=5,)#c=nh2inradius,norm=mpl.colors.LogNorm())#abundinradius
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
    figsavepath=figpath+f'real_radialavgabundiag_contsanitycheck_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    #pdb.set_trace()
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()
    '''
    savetxt=np.array([listordered_centrtopix,radialabun])
    np.savetxt(f'{source}_intstd_radialavgabun.txt',savetxt)
    '''
else:
    plt.figure(figsize=(7,5))
    plt.scatter(texinradius,abundinradius,s=5,c=nh2inradius,norm=mpl.colors.LogNorm())
    plt.yscale('log')
    plt.xlabel('$T_{rot}$ (K)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    plt.xlim(xmax=plottexmax)
    #plt.colorbar(pad=0,label='Luminosity (Lsun)')
    plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'N(CH$_3$OH) (cm$^{-2}$)')##
    figsavepath=figpath+f'texabundiag_contsanitycheck_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()

    plt.figure()
    plt.scatter(listordered_centrtopix,radialabun,s=5)#,c=nh2inradius,norm=mpl.colors.LogNorm())#abundinradius
    plt.yscale('log')
    plt.xlabel('$r$ (AU)',fontsize=14)
    plt.ylabel('X(CH$_3$OH)',fontsize=14)
    #plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
    figsavepath=figpath+f'real_radialavgabundiag_contsanitycheck_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
    plt.savefig(figsavepath,bbox_inches='tight',)#overwrite=True)
    plt.show()
    '''
    savetxt=np.array([listordered_centrtopix,radialabun])
    np.savetxt(f'{source}_intstd_radialavgabun.txt',savetxt)
    '''