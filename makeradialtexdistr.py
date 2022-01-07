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
    
source='DSi'
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'Kfield3originals_200K_trial1_noexclusions/','DSVIII':'Kfield3originals_175K_trial1_noexclusions/','DSIX':'Kfield7originals_150K_trial1_noexclusions/'}#base+'field10originals_z0_000186431_5-6mhzwidth_stdfixes/'
home=base+homedict[source]
fighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'
figpath=fighome+homedict[source]
if not os.path.exists(figpath):
    os.makedirs(figpath)
    print(f'Creating figpath {figpath}')
else:
    print(f'Figpath {figpath} already exists.')

notransmask=['DSv','DSVI','DSVII','DSVIII','DSIX']
if source in notransmask:
    texmap=home+"texmap_3sigma_allspw_withnans_weighted.fits"
else:
    texmap=home+"texmap_5transmask_3sigma_allspw_withnans_weighted.fits"

snrmap=home+"texmap_snr_allspw_weighted.fits"
abunmap=home+"ch3ohabundance_3sigma_ntotintercept_bolocamfeather.fits"
abunsnrmap=home+'ch3ohabundance_snr_ntotintercept_bolocamfeather.fits'
nh2map=home+"nh2map_3sigmacontandsurfacedensity_bolocamfeather.fits"
nh2errormap=home+"nh2map_error_bolocamfeather.fits"
lummap=home+"boltzmannlum_bolocamfeather.fits"
ntotmap=home+'ntotmap_allspw_withnans_weighted_useintercept_3sigma.fits'

texmap=fits.open(texmap)
texmapdata=texmap[0].data*u.K
snrs=fits.getdata(snrmap)
abunds=fits.getdata(abunmap)
snr_abund=fits.getdata(abunsnrmap)
nh2s=fits.getdata(nh2map)*u.cm**-2
nh2s_error=fits.getdata(nh2errormap)*u.cm**-2
lums=fits.getdata(lummap)*u.solLum
ntots=fits.getdata(ntotmap)*u.cm**-2

dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

cellsize=(np.abs(texmap[0].header['CDELT1']*u.deg)).to('arcsec')

pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')

print(pixtophysicalsize)

bmaj=texmap[0].header['BMAJ']*u.deg
bmajpix=round((bmaj/cellsize).value)
bmin=texmap[0].header['BMIN']*u.deg
bminpix=round((bmin/cellsize).value)
beamarea_sqdeg=bmaj*bmin
beamarea_sr=beamarea_sqdeg.to('sr')
bmajtophyssize=(np.tan(bmaj)*dGC).to('AU')
bmintophyssize=(np.tan(bmin)*dGC).to('AU')
'''Can probably simplify beamarea_phys to d(np.tan(bmaj)*np.tan(bmin))'''
beamarea_phys=np.pi*bmajtophyssize*bmintophyssize

pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}#y,x; DSiii was 24,24
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
rdict={'SgrB2S':12000*u.AU,'DSi':6400*u.AU,'DSii':9000*u.AU,'DSiii':6000*u.AU,'DSv':4000*u.AU,'DSIX':5000*u.AU}
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

centrtopix=(rr[mask]*pixtophysicalsize).value
yinradius=rr[0]
xinradius=rr[1]
texinradius=texmapdata[rr<r].value
snrsinradius=snrs[rr<r]/5
abundinradius=abunds[rr<r]
abundsnrinradius=snr_abund[rr<r]/5
nh2inradius=nh2s[rr<r].value
nh2errorinradius=nh2s_error[rr<r].value
lumsinradius=lums[rr<r].value
ntotsinradius=ntots[rr<r].value
massinradius=(nh2s[rr<r]*mu*beamarea_phys).to('solMass').value

lookformax=texmapdata[rr<10**0.5].value

teststack=np.stack((centrtopix,massinradius,snrsinradius,lumsinradius,nh2inradius,nh2errorinradius,texinradius),axis=1)
teststack=teststack[teststack[:,0].argsort()]
set_centrtopix=set(teststack[:,0])
listordered_centrtopix=list(set_centrtopix)
listordered_centrtopix.sort()

avglist=[]
avgtexlist=[]
avgtexerrlist=[]
avginversesigma=[]

for bin in listordered_centrtopix:
    tempmass=[]
    tempsnr=[]
    temptex=[]
    temptexerr=[]
    tempinvsig=[]
    for data in teststack:
        #print(f'Bin: {bin} AU')
        if bin == data[0]:
            if np.isnan(data[1]):
                #print(f'Nan present in bin {data[0]}')
                continue
            else:
                tempmass.append(data[1])
                tempsnr.append(data[2]*5)
                temptex.append(data[6])
                temptexerr.append(data[6]/data[2])
                tempinvsig.append(data[2]/data[6])
        else:
            pass
    if len(tempmass)==0:
        avglist.append(0)
    else:
        avg=np.average(tempmass,weights=tempsnr)
        avgtex=np.average(temptex,weights=tempsnr)
        avgtexerr=np.average(temptexerr)
        avginvsig=np.average(tempinvsig)
        avglist.append(avg)
        avgtexlist.append(avgtex)
        avgtexerrlist.append(avgtexerr)
        avginversesigma.append(avginvsig)
    #pdb.set_trace()

massderivative=list(np.diff(avglist))
edge=0
for diff in massderivative:
    if diff >= 0 and massderivative.index(diff) > 30:
        index=massderivative.index(diff)
        edge=listordered_centrtopix[index]
        break
        
fillwidth=np.copy(avgtexerrlist)#[x/2 for x in avgtexerrlist]
upperfill=list( map(add,avgtexlist,fillwidth))
lowerfill=list( map(sub,avgtexlist,fillwidth))

massestosum=[]
lumstosum=[]
nh2tomean=[]
nh2errortomean=[]
for data2 in teststack:
    if data2[0] <= edge:
        massestosum.append(data2[1])
        lumstosum.append(data2[3])
        nh2tomean.append(data2[4])
        nh2errortomean.append(data2[5])
    else:
        break
masssum=np.nansum(massestosum)
lumsum=np.nansum(lumstosum)
nh2mean=np.nanmean(nh2tomean)
nh2errormean=np.nanmean(nh2errortomean)

powerlaw_normpairs={'SgrB2S':(275,3500),'DSi':(210,3500),'DSii':(130,5000),'DSiii':(150,3000),'DSiv':(150,4500),'DSv':(160,2000),'DSVI':(130,4000),'DSVII':(75,4000),'DSVIII':(75,5000),'DSIX':(160,2500)}
powerlawpair=powerlaw_normpairs[source]
fiducial_index=0.75
fiducial_norm=powerlawpair[0]*(powerlawpair[1]/pixtophysicalsize)**fiducial_index
lowerbound_initialguess={'SgrB2S':0.25,'DSi':0.25,'DSii':0.1,'DSiii':0.25,'DSiv':0.25,'DSv':0.25,'DSVI':0.25,'DSVII':0.25,'DSVIII':0.25,'DSIX':0.25}
bpl_alpha2_initialguess={'SgrB2S':0.25,'DSi':0.5,'DSii':0.5,'DSiii':0.25,'DSiv':0.25,'DSv':0.25,'DSVI':0.25,'DSVII':0.25,'DSVIII':0.25,'DSIX':0.25}

fitter=fitting.LevMarLSQFitter()
base_bpl=powerlaws.BrokenPowerLaw1D(amplitude=1,x_break=1000,alpha_1=-1,alpha_2=bpl_alpha2_initialguess[source])#lowerbound_initialguess[source])#amp=150,xbreak=1500 for ds7

innerradius=1
radiustofit=listordered_centrtopix[innerradius:]
textofit=avgtexlist[innerradius:]
texerrtofit=avgtexerrlist[innerradius:]
weightstofit=avginversesigma[innerradius:]

popt,pcov=cf(powerlaw_profile,radiustofit,textofit,sigma=texerrtofit,bounds=([0,lowerbound_initialguess[source]],[1000,1.25]))#[1:197] for ds4 (some nans may be present further out),[5:] for ds7, [9:] for ds9
testindex=1
fit_pl=fitter(base_bpl,listordered_centrtopix[testindex:],avgtexlist[testindex:],weights=avginversesigma[testindex:])
perr=np.sqrt(np.diag(fitter.fit_info['param_cov']))

print(f'Sum: {masssum}')
print(f'Core radius: {edge} AU')
print(f'Core luminosity: {lumsum}')
print(f'Average H2 column in {edge} AU radius: {nh2mean} +/- {nh2errormean}')
plottexmax=np.nanmax(lookformax)+10
print(f'Max Tex: {np.max(lookformax)}')
#copy_centrtopix=np.copy(centrtopix)*u.AU
copy_centrtopix=np.copy(listordered_centrtopix)#np.linspace(0,r_phys.value,num=len(avgtexlist))*u.AU
#copy_centrtopix.sort()
lineartex=[]
quadrtex=[]
sublinhalftex=[]
sublinmidtex=[]
fittedtex=[]

for dist in copy_centrtopix:
    lineartex.append(linear_profile(dist).value)
    quadrtex.append(quadratic_profile(dist).value)
    sublinhalftex.append(sublinear_profile(dist).value)
    sublinmidtex.append(submid_profile(dist).value)
    fittedtex.append((popt[0]*((dist/pixtophysicalsize.value)**-popt[1])))#.value)
    
plt.rcParams["figure.dpi"]=150
'''
ax=plt.subplot(111)

vmaxdict={'DSi':1e-5}
plt.scatter(centrtopix,texinradius,s=snrsinradius,c=abundinradius,cmap='Greens',alpha=0.7)#,norm=mpl.colors.LogNorm())
plt.plot(copy_centrtopix,lineartex,color='yellow',label=r'r$^{-1}$')
plt.plot(copy_centrtopix,quadrtex,color='green',label=r'r$^{-2}$')
plt.plot(copy_centrtopix,sublinhalftex,color='purple',label=r'r$^{-0.5}$')
ax.set_xlabel('$d$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
ax.tick_params(size=14)
plt.tight_layout()
plt.colorbar(pad=0)
plt.legend()

savefigpath=home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_linear_quad_sqrt.png'
#plt.savefig(savefigpath,overwrite=True)
plt.show()

plt.close()
'''
ax=plt.subplot(111)
plt.scatter(centrtopix,texinradius,s=snrsinradius,c=abundinradius,norm=mpl.colors.LogNorm(),cmap='viridis')#,alpha=0.7)

plt.plot(copy_centrtopix,sublinhalftex,color='red',label=r'r$^{-0.5}$')
plt.plot(copy_centrtopix,sublinmidtex,color='pink',label=r'r$^{-0.75}$')
plt.plot(copy_centrtopix,lineartex,color='orange',label=r'r$^{-1}$')
plt.plot(copy_centrtopix,quadrtex,color='yellow',label=r'r$^{-1.25}$')



ax.set_xlabel('$r$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
plt.ylim(ymax=plottexmax,ymin=(np.nanmin(texinradius)-10))
ax.tick_params(size=14)
plt.tight_layout()
plt.colorbar(pad=0,label='X(CH$_3$OH)')
plt.legend(loc=1)

savefigpath=home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_interceptabundances.png'#home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_lognorm_linear_quad_sqrt_interceptabundances.png'
figsavepath=figpath+f'radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_lognorm_linear_quad_sqrt_interceptabundances_bolocamfeather.png'
#plt.savefig(savefigpath,overwrite=True)
plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)

plt.show()

plt.figure()
plt.plot(listordered_centrtopix,avglist)
plt.axvline(edge,ls='--')
plt.show()

plt.figure()
plt.scatter(texinradius,abundinradius,s=abundsnrinradius,c=nh2inradius,norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xlabel('$T_K$ (K)',fontsize=14)
plt.ylabel('X(CH$_3$OH)',fontsize=14)
plt.xlim(xmax=plottexmax)
#plt.colorbar(pad=0,label='Luminosity (Lsun)')
plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'N(CH$_3$OH) (cm$^{-2}$)')##
figsavepath=figpath+f'texabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_interceptabundances_bolocamfeather.png'
plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.scatter(centrtopix,abundinradius,s=snrsinradius,c=nh2inradius,norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xlabel('$r$ (AU)',fontsize=14)
plt.ylabel('X(CH$_3$OH)',fontsize=14)
plt.colorbar(pad=0,label='N(H$_2$) (cm$^{-2}$)')#'T$_K$ (K)')
figsavepath=figpath+f'radialavgabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_interceptabundances_bolocamfeather.png'
plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.scatter(centrtopix,nh2inradius,s=snrsinradius,c=texinradius,vmax=plottexmax,cmap='inferno')
plt.yscale('log')
plt.xlabel('$r$ (AU)',fontsize=14)
plt.ylabel('N(H$_2$) (cm$^{-2}$)',fontsize=14)
plt.colorbar(pad=0,label='T$_K$ (K)')#'T$_K$ (K)')
figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_bolocamfeather.png'
plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

fig=plt.figure()#figsize=(14,7))

spec = gridspec.GridSpec(ncols=1, nrows=2, wspace=0.5,
                         hspace=.0, height_ratios=[2, 0.5])

ax0=fig.add_subplot(spec[0])

ax1=fig.add_subplot(spec[1],sharex=ax0)

residual_bpl=list( map(sub,avgtexlist,fit_pl(listordered_centrtopix)))
residual_spl=list( map(sub,avgtexlist,fittedtex))

ax1.axhline(y=0,ls='--',color='black')
ax1.plot(copy_centrtopix,residual_spl,color='orange',label='single')
ax1.plot(copy_centrtopix,residual_bpl,color='red',label='broken')
#ax1.legend()

ax0.scatter(listordered_centrtopix,avgtexlist)#[:197] for ds4
ax0.fill_between(listordered_centrtopix,upperfill,lowerfill,alpha=0.2,color='blue')
#plt.errorbar(listordered_centrtopix,avgtexlist,yerr=avgtexerrlist,fmt='o')#[:(len(listordered_centrtopix)-2)] for ds4
ax0.plot(copy_centrtopix,fittedtex,color='orange',label=f'$\\alpha$={round(popt[1],2)} \u00B1 {round(pcov[1,1]**0.5,2)}',zorder=1)
ax0.plot(copy_centrtopix,fit_pl(copy_centrtopix),color='red',ls='-',zorder=2,label=f'$\\alpha_1$={round(fit_pl.alpha_1.value,2)} \u00B1 {round(perr[2],2)}\n$\\alpha_2$={round(fit_pl.alpha_2.value,2)} \u00B1 {round(perr[3],2)}\n$r_{{break}}$={round(fit_pl.x_break.value)} \u00B1 {round(perr[1])}')
ax1.set_xlabel('$r$ (AU)',fontsize=14)
ax0.set_ylabel('T$_K$ (K)',fontsize=14)
ax1.set_ylabel('Residuals',fontsize=10)
ax0.set_ylim(ymax=(max(avgtexlist)+30))
ax0.legend()
ax0.tick_params(direction='in')
ax1.tick_params(axis='x',top=True,direction='in')
#plt.setp(ax0.get_xticklabels(), visible=False)
#plt.setp(ax0.get_xticks(),visible=True)

plt.tight_layout()

plt.show()
print(f'Norm initial guess: {fiducial_norm}')
print(f'norm error: {pcov[0,0]**0.5}')

print(f'index initial guess: {lowerbound_initialguess[source]}')
print(f'index error: {pcov[1,1]**0.5}')


'''
plt.scatter(nh2inradius,ntotsinradius,s=snrsinradius,c=texinradius)
plt.yscale('log')
plt.xscale('log')
plt.colorbar(pad=0,label='T$_K$ (K)')
plt.xlabel('N(H$_2$) (cm$^{-2}$)')
plt.ylabel('N(CH$_3$OH) (cm$^{-2}$)')
plt.show()
'''
'''
edgepoints=circle(texmapdata,texpeakpix[0],texpeakpix[1],r)

plt.imshow(texmapdata.value,origin='lower',vmax=550,vmin=10)
plt.scatter(edgepoints[0],edgepoints[1],color='orange')
#plt.scatter(xinradius,yinradius,color='orange')
plt.show()
'''
