from astropy.io import fits
import numpy as np
import astropy.units as u
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import visualization
import os
import pdb

mpl.interactive(False)

plt.clf()

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

def subquarter_profile(r,peak):
    return peak*((r/pixtophysicalsize)**-0.25)
    
def sublinear_profile(r,peak):
    return peak*((r/pixtophysicalsize)**-0.5)

def linear_profile(r,peak):
    return peak*((r/pixtophysicalsize)**-1)
    
def quadratic_profile(r,peak):
    return peak*((r/pixtophysicalsize)**-2)
    
source='DSVIII'
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3}
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'Kfield3originals_200K_trial1_noexclusions/','DSVIII':'Kfield3originals_175K_trial1_noexclusions/'}#base+'field10originals_z0_000186431_5-6mhzwidth_stdfixes/'
home=base+homedict[source]
fighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'
figpath=fighome+homedict[source]
if not os.path.exists(figpath):
    os.makedirs(figpath)
    print(f'Creating figpath {figpath}')
else:
    print(f'Figpath {figpath} already exists.')
'''
source='SgrB2S'
fnum=1
home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/"
#home='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/OctReimage_z0_000186431_5-6mhzwidth_stdfixes/'
#home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/"
'''
notransmask=['DSv','DSVI','DSVII','DSVIII']
if source in notransmask:
    texmap=home+"texmap_3sigma_allspw_withnans_weighted.fits"
else:
    texmap=home+"texmap_5transmask_3sigma_allspw_withnans_weighted.fits"

snrmap=home+"texmap_snr_allspw_weighted.fits"
abunmap=home+"ch3ohabundance_3sigma_ntotintercept.fits"
abunsnrmap=home+'ch3ohabundance_snr_ntotintercept.fits'
nh2map=home+"nh2map_3sigmacontandsurfacedensity.fits"
lummap=home+"boltzmannlum.fits"
ntotmap=home+'ntotmap_allspw_withnans_weighted_useintercept_3sigma.fits'

texmap=fits.open(texmap)
texmapdata=texmap[0].data*u.K
snrs=fits.getdata(snrmap)
abunds=fits.getdata(abunmap)
snr_abund=fits.getdata(abunsnrmap)
nh2s=fits.getdata(nh2map)*u.cm**-2
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

pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50)}#y,x; DSiii was 24,24
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
r_phys=12000*u.AU#r*pixtophysicalsize.to('pc')
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

for y in range(np.shape(texmapdata)[0]):
    for x in range(np.shape(texmapdata)[1]):
        if (y-texpeakpix[0])**2+(x-texpeakpix[1])**2 <= r**2:
            texinradius.append(texmapdata[y,x].value)
            xinradius.append(x)
            yinradius.append(y)
            centrtopix.append((np.sqrt((y-texpeakpix[0])**2+(x-texpeakpix[1])**2)*pixtophysicalsize).value)
            snrsinradius.append(snrs[y,x]/5)#Scaling down by 5 to better be able to distinguish datapoints in plot
            abundinradius.append(abunds[y,x])
            abundsnrinradius.append(snr_abund[y,x]/5)#ditto as tex snrs
            nh2inradius.append(nh2s[y,x].value)
            massinradius.append((nh2s[y,x]*mu*beamarea_phys).to('solMass').value)#this is actually real mass, not column density
            lumsinradius.append(lums[y,x].value)
            ntotsinradius.append(ntots[y,x].value)
        '''
        if (y-nh2peakpix[0])**2+(x-nh2peakpix[1])**2 <= r**2:
            abundinradius.append(abunds[y,x])
            nh2inradius.append((nh2s[y,x]*mu*beamarea_phys).to('solMass').value)#this is actually real mass, not column density
            #pdb.set_trace()
        '''
        if (y-texpeakpix[0])**2+(x-texpeakpix[1])**2 <= 10:
            lookformax.append(texmapdata[y,x].value)
        else:
            pass
            
teststack=np.stack((centrtopix,massinradius,snrsinradius,lumsinradius),axis=1)
teststack=teststack[teststack[:,0].argsort()]
set_centrtopix=set(teststack[:,0])
listordered_centrtopix=list(set_centrtopix)
listordered_centrtopix.sort()

avglist=[]

for bin in listordered_centrtopix:
    tempmass=[]
    tempsnr=[]
    for data in teststack:
        #print(f'Bin: {bin} AU')
        if bin == data[0]:
            if np.isnan(data[1]):
                continue
            else:
                tempmass.append(data[1])
                tempsnr.append(data[2]*5)
        else:
            pass
    if len(tempmass)==0:
        avglist.append(0)
    else:
        avg=np.average(tempmass,weights=tempsnr)
        avglist.append(avg)
    #pdb.set_trace()

massderivative=list(np.diff(avglist))
edge=0
for diff in massderivative:
    if diff >= 0 and massderivative.index(diff) > 30:
        index=massderivative.index(diff)
        edge=listordered_centrtopix[index]
        break
        
massestosum=[]
lumstosum=[]
m2sum2=avglist[:index]
for data2 in teststack:
    if data2[0] <= edge:
        massestosum.append(data2[1])
        lumstosum.append(data2[3])
    else:
        break
masssum=np.nansum(massestosum)
lumsum=np.nansum(lumstosum)
msum=np.nansum(m2sum2)
print(f'Sum: {masssum}')
print(f'Sum2: {msum}')
print(f'Core radius: {edge}')
print(f'Core luminosity: {lumsum}')
plottexmax=np.max(lookformax)+10
print(f'Max Tex: {np.max(lookformax)}')
copy_centrtopix=np.copy(centrtopix)*u.AU
#copy_centrtopix.sort()
lineartex=[]
quadrtex=[]
sublinhalftex=[]
sublinquarttex=[]

for dist in copy_centrtopix:
    lineartex.append(linear_profile(dist,texmapdata[texpeakpix[0],texpeakpix[1]]).value)
    quadrtex.append(quadratic_profile(dist,texmapdata[texpeakpix[0],texpeakpix[1]]).value)
    sublinhalftex.append(sublinear_profile(dist,texmapdata[texpeakpix[0],texpeakpix[1]]).value)
    sublinquarttex.append(subquarter_profile(dist,texmapdata[texpeakpix[0],texpeakpix[1]]).value)
    
plt.rcParams["figure.dpi"]=150

ax=plt.subplot(111)
'''
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
plt.savefig(savefigpath,overwrite=True)
plt.show()

plt.close()
'''
#ax=plt.subplot(111)
#plt.scatter(centrtopix,texinradius,s=snrsinradius,c=abundinradius,norm=mpl.colors.LogNorm(vmin=1e-8),cmap='viridis',alpha=0.7)
'''
plt.plot(copy_centrtopix,lineartex,color='yellow',label=r'r$^{-1}$')
plt.plot(copy_centrtopix,quadrtex,color='green',label=r'r$^{-2}$')
plt.plot(copy_centrtopix,sublinhalftex,color='purple',label=r'r$^{-0.5}$')
plt.plot(copy_centrtopix,sublinquarttex,color='orange',label=r'r$^{-0.25}$')
'''
'''
ax.set_xlabel('$d$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
plt.ylim(ymax=plottexmax)
ax.tick_params(size=14)
plt.tight_layout()
plt.colorbar(pad=0,label='X(CH$_3$OH)')
#plt.legend(loc=1)

savefigpath=home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_interceptabundances.png'#home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_lognorm_linear_quad_sqrt_interceptabundances.png'
figsavepath=figpath+f'radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_lognorm_linear_quad_sqrt_interceptabundances.png'
#plt.savefig(savefigpath,overwrite=True)
#plt.savefig(figsavepath,overwrite=True)

plt.show()

plt.plot(listordered_centrtopix,avglist)
plt.axvline(edge,ls='--')
plt.show()

plt.scatter(texinradius,abundinradius,s=abundsnrinradius,c=ntotsinradius,norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xlabel('$T_K$ (K)',fontsize=14)
plt.ylabel('X(CH$_3$OH)',fontsize=14)
plt.xlim(xmax=plottexmax)
#plt.colorbar(pad=0,label='Luminosity (Lsun)')
plt.colorbar(pad=0,label='N(CH$_3$OH) (cm$^{-2}$)')#'N(H2) (cm-2)')#
figsavepath=figpath+f'radialabundiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_interceptabundances.png'
#plt.savefig(figsavepath,overwrite=True)
plt.show()
'''
plt.scatter(nh2inradius,ntotsinradius,s=snrsinradius,c=texinradius)
plt.yscale('log')
plt.xscale('log')
plt.colorbar(pad=0,label='T$_K$ (K)')
plt.xlabel('N(H$_2$) (cm$^{-2}$)')
plt.ylabel('N(CH$_3$OH) (cm$^{-2}$)')
plt.show()

'''
edgepoints=circle(texmapdata,texpeakpix[0],texpeakpix[1],r)

plt.imshow(texmapdata.value,origin='lower',vmax=550,vmin=10)
plt.scatter(edgepoints[0],edgepoints[1],color='orange')
#plt.scatter(xinradius,yinradius,color='orange')
plt.show()
'''

