from astropy.io import fits
import numpy as np
import astropy.units as u
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import visualization

mpl.interactive(False)

plt.clf()

def circle(data,ycenter,xcenter,rad):
    edgex=[]
    edgey=[]
    for i in range(np.shape(data)[0]):
        for j in range(np.shape(data)[1]):
            if (j-xcenter)**2+(i-ycenter)**2==rad**2:
                edgex.append(j)
                edgey.append(i)
    return np.vstack((edgex,edgey))
    
    
source='DSv'
fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10}
fnum=fielddict[source]
print(f'Source: {source}')
base=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'
homedict={'SgrB2S':"new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"Kfield10originals_noexclusions/",'DSiii':"Kfield10originals_noexclusions/",'DSiv':"Kfield10originals_noexclusions/",'DSv':"Kfield10originals_noexclusions_include4-3_150K_trial2/"}#base+'field10originals_z0_000186431_5-6mhzwidth_stdfixes/'
home=base+homedict[source]
'''
source='SgrB2S'
fnum=1
home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/"
#home='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/OctReimage_z0_000186431_5-6mhzwidth_stdfixes/'
#home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/"
'''
if source == 'DSv':
    texmap=home+"texmap_0transmask_3sigma_allspw_withnans_weighted.fits"
else:
    texmap=home+"texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
snrmap=home+"texmap_snr_allspw_weighted.fits"
abunmap=home+"ch3ohabundance.fits"

texmap=fits.open(texmap)
texmapdata=texmap[0].data*u.K
snrs=fits.getdata(snrmap)
abunds=fits.getdata(abunmap)

dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

cellsize=(np.abs(texmap[0].header['CDELT1']*u.deg)).to('arcsec')

pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')

print(pixtophysicalsize)

pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19)}#y,x

texpeakpix=pixdict[source]
#(36,43)#DSi hotspot
#texpeakpix=(73,56)#SgrB2S hotspot
#x-1, y-1 from DS9

print(f'Center p: {texmapdata[texpeakpix[0],texpeakpix[1]]}')

r=35
#pixradius=math.ceil((0.08*u.pc/pixtophysicalsize).to(''))
r_phys=r*pixtophysicalsize.to('pc')
print(f'physical radius: {r_phys}')

xpixs=np.arange(texpeakpix[0],(texpeakpix[0]+r))

texinradius=[]
xinradius=[]
yinradius=[]
centrtopix=[]
snrsinradius=[]
abundinradius=[]

for y in range(np.shape(texmapdata)[0]):
    for x in range(np.shape(texmapdata)[1]):
        if (y-texpeakpix[0])**2+(x-texpeakpix[1])**2 <= r**2:
            texinradius.append(texmapdata[y,x].value)
            xinradius.append(x)
            yinradius.append(y)
            centrtopix.append((np.sqrt((y-texpeakpix[0])**2+(x-texpeakpix[1])**2)*pixtophysicalsize).value)
            snrsinradius.append(snrs[y,x]/5)#Scaling down by 5 to better be able to distinguish datapoints in plot
            abundinradius.append(abunds[y,x])
        else:
            pass
            
plt.rcParams["figure.dpi"]=150

ax=plt.subplot(111)

vmaxdict={'DSi':1e-5}
plt.scatter(centrtopix,texinradius,s=snrsinradius,c=abundinradius,cmap='rainbow')#,norm=mpl.colors.LogNorm())
ax.set_xlabel('$d$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
ax.tick_params(size=14)
plt.tight_layout()
plt.colorbar(pad=0)

savefigpath=home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU.png'
plt.savefig(savefigpath,overwrite=True)
plt.show()

plt.close()

ax=plt.subplot(111)
plt.scatter(centrtopix,texinradius,s=snrsinradius,c=abundinradius,norm=mpl.colors.LogNorm(),cmap='rainbow')
ax.set_xlabel('$d$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
ax.tick_params(size=14)
plt.tight_layout()
plt.colorbar(pad=0)

savefigpath=home+f'figures/radialtexdiag_r{r}px_rphys{int(pixtophysicalsize.value)}AU_lognorm.png'
plt.savefig(savefigpath,overwrite=True)
plt.show()


'''
edgepoints=circle(texmapdata,texpeakpix[0],texpeakpix[1],r)

plt.imshow(texmapdata.value,origin='lower',vmax=550,vmin=10)
plt.scatter(edgepoints[0],edgepoints[1],color='orange')
#plt.scatter(xinradius,yinradius,color='orange')
plt.show()
'''

