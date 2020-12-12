from astropy.io import fits
import numpy as np
import astropy.units as u
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.interactive(True)

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
    
home='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/OctReimage_z0_000186431_5-6mhzwidth_stdfixes/'
#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/"
texmap=home+"texmap_3sigma_allspw_withnans_weighted.fits"

texmapdata=fits.getdata(texmap)*u.K

dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

cellsize=0.05*u.arcsec

pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')

print(pixtophysicalsize)

texpeakpix=(36,41)
#(73,56)SgrB2S hotspot
#x-1, y-1 from DS9

print(f'max: {texmapdata[texpeakpix[0],texpeakpix[1]]}')

r=12
#pixradius=math.ceil((0.08*u.pc/pixtophysicalsize).to(''))
r_phys=r*pixtophysicalsize.to('pc')
print(f'physical radius: {r_phys}')

xpixs=np.arange(texpeakpix[0],(texpeakpix[0]+r))

texinradius=[]
xinradius=[]
yinradius=[]
centrtopix=[]

for y in range(np.shape(texmapdata)[0]):
    for x in range(np.shape(texmapdata)[1]):
        if (y-texpeakpix[0])**2+(x-texpeakpix[1])**2 <= r**2:
            texinradius.append(texmapdata[y,x].value)
            xinradius.append(x)
            yinradius.append(y)
            centrtopix.append((np.sqrt((y-texpeakpix[0])**2+(x-texpeakpix[1])**2)*pixtophysicalsize).value)
        else:
            pass

ax=plt.subplot(111)
plt.scatter(centrtopix,texinradius)
ax.set_xlabel('$d$ (AU)',fontsize=14)
ax.set_ylabel('$T_K$ (K)',fontsize=14)
ax.tick_params(size=14)
plt.tight_layout()

plt.show()
'''
edgepoints=circle(texmapdata,texpeakpix[0],texpeakpix[1],r)

plt.imshow(texmapdata.value,origin='lower',vmax=550,vmin=10)
plt.scatter(edgepoints[0],edgepoints[1],color='orange')
#plt.scatter(xinradius,yinradius,color='orange')
plt.show()
'''

