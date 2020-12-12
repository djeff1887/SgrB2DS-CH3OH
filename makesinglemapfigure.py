import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import glob
import numpy as np
import copy
import matplotlib as mpl


cm= copy.copy(mpl.cm.get_cmap("inferno"))
cm.set_bad('black')

home='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/'
sourcepath=home+'SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/'+'figures/'
saveimgpath=sourcepath+'texmap_ntotcontours-2nsig.png'

cntrfile="/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/SgrB2DS_field1_spw2_cube_robust0_niter1e6_chunks16_minimize.image.pbcor_continuum.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#input('Map fits file: ')

infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
'''
if '"' in infile:
    infile.replace('"','')
    print('" Replaced')
else:
    print('No quotations found.')
    pass
'''
hdu=fits.open(infile)[0]
cntrhdu=fits.open(cntrfile)[0]

'''
Stopgap masking values that are unphysically large
cntrshape=np.shape(cntrhdu.data)
upperntot=1e20
for y in range(cntrshape[0]):
    for x in range(cntrshape[1]):
        if cntrhdu.data[y,x] > upperntot:
            cntrhdu.data[y,x]=np.nan
        else:
            pass
assert upperntot not in cntrhdu.data, 'Unphysical values in Ntot image'
'''
cntrrms=np.nanstd(cntrhdu.data)
cntrlist=cntrrms*np.array([1,2,4,8,16,32,64])

cntrdata=cntrhdu.data

hduwcs=WCS(hdu)

sliced=['y','x']
ax=plt.subplot(projection=hduwcs,slices=sliced)

ra=ax.coords[1]
dec=ax.coords[0]

img=ax.imshow(hdu.data,vmax=550,cmap=cm)#tmax=550, tmin=10, ntotmax=6.26e17, dsintotmax=2.21e17

ax.contour(cntrhdu.data, levels=cntrlist, colors='white')#, alpha=0.5)#ax.contour(data=hdu.data)#, colors='black')#binary/Greys are other good cmaps
#ax.contour(cntrhdu.data,levels=(-1*cntrlist),colors='white',linestyles='dashed')#YlGn is alt for both
    
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=0.1)
ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.5)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(img,pad=0)#plt.colorbar()
#plt.savefig(saveimgpath)
plt.show()

