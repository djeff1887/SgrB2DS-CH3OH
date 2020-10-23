import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import glob


cm=plt.cm.get_cmap('Blues_r')
cm.set_bad('black')

home='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/'
sourcepath=home+'DSi/z0_000186407_box1_5-6mhzwidth/'+'figures/'
saveimgpath=sourcepath+'ntotmap_allspw_withnans_weighted_sgrb2dsscale.png'

infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#input('Map fits file: ')

'''
if '"' in infile:
    infile.replace('"','')
    print('" Replaced')
else:
    print('No quotations found.')
    pass
'''
hdu=fits.open(infile)[0]

hduwcs=WCS(hdu)

sliced=['y','x']
ax=plt.subplot(projection=hduwcs,slices=sliced)

ra=ax.coords[1]
dec=ax.coords[0]

img=ax.imshow(hdu.data,vmax=6.26e17,cmap=cm)#tmax=550, tmin=10, ntotmax=6.26e17, dsintotmax=2.21e17

dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=0.1)
ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.5)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(img,pad=0)#plt.colorbar()
plt.savefig(saveimgpath)
plt.show()

