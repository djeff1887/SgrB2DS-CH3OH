import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import glob
import numpy as np
import copy
import matplotlib as mpl
import astropy.units as u
import radio_beam
from astropy import wcs
import astropy.units as u
from astropy import coordinates
from matplotlib.patches import Rectangle
#from astropy.visualization.wcsaxes import Quadrangle

plt.close()

mpl.interactive(True)

def make_scalebar(ax, left_side, length, color='black', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec, coordsys='icrs'):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform(coordsys),
                   zorder=3)
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform(coordsys),
                  color=color,
                  fontsize=fontsize,
                 zorder=2,bbox=dict(facecolor='white', alpha=0.6))
    ax.axis(axlims)
    return lines,txt

cm= copy.copy(mpl.cm.get_cmap("bone"))#mom0 bone, temp inferno, nupper Blues_r, detections CMRmap
cm.set_bad('black')
dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
source='SgrB2'
fnum=1

home=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/'
sourcepath=home+f'{source}/field10originals_z0_000186431_5-6mhzwidth_stdfixes/'+'figures/'
saveimgpath=sourcepath+'texmap_ntotcontours-3-6-8to48.png'

cntrfile="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.pbcor.fits"
#cntrfile="/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/SgrB2DS_field1_spw2_cube_robust0_niter1e6_chunks16_minimize.image.pbcor_continuum.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#input('Map fits file: ')

infile="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_robust0.image.tt0.pbcor.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/ch3ohdetections3_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/texmap_3transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/ch3ohdetections3_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/texmap_3transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
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
cntrlist=cntrrms*np.array([3,6,8,16,32,48])

cntrdata=np.squeeze(cntrhdu.data)

hduwcs=WCS(hdu)
cntrwcs=WCS(cntrhdu).celestial
texmaxpix=hduwcs.array_shape
hdubeam=radio_beam.Beam.from_fits_header(hdu.header)

sliced=['x','y']#Must be in same order as axes are given in fits header, at least. 
ax=plt.subplot(projection=hduwcs,slices=sliced)
#ax.set_figheight(15)
#ax.set_figwidth(15)
plt.rcParams['figure.dpi'] = 150

ra=ax.coords[0]
dec=ax.coords[1]

img=ax.imshow(hdu.data,vmax=0.005,interpolation=None, cmap=cm)#vmaxdsi=300,vmaxsgrb2s=605 tmin=10, ntotmax=6.26e17, dsintotmax=2.21e17
lims=ax.axis()
#ax.contour(cntrdata, levels=cntrlist, colors='white',transform=ax.get_transform(cntrwcs),linewidths=1)#, alpha=0.5)#ax.contour(data=hdu.data)#, colors='black')#binary/Greys are other good cmaps
#ax.contour(cntrhdu.data,levels=(-1*cntrlist),colors='white',linestyles='dashed')#YlGn is alt for both

scale=5000*u.AU
lenn=np.arctan(scale/dGC)
#SgrB2S
if source=='SgrB2S':
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.6618 -28:23:48.734', unit=(u.hour, u.deg), 
                                           frame='icrs'),
                  length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
                  text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.6590 -28:23:48.772', unit=(u.hour,u.deg), frame='icrs')


#DSi
elif source=='DSi':
    make_scalebar(ax, coordinates.SkyCoord('17:47:19.5180 -28:23:51.359', unit=(u.hour, u.deg), 
                                           frame='icrs'),
                  length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
                  text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:19.4976 -28:23:51.384', unit=(u.hour,u.deg), frame='icrs')
else:
    print('No scalebar available for this object')

pixscale=np.mean(wcs.utils.proj_plane_pixel_scales(hduwcs))*u.deg
hdubeamplot=hdubeam.ellipse_to_plot(10,10,pixscale)
hdubeamplot.set_facecolor('none')
hdubeamplot.set_edgecolor('cyan')
    
ax.axis(lims)
ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(img,pad=0,label='mJy')#'n$_{transition}$')#plt.colorbar()
#plt.savefig(saveimgpath)

'''Plotting shapes
r = Rectangle((2000,2000),1000,1000,facecolor='red')
s = Rectangle((2000,2000),500,500,facecolor='green')
t = Rectangle((2000,2000),100,100,facecolor='blue')
f=Rectangle((266.8350000, 28.3958333), 0.1111667, 0.1111667, label='EMIR FOV',edgecolor='orange', facecolor='orange', linewidth=100,linestyle='-',transform=ax.get_transform('icrs'))

ax.add_patch(r)
ax.add_patch(s)
ax.add_patch(t)
ax.add_patch(f)

plt.legend(loc='upper right')
'''
              
plt.show()

