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
import matplotlib as mpl
import os
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
set=0
colordict={0:('trot','inferno','texmap_3sigma_allspw_withnans_weighted.fits'),1:('mom0','bone'),2:('nupper','Blues_r'),3:('detections','CMRmap',"ch3ohdetections5_3sigma_allspw_withnans_weighted.fits"),4:('abundance','viridis','ch3ohabundance_3sigma_ntotintercept.fits')}
mode=colordict[set][0]
color=colordict[set][1]
pathsuffix=colordict[set][2]

print(f'Mode: {mode}')
cm= copy.copy(mpl.cm.get_cmap(color))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap, abundance cividis
cm.set_bad('black')
dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
source='SgrB2S'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3}
fnum=fields[source]

home=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/'
sourcepath=home+f'{source}/field10originals_z0_000186431_5-6mhzwidth_stdfixes/'+'figures/'
#saveimgpath=sourcepath+'texmap_ntotcontours-3-6-8to48.png'

cntrfile="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.pbcor.fits"
#cntrfile="/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/SgrB2DS_field1_spw2_cube_robust0_niter1e6_chunks16_minimize.image.pbcor_continuum.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/ntotmap_allspw_withnans_weighted.fits"#input('Map fits file: ')
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/200K_field10originals_z0_000190713_5-6mhzwidth_stdfixes/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/200K_field10originals_z0_000190713_5-6mhzwidth_stdfixes/mom0/CH3OH~25_3-24_4E1vt0.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/200K_field10originals_z0_000190713_5-6mhzwidth_stdfixes/ch3ohdetections3_3sigma_allspw_withnans_weighted.fits"
#infile='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/200K_field10originals_z0_000190713_5-6mhzwidth_stdfixes/texmap_5transmask_3sigma_allspw_withnans_weighted.fits'

#infile=cntrfile

sourcedict={'SgrB2S':'/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/Kfield10originals_noexclusions/','DSiii':'/Kfield10originals_noexclusions/','DSiv':'/Kfield10originals_noexclusions/','DSv':f'/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':'/Kfield2originals_trial3_8_6-8_7excluded/','DSVII':'/Kfield3originals_200K_trial1_noexclusions/','DSVIII':'/Kfield3originals_175K_trial1_noexclusions/'}#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial4_3kmsslab_newexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/field10/CH3OH/DSi/field10originals_spatialandvelocitymaskingtrial4_3kmsslab_newexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_large_z0_0002306756533745274_5-6mhzwidth_stdfixes/ch3ohdetections5_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/employingrepresentativelineprocedure_K_OctReimage_restfreqfix_newvelmask_newpeakamp/mom0/CH3OH~10_2+-9_3+vt0_masked.fits"
#infile="/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/peakintensity6-7_octreimageKlarge.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_large_z0_0002306756533745274_5-6mhzwidth_stdfixes/mom0/CH3OH~6_1--7_2-vt1.fits"
#infile="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_robust0.image.tt0.pbcor.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/ch3ohdetections3_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/texmap_3transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/ch3ohdetections3_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#infile="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"

infile=home+source+'/'+sourcedict[source]+pathsuffix

savefigbase=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}'
savefighome=savefigbase+sourcedict[source]

if not os.path.exists(savefighome):
    makedir=input(f'Create savefigpath {savefighome}?  (y/n)')
    if makedir=='y':
        os.makedirs(savefighome)
    elif makedir=='n':
        print('Edit script to have correct savefighome')
        pdb.set_trace()
else:
    print(f'Savefigpath {savefighome} already exists')
    pass 

pathsuffix2=pathsuffix.replace('fits','png')
savefigpath=savefighome+pathsuffix2

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
cntrlist=cntrrms*np.array([3,6,8,16,32])

cntrdata=np.squeeze(cntrhdu.data)

hduwcs=WCS(hdu)
cntrwcs=WCS(cntrhdu).celestial
texmaxpix=hduwcs.array_shape
hdubeam=radio_beam.Beam.from_fits_header(hdu.header)

sliced=['x','y']#,0,0]#Must be in same order as axes are given in fits header, at least. 
ax=plt.subplot(projection=hduwcs,slices=sliced)
#ax.set_figheight(15)
#ax.set_figwidth(15)
plt.rcParams['figure.dpi'] = 150

ra=ax.coords[0]
dec=ax.coords[1]
if mode == 'trot':
    vmaxdict={'SgrB2S':520,'DSi':320,'DSii':225,'DSiii':300,'DSiv':313,'DSv':281,'DSVI':378,'DSVII':250,'DSVIII':225}
    img=ax.imshow(np.squeeze(hdu.data),vmax=vmaxdict[source],interpolation=None, cmap=cm)#, norm=mpl.colors.LogNorm())#vmaxcntm=0.005, vmaxdsi=300 (no min value),vmaxsgrb2s=605 tmin=10 (try no min value), ntotmax=6.26e17, dsintotmax=2.21e17
elif mode == 'abundance':
    img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm, norm=mpl.colors.LogNorm(vmax=5e-5,vmin=1e-7))
else:
    img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm)
lims=ax.axis()
ax.contour(cntrdata, levels=cntrlist, colors='white',transform=ax.get_transform(cntrwcs),linewidths=0.5)#, alpha=0.5)#ax.contour(data=hdu.data)#, colors='black')#binary/Greys are other good cmaps
#ax.contour(cntrdata,levels=(-1*cntrlist),colors='white',linestyles='dashed')#YlGn is alt for both

scaledict={'SgrB2S':5000*u.AU,'DSi':5000*u.AU,'DSii':2000*u.AU,'DSiii':2000*u.AU,'DSiv':2000*u.AU,'DSv':2000*u.AU,'DSVI':5000*u.AU,'DSVII':5000*u.AU,'DSVIII':5000*u.AU}
scale=scaledict[source]
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
    
elif source == 'DSii':
    print(f'Scalebar source: DSii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.1115 -28:23:47.596', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.1087 -28:23:48.634', unit=(u.hour,u.deg), frame='icrs')
    pass
elif source == 'DSiii':
    print(f'Scalebar source: DSiii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.0577 -28:23:49.912', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.0549 -28:23:49.950', unit=(u.hour,u.deg), frame='icrs')
    
elif source=='DSall':
    make_scalebar(ax, coordinates.SkyCoord('17:47:24.0 -28:26:25', unit=(u.hour, u.deg), 
                                           frame='icrs'),
                  length=lenn, label=f'{scale.value} pc', fontsize=12, 
                  text_offset=0.02*u.arcsec)
else:
    print('No scalebar available for this object')

pixscale=np.mean(wcs.utils.proj_plane_pixel_scales(hduwcs))*u.deg
hdubeamplot=hdubeam.ellipse_to_plot(10,10,pixscale)
hdubeamplot.set_facecolor('none')
hdubeamplot.set_edgecolor('cyan')

labeldict={0:'T$_{K}$ (K)',1:r'S$_\nu$ (K)',2:r'N$_{upper}$ (cm$^{-2}$)',3:'$n_{transition}$',4:'X(CH$_3$OH)'}
ax.axis(lims)
ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(img,pad=0,label=labeldict[set])#r'S$_\nu$ (Jy)')#'$n_{transition}$')#plt.colorbar()
plt.savefig(savefigpath)

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

