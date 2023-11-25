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
import math
import pdb
import matplotlib.ticker as mticker
from utilities import *
import sys
from astropy.table import QTable
#from astropy.visualization.wcsaxes import Quadrangle

#plt.close()

#mpl.interactive(True)

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

def make_tickstrings(list_of_float):
    list_of_strings=[]
    for f in list_of_float:
        order_of_mag=math.floor(np.log10(f))
        scaled=round(f/10**order_of_mag)#To revert back to the previous version, use math.floor here
        #pdb.set_trace()
        if scaled==1:
            scaled=10
            s1=r'$10^{'
            s2=fr'{order_of_mag}'
            s3=r'}$'
            string=s1+s2+s3
        else:
            end1=r'$\times10^{'
            end2=fr'{order_of_mag}'
            end3=r'}$'
            end=end1+end2+end3
            string=str(scaled)+end
        list_of_strings.append(string)
    return list_of_strings

set=3
colordict={0:('trot','inferno','bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'),1:('mom0','bone',"CH3OH~5_1-4_2E1vt0_masked.fits"),2:('nupper','Blues_r'),3:('detections','CMRmap',"ch3ohdetections0_3sigma_allspw_withnans_weighted.fits"),4:('abundance','viridis','bootstrap_ch3ohabundance_3sigma_ntotintercept_bolocamfeather_smoothedtobolocam.fits'),5:('nh2','Greys_r','bootstrap_nh2map_3sigma_bolocamfeather_smoothedtobolocam.fits')}
mode=colordict[set][0]
color=colordict[set][1]
pathsuffix=colordict[set][2]

print(f'Mode: {mode}')
cm= copy.copy(mpl.cm.get_cmap(color))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap, abundance cividis
cm.set_bad('black')
dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf
source=os.getenv('source')#
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
fnum=fields[source]

home=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/'

cntrfile='/orange/adamginsburg/sgrb2/2017.1.00114.S/imaging_results/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter_feathered_with_bolocam.fits'
tblfile='pacman_sep2023revolution/nov92023_allpropertytable.fits'

infile=home+source+sourcedict[source]+pathsuffix

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

hdu=fits.open(infile)[0]
cntrhdu=fits.open(cntrfile)[0]

cellsize=(np.abs(hdu.header['CDELT1']*u.deg)).to('arcsec')
pixtophysicalsize=(np.tan(cellsize)*dGC).to('AU')

cntrrms=0.0002#Jy
cntrlist=cntrrms*np.array([5,9,27,81,128,281])

cntrdata=np.squeeze(cntrhdu.data)

hduwcs=WCS(hdu)
cntrwcs=WCS(cntrhdu).celestial
texmaxpix=hduwcs.array_shape
hdubeam=radio_beam.Beam.from_fits_header(hdu.header)

plt.rcParams['figure.dpi'] = 300
sliced=['x','y']#,0,0]#Must be in same order as axes are given in fits header, at least. 
plt.figure()
ax=plt.subplot(projection=hduwcs,slices=sliced)
#ax.set_figheight(15)
#ax.set_figwidth(15)

ra=ax.coords[0]
dec=ax.coords[1]
if mode == 'trot':
    vmaxdict={'SgrB2S':352,'DSi':306,'DSii':239,'DSiii':351,'DSiv':398,'DSv':357,'DSVI':427,'DSVII':282,'DSVIII':300,'DSIX':252}
    img=ax.imshow(np.squeeze(hdu.data),vmax=vmaxdict[source],vmin=25,interpolation=None, cmap=cm)#, norm=mpl.colors.LogNorm())#vmaxcntm=0.005, vmaxdsi=300 (no min value),vmaxsgrb2s=605 tmin=10 (try no min value), ntotmax=6.26e17, dsintotmax=2.21e17
elif mode == 'abundance':
    abundadjust={'SgrB2S':5e-8,'DSiii':5e-9,'DSiv':9e-8,'DSv':9e-9,'DSVI':6e-8,'DSVII':2e-8,'DSVIII':3e-8,'DSIX':5e-8}
    maxadjust={'SgrB2S':1.13e-6,'DSi':9.3e-7,'DSiv':7.9e-7,'DSVI':6.1e-7,'DSVII':3.88e-7,'DSVIII':5.9e-7}#'DSiii':7e-8,
    if source in list(abundadjust.keys()) and source in list(maxadjust.keys()):
        print('bothadjust')
        maxfix=maxadjust[source]
        minfix=abundadjust[source]
        img=ax.imshow(np.squeeze(hdu.data),interpolation=None,cmap=cm, norm=mpl.colors.LogNorm(vmin=minfix,vmax=maxfix))
    elif source in list(abundadjust.keys()) or source in list(maxadjust.keys()):
        if source in list(abundadjust.keys()):
            print('abundadjust')
            minfix=abundadjust[source]
            img=ax.imshow(np.squeeze(hdu.data),interpolation=None,cmap=cm, norm=mpl.colors.LogNorm(vmin=minfix))#(vmax=5e-5,vmin=1e-7))
            maxfix=False
        if source in list(maxadjust.keys()):
            print('maxadjust')
            maxfix=maxadjust[source]
            img=ax.imshow(np.squeeze(hdu.data),interpolation=None,cmap=cm, norm=mpl.colors.LogNorm(vmax=maxfix))
    else:
        img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm, norm=mpl.colors.LogNorm())
        maxfix=False

elif mode == 'nh2':
    minnh2adjust={'SgrB2S':9e22,'DSVI':9e22,'DSVII':9e22}#,'DSi':9e22}
    maxnh2adjust={'DSii':9e23,'DSv':5e23,'DSIX':3e23}
    if source in minnh2adjust.keys():
        minfix=minnh2adjust[source]
        img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm, norm=mpl.colors.LogNorm(vmin=minfix,))
    elif source in maxnh2adjust.keys():
        maxfix=maxnh2adjust[source]
        img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm, norm=mpl.colors.LogNorm(vmax=maxfix,))
    else:
        img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm, norm=mpl.colors.LogNorm())
else:
    img=ax.imshow(np.squeeze(hdu.data),interpolation=None, cmap=cm)
lims=ax.axis()
if mode == 'detections':
    ax.contour(cntrdata, levels=cntrlist, colors='black',transform=ax.get_transform(cntrwcs),linewidths=0.5,zorder=1)
    
    tbl=QTable.read(tblfile)
    tblindex=np.where(tbl['Source']==tblconversion[source])[0]
    srcradius=tbl['Radius'][tblindex]
    radius_pixels=math.floor((srcradius/pixtophysicalsize).to('').value)
    if source == 'SgrB2S':
        circle=mpl.patches.Circle(xy=(65,70),radius=radius_pixels,facecolor='none',edgecolor='cyan',linewidth=1)
    else:
        circle=mpl.patches.Circle(xy=(pixdict[source][1],pixdict[source][0]),radius=radius_pixels,facecolor='none',edgecolor='cyan',linewidth=1)
    ax.scatter(pixdict[source][1],pixdict[source][0],marker='*',color='cyan')
    ax.add_patch(circle)
else:
    ax.contour(cntrdata, levels=cntrlist, colors='white',transform=ax.get_transform(cntrwcs),linewidths=0.5,zorder=1)

scaledict={'SgrB2S':5000*u.AU,'DSi':5000*u.AU,'DSii':2000*u.AU,'DSiii':2000*u.AU,'DSiv':2000*u.AU,'DSv':2000*u.AU,'DSVI':5000*u.AU,'DSVII':5000*u.AU,'DSVIII':5000*u.AU,'DSIX':5000*u.AU}
scale=scaledict[source]
lenn=np.arctan(scale/dGC)
h_img=hdu.header['NAXIS1']
l_img=hdu.header['NAXIS2']

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
              text_offset=0.015*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.0549 -28:23:49.950', unit=(u.hour,u.deg), frame='icrs')
    
elif source == 'DSv':
    scalebar_vert=math.ceil(h_img/10)
    scalebar_hor=math.ceil(l_img/10)
    scalebar_position=hduwcs.wcs_pix2world([[scalebar_vert,scalebar_hor],[0,0]],0)

    make_scalebar(ax, coordinates.SkyCoord(scalebar_position[0,0],scalebar_position[0,1], unit=u.deg, 
                                           frame='icrs'),
                  length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
                  text_offset=0.015*u.arcsec)
    
else:#Find the ideal corner position for the scalebar automatically
    scalebar_vert=math.ceil(h_img/10)
    scalebar_hor=math.ceil(l_img/10)
    scalebar_position=hduwcs.wcs_pix2world([[scalebar_vert,scalebar_hor],[0,0]],0)

    make_scalebar(ax, coordinates.SkyCoord(scalebar_position[0,0],scalebar_position[0,1], unit=u.deg, 
                                           frame='icrs'),
                  length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
                  text_offset=0.02*u.arcsec)
'''
else:
    print('No scalebar available for this object')
'''
pixscale=np.mean(wcs.utils.proj_plane_pixel_scales(hduwcs))*u.deg
if source == 'SgrB2S':
    beam_vert=math.ceil(h_img/20)
    beam_hor=math.ceil(l_img/17)
elif source == 'DSi':
    beam_vert=math.ceil(h_img/20)
    beam_hor=math.ceil(l_img/1.15)
else:
    beam_vert=math.ceil(h_img/7.78)#20 for DSi
    beam_hor=math.ceil(l_img/1.17)
hdubeamplot=hdubeam.ellipse_to_plot(beam_hor,beam_vert,pixscale)
ax.add_patch(hdubeamplot)

hdubeamplot.set_facecolor('gray')
hdubeamplot.set_edgecolor('white')

labeldict={0:'T$_{K}$ (K)',1:r'Intensity (K km s$^{-1}$)',2:r'N$_{upper}$ (cm$^{-2}$)',3:'$n_{transition}$',4:'X(CH$_3$OH)',5:'$N$(H$_2$) (cm$^{-2}$)'}
ax.axis(lims)
ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(img,label=labeldict[set],pad=0)
docbartickfix=['DSii','DSIX']
skipcbartickfix=['SgrB2S']

if mode == 'abundance':
    if maxfix:
        #ticklimit=False
        #cbarticktext=cbar1.ax.get_yticklabels()
        cbarticks=cbar1.get_ticks()
        tickbools=maxfix<cbarticks
        index=-1
        for boo in tickbools:
            #endloop=False
            if not boo:
                index+=1
            elif source in skipcbartickfix:
                print(f'Manually skipping {source} cbar correction.')
                break
            else:
                #lowertick=cbarticks[index]
                #tickconvert=maxfix/lowertick
                #cbarticks[index+1]=maxfix
                cbarticks=np.insert(cbarticks,(index+1),maxfix)
                ticks_loc=cbarticks.tolist()
                cbar1.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
                strcbarticks=make_tickstrings(cbarticks)
                cbar1.ax.set_yticklabels(strcbarticks)#print(cbarticks)
                break
                
    elif source in docbartickfix:
        print(f'Manual cbar correciton for {source} yay')
        maxfix=np.nanmax(np.squeeze(hdu.data))
        cbarticks=cbar1.get_ticks()
        tickbools=maxfix<cbarticks
        index=-1
        for boo in tickbools:
            #endloop=False
            if not boo:
                index+=1
            elif source in skipcbartickfix:
                print(f'Manually skipping {source} cbar correction.')
                break
            else:
                #lowertick=cbarticks[index]
                #tickconvert=maxfix/lowertick
                #cbarticks[index+1]=maxfix
                cbarticks=np.insert(cbarticks,(index+1),maxfix)
                ticks_loc=cbarticks.tolist()
                cbar1.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
                strcbarticks=make_tickstrings(cbarticks)
                cbar1.ax.set_yticklabels(strcbarticks)#print(cbarticks)
                break

#pdb.set_trace()
#sys.exit()
print(f'Saving figure at {savefigpath}')
plt.savefig(savefigpath)
print('Saved')

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
              
#plt.show()

