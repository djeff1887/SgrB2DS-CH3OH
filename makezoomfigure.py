import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import glob

cm=plt.cm.get_cmap('inferno')
cm.set_bad('black')

#sourcepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/'
sgrb2dspath="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.pbcor.fits"
sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/field10originals_z0_000186431_5-6mhzwidth_stdfixes/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"

sgrb2dshdu=fits.open(sgrb2dspath)[0]
sgrb2dswcs=WCS(sgrb2dshdu)

sgrb2shdu=fits.open(sgrb2stexmap)[0]
sgrb2dsihdu=fits.open(sgrb2dsitexmap)[0]


centerx=1450
centery=3330
width=65

centerx2=1690
centery2=3245
width2=25

sliced=["x","y",0,0]#[0,0,"y","x"]
ax=plt.subplot(projection=sgrb2dswcs,slices=sliced)
ra=ax.coords[0]
dec=ax.coords[1]
axins=ax.inset_axes([-0.75,0.33,0.5,0.5])
axins.imshow(sgrb2dshdu.data[0,0,:,:],origin='lower',cmap='gray_r')
axins.set_xlim((centerx-width),(centerx+width))
axins.set_ylim((centery-width),(centery+width))
ax.imshow(sgrb2dshdu.data[0,0,:,:], origin='lower',cmap='gray_r')
axins2=axins.inset_axes([-1.25,0,1,1])
axins2.imshow(sgrb2shdu.data,vmin=73,vmax=605,origin='lower',cmap='inferno')

axins3=ax.inset_axes([1.25,0.33,0.5,0.5])
axins3.imshow(sgrb2dshdu.data[0,0,:,:],origin='lower',cmap='gray_r')
axins3.set_xlim((centerx2-width2),(centerx2+width2))
axins3.set_ylim((centery2-width2),(centery2+width2))
axins4=axins3.inset_axes([1.25,0,1,1])
axins4.imshow(sgrb2dsihdu.data,vmin=73,vmax=605,origin='lower',cmap='inferno')
#plt.grid(color='white', ls='solid')
dec.set_axislabel('Dec')
ra.set_axislabel('RA')
#ra.set_ticks_visible(False)
#dec.set_ticks_visible(False)
ax.tick_params(direction='in')
axins.tick_params(direction='in')
axins2.tick_params(direction='in')
axins3.tick_params(direction='in')
axins4.tick_params(direction='in')
ra.set_ticklabel_visible(False)
dec.set_ticklabel_visible(False)
axins.xaxis.set_ticklabels([])
axins.yaxis.set_ticklabels([])
axins2.xaxis.set_ticklabels([])
axins2.yaxis.set_ticklabels([])
axins3.xaxis.set_ticklabels([])
axins3.yaxis.set_ticklabels([])
axins4.xaxis.set_ticklabels([])
axins4.yaxis.set_ticklabels([])
ax.indicate_inset_zoom(axins)
ax.indicate_inset_zoom(axins3)
#axins.indicate_inset_zoom(axins2)
plt.show()