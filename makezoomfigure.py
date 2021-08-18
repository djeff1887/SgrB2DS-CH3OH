import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import visualization
import glob
import math

cm=plt.cm.get_cmap('inferno')
cm.set_bad('black')

#sourcepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/'
sgrb2dspath="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.pbcor.fits"
sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
#sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsiitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/Kfield10originals_noexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
sgrb2dsiiitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
sgrb2dsivtexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
dsvtexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_150K_trial2/texmap_0transmask_3sigma_allspw_withnans_weighted.fits"
ds6texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial2_16_6-16_7excluded/texmap_5transmask_3sigma_allspw_withnans_weighted.fits"
ds7texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_trial1_noexclusions/texmap_3sigma_allspw_withnans_weighted.fits"

sgrb2dshdu=fits.open(sgrb2dspath)[0]
sgrb2dsdata=sgrb2dshdu.data.squeeze()
sgrb2dswcs=WCS(sgrb2dshdu.header).celestial

sgrb2shdu=fits.open(sgrb2stexmap)[0]
sgrb2dsihdu=fits.open(sgrb2dsitexmap)[0]
dsiihdu=fits.open(sgrb2dsiitexmap)[0]
dsiiihdu=fits.open(sgrb2dsiiitexmap)[0]
dsivhdu=fits.open(sgrb2dsivtexmap)[0]
dsvhdu=fits.open(dsvtexmap)[0]
ds6hdu=fits.open(ds6texmap)[0]
ds7hdu=fits.open(ds7texmap)[0]

tmax=520
jymax=0.050831314
jymaxfull=0.01
axins_dims=0.45

centerx=1453
centery=3317
width=(122/2)

centerx2=1687
centery2=3246
width2=(74/2)

centerx3=1567
centery3=3309
width3=(46/2)

centerx4=1583
centery4=3265
width4=(50/2)

centerx5=1640
centery5=3371
width5=(65/2)

centerx6=1655
centery6=3211
width6=(40/2)

centerx7=1282
centery7=2677
width7=(125/2)

centerx8=992
centery8=2363
width8=(150/2)

sliced=["x","y"]#[0,0,"y","x"]
ax=plt.subplot(projection=sgrb2dswcs,slices=sliced)
ra=ax.coords[0]
dec=ax.coords[1]
axins=ax.inset_axes([-0.75,0.75,axins_dims,axins_dims])
axins.imshow(sgrb2dsdata,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray')
axins.set_xlim((centerx-width),(centerx+width))
axins.set_ylim((centery-width),(centery+width))
ax.imshow(sgrb2dsdata, origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymaxfull, min_cut=0),cmap='gray_r')
axins2=axins.inset_axes([-1.25,0,1,1])
axins2.imshow(sgrb2shdu.data,vmax=tmax,origin='lower',cmap='inferno')

axins3=ax.inset_axes([1.25,0.15,axins_dims,axins_dims])
axins3.imshow(sgrb2dsdata,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray')
axins3.set_xlim((centerx2-width2),(centerx2+width2))
axins3.set_ylim((centery2-width2),(centery2+width2))
axins4=axins3.inset_axes([1.25,0,1,1])
axins4.imshow(sgrb2dsihdu.data, vmax=tmax,origin='lower',cmap='inferno')

axins5=ax.inset_axes([1.25,0.45,axins_dims,axins_dims])
axins5.imshow(sgrb2dsdata,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray')
axins5.set_xlim((centerx3-width3),(centerx3+width3))
axins5.set_ylim((centery3-width3),(centery3+width3))
axins6=axins5.inset_axes([1.25,0,1,1])
axins6.imshow(dsiihdu.data, vmax=tmax,origin='lower',cmap='inferno')

axins7=ax.inset_axes([1.25,-0.15,axins_dims,axins_dims])
axins7.imshow(sgrb2dsdata,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray')
axins7.set_xlim((centerx4-width4),(centerx4+width4))
axins7.set_ylim((centery4-width4),(centery4+width4))
axins8=axins7.inset_axes([1.25,0,1,1])
axins8.imshow(dsiiihdu.data, vmax=tmax,origin='lower',cmap='inferno')

axins9=ax.inset_axes([1.25,0.75,axins_dims,axins_dims])
axins9.imshow(sgrb2dsdata,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray')
axins9.set_xlim((centerx5-width5),(centerx5+width5))
axins9.set_ylim((centery5-width5),(centery5+width5))
axins10=axins9.inset_axes([1.25,0,1,1])
axins10.imshow(dsivhdu.data, vmax=tmax,origin='lower',cmap='inferno')

axins11=ax.inset_axes([-0.75,0.45,axins_dims,axins_dims])
axins11.imshow(sgrb2dsdata,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray')
axins11.set_xlim((centerx4-width4),(centerx4+width4))
axins11.set_ylim((centery4-width4),(centery4+width4))
axins12=axins11.inset_axes([-1.25,0,1,1])
axins12.imshow(dsiiihdu.data,vmax=tmax,origin='lower',cmap='inferno')

axins13=ax.inset_axes([-0.75,0.15,axins_dims,axins_dims])
axins13.imshow(sgrb2dsdata,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray')
axins13.set_xlim((centerx7-width7),(centerx7+width7))
axins13.set_ylim((centery7-width7),(centery7+width7))
axins14=axins13.inset_axes([-1.25,0,1,1])
axins14.imshow(ds6hdu.data,vmax=tmax,origin='lower',cmap='inferno')

axins15=ax.inset_axes([-0.75,-0.15,axins_dims,axins_dims])
axins15.imshow(sgrb2dsdata,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray')
axins15.set_xlim((centerx8-width8),(centerx8+width8))
axins15.set_ylim((centery8-width8),(centery8+width8))
axins16=axins15.inset_axes([-1.25,0,1,1])
axins16.imshow(ds7hdu.data,vmax=tmax,origin='lower',cmap='inferno')

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
axins5.tick_params(direction='in')
axins6.tick_params(direction='in')
axins7.tick_params(direction='in')
axins8.tick_params(direction='in')
axins9.tick_params(direction='in')
axins10.tick_params(direction='in')
axins11.tick_params(direction='in')
axins12.tick_params(direction='in')
axins13.tick_params(direction='in')
axins14.tick_params(direction='in')
axins15.tick_params(direction='in')
axins16.tick_params(direction='in')

ra.set_ticklabel_visible(True)
dec.set_ticklabel_visible(True)
axins.xaxis.set_ticklabels([])
axins.yaxis.set_ticklabels([])
axins2.xaxis.set_ticklabels([])
axins2.yaxis.set_ticklabels([])
axins3.xaxis.set_ticklabels([])
axins3.yaxis.set_ticklabels([])
axins4.xaxis.set_ticklabels([])
axins4.yaxis.set_ticklabels([])
axins5.xaxis.set_ticklabels([])
axins5.yaxis.set_ticklabels([])
axins6.xaxis.set_ticklabels([])
axins6.yaxis.set_ticklabels([])
axins7.xaxis.set_ticklabels([])
axins7.yaxis.set_ticklabels([])
axins8.xaxis.set_ticklabels([])
axins8.yaxis.set_ticklabels([])
axins9.xaxis.set_ticklabels([])
axins9.yaxis.set_ticklabels([])
axins10.xaxis.set_ticklabels([])
axins10.yaxis.set_ticklabels([])
axins11.xaxis.set_ticklabels([])
axins11.yaxis.set_ticklabels([])
axins12.xaxis.set_ticklabels([])
axins12.yaxis.set_ticklabels([])
axins13.xaxis.set_ticklabels([])
axins13.yaxis.set_ticklabels([])
axins14.xaxis.set_ticklabels([])
axins14.yaxis.set_ticklabels([])
axins15.xaxis.set_ticklabels([])
axins15.yaxis.set_ticklabels([])
axins16.xaxis.set_ticklabels([])
axins16.yaxis.set_ticklabels([])

ax.indicate_inset_zoom(axins)
ax.indicate_inset_zoom(axins3)
ax.indicate_inset_zoom(axins5)
ax.indicate_inset_zoom(axins7)
ax.indicate_inset_zoom(axins9)
ax.indicate_inset_zoom(axins11)
ax.indicate_inset_zoom(axins13)
ax.indicate_inset_zoom(axins15)
#axins.indicate_inset_zoom(axins2)
plt.show()