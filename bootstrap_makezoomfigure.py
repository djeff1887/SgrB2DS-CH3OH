import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy import visualization
import glob
import math
import numpy as np
import astropy.units as u
import matplotlib.patheffects as pe
import matplotlib as mpl
import radio_beam
import pdb
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


plt.close('all')
mpl.interactive(True)

cm=plt.cm.get_cmap('inferno')
cm.set_bad('black')

#sourcepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/'
sgrb2dspath='/orange/adamginsburg/sgrb2/2017.1.00114.S/imaging_results/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter_feathered_with_bolocam.fits'
notfeath_sgrb2dspath="/blue/adamginsburg/d.jeff/imaging_results/adamcleancontinuum/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust0_selfcal4_finaliter.image.tt0.pbcor.fits"# Using this one b/c the feathered one has that weird background
sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
#sgrb2stexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/nov2022continuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
#sgrb2dsitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/DSi/z0_000186407_box1_5-6mhzwidth/texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsiitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/nov2022continuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsiiitexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/nov2022continuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
sgrb2dsivtexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/nov2022contniuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
dsvtexmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/nov2022contniuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
ds6texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/nov2022continuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
ds7texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/nov2022contniuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
ds8texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVIII/nov2022contniuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"
ds9texmap="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/nov2022contniuumsanitycheck/bootstrap_texmap_3sigma_allspw_withnans_weighted.fits"

sgrb2dshdu=fits.open(sgrb2dspath)[0]
sgrb2dsdata=sgrb2dshdu.data.squeeze()
sgrb2dshdu2=fits.open(notfeath_sgrb2dspath)[0]
sgrb2dsdata2=sgrb2dshdu2.data.squeeze()
sgrb2dswcs=WCS(sgrb2dshdu.header).celestial

contrestfreq=sgrb2dshdu.header['RESTFRQ']*u.Hz
equiv=u.brightness_temperature(contrestfreq)
beam_sgrb2dshdu=radio_beam.Beam.from_fits_header(sgrb2dshdu.header)
beam_sgrb2dshdu=beam_sgrb2dshdu.sr
sgrb2dsdata=((sgrb2dsdata/beam_sgrb2dshdu)*u.Jy).to('K',equivalencies=equiv).value

sgrb2shdu=fits.open(sgrb2stexmap)[0]
sgrb2dsihdu=fits.open(sgrb2dsitexmap)[0]
dsiihdu=fits.open(sgrb2dsiitexmap)[0]
dsiiihdu=fits.open(sgrb2dsiiitexmap)[0]
dsivhdu=fits.open(sgrb2dsivtexmap)[0]
dsvhdu=fits.open(dsvtexmap)[0]
ds6hdu=fits.open(ds6texmap)[0]
ds7hdu=fits.open(ds7texmap)[0]
ds8hdu=fits.open(ds8texmap)[0]
ds9hdu=fits.open(ds9texmap)[0]

'''for alma'''
pixscale=0.05*u.arcsec
hpbw=42.9*u.arcsec
pbdiam=math.ceil(hpbw/pixscale)
pbrad=pbdiam/2
theta=np.linspace(0,2*np.pi,150)
a=pbrad*np.cos(theta)
b=pbrad*np.sin(theta)

fluxscalefactor=1#000
tmax=575
tmin=20
jymax=((0.05*fluxscalefactor/beam_sgrb2dshdu)*u.Jy).to('K',equivalencies=equiv).value#Was 0.050831314 Jy
jymaxfull=((0.05*fluxscalefactor/beam_sgrb2dshdu)*u.Jy).to('K',equivalencies=equiv).value#was 0.01 Jy
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

centerx9=1040
centery9=2191
width9=(101/2)

centerx10=667
centery10=800
width10=(69/2)

lefthoriz=-0.75
righthoriz=1.25

cntrlist=jymax*np.array([1,1.5,2,2.5,3])

ok=np.isfinite(sgrb2dsdata)*np.isfinite(sgrb2dsdata2)
maskfeather=np.ma.masked_where(ok==False,sgrb2dsdata)
sgrb2dsdata=maskfeather.filled(fill_value=np.nan)

annote_size=15
plt.rcParams["figure.dpi"]=150
sliced=["x","y"]#[0,0,"y","x"]
fig=plt.figure(dpi=150,figsize=(8,7))#tight_layout=True) 16,7
ax=plt.subplot(projection=sgrb2dswcs,slices=sliced)
ax.axis([348,2164,448,3735])#This sets the borders of the data in the continuum
ra=ax.coords[0]
dec=ax.coords[1]
axins=ax.inset_axes([lefthoriz,0.864,axins_dims,axins_dims])
axins.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins.set_xlim((centerx-width),(centerx+width))
axins.set_ylim((centery-width),(centery+width))
axins.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
continuum=ax.imshow((sgrb2dsdata*fluxscalefactor), origin='lower',norm=visualization.simple_norm((sgrb2dsdata*fluxscalefactor), stretch='sqrt', max_cut=jymaxfull, min_cut=0),cmap='gray_r')
axins2=axins.inset_axes([-1.25,0,1,1])
axins2show=axins2.imshow(sgrb2shdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
s_txt=axins2.text(3,102,'Sgr B2(S)',size=annote_size,color='black')
s_txt.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

shape=sgrb2shdu.data.shape
xratio=3/shape[1]
yratio=102/shape[0]

axins3=ax.inset_axes([righthoriz,0.264,axins_dims,axins_dims])
axins3.imshow(sgrb2dsdata*fluxscalefactor,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray_r')
axins3.set_xlim((centerx2-width2),(centerx2+width2))
axins3.set_ylim((centery2-width2),(centery2+width2))
axins3.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins4=axins3.inset_axes([1.25,0,1,1])
axins4.imshow(sgrb2dsihdu.data, vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape1=sgrb2dsihdu.data.shape
xpos1=round(shape1[1]*xratio)
ypos1=round(shape1[0]*yratio)
txt1=axins4.text(xpos1,ypos1,'DS1',size=annote_size,color='black')
txt1.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins5=ax.inset_axes([righthoriz,0.564,axins_dims,axins_dims])
axins5.imshow(sgrb2dsdata*fluxscalefactor,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray_r')
axins5.set_xlim((centerx3-width3),(centerx3+width3))
axins5.set_ylim((centery3-width3),(centery3+width3))
axins5.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins6=axins5.inset_axes([1.25,0,1,1])
axins6.imshow(dsiihdu.data, vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape2=dsiihdu.data.shape
xpos2=round(shape2[1]*xratio)
ypos2=round(shape2[0]*yratio)
txt2=axins6.text(xpos2,ypos2,'DS2',size=annote_size,color='black')
txt2.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins7=ax.inset_axes([righthoriz,-0.036,axins_dims,axins_dims])
axins7.imshow(sgrb2dsdata*fluxscalefactor,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray_r')
axins7.set_xlim((centerx6-width6),(centerx6+width6))
axins7.set_ylim((centery6-width6),(centery6+width6))
axins7.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins8=axins7.inset_axes([1.25,0,1,1])
axins8.imshow(dsvhdu.data, vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape5=dsvhdu.data.shape
xpos5=round(shape5[1]*xratio)
ypos5=round(shape5[0]*yratio)
txt5=axins8.text(xpos5,ypos5,'DS5',size=annote_size,color='black')
txt5.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins9=ax.inset_axes([righthoriz,0.864,axins_dims,axins_dims])
axins9.imshow(sgrb2dsdata*fluxscalefactor,origin='lower',norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt', max_cut=jymax),cmap='gray_r')
axins9.set_xlim((centerx5-width5),(centerx5+width5))
axins9.set_ylim((centery5-width5),(centery5+width5))
axins9.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins10=axins9.inset_axes([1.25,0,1,1])
axins10.imshow(dsivhdu.data, vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape4=dsivhdu.data.shape
xpos4=round(shape4[1]*xratio)
ypos4=round(shape4[0]*yratio)
txt4=axins10.text(xpos4,ypos4,'DS4',size=annote_size,color='black')
txt4.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins11=ax.inset_axes([lefthoriz,0.564,axins_dims,axins_dims])
axins11.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins11.set_xlim((centerx4-width4),(centerx4+width4))
axins11.set_ylim((centery4-width4),(centery4+width4))
axins11.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins12=axins11.inset_axes([-1.25,0,1,1])
axins12.imshow(dsiiihdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape3=dsiiihdu.data.shape
xpos3=round(shape3[1]*xratio)
ypos3=round(shape3[0]*yratio)
txt3=axins12.text(xpos3,ypos3,'DS3',size=annote_size,color='black')
txt3.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins13=ax.inset_axes([lefthoriz,0.264,axins_dims,axins_dims])
axins13.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins13.set_xlim((centerx7-width7),(centerx7+width7))
axins13.set_ylim((centery7-width7),(centery7+width7))
axins13.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins14=axins13.inset_axes([-1.25,0,1,1])
axins14.imshow(ds6hdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape6=ds6hdu.data.shape
xpos6=round(shape6[1]*xratio)
ypos6=round(shape6[0]*yratio)
txt6=axins14.text(xpos6,ypos6,'DS6',size=annote_size,color='black')
txt6.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins15=ax.inset_axes([lefthoriz,-0.036,axins_dims,axins_dims])
axins15.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins15.set_xlim((centerx8-width8),(centerx8+width8))
axins15.set_ylim((centery8-width8),(centery8+width8))
axins15.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins16=axins15.inset_axes([-1.25,0,1,1])
axins16.imshow(ds7hdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape7=ds7hdu.data.shape
xpos7=round(shape7[1]*xratio)
ypos7=round(shape7[0]*yratio)
txt7=axins16.text(xpos7,ypos7,'DS7',size=annote_size,color='black')
txt7.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins17=ax.inset_axes([lefthoriz,-0.336,axins_dims,axins_dims])
axins17.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins17.set_xlim((centerx10-width10),(centerx10+width10))
axins17.set_ylim((centery10-width10),(centery10+width10))
axins17.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins18=axins17.inset_axes([-1.25,0,1,1])
axins18.imshow(ds9hdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape9=ds9hdu.data.shape
xpos9=round(shape9[1]*xratio)
ypos9=round(shape9[0]*yratio)
txt9=axins18.text(xpos9,ypos9,'DS9',size=annote_size,color='black')
txt9.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

axins19=ax.inset_axes([righthoriz,-0.336,axins_dims,axins_dims])
axins19.imshow(sgrb2dsdata*fluxscalefactor,origin='lower', norm=visualization.simple_norm(sgrb2dsdata, stretch='sqrt',max_cut=jymax),cmap='gray_r')
axins19.set_xlim((centerx9-width9),(centerx9+width9))
axins19.set_ylim((centery9-width9),(centery9+width9))
axins19.contour((sgrb2dsdata*fluxscalefactor),levels=cntrlist,colors='red',linewidths=0.7)
axins20=axins19.inset_axes([1.25,0,1,1])
im=axins20.imshow(ds8hdu.data,vmax=tmax,vmin=tmin,origin='lower',cmap='inferno')
shape8=ds8hdu.data.shape
xpos8=round(shape8[1]*xratio)
ypos8=round(shape8[0]*yratio)
txt8=axins20.text(xpos8,ypos8,'DS8',size=annote_size,color='black')
txt8.set_path_effects([pe.withStroke(linewidth=3,foreground='w')])

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
axins17.tick_params(direction='in')
axins18.tick_params(direction='in')
axins19.tick_params(direction='in')
axins20.tick_params(direction='in')

ra.set_ticklabel_visible(True)
ra.set_axislabel('RA (J2000)')
dec.set_ticklabel_visible(True)
dec.set_axislabel('Dec (J2000)',minpad=-0.7)
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
axins17.xaxis.set_ticklabels([])
axins17.yaxis.set_ticklabels([])
axins18.xaxis.set_ticklabels([])
axins18.yaxis.set_ticklabels([])
axins19.xaxis.set_ticklabels([])
axins19.yaxis.set_ticklabels([])
axins20.xaxis.set_ticklabels([])
axins20.yaxis.set_ticklabels([])

ax.indicate_inset_zoom(axins)
ax.indicate_inset_zoom(axins3)
ax.indicate_inset_zoom(axins5)
ax.indicate_inset_zoom(axins7)
ax.indicate_inset_zoom(axins9)
ax.indicate_inset_zoom(axins11)
ax.indicate_inset_zoom(axins13)
ax.indicate_inset_zoom(axins15)
ax.indicate_inset_zoom(axins17)
ax.indicate_inset_zoom(axins19)
#axins.indicate_inset_zoom(axins2)

#plt.Circle((centerx3,centery3),15)
savefigpath='/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/zoominfigs/contfix_nograd_contours_bolocam_zoomin.png'
cbax=axins20.inset_axes([1.25,0,1,1])#inset_axes(axins2,width='5%',height='100%',loc='lower right', bbox_to_anchor=a,borderpad=0,bbox_transform=ax.transAxes)#,bbox_transform=axins20.transAxes,borderpad=0)
cb=plt.colorbar(shrink=2,pad=0, mappable=im,cax=cbax)
cb2=plt.colorbar(continuum,shrink=1.2,pad=0.5)#,anchor=(0.0,1.0))
#cb.set_label(label='T$_{rot}$ (K)',size=15)
cb2.set_label(label=r'$T_{b}$ (K)',size=15) #r'$S_{\nu}$ (mJy/beam)'
#axins2show.figure.axes[1].tick_params(axis='y',labelsize=13)
cb2.ax.tick_params(labelsize=13)#continuum.figure.axes[1].tick_params(axis='y',labelsize=13)
#plt.savefig(savefigpath,bbox_inches='tight',overwrite=True)
#plt.tight_layout()
plt.show()
