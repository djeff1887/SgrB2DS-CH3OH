import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.ticker as mtick
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import regions
import math
import matplotlib as mpl
import pdb
from utilities import *
import copy
from astropy.wcs import WCS
import matplotlib.gridspec as gridspec
import os
#from utilities import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.table import Table

mpl.interactive(True)

#plt.rcParams['figure.dpi']=300

source='DSi'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
romannums={'DSi':'DSI','DSii':'DSii'}
fnum=fields[source]

immode='mom0'
colormap={'mom0':'bone_r','tex':'inferno','numtrans':'CMRmap','nupper':'Blues_r','nh2':'cividis'}
cm = copy.copy(mpl.cm.get_cmap(colormap[immode]))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap
cm.set_bad('white')

sourcelocs={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
mom0path=sourcepath+'mom0/*_masked.fits'
cntrfile='/orange/adamginsburg/sgrb2/2017.1.00114.S/imaging_results/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter_feathered_with_bolocam.fits'

savefighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'+sourcelocs[source]

if os.path.exists(savefighome):
    print(f'Figure directory {savefighome} already exists.')
    pass
else:
    print(f'Making figure directory {savefighome}')
    os.makedirs(savefighome)
    
savefigpath=savefighome+f'pacman_multi{immode}fig.png'

#Gathers a core's mom0 images from the relevant data version
mom0images=glob.glob(mom0path)

samplefits=fits.open(mom0images[0])
samplewcs=WCS(samplefits[0])

cntrhdu=fits.open(cntrfile)[0]
cntrrms=0.0002#Jy
cntrlist=cntrrms*np.array([5,9,27,81,128,281])
cntrdata=np.squeeze(cntrhdu.data)
hduwcs=WCS(samplefits[0])#of sample line data
cntrwcs=WCS(cntrhdu).celestial
hdubeam=radio_beam.Beam.from_fits_header(samplefits[0].header)#sample line beam

#These are the xy positions of hot cores in the 1mm continuum data, along with the dimensions of their images
corecntmpositions={'SgrB2S':(1453,3317,122/2),'DSi':(1687,3246,74/2),'DSii':(1567,3309,46/2),'DSiii':(1583,3265,50/2),
                   'DSiv':(1640,3371,65/2),'DSv':(1655,3211,40/2),'DSVI':(1282,2677,125/2),
                   'DSVII':(992,2363,150/2),'DSVIII':(1040,2191,101/2),'DSIX':(667,800,69/2)}
corecntmx=corecntmpositions[source][0]
corecntmy=corecntmpositions[source][1]
corecntmhalf=int(corecntmpositions[source][2])
cntmcutout=cntrdata[corecntmy-corecntmhalf:corecntmy+corecntmhalf,corecntmx-corecntmhalf:corecntmx+corecntmhalf]

mastereuks=[]
masterqns=[]
masteroldqns=[]
mastertable=Table.read('multiplot_methanoltransitiontable.fits')#np.genfromtxt(sourcepath+'mastereuksqnsfreqsdegens.txt',dtype=str)
for master in mastertable:
    mastereuks.append(float(master['$E_U$']))
    masterqns.append(master['Transition'])
    masteroldqns.append(master['OldTransition'])
'''
for master in range(len(masterlist[:,0])):
    mastereuks.append(float(masterlist[master,0]))
    if masterlist[master,1] not in excludedlines[source]:
        masterqns.append(masterlist[master,1])
'''
mastereuks.sort()

sortedqns=[]
sortedoldqns=[]
sortedmom0str=[]
sortedmom0=[]
sortedeuks=[]
for euk in mastereuks:
    for master in range(len(mastereuks)):
        if euk == float(mastertable['$E_U$'][master]) and masterqns[master] not in sortedqns:
            sortedqns.append(masterqns[master])
            sortedoldqns.append(masteroldqns[master])
            sortedeuks.append(euk)
        else:
            continue
#print(sortedqns)

qnstoplot=[]
eukstoplot=[]
for (qn,plotqn,eupper) in zip(sortedoldqns,sortedqns,sortedeuks):
    assert len(sortedqns)==len(sortedeuks), 'QN and Eupper list lengths don\'t match'
    for mom0 in mom0images:
        if qn_replace(qn) in mom0 and qn_replace(qn) not in excludedlines[source] and qn_replace(qn) not in qnstoplot:
            sortedmom0str.append(mom0)
            sortedmom0.append(fits.getdata(mom0))
            qnstoplot.append(plotqn)
            eukstoplot.append(int(eupper))
        else:
            #print(qn_replace(qn))
            continue

#pdb.set_trace()

figsizes={'DSi':(8,4),'DSii':(12,6),'DSiii':(11,5),'DSiv':(10,6),'DSv':(8,3),'DSVI':(12,7),'DSVII':(12,10),'DSVIII':(10,6),'DSIX':(12,11)}
if source in figsizes.keys():
    fs=figsizes[source]
else:
    fs=(12,13)
numcols=5
numrows=math.ceil(len(sortedmom0)/numcols)
fig,ax=plt.subplots(numrows,numcols,sharey=True,figsize=fs)#,use_gridspec=True)
plt.rcParams['figure.dpi'] = 300
print(f'Number of rows: {numrows}')
print(f'Number of columns: {numcols}')
maxtb=np.nanmax(sortedmom0)
brightestimage=int(np.where(sortedmom0==maxtb)[0])
mintb=np.nanmin(sortedmom0)
brightestline=sortedmom0str[brightestimage]
firstpanel=True

imgdims=np.shape(sortedmom0[0])

#gs1 = gridspec.GridSpec(numrows, numcols)
#gs1.update(wspace=0.025, hspace=0.05)

#brightestline={'SgrB2S':600,'DSi':450,'DSii':0,'DSiii':0,'DSiv':0,'DSv':0,'DSVI':0,'DSVII':0,'DSVIII':0}
sourcewspace={'SgrB2S':0,'DSi':-0.54,'DSii':-0.54,'DSiii':-0.64,'DSiv':-0.12,'DSv':-0.32,'DSVI':-0.22,'DSVII':0,'DSVIII':-0.12,'DSIX':0}
sourcehspace={'SgrB2S':-0.72,'DSi':0.0,'DSii':0.0,'DSiii':0,'DSiv':0.0,
              'DSv':0,'DSVI':0,'DSVII':-0.6825,'DSVIII':0.0,'DSIX':-0.72}#0.685
fontsizes={'DSi':7,'DSiii':9}
eupperfonts={'DSi':6}
if source in fontsizes:
    fs=fontsizes[source]
else:
    fs=10
if source in eupperfonts:
    fs2=eupperfonts[source]
else:
    fs2=fontsizes[source]
i=0
for row in np.arange(numrows):
    for col in np.arange(numcols):
        if i >= len(sortedmom0):
            ax[row,col].set_axis_off()
            #break
        else:
            transition=sortedmom0[i]
            if firstpanel:
                im=ax[row,col].imshow(transition,origin='lower',cmap=cm)
                fmax=maxtb#np.nanmax(transition)
                fmin=mintb#np.nanmin(transition)
                ax[row,col].set_xticklabels([])
                ax[row,col].set_yticklabels([])
                ax[row,col].tick_params(direction='in')
                ax[row,col].set_xticks([])
                ax[row,col].set_yticks([])
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=fs)#12
                ax[row,col].annotate(fr'$E_u$={eukstoplot[i]} K',xy=(0,0),xytext=((imgdims[0]/1.8),(imgdims[1]/20)),fontsize=fs2)
                
                ax[row,col].contour(cntmcutout, levels=cntrlist,
                                    colors='red',linewidths=0.5,zorder=1)
                firstpanel=False
                i+=1
            elif row==(numrows-1) and col==0:
                ax[row,col].imshow(transition,vmax=fmax,vmin=fmin,origin='lower',cmap=cm)
                ax[row,col].set_xticklabels([])
                ax[row,col].set_yticklabels([])
                ax[row,col].tick_params(direction='in')
                ax[row,col].set_xticks([])
                ax[row,col].set_yticks([])
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=fs)#12
                ax[row,col].annotate(fr'$E_u$={eukstoplot[i]} K',xy=(0,0),xytext=((imgdims[0]/1.8),(imgdims[1]/20)),fontsize=fs2)#12
                i+=1
            else:
                ax[row,col].imshow(transition,vmax=fmax,vmin=fmin,origin='lower',cmap=cm)
                ax[row,col].set_xticklabels([])
                ax[row,col].set_yticklabels([])
                ax[row,col].tick_params(direction='in')
                ax[row,col].set_xticks([])
                ax[row,col].set_yticks([])
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=fs)#12
                ax[row,col].annotate(fr'$E_u$={eukstoplot[i]} K',xy=(0,0),xytext=((imgdims[0]/1.8),(imgdims[1]/20)),fontsize=fs2)
                i+=1

heights={'SgrB2S':'290%','DSi':'300%','DSii':'300%','DSiii':'300%','DSiv':'300%','DSv':'200%','DSVI':'300%','DSVII':'198%','DSVIII':'300%','DSIX':'198%'}
bboxparams={'SgrB2S':(1,0,1,1),'DSi':(1,0,1,1),'DSiii':(0.69,0,1,1),'DSiv':(1,0,1,1),'DSv':(0.845,0,1,1),'DSVI':(0.9,0,1,1),'DSVII':(1,0,1,1)}
if source in heights:
    strheight=heights[source]
else:
    strheight='225%'
    
if source in bboxparams.keys():
    a=bboxparams[source]#
    lowercorner=ax[(numrows-1),(numcols-1)]#[1,4]
else:
    a=(1,0,1,1)#(1.05,-0.1,1,1)
    lowercorner=ax[(numrows-1),(numcols-1)]
axins=inset_axes(lowercorner,width='5%',height=strheight,loc='lower left', bbox_to_anchor=a,bbox_transform=lowercorner.transAxes,borderpad=0)#div=make_axes_locatable(ax[(numrows-1),(numcols-1)])
#cax=div.append_axes("right",size='10%',pad=0.1)#plt.axes([0.78, 0.11, 0.05, 0.77])
plt.colorbar(mappable=im,cax=axins,label=r'$I_{\nu}$ (K km s$^{-1}$)')#shrink=5             
fig.subplots_adjust(wspace=sourcewspace[source],hspace=sourcehspace[source])
'''
for ax in fig.get_axes():
    ss = ax.get_subplotspec()
    ax.spines.top.set_visible(ss.is_first_row())
    ax.spines.bottom.set_visible(ss.is_last_row())
    ax.spines.left.set_visible(ss.is_first_col())
    ax.spines.right.set_visible(ss.is_last_col())
'''
#plt.savefig(savefigpath,dpi=300)
#plt.tight_layout()
#gs1.tight_layout(fig,pad=1)
plt.show()
