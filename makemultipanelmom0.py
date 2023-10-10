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

savefighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}/'+sourcelocs[source]

if os.path.exists(savefighome):
    print(f'Figure directory {savefighome} already exists.')
    pass
else:
    print(f'Making figure directory {savefighome}')
    os.makedirs(savefighome)
    
savefigpath=savefighome+f'contfix_multi{immode}fig.png'

mom0images=glob.glob(mom0path)

samplefits=fits.open(mom0images[0])
samplewcs=WCS(samplefits[0])

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
for euk in mastereuks:
    for master in range(len(mastereuks)):
        if euk == float(mastertable['$E_U$'][master]) and masterqns[master] not in sortedqns:
            sortedqns.append(masterqns[master])
            sortedoldqns.append(masteroldqns[master])
        else:
            continue
#print(sortedqns)

qnstoplot=[]
for (qn,plotqn) in zip(sortedoldqns,sortedqns):
    for mom0 in mom0images:
        if qn_replace(qn) in mom0 and qn_replace(qn) not in excludedlines[source] and qn_replace(qn) not in qnstoplot:
            sortedmom0str.append(mom0)
            sortedmom0.append(fits.getdata(mom0))
            qnstoplot.append(plotqn)
            #pdb.set_trace()
        else:
            #print(qn_replace(qn))
            continue
            
figsizes={'DSi':(11,9),'DSiii':(16,5),'DSv':(8,3),'DSVI':(12,7),'DSVII':(12,10),'DSVIII':(12,11),'DSIX':(12,11)}
if source in figsizes.keys():
    fs=figsizes[source]
else:
    fs=(12,13)
numcols=5
numrows=math.ceil(len(sortedmom0)/numcols)
fig,ax=plt.subplots(numrows,numcols,sharey=True,figsize=fs)#,use_gridspec=True)
plt.rcParams['figure.dpi'] = 150
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
sourcewspace={'SgrB2S':0,'DSi':0,'DSii':0,'DSiii':-0.59,'DSiv':0,'DSv':-0.28,'DSVI':-0.16,'DSVII':0,'DSVIII':0,'DSIX':0}
sourcehspace={'SgrB2S':-0.72,'DSi':-0.05,'DSii':-0.77,'DSiii':0,'DSiv':-0.77,
              'DSv':0,'DSVI':0,'DSVII':-0.6825,'DSVIII':-0.72,'DSIX':-0.72}
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
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=10)#12
                #print(f'{sortedmom0str[i]},{qnstoplot[i]}')
                firstpanel=False
                i+=1
            else:
                ax[row,col].imshow(transition,vmax=fmax,vmin=fmin,origin='lower',cmap=cm)
                ax[row,col].set_xticklabels([])
                ax[row,col].set_yticklabels([])
                ax[row,col].tick_params(direction='in')
                ax[row,col].set_xticks([])
                ax[row,col].set_yticks([])
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=10)#12
                #print(f'{sortedmom0str[i]}, {qnstoplot[i]}')
                i+=1
                #if row == (numrows-1) and col == (numcols-1):

heights={'SgrB2S':'293.5%','DSi':'383%','DSii':'200.5%','DSiii':'200%','DSiv':'200.5%','DSv':'200%','DSVI':'300%','DSVII':'200%','DSVIII':'199.5%','DSIX':'199.5%'}
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
plt.savefig(savefigpath,dpi=300)
#plt.tight_layout()
#gs1.tight_layout(fig,pad=1)
plt.show()
