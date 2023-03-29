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
from utilities import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.table import Table

mpl.interactive(True)

source='DSiii'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
romannums={'DSi':'DSI','DSii':'DSii'}
fnum=fields[source]

immode='mom0'
colormap={'mom0':'bone_r','tex':'inferno','numtrans':'CMRmap','nupper':'Blues_r'}
cm = copy.copy(mpl.cm.get_cmap(colormap[immode]))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap
cm.set_bad('white')

sourcelocs={'SgrB2S':'/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/','DSi':'/nov2022continuumsanitycheck/','DSii':'/nov2022continuumsanitycheck/','DSiii':'/nov2022continuumsanitycheck/','DSiv':'/nov2022contniuumsanitycheck/','DSv':f'/nov2022contniuumsanitycheck/','DSVI':'/nov2022continuumsanitycheck/','DSVII':f'/nov2022contniuumsanitycheck/','DSVIII':f'/nov2022contniuumsanitycheck/','DSIX':f'/nov2022contniuumsanitycheck/'}#{'SgrB2S':'new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'Kfield10originals_noexclusions/','DSiii':'Kfield10originals_noexclusions/','DSiv':'Kfield10originals_noexclusions/','DSv':f'Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':'Kfield2originals_trial3_8_6-8_7excluded/','DSVII':f'Kfield3originals_trial1_noexclusions/','DSVIII':f'Kfield3originals_175K_trial1_noexclusions/','DSIX':'Kfield7originals_150K_trial1_noexclusions/'}

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
mom0path=sourcepath+'mom0/*_masked.fits'

if source == 'DSi':
    savefighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{romannums[source]}/'+sourcelocs[source]
else:
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
excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1','15_6-15_7E1vt1','9_6-9_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','8_6-8_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','15_6-15_7E1vt1'],'DSii':'','DSiii':'','DSiv':'','DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1'],'DSVII':'','DSVIII':'','DSIX':'','DSX':''}
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
            
figsizes={'DSi':(11,9),'DSiii':(16,5)}
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
sourcewspace={'SgrB2S':0,'DSi':0,'DSii':0,'DSiii':-0.61,'DSiv':0,'DSv':0,'DSVI':0,'DSVII':0,'DSVIII':0,'DSIX':0}
sourcehspace={'SgrB2S':-0.72,'DSi':0,'DSii':-0.77,'DSiii':0,'DSiv':-0.72,
              'DSv':-0.67,'DSVI':-0.602,'DSVII':-0.72,'DSVIII':-0.72,'DSIX':-0.72}
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
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=12)
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
                ax[row,col].annotate(f'{qnstoplot[i]}',xy=(0,0),xytext=((imgdims[0]/20),(imgdims[1]/1.2)),fontsize=12)
                #print(f'{sortedmom0str[i]}, {qnstoplot[i]}')
                i+=1
                #if row == (numrows-1) and col == (numcols-1):

heights={'SgrB2S':'290%','DSi':'400%','DSii':'200%','DSiii':'100%','DSiv':'225%','DSv':'70%','DSVI':'132%','DSVII':'225%','DSVIII':'220%'}
bboxparams={'SgrB2S':(1.05,0,1,1),'DSi':(1.05,0,1,1),'DSiii':(1.05,0.26,1,1),'DSv':(1.05,0.34,1,1),'DSVI':(1.05,0.25,1,1),'DSVII':(1.05,-0.1,1,1)}
if source in heights:
    strheight=heights[source]
else:
    strheight='225%'
    
if source in bboxparams.keys():
    a=bboxparams[source]#
    lowercorner=ax[(numrows-1),(numcols-1)]#[1,4]
else:
    a=(1.05,0,1,1)#(1.05,-0.1,1,1)
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
plt.savefig(savefigpath,dpi=150)
#plt.tight_layout()
#gs1.tight_layout(fig,pad=1)
plt.show()
