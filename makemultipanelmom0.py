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

source='SgrB2S'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10}
romannums={'DSi':'DSI','DSii':'DSii'}
fnum=fields[source]

immode='mom0'
colormap={'mom0':'bone','tex':'inferno','numtrans':'CMRmap','nupper':'Blues_r'}
cm = copy.copy(mpl.cm.get_cmap(colormap[immode]))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap
cm.set_bad('black')

sourcelocs={'SgrB2S':'new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/Kfield10originals_noexclusions/','DSiii':'/Kfield10originals_noexclusions/','DSiv':'/Kfield10originals_noexclusions/','DSv':f'/Kfield10originals_noexclusions_include4-3_K_trial1/','DSVI':'/Kfield2originals_trial3_8_6-8_7excluded/','DSVII':f'/Kfield3originals_trial1_noexclusions/','DSVIII':f'/Kfield3originals_150K_trial1_noexclusions/'}

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
mom0path=sourcepath+'mom0/*_masked.fits'

if source == 'DSi':
    savefighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{romannums[source]}'+sourcelocs[source]
else:
    savefighome=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}'+sourcelocs[source]

if os.path.exists(savefighome):
    print(f'Figure directory {savefighome} already exists.')
    pass
else:
    print(f'Making figure directory {savefighome}')
    os.makedirs(savefighome)
    
savefigpath=savefighome+f'firstmulti{immode}fig.png'

mom0images=glob.glob(mom0path)

samplefits=fits.open(mom0images[0])
samplewcs=WCS(samplefits[0])

mastereuks=[]
masterqns=[]
excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','15_6-15_7E1vt1'],'DSii':'','DSiii':'','DSiv':'','DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1'],'DSVII':'','DSVIII':''}
masterlist=np.genfromtxt(sourcepath+'mastereuksqnsfreqsdegens.txt',dtype=str)
for master in range(len(masterlist[:,0])):
    mastereuks.append(float(masterlist[master,0]))
    masterqns.append(masterlist[master,1])
mastereuks.sort()

sortedqns=[]
sortedmom0=[]
for euk in mastereuks:
    for master in range(len(mastereuks)):
        if euk == float(masterlist[master,0]):
            sortedqns.append(masterqns[master])
        else:
            continue
#print(sortedqns)

for qn in sortedqns:
    for mom0 in mom0images:
        if qn_replace(qn) in mom0:
            sortedmom0.append(mom0)
        else:
            continue
        
numcols=5
numrows=math.ceil(len(sortedqns)/numcols)
fig,ax=plt.subplots(numrows,numcols,sharey=True,figsize=(14,13))#,use_gridspec=True)
plt.rcParams['figure.dpi'] = 150
print(f'Number of rows: {numrows}')
print(f'Number of columns: {numcols}')
firstpanel=True

gs1 = gridspec.GridSpec(numrows, numcols)
gs1.update(wspace=0.025, hspace=0.05)

brightestline={'SgrB2S':600,'DSi':450,'DSii':0,'DSiii':0,'DSiv':0,'DSv':0,'DSVI':0,'DSVII':0,'DSVIII':0}

for row in np.arange(numrows):
    for col in np.arange(numcols):
        transition=fits.getdata(sortedmom0[(col+(row*numcols))-1])
        whichline=(col+(row*numcols))
        if firstpanel:
            im=ax[row,col].imshow(transition,origin='lower',cmap=cm)
            fmax=np.nanmax(transition)
            fmin=np.nanmin(transition)
            ax[row,col].set_xticklabels([])
            ax[row,col].set_yticklabels([])
            ax[row,col].tick_params(direction='in')
            firstpanel=False
        else:
            ax[row,col].imshow(transition,vmax=fmax,vmin=fmin,origin='lower',cmap=cm)
            ax[row,col].set_xticklabels([])
            ax[row,col].set_yticklabels([])
            ax[row,col].tick_params(direction='in')
        #if row == (numrows-1) and col == (numcols-1):

cax=plt.axes([0.78, 0.11, 0.05, 0.77])
plt.colorbar(mappable=im,shrink=5,cax=cax,label=r'I$_{\nu}$ (K km s$^{-1}$)')              
fig.subplots_adjust(wspace=-0.75,hspace=0)
plt.savefig(savefigpath,dpi=150)
plt.show()