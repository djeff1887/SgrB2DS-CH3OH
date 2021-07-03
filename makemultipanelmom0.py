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

source='SgrB2S'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10}
fnum=fields[source]

immode='mom0'
colormap={'mom0':'bone','tex':'inferno','numtrans':'CMRmap','nupper':'Blues_r'}
cm = copy.copy(mpl.cm.get_cmap(colormap[immode]))#mom0 bone, temperature inferno, nupper Blues_r, detections CMRmap
cm.set_bad('black')

sourcelocs={'SgrB2S':'new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/field10originals_spatialandvelocitymaskingtrial5_newexclusions3andexclusionsnotinfit/','DSii':'/field10originals_noexclusions/'}

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
mom0path=sourcepath+'mom0/*_masked.fits'

mom0images=glob.glob(mom0path)

samplefits=fits.open(mom0images[0])
samplewcs=WCS(samplefits[0])

mastereuks=[]
masterqns=[]
excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2'],'DSii':''}
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
numrows=math.ceil(len(mastereuks)/numcols)
fig,ax=plt.subplots(numrows,numcols,sharey=True)#,use_gridspec=True)
print('Number of rows: ', numrows)
firstpanel=True

gs1 = gridspec.GridSpec(numrows, numcols)
gs1.update(wspace=0.025, hspace=0.05)

for row in np.arange(numrows):
    for col in np.arange(numcols):
        transition=fits.getdata(sortedmom0[col+(row*numcols)])
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

plt.colorbar(mappable=im,pad=0.05,shrink=5)
fig.subplots_adjust(wspace=0,hspace=0)
plt.show()