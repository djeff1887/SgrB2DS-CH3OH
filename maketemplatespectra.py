import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import regions
import math

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

source='SgrB2S'
z=0.0002306756533745274

inpath=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/OctReimage/'
path=glob.glob(inpath+'*spw3*')

home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/mastereuksqnsfreqs.txt"

detections=np.loadtxt(home,dtype=str)
detshape=np.shape(detections)

cube=sc.read(path[0])

cubewcs=cube.wcs
targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]] #SgrB2S
#[[0,0,0],[266.8332569, -28.3969, 0]] #DSii/iii
targetpixcrd=cubewcs.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    
pixxcrd,pixycrd=int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1]))
print(f'x: {pixxcrd}/y: {pixycrd}')
assert pixxcrd >= 0 and pixycrd >= 0, 'Negative pixel coords'
spw=cube[:,pixycrd,pixxcrd]

freqs=cube.spectral_axis#Hz
freqflip=False
if freqs[1] < freqs[0]:
    freqs=freqs[::-1]
    freqflip=True
    print('Corrected decreasing frequency axis')
else:
    pass

freq_min=freqs[0]
#print(freq_max)
freq_max=freqs[(len(freqs)-1)]

assert freq_max > freq_min, 'Inverted spectral axis'
print('Passed increasing spectral axis check')

print('Plotting spectra')

xlims=(232707306565.9205,232919343267.9113)
plt.plot(freqs,spw.value,drawstyle='steps',color='black')
plt.xlim(xlims[0],xlims[1])
plt.ylim(ymax=100)

firstline=True
print('Begin detected line vlines')
for row in range(detshape[0]):
    line=float(detections[row,2])
    if line >= xlims[0] and line <= xlims[1]:
        if firstline:
            plt.axvline(x=float(detections[row,2]),linestyle='--',color='green', label='CH$_3$OH')
            firstline=False
        elif not firstline:
            plt.axvline(x=float(detections[row,2]),linestyle='--',color='green')
    else:
        print(f'Line {detections[row,1]} at {(line*u.Hz).to("GHz")} skipped; out of frequency range')


plt.xlabel(r'$\nu$ (Hz)')
plt.ylabel('T$_b$ (K)')
plt.rcParams['figure.dpi'] = 150
plt.legend()
plt.show()