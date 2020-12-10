import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import matplotlib as mpl
import os

mpl.interactive(True)

home="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/"
slabpath=home+"spectralslabs/km_s/*.fits"
stackeddir=home+"3dmasktesting/"

files=glob.glob(slabpath)
'''Initialize the sum'''
sum=sc.read(files[1])
brightlines=[1,2,3,5,9]

'''If line is the same as the initialized line, skip. Otherwise, add'''
for line in brightlines:
    if files[line] == files[1]:
        print(f'{files[line]}=={files[1]}')
        continue
    else:
        print(f'Loading {files[line]}')
        slab=sc.read(files[line])
        plt.plot(slab.spectral_axis,slab[:,75,60].value,color='gray',drawstyle='steps')
        sum=sum+slab
    
'''Normalize by number of lines'''
sum=sum/len(brightlines)

if not os.path.isdir(stackeddir):
    print(f'Creating 3D mask testing directory {stackeddir}')
    os.mkdir(stackeddir)
else:
    print(f'3D mask testing directory {stackeddir} already exists')
    pass
    
sum.write(stackeddir+'stackedlines1-3-4-6-11.fits')

plt.plot(sum.spectral_axis,sum[:,75,60].value,color='red',drawstyle='steps')



