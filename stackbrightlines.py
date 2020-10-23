import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import matplotlib as mpl

mpl.interactive(True)

files=glob.glob("/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/spectralslabs/km_s/*.fits")
sum=sc.read(files[1])
brightlines=[1,3,4,6,11]

for line in brightlines:
    if files[line] == files[1]:
        print(f'{files[line]}=={files[1]}')
        continue
    else:
        print(f'Loading {files[line]}')
        slab=sc.read(files[line])
        plt.plot(slab.spectral_axis,slab[:,75,60].value,color='gray',drawstyle='steps')
        sum=sum+slab
    
sum=sum/len(brightlines)

plt.plot(sum.spectral_axis,sum[:,75,60].value,color='red',drawstyle='steps')



