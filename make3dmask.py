import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import matplotlib as mpl
import pickle as pkl

mpl.interactive(True)

images=['spw0','spw1','spw2','spw3']

inpkl=open("/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/testbox2dict.obj","rb")
indict=pkl.load(inpkl)

cubefn="/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage/spw0minimize.image.pbcor_line.fits"

for img in images:
    if img in cubefn
        spw=img
        break
    else:
        pass

transition_keys=list(indict[spw].keys())
std=sc.read("/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/errorimgs/{spw}"#indict[spw][transition_keys[0]]['stddev']

cube=sc.read(cubefn)
brightest_cube=sc.read("/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/3dmasktesting/stackedlines1-3-4-6-11.fits")

boolmask=cube > (std)

testboolmask2=brightest_cube > 8*u.K

testmasked2=brightest_cube.with_mask(testboolmask2)

testmasked2.write('/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/z0_0002306756533745274_testbox2_5-6mhzwidth/3dmasktesting/test8Kmask.fits')



'''
# compute various moments & statistics along the spcetral dimension
peak_velocity = brightest_cube.spectral_axis[brightest_cube.argmax(axis=0)]
max_map = peak_amplitude = brightest_cube.max(axis=0)
width_map = brightest_cube.linewidth_sigma() # or vcube.moment2(axis=0)**0.5
fwhm_map = brightest_cube.linewidth_fwhm() # FOR TESTING
sqrtmom2_map = brightest_cube.moment2(axis=0)**0.5 # FOR TESTING
centroid_map = brightest_cube.moment1(axis=0)
'''