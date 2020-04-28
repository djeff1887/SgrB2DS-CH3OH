import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import regions

home='/home/d.jeff/SgrB2DS_field1/VelocityMoms3/'#Make sure to include slash after path
fname='/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'
cube=sc.read(fname)
header=fits.getheader(fname)

ds9region=regions.read_ds9('ufrc/adamginsburg/d.jeff/imaging_results/383-649.reg')

subcube=cube.subcube_from_regions(ds9region)
print(subcube.spectral_axis)


