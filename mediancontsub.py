import numpy as np
from spectral_cube import SpectralCube as sc
import glob
import astropy.units as u

spw3=sc.read("/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw3_cube.image.fits")
#spw0=sc.read("/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw0_cube.image.fits")

cube_w=spw3.wcs

targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)

testonpix=spw3[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]

testonreg=spw3[:,762:772,496:506]
testmed=testonreg.median(axis=0)
print(testmed)
