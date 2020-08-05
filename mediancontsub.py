import numpy as np
from spectral_cube import SpectralCube as sc
import glob
import astropy.units as u
import matplotlib.pyplot as plt
import time

#plt.close('all')
#spw3=sc.read("/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw3_cube.image.fits")
spw0=sc.read("/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw0_cube.image.fits")

cube_w3=spw0.wcs

targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
targetpixcrd=cube_w3.all_world2pix(targetworldcrd,1,ra_dec_order=True)

testonpix=spw0[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]

testonreg=spw0[:,700:900,700:900]
testonreg.allow_huge_operations=True
testonreg=testonreg.mask_out_bad_beams(threshold=0.01)
#print(testonreg[1:(len(testonreg)-2)].beams)
print('Begin median calc...')
starttime=time.time()
testmed=testonreg.median(axis=0)
elapsed=time.time()-starttime
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
#print(f'Median calc elapsed time: {elapsed}')
testmedvalues=testmed.value
plt.imshow(testmedvalues)
plt.show()
testmedsub=testonreg-testmed



