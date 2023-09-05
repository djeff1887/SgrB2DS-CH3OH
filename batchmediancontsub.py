import numpy as np
from spectral_cube import SpectralCube as sc
import glob
import astropy.units as u
import matplotlib.pyplot as plt
import time
import tempfile

home="/blue/adamginsburg/d.jeff/imaging_results/usedasktests/"
#plt.close('all')
print('Loading cube')
starttime=time.time()
#spw3=sc.read("/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw3_cube.image.fits")
spw0=sc.read("/blue/adamginsburg/d.jeff/imaging_results/noncontsubcubes/SgrB2DS_field1_spw0_cube.image.fits",use_dask=True)
elapsed=time.time()-starttime
print('Cube loaded in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))

tempfile.tempdir='/blue/adamginsburg/d.jeff/tmp'
print('Rechunking cube, save_to_tmp_dir=True, default chunks')
starttime=time.time()
spw0=spw0.rechunk(save_to_tmp_dir=True)
elapsed=time.time()-starttime
print(f'Rechunking done in {time.strftime("%H:%M:%S", time.gmtime(elapsed))}')
'''
cube_w3=spw3.wcs

targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)


testonpix=spw3[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]

testonreg=spw3[1:(len(spw3)-2),400:900,400:900]
'''
#spw3.allow_huge_operations=True
#spw3medsub=spw3[1:(len(spw3)-2)]
spw0.allow_huge_operations=True
print('Begin bad beam masking')
starttime=time.time()
spw0=spw0.mask_out_bad_beams(threshold=0.01)
elapsed=time.time()-starttime
print('Bad beam done in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
#print(testonreg[1:(len(testonreg)-2)].beams)

print('Begin spw0 median calc...')
starttime=time.time()
testmed0=spw0.median(axis=0)
elapsed=time.time()-starttime
print('Done in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
'''
print('Begin spw3 median calc...')
starttime=time.time()
testmed3=spw3medsub.median(axis=0)
elapsed=time.time()-starttime
print('Done in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
'''
print('Writing median cont image')
starttime=time.time()
testmed0.write(home+'field1spw0mediancont_fullsize_usedask_rechunkdefault.fits',format='fits',overwrite=True)
elapsed=time.time()-starttime
print('median cont image write done in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
'''
testmed3.write('/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/field1spw3mediancont.fits',format='fits')
'''
print('Performing median subtraction')
starttime=time.time()
medsub0=spw0-testmed0
elapsed=time.time()-starttime
print('medsub done in')
print(time.strftime("%H:%M:%S", time.gmtime(elapsed)))
print('Writing medcontsub cube')
medsub0.write(home+'SgrB2DS_field1_spw0_cube_fullsize_usedask_rechunkdefault.image.fits',format='fits')
print('Done')
'''
spw3medsub=spw3-testmed3
spw3medsub.write('/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw3_cube_medsub.image.fits',format='fits')
print('Finished')
'''


