import glob
from spectral_cube import SpectralCube as sc
import os
import astropy.units as u
from astropy.io import fits
import pdb
import numpy as np

fieldnum=os.getenv('FIELDNUM')

inpath=f"/orange/adamginsburg/sgrb2/d.jeff/field{fieldnum}/"

print(f'Grabbing cubes from {inpath}')
cubes=glob.glob(inpath+f"SgrB2DS_field{fieldnum}_spw*_cube.image.pbcor.fits")

images=['spw0','spw1','spw2','spw3']

orderedcubes=[]

print('Begin cube ordering procedure')
for spew in images:
    for f1 in cubes:
        if spew in f1:
            orderedcubes.append(f1)
    #images.append(files[77:81])#[57:61])
    
assert 'spw0' in orderedcubes[0], f'Cube list out of order'

print('Cubes in numerical order')

outpath=f'/orange/adamginsburg/sgrb2/d.jeff/data/field{fieldnum}originals_K/'

if os.path.exists(outpath):
    pass
else:
    print(f'Creating outpath {outpath}')
    os.mkdir(outpath)


for spw, cube in zip(images, orderedcubes):
    kcubename=outpath+spw+'minimize.image.pbcor.fits'
    if os.path.isfile(kcubename):
        print(f'{kcubename} already exists. Skipping...\n')
        continue
    else:
        print(f'Grabbing datacube from {cube}')
        fullsizecube=sc.read(cube)#,use_dask=True)
        fullsizecube.allow_huge_operations=True
        print('Minimizing Jy/beam cube')
        cubemin=fullsizecube.minimal_subcube()
        print('Converting from Jy/beam to K')
        cubeKmin=cubemin.to(u.K)
        print(f'Saving to {kcubename}')
        #hdulist.writeto(kcubename,overwrite=True)
        cubeKmin.write(kcubename,format='fits',overwrite=True)
        print('Finished\n')