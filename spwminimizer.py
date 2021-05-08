from spectral_cube import SpectralCube as sc
from astropy.io import fits
import astropy.units as u
import time
import glob
import os

home='/orange/adamginsburg/sgrb2/d.jeff/field10pbcors/'
print('Get cubes')
cubes=glob.glob(home+'*pbcor.fits')
fieldnum=10
out='/orange/adamginsburg/sgrb2/d.jeff/data/field10originalimages/'

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in cubes:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

for cube in datacubes:
    for i in range(4):
        if f'spw{i}' not in cube:
            print(f'spw{i} not in {cube}')
            continue
        else:
            outfile=out+f'SgrB2DS_field{fieldnum}_spw{i}_cube_minimize.image.pbcor.fits'
            if os.path.isfile(outfile):
                print(f'Outfile {outfile} already exists.')
                continue
            else:
                data=sc.read(cube)
                data.allow_huge_operations=True
                print(f'Minimizing {cube}')
                starttime=time.time()
                mincube=data.minimal_subcube()
                elapsed=time.time()-starttime
                print(f'Done in {time.strftime("%H:%M:%S", time.gmtime(elapsed))}')
                mincube.write(outfile)
'''
spw2.allow_huge_operations=True

print('Minimizing spw1')
starttime=time.time()
spw1min=spw1.minimal_subcube()
elapsed=time.time()-starttime
print(f'Done in {time.strftime("%H:%M:%S", time.gmtime(elapsed))}')
print('Minimizing spw2')
starttime=time.time()
spw2min=spw2.minimal_subcube()
elapsed=time.time()-starttime
print(f'Done in {time.strftime("%H:%M:%S", time.gmtime(elapsed))}')
spw1min.write("/blue/adamginsburg/d.jeff/imaging_results/noncontsubcubes/SgrB2DS_field1_spw1_cube_minimize.image.fits",format="fits")
spw2min.write("/blue/adamginsburg/d.jeff/imaging_results/noncontsubcubes/SgrB2DS_field1_spw2_cube_minimize.image.fits",format="fits")
'''