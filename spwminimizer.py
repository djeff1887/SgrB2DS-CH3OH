from spectral_cube import SpectralCube as sc
from astropy.io import fits
import astropy.units as u
import time
import glob
import os

home='/blue/adamginsburg/d.jeff/1e6-16_results/'
print('Get cubes')
cubes=glob.glob(home+'*pbcor*')
#spw1=sc.read("/blue/adamginsburg/d.jeff/imaging_results/noncontsubcubes/SgrB2DS_field1_spw1_cube.image.fits",use_dask=True)
#spw2=sc.read("/blue/adamginsburg/d.jeff/imaging_results/noncontsubcubes/SgrB2DS_field1_spw2_cube.image.fits",use_dask=True)

for cube in cubes:
    for i in range(4):
        if f'spw{i}' not in cube:
            print(f'spw{i} not in {cube}')
            continue
        else:
            outfile=home+f'SgrB2DS_field1_spw{i}_cube_robust0_niter1e6_chunks16_minimize.image.pbcor.fits'
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