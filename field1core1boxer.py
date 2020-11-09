import glob
from spectral_cube import SpectralCube as sc
import os
import astropy.units as u

home='/blue/adamginsburg/d.jeff/imaging_results/'
cubes=glob.glob(home+'*.fits')
region='fk5; box(266.8333438, -28.3966103, 0.0014028, 0.0014028)'#box(266.8315833, -28.3971867, 0.0006528, 0.0006528)' #box(266.8353410,-28.3962005,0.0016806,0.0016806)'
filepath='/blue/adamginsburg/d.jeff/imaging_results/DSii_iiibox1/'

images=[]#'spw0','spw2','spw1','spw3']
for files in cubes:
    images.append(files[57:61])
assert 'spw0' in images, f'image name list does not match spw# format'

for spw, cube in zip(images, cubes):
    boxcubename=filepath+spw+'minimize_hasnoncontbeams.image_line.fits'
    if os.path.isfile(boxcubename):
        print(f'{boxcubename} already exists. Skipping...\n')
        continue
    else:
        print(f'Grabbing datacube from {cube}')
        fullsizecube=sc.read(cube,use_dask=True)
        spwrestfreq=fullsizecube.header['RESTFRQ']*u.Hz
        print('Creating subcube and converting from Jy/beam to K')
        boxedsubcubeK=fullsizecube.subcube_from_ds9region(region).to(u.K)
        #print('Converting spectral axis to km/s')
        #boxedsubcubeKkms=boxedsubcubeK.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
        print(f'Saving to {boxcubename}')
        boxedsubcubeK.write(boxcubename,format='fits',overwrite=True)
        print('Finished\n')
    