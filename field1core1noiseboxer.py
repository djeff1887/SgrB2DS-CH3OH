import glob
from spectral_cube import SpectralCube as sc
import os
import astropy.units as u

home='/blue/adamginsburg/d.jeff/imaging_results/'
imglist=glob.glob(home+'*.fits')
region='fk5; box(266.8351754,-28.3961033,0.0008333,0.0008333)'#'/blue/adamginsburg/d.jeff/imaging_results/field1core1box.reg'
filepath='/blue/adamginsburg/d.jeff/imaging_results/field1core1statcontnoiseimages/'

images=[]#'spw0','spw2','spw1','spw3']
for files in cubes:
    images.append(files[57:61])
assert 'spw0' in images, f'image name list does not match spw# format'

for spw, image in zip(images, imglist):
    boxcubename=filepath+spw+'minimize_hasnoncontbeams.image_line.fits'
    if os.path.isfile(boxcubename):
        print(f'{boxcubename} already exists. Skipping...\n')
        continue
    else:
        print(f'Grabbing datacube from {cube}')
        fullsizecube=fits.getdata(image)*u.Jy/u.beam
        print('Creating subcube and converting from Jy/beam to K')
        boxedsubcubeK=#fullsizecube.subcube_from_ds9region(region).to(u.K)
        #print('Converting spectral axis to km/s')
        #boxedsubcubeKkms=boxedsubcubeK.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
        print(f'Saving to {boxcubename}')
        boxedsubcubeK.write(boxcubename,format='fits',overwrite=True)
        print('Finished\n')
    