import glob
from spectral_cube import SpectralCube as sc
import os
import astropy.units as u
from astropy.io import fits
import pdb

inpath='/orange/adamginsburg/sgrb2/d.jeff/data/field10originalimages/'
#inpath='/blue/adamginsburg/d.jeff/imaging_results/data/OctReimage/'
beamcubes=glob.glob(inpath+'*.fits')
home='/orange/adamginsburg/sgrb2/d.jeff/products/field10originalimages/'
#home='/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/'
cubes=glob.glob(home+'*pbcor_line.fits')
region='fk5; box(266.8316387, -28.3971867, 0.0010556, 0.0010556)'#DSi-large
#region='fk5; box(266.8353410,-28.3962005,0.0016806,0.0016806)'#box(266.8333438, -28.3966103, 0.0014028, 0.0014028)'#box(266.8315833, -28.3971867, 0.0006528, 0.0006528)' 
outpath='/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSi/field10originals/'#imaging_results/DSii_iiibox1/'
statfixpath='/blue/adamginsburg/d.jeff/SgrB2DSstatcontfix/field10originals/'

cubestobox=[]

images=['spw0','spw1','spw2','spw3']

orderedcubes=[]
orderedbeamcubes=[]

for spew in images:
    for f1,f2 in zip(cubes,beamcubes):
        if spew in f1:
            orderedcubes.append(f1)
        if spew in f2:
            orderedbeamcubes.append(f2)
    #images.append(files[77:81])#[57:61])
    
assert 'spw0' in orderedcubes[0] and 'spw0' in orderedbeamcubes[0], f'Cube list out of order'

print('Cube lists successfully reordered')

if not os.path.isdir(outpath):
    print(f'Creating filepath {outpath}')
    os.makedirs(outpath)
else:
    print(f'{outpath} already exists. Proceeding...\n')
    
if 'products' in home:
    print('STATCONT products detected\n')
    if not os.path.isdir(statfixpath):
        print(f'Creating beamfix directory {statfixpath}')
        os.makedirs(statfixpath)
    else:
        print(f'{statfixpath} already exists. Proceeding...')
    for beamcube,statcube in zip(orderedbeamcubes,orderedcubes):
        print(f'Extracting beams from {beamcube}')
        beamfits=fits.open(beamcube)
        cubefits=fits.open(statcube)
        beams=beamfits[1]
        cubedata=cubefits[0]
        newhdulist=fits.HDUList([cubedata,beams])
        print(f'Beamlist merged with Primary HDU in {statcube}')
        cubewithbeampath=statcube.replace(home,statfixpath)
        print(f'Saving new fits file {cubewithbeampath}\n')
        newhdulist.writeto(cubewithbeampath)
        cubestobox.append(cubewithbeampath)
else:
    cubestobox=cubes
        
        
for spw, cube in zip(images, cubestobox):
    boxcubename=outpath+spw+'minimize.image.pbcor_line.fits'
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
    