import numpy as np
from astropy.io import fits

source='DSv'
fnum=10
detectnum=5

home=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/'
sourcepath=home+f'{source}/200K_field10originals_z0_000190713_5-6mhzwidth_stdfixes/'
print(f'Accessing images in {sourcepath}')
transitionimg=fits.open(sourcepath+"texmap_3transmask_3sigma_allspw_withnans_weighted.fits")[0]
teximg=fits.open(sourcepath+"texmap_3sigma_allspw_withnans_weighted.fits")[0]

transitionmap=transitionimg.data
texmap=teximg.data

print('Masking temperature map')
maskedarr=np.ma.masked_where(transitionmap<detectnum,texmap)
maskedmap=maskedarr.filled(fill_value=np.nan)

print('Creating new fits file')
maskedhdu=fits.PrimaryHDU(maskedmap)
maskedhdu.header=teximg.header
maskedfits=fits.HDUList(maskedhdu)
maskedoutpath=sourcepath+f"texmap_{detectnum}transmask_3sigma_allspw_withnans_weighted.fits"
print(f'Saving at {maskedoutpath}')
maskedfits.writeto(maskedoutpath)
print('Done')

