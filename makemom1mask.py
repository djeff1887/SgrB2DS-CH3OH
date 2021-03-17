import numpy as np
from spectral_cube import SpectralCube as sc
from astropy.io import fits
import astropy.units as u
import astropy.constants as const
from spectral_cube import BooleanArrayMask

c=const.c

def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('km s-1')
    
def ztokms(z):
    return c*z

infile="/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage/spw1minimize.image.pbcor_line.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_large_z0_0002306756533745274_5-6mhzwidth_stdfixes/spectralslabs/km_s/CH3OH~5_1-4_2E1vt0_slab.fits"
infile2="/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/mom1/CH3OH~8_0-7_1E1vt0.fits"#"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/OctReimage_large_z0_0002306756533745274_5-6mhzwidth_stdfixes/mom1/CH3OH~5_1-4_2E1vt0.fits"
print(f'Grabbing mom1 from {infile2}')

reffreq=220027805942.10373*u.Hz#216895488491.4782*u.Hz#redshifted
z0=0.0002306756533745274
v0=c*z0

restfreq=reffreq*(1+z0)

rawcube=sc.read(infile)
datacube=rawcube.spectral_slab(reffreq-15*u.MHz, reffreq+15*u.MHz)
wcs=datacube.wcs
mom1=fits.getdata(infile2)*u.km/u.s

spatdims=datacube.shape[1:]

spectraxis=datacube.with_spectral_unit('Hz').spectral_axis

shell=np.empty([len(spectraxis),spatdims[0],spatdims[1]],dtype=bool)

width=7*u.MHz
print(f'Per-pixel mask width: {width}')
#shellcube=sc(data=shell,wcs=wcs).with_spectral_axis('Hz')

print('Begin masking procedure')
for y in range(spatdims[0]):
    print(f'Start Row {y}')
    for x in range(spatdims[1]):
        dopplershift=mom1[y,x]
        #print(f'dopplershift from reference: {dopplershift}')
        v_shifted=v0+dopplershift
        z_shifted=v_shifted/c
        freq_shifted=restfreq/(1+z_shifted)
        #print(f'shifted frequency: {freq_shifted}')
        tempmask=np.logical_and(spectraxis < (freq_shifted+width), spectraxis > (freq_shifted-width))
        shell[:,y,x]=tempmask
        

shellmaskcube=BooleanArrayMask(mask=shell,wcs=wcs)

maskedcube=datacube.with_mask(shellmaskcube)

print('Saving')
datacube.write('unmaskedspw08-7slab.fits',overwrite=True)
maskedcube.write(f'{width.value}mhzmasked_spw08-7slab.fits',overwrite=True)
shellmaskcube.write(f'{width.value}mhz8-7mask.fits',overwrite=True)
print('Done.')
#maskedcube.to_ds9('7f000001:41898')
#condition=restfreq+(0.77*u.km/u.s)#and minus
