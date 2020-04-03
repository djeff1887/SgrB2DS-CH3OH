import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob

table = utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=[linelist],
                                show_upper_degeneracy=True))
lines=table['Freq']*10**9
qns=table['QNs']
print(qns)
files=glob.glob('/home/d.jeff/SgrB2DS_field1/AutoMoments/*')
#print(files)

def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix-1),(xpix-1)])
    return vals
'''Rough estimate: 0.75 arcsec/15pixels >>> 0.05 arsec/pixel >>> 2.424e-7arcsec/pixel
Taken from DS9 tradiation region'''    
def solid_angle(rad):
    return np.pi*np.sin(rad.to(u.deg))**2*u.sr
    
    '''
temp=fits.getdata(files[0])
print(temp[648,382])
plt.imshow(temp)
plt.scatter(383,649,color='purple')  
plt.show() ''' 
fluxes=fluxvalues(383,649,files)*u.Jy*u.km/u.s
Omega=solid_angle(2.424e-7*u.rad)
vint_intensities=fluxes/Omega
#vint_trads=vint_intensities*((cnst.c)**2
#print(vint_intensities)


