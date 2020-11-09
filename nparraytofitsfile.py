import numpy as np
from astropy.io import fits
import astropy.units as u
from spectral_cube import SpectralCube as sc

swapaxis2toaxis0=np.swapaxes(nugserrormap,0,2)

swapxwithy=np.swapaxes(swapaxis2toaxis0,1,2)

nugscube=fits.PrimaryHDU(swapxwithy)

nugscube.header=transmom0header

nugscube.header['BUNIT']='cm-2'

nugscube.header['BTYPE']='Upper-state column density error'

nugshdul=fits.HDUList([nugscube])

nugshdul.writeto('NupperColDens/field1/testcore1/debug/alltransitions_nupper_error.fits',overwrite=True)
