#from spectral_cube import SpectralCube as sc
import numpy as np
from astropy.io import fits
import astropy.units as u

NumIterations=25

datatxt='/blue/adamginsburg/d.jeff/imaging_results/XCLASSstuff/MHzKfiles/649_383_spw2mhzkspec.txt'

cube='/blue/adamginsburg/d.jeff/imaging_results/XCLASSstuff/sgrb2ds_649_383_spw2subcubemhz.fits'

header=fits.getheader(cube)

#RestFreq=(header['RESTFRQ']*u.Hz).to('MHz').value

vLSR=69.959#255

expdata=LoadASCIIFile(datatxt)

#ASCIIDataFileName=FileName

TelescopeSize=0.285

Inter_Flag=True

t_back_flag=False

nH_flag=False

MolfitsFileName='/blue/adamginsburg/d.jeff/imaging_results/XCLASSstuff/mols.molfit'

iso_flag=False

FreqStep=(expdata[1][0]-expdata[0][0])

if FreqStep < 0:
    expdata=expdata[::-1]
    FrequencyWidth=(expdata[1][0]-expdata[0][0])
    print('Flipping spectral axis')
else:
    pass

FreqMax=expdata[len(expdata)-1][0]

FreqMin=expdata[0][0]

FreqStep=(expdata[1][0]-expdata[0][0])

assert FreqMax > FreqMin, 'Decreasing frequency/velocity axis'

RestFreq=0

modeldata, log, TransEnergies, IntOpt, jobDir=myXCLASS()

PlotTitle="field1 spw2 CH3OH Synthetic Spectra"

#TransEnergies=""

myXCLASSPlot()


