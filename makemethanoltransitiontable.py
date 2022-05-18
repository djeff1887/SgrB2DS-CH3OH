from utilities import *
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import astropy.units as u
import astropy.constants as cnst

c=cnst.c

source = 'SgrB2S'
sourcepath=sourcedict[source]

sgrb2sz=0.000234806

mom0home=sourcepath+'mom0/*masked.fits'

txtfile=sourcepath+'mastereuksqnsfreqsdegens.txt'
array=np.genfromtxt(txtfile,dtype=str)
transitions=array[:,1]
freqs=(array[:,2].astype(np.float)*u.Hz).to('GHz')
freqs=[str(round((x.value*(1+sgrb2sz)),5)) for x in freqs]

stack=hstack([transitions,freqs])
table=Table(stack,names=['Transition','Frequency'])
table.write('methanoltransitiontable.fits',overwrite=True)
table.write('methanoltransitiontable.csv',overwrite=True)
#mom0s=glob.glob(mom0home)

print(freqs)

