import numpy as np
from astropy.table import QTable, Table, vstack
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib as mpl
import pdb
import glob
import sys
import scipy.constants as cnst
import math
from utilities import *

def reallum(T,Terr):
    l=(4*np.pi*geomeanbeamarea**2)*sigmasb*T**4
    err=np.sqrt(((16*np.pi*(geomeanbeamarea)**2)*sigmasb*(T**3))**2*(Terr)**2)
    return l.to('solLum'), err.to('solLum')

sigmasb=cnst.sigma*u.W/(u.m**2*u.K**4)

geomeanbeamarea=1908.50019812*u.AU

home='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/'

comptable=QTable.read('/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/sep2023revolution/pacman_densabunpowersummarytable.fits')
lums=np.array(comptable['Luminosity'])#list(comptable[6])[1:]
errorlum=np.array(comptable['Luminosity_error'])#list(comptable[7])[1:])#list(comptable[7])[1:]
temps=np.array(comptable['T_max'])*u.K#list(comptable[0])[1:]
errortemps=np.array(comptable['T_max_error'])*u.K#list(comptable[1])[1:]

methlums=[]
errmethlums=[]

print('Compute upper limit luminosity from methanol')
for trot,errtrot in zip(temps,errortemps):
    newlum,newlumerr=reallum(trot,errtrot)
    upper=newlum+newlumerr#Computes the upper limit based on the 1 sigma error, used to get the error in dex
    dex=np.log10((upper/newlum).value)
    #lower=newlum-(1-(newlumerr
    oom=5#math.floor(np.log10(newlum.value))
    scaled=np.log10(newlum.value)#upper/10**oom
    methlums.append(scaled)
    errmethlums.append(dex)
    print(f'Dex: {round(dex,2)}')
    print(f'log10(solLum):{round(scaled,2)}')
    #print(f'{round(scaled,2)}*10^{oom} solLum')
    
comptable.add_columns([methlums,errmethlums],names=['Lstar_ch3oh','Lstar_ch3oh error'])
comptable.write(datadir+'pacman_lstarch3oh_densabunpowersummarytable.fits')


#print(comptable)
