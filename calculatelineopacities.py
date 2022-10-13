from utilities import *
import astropy.units as u
import astropy.constants as cnst
from pyspeckit.spectrum.models import lte_molecule
import numpy as np
from utilities import velocitytofreq, Q_rot_asym
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from astropy.table import QTable

mpl.interactive(True)
plt.close('all')

def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def N_u(nu,Aij,velocityintegrated_intensity_K):#,velint_intK_err):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    #nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err
    return nuppercalc#,nuppererr

k=cnst.k_B
h=cnst.h
c=cnst.c
sources=sourcedict.keys()
alltablepaths=glob.glob('OpticalDepthTables/*.fits')
testT=[150,300,500]*u.K

for tblpath in alltablepaths:
    s=tblpath.replace('OpticalDepthTables/','')
    s=s.replace('.fits','')
    qtable=QTable.read(tblpath)
    opacities=[]
    masternuppers=[]
    masterngs=[]
    for line in qtable:
        peakflux=line['Peak Velocity Integrated Flux']
        errpeakflux=line['Peak Flux Error']
        restfreq=line['Freq']
        euj=line['EU(J)']
        euk=line['EU(K)']
        aij=line['Aij']
        degen=line['g']
        fwhm=line['Line Width']
        peaktb=(peakflux/fwhm).to('K')
        fwhm_Hz=velocitytofreq(fwhm,restfreq)

        qrot=Q_rot_asym(testT[1])
        nupper=N_u(restfreq,aij,peakflux)
        masternuppers.append(nupper.to('cm-2').value)
        masterngs.append(nupper.to('cm-2').value/degen)
        ntot=lte_molecule.ntot_of_nupper(nupper,euj,peaktb,qrot,degen)
        phi_nu=lineprofile(fwhm_Hz,restfreq,restfreq)

        intertau=lte_molecule.line_tau(peaktb,ntot,qrot,degen,restfreq,euj,aij)
        tau=(intertau*phi_nu).to('')
        opacities.append(tau)

    plt.figure()
    plt.scatter(masterngs,opacities,c=qtable['EU(K)'].value)
    plt.xlabel(r'$N_{u}$ cm$^{-2}$')
    plt.ylabel(r'$\tau$')
    plt.title(f'Qrot({testT[1]})')
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.tight_layout()
    plt.show()
'''    
sgrb2scontsource_43tb=100*u.K#61.9*u.K#peak taken from CARTA
qrot=Q_rot_asym(testT)#3203.211430184744
fourthreeaij=4.686514876992861e-05*u.Hz
restfreq=218.44005*u.GHz
sigma_vel=5*u.km/u.s
sigma_freq=velocitytofreq(sigma_vel,restfreq)
freqrange=np.linspace(218.4*u.GHz,218.5*u.GHz,100)

nupper_pix=N_u(restfreq,fourthreeaij,713*u.K*u.km*u.s**-1)#9.44e13*u.cm**-2
eupper43=45.45988*u.K
eupperJ=(3/2)*k*eupper43
degen=9.0
ntot=lte_molecule.ntot_of_nupper(nupper_pix,eupperJ,sgrb2scontsource_43tb,qrot,degen)


phi_nu=[]
for f in freqrange:
    phi_nu.append(lineprofile(sigma_freq,restfreq,f).value)

#pdb.set_trace()
phi_nu=np.array(phi_nu)*u.MHz**-1

intertau=lte_molecule.line_tau(sgrb2scontsource_43tb,ntot,qrot,degen,restfreq,eupperJ,fourthreeaij)
tau=(intertau*phi_nu).to('')

plt.plot(freqrange,tau)
plt.xlabel('Freq (GHz)')
plt.ylabel(r'{\tau}')
plt.title(f'{testT}, {sigma_vel} line width')
plt.show()

#print(f'Tau = {tau}')
'''
