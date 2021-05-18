import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from scipy.optimize import curve_fit as cf
from astropy.io import fits
import glob
import radio_beam
from astropy.modeling import models, fitting#Fittable1DModel, Parameter, fitting
from utilities import *#Q_rot_asym,mulu,vradio,t_rad,nupper_estimated,opticaldepth,qngrabber
import matplotlib as mpl
#from cubecoretexmap import t_rad,nupper_estimated,opticaldepth
Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

mpl.interactive(True)

plt.close('all')
linelist='JPL'

'''Collect constants for N_tot and N_upper calculations'''

source='SgrB2S'

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K
testT=300*u.K
qrot_partfunc=Q_rot_asym(testT).to('')
testntot=1e17*u.cm**-2

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

dopplershifts={'SgrB2S':0.000234806,'DSi':0.000186431,'DSv':0.000186431}#:0.000190713}/old doppler S: 0.0002306756533745274

z=dopplershifts[source]

files=glob.glob('/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage_K/*.fits')

imgnames=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in imgnames:
    for f1 in files:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

#imgnum=0
testline=0


#cube=sc.read(datacubes[imgnum])



#targetworldcrd=[[0,0,0],[266.8316149,-28.3972040,0]]#DSi
#targetworldcrd=[[0,0,0],[266.8321311,-28.3976633,0]]#DSv
#targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]#SgrB2S




for datacube, img in zip(datacubes,imgnames):
    print('Getting ready - '+img)
    cube=sc.read(datacube)
    cube.allow_huge_operations=True
    cube_w=cube.wcs#WCS(files[imgnum])
    targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]#SgrB2S
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    
    freqs=cube.spectral_axis
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    linewidth=0.00485*u.GHz#Half of original 0.0097GHz
    lw2=linewidth/8
    originallinewidth=(11231152.36688232*u.Hz/2)
    
            
    '''Generate methanol table for contaminant search'''    
    methanol_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=['JPL'], show_upper_degeneracy=True)
    
    minmethtable=utils.minimize_table(methanol_table)
    
    mlines=(minmethtable['Freq']*10**9)/(1+z)*u.Hz
    mqns=minmethtable['QNs']
    meuks=minmethtable['EU_K']*u.K
    meujs=[]
    for euk in meuks:
        meujs.append(KtoJ(euk))
    mdegs=methanol_table['Upper State Degeneracy']
    mlog10aijs=minmethtable['log10_Aij']
    maijs=10**mlog10aijs*u.s**-1
    
    #zeros=np.zeros(cube.shape[0])
    #baseline=np.array(list(zip(plot.value,zeros)))#np.vstack((plot,zeros))
    
    baseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    baseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    modelspec=baseline
    print('Begin model loops')
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],cube.shape[0]).to('GHz')
    modelgaus=models.Gaussian1D(mean=freqs[0], stddev=11 * u.MHz, amplitude=0*u.K)
    for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
        print(f'Transition: {qn} @ {line.to("GHz")}')
        restline=line*(1+z)
        
        est_nupper=nupper_estimated(testntot,deg,qrot_partfunc,euj,testT).to('cm-2')
        est_tau=opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        
        modelline=models.Gaussian1D(mean=line, stddev=1 * u.MHz, amplitude=trad)
        #modelgaus+=modelline
        modelspec+=modelline

  
    print('Plotting model spectra')
    #plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],cube.shape[0])
    cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].quicklook()
    plt.plot(freqs,modelspec(freqs),drawstyle='steps',color='red')
    #plt.plot(plot,cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].value,color='black',drawstyle='steps')
    plt.xlabel(r'$\nu$ (Hz)')
    plt.ylabel('T$_b$ (K)')
    print('Overplotting axvlines and transition annotations')
    for line,qn in zip(mlines,mqns):
        plt.axvline(x=line.value,linestyle='--',color='yellow',ymin=0.25)
        plt.annotate(qn, (line.value, 0), (line.value-0.002,40),rotation=90)
    plt.show()
