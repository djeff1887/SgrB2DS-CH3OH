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

source='DSVI'
sourceisnew=True

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K

trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K}
testT=trotdict[source]
qrot_partfunc=Q_rot_asym(testT).to('')
testntot=1e16*u.cm**-2

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

dopplershifts={'SgrB2S':0.000234806,'DSi':0.000186431,'DSii':0.00015954965399894244,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407}#:0.000190713}/old doppler S: 0.0002306756533745274/0.00015954965399894244/0.00016236367659115043

z=dopplershifts[source]
z_vel=z*c

sourcelocs={'SgrB2S':'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage_K/*.fits','DSi':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSi/field10originals_K/*.fits",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiii/field10originals/*.fits",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiv/field10originals_K/*.fits",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSv/field10originals_K/*.fits",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSVI/field2originals_K/*.fits"}
dopplershiftimg={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/4_2-3_1vt=0repline_mom1.fits",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial5carryover_field10errors/8_1-7_0vt=0repline_mom1.fits",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/mom1/CH3OH~8_0-7_1E1vt0.fits",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/field10originals_noexclusions/mom1/CH3OH~10_2--9_3-vt0.fits"}

sourcelocs={'SgrB2S': r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/SgrB2S/OctReimage_K','DSi':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSi/field10originals_K','DSii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSii/field10originals_K','DSiii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiii/field10originals_K','DSiv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiv/field10originals_K','DSv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSv/field10originals_K','DSVI':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVI/field2originals_K'}

sourcepath=sourcelocs[source]

print(f'Collecting spectra from {sourcepath}')
inspecs=glob.glob(sourcepath+'/*.txt')

images=['spw0','spw1','spw2','spw3']

spectra=[]

for spew in images:
    for f1 in inspecs:
        if spew in f1:
            spectra.append(f1)
            continue
    
assert 'spw0' in spectra[0], 'Cube list out of order'

print('Spectra in sequential order')

#imgnum=0
testline=0

linewidth=3*u.km/u.s#2.5 km/s is ideal for DSVI
print(f'Absolute model line width: {linewidth}\n')

#cube=sc.read(spectrums[imgnum])

#targetworldcrd=[[0,0,0],[266.8321311,-28.3976633,0]]#DSv
#targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]#SgrB2S

for spectrum, img in zip(spectra,images):
    print('Getting ready - '+img)
    spec=np.genfromtxt(spectrum)
    '''
    cube.allow_huge_operations=True
    cube_w=cube.wcs#WCS(files[imgnum])
    targetworldcrd=[[0,0,0],[266.8380037, -28.4050741,0]]#DSVI cont peak
    #targetworldcrd=[[0,0,0],[266.8321331, -28.3976585, 0]]#DSv cont peak
    #targetworldcrd=[[0,0,0],[266.8323834, -28.3954424,0]]#DSiv cont peak
    #targetworldcrd=[[0,0,0],[266.8331150,-28.3968609,0]]#DSiii lowest S/N
    #targetworldcrd=[[0,0,0],[266.8332758,-28.3969269,0]]#DSiii cont peak
    #targetworldcrd=[[0,0,0],[266.83183633,-28.3971697,0]]#Near center of new DSi offset hotspot
    #targetworldcrd=[[0,0,0],[266.8347475,-28.3967675,0]]#Continuum bump south from main SgrB2S body (in "extended" region)
    #targetworldcrd=[[0,0,0],[266.8352867,-28.3958092,0]]#Lowest absorption in 6-7 torsional line (very line poor)
    #targetworldcrd=[[0,0,0],[266.8316149,-28.3972040,0]]#DSi sample pixel
    #targetworldcrd=[[0,0,0],[266.8352814,-28.3959488,0]]#SgrB2S HII region continuum peak
    #targetworldcrd=[[0,0,0],[266.8353290,-28.3962279,0]]#SgrB2S southernmost molecular ridge continuum peak
    #targetworldcrd=[[0,0,0],[266.8354136,-28.3960233,0]]#SgrB2S CH3OH hotspot (pre STATCONT catastrophe)
    #targetworldcrd=[[0,0,0],[266.8353844,-28.3960720,0]]#SgrB2S hotspot-adjacent
    #targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]#SgrB2S sample pixel
    #targetworldcrd=[[0,0,0],[266.8335363,-28.3963158,0]]#DSii sample pixel
    refpix={'SgrB2S':[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]],'DSi':[[0,0,0],[266.8316149,-28.3972040,0]],'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]]}
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    targetypix=int(round(targetpixcrd[1][1]))
    targetxpix=int(round(targetpixcrd[1][0]))
    
    if targetworldcrd != refpix[source]:
        print('Non-reference pixel detected\nAdjusting dopplershift\n')
        vel=z*c
        dif=zconvimg[targetypix,targetxpix]
        velcorr=vel+dif
        z=velcorr/c
    '''
    freqs=spec[:,0]#cube.spectral_axis
    data=spec[:,1]
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        data=data[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    #linewidth=0.00485*u.GHz#Half of original 0.0097GHz
    #lw2=linewidth/8
    originallinewidth=(11231152.36688232*u.Hz/2)
    
            
    '''Generate methanol table for contaminant search'''    
    methanol_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=['JPL','SLAIM','CDMS'], show_upper_degeneracy=True)
    
    if len(methanol_table['Chemical Name']) == 0:
        print(f'No transitions in {img}. Continue')
        continue
    else:
        pass
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
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],np.shape(spec)[0]).to('GHz')
    modelgaus=models.Gaussian1D(mean=freqs[0], stddev=11 * u.MHz, amplitude=0*u.K)
    
    for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
        print(f'Transition: {qn} @ {line.to("GHz")}')
        restline=line*(1+z)
        
        est_nupper=nupper_estimated(testntot,deg,qrot_partfunc,euj,testT).to('cm-2')
        est_tau=opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
        print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        print(f'Estimated brightness: {"{:.3f}".format(trad)}')
        
        modlinewidth=velocitytofreq(linewidth,line)
        print(f'Model linewidth (Hz): {modlinewidth}')
        modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
        #modelgaus+=modelline
        modelspec+=modelline

  
    print('Plotting model spectra')
    #plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],cube.shape[0])
    #cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].quicklook()
    plt.plot(freqs,data,drawstyle='steps-mid')
    plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color='black')
    #plt.plot(plot,cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].value,color='black',drawstyle='steps')
    plt.xlabel(r'$\nu$ (Hz)')
    plt.ylabel('T$_b$ (K)')
    plt.title(f'modelTex: {testT}, modelntot: {testntot}')
    print('Overplotting axvlines and transition annotations')
    for line,qn in zip(mlines,mqns):
        plt.axvline(x=line.value,linestyle='--',color='yellow',ymin=0.25)
        plt.annotate(qn, (line.value, 0), (line.value-0.002,40),rotation=90)
    plt.show()
