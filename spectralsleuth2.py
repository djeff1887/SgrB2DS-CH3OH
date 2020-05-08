import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import regions

fname='/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'#glob.glob('/ufrc/adamginsburg/d.jeff/imaging_results/*.fits')
z=0.0002333587
#chem= input('Molecule?: ')
#chem=(' '+chem+' ')
linelist=input('Linelist? (Lovas, SLAIM, JPL, CDMS, ToyoMA, OSU, Recomb, Lisa, RFI): ')

speciesdata={}
imgnames=['spw1']

#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw0_cube.image.fits'
cube=sc.read(fname)
header=fits.getheader(fname)

freqs=cube.spectral_axis

freq_max=freqs[0]*(1+z)#215*u.GHz
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
linewidth=0.5*0.0097*u.GHz
    
methanol_table= utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True))

mlines=(methanol_table['Freq']*10**9)/(1+z)
mqns=methanol_table['QNs']

for i in range(len(methanol_table['Freq'])):

    mfreqmin=methanol_table['Freq'][i]*u.GHz-(linewidth/2)
    mfreqmax=methanol_table['Freq'][i]*u.GHz+(linewidth/2)
    
    temptable=Splatalogue.query_lines(mfreqmin, mfreqmax,energy_max=1840, energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True)#chemical_name=chem,
    if len(temptable)==0:
        print('No  lines in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
        continue
    else:
        print('Lines identified for CH3OH '+mqns[i]+'.')
        table = utils.minimize_table(temptable)

        '''
        table2=Splatalogue.query_lines(freq_min, freq_max,
                                energy_max=1840, chemical_name=chem, energy_type='eu_k',
                                line_lists=[linelist],
                                show_upper_degeneracy=True, only_NRAO_recommended=True)
        '''
    species=table['Species']
    lines=(table['Freq']*10**9)/(1+z)#Redshifted
    qns=table['QNs']
    print('Plotting spectra')
    spw=cube[:,649,383]
    fig=plt.figure()
    ax=plt.subplot(111)
    plt.plot(freqs,spw.value,drawstyle='steps')
    
    for stuff in range(len(species)):
        print('Querying lines for '+species[stuff]+' in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
        contam_table=utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max,energy_max=1840, chemical_name=species[stuff], energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True))
        contam_lines=(contam_table['Freq']*10**9)/(1+z)
        contam_qns=contam_table['QNs']
        
        speciesdata[stuff]={'name':species[stuff],'freqs':(contam_lines*u.Hz),'qns':contam_qns,'E_u(K)':(contam_table['EU_K'])}#,'methanoltable':methanol_table}#'lines':lines,'qns':qns}
        
        print('Plotting lines...')
        for j in range(len(mlines)):
            ax.axvline(x=mlines[j],color='green')
        for k in range(len(contam_lines)):
            ax.axvline(x=contam_lines[k],color='red')
            
        ax.set_title((imgnames[0]+' '+species[stuff]+' '+linelist+' Spectral Sleuthing (2)'))
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Jy/beam')
        plt.show()
        continue
    '''
    
    '''
    '''
    if i == 2:
        
        
        
            #plt.annotate((lines[i]),xy=(lines[i],0),xytext=(lines[i],(spw1[i].value+0.01)),rotation=90)
        ax.set_title((imgnames[0]+' '+linelist+' Spectral Sleuthing (2)'))
        ax.set_ylabel('Jy/beam')
        ax.set_xlabel('Frequency (Hz)')
        plt.show()
        continue
     '''  
''' 
     spw=cube[:,649,383]
     fig=plt.figure()
     ax=plt.subplot(111)
     plt.plot(freqs,spw.value,drawstyle='steps')
     for k in range(len(mlines)):
       ax.axvline(x=mlines[k],color='green')
     for l in range(len(lines)):
       ax.axvline(x=lines[l],color='red')
            #plt.annotate((lines[i]),xy=(lines[i],0),xytext=(lines[i],(spw1[i].value+0.01)),rotation=90)
        continue
'''




print('Red lines are chems, green lines are CH3OH.')