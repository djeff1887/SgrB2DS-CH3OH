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

fname='/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw2_cube_medsub.image.fits'#glob.glob('/ufrc/adamginsburg/d.jeff/imaging_results/*.fits')
z=0.0002333587
#chem= input('Molecule?: ')
#chem=(' '+chem+' ')
linelist=input('Linelist? (Lovas, SLAIM, JPL, CDMS, ToyoMA, OSU, Recomb, Lisa, RFI): ')

speciesdata={}
imgnames=['spw2']

#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw0_cube.image.fits'
cube=sc.read(fname)
header=fits.getheader(fname)

freqs=cube.spectral_axis

freq_max=freqs[0]*(1+z)#215*u.GHz
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
linewidth=0.5*0.0097*u.GHz
    
'''Generate methanol table for contaminant search'''    
methanol_table= utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True))

mlines=(methanol_table['Freq']*10**9)/(1+z)
mqns=methanol_table['QNs']

skips=['CH2Cl2','(CH3)2CO']
knowns=['CH3OH','CH3OCHO','HOONO']

for i in range(len(methanol_table['Freq'])):
    '''Check for contaminants within linewidth around methanol lines'''
    mfreqmin=methanol_table['Freq'][i]*u.GHz-(linewidth/2)
    mfreqmax=methanol_table['Freq'][i]*u.GHz+(linewidth/2)
    
    temptable=Splatalogue.query_lines(mfreqmin, mfreqmax,energy_max=1840, energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True)#chemical_name=chem,
    if len(temptable)==0:
        print('No  lines in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
        continue
    else:
        print('Possible contaminants identified for CH3OH '+mqns[i]+'.')
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
    print('Plotting background spectra')
    if imgnames[0]=='spw0':
        spw=cube[:,762,496]
    else:
        spw=cube[:,649,383]
        
    fig=plt.figure()
    ax=plt.subplot(111)
    plt.plot(freqs,spw.value,drawstyle='steps')
    print('Begin contaminant plotting')
    
    for stuff in range(len(species)):
        for much in range(len(skips)):
            if skips[much] in species[stuff]:
                print(skips[much])
                print(species[stuff])
                print('Skipping problem '+skips[much]+'...')
                continue

        '''
        elif species[stuff]=='CH2Cl2':
            print('Skipping CH2Cl2...')#Not sure why, but this molecule doesn't have a Splatalogue entry/isn't queried properly in spw2
            continue
            '''
        
        for more in range(len(knowns)):#Skips known contaminants
            if knowns[more] in species[stuff]:
                print('Skipping known '+knowns[more]+'...')
                continue   

        if stuff > 0:
            if species[stuff] in species[stuff-1]:#Prevents repeats of same molecule plotting, for efficiency
                print('Skipping duplicate '+species[stuff]+'...')
                continue
          
        else:
            print('Querying lines for '+species[stuff]+' in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
            if 'v' in species[stuff]:#Removes v=___ from chem name for queries
                vequals=species[stuff]
                species[stuff]=vequals.replace(vequals[vequals.index('v'):],'')
                contam_table=utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max,energy_max=1840, chemical_name=species[stuff], energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True))
                contam_lines=(contam_table['Freq']*10**9)/(1+z)
                contam_qns=contam_table['QNs']
            else:
                contam_table=utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max,energy_max=1840, chemical_name=species[stuff], energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True))
                contam_lines=(contam_table['Freq']*10**9)/(1+z)
                contam_qns=contam_table['QNs']
        
            speciesdata[species[stuff]]={'freqs':(contam_lines*u.Hz),'qns':contam_qns,'E_u(K)':(contam_table['EU_K'])}#,'methanoltable':methanol_table}#'lines':lines,'qns':qns}
        
            print('Plotting lines...')
            for j in range(len(mlines)):
                ax.axvline(x=mlines[j],color='green')
            for k in range(len(contam_lines)):
                ax.axvline(x=contam_lines[k],color='red')
            
            ax.set_title((imgnames[0]+' '+species[stuff]+' '+linelist+' Precision Sleuthing'))
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Jy/beam')
            plt.show()



print('Red lines are chems, green lines are CH3OH.')