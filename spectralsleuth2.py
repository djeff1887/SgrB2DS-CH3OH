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

skips=['CH2Cl2','(CH3)2CO', 'c-C6H5CCH','c-C4H4O2','HC7O','H13C(O)NH2','a-(CH3)2CHCH2CN','H2NCH2COOH-Iv=0'
, 'g-(CH3)2CHCH2CN','t-HC(O)SH','CH3CHNH2COOH-I','C2H5CHCNCH3','NCC(O)NH2','c-C6H5OH','(Z)-HC2CHCHCN', 'c-C6H5CN','g\'Ga-(CH2OH)2','H2C(CN)2','AA-n-C4H9CN','GA-n-C4H9CN','OC(CN)2','CH2(OD)CHO','g\'Gg-(CH2OH)2','aa-(C2H5)2O','C2H3NH20&#150;&larr;0+','CH3CHNH2COOH-II','c-HC(O)SH']
knowns=['CH3OH','CH3OCHO','HOONO']

for i in range(len(methanol_table['Freq'])):
    '''Check for contaminants within linewidth around methanol lines'''
    mfreqmin=methanol_table['Freq'][i]*u.GHz-(linewidth/2)
    mfreqmax=methanol_table['Freq'][i]*u.GHz+(linewidth/2)
    
    temptable=Splatalogue.query_lines(mfreqmin, mfreqmax,energy_max=1840, energy_type='eu_k', line_lists=[linelist],show_upper_degeneracy=True,only_NRAO_recommended=True)#chemical_name=chem,
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
        '''
        if any(species[stuff] in knowns):
            print('Skipping known '+species[stuff]+'...')
            continue
            
        '''
        
        if stuff > 0 and species[stuff] == species[stuff-1]:#Prevents repeats of same molecule plotting, for efficiency
            print('Skipping duplicate '+species[stuff]+'...')
            continue
            
        if species[stuff] in skips:
            print('Skipping problem '+species[stuff]+'...')
            continue
        
        if species[stuff] in knowns:
            print('Skipping known '+knowns[more]+'...')
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