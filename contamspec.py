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

files=glob.glob('/ufrc/adamginsburg/d.jeff/imaging_results/*.fits')
z=0.0002333587
#chem= input('Molecule?: ')
#chem=(' '+chem+' ')
contaminants=['CH3OCHO','HOONO']

linelist='JPL'#input('Linelist? (Lovas, SLAIM, JPL, CDMS, ToyoMA, OSU, Recomb, Lisa, RFI): ')

contamdata={}
imgnames=['spw2','spw1','spw0']

for i in range(len(files)-2):
    print('Getting ready')
    cube=sc.read(files[i])
    header=fits.getheader(files[i])
    
    freqs=cube.spectral_axis
    
    freq_max=freqs[0]*(1+z)#215*u.GHz
    freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    linewidth=0.00485*u.GHz#Half of original 0.0097GHz
            
    '''Generate methanol table for contaminant search'''    
    methanol_table= utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True))
        
    mlines=(methanol_table['Freq']*10**9)/(1+z)
    mqns=methanol_table['QNs']
    print('Setting figure variable')
    fig=plt.figure()
    if i == 2:
        print('Plotting backdrop spectrum')
        spw=cube[:,762,496]
        plt.plot(freqs,spw.value,drawstyle='steps')
        plt.ylabel('Jy/beam')
        plt.xlabel('Frequency (Hz)')
        plt.title((imgnames[i]+' '+'Contaminant-labeled Spectra'))
        ax=plt.subplot(111)
        print('Plotting mlines')
        for a in range(len(mlines)):
            centroid=mlines[a]*u.Hz
            minfreq=centroid-linewidth
            maxfreq=centroid+linewidth
            contamtable=Splatalogue.query_lines(minfreq, maxfreq,energy_max=1840, energy_type='eu_k', chemical_name=chem, line_lists=[linelist],show_upper_degeneracy=True)
            if len(contamtable)==0:
                print('No '+chem+' lines in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
            else:
                print(chem+' contaminants identified for CH3OH '+mqns[a]+'.')
                table = utils.minimize_table(contamtable)
            if a == 0:
                ax.axvline(x=centroid.value,color='green',label='CH3OH')
            else:
                ax.axvline(x=centroid.value,color='green')
            ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
        plt.show()
        
    elif i != 2:
        print('Plotting backdrop spectrum')
        spw=cube[:,649,383]
        plt.plot(freqs,spw.value,drawstyle='steps')
        plt.ylabel('Jy/beam')
        plt.xlabel('Frequency (Hz)')
        plt.title((imgnames[i]+' '+'Contaminant-labeled Spectra'))
        ax=plt.subplot(111)
        print('Plotting mlines')
        for b in range(len(mlines)):
            centroid=mlines[b]*u.Hz
            minfreq=centroid-linewidth
            maxfreq=centroid+linewidth
            for k in range(len(contaminants)):
                contamtable=Splatalogue.query_lines(minfreq, maxfreq,energy_max=1840, energy_type='eu_k', chemical_name=contaminants[k], line_lists=[linelist],show_upper_degeneracy=True)
                if len(contamtable)==0:
                    print('No '+contaminants[k]+' lines in frequency range '+str(minfreq)+'-'+str(maxfreq)+' GHz.')
                else:
                    print(contaminants[k]+' contaminants identified for CH3OH '+mqns[b]+'.')
                    table = utils.minimize_table(contamtable)
                    line=(table['Freq']*10**9)/(1+z)#Redshifted
                    qns=table['QNs']
                    for f in range(len(table)):
                        ax.axvline(x=line[f],color='red')
            if b == 0:
                ax.axvline(x=centroid.value,color='green',label='CH3OH')
            else:
                ax.axvline(x=centroid.value,color='green')
            if (freqs[0]-freqs[1])<0:
                ax.plot(freqs[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],spw.value[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],drawstyle='steps',color='orange')
            else:
                ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
        plt.legend()
        plt.show()
            #ax.plot(interval,spw.value[(centrchan-numchans):(centrchan+numchans)],drawstyle='steps')
    '''    
    for j in range(len(contaminants)):
        contamtable=Splatalogue.query_lines(m, mfreqmax,energy_max=1840, energy_type='eu_k', chem_name=chem, line_lists=[linelist],show_upper_degeneracy=True)
        if len(temptable)==0:
            print('No '+chem+' lines in frequency range '+str(mfreqmin)+'-'+str(mfreqmax)+' GHz.')
            continue
        else:
            print(chem+' contaminants identified for CH3OH '+mqns[i]+'.')
                table = utils.minimize_table(temptable)
        
                '''
    '''
                table2=Splatalogue.query_lines(freq_min, freq_max,
                                        energy_max=1840, chemical_name=chem, energy_type='eu_k',
                                        line_lists=[linelist],
                                        show_upper_degeneracy=True, only_NRAO_recommended=True)
    '''
    '''
            species=table['Species']
            lines=(table['Freq']*10**9)/(1+z)#Redshifted
            qns=table['QNs']
            print('Plotting spectra')
            spw=cube[:,762,496]#[:,649,383]-spw1&2
            
            ax=plt.subplot(111)
            
            
            for stuff in range(len(species)):
                if stuff > 0:
                    if species[stuff] == species[stuff-1]:#Prevents repeats of same molecule plotting, for efficiency
                        continue
                elif species[stuff]=='CH2Cl2':
                    print('Skipping CH2Cl2...')#Not sure why, but this molecule doesn't have a Splatalogue entry/isn't queried properly in spw2
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
        
                    for k in range(len(contam_lines)):
                        ax.axvline(x=contam_lines[k],color='red')
                    
                    ax.set_title((imgnames[0]+' '+species[stuff]+' '+linelist+' Precision Sleuthing'))
                    ax.set_xlabel('Frequency (Hz)')
                    ax.set_ylabel('Jy/beam')
                    plt.show()
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