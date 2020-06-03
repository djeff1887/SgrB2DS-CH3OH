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
chem= input('Molecule?: ')
chem=(' '+chem+' ')
linelist=input('Linelist? (Lovas, SLAIM, JPL, CDMS, ToyoMA, OSU, Recomb, Lisa, RFI): ')
linewidth=0.00485*u.GHz


speciesdata={}
imgnames=['spw2','spw1','spw0']

for i in range(len(files)):
    fname=files[i]#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'#'/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw0_cube.image.fits'
    cube=sc.read(fname)
    header=fits.getheader(fname)

    freqs=cube.spectral_axis
    
    numchans=int(round(np.abs((linewidth.to('Hz')).value/(freqs[1].value-freqs[0].value))))

    freq_max=freqs[np.argmin(freqs)]*(1+z)#215*u.GHz
    freq_min=freqs[np.argmax(freqs)]*(1+z)#235*u.GHz
    
    methanol_table= utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True))

    mlines=(methanol_table['Freq']*10**9)/(1+z)
    mqns=methanol_table['QNs']
    
    temptable=Splatalogue.query_lines(freq_min, freq_max, chemical_name=chem,
                                energy_max=1840, energy_type='eu_k',
                                line_lists=[linelist],show_upper_degeneracy=True)
    if len(temptable)==0:
        print('No '+chem+' lines in frequency range '+str(freq_min)+'-'+str(freq_max)+' GHz.')
        continue
    else:
        print('Lines identified in '+imgnames[i]+'.')
        table = utils.minimize_table(temptable)

        table2=Splatalogue.query_lines(freq_min, freq_max,
                                energy_max=1840, chemical_name=chem, energy_type='eu_k',
                                line_lists=[linelist],
                                show_upper_degeneracy=True, only_NRAO_recommended=True)
    lines=(table['Freq']*10**9)/(1+z)#Redshifted
    qns=table['QNs']
    speciesdata[imgnames[i]]={'freqs':(lines*u.Hz),'qns':qns,'EU_K':(table['EU_K']),'methanoltable':methanol_table}#'lines':lines,'qns':qns}
    
    print('Plotting lines...')
    if i == 2:
        spw=cube[:,762,496]
        fig=plt.figure()
        ax=plt.subplot(111)
        plt.plot(freqs,spw.value,drawstyle='steps')
        for j in range(len(mlines)):
            centroid=mlines[j]*u.Hz
            minfreq=centroid-linewidth
            maxfreq=centroid+linewidth
            centrchan=int(cube.closest_spectral_channel(centroid))
            #interval=np.linspace(cube.closest_spectral_channel(maxfreq.value),cube.closest_spectral_channel(minfreq.value),numchans*2)
            ax.axvline(x=centroid.value,color='green')
            ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
        for k in range(len(lines)):
            ax.axvline(x=lines[k],color='red')
            #plt.annotate((lines[i]),xy=(lines[i],0),xytext=(lines[i],(spw1[i].value+0.01)),rotation=90)
        ax.set_title((imgnames[i]+chem+linelist+' '+'Spectral Sleuthing'))
        ax.set_ylabel('Jy/beam')
        ax.set_xlabel('Frequency (Hz)')
        plt.show()
        continue
    elif i != 2:
        spw=cube[:,649,383]
        '''
        if (freqs[0]-freqs[1])<0:
            freqs=freqs[::-1]
            pass
        else:
            pass
        '''
        fig=plt.figure()
        ax=plt.subplot(111)
        plt.plot(freqs,spw.value,drawstyle='steps')
        for k in range(len(mlines)):
            centroid=mlines[k]*u.Hz
            minfreq=centroid-linewidth
            maxfreq=centroid+linewidth
            centrchan=int(cube.closest_spectral_channel(centroid))
            #interval=np.linspace((maxfreq.value),(minfreq.value),numchans*2)
            ax.axvline(x=centroid.value,color='green')
            if (freqs[0]-freqs[1])<0:
                ax.plot(freqs[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],spw.value[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],drawstyle='steps',color='orange')
            else:
                ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
            #ax.plot(interval,spw.value[(centrchan-numchans):(centrchan+numchans)],drawstyle='steps')
        for l in range(len(lines)):
            ax.axvline(x=lines[l],color='red')
            #plt.annotate((lines[i]),xy=(lines[i],0),xytext=(lines[i],(spw1[i].value+0.01)),rotation=90)
        ax.set_title((imgnames[i]+chem+linelist+' '+'Spectral Sleuthing'))
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Jy/beam')
        plt.show()
        continue

print('Red lines are'+chem+', green lines are CH3OH.')