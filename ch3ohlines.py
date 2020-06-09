import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
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
contaminants=[' CH3OCHO ',' HOONO ',' C3H6O2 ',' g-CH3CH2OH ']
colors=cm.rainbow(np.linspace(0,1,len(contaminants)))

linelist='JPL'#input('Linelist? (Lovas, SLAIM, JPL, CDMS, ToyoMA, OSU, Recomb, Lisa, RFI): ')

mdict={}
contamdata={}
imgnames=['spw2','spw1','spw0']

def specmaker(plot,x,y,xmin,xmax,center,trans):
    plot.set_xlim(xmin.value,xmax.value)
    plot.axvline(x=center,color='green',linestyle='--')
    print(len(y.value[cube.closest_spectral_channel(xmax):cube.closest_spectral_channel(xmin)]))
    left=cube.closest_spectral_channel(xmin)
    right=cube.closest_spectral_channel(xmax)
    tempx=x[0:90]
    tempy=y.value[0:90]
    if left > right:
        plot.set_ylim(-100,100)
        print(tempx[0:5])
        print(tempy[0])
        print(type(right))
        print('Inverted cube')
        plot.plot(tempx,tempy,drawstyle='steps')
        #plot.axhline(y=0,color='red',linestyle='--')
    else:
        plot.set_ylim(-1,10)
        plot.plot(x[left:right],y.value[left:right],drawstyle='steps')
    

for i in range(len(files)-1):
    print('Getting ready - '+imgnames[i])
    cube=sc.read(files[i])
    header=fits.getheader(files[i])
    
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
    linewidth=0.00485*u.GHz#Half of original 0.0097GHz
            
    '''Generate methanol table for contaminant search'''    
    methanol_table= utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True))
        
    mlines=(methanol_table['Freq']*10**9)/(1+z)
    mqns=methanol_table['QNs']
    
    mins=[]
    maxs=[]
    if i == 2:
        print('Plotting backdrop spectrum')
        spw=cube[:,762,496]
        plt.plot(freqs,spw.value,drawstyle='steps')
        plt.ylabel('Jy/beam')
        plt.xlabel('Frequency (Hz)')
        plt.title((imgnames[i]+' '+'Contaminant-labeled Spectra'))
        ax=plt.subplot(111)
        print('Plotting mlines and ROIs')
        for a in range(len(mlines)):
            centroid=mlines[a]*u.Hz
            minfreq=centroid-linewidth
            maxfreq=centroid+linewidth
            mins.append(minfreq)
            maxs.append(maxfreq)
            if a == 0:
                ax.axvline(x=centroid.value,color='green',label='CH3OH')
            else:
                ax.axvline(x=centroid.value,color='green')
            ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
        print('Begin plotting contaminant lines')
        for j in range(len(contaminants)):
            print('Checking'+contaminants[j])
            dum=0
            for d in range(len(mins)):
                contamtable=Splatalogue.query_lines((mins[d]*(1+z)), (maxs[d]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[j], line_lists=[linelist],show_upper_degeneracy=True)
                if len(contamtable)==0:
                    print('No '+contaminants[j]+' lines in frequency range '+str(mins[d])+'-'+str(maxs[d])+'.')
                else:
                    print(contaminants[j]+' contaminants identified for CH3OH '+mqns[d]+' at '+str(mins[d]+linewidth)+' GHz.')
                    table = utils.minimize_table(contamtable)
                    line=(table['Freq']*10**9)/(1+z)#Redshifted
                    qns=table['QNs']
                    for g in range(len(table)):
                        if g==0 and dum==0:
                            ax.axvline(x=line[g],color=colors[j],label=contaminants[j])
                            print('hiii')
                            dum+=1
                        else:
                            ax.axvline(x=line[g],color=colors[j])
        plt.legend()
        plt.show()
        
    elif i != 2:
        print('Setting figure and ax variables')
        numrows=3
        numcols=5
        fig,ax=plt.subplots(numrows,numcols,sharey=True)
        print('Gathering intensity values')
        spw=cube[:,649,383]
        if freqflip==True:
            spw=spw[::-1]
        else:
            pass

        
        '''        
        plt.plot(freqs,spw.value,drawstyle='steps')
        plt.ylabel('Jy/beam')
        plt.xlabel('Frequency (Hz)')
        plt.title((imgnames[i]+' '+'Contaminant-labeled Spectra'))
        ax=plt.subplot(111)
        '''
        print('Generating mline data')
        for b in range(len(mlines)):
            centroid=mlines[b]*u.Hz
            minfreq=centroid-(linewidth*1.5)
            maxfreq=centroid+(linewidth*1.5)
            mins.append(minfreq)
            maxs.append(maxfreq)
            
        print('Begin figure plot loops')
        rowoffset=0
        for row in range(numrows):
            print('Start Row '+str(row)+'.')
            for col in range(numcols):
                if row == (numrows-1):
                    if col >= int(len(mlines)/numrows):
                        break
                    else:
                        specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+rowoffset],maxs[col+rowoffset], mlines[col+rowoffset], mqns[col+rowoffset])
                else:
                    specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+rowoffset],maxs[col+rowoffset], mlines[col+rowoffset], mqns[col+rowoffset])
                '''
                if row == 0:
                    specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col],maxs[col],mlines[col],mqns[col])
                    continue
                if row == 1:
                    specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+5],maxs[col+5], mlines[col+5],mqns[col+5])
                    continue
                if row == 2:
                    if col >= int(len(mlines)/numrows):
                        break
                    else:
                        specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+10],maxs[col+10], mlines[col+10], mqns[col+10])
                        continue
                '''
            rowoffset+=5
        plt.title(imgnames[i]+' CH3OH Lines')
        print('Plotting complete. plt.show()')
        plt.show()
        '''
            if b == 0:
                ax.axvline(x=centroid.value,color='green',label='CH3OH')
            else:
                ax.axvline(x=centroid.value,color='green')
            if (freqs[0]-freqs[1])<0:
                ax.plot(freqs[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],spw.value[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],drawstyle='steps',color='orange')
            else:
                ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
        '''
        '''
        print('Begin plotting contaminant lines')
        for k in range(len(contaminants)):
            print('Checking'+contaminants[k]+'...')
            dummy=0
            for c in range(len(mins)):
                contamtable=Splatalogue.query_lines((mins[c]*(1+z)), (maxs[c]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[k], line_lists=[linelist],show_upper_degeneracy=True)
                if len(contamtable)==0:
                    print('No '+contaminants[k]+' lines in frequency range '+str(mins[c])+'-'+str(maxs[c])+'.')
                    continue
                else:
                    print(contaminants[k]+' contaminants identified for CH3OH '+mqns[c]+' in frequency range '+str(mins[c])+'-'+str(maxs[c])+'.')
                    table = utils.minimize_table(contamtable)
                    line=(table['Freq']*10**9)/(1+z)#Redshifted
                    qns=table['QNs']
                    dummy+=1
                    for f in range(len(table)):
                        if f == 0 and dummy == 1:
                            ax.axvline(x=line[f],color=colors[k],label=contaminants[k])
                        else:
                            ax.axvline(x=line[f],color=colors[k])
        plt.legend()
        plt.show()
        '''
            #ax.plot(interval,spw.value[(centrchan-numchans):(centrchan+numchans)],drawstyle='steps')