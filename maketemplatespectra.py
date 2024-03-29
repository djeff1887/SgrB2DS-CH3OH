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
import math

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

source='SgrB2S'
dopplershifts={'SgrB2S':0.000234806,'DSi':0.000186431}#:0.000190713}

otherspecies=[' CH3OCHO ',' HOONO ',' HNCO ',' H2CO ',' C2H5OH ',' CH3CHO ',' CH3NH2 ', ' DCN ',' OCS ',' CH3OCH3 ',' SO ',' NH2CHO ']#' g-CH3CH2OH ' ' C3H6O2 ' skeptical looking at full spectrum, removing DCN b/c it causes issues here
specplotnames=['CH$_3$OCHO','HOONO','HNCO','H$_2$CO','C$_2$H$_5$OH','CH$_3$CHO','CH_3NH$_2$','DCN','OCS','CH$_3$OCH$_3$', 'SO', 'NH$_2$CHO']#' C$_3$H$_6$O$_2$ ',' g-CH$_3$CH$_2$OH '
colors=cm.rainbow(np.linspace(0,1,len(otherspecies)))

z=dopplershifts[source]

inpathdict={'SgrB2S':'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage_K/','DSi':'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSi/field10originals/'}
inpath=inpathdict[source]
path=glob.glob(inpath+'*spw3*')

homedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/"}
home=homedict[source]+'mastereuksqnsfreqsdegens.txt'
stdhome=homedict[source]+'errorimgs/std/spw3fluxstd.fits'

detections=np.loadtxt(home,dtype=str)
detshape=np.shape(detections)

cube=sc.read(path[0])
stdimage=fits.getdata(stdhome)*u.K

cubewcs=cube.wcs
pixdict={'SgrB2S':[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]],'DSi':[[0,0,0],[266.8316149,-28.3972040,0]],'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]]}
#targetworldcrd=[[0,0,0],[266.8316149,-28.3972040,0]]#DSi
targetworldcrd=pixdict[source] #SgrB2S
#[[0,0,0],[266.8332569, -28.3969, 0]] #DSii/iii
targetpixcrd=cubewcs.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    
pixxcrd,pixycrd=int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1]))
print(f'x: {pixxcrd}/y: {pixycrd}')
assert pixxcrd >= 0 and pixycrd >= 0, 'Negative pixel coords'
spw=cube[:,pixycrd,pixxcrd]
stdvalue=stdimage[pixycrd,pixxcrd]

freqs=cube.spectral_axis#Hz
freqflip=False
if freqs[1] < freqs[0]:
    freqs=freqs[::-1]
    freqflip=True
    print('Corrected decreasing frequency axis')
else:
    pass

freq_min=freqs[0]
#print(freq_max)
freq_max=freqs[(len(freqs)-1)]

assert freq_max > freq_min, 'Inverted spectral axis'
print('Passed increasing spectral axis check')

print('Plotting spectra')

xlims=(232707306565.9205,232919343267.9113)
emitxlims=(xlims[0]*(1+z),xlims[1]*(1+z))
plt.plot(freqs,spw.value,drawstyle='steps',color='black')
plt.xlim(xlims[0],xlims[1])
plt.ylim(ymax=100)#100 for SgrB2S

firstline=True
print('Begin detected line vlines')
for row in range(detshape[0]):
    line=float(detections[row,2])
    if line >= xlims[0] and line <= xlims[1]:
        if firstline:
            print('Begin detected line vlines')
            plt.axvline(x=float(detections[row,2]),linestyle='--',color='blue', label='CH$_3$OH')
            firstline=False
        elif not firstline:
            plt.axvline(x=float(detections[row,2]),linestyle='--',color='blue')
    else:
        print(f'Line {detections[row,1]} at {(line*u.Hz).to("GHz")} skipped; out of frequency range')
        
for mol,molname,linecolor in zip(otherspecies,specplotnames,colors):
    firstline=True
    temptable=Splatalogue.query_lines(emitxlims[0]*u.Hz, emitxlims[1]*u.Hz,energy_max=1840, energy_type='eu_k', chemical_name=mol, line_lists=['CDMS','JPL','SLAIM'],show_upper_degeneracy=True,only_NRAO_recommended=True)
    if len(temptable)!=0:
        print(f'{mol} lines available in range.')
        table = utils.minimize_table(temptable)
        otherlines=(table['Freq']*10**9*u.Hz)/(1+z)#Redshifted
        euks=table['EU_K']
        qns=table['QNs']
        for line,euk,qn in zip(otherlines,euks,qns):
            if euk < 1000:
                closest_channel=cube.closest_spectral_channel(line)
                if spw[closest_channel] >= 3*stdvalue:
                    if firstline:
                        print('Begin plotting')
                        plt.axvline(x=line.value,linestyle='--',color=linecolor, label=molname)
                        print(f'{qn}; {euk} K; {line} Hz')
                        firstline=False
                    elif not firstline:
                        plt.axvline(x=line.value,linestyle='--',color=linecolor)
                        print(f'{qn}; {euk} K; {line} Hz')
                else:
                    print(f'Line {qn} at {line} below 3 sigma level ({3*stdvalue})')
            else:
                print(f'Line {qn} with Eupper {euk} at {line} Hz above 1000 K threshold.')
    else:
        print(f'No {mol} lines available in range.')


plt.xlabel(r'$\nu$ (Hz)',fontsize=14)
plt.ylabel('T$_b$ (K)',fontsize=14)
plt.rcParams['figure.dpi'] = 150
plt.legend()
plt.show()
