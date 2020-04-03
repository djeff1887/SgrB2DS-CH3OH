import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
        
home='/home/d.jeff/SgrB2DS_field1/AutoMoments/'
fname='/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'
cube=sc.read(fname)

freqs=cube.spectral_axis#Hz
cube_unmasked=cube.unmasked_data
data=cube_unmasked[:,368,628]#[:,383,649]#Jy/beam
#test=cube_unmasked[:,383:413,649:679]

print(np.shape(data))

#sum=np.sum(test[:,0:,0:])
#print(np.shape(sum))

spec=np.stack((freqs,data),axis=1)
print(np.shape(spec))

plt.plot(spec[:,0],spec[:,1])
#plt.show()

z=0.0002333587
freq_max=freqs[0]*(1+z)#215*u.GHz
print(freq_max)
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
print(freq_min)
linelist='JPL'

table = utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=[linelist],
                                show_upper_degeneracy=True))
                                
def lineloop(line_list,line_width,iterations,quantum_numbers):
    for i in range(iterations):
        line=line_list[i]*u.Hz
        #print('Computing frequency limits')
        nu_upper=line+line_width
        nu_lower=line-line_width
        #print('Done')
        #print('Make spectral slab')
        slab=cube.spectral_slab(nu_upper,nu_lower)
        #print('Done')
        #print('Moment 0')
        slabmom0=slab.moment0()
        #print('Done')
        #print('Saving...')
        transition=qn_replace(quantum_numbers[i])
        #name='test'+str(i)
        slabmom0.write(home+'CH3OH~'+transition+'.fits')
        print('Done')
        
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    return string

#print(np.shape(lines))
lines=table['Freq']*10**9
qns=table['QNs']
#testqn=qns[0:50]
#print(np.size(qns))
#print(testqn)


#print(testqn)
'''
for i in range(len(test)):
    plt.axvline(x=test[i],color='red')
plt.show()
'''
linewidth=0.0097*u.GHz#from small line @ 219.9808GHz
lineloop(lines,linewidth,len(lines),qns)