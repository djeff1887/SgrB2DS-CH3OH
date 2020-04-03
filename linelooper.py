import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
        
home='/home/d.jeff/SgrB2DS_field1/AutoMoments/'
fname='/ufrc/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw1_cube_medsub.image.fits'
cube=sc.read(fname)

freqs=cube.spectral_axis#Hz
cube_unmasked=cube.unmasked_data
data=cube_unmasked[:,368,628]#[:,383,649]#Jy/beam
#test=cube_unmasked[:,383:413,649:679]

#print(np.shape(data))

#sum=np.sum(test[:,0:,0:])
#print(np.shape(sum))

spec=np.stack((freqs,data),axis=1)
#print(np.shape(spec))

plt.plot(spec[:,0],spec[:,1])
#plt.show()

z=0.0002333587
freq_max=freqs[0]*(1+z)#215*u.GHz
#print(freq_max)
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
#print(freq_min)
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
lines=table['Freq']*10**9*u.Hz
qns=table['QNs']
tex=table['EU_K']
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
#lineloop(lines,linewidth,len(lines),qns)
######

#print(qns)
files=glob.glob('/home/d.jeff/SgrB2DS_field1/AutoMoments/*')
#print(files)

def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix-1),(xpix-1)])
    return vals
'''Rough estimate: 0.75 arcsec/15pixels >>> 0.05 arsec/pixel >>> 2.424e-7rad/pixel
Taken from DS9 tradiation region'''    
def solid_angle(angle):
    return np.pi*np.sin(angle)**2*u.sr
    
def KtoJ(T):
    return (3/2)*k*T
    
def Ntot_rj_thin_nobg(nu,linewidth,s,g,q,eu_J,T_ex,vint_intensity):
    #nu=nu
    #T_ex=T_ex
    #T_r=T_r
    return ((3*k)/(8*np.pi**3*nu*mu_a**2*s*R_i))*(q/g)*np.exp((eu_J/(k*T_ex)))*((vint_intensity)#((nu+templatewidth)-(nu-templatewidth)))

def N_u(ntot,qrot,gu,eu_J,T_ex):
    return ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def Q_rot_asym(T):
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
fluxes=fluxvalues(383,649,files)*u.Jy*u.Hz#/u.sr
howmanybeams=2.424e-7/8.57915480931599e-5#Omega=solid_angle(8.57915480931599e-5*u.deg)#BMAJ
jyhz=fluxes*howmanybeams
#print(jyhz.to('Jy Hz'))
nohzflux=(jyhz[7]*c/lines[0]).to('Jy km s-1')
#print(nohzflux)
#print(lines)#vint_intensities.to('erg s-1 cm-2 sr-1 Hz-1 km s-1'))
vint_trads=nohzflux*((c)**2/(2*k*lines[0]**2))
vint_trads=vint_trads.to('K km s-1')
#print(vint_trads)
qrot=Q_rot_asym(tex[0])
s_j=(4**2-2**2)/(4*(2*4+1))
ntot=Ntot_rj_thin_nobg(lines[0],linewidth,s,9,KtoJ(tex[0])


