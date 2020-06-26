import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam

plt.close('all')
linelist='JPL'

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
mu_a=(0.896e-18*u.statC*u.cm).to('cm(3/2) g(1/2) s-1 cm')

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

files=glob.glob('/ufrc/adamginsburg/d.jeff/imaging_results/*.fits')
z=0.0002333587
imgnames=['spw2','spw1','spw0']



def gauss(x,A,mu,sig):
    return A*np.exp((-1/2)*((x-mu)/sig)**2)

def Q_rot_asym(T):
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('cm s-1')
    
def KtoJ(T):
    return (3/2)*k*T
    
def qngrabber(nums):
    temp=nums.split('(')
    temp2=temp[1].split(',')
    jupper=int(temp[0])
    if linelist == 'JPL':
        temp3=temp2[0].split(')')
        kupper=temp3[0]
        if 'a' in kupper:#What are these things?
            kupper=0
        else:
            kupper=int(temp3[0])
    else:
        kupper=int(temp2[0])
    
    return jupper, kupper
    
def Ntot_rj_thin_nobg(nu,line_width,s,g,q,eu_J,T_ex,T_r):
    #nu=nu
    #T_ex=T_ex
    #T_r=T_r
    return ((3*k)/(8*np.pi**3*nu*mu_a**2*s*R_i))*(q/g)*np.exp((eu_J/(k*T_ex)))*((T_r/f)*line_width)
    
def Tb2(nu,line_width,s,n_u):
    return (8*np.pi**3*nu*mu_a**2*s*R_i*f*n_u)/(3*k*line_width)

def Tb3(ntot,nu,line_width,s,g,q,eu_J,T_ex):
    return ((8*np.pi**3*nu*mu_a**2*s*R_i*g*f)/(3*k*q*np.exp(eu_J/(k*T_ex))*line_width))*ntot
    
def N_u(ntot,qrot,gu,eu_J,T_ex):
    return ntot/((qrot/gu)*np.exp(eu_J/(k*T_ex)))
    
def t_brightness(eu_J,gu,qrot,ntot,n_u):
    return (eu_J/k)*(((ntot*gu)/(qrot*n_u))-1)**-1

def kkms(beams,data):
    intensitylist=[]
    t_bright=[]
    for i in range(len(data)):
        temp=(data[i]).to('Jy/beam')
        #print(temp)
        equiv=u.brightness_temperature(data.spectral_axis[i])
        #print(equiv)
        jy_sr=temp/beams[i]
        #print(jy_sr)
        conversion=jy_sr.to(u.K,equivalencies=equiv)
        t_bright.append(conversion.value)
        #print(conversion)
        #velflux_T=conversion*lwvel
        #print(velflux_T)
        #print('\n')
        #intensitylist.append(velflux_T)
    return t_bright

print('Getting ready - '+imgnames[0])
cube=sc.read(files[0])
#header=fits.getheader(files[0])
#beamer=radio_beam.Beam.from_fits_header(hdu).value

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

assert freq_max > freq_min, 'Decreasing frequency axis'


linewidth=0.00485*u.GHz#Half of original 0.0097GHz
lw2=linewidth/4

        
'''Generate methanol table for contaminant search'''    
methanol_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=['JPL'], show_upper_degeneracy=True)

minmethtable=utils.minimize_table(methanol_table)

mlines=(minmethtable['Freq']*10**9)/(1+z)*u.Hz
mqns=minmethtable['QNs']
meuks=minmethtable['EU_K']*u.K
meujs=[]
for euk in meuks:
    meujs.append(KtoJ(euk))
mdegs=methanol_table['Upper State Degeneracy']

testline=8
plotwidth=linewidth*1.5
lwvel=vradio(lw2,mlines[testline])
print(f'Transition: {mqns[testline]}\nEU_K: {meuks[testline]}')

spwwindow=cube.spectral_slab((mlines[testline]-plotwidth),(mlines[testline]+plotwidth))[:,649,383]
beamlist=spwwindow.beams
beamlist=(beamlist.value)*u.sr/u.beam
t_brights=kkms(beamlist,spwwindow)

#print(t_brights)


Tphys=np.linspace(1,1000,100)*u.K
testT=165*u.K
n_totes=[1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23]*u.cm**-2
plot=np.linspace((mlines[testline]-plotwidth),(mlines[testline]+plotwidth),30)

b1=[]
b2=[]
qs=[]

q=Q_rot_asym(testT).to('')
J,K=qngrabber(mqns[testline])
s_j=(J**2-K**2)/(J*(2*J+1))
n_total=1e18*u.cm**-2
n_upper=N_u(n_total,q,mdegs[testline],meujs[testline],meuks[testline]).to('cm-2')
#testtbright=t_brightness(meujs[testline],mdegs[testline],q,n_total,n_upper)
testtbright2=Tb3(n_total,mlines[testline],lwvel,s_j,mdegs[testline],q,meujs[testline],meuks[testline]).to('K')#(mlines[testline],lwvel,s_j,n_upper).to('K')

plotprofile=[]
for i in range(len(plot)):
    temp=gauss(plot[i],testtbright2,mlines[testline],lw2)
    plotprofile.append(temp/u.K)

'''
b1.append(testtbright/u.K)
b2.append(testtbright2/u.K)
qs.append(q)
'''
        
print(f'q: {q} n_upper: {n_upper} nu/g: {n_upper/mdegs[testline]} Tb2: {testtbright2}')
plt.plot(spwwindow.spectral_axis,t_brights,drawstyle='steps')
plt.plot(plot,plotprofile)
plt.axvline(x=mlines[testline].value,ls='--')
plt.title(f'Transition: {mqns[testline]} EU_K: {meuks[testline]} Tphys: {testT}')
#plt.plot(Tphys.value,b1)
#plt.plot(Tphys.value,b2)
#plt.legend()
plt.show()

#plt.plot(Tphys.value,qs)
#plt.show()
