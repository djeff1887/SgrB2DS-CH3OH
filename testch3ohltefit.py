import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from scipy.optimize import curve_fit as cf
from astropy.io import fits
import glob
import radio_beam

from astropy.modeling import models, fitting#Fittable1DModel, Parameter, fitting

plt.close('all')
linelist='JPL'

'''Collect constants for N_tot and N_upper calculations'''

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

files=glob.glob('/blue/adamginsburg/d.jeff/imaging_results/*.fits')
z=0.000234806#<<<avg of the two to the left#0.000236254#0.0002333587
imgnames=['spw3','spw0','spw2','spw1']

def gauss(x,A,mu,sig):#Standard Gaussian equation
    return A*np.exp((-1/2)*((x-mu)/sig)**2)

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
def mulu(aij,nu):#Rearranged from Eq 11 (Magnum & Shirley 2015), returns product in units of cm5 g s-2
    return (3*h*c**3*aij)/(64*np.pi**4*nu**3)
    
def vradio(frequency,rest_freq):#Converts between frequency to radio velocity (used to convert the model line width to velocity for column density calculations
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('cm s-1')
    
def component_restfrequency(cmpntvel,rest_freq):#Computes new rest velocities of new velocity components for methanol
    ref_freq=minmethtable['Freq'][testline]*10**9*u.Hz
    c_kms=c.to('km s-1')
    z_vel=z*c_kms
    cmpntshift_vel=z_vel+cmpntvel
    cmpnt_zvel=cmpntshift_vel/c_kms
    cmpnt_restfreq=ref_freq/(1+cmpnt_zvel)
    
    return cmpnt_restfreq
    
def KtoJ(T):#Convert from excitation temperature (Kelvin) to energy (Joules)
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
    
def rjequivtemp(nu,T_ex):
    return ((h*nu)/k)/(np.exp((h*nu)/(k*T_ex))-1)

def Tb3(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):#Rearranged from Eq 82, M&S 2015
    print('In function variables:')
    print(f'ntot: {ntot} nu: {nu} line_width: {line_width} mulu_2: {mulu_2} g: {g} q: {q} eu_J: {eu_J} Tex: {T_ex}')
    return ((8*np.pi**3*nu*mulu_2*R_i*g*f)/(3*k*q*np.exp(eu_J/(k*T_ex))*line_width))*ntot
    
def Tbthick(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):
    return (1-np.exp(((-8*np.pi**3*mulu_2*R_i*g)/(3*h*q*line_width))*((np.exp((h*nu)/(k*T_ex))-1)/np.exp((eu_J)/(k*T_ex)))*ntot))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))

'''#Used pyspeckit equations to get brightness temperature with optical depth, does not work though   
def Tbthick2(nupper,nu,line_width,aij,eu_J,T_ex):
    return (1-np.exp(((-c**2*aij)/(8*np.pi*nu**2))*(np.exp((h*nu)/(k*T_ex)))*nupper))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))
'''

def opticaldepth(Tr,nu,T_ex):
    return -np.log(1-(Tr/(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))))
    
def opticaldepth2(mulu_2,nu,line_width,T_ex,nupper):
    return ((8*np.pi**3*nu*mulu_2)/(3*h*c*line_width))*(np.exp((h*nu)/(k*T_ex))-1)*nupper
    
def N_u(ntot,qrot,gu,eu_J,T_ex):#Computes upper-state column density for a given molecule (methanol in this case)
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
    
imgnum=3
testline=0
print('Getting ready - '+imgnames[imgnum])
cube=sc.read(files[imgnum])

'''Set WCS params [RA,DEC,FRQ].(For non-contsub cubes not using spectral-cube, the Stokes axis hasn't been removed, so the input array must be [RA,DEC,STOKES,FRQ]. '''
cube_w=cube.wcs#WCS(files[imgnum])
if 'medsub' in files[imgnum]:
    contsub=True
else:
    contsub=False

targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)

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
lw2=linewidth/8

        
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
mlog10aijs=minmethtable['log10_Aij']
maijs=10**mlog10aijs*u.s**-1

plotwidth=linewidth*1.5
lwvel=vradio(lw2,mlines[testline])
print(f'Transition: {mqns[testline]}\nEU_K: {meuks[testline]}')

cmpnt2veloffset=1*u.km/u.s
cmpnt2frq=component_restfrequency(cmpnt2veloffset,mlines[testline])
print(f'2nd component: {cmpnt2frq}')

spwwindow=cube.spectral_slab((mlines[testline]-plotwidth),(mlines[testline]+plotwidth))[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]
beamlist=spwwindow.beams
beamlist=(beamlist.value)*u.sr/u.beam
t_brights=kkms(beamlist,spwwindow)
        
#print(t_brights)
peakK=spwwindow[np.argmax(spwwindow)]

Tphys=np.linspace(1,1000,100)*u.K
testT=550*u.K
n_totes=[1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23]*u.cm**-2
plot=np.linspace((mlines[testline]-plotwidth),(mlines[testline]+plotwidth),30)
#plotcmpnt=

b1=[]
b2=[]
qs=[]

q=Q_rot_asym(testT).to('')
#J,K=qngrabber(mqns[testline])
#s_j=(J**2-K**2)/(J*(2*J+1))
mulu2=mulu(maijs[testline],mlines[testline]).to('cm5 g s-2')#(0.896e-18*u.statC*u.cm).to('cm(3/2) g(1/2) s-1 cm')
n_total=1e17*u.cm**-2
n_upper=N_u(n_total,q,mdegs[testline],meujs[testline],testT).to('cm-2')
print(f'pre-function n_upper: {n_upper}')
#testtbright=t_brightness(meujs[testline],mdegs[testline],q,n_total,n_upper)
testtbright=Tb3(n_total,mlines[testline],lwvel,mulu2,mdegs[testline],q,meujs[testline],testT).to('K')#(mlines[testline],lwvel,s_j,n_upper).to('K')
testtbthick=Tbthick(n_total,mlines[testline],lwvel,mulu2,mdegs[testline],q,meujs[testline],testT).to('K')
testtbcmpnt2=Tbthick(n_total,mlines[testline],lwvel,mulu2,mdegs[testline],q,meujs[testline],testT).to('K')
testtau=opticaldepth(testtbthick,mlines[testline],testT)
#testtau2=opticaldepth(testtbthick2,mlines[testline],testT)
testtau3=opticaldepth2(mulu2,mlines[testline],lw2,testT,n_upper).to('')

plotprofilethin=[]
plotprofilethick=[]
plotprofilecmpnt2=[]
if contsub:
    contoffset=0
else:
    contoffset=20

for i in range(len(plot)):
    temp=gauss(plot[i],testtbright,mlines[testline],lw2)
    temp2=gauss(plot[i],testtbthick,mlines[testline],lw2)
    temp3=gauss(plot[i],testtbthick,cmpnt2frq,lw2)
    plotprofilethin.append((temp/u.K)+contoffset)
    plotprofilethick.append((temp2/u.K)+contoffset)
    plotprofilecmpnt2.append((temp3/u.K)+contoffset)

comboplotprofile=[]
for j in range(len(plot)):
    gaussadd=(plotprofilethick[j]-contoffset)+(plotprofilecmpnt2[j]-contoffset)
    comboplotprofile.append(gaussadd+contoffset)
'''    
class Gaussian1D(Fittable1DModel):
    amplitude = Parameter(default=0, unit=u.K)
    mean = Parameter(default=np.median(spwwindow.spectral_axis).value,unit=u.Hz)
    stddev = Parameter(default=1, unit=u.MHz)

    @staticmethod
    def evaluate(x, amplitude, mean, stddev):
        return (amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2)))#+contoffset*u.K

    @staticmethod
    def fit_deriv(x, amplitude, mean, stddev):
        d_amplitude = np.exp((-(1 / (stddev**2)) * (x - mean)**2))
        d_mean = (2 * amplitude *
                  np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                  (x - mean) / (stddev**2))
        d_stddev = (2 * amplitude *
                    np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                    ((x - mean)**2) / (stddev**3))
        return [d_amplitude, d_mean, d_stddev]
        
'''
        
fit_g=fitting.LevMarLSQFitter()
linemin=0
linemax=29
testclass=models.Gaussian1D(mean=np.median(spwwindow.spectral_axis), stddev=1 * u.MHz, amplitude=1 * u.K)#testtbthick,mlines[testline],lw2)
testfit=fit_g(testclass,spwwindow.spectral_axis[linemin:linemax],t_brights[linemin:linemax]*u.K)
        
print(f'q: {q} n_upper: {n_upper} nu/g: {n_upper/mdegs[testline]} Tb: {testtbright}')
print(f'Tbthick: {testtbthick}')
#print(f'Tbthick-2: {testtbthick2}')
print(f'thin-thick: {testtbright-testtbthick}')
print(f'tau: {testtau}')
#print(f'tau2: {testtau2}')
print(f'tau3: {testtau3}')
print(f'aij: {maijs[testline]} lines: {mlines[testline]}')
print(f'Model fit params: {testfit}')
print(f'Model fit info: {fit_g.fit_info}')
plt.plot(spwwindow.spectral_axis,t_brights,drawstyle='steps')
plt.plot(plot,plotprofilethin,label=(r'$\tau << 1$'))
plt.plot(plot,plotprofilethick,label=(r'$\tau \geq 1$'))
plt.plot(spwwindow.spectral_axis[linemin:linemax],t_brights[linemin:linemax],drawstyle='steps',color='orange')
'''
plt.plot(plot,plotprofilecmpnt2,label=(r'cmpnt2'))
plt.plot(plot,comboplotprofile,label=('combo'))
'''
plt.plot(plot,testfit(plot),label='LMLSQ fit')
plt.axvline(x=mlines[testline].value,ls='--')
plt.title(f'Transition: {mqns[testline]} EU_K: {meuks[testline]} Tphys: {testT}')
#plt.plot(Tphys.value,b1)
#plt.plot(Tphys.value,b2)
plt.legend()
plt.show()

#plt.plot(Tphys.value,qs)
#plt.show()

