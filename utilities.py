import astropy.units as u
import scipy.constants as cnst
import numpy as np

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
R_i=1
f=1
Tbg=2.7355*u.K

'''LTE Analysis'''
def Tbthick(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):
    print(f'ntot: {ntot}, nu: {nu},line_width: {line_width},mulu_2: {mulu_2},g: {g},q: {q},eu_J: {eu_J},T_ex: {T_ex}')
    return (1-np.exp(((-8*np.pi**3*mulu_2*R_i*g)/(3*h*q*line_width))*((np.exp((h*nu)/(k*T_ex))-1)/np.exp((eu_J)/(k*T_ex)))*ntot))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
def mulu(aij,nu):#Rearranged from Eq 11 (Magnum & Shirley 2015), returns product in units of cm5 g s-2
    return (3*h*c**3*aij)/(64*np.pi**4*nu**3)

'''Converts given line list in frequency to radio velocity'''
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('km s-1')
    
'''Converts given velocity width to frequency width'''
def velocitytofreq(velocity,ref_freq):
    frequency=((velocity)/c.to(u.km/u.s))*ref_freq
    return frequency.to('MHz')
    
def freqtovelocity(freq,ref_freq):
    velocity=(freq/ref_freq)*c.to('km s-1')
    return velocity.to('km s-1')

'''Compute Rayleigh-Jean equivalent brightness temperature per M&S 2015 Eq. 24'''
def rjequivtemp(nu,T_ex):
    return ((h*nu)/k)/(np.exp((h*nu)/(k*T_ex))-1)
    
'''Compute upper-state column density and associated error of transition of interest'''      
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err#(velint_intK_err.to('K km s-1')/velocityintegrated_intensity_K.to('K km s-1'))*nuppercalc
    return nuppercalc,nuppererr#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
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
        
def t_rad(tau_nu, ff, nu, T_ex):#Radiation temperature per M&S (2015) Eq 27
    return ff*(1-np.exp(-tau_nu))*(rjequivtemp(nu, T_ex)-rjequivtemp(nu,Tbg))
    
def nupper_estimated(n_tot,g,q,euj,tex):#LTE model upper-state column density (nupper) using a fiducial/test total column density (ntot) and excitation temperature (Tex)
    return n_tot*(g/q)*np.exp(-euj/(k*tex))
    
def opticaldepth(aij,nu,T_ex,nupper,lw):#LTE model optical depth using nupper from the above equation
    return (c**2/(8*np.pi*nu**2*lw))*aij*nupper*np.exp((h*nu)/(k*T_ex))
    
'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    return string
    