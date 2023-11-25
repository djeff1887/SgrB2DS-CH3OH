import astropy.units as u
import scipy.constants as cnst
import numpy as np
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
import os

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
'''
Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=150*u.GHz,
                                                         fmax=300*u.GHz)
'''
'''Hot core locations'''

fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}

dataversion='pacman_sep2023revolution'#homedict[source].replace('/','')
datadir=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/{dataversion}/'

sourcedict={'SgrB2S':'/sep2023-5removelasttorsionalline/','DSi':'/sep2023-5addvt2linesbackin/','DSii':'/sep2023-2widerrefslab/','DSiii':'/sep2023-3vt2doublet/','DSiv':'/sep2023-4nextinline/','DSv':f'/sep2023phi_nu&doublet/','DSVI':'/sep2023-2removenewvt1line/','DSVII':f'/sep2023phi_nu&doublet/','DSVIII':f'/sep2023phi_nu&doublet/','DSIX':f'/sep2023phi_nu&doublet/'}
excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1','15_6-15_7E1vt1','9_6-9_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','8_6-8_7E1vt1','16_6-16_7E1vt1','10_6-10_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','23_5-22_6E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','15_6-15_7E1vt1','16_6-16_7E1vt1','9_6-9_7E1vt1','10_6-10_7E1vt1','11_6-11_7E1vt1','12_6-12_7E1vt1','13_6-13_7E1vt1'],'DSii':['7_6-7_7E1vt1','9_6-9_7E1vt1','14_6-14_7E1vt1','10_6-10_7E1vt1','13_6-13_7E1vt1','11_6-11_7E1vt1','23_5-22_6E1vt0'],'DSiii':'','DSiv':['8_6-8_7E1vt1','7_6-7_7E1vt1','9_6-9_7E1vt1','10_6-10_7E1vt1','11_6-11_7E1vt1','12_6-12_7E1vt1','13_6-13_7E1vt1','14_6-14_7E1vt1','6_1--7_2-vt1'],'DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1','17_6-17_7E1vt1'],'DSVII':["6_1--7_2-vt1"],'DSVIII':'','DSIX':''}

targetworldcrds={'SgrB2S':[[0,0,0],[266.8351718,-28.3961210, 0]], 'DSi':[[0,0,0],[266.8316149,-28.3972040,0]], 'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]],'DSVII':[[0,0,0],[266.8426074, -28.4094401,0]],'DSVIII':[[0,0,0],[266.8418408, -28.4118242, 0]],'DSIX':[[0,0,0],[266.8477371, -28.4311386,0]],'DS10':[[0,0,0],[266.8373798, -28.4009340,0]],'DS11':[[0,0,0],[266.8374572, -28.3996894, 0]],'DSX':[[0,0,0],[266.8452950, -28.4282608,0]]}
pixdict={'SgrB2S':(73,54),'DSi':(36,40),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}
tblconversion={'SgrB2S':'SgrB2S','DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9'}
#sourcedict={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/Kfield10originals_noexclusions/",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_150K_trial2/",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial3_8_6-8_7excluded/",'DSVII':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_200K_trial1_noexclusions/','DSVIII':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVIII/Kfield3originals_175K_trial1_noexclusions/','DSIX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/Kfield7originals_150K_trial1_noexclusions/'}

dopplershifts={'SgrB2S':0.00023099669803283718,'DSi':0.00018761288466593936,'DSii':0.00016236367659115043,'DSiii':0.000176,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016236727257136008,'DSVIII':0.0001661546432045067,'DSIX':0.00015787296484373237}

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
    return nuppercalc.to('cm-2'),nuppererr.to('cm-2')#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def KtoJ(T):#Convert from excitation temperature (Kelvin) to energy (Joules)
    return k*T

def JtoK(E):
    return (E/k).to('K')
    
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
    string=string.replace(',','&')
    return string
    
def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def Ctau(tau):
    return tau/(1-np.exp(-tau))
