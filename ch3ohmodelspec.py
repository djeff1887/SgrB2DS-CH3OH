# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:51:30 2020

@author: zen83
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cnst
import astropy.units as u
from astroquery.splatalogue import Splatalogue as splat

'''for the 10_2-9_3 state'''
testrestfreq=232.41852100*u.GHz

dnu_spw2=488310.9724731*u.Hz

templatewidth=8*dnu_spw2
templatewidthkms=lprime.hztokms(templatewidth, testrestfreq)

freq_min=215*u.GHz
freq_max=235*u.GHz
linelist='JPL'

tau=0.001

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

slaim = splat.query_lines(freq_min, freq_max, chemical_name=' CH3OH',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=[linelist],
                                show_upper_degeneracy=True)
#print(tbl.keys())
frqs=slaim['Freq-GHz(rest frame,redshifted)']*u.GHz
aij=slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
eu_K=slaim['E_U (K)']*u.K
degs=slaim['Upper State Degeneracy']
qns=slaim['Resolved QNs']

def qngrabber(nums,state):
    temp=nums[state].split('(')
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
    
def Q_rot_asym(T):
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)

def Ntot_Nu(qrot,gu,E_u,T):
    val=(qrot/gu)*np.exp(E_u/(k*T))
    return val

def KtoJ(T):
    return (3/2)*k*T

def Ntot(nu,gu,s,q,E,T_ex,T_b):
    nu=nu
    T_ex=T_ex
    T_b=T_b
    val=((3*k)/(8*np.pi*nu**2*s*R_i))*(q/gu)*np.exp(E/(k*T_ex))
    return val

def Ntot_rj_thin_nobg(nu,linewidth,s,g,q,eu_J,T_ex,T_r):
    #nu=nu
    #T_ex=T_ex
    #T_r=T_r
    return ((3*k)/(8*np.pi**3*nu*mu_a**2*s*R_i))*(q/g)*np.exp((eu_J/(k*T_ex)))*((T_r/f)*linewidth)#((nu+templatewidth)-(nu-templatewidth)))

def N_u(ntot,qrot,gu,eu_J,T_ex):
    return ntot/(qrot*np.exp(eu_J/(k*T_ex)))

def T_brightness(nu,Bnu):
    return (((c**2/(2*k*nu**2))*Bnu).to('K'))

def Brj_nu(nu,T):
    return ((2*(nu)**2*(k)*(T)/(c**2)))#*(u.W/(u.m**2*u.Hz*u.sr))

def Jrj_nu(nu,Bnu):
    return (c**2/(2*k*nu**2))*Bnu

#tbl=splat.query_lines(min_frequency=231*u.GHz,max_frequency=235*u.GHz,chemical_name='CH3OH',energy_max=1840,energy_type='eu_k')
eus=[]
nugs=[]

g_k=1
g_i=1
mu_a=(0.896e-18*u.statC*u.cm).to('cm(3/2) g(1/2) s-1 cm')

testTex=100*u.K

brj=Brj_nu(testrestfreq.to('Hz'),testTex)
tr=T_brightness(testrestfreq,brj)
qrot1=Q_rot_asym(testTex)

def main():
    for i in range(len(degs)):
        statenum=i
        J,K=qngrabber(qns,statenum)

        #g_j=2*J+1
        s_j=(J**2-K**2)/(J*(2*J+1))
        #g_u=g_j*g_k*g_i

        testE=eu_K[statenum]
        eus.append(testE/u.K)
        testg=degs[statenum]
        
        testn_tot=Ntot_rj_thin_nobg(testrestfreq,templatewidthkms,s_j,testg,qrot1,KtoJ(testE),testTex,tr)
        n_u_overg=N_u(testn_tot,qrot1,testg,KtoJ(testE),testTex)
        nugs.append((n_u_overg.to('cm-2'))*u.cm**2)
main()
#nugs=nugs*(u.cm)**2
plt.scatter(np.log10(eus),np.log10(nugs))
plt.title(r'Full CH$_3$OH Rotational Diagram @ '+str(testTex))
plt.xlim(xmin=1.5,xmax=3.3)
plt.ylim(ymax=14.5)
plt.xlabel(r'log$_{10}$(E$_u$) (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g) (cm$^{-2}$)')
plt.show()

#spw0=np.ones((1,2))

for i in range(len(frqs)):
    if frqs[i]>=216.468*u.GHz:
        if frqs[i]<=218.343*u.GHz:
            plt.scatter(np.log10(eus[i]),np.log10(nugs[i]),c='red')
    if frqs[i]>=218.303*u.GHz:
        if frqs[i]<=220.176*u.GHz:
            plt.scatter(np.log10(eus[i]),np.log10(nugs[i]),c='blue')
    if frqs[i]>=230.417*u.GHz:
        if frqs[i]<=232.292*u.GHz:
            plt.scatter(np.log10(eus[i]),np.log10(nugs[i]),c='green')
        
#plt.scatter(spw0[0],spw0[1],c='red')
plt.title(linelist+r' CH$_3$OH Rotational Diagram @ '+str(testTex))
plt.xlim(xmin=1.5,xmax=3.3)
plt.ylim(ymax=14.5)
plt.xlabel(r'log$_{10}$(E$_u$) (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g) (cm$^{-2}$)')
plt.show()
#nratio=Ntot_Nu(qrot1,g_u,KtoJ(testE),testT)
#plt.plot(qrot1[1:49],nratio[1:49])
#print(nratio.to(''))
