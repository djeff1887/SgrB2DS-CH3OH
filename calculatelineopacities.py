from utilities import *
import astropy.units as u
import astropy.constants as cnst
from pyspeckit.spectrum.models import lte_molecule
import numpy as np
from utilities import velocitytofreq, Q_rot_asym
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
from astropy.table import QTable
import sys
import pickle

mpl.interactive(True)
plt.close('all')
plt.rcParams["figure.dpi"]=150

def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def N_u(nu,Aij,velocityintegrated_intensity_K):#,velint_intK_err):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    #nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err
    return nuppercalc#,nuppererr

repodict={'SgrB2S':'/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/Kfield10originals_noexclusions/','DSiii':'/Kfield10originals_noexclusions/','DSiv':'/Kfield10originals_noexclusions/','DSv':f'/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':'/Kfield2originals_trial3_8_6-8_7excluded/','DSVII':'/Kfield3originals_200K_trial1_noexclusions/','DSVIII':'/Kfield3originals_175K_trial1_noexclusions/','DSIX':f'/Kfield7originals_150K_trial1_noexclusions/'}

k=cnst.k_B
h=cnst.h
c=cnst.c
sources=sourcedict.keys()
alltablepaths=glob.glob('OpticalDepthTables/*ntot_4-3peak.fits')#('OpticalDepthTables/*_contpeak_nothiiregion.fits')#('_4-3peak.fits')
testT=[150,300,500]*u.K

for tblpath in alltablepaths:
    s=tblpath.replace('OpticalDepthTables/','')
    s=s.replace('_ntot_4-3peak.fits','')#('_contpeak_nothiiregion.fits','')#('_4-3peak.fits','')
    print(f'Source: {s}')
    qtable=QTable.read(tblpath)
    tau150=[]
    tau300=[]
    tau500=[]
    masternuppers=[]
    masterngs=[]
    cnuppers05=[]
    zipped={}

    for line in qtable:
        peakflux=line['Peak Velocity Integrated Flux']
        errpeakflux=line['Peak Flux Error']
        restfreq=line['Freq']
        euj=line['EU(J)']
        euk=line['EU(K)']
        aij=line['Aij']
        degen=line['g']
        fwhm=line['Line Width']
        peaktb=(peakflux/fwhm).to('K')
        fwhm_Hz=velocitytofreq(fwhm,restfreq)
        
        nupper=N_u(restfreq,aij,peakflux)
        ntot=line['Ntot']#10**(17.75)*u.cm**-2#lte_molecule.ntot_of_nupper(nupper,euj,t,qrot,degen)
        masternuppers.append(nupper.to('cm-2').value)
        masterngs.append(nupper.to('cm-2').value/degen)
        for t in testT:
            qrot=Q_rot_asym(t)
            phi_nu=lineprofile(fwhm_Hz,restfreq,restfreq)
            intertau=lte_molecule.line_tau(t,ntot,qrot,degen,restfreq,euj,aij)
            tau=(intertau*phi_nu).to('')
            
            tauC=0.3
            cnupper05=nupper*(tauC/(1-np.exp(-tauC)))
            cnuppers05.append(cnupper05)
            pctdiff=(cnupper05-nupper)/((cnupper05+nupper)/2)*100
            print(f'QN:{line["QNs"]} Tex:{t} Qrot:{qrot.to("")} Tau:{tau} Ntot:{ntot.to("cm-2")}, Nupper:{nupper.to("cm-2")} Underestimated by: {round(pctdiff.value,2)}%')
            if t == testT[0]:
                tau150.append(tau.value)
            elif t == testT[1]:
                tau300.append(tau.value)
            elif t == testT[2]:
                tau500.append(tau.value)

    savefigbase=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{s}'
    savefighome=savefigbase+repodict[s]
    savefigpath1=savefighome+'qrotfix_lineopacities.png'
    savefigpath2=savefighome+'qrotfix_linear_lineopacities.png'
    savefigpath3=savefighome+'qrotfix_tauvseupper.png'
    savefigpath4=savefighome+'qrotfix_linear_tauvseupper.png'

    plt.figure()
    plt.scatter(masterngs,tau150,c=qtable['EU(K)'].value,label='$Q_{rot}$(150 K)')
    plt.scatter(masterngs,tau300,c=qtable['EU(K)'].value,marker='*',label='$Q_{rot}$(300 K)')
    plt.scatter(masterngs,tau500,c=qtable['EU(K)'].value,marker='s',label='$Q_{rot}$(500 K)')
    plt.xlabel(r'log$_{10}$($N_{u}/g_u$) (cm$^{-2}$)',fontsize=15)
    plt.ylabel(r'$\tau$',fontsize=15)
    #plt.title(f'{s}, Qrot({testT[1]})')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.colorbar(label='$E_U$ (K)')
    plt.tight_layout()
    plt.savefig(savefigpath1)
    plt.show()

    plt.figure()
    plt.scatter(masterngs,tau150,c=qtable['EU(K)'].value,label='$Q_{rot}$(150 K)')
    plt.scatter(masterngs,tau300,c=qtable['EU(K)'].value,marker='*',label='$Q_{rot}$(300 K)')
    plt.scatter(masterngs,tau500,c=qtable['EU(K)'].value,marker='s',label='$Q_{rot}$(500 K)')
    plt.xlabel(r'log$_{10}$($N_{u}/g_u$) (cm$^{-2}$)',fontsize=15)
    plt.ylabel(r'$\tau$',fontsize=15)
    #plt.title(f'{s}, Qrot({testT[1]})')
    plt.legend()
    plt.colorbar(label='$E_U$ (K)')
    plt.tight_layout()
    plt.savefig(savefigpath2)
    plt.show()
    
    plt.figure()
    plt.scatter(qtable['EU(K)'].value,tau150,label='$Q_{rot}$(150 K)')
    plt.scatter(qtable['EU(K)'].value,tau300,marker='*',label='$Q_{rot}$(300 K)')
    plt.scatter(qtable['EU(K)'].value,tau500,marker='s',label='$Q_{rot}$(500 K)')
    plt.xlabel(r'$E_U$ (K)',fontsize=15)
    plt.ylabel(r'$\tau$',fontsize=15)
    #plt.title(f'{s}, Qrot({testT[1]})')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    #plt.colorbar(label='$E_U$ (K)')
    plt.tight_layout()
    plt.savefig(savefigpath3)
    plt.show()

    plt.figure()
    plt.scatter(qtable['EU(K)'].value,tau150,label='$Q_{rot}$(150 K)')
    plt.scatter(qtable['EU(K)'].value,tau300,marker='*',label='$Q_{rot}$(300 K)')
    plt.scatter(qtable['EU(K)'].value,tau500,marker='s',label='$Q_{rot}$(500 K)')
    plt.xlabel(r'$E_U$ (K)',fontsize=15)
    plt.ylabel(r'$\tau$',fontsize=15)
    #plt.title(f'{s}, Qrot({testT[1]})')
    plt.legend()
    #plt.colorbar(label='$E_U$ (K)')
    plt.tight_layout()
    plt.savefig(savefigpath4)
    plt.show()
    
    if s == 'DSi':
        for i,j in zip(qtable['EU(K)'].value,np.squeeze(tau300)):
            zipped.update({i:j})
        myFile=open(f'{s}_tau300.obj','wb')
        pickle.dump(zipped,myFile)
        myFile.close
        print(f'Saved pickle of taus to this folder')
    #sys.exit()
'''    
sgrb2scontsource_43tb=100*u.K#61.9*u.K#peak taken from CARTA
qrot=Q_rot_asym(testT)#3203.211430184744
fourthreeaij=4.686514876992861e-05*u.Hz
restfreq=218.44005*u.GHz
sigma_vel=5*u.km/u.s
sigma_freq=velocitytofreq(sigma_vel,restfreq)
freqrange=np.linspace(218.4*u.GHz,218.5*u.GHz,100)

nupper_pix=N_u(restfreq,fourthreeaij,713*u.K*u.km*u.s**-1)#9.44e13*u.cm**-2
eupper43=45.45988*u.K
eupperJ=(3/2)*k*eupper43
degen=9.0
ntot=lte_molecule.ntot_of_nupper(nupper_pix,eupperJ,sgrb2scontsource_43tb,qrot,degen)


phi_nu=[]
for f in freqrange:
    phi_nu.append(lineprofile(sigma_freq,restfreq,f).value)

#pdb.set_trace()
phi_nu=np.array(phi_nu)*u.MHz**-1

intertau=lte_molecule.line_tau(sgrb2scontsource_43tb,ntot,qrot,degen,restfreq,eupperJ,fourthreeaij)
tau=(intertau*phi_nu).to('')

plt.plot(freqrange,tau)
plt.xlabel('Freq (GHz)')
plt.ylabel(r'{\tau}')
plt.title(f'{testT}, {sigma_vel} line width')
plt.show()

#print(f'Tau = {tau}')
'''
