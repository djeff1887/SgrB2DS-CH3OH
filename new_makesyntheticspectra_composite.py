
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
from utilities import *#Q_rot_asym,mulu,vradio,t_rad,nupper_estimated,opticaldepth,qngrabber
import matplotlib as mpl
#from cubecoretexmap import t_rad,nupper_estimated,opticaldepth
import pdb
from pyspeckit.spectrum.models import lte_molecule

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

mpl.interactive(True)

plt.close('all')
linelist='JPL'

def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

'''Collect constants for N_tot and N_upper calculations'''

source='SgrB2S'

sourceisnew=True

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K


trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':215*u.K,'DSIX':150*u.K}

testT=trotdict[source]
qrot_partfunc=Q_rot_asym(testT).to('')


R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1


dopplershifts={'SgrB2S':0.000228,'DSi':0.0001865,'DSii':0.000163,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001662062062062062,'DSIX':0.00015453732389175085}#:0.000190713}/old doppler S: 0.0002306756533745274/0.00015954965399894244/0.00016236367659115043

s_othermol_dshift_v={' CH3CHO ':67.45330305*u.km/u.s,' C2H5OH ':67.45330305*u.km/u.s,' CH3OCHO ':67.45330305*u.km/u.s,' C(18)O ':69.551850256*u.km/u.s,' 13CH3OH ':67.5*u.km/u.s,' SO ':70.5*u.km/u.s}#' CH3OH ':68352.680424
ds2_othermol_dshift_v={' CH3OCHO ':49*u.km/u.s,' CH3CHO ':49*u.km/u.s,' C2H5OH ':49.3*u.km/u.s}#47831.782945392486 m / s
othermol_dopplershift={' CH3CHO ':0.000225,' C2H5OH ':0.000225,' CH3OCHO ':0.000225,' C(18)O ':0.000232}

sourceothers={'SgrB2S':s_othermol_dshift_v,'DSii':ds2_othermol_dshift_v}
othermol_dshift_v=sourceothers[source]

z=dopplershifts[source]
z_vel=z*c

sourcelocs={'SgrB2S': r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/SgrB2S/OctReimage_K','DSi':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSi/field10originals_K','DSii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSii/field10originals_K','DSiii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiii/field10originals_K','DSiv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiv/field10originals_K','DSv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSv/field10originals_K','DSVI':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVI/field2originals_K','DSVII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVII/field3originals_K','DSIX':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSIX/field7originals_K'}

texlocs={'SgrB2S':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS1/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS2/Kfield10originals_noexclusions/','DSiii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DSiii/Kfield10originals_noexclusions/','DSiv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DSiv/Kfield10originals_noexclusions/','DSv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS5/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS6/DSVI/Kfield2originals_trial2_16_6-16_7excluded/','DSVII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS7/Kfield3originals_trial1_noexclusions/','DSIX':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS9/Kfield7originals_150K_trial1_noexclusions/'}

arabicswitch={'DSiii':'DS3','DSiv':'DS4','DSv':'DS5'}

if source in arabicswitch.keys():
    texlocs[source]=texlocs[source].replace(source,arabicswitch[source])
    texmappath=texlocs[source]+'texmap_3sigma_allspw_withnans_weighted.fits'
else:
    texmappath=texlocs[source]+'texmap_3sigma_allspw_withnans_weighted.fits'

texmapdata=fits.getdata(texmappath)*u.K

pixdict={'SgrB2S':(70,59),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSIX':(34,35)}#61,64

targetpix=pixdict[source]

testT=texmapdata[targetpix[0],targetpix[1]]

sourcepath=sourcelocs[source]

print(f'Collecting spectra from {sourcepath}')
incubes=glob.glob(sourcepath)
inspecs=glob.glob(sourcepath+'/*.txt')
instds=glob.glob(texlocs[source]+'errorimgs/std/*.fits')
images=['spw0','spw1','spw2','spw3']

spectra=[]

stds=[]

for spew in images:
    for f1 in inspecs:
        if spew in f1:
            if str(targetpix[0]) in f1 and str(targetpix[1]) in f1:
                spectra.append(f1)
                continue
            else:
                continue
    for f2 in instds:
        if spew in f2:
            stds.append(f2)
            continue
        else:
            continue
    
assert 'spw0' in spectra[0] and 'spw0' in stds[0], 'List out of order'

print('Spectra and stds are sequential order')

#imgnum=0
testline=0

linewidth=2.5*u.km/u.s#2.5 km/s is ideal for DSVI
print(f'Absolute model line width: {linewidth}\n')

specieslist=[' CH3OH ',' CH3OCHO ',' HOONO ',' HNCO ',' DCN ',' H2CO ',' C2H5OH ',' CH3CHO ',' CH3COOH ',' CH3NH2 ', ' CH3OCH3 ', ' HC3N ',' NH2CHO ', ' NH2CN ',' NH2D ',' SO2 ',' SO ',' t-HCOOH ',' a-H2CCHOH ',' s-H2CCHOH ',' H2CCNH ',' CH3CH2CHO ',' HCCCHO ',' SiS ',' CH2DOH ',' C(18)O ',' HDO ',' CH2CHCN ',' CH3CH2CN ',' c-H2COCH2 ', ' c-HCCCH ',' CCS ',' CH2NH ',"Ethylene Glycol",' cis-CH2OHCHO ','Acetone', ' CH3CN ',' CH2CHCN ',' CH3O13CHO ', ' SiO ', ' OCS ', 'N2D+', ' CH3C(15)N ',' CH3CCH ',' CH3SH ',' 13CS ', ' H2S ', ' SO ',' CH3(18)OH ', ' 13CH3OH ']

speciesname=' CH3OCHO '
molnum=specieslist.index(speciesname)#42

#maskcolumnissuemols=['OCS', 'SiO','C(18)O']

testntot=1e14*u.cm**-2

sgrb2scolumns={' CH3OH ':2e17*u.cm**-2,' CH3OCHO ':7e15*u.cm**-2, ' CH3CHO ':5e15*u.cm**-2,' C2H5OH ':6e16*u.cm**-2,' CH3OCH3 ':1.5e15*u.cm**-2,' DCN ':5e15*u.cm**-2, ' OCS ':8e16*u.cm**-2,' 13CH3OH ':1.5e17*u.cm**-2,' H2CO ':5e16*u.cm**-2,' HC3N ':1.8e15*u.cm**-2, ' C(18)O ':1.3e19*u.cm**-2,' 13CS ':8e15*u.cm**-2,' SO2 ':3e16*u.cm**-2,' NH2CHO ':5e13*u.cm**-2,' HNCO ':2e16*u.cm**-2,' SO ':2e13*u.cm**-2,' SiO ':1e15*u.cm**-2,' H2S ':1e18*u.cm**-2,' H2CCO ':1e16*u.cm**-2,}#' H2CS ':1e18*u.cm**-2,' CH3(18)OH ':2.5e16*u.cm**-2,' CH3COOH ':2e15*u.cm**-2,,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' c-HCCCH ':2.5e15*u.cm**-2, 'Acetone':6e13*u.cm**-2,' CH3C(15)N ':3e13*u.cm**-2,' SiN ':2e15*u.cm**-2, ' CH3NH2 ':9e15*u.cm**-2,}#' HOONO ':5e15*u.cm**-2,

dsicolumns={' CH3OH ':1e17*u.cm**-2,' CH3OCHO ':3e14*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':0.7e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':2e14*u.cm**-2,' CH3COOH ':0.7e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':6e13*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':0.8e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':8e14*u.cm**-2, ' SO ':7e11*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':5e13*u.cm**-2,' C(18)O ':5e17*u.cm**-2,' CH2DOH ':5e15*u.cm**-2, ' CH2NH ':1e13*u.cm**-2,"Ethylene Glycol":1e13*u.cm**-2,'Acetone':0.25e13*u.cm**-2, ' SiO ':7e13*u.cm**-2, ' OCS ':5.25e15*u.cm**-2, ' 13CS ':5e14*u.cm**-2, ' 13CH3OH ':2e15*u.cm**-2}

ds2columns={' CH3OH ':6e16*u.cm**-2,' CH3OCHO ':5e14*u.cm**-2,' CH3CHO ':6e14*u.cm**-2,' C2H5OH ':7e15*u.cm**-2,' CH3OCH3 ':2.8e14*u.cm**-2,' DCN ':7e14*u.cm**-2,' OCS ':3.5e16*u.cm**-2,' 13CH3OH ':1.5e16*u.cm**-2,' H2CO ':4e16*u.cm**-2,' HC3N ':5e14*u.cm**-2,' C(18)O ':5e18*u.cm**-2,' 13CS ':2e15*u.cm**-2,' SO2 ':2.5e15*u.cm**-2,' NH2CHO ':1e13*u.cm**-2,' HNCO ':4e15*u.cm**-2,' SO ':4e12*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e17*u.cm**-2,}#' CH3(18)OH ':2.5e15*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' (13)CN ':1.5e15*u.cm**-2  NH2CN ':1e13*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,  ' SiS ':4e15*u.cm**-2,}' HOONO ':1e16*u.cm**-2,'Carbon Dioxide':5e16*u.cm**-2,' NaCl ':1e16*u.cm**-2,' CH3COOH ':7e14*u.cm**-2,' HDCO ':1e18*u.cm**-2,

ds5columns={' CH3OH ':1e15*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1.75e15*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e13*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':0.3e14*u.cm**-2}


ds6columns={' CH3OH ':1e17*u.cm**-2,' CH3OCHO ':5e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1.3e16*u.cm**-2,' DCN ':3e15*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':5e16*u.cm**-2, ' CH3CHO ':3e15*u.cm**-2,' CH3COOH ':5e15*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':1e13*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' SiO ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' 13CH3OH ':1e17*u.cm**-2}

ds7columns={' CH3OH ':1e16*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':1e14*u.cm**-2}

ds9columns={' CH3OH ':1e15*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':1e14*u.cm**-2}

weeds=[' CH3OCHO ', ' CH3CHO ']

sourcecolumns={'SgrB2S':sgrb2scolumns,'DSi':dsicolumns, 'DSii':ds2columns,'DSv':ds5columns,'DSVI':ds6columns,'DSVII':ds7columns,'DSIX':ds9columns}

columndict=sourcecolumns[source]
if specieslist[molnum] in columndict:
    testntot=columndict[specieslist[molnum]]
else:
    pass

methntot=columndict[' CH3OH ']
plt.rcParams['figure.dpi'] = 150
plt.figure(1, figsize=(20,10))
molcolors=['red','cyan','orange','brown','deepskyblue','darkviolet','yellow','pink','darkviolet','darkkhaki','silver','blue','lime','magenta','grey','plum','fuchsia','darkcyan']
spwmoldict={}
dummylist=[]
firstmolline={}#list(np.ones(len(columndict.keys())))
for m in columndict.keys():
    firstmolline.update({m:1})

for spectrum, img, stdimage in zip(spectra,images,stds):
    print('Getting ready - '+img)
    spec=np.genfromtxt(spectrum)
    error=fits.getdata(stdimage)[targetpix[0],targetpix[1]]*u.K

    freqs=spec[:,0]*u.MHz#cube.spectral_axis
    data=spec[:,1]*u.K
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        data=data[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    print('Plotting model spectra')
    plt.plot(freqs.value,data.value,drawstyle='steps-mid',color='black')
    
    '''Generate methanol table for use during contaminant search'''
    methanol_table=Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=['JPL','CDMS','SLAIM','ToyoMA','OSU','RFI','Lisa'], show_upper_degeneracy=True)
    minmethtable=utils.minimize_table(methanol_table)
    mlines=((minmethtable['Freq']*10**9)/(1+z)*u.Hz).to('MHz')
    mqns=minmethtable['QNs']
    meuks=minmethtable['EU_K']*u.K
    meujs=[]
    for euk in meuks:
        meujs.append(KtoJ(euk))
    mdegs=methanol_table['Upper State Degeneracy']
    mlog10aijs=minmethtable['log10_Aij']
    maijs=10**mlog10aijs*u.s**-1
    
    '''Create background model for the methanol lines and other species'''
    baseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    #mbaseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    baseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    #mbaseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    #modelspec=baseline
    methmodelspec=baseline
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],np.shape(spec)[0]).to('MHz')
    modeldict={}
       
    for molecule,hue,first in zip(list(columndict.keys())[1:], molcolors,firstmolline):     
        '''Generate species table for contaminant search'''
        modelspec=baseline
        species_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=molecule, energy_max=1840, energy_type='eu_k', line_lists=['JPL','CDMS','SLAIM','ToyoMA','OSU','RFI','Lisa'], show_upper_degeneracy=True)
        if len(species_table['Chemical Name']) == 0:
            print(f'No transitions for {molecule} in {img}. Continue')
            continue
        else:
            pass
        
        minchemtable=utils.minimize_table(species_table)
        if molecule in othermol_dshift_v.keys():
            otherz=othermol_dshift_v[molecule]/c
            clines=((minchemtable['Freq']*10**9)/(1+otherz)*u.Hz).to('MHz')
        else:
            clines=((minchemtable['Freq']*10**9)/(1+z)*u.Hz).to('MHz')
        
        cqns=minchemtable['QNs']
        ceuks=minchemtable['EU_K']*u.K
        ceujs=[]
        for euk in ceuks:
            ceujs.append(KtoJ(euk))
        cdegs=species_table['Upper State Degeneracy']
        clog10aijs=minchemtable['log10_Aij']
        caijs=10**clog10aijs*u.s**-1
        cntot=columndict[molecule]
        print(f'Begin model loops for {molecule}')
        
        linedetections=[]
        for line,deg,euj,aij,qn in zip(clines,cdegs,ceujs,caijs,cqns):
            #print(f'Transition: {qn} @ {line.to("GHz")}')
            if molecule in othermol_dshift_v.keys():
                restline=line*(1+otherz)
            else:
                restline=line*(1+z)
            est_nupper=nupper_estimated(cntot,deg,qrot_partfunc,euj,testT).to('cm-2')
            modlinewidth=velocitytofreq(linewidth,line)
            phi_nu=lineprofile(modlinewidth,restline,restline)
            intertau=lte_molecule.line_tau(testT, cntot, qrot_partfunc, deg, restline, euj, aij) #opticaldepth(aij,restline,testT,est_nupper,modlinewidth).to('')
            est_tau=(intertau*phi_nu).to('')
            #print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
            trad=t_rad(f,est_tau,restline,testT).to('K')
            if trad >= 3*error:
                #print(f'Estimated brightness: {"{:.3f}".format(trad)}')
                #modlinewidth=velocitytofreq(linewidth,line)
                #print(f'Model linewidth (Hz): {modlinewidth}')
                modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
                #modelgaus+=modelline
                modelspec+=modelline
                linedetections.append(True)
            else:
                #print(f'{qn} line brightness ({trad}) below 3sigma threshold ({3*error})')
                linedetections.append(False)
                continue
        if molecule == ' CH3CHO ':
            dummylist.append((freqs,modelspec(freqs)))#spwmoldict.update({img:(freqs,modelspec(freqs))})
        if firstmolline[first]:
            plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color=hue,label=molecule)
            firstmolline[first]=0
        else:
            plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color=hue)
    
    print('Begin CH3OH modeling\n')
    mdetections=[]
    for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
        #print(f'Transition: {qn} @ {line.to("GHz")}')
        restline=line*(1+z)
        modlinewidth=velocitytofreq(linewidth,line)
        phi_nu=lineprofile(modlinewidth,restline,restline)
        
        est_nupper=nupper_estimated(methntot,deg,qrot_partfunc,euj,testT).to('cm-2')
        intertau=lte_molecule.line_tau(testT, methntot, qrot_partfunc, deg, restline, euj, aij)#opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
        est_tau=(intertau*phi_nu).to('')
        #print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        if trad >= 3*error:
            #print(f'Estimated brightness: {"{:.3f}".format(trad)}')
            #print(f'Model linewidth (Hz): {modlinewidth}')
            modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
            #modelgaus+=modelline
            methmodelspec+=modelline
            mdetections.append(True)
        else:
            #print(f'{qn} line brightness ({"{:.3f}".format(trad)}) below 3sigma threshold ({3*error})')
            mdetections.append(False)
            continue

    #compositespec=modelspec+methmodelspec
    plt.plot(freqs,methmodelspec(freqs),drawstyle='steps-mid',color='green',linestyle='--')#plt.plot(freqs,compositespec(freqs),drawstyle='steps-mid',color='green')
    #plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color='orange')
    #plt.plot(plot,cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].value,color='black',drawstyle='steps')
    '''
    print('Overplotting axvlines and transition annotations')
    for line,qn,detected in zip(clines,cqns,linedetections):
        if detected:
            plt.axvline(x=line.value,linestyle='--',color='yellow',ymin=0.25)
            plt.annotate(qn, (line.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    
    for mline,mqn,detected in zip(mlines,mqns,mdetections):
        if detected:
            plt.axvline(x=mline.value,linestyle='--',color='pink',ymin=0.25)
            plt.annotate(mqn, (mline.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    '''
plt.xlabel(r'$\nu$ (Hz)')
plt.ylabel('T$_b$ (K)')
plt.ylim(ymax=100)
#plt.title(f'source: {source}, modelTex: {testT}, modelntot: {testntot}, {specieslist[molnum]}')
plt.legend()
plt.show()
