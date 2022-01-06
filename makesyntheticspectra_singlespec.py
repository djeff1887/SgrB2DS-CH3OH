
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
Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

mpl.interactive(True)

plt.close('all')
linelist='JPL'

'''Collect constants for N_tot and N_upper calculations'''

source='DSVIII'
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

trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':150*u.K,'DSIX':150*u.K}
testT=trotdict[source]
qrot_partfunc=Q_rot_asym(testT).to('')


R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

dopplershifts={'SgrB2S':0.000234806,'DSi':0.000186431,'DSii':0.00015954965399894244,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001661546432045067,'DSIX':0.00015453732389175085}#:0.000190713}/old doppler S: 0.0002306756533745274/0.00015954965399894244/0.00016236367659115043

z=dopplershifts[source]
z_vel=z*c

sourcelocs={'SgrB2S':'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/SgrB2S/OctReimage_K/*.fits','DSi':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSi/field10originals_K/*.fits",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiii/field10originals/*.fits",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiv/field10originals_K/*.fits",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSv/field10originals_K/*.fits",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSVI/field2originals_K/*.fits",'DSIX':"/blue/adamginsburg/d.jeff/XCLASS2021/files/DSIX/field7originals_K/"}
dopplershiftimg={'SgrB2S':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/4_2-3_1vt=0repline_mom1.fits",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial5carryover_field10errors/8_1-7_0vt=0repline_mom1.fits",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/field10originals_noexclusions/mom1/CH3OH~8_0-7_1E1vt0.fits",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/field10originals_noexclusions/mom1/CH3OH~10_2--9_3-vt0.fits"}

texlocs={'SgrB2S':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSi/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSii/Kfield10originals_noexclusions/','DSiii':r'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiii/Kfield10originals_noexclusions/','DSiv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSiv/Kfield10originals_noexclusions/','DSv':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field10/CH3OH/DSv/Kfield10originals_noexclusions_include4-3_trial1/','DSVI':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field2/CH3OH/DSVI/Kfield2originals_trial2_16_6-16_7excluded/','DSVII':"/blue/adamginsburg/d.jeff/SgrB2DSreorg/field3/CH3OH/DSVII/Kfield3originals_trial1_noexclusions/",'DSIX':'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field7/CH3OH/DSIX/Kfield7originals_150K_trial1_noexclusions/'}

texmappath=texlocs[source]+'texmap_3sigma_allspw_withnans_weighted.fits'

texmapdata=fits.getdata(texmappath)*u.K

pixdict={'SgrB2S':(61,64),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSIX':(34,35)}

targetpix=pixdict[source]

testT=texmapdata[targetpix[0],targetpix[1]]

sourcepath=sourcelocs[source]

print(f'Collecting spectra from {sourcepath}')
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


specieslist=[' CH3OH ',' CH3OCHO ',' HOONO ',' HNCO ',' DCN ',' H2CO ',' C2H5OH ',' CH3CHO ',' CH3COOH ',' CH3NH2 ', ' CH3OCH3 ', ' HC3N ',' NH2CHO ', ' NH2CN ',' NH2D ',' SO2 ',' SO ',' t-HCOOH ',' a-H2CCHOH ',' s-H2CCHOH ',' H2CCNH ',' CH3CH2CHO ',' HCCCHO ',' SiS ',' CH2DOH ',' C(18)O ',' HDO ',' CH2CHCN ',' CH3CH2CN ',' c-H2COCH2 ', ' c-HCCCH ',' CCS ',' CH2NH ',"Ethylene Glycol",' cis-CH2OHCHO ','Acetone', ' CH3CN ',' CH2CHCN ']

#molnum=len(specieslist)-38#38
molnum=specieslist.index(' CH3OCHO ')

maskcolumnissuemols=['OCS', 'SiO','C(18)O']

print(f'Species: {specieslist[molnum]}')

testntot=3e10*u.cm**-2

sgrb2scolumns={' CH3OH ':5e16*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':2e15*u.cm**-2, ' SO ':1.25e12*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':1e14*u.cm**-2,' C(18)O ':1e18*u.cm**-2, ' c-HCCCH ':2.5e15*u.cm**-2, 'Acetone':6e13*u.cm**-2}

dsicolumns={' CH3OH ':1e16*u.cm**-2,' CH3OCHO ':3e14*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':0.7e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':2e14*u.cm**-2,' CH3COOH ':0.7e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':6e13*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':0.8e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':8e14*u.cm**-2, ' SO ':7e11*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':5e13*u.cm**-2,' C(18)O ':5e17*u.cm**-2,' CH2DOH ':5e15*u.cm**-2, ' CH2NH ':1e13*u.cm**-2,"Ethylene Glycol":1e13*u.cm**-2,'Acetone':0.25e13*u.cm**-2}

ds6columns={' CH3OH ':1e16*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2}

ds7columns={' CH3OH ':1e16*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2}

ds9columns={' CH3OH ':1e15*u.cm**-2,' CH3OCHO ':2e13*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2}

sourcecolumns={'SgrB2S':sgrb2scolumns,'DSi':dsicolumns,'DSVI':ds6columns,'DSVII':ds7columns,'DSIX':ds9columns}

columndict=sourcecolumns[source]
if specieslist[molnum] in columndict:
    testntot=columndict[specieslist[molnum]]
else:
    pass

methntot=columndict[' CH3OH ']
for spectrum, img, stdimage in zip(spectra,images,stds):
    print('Getting ready - '+img)
    spec=np.genfromtxt(spectrum)
    error=fits.getdata(stdimage)[targetpix[0],targetpix[1]]*u.K
    '''
    refpix={'SgrB2S':[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]],'DSi':[[0,0,0],[266.8316149,-28.3972040,0]],'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]]}
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    targetypix=int(round(targetpixcrd[1][1]))
    targetxpix=int(round(targetpixcrd[1][0]))
    
    if targetworldcrd != refpix[source]:
        print('Non-reference pixel detected\nAdjusting dopplershift\n')
        vel=z*c
        dif=zconvimg[targetypix,targetxpix]
        velcorr=vel+dif
        z=velcorr/c
    '''
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
    
    #linewidth=0.00485*u.GHz#Half of original 0.0097GHz
    #lw2=linewidth/8
    originallinewidth=(11231152.36688232*u.Hz/2)
    
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
            
    '''Generate species table for contaminant search'''    
    species_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=specieslist[molnum], energy_max=1840, energy_type='eu_k', line_lists=['JPL','CDMS','SLAIM','ToyoMA','OSU','RFI','Lisa'], show_upper_degeneracy=True)
    
    if len(species_table['Chemical Name']) == 0:
        print(f'No transitions in {img}. Continue')
        continue
    else:
        pass
    '''
    for rawfreq in species_table['Freq-GHz(rest frame,redshifted)']:
        row=0
        skipmergefreqs=False
        if np.isnan(float(rawfreq)):
            print('Nan detected in query table\nSkipping merge_frequencies\n')
            minchemtable=utils.minimize_table(species_table,merge=False)
            clines=((minchemtable['MeasFreqGHz']*10**9)/(1+z)*u.Hz).to('MHz')
            #skipmergefreqs=True
            break
            #species_table.remove_row(row)
            #row+=1
        else:
            row+=1
            pass
        
    if skipmergefreqs:
    '''
    minchemtable=utils.minimize_table(species_table)
    clines=((minchemtable['Freq']*10**9)/(1+z)*u.Hz).to('MHz')
    
    cqns=minchemtable['QNs']
    ceuks=minchemtable['EU_K']*u.K
    ceujs=[]
    for euk in ceuks:
        ceujs.append(KtoJ(euk))
    cdegs=species_table['Upper State Degeneracy']
    clog10aijs=minchemtable['log10_Aij']
    caijs=10**clog10aijs*u.s**-1
    
    #zeros=np.zeros(cube.shape[0])
    #baseline=np.array(list(zip(plot.value,zeros)))#np.vstack((plot,zeros))
    
    baseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    #mbaseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    baseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    #mbaseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    modelspec=baseline
    methmodelspec=baseline
    print('Begin model loops')
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],np.shape(spec)[0]).to('MHz')
    #modelgaus=models.Gaussian1D(mean=freqs[0], stddev=11 * u.MHz, amplitude=0*u.K)
    #methmodelgaus=models.Gaussian1D(mean=freqs[0], stddev=11 * u.MHz, amplitude=0*u.K)
    linedetections=[]
    for line,deg,euj,aij,qn in zip(clines,cdegs,ceujs,caijs,cqns):
        print(f'Transition: {qn} @ {line.to("GHz")}')
        restline=line*(1+z)
        
        est_nupper=nupper_estimated(testntot,deg,qrot_partfunc,euj,testT).to('cm-2')
        est_tau=opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
        print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        if trad >= 3*error:
            print(f'Estimated brightness: {"{:.3f}".format(trad)}')
            modlinewidth=velocitytofreq(linewidth,line)
            print(f'Model linewidth (Hz): {modlinewidth}')
            modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
            #modelgaus+=modelline
            modelspec+=modelline
            linedetections.append(True)
        else:
            print(f'{qn} line brightness ({trad}) below 3sigma threshold ({3*error})')
            linedetections.append(False)
            continue
        
    print('\nBegin CH3OH modeling')
    mdetections=[]
    for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
        print(f'Transition: {qn} @ {line.to("GHz")}')
        restline=line*(1+z)
        
        est_nupper=nupper_estimated(methntot,deg,qrot_partfunc,euj,testT).to('cm-2')
        est_tau=opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
        print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        if trad >= 3*error:
            print(f'Estimated brightness: {"{:.3f}".format(trad)}')
            modlinewidth=velocitytofreq(linewidth,line)
            print(f'Model linewidth (Hz): {modlinewidth}')
            modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
            #modelgaus+=modelline
            methmodelspec+=modelline
            mdetections.append(True)
        else:
            print(f'{qn} line brightness ({"{:.3f}".format(trad)}) below 3sigma threshold ({3*error})')
            mdetections.append(False)
            continue

    compositespec=modelspec+methmodelspec
    print('Plotting model spectra')
    #plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],cube.shape[0])
    #cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].quicklook()
    plt.rcParams['figure.dpi'] = 150
    plt.figure(1, figsize=(20,10))
    plt.plot(freqs.value,data.value,drawstyle='steps-mid',color='blue')
    plt.plot(freqs,compositespec(freqs),drawstyle='steps-mid',color='green')
    plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color='orange')
    #plt.plot(plot,cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))].value,color='black',drawstyle='steps')
    plt.xlabel(r'$\nu$ (Hz)')
    plt.ylabel('T$_b$ (K)')
    plt.title(f'source: {source}, modelTex: {testT}, modelntot: {testntot}, {specieslist[molnum]}')
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
    
    plt.show()
