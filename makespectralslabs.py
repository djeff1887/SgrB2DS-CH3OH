import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import os
from astropy.modeling import models, fitting
import time
import pdb
import pickle
from astropy.wcs import WCS
import matplotlib as mpl
import copy
from astropy import coordinates
from spectral_cube import BooleanArrayMask

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

'''Collect constants and CH3OH-specific quantum parameters'''
print('Setting constants')
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
f=1
Tbg=2.7355*u.K
nu_bg=282*u.GHz

'''Loop through a given list of lines (in Hz), computing and saving moment0 maps of the entered data cube'''
def linelooplte(line_list,line_width,iterations,quantum_numbers):
    print('\ncubelooperLTE...')
    print('Grab cube and reference pixel')
    targetspec_K=cube[:,pixycrd,pixxcrd]
    #targetpixjybeam=targetpixjybeam.mask_channels(np.isnan(targetpixjybeam)==False)
    cubebeams=(cube.beams.value)*u.sr/u.beam
    #print('Convert from Jy/beam to K')
    #targetspec_K=targetpixjybeam.to(u.K)#JybeamtoK(cubebeams,targetpixjybeam)
    print('Compute cube brightness temperature stddev')
    targetspecK_stddev=stddata[stdpixycrd,stdpixxcrd]#np.nanstd(targetspec_K)
    transitionbeamlist=[]
    transitionfluxlist=[]
    for i in range(iterations):
        print(f'Start {quantum_numbers[i]} moment0 procedure')
        temptransdict={}
        line=line_list[i]#*u.Hz
        restline=line*(1+z)
        nu_upper=line+line_width
        nu_lower=line-line_width
        if nu_upper > max(cube.spectral_axis):
            print(f'Spectral slab exceeds upper bound of spectral axis ({nu_upper} > {max(cube.spectral_axis)})')
            print(f'Skipping {quantum_numbers[i]}\n')
            continue
        if nu_lower < min(cube.spectral_axis):
            print(f'Spectral slab exceeds lower bound of spectral axis ({nu_lower} < {min(cube.spectral_axis)})')
            print(f'Skipping {quantum_numbers[i]}\n')
            continue
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        slabstart=time.time()
        slab=cube.spectral_slab(nu_upper,nu_lower)
        oldstyleslab=cube.spectral_slab((nu_upper-nu_offset),(nu_lower+nu_offset))
        slabend=time.time()-slabstart
        print(f'{quantum_numbers[i]} spectral slab done in {time.strftime("%H:%M:%S", time.gmtime(slabend))}')
        slabbeams=(slab.beams.value)*u.sr/u.beam
        #print(f'slabbeams: {slabbeams}')
        slab_K=slab[:,pixycrd,pixxcrd]#JybeamtoK(slabbeams,)
        #print(f'slab_K: {slab_K}')
        mulu2=(mulu(aijs[i],restline)).to('cm5 g s-2')
        linewidth_vel=vradio(singlecmpntwidth,line)
        tbthick=Tbthick(testntot,restline,linewidth_vel,mulu2,degeneracies[i],qrot_partfunc,eujs[i],testT).to('K')
        peak_amplitude=slab_K.max(axis=0)
        
        est_nupper=nupper_estimated(testntot,degeneracies[i],qrot_partfunc,eujs[i],testT).to('cm-2')
        est_tau=opticaldepth(aijs[i],restline,testT,est_nupper,originallinewidth).to('')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        print('LTE params calculated')
        print(f'tbthick: {tbthick}\n targetspecK_stddev: {targetspecK_stddev}\n peak_amplitude: {peak_amplitude}')
        print(f'est_nupper: {est_nupper}\n est_tau: {est_tau}\n trad: {trad}')
        pdb.set_trace()
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        moment0filename=home+f'{chem}~'+transition+'.fits'
        maskedmom0fn=home+f'{chem}~'+transition+'_masked.fits'
        maskresidualfn=home+f'{chem}~'+transition+'_residual.fits'
        slabfilename=slabpath+f'{chem}~'+transition+'_slab.fits'
        #maskedslabfn=slabpath+'CH3OH~'+transition+'_maskedslab.fits'
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(maskedmom0fn):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(maskedmom0fn)*u.K*u.km/u.s
            #isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(maskedmom0fn)
            isfilestdflux=fits.getdata(f'{stdpath}{images[imgnum]}fluxstd.fits')*u.K#This is confusing, notation-wise, but I'm leaving it this way for now since it's consistent between the two forks in the loop. For future reference: isfilestdflux is the error on the measured brightnesstemp in K, whereas isfilemom0 pulls from the moment0 maps and is in K km/s
            temptransdict.update([('freq',line),('flux',isfilemom0),('stddev',isfilestdflux),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename)])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            masterdegens.append(degeneracies[i])
            print('\nDictionaries populated for this transition.')
            if os.path.isfile(slabfilename):
                print('Proceeding...\n')
                pass
            else:
                slab.write(slabfilename)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            for moment in [1,2]:
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                momentnfilename=sourcepath+f'mom{moment}/'+'CH3OH~'+transition+'.fits'
                if os.path.isfile(momentnfilename):
                    print(f'{transition} moment{moment} file already exists.')
                    continue
                elif moment == 1:
                    print(f'Computing moment 1 and saving to {momentnfilename}\n')
                    slabmom1=slab.moment1()
                    slabmom1.write(momentnfilename)
                elif moment == 2:
                    print(f'Computing moment 2 and saving to {momentnfilename}\n')
                    slabmom2=slab.moment2()
                    slabmom2.write(momentnfilename)
            pass
        elif tbthick >= targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:#*u.K:
            print('Commence moment0')
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])#spwrestfreq)
            #cubemask=BooleanArrayMask(mask=cubemaskarray,wcs=slab.wcs)
            #pdb.set_trace()
            oldstyleslab=oldstyleslab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
            #maskedslab=slab.with_mask(cubemask)
            momstart=time.time()
            print('Unmasked moment0 computing...')
            slabmom0=oldstyleslab.moment0()
            #print('Masked moment0 computing...')
            #maskslabmom0=maskedslab.moment0()
            momend=time.time()-momstart
            print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
            
            #print('\nComputing masking residuals')
            #mom0maskresiduals=maskslabmom0-slabmom0
            print('\nSaving...')
            #name='test'+str(i)
            slabmom0.write((moment0filename),overwrite=True)
            #maskslabmom0.write((maskedmom0fn))
            #mom0maskresiduals.write((maskresidualfn))
            
            moment0beam=slabmom0.beam.value*u.sr
            targetpixflux=slabmom0[pixycrd,pixxcrd]
            temptransdict.update([('freq',line),('flux',slabmom0),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename)])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            masterdegens.append(degeneracies[i])
            print(f'{quantum_numbers[i]} calculations complete.\n')
            if os.path.isfile(slabfilename):
                print(f'Spectral slab {slabfilename} already exists.\nProceeding...\n')
                pass
            else:
                slab.write(slabfilename)
                print(f'Slab written to {slabfilename}.')
                #maskedslab.write(maskedslabfn)
                #print(f'Masked slab written to {maskedslabfn}. Proceeding...\n')
            for moment in [1,2]:
                momentnfilename=sourcepath+f'mom{moment}/CH3OH~'+transition+'.fits'
                if moment == 1:
                    print(f'Computing moment 1 and saving to {momentnfilename}\n')
                    slabmom1=slab.moment1()
                    slabmom1.write(momentnfilename,overwrite=True)
                elif moment == 2:
                    print(f'Computing moment 2 and saving to {momentnfilename}\n')
                    slabmom2=slab.moment2()
                    slabmom2.write(momentnfilename,overwrite=True)
            pass
        else:
            if not tbthick >= targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:
                print('LTE Model max brightnessTemp below 1sigma threshold')
                print(f'{quantum_numbers[i]} skipped, possible contamination\n')
                pass
            elif tbthick >= targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print(f'Line amplitude ({peak_amplitude}) less than 3 sigma criterion ({3*targetspecK_stddev})')
                print(f'{quantum_numbers[i]} skipped\n')
            elif not tbthick >= targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print('1 sigma LTE model and 3 sigma amplitude criteria not met')
                print(f'{quantum_numbers[i]} skipped\n')
                pass
    spectraKdict.update({images[imgnum]:targetspec_K})
    print('lines looped.\n')
    
def KtoJ(T):
    return (3/2)*k*T
    
'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    string=string.replace(',','&')
    return string
    
'''Converts given line list in frequency to radio velocity'''
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('km s-1')

print('Setting input LTE parameters')
testT=300*u.K
testtau=0.2
testntot=1e17*u.cm**-2
print(f'Input Tex: {testT}\nInput Ntot: {testntot}')

def Tbthick(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):
    print(f'ntot: {ntot}, nu: {nu},line_width: {line_width},mulu_2: {mulu_2},g: {g},q: {q},eu_J: {eu_J},T_ex: {T_ex}')
    return (1-np.exp(((-8*np.pi**3*mulu_2*R_i*g)/(3*h*q*line_width))*((np.exp((h*nu)/(k*T_ex))-1)/np.exp((eu_J)/(k*T_ex)))*ntot))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))
    
def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)

def mulu(aij,nu):#Rearranged from Eq 11 (Magnum & Shirley 2015), returns product in units of cm5 g s-2
    return (3*h*c**3*aij)/(64*np.pi**4*nu**3)
    
def rjequivtemp(nu,T_ex):
    return ((h*nu)/k)/(np.exp((h*nu)/(k*T_ex))-1)
    
def t_rad(tau_nu, ff, nu, T_ex):
    return ff*(1-np.exp(-tau_nu))*(rjequivtemp(nu, T_ex)-rjequivtemp(nu,Tbg))
    
def nupper_estimated(n_tot,g,q,euj,tex):
    return n_tot*(g/q)*np.exp(-euj/(k*tex))
    
def opticaldepth(aij,nu,T_ex,nupper,lw):
    return (c**2/(8*np.pi*nu**2*lw))*aij*nupper*np.exp((h*nu)/(k*T_ex))
    
qrot_partfunc=Q_rot_asym(testT).to('')

source='SgrB2S'
chem='H2CO'
fnum=1

dopplershifts={'SgrB2S':0.0002306756533745274,'DSi':0.000186431,'DSv':0.000186431}#:0.000190713}

z=dopplershifts[source]

outpath=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/OctReimage_K/'

incubes=glob.glob(outpath+"*pbcor_line.fits")

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in incubes:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'
    
sourcelocs={'SgrB2S':'K_OctReimage_z0_0002306756533745274_5-6mhzwidth_stdfixes/','DSi':'field10originals_z0_000186431_5-6mhzwidth_stdfixes/','DSv':f'{int(testT.value)}K_field10originals_z0_00186431_5-6mhzwidth_stdfixes_test/'}
stdhome='/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/'


sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/{chem}/{source}/'+sourcelocs[source]
nupperpath=sourcepath+'nuppers/'
stdpath=sourcepath+'errorimgs/std/'
slabpath=sourcepath+'spectralslabs/km_s/'
mom0path=sourcepath+'mom0/'
rotdiagpath=sourcepath+'pixelwiserotationaldiagrams/'
figpath=sourcepath+'figures/'
picklepath=sourcepath+f'{chem}linesdict.obj'

spwdict={}
kstddict={}
kkmsstddict={}
spectraKdict={}

masterlines=[]
masterqns=[]
mastereuks=[]
mastereujs=[]
masterdegens=[]
masterlog10aijs=[]
masteraijs=[]
masterslicedqns=[]
masterrestfreqs=[]

masterfluxes=[]
masterbeams=[]
masterstddevs=[]

if os.path.isdir(slabpath):
    print(f'Source path directory tree {sourcepath} already exists.\n')
    if os.path.isdir(sourcepath+'mom1/'):
        print('Moment 1/2 directories already exist.')
    else:
        for moment in [1,2]:
            momnpath=sourcepath+f'mom{moment}/'
            print(f'Creating moment {moment} directory at {momnpath}')
            os.mkdir(momnpath)
        
    pass
else:
    print(f'Making source path {sourcepath}')
    os.makedirs(sourcepath)
    print(f'Making nupper folder {nupperpath}')
    os.mkdir(nupperpath)
    print(f'Making error folder {stdpath}')
    os.makedirs(stdpath)
    print(f'Making spectral slab folder {slabpath}\n')
    os.makedirs(slabpath)
    for moment in [0,1,2]:
            momnpath=sourcepath+f'mom{moment}/'
            print(f'Creating moment {moment} directory at {momnpath}')
            os.mkdir(momnpath)
    #print(f'Making mom0 folder {mom0path}')
    #os.mkdir(mom0path)
    print(f'Making rotational diagram folder')
    os.mkdir(rotdiagpath)
    print(f'Making figures folder')
    os.mkdir(figpath)
    
for imgnum in range(len(datacubes)):
    print(f'Accessing data cube {datacubes[imgnum]}')
    assert images[imgnum] in datacubes[imgnum], f'{images[imgnum]} not in filename {datacubes[imgnum]}'
    home=sourcepath+'mom0/'#f'{images[imgnum]}/'#Make sure to include slash after path
    readstart=time.time()
    cube=sc.read(datacubes[imgnum])
    readelapsed=time.time()-readstart
    print(f'Cube read in {time.strftime("%H:%M:%S", time.gmtime(readelapsed))}')
    
    #cube=cube.rechunk(save_to_tmp_dir=True)
    header=fits.getheader(datacubes[imgnum])
    
    stdimage=fits.open(stdhome+images[imgnum]+'minimize.image.pbcor_noise.fits')
    stddata=stdimage[0].data*u.K
    stdwcs=WCS(stdimage[0].header)
    
    print('Acquiring cube rest frequency and computing target pixel coordinates')
    spwrestfreq=header['RESTFRQ']*u.Hz
    masterrestfreqs.append(spwrestfreq)
    
    freqs=cube.spectral_axis#Hz
    freqflip=False
    if freqs[1] < freqs[0]:
        freqs=freqs[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
        
    velcube=cube.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
    #print(velcube.spectral_axis)
    cube_unmasked=velcube.unmasked_data
    
    targetworldcrds={'SgrB2S':[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]], 'DSi':[[0,0,0],[266.8316149,-28.3972040,0]], 'DSv':[[0,0,0],[266.8321311,-28.3976633,0]]}
    cube_w=cube.wcs
    #targetworldcrd=[[0,0,0],[266.8324225,-28.3954419,0]]#DSiv
    targetworldcrd=targetworldcrds[source]#[[0,0,0],[266.8316149,-28.3972040,0]] #DSi
    #targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]] #SgrB2S
    #[[0,0,0],[266.8332569, -28.3969, 0]] #DSii/iii
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    fullsize_targetpixcrd=stdwcs.wcs_world2pix(targetworldcrd,1,ra_dec_order=True)
    stdpixxcrd,stdpixycrd=int(round(fullsize_targetpixcrd[1][0])),int(round(fullsize_targetpixcrd[1][1]))
    print(f'Stddev position - x: {stdpixxcrd}/y: {stdpixycrd}')
    
    assert stdpixxcrd >= 0 and stdpixycrd >= 0, 'Negative std pixel coords'
    
    
    pixxcrd,pixycrd=int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1]))
    print(f'Flux position - x: {pixxcrd}/y: {pixycrd}')
    
    assert pixxcrd >= 0 and pixycrd >= 0, 'Negative pixel coords'
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    #print(freq_max)
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Inverted spectral axis'
    print('Passed increasing spectral axis check')
    #print(freq_min)
    linelists=['JPL','SLAIM','CDMS']
    
    print('Peforming Splatalogue queries')
    
    '''Needed for upper state degeneracies'''                                
    for database in linelists:
        sparetable=Splatalogue.query_lines(freq_min, freq_max, chemical_name=f' {chem} ', 
            energy_max=1840, energy_type='eu_k',line_lists=[database],show_upper_degeneracy=True)
        if len(sparetable['Freq-GHz(rest frame,redshifted)'])==0:
            print(f'No {chem} lines found in {database} catalogue')
            continue
        else:
            print(f'{chem} lines identified in {database} catalogue')
            if str(sparetable['Freq-GHz(rest frame,redshifted)'][0])=='--':
                print(f'{database} catalogue does not contain theoretical values. Switching to SLAIM')
                sparetable=Splatalogue.query_lines(freq_min, freq_max, chemical_name=f' {chem} ', 
            energy_max=1840, energy_type='eu_k',line_lists=['SLAIM'],show_upper_degeneracy=True)
                break
            break
                                
    maintable = sparetable
    '''
    utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=f' {chem} ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=[database],
                                    show_upper_degeneracy=True))              
    '''
    #pdb.set_trace()
    print('Gathering Splatalogue table parameters')    
    lines=maintable['Freq-GHz(rest frame,redshifted)']*u.GHz/(1+z)#Redshifted to source
    #masterlines.append(lines)
    #vel_lines=vradio(lines,spw1restfreq)
    qns=maintable['Resolved QNs']
    euks=maintable['E_U (K)']*u.K
    eujs=[]
    for eupper_K in euks:
        eujs.append(KtoJ(eupper_K))
    degeneracies=sparetable['Upper State Degeneracy']
    log10aijs=maintable['Log<sub>10</sub> (A<sub>ij</sub>)']
    aijs=10**log10aijs*u.Hz
    
    '''
    for i in range(len(test)):
        plt.axvline(x=test[i],color='red')
    plt.show()
    '''
    singlecmpntwidth=(0.00485/8)*u.GHz
    linewidth=(15.15*u.MHz)
    originallinewidth=(11231152.36688232*u.Hz/2)#0.005*u.GHz####0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
    nu_offset=linewidth-originallinewidth
    linewidth_vel=vradio(singlecmpntwidth,spwrestfreq)#(singlecmpntwidth*c.to(u.km/u.s)/spwrestfreq).to('km s-1')#vradio(linewidth,spw1restfreq)
    #slicedqns=[]
    
    pixeldict={}
    transitiondict={}
    linelooplte(lines,linewidth,len(lines),qns)
    spwdict.update([(images[imgnum],transitiondict)])
    tempkeys=list(spwdict[images[imgnum]].keys())
    
    '''
    kstdimgpath=stdpath+f'{images[imgnum]}fluxstd.fits'
    kkmsstdimgpath=stdpath+f'{images[imgnum]}intensitystd.fits'
    if os.path.isfile(kkmsstdimgpath):
        print(f'{images[imgnum]} brightness std image already exists')
        spwstdarray=fits.getdata(kstdimgpath)*u.K
        kkmsstdarray=fits.getdata(kkmsstdimgpath)*u.K*u.km/u.s
        print(f'Retrieved brightness std data from {kstdimgpath} and {kkmsstdimgpath}\n')
    else:
        print(f'Start {images[imgnum]} std calculations')
        spwstdarray,kkmsstdarray=pixelwisestd(cube)
        for stdarray, imgpath in zip([spwstdarray,kkmsstdarray],[kstdimgpath,kkmsstdimgpath]):
            print('Set Primary HDU')
            hdu=fits.PrimaryHDU(stdarray.value)
            This transmoment0 file has intensity (K km/s) units
            if len(tempkeys) == 0:
                print(f'No transitions detected in this spw ({images[imgnum]})')
                transmomslab=cube.spectral_slab((lines[0]-linewidth),(lines[0]+linewidth))
                transmoment0=transmomslab.moment0()
                transmom0header=transmoment0.header
                print(f'Set transmoment0 to moment0 from {(lines[0]+linewidth).to("GHz")} to {(lines[0]-linewidth).to("GHz")}')
            else:
                transmoment0=fits.open(spwdict[images[imgnum]][tempkeys[0]]['filename'])
                transmom0header=transmoment0[0].header
                print(f'Set header from {spwdict[images[imgnum]][tempkeys[0]]["filename"]}')
            hdu.header=transmom0header
            if np.all(stdarray==spwstdarray):
                hdu.header['BUNIT']='K'
            else:
                hdu.header['BUNIT']='K km s-1'
            print('Wrapping Primary HDU in HDUList')
            hdul=fits.HDUList([hdu])
            print(f'Writing to {imgpath}')
            hdul.writeto(imgpath)
        print(f'{images[imgnum]} std calculations complete.\n')
    
    #transitiondict.update({'restfreq':spwrestfreq})
    #,('pixel_0',(pixycrd,pixxcrd))])
    kstddict.update([(images[imgnum],spwstdarray)])
    kkmsstddict.update([(images[imgnum],kkmsstdarray)])
    '''
    print(f'Finished loop for {images[imgnum]}')

if os.path.isfile(picklepath):
    print(f'pickle {picklepath} already exists.')
else:
    print('Saving dictionary pickle...')
    f=open(picklepath,'wb')
    pickle.dump(spwdict,f)
    f.close()
    print(f'Dictionary pickle saved at {picklepath}')

print('Saving EU_K/QN/etc txt file')
eukqns=np.column_stack((mastereuks,masterqns,masterlines,masterdegens))
np.savetxt(sourcepath+'mastereuksqnsfreqsdegens.txt',eukqns,fmt='%s',header=f'{chem} transitions, excitation temperatures, and degeneracies used in this folder. Temperatures in units of K, frequencies are redshifted ({z}/{(z*c).to("km s-1")}) and in Hz.')
