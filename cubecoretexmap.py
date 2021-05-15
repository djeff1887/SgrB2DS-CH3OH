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
from astropy.nddata import Cutout2D
from spectral_cube.io.casa_masks import make_casa_mask

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

'''This wing of the script takes in continuum-subtracted cubes, cuts out a subcube around a region of interest based on a DS9 region, and converts the subcubes into brightness temperature (K) units'''

print('Begin Jy/beam-to-K and region subcube conversion\n')

#source='DSv'
#source='DSi'
source='SgrB2S'
fnum=1

#inpath="/orange/adamginsburg/sgrb2/d.jeff/data/field10originalimages/"
inpath='/blue/adamginsburg/d.jeff/imaging_results/data/OctReimage/'
beamcubes=glob.glob(inpath+'*.fits')
#home="/orange/adamginsburg/sgrb2/d.jeff/products/field10originalimages/"
home='/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/'
cubes=glob.glob(home+'*pbcor_line.fits')
#region='fk5; box(266.8321311,-28.3976633, 0.0010833, 0.0010833)'#DSv
#region='fk5; box(266.8324225,-28.3954419, 0.0010417, 0.0010417)'#DSiv
#region='fk5; box(266.8316387, -28.3971867, 0.0010556, 0.0010556)'#DSi-large
region='fk5; box(266.8353410,-28.3962005,0.0016806,0.0016806)'#SgrB2S-box2
#region='fk5; box(266.8350804, -28.3963256, 0.0023889, 0.0023889'#SgrB2S-large
#box(266.8333438, -28.3966103, 0.0014028, 0.0014028)' #DSii/iii
#box(266.8315833, -28.3971867, 0.0006528, 0.0006528)' #DSi-small
#outpath=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/OctReimage/'
outpath=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/OctReimage_K/'#imaging_results/DSii_iiibox1/'
#statfixpath=f'/blue/adamginsburg/d.jeff/SgrB2DSstatcontfix/field10originals/'
statfixpath=f'/blue/adamginsburg/d.jeff/SgrB2DSstatcontfix/OctReimage_K/'

if os.path.isdir(outpath):
    print(f'Minicube directory "{outpath}" already exists. Proceeding to line loops/LTE fitting procedure.\n')
    pass
else:
    cubestobox=[]
    
    images=['spw0','spw1','spw2','spw3']
    
    orderedcubes=[]
    orderedbeamcubes=[]
    
    for spew in images:
        for f1,f2 in zip(cubes,beamcubes):
            if spew in f1:
                orderedcubes.append(f1)
            if spew in f2:
                orderedbeamcubes.append(f2)
        #images.append(files[77:81])#[57:61])
        
    assert 'spw0' in orderedcubes[0] and 'spw0' in orderedbeamcubes[0], f'Cube list out of order'
    
    print('Cube lists successfully reordered')
    
    if not os.path.isdir(outpath):
        print(f'Creating filepath {outpath}')
        os.makedirs(outpath)
    else:
        print(f'{outpath} already exists. Proceeding...\n')
        
    if 'products' in home:
        print('***STATCONT products detected***\n')
        if not os.path.isdir(statfixpath):
            print(f'Creating beamfix directory {statfixpath}')
            os.makedirs(statfixpath)
            for beamcube,statcube in zip(orderedbeamcubes,orderedcubes):
                print(f'Extracting beams from {beamcube}')
                beamfits=fits.open(beamcube)
                cubefits=fits.open(statcube)
                beams=beamfits[1]
                cubedata=cubefits[0]
                newhdulist=fits.HDUList([cubedata,beams])
                print(f'Beamlist merged with Primary HDU in {statcube}')
                cubewithbeampath=statcube.replace(home,statfixpath)
                print(f'Saving new fits file {cubewithbeampath}\n')
                newhdulist.writeto(cubewithbeampath)
                cubestobox.append(cubewithbeampath)
        else:
            print(f'{statfixpath} already exists. Grabbing and reordering statfix cubes')
            statfixcubes=glob.glob(statfixpath+'*.fits')
            for spew in images:
                for statcube in statfixcubes:
                    if spew in statcube:
                        cubestobox.append(statcube)
                        continue
            print('Statfix cubes reordered.\n')
            
    else:
        cubestobox=orderedcubes
            
    pdb.set_trace()
            
    for sppw, cub in zip(images, cubestobox):
        boxcubename=outpath+sppw+'minimize.image.pbcor_line.fits'
        if os.path.isfile(boxcubename):
            print(f'{boxcubename} already exists. Skipping...\n')
            continue
        else:
            print(f'Grabbing datacube from {cub}')
            fullsizecube=sc.read(cub,use_dask=True)
            fullsizecube.allow_huge_operations=True
            spwrestfreq=fullsizecube.header['RESTFRQ']*u.Hz
            print('Creating subcube and converting from Jy/beam to K')
            boxedsubcubeK=fullsizecube.subcube_from_ds9region(region).to(u.K)
            #print('Converting spectral axis to km/s')
            #boxedsubcubeKkms=boxedsubcubeK.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
            print(f'Saving to {boxcubename}')
            boxedsubcubeK.write(boxcubename,format='fits',overwrite=True)
            print('Finished\n')
        
       
'''This wing of the code runs the linelooper LTE modeling and kinetic temperature determination on the newly created region-specific subcubes'''

print('Begin core cube to Tex map process\n') 

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

dopplershifts={'SgrB2S':0.000234806,'DSi':0.000186431,'DSv':0.000186431}#:0.000190713}/old doppler S: 0.0002306756533745274

z=dopplershifts[source]
#z=0.00017594380066803095 #SgrB2DSII?
#z=0.000186431 #SgrB2DSi/DSiv(?)
#z=0.0002306756533745274#<<average of 2 components of 5_2-4_1 transition using old redshift(0.000236254)#0.000234806#0.0002333587 SgrB2S
print(f'Doppler shift: {z} / {(z*c).to("km s-1")}\n')

print('Setting input LTE parameters')
testT=300*u.K#500*u.K
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
    
def KtoJ(T):
    return (3/2)*k*T
    
def JybeamtoK(beams,data):
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
    
'''Converts given line list in frequency to radio velocity'''
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('km s-1')
    
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
        print(f'\nStart {quantum_numbers[i]} moment0 procedure')
        temptransdict={}
        line=line_list[i]#*u.Hz
        restline=line*(1+z)
        nu_upper=line+line_width
        nu_lower=line-line_width
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        slabstart=time.time()
        slab=cube.spectral_slab(nu_upper,nu_lower)
        oldstyleslab=cube.spectral_slab((nu_upper-nu_offset),(nu_lower+nu_offset))
        slabend=time.time()-slabstart
        print(f'{quantum_numbers[i]} spectral slab done in {time.strftime("%H:%M:%S", time.gmtime(slabend))}')
        #pdb.set_trace()
        peakchannel=slab.closest_spectral_channel(line)
        print(f'Peak channel: {peakchannel}')
        slabbeams=(slab.beams.value)*u.sr/u.beam
        #print(f'slabbeams: {slabbeams}')
        slab_K=slab[:,pixycrd,pixxcrd]#JybeamtoK(slabbeams,)
        #print(f'slab_K: {slab_K}')
        mulu2=(mulu(aijs[i],line)).to('cm5 g s-2')
        linewidth_vel=vradio(singlecmpntwidth,line)
        tbthick=Tbthick(testntot,restline,linewidth_vel,mulu2,degeneracies[i],qrot_partfunc,eujs[i],testT).to('K')
        peak_amplitude=slab_K[peakchannel]#slab_K.max(axis=0)
        if peak_amplitude == slab_K.max(axis=0):
            print('Peak amplitude == Max amplitude in slab')
        else:
            print('Other bright line in slab')
            print(f'Max brightness in slab: {slab_K.max(axis=0)}\n')
        
        est_nupper=nupper_estimated(testntot,degeneracies[i],qrot_partfunc,eujs[i],testT).to('cm-2')
        est_tau=opticaldepth(aijs[i],restline,testT,est_nupper,originallinewidth).to('')
        trad=t_rad(f,est_tau,restline,testT).to('K')
        
        print('LTE params calculated')
        print(f'tbthick: {tbthick}\n targetspecK_stddev: {targetspecK_stddev}\n peak_amplitude: {peak_amplitude}')
        print(f'est_nupper: {est_nupper}\n est_tau: {est_tau}\n trad: {trad}')
        #if tbthick >= targetspecK_stddev:
        #    print(f'\n est tau from data at {testT}: {(peak_amplitude/trad).to("")}')
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        moment0filename=home+'CH3OH~'+transition+'_raw.fits'
        maskedmom0fn=home+'CH3OH~'+transition+'_masked.fits'
        maskresidualfn=home+'CH3OH~'+transition+'_residual.fits'
        slabfilename=slabpath+'CH3OH~'+transition+'_slab.fits'
        maskedslabfn=slabpath+'CH3OH~'+transition+'_maskedslab.fits'
        maskfn=slabpath+'CH3OH~'+transition+'_mask.fits'
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(maskedmom0fn):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(maskedmom0fn)*u.K*u.km/u.s
            #isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(maskedmom0fn)
            isfilestdflux=stddata#fits.getdata(f'{stdpath}{images[imgnum]}fluxstd.fits')*u.K#This is confusing, notation-wise, but I'm leaving it this way for now since it's consistent between the two forks in the loop. For future reference: isfilestdflux is the error on the measured brightnesstemp in K, whereas isfilemom0 pulls from the moment0 maps and is in K km/s
            temptransdict.update([('freq',restline),('flux',isfilemom0),('stddev',isfilestdflux),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line)])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print('\nDictionaries populated for this transition.')
            if os.path.isfile(maskedslabfn):
                print('Proceeding...\n')
                pass
            else:
                slab.write(maskedslabfn)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            for moment in [1,2]:
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                momentnfilename=sourcepath+f'mom{moment}/'+'CH3OH~'+transition+'.fits'
                if os.path.isfile(momentnfilename):
                    print(f'{transition} moment{moment} file already exists.\n')
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
        elif trad >= 3*targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:#*u.K:
            print('Commence moment0 procedure\n')
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])#spwrestfreq)
            #cubemask=BooleanArrayMask(mask=cubemaskarray,wcs=slab.wcs)
            print(f'Create {quantum_numbers[i]} spatial-velocity mask')
            slabspecax=slab.spectral_axis
            slabmom1=slab.moment1()
            slabfwhm=(7*u.MHz/line)*c.to('km s-1')#slab.linewidth_fwhm()
            cubemask=(slabspecax[:,None,None] < (slabmom1 + slabfwhm)[None,:,:]) & (slabspecax[:,None,None] > (slabmom1 - slabfwhm)[None,:,:])
            oldstyleslab=oldstyleslab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
            #if imgnum > 0:
            #    pdb.set_trace()
            print('Masking spectral slab')
            maskedslab=slab.with_mask(cubemask)
            momstart=time.time()
            print('Unmasked moment0 computing...\n')
            slabmom0=oldstyleslab.moment0()
            print('Masked moment0 computing...\n')
            maskslabmom0=maskedslab.moment0()
            momend=time.time()-momstart
            #print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
            
            print('\nComputing masking residuals')
            mom0maskresiduals=maskslabmom0-slabmom0
            print('\nSaving...')
            #name='test'+str(i)
            slabmom0.write((moment0filename),overwrite=True)
            maskslabmom0.write((maskedmom0fn))
            mom0maskresiduals.write((maskresidualfn))
            #maskboolarr=BooleanArrayMask(mask=cubemask,wcs=slab.wcs)
            #make_casa_mask(maskedslab,maskfn,append_to_image=False,add_stokes=False)
            moment0beam=slabmom0.beam.value*u.sr
            targetpixflux=slabmom0[pixycrd,pixxcrd]
            temptransdict.update([('freq',restline),('flux',maskslabmom0),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line)])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i])
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print(f'{quantum_numbers[i]} calculations complete.\n')
            if os.path.isfile(slabfilename):
                print(f'Spectral slab {slabfilename} already exists.\nProceeding...\n')
                pass
            else:
                slab.write(slabfilename)
                print(f'Slab written to {slabfilename}.')
                maskedslab.write(maskedslabfn)
                print(f'Masked slab written to {maskedslabfn}. Proceeding...\n')
            for moment in [1,2]:
                momentnfilename=sourcepath+f'mom{moment}/CH3OH~'+transition+'.fits'
                if moment == 1:
                    print(f'Computing moment 1 and saving to {momentnfilename}\n')
                    #slabmom1=slab.moment1()
                    slabmom1.write(momentnfilename)
                elif moment == 2:
                    print(f'Computing moment 2 and saving to {momentnfilename}\n')
                    slabmom2=slab.moment2()
                    slabmom2.write(momentnfilename)
            pass
        else:
            if not trad >= 3*targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:
                print('LTE Model max brightnessTemp below 3sigma threshold')
                print(f'{quantum_numbers[i]} skipped, possible contamination\n')
                pass
            elif trad >= 3*targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print(f'Line amplitude ({peak_amplitude}) less than 3 sigma criterion ({3*targetspecK_stddev})')
                print(f'{quantum_numbers[i]} skipped\n')
            elif not trad >= 3*targetspecK_stddev and not peak_amplitude >= 3* targetspecK_stddev:
                print('3 sigma LTE model and 3 sigma amplitude criteria not met')
                print(f'{quantum_numbers[i]} skipped\n')
                pass
    spectraKdict.update({images[imgnum]:targetspec_K})
    print('lines looped.\n')

'''Compute the pixelwise standard deviation for error calculations'''
def pixelwisestd(datacube):
    rowdims=len(cube[1])
    coldims=len(cube[2])
    stdarray=np.empty((rowdims,coldims))
    for row in range(rowdims):
        print(f'Start Row {row} std calcs')
        for col in range(coldims):
            targetpixspecstd=cube[:,row,col].std()
            stdarray[row,col]=targetpixspecstd.value
    print('Compute intensity std')
    intensitystds=(stdarray*u.K)*linewidth_vel
    return stdarray*u.K,intensitystds
    
'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    return string
    
'''Gathers beam data from moment map headers'''    
def beamer(momentmap):
    print(f'in beamer function: {momentmap}')
    hdu=fits.getheader(momentmap)
    momentbeam=radio_beam.Beam.from_fits_header(hdu).value
    return momentbeam*u.sr

'''Reorders Splatalogue table parameters to match the glob.glob filename order'''
def unscrambler(filenames,sliced_qns,linelist):
    #print('Start unscrambler')
    unscrambled_qns=[]
    unscrambled_freqs=[]
    unscrambled_eus=[]
    unscrambled_degs=[]
    unscrambled_aijs=[]
    tempfiles=np.copy(filenames)
    for i in range(len(filenames)):
        print(f'filename: {filenames[i]}')
        tempfiles[i]=tempfiles[i].replace('.fits','')
        for j in range(len(sliced_qns)):
            #print(f'sliced_qns: {sliced_qns[j]}')
            #print(f'comparison qns: {tempfiles[i][55:]}')
            comp=(sliced_qns[j]==tempfiles[i][79:])
            if comp==True:
                print(f'{sliced_qns[j]} == {tempfiles[i][79:]}\comp==True')
                unscrambled_qns.append(sliced_qns[j])
                unscrambled_freqs.append(linelist[j])
                unscrambled_eus.append(mastereuks[j]/u.K)
                unscrambled_degs.append(masterdegens[j])
                unscrambled_aijs.append(masteraijs[j])
                break
            else: 
                print(f'{sliced_qns[j]} != {tempfiles[i][79:]}')
    return unscrambled_qns,unscrambled_freqs,unscrambled_eus,unscrambled_degs,unscrambled_aijs
    
def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix-1),(xpix-1)])#Corrects for different pixel counting procedures
    return vals
    
'''Compute Kkm/s intensity from datadict'''
def JybeamtoKkms(fluxdict):
    intensitylist={}
    t_bright={}
    dictkeys=fluxdict.keys()
    for key in dictkeys:
        temptransdict=fluxdict[key]
        temptransdictkeys=list(temptransdict.keys())
        print(temptransdictkeys)
        for i in range(len(temptransdictkeys)):
            if 'restfreq' in temptransdictkeys[i]:
                continue
            else:
                temp=(temptransdict[temptransdictkeys[i]]['flux']/linewidth_vel).to('Jy')
                #print(temp)
                equiv=u.brightness_temperature(temptransdict[temptransdictkeys[i]]['freq'])
                #print(equiv)
                jy_sr=temp/temptransdict[temptransdictkeys[i]]['beam']
                #print(jy_sr)
                conversion=jy_sr.to(u.K,equivalencies=equiv)
                t_bright.update({temptransdictkeys[i]:conversion})
                #print(conversion)
                velflux_T=conversion*linewidth_vel
                d_velfluxT=(temptransdict[temptransdictkeys[i]]['stddev']/conversion)*velflux_T
                intensityerror.append(d_velfluxT)
                #print(velflux_T)
                #print('\n')
                intensitylist.update({temptransdictkeys[i]:velflux_T})
    return intensitylist,t_bright
    
def brightnessTandintensities(fluxdict):
    intensitydict={}
    t_bright={}
    dictkeys=fluxdict.keys()
    for key in dictkeys:
        temptransdict=fluxdict[key]
        temptransdictkeys=list(temptransdict.keys())
        print(f'Transition keys in brightnessTandintensities: {temptransdictkeys}')
        for i in range(len(temptransdictkeys)):
            if 'restfreq' in temptransdictkeys[i]:
                continue
            else:
                velflux_T=temptransdict[temptransdictkeys[i]]['flux']
                intensitydict.update({temptransdictkeys[i]:velflux_T})
                temp=velflux_T/linewidth_vel
                t_bright.update({temptransdictkeys[i]:temp})
                d_velfluxT=(temptransdict[temptransdictkeys[i]]['stddev'])#/temp)*velflux_T
                intensityerror.append(d_velfluxT)
                
    return intensitydict,t_bright
                
def jupperfinder(quan_nums):
    j_upper=[]
    k_upper=[]
    for i in range(len(quan_nums)):
        for j in range(len(quan_nums[i])):
            comp=quan_nums[i][j].isdigit()
            if comp==False:
                appendage=quan_nums[i][:(j)]
                j_upper.append(int(appendage))
                for k in range(1,len(quan_nums[i][j:])):
                    secondary=quan_nums[i][j+k]
                    if k == 1:
                        if secondary=='-':
                            continue
                    elif secondary.isdigit()==False:
                        appendage=quan_nums[i][(j+1):(j+k)]
                        k_upper.append(int(appendage))
                        break
                break
    
    return j_upper,k_upper
    
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err#(velint_intK_err.to('K km s-1')/velocityintegrated_intensity_K.to('K km s-1'))*nuppercalc
    return nuppercalc,nuppererr#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def S_j(j_upper,k_upper):#Works for symmetric tops
    return (j_upper**2-k_upper**2)/(j_upper*(2*j_upper+1))
    
def KtoJ(T):
    return (3/2)*k*T
    
def Ntot_rj_thin_nobg(nu,s,g,q,eu_J,T_ex,vint_intensity):
    #nu=nu
    #T_ex=T_ex
    #T_r=T_r
    return ((3*k)/(8*np.pi**3*nu*mu_a**2*s*R_i))*(q/g)*np.exp((eu_J/(k*T_ex)))*((vint_intensity))#((nu+templatewidth)-(nu-templatewidth)))
    
def t_rad(tau_nu, ff, nu, T_ex):
    return ff*(1-np.exp(-tau_nu))*(rjequivtemp(nu, T_ex)-rjequivtemp(nu,Tbg))
    
def nupper_estimated(n_tot,g,q,euj,tex):
    return n_tot*(g/q)*np.exp(-euj/(k*tex))
    
def opticaldepth(aij,nu,T_ex,nupper,lw):
    return (c**2/(8*np.pi*nu**2*lw))*aij*nupper*np.exp((h*nu)/(k*T_ex))
     
qrot_partfunc=Q_rot_asym(testT).to('')

sanitytable2 = utils.minimize_table(Splatalogue.query_lines(200*u.GHz, 300*u.GHz, chemical_name=' HCN ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=['JPL'],
                                    show_upper_degeneracy=True))

incubes=glob.glob(outpath+"*pbcor_line.fits")#'/blue/adamginsburg/d.jeff/imaging_results/field1core1box2/*.fits')

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in incubes:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

maskname="/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/7.0mhzmasked_spw08-7slab.fits"
maskeddatacube=sc.read(maskname)
maskeddatacube=maskeddatacube.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=220027805942.10373*u.Hz)
stdhome='/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/'

cubemaskarray=maskeddatacube.get_mask_array()

sourcelocs={'SgrB2S':'K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/field10originals_z0_000186431_5-6mhzwidth_stdfixes/','DSv':f'/{int(testT.value)}K_field10originals_z0_00186431_5-6mhzwidth_stdfixes_test/'}

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
nupperpath=sourcepath+'nuppers/'
stdpath=sourcepath+'errorimgs/std/'
slabpath=sourcepath+'spectralslabs/km_s/'
mom0path=sourcepath+'mom0/'
rotdiagpath=sourcepath+'pixelwiserotationaldiagrams/'
figpath=sourcepath+'figures/'

picklepath=sourcepath+'ch3ohlinesdict.obj'

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
    

#pdb.set_trace()

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
    stdcellsize=(np.abs(stdimage[0].header['CDELT1']*u.deg)).to('arcsec')
    stdcutoutsize=round(((float(region[43:52])*u.deg)/stdcellsize).to('').value)
    stddata=stdimage[0].data*u.K
    
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
    stdwcs=WCS(stdimage[0].header)#WCS(stdimage[0].header)
    
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
    
    stdcutout=Cutout2D(stddata,(stdpixxcrd,stdpixycrd),(stdcutoutsize))
    assert np.shape(stdcutout)[0]==cube.shape[1], 'Standard deviation cutout size mismatch'
    
    #pdb.set_trace()
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    #print(freq_max)
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Inverted spectral axis'
    print('Passed increasing spectral axis check')
    #print(freq_min)
    linelist='JPL'
    
    print('Peforming Splatalogue queries')
    maintable = utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=[linelist],
                                    show_upper_degeneracy=True))
    '''Needed for upper state degeneracies'''                                
    sparetable=Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=[linelist],
                                    show_upper_degeneracy=True)
                                    
    
    print('Gathering Splatalogue table parameters')    
    lines=maintable['Freq']*10**9*u.Hz/(1+z)#Redshifted to source
    #masterlines.append(lines)
    #vel_lines=vradio(lines,spw1restfreq)
    qns=maintable['QNs']
    euks=maintable['EU_K']*u.K
    eujs=[]
    for eupper_K in euks:
        eujs.append(KtoJ(eupper_K))
    degeneracies=sparetable['Upper State Degeneracy']
    log10aijs=maintable['log10_Aij']
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
    
    #kstdimgpath=stdpath+f'{images[imgnum]}fluxstd.fits'
    kkmsstdimgpath=stdpath+f'{images[imgnum]}intensitystd.fits'
    if os.path.isfile(kkmsstdimgpath):
        print(f'{images[imgnum]} brightness std image already exists')
        spwstdarray=stdcutout.data#fits.getdata(kstdimgpath)*u.K
        kkmsstdarray=fits.getdata(kkmsstdimgpath)*u.K*u.km/u.s
        print(f'Retrieved integrated intensity std data from {kkmsstdimgpath}\n')
    else:
        print(f'Start {images[imgnum]} std calculations')
        spwstdarray=stdcutout.data
        kkmsstdarray=stdcutout.data*linewidth_vel#Stopgap until figure out how to do this on per-transition basis around cores only#pixelwisestd(cube)#spwstdarray,
        #for stdarray, imgpath in zip([spwstdarray,kkmsstdarray],[kstdimgpath,kkmsstdimgpath]):
        print('Set Primary HDU')
        hdu=fits.PrimaryHDU(spwstdarray.value)
        '''This transmoment0 file has intensity (K km/s) units'''
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
        #if np.all(stdarray==spwstdarray):
        #    hdu.header['BUNIT']='K'
        #else:
        hdu.header['BUNIT']='K km s-1'
        print('Wrapping Primary HDU in HDUList')
        hdul=fits.HDUList([hdu])
        print(f'Writing to {kkmsstdimgpath}')
        hdul.writeto(kkmsstdimgpath,overwrite=True)
        print(f'{images[imgnum]} std calculations complete.\n')
    #transitiondict.update({'restfreq':spwrestfreq})
    #,('pixel_0',(pixycrd,pixxcrd))])
    kstddict.update([(images[imgnum],spwstdarray)])
    kkmsstddict.update([(images[imgnum],kkmsstdarray)])
    print(f'Finished loop for {images[imgnum]}')
    #masterqns.append(slicedqns)
    #pdb.set_trace()

pdb.set_trace()

if os.path.isfile(picklepath):
    print(f'pickle {picklepath} already exists.')
else:
    print('Saving dictionary pickle...')
    f=open(picklepath,'wb')
    pickle.dump(spwdict,f)
    f.close()
    print(f'Dictionary pickle saved at {picklepath}')

print('Computing K km/s intensities and K brightness temperatures')
intensityerror=[]
intensities,t_brights=brightnessTandintensities(spwdict)#JybeamtoKkms(spwdict)

#print(t_brights)
print(intensityerror)

'''    
def T_ex(tb,datums):
    ts=[]
    for i in range(len(datums.keys())):
        insert=tb[i]
        nu=datums[i]['freq']
        if tb[i]>0:
            tex=(h*nu)/(k*np.log(((h*nu)/(insert*k))+1))
            ts.append(tex)
        else:
            continue
    
    return ts'''

    
#howmanybeams=2.424e-7/8.57915480931599e-5#Omega=solid_angle(8.57915480931599e-5*u.deg)#BMAJ
#jyhz=fluxes*howmanybeams
#print(jyhz.to('Jy Hz'))
#nohzflux=(jyhz[7]*c/lines[0]).to('Jy km s-1')
#print(nohzflux)
#print(lines)#vint_intensities.to('erg s-1 cm-2 sr-1 Hz-1 km s-1'))
#vint_trads=nohzflux*((c)**2/(2*k*lines[0]**2))
#vint_trads=vint_trads.to('K km s-1')
#print(vint_trads)
#texs=T_ex(t_brights,datadict)
print('Begin fitting procedure\nCompute N_uppers')
spwdictkeys=spwdict.keys()
print(f'spwdictkeys: {spwdictkeys}')
testyshape=np.shape(cube)[1]#60
testxshape=np.shape(cube)[2]#60
testzshape=len(mastereuks)
nugsmap=np.empty(shape=(testyshape,testxshape,testzshape))
nugserrormap=np.empty(shape=(testyshape,testxshape,testzshape))
orderedeuks=[]
ordereddegens=[]
print(f'Begin pixel loops of shape ({testyshape},{testxshape})')
pixelzcoord_nupper=0
pixelzcoord_nuperr=0
for key in spwdictkeys:
    transdict=spwdict[key]
    #print(f'transdict: {transdict}')
    transitionkeys=list(spwdict[key])
    #print(f'transitionkeys: {transitionkeys}')
    for transkey in range(len(transitionkeys)):#Need to figure out way to store the n_us per pixel, per moment map. possibly append in 3D array
        print(f'Transition: {transitionkeys[transkey]}/Nupper array z-coord: {pixelzcoord_nupper}')
        nupperimage_filepath=nupperpath+'CH3OH~'+transitionkeys[transkey]+'.fits'
        nuperrorimage_filepath=nupperpath+'CH3OH~'+transitionkeys[transkey]+'error.fits'
        orderedeuks.append(transdict[transitionkeys[transkey]]['euk'])
        ordereddegens.append(transdict[transitionkeys[transkey]]['degen'])
        
        nupperimgexists=False
        nuperrorimgexists=False
        if os.path.isfile(nupperimage_filepath):
            print(f'{nupperimage_filepath} already exists.\nPopulating nuppers array...\n')
            tempnupper=fits.getdata(nupperimage_filepath)
            nupperimgexists=True
            nugsmap[:,:,pixelzcoord_nupper]=tempnupper
            pixelzcoord_nupper+=1
        if os.path.isfile(nuperrorimage_filepath):
            print(f'{nuperrorimage_filepath} already exists\nPopulating nupper error array...\n')
            tempnuerr=fits.getdata(nuperrorimage_filepath)
            nuperrorimgexists=True
            nugserrormap[:,:,pixelzcoord_nuperr]=tempnuerr
            pixelzcoord_nuperr+=1
        elif not nupperimgexists or not nuperrorimgexists:
            for y in range(testyshape):
                print(f'Row {y} Looping')
                n_us=[]#np.empty(np.shape(intensityerror))
                n_uerr=[]#np.empty(np.shape(intensityerror))
                for x in range(testxshape):
                    if nupperimgexists:
                        n_us.append((tempnupper[y,x])/transdict[transitionkeys[transkey]]['degen'])
                    if nuperrorimgexists:
                        n_uerr.append((tempnuerr[y,x])/transdict[transitionkeys[transkey]]['degen'])
                    else:
                        tempnupper,tempnuerr=N_u(transdict[transitionkeys[transkey]]['freq'],transdict[transitionkeys[transkey]]['aij'],intensities[transitionkeys[transkey]][y,x],kkmsstddict[key][y,x])
                        n_us.append((tempnupper.to('cm-2')*u.cm**2)/transdict[transitionkeys[transkey]]['degen'])
                        n_uerr.append((tempnuerr.to('cm-2')*u.cm**2)/transdict[transitionkeys[transkey]]['degen'])  
                nugsmap[y,:,pixelzcoord_nupper]=n_us
                nugserrormap[y,:,pixelzcoord_nuperr]=n_uerr
            if not nupperimgexists:
                nupperimgdata=nugsmap[:,:,pixelzcoord_nupper]
                primaryhdu=fits.PrimaryHDU(nupperimgdata)
                transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
                transmom0header=transmoment0[0].header
                primaryhdu.header=transmom0header
                primaryhdu.header['BTYPE']='Upper-state column density'
                primaryhdu.header['BUNIT']='cm-2'
                hdul=fits.HDUList([primaryhdu])
                hdul.writeto(nupperimage_filepath,overwrite=True)
            if not nuperrorimgexists:
                nuperrorimgdata=nugserrormap[:,:,pixelzcoord_nuperr]
                primaryhduerr=fits.PrimaryHDU(nuperrorimgdata)
                transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
                transmom0header=transmoment0[0].header
                primaryhduerr.header=transmom0header
                primaryhduerr.header['BTYPE']='Upper-state column density'
                primaryhduerr.header['BUNIT']='cm-2'
                hdulerr=fits.HDUList([primaryhduerr])
                hdulerr.writeto(nuperrorimage_filepath)
            pixelzcoord_nupper+=1
            pixelzcoord_nuperr+=1
            #pdb.set_trace()
        #print(n_us)
        #print(n_uerr)
      
print('pixels looped, Nupper calcs complete\n')  

'''
for key in spwdictkeys:
    transitionkeys=list(spwdict[key])
    for transkey in range(len(transitionkeys)):
        nupperimage_filepath=filepath2+'CH3OH~'+transitionkeys[transkey]+'.fits'
        nuperrorimage_filepath=filepath2+'CH3OH~'+transitionkeys[transkey]+'error.fits'
        if os.path.isfile(nupperimage_filepath):
            print(f'{nupperimage_filepath} already exists')
        elif os.path.isfile(nupperimage_filepath)==False:
            nupperimgdata=nugsmap[:,:,transkey]
            primaryhdu=fits.PrimaryHDU(nupperimgdata)
            primaryhdu.header['UNIT']='cm-2'
            hdul.writeto(nupperimage_filepath,overwrite=True)
            hdul=fits.HDUList([primaryhdu])
        elif os.path.isfile(nuperrorimage_filepath):
            print(f'{nuperrorimage_filepath} already exists')
        elif os.path.isfile(nuperrorimage_filepath)==False:
            primaryhduerr=fits.PrimaryHDU(nuperrorimgdata)
            primaryhduerr.header['UNIT']='cm-2'
            hdulerr=fits.HDUList([primaryhduerr])
            hdulerr.writeto(nuperrorimage_filepath)
'''

print('Setting up and executing model fit')
texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))
texsigclipmap=np.empty((testyshape,testxshape))
texsnrmap=np.empty((testyshape,testxshape))
numtransmap=np.empty((testyshape,testxshape))
degensforfit=[]
snr=3

fitdict={}
#pdb.set_trace()
for y in range(testyshape):
    print(f'Start Row {y} fitting')
    for x in range(testxshape):
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        for zed in range(testzshape):
            if nugsmap[y,x,zed] <= 0 or np.isnan(nugsmap[y,x,zed]):
                continue
            else:
                nupperstofit.append(nugsmap[y,x,zed])
                eukstofit.append(mastereuks[zed])
                nuperrors.append(nugserrormap[y,x,zed])
                degensforfit.append(ordereddegens[zed])
        #pdb.set_trace()
        numtransmap[y,x]=len(nupperstofit)        
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texsnrmap[y,x]=np.nan
            texsigclipmap[y,x]=obsTex
            texerrormap[y,x]=np.nan
        else:
            #log10nuerr=[]
            errstofit=[]
            
            for num in range(len(nupperstofit)):
                templog10=(1/nupperstofit[num])*nuperrors[num]
                temperrfit=1/templog10
                #log10nuerr.append(templog10)
                errstofit.append(temperrfit)
            
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit),weights=errstofit)
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTex=-np.log10(np.e)/(fit_lin.slope)
            obsNtot=qrot_partfunc*10**(np.log10(nupperstofit[0])+fit_lin.slope*eukstofit[0])
            dobsTex=(eukstofit[0]*u.K*np.log(10)*np.log(np.e))/(np.log(nupperstofit[0]/degensforfit[0])-np.log(obsNtot/qrot_partfunc))**2
            sigTex=(obsTex*u.K/dobsTex).to('')
            
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTex.to('K').value
            texsnrmap[y,x]=sigTex
            if sigTex >= snr:
                texsigclipmap[y,x]=obsTex
            else:
                texsigclipmap[y,x]=np.nan

detectnum=5
transmaskarr=np.ma.masked_where(numtransmap<detectnum,texsigclipmap)
transmasktexmap=transmaskarr.filled(fill_value=np.nan)#np.array(np.ma.masked_where(numtransmap<detectnum,texsigclipmap))
            
transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
transmom0header=transmoment0[0].header

primaryhdutex=fits.PrimaryHDU(texmap)
primaryhdutex.header=transmom0header
primaryhdutex.header['BTYPE']='Excitation temperature'
primaryhdutex.header['BUNIT']='K'
hdultex=fits.HDUList([primaryhdutex])
print(f'Saving raw temperature map at {sourcepath+"texmap_allspw_withnans_weighted.fits"}\n')
hdultex.writeto(sourcepath+'texmap_allspw_withnans_weighted.fits',overwrite=True)

primaryhduntot=fits.PrimaryHDU(ntotmap)
primaryhduntot.header=transmom0header
primaryhduntot.header['BTYPE']='Total column density'
primaryhduntot.header['BUNIT']='cm-2'
hdulntot=fits.HDUList([primaryhduntot])
print(f'Saving raw ntot map at {sourcepath+"ntotmap_allspw_withnans_weighted.fits"}\n')
hdulntot.writeto(sourcepath+'ntotmap_allspw_withnans_weighted.fits',overwrite=True)

primaryhdutexerr=fits.PrimaryHDU(texerrormap)
primaryhdutexerr.header=transmom0header
primaryhdutexerr.header['BTYPE']='Excitation temperature'
primaryhdutexerr.header['BUNIT']='K'
hdultexerror=fits.HDUList([primaryhdutexerr])
print(f'Saving temperature error map at {sourcepath+"texmap_error_allspw_withnans_weighted.fits"}\n')
hdultexerror.writeto(sourcepath+'texmap_error_allspw_withnans_weighted.fits',overwrite=True)

primaryhdutexclip=fits.PrimaryHDU(texsigclipmap)
primaryhdutexclip.header=transmom0header
primaryhdutexclip.header['BTYPE']='Excitation temperature'
primaryhdutexclip.header['BUNIT']='K'
hdultexclip=fits.HDUList([primaryhdutexclip])
nsigmatexpath=sourcepath+f'texmap_{snr}sigma_allspw_withnans_weighted.fits'
print(f'Saving {snr}sigma temperature map at {nsigmatexpath}')
hdultexclip.writeto(nsigmatexpath,overwrite=True)

primaryhdutexsnr=fits.PrimaryHDU(texsnrmap)
primaryhdutexsnr.header=transmom0header
primaryhdutexsnr.header['BTYPE']='Excitation temperature SNR'
primaryhdutexsnr.header['BUNIT']=''
hdultexsnr=fits.HDUList([primaryhdutexsnr])
print(f'Saving {snr}sigma temperature map at {sourcepath+"texmap_{snr}sigma_allspw_withnans_weighted.fits"}\n')
hdultexsnr.writeto(sourcepath+'texmap_snr_allspw_weighted.fits',overwrite=True)

primaryhdunumtrans=fits.PrimaryHDU(numtransmap)
primaryhdunumtrans.header=transmom0header
primaryhdunumtrans.header['BTYPE']='Number CH3OH 3sigma Detected Transitions'
primaryhdunumtrans.header['BUNIT']=''
hdulnumtrans=fits.HDUList([primaryhdunumtrans])
print(f'Saving number of {snr}sigma detected CH3OH lines map at {sourcepath+"ch3ohdetections_{snr}sigma_allspw_withnans_weighted.fits"}\n')
hdulnumtrans.writeto(sourcepath+f"ch3ohdetections{detectnum}_{snr}sigma_allspw_withnans_weighted.fits",overwrite=True)

primaryhdutransmasktex=fits.PrimaryHDU(transmasktexmap)
primaryhdutransmasktex.header=transmom0header
primaryhdutransmasktex.header['BTYPE']='Excitation temperature'
primaryhdutransmasktex.header['BUNIT']='K'
hdultransmasktex=fits.HDUList([primaryhdutransmasktex])
nsigmatransmaskedpath=sourcepath+f"texmap_{detectnum}transmask_{snr}sigma_allspw_withnans_weighted.fits"
print(f'Saving {detectnum} transition masked temperature map at {nsigmatransmaskedpath}\n')
hdultransmasktex.writeto(nsigmatransmaskedpath,overwrite=True)

nugs_swapaxis2toaxis0=np.swapaxes(nugsmap,0,2)
nugserr_swapaxis2toaxis0=np.swapaxes(nugserrormap,0,2)

nugs_swapxwithy=np.swapaxes(nugs_swapaxis2toaxis0,1,2)
nugserr_swapxwithy=np.swapaxes(nugserr_swapaxis2toaxis0,1,2)

nugscube=fits.PrimaryHDU(nugs_swapxwithy)
nugserrcube=fits.PrimaryHDU(nugserr_swapxwithy)

nugscube.header=transmom0header
nugserrcube.header=transmom0header

nugscube.header['BUNIT']='cm-2'
nugserrcube.header['BUNIT']='cm-2'

nugserrcube.header['BTYPE']='Upper-state column density error'
nugscube.header['BTYPE']='Upper-state column density'

nugshdul=fits.HDUList([nugscube])
nugserrhdul=fits.HDUList([nugserrcube])

print('Saving alltransition and E_U(K) lists\n')
nugshdul.writeto(sourcepath+'alltransitions_nuppers.fits',overwrite=True)
nugserrhdul.writeto(sourcepath+'alltransitions_nupper_error.fits',overwrite=True)
eukqns=np.column_stack((mastereuks,masterqns,masterlines,ordereddegens))
np.savetxt(sourcepath+'mastereuksqnsfreqsdegens.txt',eukqns,fmt='%s',header=f'Methanol transitions, excitation temperatures, and degeneracies used in this folder. Temperatures in units of K, frequencies are redshifted ({z}/{(z*c).to("km s-1")}) and in Hz.')

'''This wing of the code plots up the temperature and ntot maps'''

print('Begin plotting procedure.\n')

def make_scalebar(ax, left_side, length, color='black', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec, coordsys='icrs'):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform(coordsys),
                   zorder=3)
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform(coordsys),
                  color=color,
                  fontsize=fontsize,
                 zorder=2,bbox=dict(facecolor='white', alpha=0.6))
    ax.axis(axlims)
    return lines,txt

colormap= copy.copy(mpl.cm.get_cmap("inferno"))
colormap.set_bad('black')
dGC=8.34*u.kpc#per Meng et al. 2019 https://www.aanda.org/articles/aa/pdf/2019/10/aa35920-19.pdf

plottexhdu=fits.open(nsigmatransmaskedpath)[0]

plottexwcs=WCS(plottexhdu)

sliced=['x','y']
ax=plt.subplot(projection=plottexwcs,slices=sliced)
plt.rcParams['figure.dpi'] = 150

ra=ax.coords[0]
dec=ax.coords[1]

plottedtex=ax.imshow(plottexhdu.data,cmap=colormap,vmax=1000,vmin=10)

scale=5000*u.AU
lenn=np.arctan(scale/dGC)

if source == 'DSi':
    print(f'Scalebar source : DSi')
    make_scalebar(ax, coordinates.SkyCoord('17:47:19.5180 -28:23:51.359', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:19.4976 -28:23:51.384', unit=(u.hour,u.deg), frame='icrs')

elif source == 'SgrB2S':
    print(f'Scalebar source: SgrB2S')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.6618 -28:23:48.734', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.6590 -28:23:48.772', unit=(u.hour,u.deg), frame='icrs')

ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(plottedtex,pad=0,label='K')#plt.colorbar()
#plt.savefig(saveimgpath)
plt.show()

print('Cube-->core-->Texmap complete.')
