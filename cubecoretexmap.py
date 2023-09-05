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
import sys
import reproject
from utilities import *
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from pyspeckit.spectrum.models import lte_molecule

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

'''This wing of the script takes in continuum-subtracted cubes, cuts out a subcube around a region of interest based on a DS9 region, and converts the subcubes into brightness temperature (K) units'''
print('Cube-->Core-->Tex start\n')
print('Begin Jy/beam-to-K and region subcube conversion\n')

source='DSiii'
print(f'Source: {source}\n')
fields={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7,'DS10':1,'DS11':1,'DSXI':8}
fnum=fields[source]

#inpath="/orange/adamginsburg/sgrb2/d.jeff/data/field10originalimages/"
inpaths={1:'/orange/adamginsburg/sgrb2/d.jeff/data/OctReimage_K/',10:"/orange/adamginsburg/sgrb2/d.jeff/data/field10originals_K/",2:"/orange/adamginsburg/sgrb2/d.jeff/data/field2originals_K/",3:"/orange/adamginsburg/sgrb2/d.jeff/data/field3originals_K/",7:"/orange/adamginsburg/sgrb2/d.jeff/data/field7originals_K/",8:"/orange/adamginsburg/sgrb2/d.jeff/data/field8originals_K/"}
inpath=inpaths[fnum]#'/blue/adamginsburg/d.jeff/imaging_results/data/OctReimage/'
beamcubes=glob.glob(inpath+'*.fits')
homes={1:'/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/',10:"/orange/adamginsburg/sgrb2/d.jeff/products/field10originals_K/",2:"/orange/adamginsburg/sgrb2/d.jeff/products/field2originals_K/",3:"/orange/adamginsburg/sgrb2/d.jeff/products/field3originals_K/",7:"/orange/adamginsburg/sgrb2/d.jeff/products/field7originals_K/",8:"/orange/adamginsburg/sgrb2/d.jeff/products/field8originals_K/"}
home=homes[fnum]#'/blue/adamginsburg/d.jeff/imaging_results/products/OctReimage/'
cubes=glob.glob(home+'*pbcor_line.fits')
sourceregs={'SgrB2S':'fk5; box(266.8353410, -28.3962005, 0.0016806, 0.0016806)','DSi':'fk5; box(266.8316387, -28.3971867, 0.0010556, 0.0010556)','DSii':'fk5; box(266.8335363, -28.3963159, 0.0006389, 0.0006389)','DSiii':'fk5; box(266.8332758, -28.3969270, 0.0006944, 0.0006944)','DSiv':'fk5; box(266.8323834, -28.3954424, 0.0009000, 0.0009000)','DSv':'fk5; box(266.8321331, -28.3976585, 0.0005556, 0.0005556)','DSVI':'fk5; box(266.8380037, -28.4050741, 0.0017361, 0.0017361)','DSVII':'fk5; box(266.8426074, -28.4094401, 0.0020833, 0.0020833)', 'DSVIII':'fk5; box(266.8418408, -28.4118242, 0.0014028, 0.0014028)','DSIX':'fk5; box(266.8477371, -28.4311386, 0.0009583, 0.0009583)','DS10':'fk5; box(266.8373798, -28.4009340, 0.001, 0.001)','DS11':'fk5; box(266.8374572, -28.3996894, 0.0009, 0.0009)','DSXI':'fk5; box(266.8404733, -28.4286378, 0.0013194, 0.0013194)'}
#region='fk5; box(266.8321311,-28.3976633, 0.0010833, 0.0010833)'#DSv
#region='fk5; box(266.8323834,-28.39544244, 0.0009000, 0.0009000)'#DSiv
region=sourceregs[source]#'fk5; box(266.8316387, -28.3971867, 0.0010556, 0.0010556)'#DSi-large
#region='fk5; box(266.8353410, -28.3962005, 0.0016806, 0.0016806)'#SgrB2S-box2
#region='fk5; box(266.8350804, -28.3963256, 0.0023889, 0.0023889)'#SgrB2S-large
#box(266.8335363, -28.3963159, 0.0006389, 0.0006389)' #DSii
#/iii
#box(266.8315833, -28.3971867, 0.0006528, 0.0006528)' #DSi-small
outpath_base=f'/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/{source}/'
outstatpath_end={1:'OctReimage_K/',10:'field10originals_K/',2:'field2originals_K/',3:'field3originals_K/',7:'field7originals_K/',8:'field8originals_K/'}
outpath=outpath_base+outstatpath_end[fnum]#f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/field10originals/'
#outpath=f'/blue/adamginsburg/d.jeff/SgrB2DSminicubes/{source}/OctReimage_K/'#imaging_results/DSii_iiibox1/'
statfixpath_base='/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSstatcontfix/'
statfixpath=statfixpath_base+outstatpath_end[fnum]#f'/blue/adamginsburg/d.jeff/SgrB2DSstatcontfix/OctReimage_K/'

regionparams=[float(val) for val in region[9:(len(region)-1)].split(', ')]

if os.path.isdir(outpath):
    print(f'Minicube directory "{outpath}" already exists. Proceeding to line loops/LTE fitting procedure.\n')
    sourceisnew=False
    pass
else:
    sourceisnew=True
  
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
            
    #pdb.set_trace()
            
    for sppw, cub in zip(images, cubestobox):
        boxcubename=outpath+sppw+'minimize.image.pbcor_line.fits'
        fullsizecube=sc.read(cub,use_dask=True)
        if os.path.isfile(boxcubename):
            print(f'{boxcubename} already exists. Skipping...\n')
            continue
        elif fullsizecube.unit=='K':
            print(f'Grabbing datacube from {cub}')
            fullsizecube.allow_huge_operations=True
            spwrestfreq=fullsizecube.header['RESTFRQ']*u.Hz
            print('Boxed cube already has brightness (K) units.')
            print(f'Creating subcube containing region around core {source}')
            boxedsubcubeK=fullsizecube.subcube_from_ds9region(region)
            print(f'Saving to {boxcubename}')
            boxedsubcubeK.write(boxcubename,format='fits',overwrite=True)
            print(f'Done. Continuing to next spw')
            continue
        else:
            print(f'Grabbing datacube from {cub}')
            fullsizecube.allow_huge_operations=True
            spwrestfreq=fullsizecube.header['RESTFRQ']*u.Hz
            print('Creating subcube and converting from Jy/beam to K')
            boxedsubcubeK=fullsizecube.subcube_from_ds9region(region).to(u.K)
            #print('Converting spectral axis to km/s')
            #boxedsubcubeKkms=boxedsubcubeK.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
            print(f'Saving to {boxcubename}')
            boxedsubcubeK.write(boxcubename,format='fits',overwrite=True)
            print('Finished\n')
        
if sourceisnew:
    print(f'Since source {source} is new, pausing script for staging.')
    print('Enter "exit" to stage the core.\n')
    sys.exit()#pdb.set_trace() 
else:
    pass      
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

Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=150*u.GHz,
                                                         fmax=300*u.GHz)

dopplershifts={'SgrB2S':0.000234806,'DSi':0.0001842772437139578,'DSii':0.00016236367659115043,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001661546432045067,'DSIX':0.00015453732389175085,'DS10':0.00015794099431186572,'DS11':0.00016667684485037162,'DSX':0.00016375916278648755}#:0.000190713}/old doppler S: 0.0002306756533745274/old doppler I: 0.000186431

z=dopplershifts[source]
#z=0.00017594380066803095 #SgrB2DSII?
#z=0.000186431 #SgrB2DSi/DSiv(?)
#z=0.0002306756533745274#<<average of 2 components of 5_2-4_1 transition using old redshift(0.000236254)#0.000234806#0.0002333587 SgrB2S
print(f'Doppler shift: {z} / {(z*c).to("km s-1")}\n')

print('Setting input LTE parameters')
trotdict={'SgrB2S':300*u.K,'DSi':400*u.K,'DSii':200*u.K,'DSiii':300*u.K,'DSiv':150*u.K,'DSv':150*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':175*u.K,'DSIX':150*u.K,'DS10':150*u.K,'DS11':100*u.K,'DSX':100*u.K}#old dsv - 1e16, 150K
testT=trotdict[source]#500*u.K
ntotdict={'SgrB2S':1e18*u.cm**-2,'DSi':1e18*u.cm**-2,'DSii':1e17*u.cm**-2,'DSiii':1e17*u.cm**-2,'DSiv':1e18*u.cm**-2,'DSv':1e17*u.cm**-2,'DSVI':1e17*u.cm**-2,'DSVII':1e16*u.cm**-2,'DSVIII':1e16*u.cm**-2,'DSIX':1e16*u.cm**-2,'DS10':1e16*u.cm**-2,'DS11':1e16*u.cm**-2,'DSX':1e15*u.cm**-2}
testntot=ntotdict[source]
print(f'Input Tex: {testT}\nInput Ntot: {testntot}')
    
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

'''    
def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix-1),(xpix-1)])#Corrects for different pixel counting procedures
    return vals
'''
 
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
                temp=velflux_T/linewidth_vel#***What's going on here? Should this be fixed because it uses one uniform linewidth
                t_bright.update({temptransdictkeys[i]:temp})
                d_velfluxT=(temptransdict[temptransdictkeys[i]]['stddev'])#/temp)*velflux_T
                intensityerror.append(d_velfluxT)
                #pdb.set_trace()
                
    return intensitydict,t_bright
'''
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):#taken from pyspeckit documentation https://pyspeckit.readthedocs.io/en/latest/lte_molecule_model.html?highlight=Aij#lte-molecule-model
    nuppercalc=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K
    nuppererr=((8*np.pi*k*nu**2)/(h*c**3*Aij))*velint_intK_err#(velint_intK_err.to('K km s-1')/velocityintegrated_intensity_K.to('K km s-1'))*nuppercalc
    return nuppercalc,nuppererr#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
'''
def S_j(j_upper,k_upper):#Works for symmetric tops
    return (j_upper**2-k_upper**2)/(j_upper*(2*j_upper+1))
    
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
        line_width_freq=velocitytofreq(line_width,line)
        #pdb.set_trace()
        nu_upper=line+line_width_freq
        nu_lower=line-line_width_freq
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        #slabstart=time.time()
        slab=cube.spectral_slab(nu_upper,nu_lower)
        oldstyleslab=cube.spectral_slab((nu_upper-nu_offset),(nu_lower+nu_offset))
        peakchannel=slab.closest_spectral_channel(line)
        print(f'Peak channel: {peakchannel}')
        slabbeams=(slab.beams.value)*u.sr/u.beam
        slab_K=slab[:,pixycrd,pixxcrd]
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
        intertau=lte_molecule.line_tau(testT,testntot,qrot_partfunc,degeneracies[i],restline,eujs[i],aijs[i])
        phi_nu=lineprofile(linewidth_sigma,restline,restline)
        est_tau=(intertau*phi_nu).to('')
        trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=testT).to('K')
        
        print('LTE params calculated')
        print(f'tbthick: {tbthick}\n targetspecK_stddev: {targetspecK_stddev}\n peak_amplitude: {peak_amplitude}')
        print(f'est_nupper: {est_nupper}\n est_tau: {est_tau}\n trad: {trad}')
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        moment0filename=home+'CH3OH~'+transition+'_raw.fits'
        maskedmom0fn=home+'CH3OH~'+transition+'_masked.fits'
        maskresidualfn=home+'CH3OH~'+transition+'_residual.fits'
        maskedmom0errfn=home+'CH3OH~'+transition+'_error.fits'
        slabfilename=slabpath+'CH3OH~'+transition+'_slab.fits'
        maskedslabfn=slabpath+'CH3OH~'+transition+'_maskedslab.fits'
        maskfn=slabpath+'CH3OH~'+transition+'_mask.fits'
        peakintfn=home+'CH3OH~'+transition+'_peakint.fits'
        if os.path.isfile(maskedmom0fn):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(maskedmom0fn)*u.K*u.km/u.s
            #isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(maskedmom0fn)
            isfilestdflux=stddata#fits.getdata(f'{stdpath}{images[imgnum]}fluxstd.fits')*u.K#This is confusing, notation-wise, but I'm leaving it this way for now since it's consistent between the two forks in the loop. For future reference: isfilestdflux is the error on the measured brightnesstemp in K, whereas isfilemom0 pulls from the moment0 maps and is in K km/s
            isfilemom0err=fits.getdata(maskedmom0errfn)*u.K*u.km/u.s
            temptransdict.update([('freq',restline),('flux',isfilemom0),('stddev',isfilestdflux),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line),('mom0err',isfilemom0err)])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print('\nDictionaries populated for this transition.')
            if os.path.isfile(maskedslabfn):
                print('Masked slab already exists...\n')
                pass
            else:
                slab.write(maskedslabfn)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            if os.path.isfile(peakintfn):
                print('Peak intensity file already exists.')
                pass
            else:
                print('Peak intensity procedure starting')
                print('Creating masked slab')
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                slab3sigmamask=slab > (3*stdcutout.data)
                slab=slab.with_mask(slab3sigmamask)
                slabspecax=slab.spectral_axis
                slabmom1=slab.moment1()
                slabfwhm=slab.linewidth_fwhm()#(7*u.MHz/line)*c.to('km s-1')#
                cubemask=(slabspecax[:,None,None] < (velocityfield_representative + fwhm_representative)[None,:,:]) & (slabspecax[:,None,None] > (velocityfield_representative - fwhm_representative)[None,:,:])
                maskedslab=slab.with_mask(cubemask)
                print('Computing peak intensity image')
                maskedpeakint=maskedslab.max(axis=0)
                print(f'Saving to {peakintfn}')
                maskedpeakint.write(peakintfn)
                print('Peak intensity image saved.\n')
            for moment in [1,2]:
                slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                momentnfilenames=[(sourcepath+'mom1/CH3OH~'+transition+'.fits'),(sourcepath+'mom2/CH3OH~'+transition+'_var.fits')]
                fwhmfilename=sourcepath+'mom2/CH3OH~'+transition+'_fwhm.fits'
                if os.path.isfile(momentnfilenames[moment-1]):
                    print(f'{transition} moment{moment} file already exists.\n')
                    continue
                elif moment == 1:
                    print(f'Computing moment 1 and saving to {momentnfilenames[moment-1]}\n')
                    slabmom1=slab.moment1()
                    slabmom1.write(momentnfilenames[moment-1])
                elif moment == 2:
                    print(f'Computing moment 2 and saving to {momentnfilenames[moment-1]}\n')
                    slabmom2=slab.moment2()
                    slabfwhm=slab.linewidth_fwhm()
                    slabmom2.write(momentnfilenames[moment-1])
                    slabfwhm.write(fwhmfilename)
            pass
        elif trad >= 3*targetspecK_stddev and peak_amplitude >= 3* targetspecK_stddev:#*u.K:
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])#spwrestfreq)
            if transition in excludedlines[source]:
                print(f'\nExcluded line detected: {quantum_numbers[i]}, E_U: {euks[i]}, Freq: {line.to("GHz")}')
                #sigma1mom0=np.empty((cube.shape[1],cube.shape[2]))
                sigma1mom0=stdcutout.data*linewidth_vel
                sigma1hdu=fits.PrimaryHDU(sigma1mom0.value)
                sigma1hdu.header=fits.open(mom0path+'CH3OH~5_1-4_2E1vt0_raw.fits')[0].header
                sigma1hdu.header['RESTFRQ']=line.value
                sigma1hdul=fits.HDUList([sigma1hdu])
                print(f'1sigma moment0 of shape ({cube.shape[1]},{cube.shape[2]}) populated\nTarget pix flux value: {stddata[pixycrd,pixxcrd]*linewidth_vel}')
                moment0beam=slab.beams
                maskslabmom0=sigma1mom0
                #maskslabmom0.write((maskedmom0fn))
                #pdb.set_trace()
                sigma1hdul.writeto(moment0filename,overwrite=True)
                print(f'Saved to {moment0filename}')
                kkmsstdarray=maskslabmom0
            else:
                print('Commence moment0 procedure\n')
                #cubemask=BooleanArrayMask(mask=cubemaskarray,wcs=slab.wcs)
                print(f'Create {quantum_numbers[i]} spatial-velocity mask')
                slab3sigmamask=slab > (3*stdcutout.data)
                slab=slab.with_mask(slab3sigmamask)
                slabspecax=slab.spectral_axis
                slabmom1=slab.moment1()
                slabfwhm=slab.linewidth_fwhm()#(7*u.MHz/line)*c.to('km s-1')#
                cubemask=(slabspecax[:,None,None] < (velocityfield_representative + fwhm_representative)[None,:,:]) & (slabspecax[:,None,None] > (velocityfield_representative - fwhm_representative)[None,:,:])
                oldstyleslab=oldstyleslab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])
                #if imgnum > 0:
                #    pdb.set_trace()
                print('Masking spectral slab')
                maskedslab=slab.with_mask(cubemask)
                #momstart=time.time()
                print('Unmasked moment0 computing...\n')
                slabmom0=oldstyleslab.moment0()
                print('Masked moment0 computing...\n')

                contmom0=reprojcont_K*slabfwhm#continuum sanity check

                if transition in doublet:
                    maskslabmom0=(maskedslab.moment0()+contmom0)/2
                    print(f'\nDoublet line identified: {quantum_numbers[i]}')
                    print(f'Value divided in half to compensate for line blending.\n')
                else:
                    maskslabmom0=maskedslab.moment0()+contmom0
                #momend=time.time()-momstart
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
                kkmsstdarray=stdcutout.data*fwhm_representative
                if os.path.isfile(slabfilename):
                    print(f'Spectral slab {slabfilename} already exists.\nProceeding...\n')
                    pass
                else:
                    slab.write(slabfilename)
                    print(f'Slab written to {slabfilename}.')
                    maskedslab.write(maskedslabfn)
                    print(f'Masked slab written to {maskedslabfn}.\n')
                    print('Creating moment0 error Primary HDU')
                    kkmshdu=fits.PrimaryHDU(kkmsstdarray.value)
                    kkmshdu.header=maskedslab.header
                    kkmshdu.header['BUNIT']='K km s-1'
                    print('Wrapping moment0 error Primary HDU in HDUList')
                    kkmshdul=fits.HDUList([kkmshdu])
                    print(f'Writing moment0 error to {maskedmom0errfn}')
                    kkmshdul.writeto(maskedmom0errfn,overwrite=True)
                if os.path.isfile(peakintfn):
                    print('Peak intensity file already exists.')
                    pass
                else:
                    print('Peak intensity procedure starting')
                    print('Computing peak intensity image')
                    maskedpeakint=maskedslab.max(axis=0)
                    print(f'Saving to {peakintfn}')
                    maskedpeakint.write(peakintfn)
                    print('Peak intensity image saved.\n')
                for moment in [1,2]:
                    moment1filename=sourcepath+'mom1/CH3OH~'+transition+'.fits'
                    moment2filename=sourcepath+'mom2/CH3OH~'+transition+'_var.fits'
                    fwhmfilename=sourcepath+'mom2/CH3OH~'+transition+'_fwhm.fits'
                    if moment == 1:
                        print(f'Computing moment 1 and saving to {moment1filename}\n')
                        #slabmom1=slab.moment1()
                        slabmom1.write(moment1filename)
                    elif moment == 2:
                        print(f'Computing moment 2s and saving to {moment2filename} and {fwhmfilename}\n')
                        slabmom2=slab.moment2()
                        slabfwhm=slab.linewidth_fwhm()
                        slabmom2.write(moment2filename)
                        slabfwhm.write(fwhmfilename)
                pass
            temptransdict.update([('freq',restline),('flux',maskslabmom0),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename),('shift_freq',line),('mom0err',kkmsstdarray)])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            masterqns.append(quantum_numbers[i])
            masterlines.append(line_list[i].value)
            print(f'{quantum_numbers[i]} calculations complete.\n')
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
     
qrot_partfunc=qrot(testT)#Q_rot_asym(testT).to('')
initialguess_partfun=np.copy(qrot_partfunc)

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

#maskname="/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/7.0mhzmasked_spw08-7slab.fits"
#maskeddatacube=sc.read(maskname)
#maskeddatacube=maskeddatacube.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=220027805942.10373*u.Hz)

stdhomedict={1:'/orange/adamginsburg/sgrb2/d.jeff/products/OctReimage_K/',10:'/orange/adamginsburg/sgrb2/d.jeff/products/field10originals_K/',2:'/orange/adamginsburg/sgrb2/d.jeff/products/field2originals_K/',3:'/orange/adamginsburg/sgrb2/d.jeff/products/field3originals_K/',7:'/orange/adamginsburg/sgrb2/d.jeff/products/field7originals_K/'}
stdhome=stdhomedict[fnum]

#cubemaskarray=maskeddatacube.get_mask_array()

sourcelocs={'SgrB2S':'/lateaug2023corrected/','DSi':'/lateaug2023corrected/','DSii':'/lateaug2023corrected_real_removecontaminants/','DSiii':'/sep2023phi_nu&doublet/','DSiv':'/lateaug2023corrected/','DSv':f'/lateaug2023corrected/','DSVI':'/aug2023qrotfix/','DSVII':f'/aug2023qrotfix/','DSVIII':f'/aug2023qrotfix/','DSIX':f'/aug2023qrotfix/','DS10':'/march2023discovery_5kmslw/','DS11':f'/march2023discovery_5kmslw_{int(testT.value)}K/','DSX':f'/Kfield7originals_{int(testT.value)}K_trial1_noexclusions/'}#'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/'

origsourcelocs={'SgrB2S':'/new_testingstdfixandontheflyrepstuff_K_OctReimage_restfreqfix_newvelmask_newpeakamp/','DSi':'/Kfield10originals_trial7_field10errors_newexclusion_matchslabwidthtorep/','DSii':'/Kfield10originals_noexclusions/','DSiii':'/Kfield10originals_noexclusions/','DSiv':'/Kfield10originals_noexclusions/','DSv':f'/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':'/Kfield2originals_trial3_8_6-8_7excluded/','DSVII':f'/Kfield3originals_{int(testT.value)}K_trial1_noexclusions/','DSVIII':f'/Kfield3originals_{int(testT.value)}K_trial1_noexclusions/','DSIX':f'/Kfield7originals_{int(testT.value)}K_trial1_noexclusions/','DSX':f'/Kfield7originals_{int(testT.value)}K_trial1_noexclusions/'}#

representativelines={'SgrB2S':'20_1-20_0vt=0','DSi':'8_1-7_0vt=0','DSii':'8_1-7_0vt=0','DSiii':'10_2--9_3-vt=0','DSiv':'20_1-20_0vt=0','DSv':'8_1-7_0vt=0','DSVI':'8_1-7_0vt=0','DSVII':'8_1-7_0vt=0','DSVIII':'8_1-7_0vt=0','DSIX':'8_1-7_0vt=0','DS10':'10_2--9_3-vt=0','DS11':'8_1-7_0vt=0','DSX':'8_1-7_0vt=0'}#oldS 4_2-3_1vt=0
representativelws={'SgrB2S':(10*u.km/u.s),'DSi':(3*u.km/u.s),'DSii':(3*u.km/u.s),'DSiii':(3*u.km/u.s),'DSiv':(4*u.km/u.s),'DSv':(4*u.km/u.s),'DSVI':(3*u.km/u.s),'DSVII':(2.5*u.km/u.s),'DSVIII':(2.5*u.km/u.s),'DSIX':(5*u.km/u.s),'DS10':(5*u.km/u.s),'DS11':(5*u.km/u.s),'DSX':(4*u.km/u.s)}#{'SgrB2S':8*u.MHz,'DSi':3.6*u.MHz}#11MHz for ~10 km/s
representativecubes={'SgrB2S':0,'DSi':1,'DSii':1,'DSiii':2,'DSiv':0,'DSv':1,'DSVI':1,'DSVII':1,'DSVIII':1,'DSIX':1,'DS10':2,'DS11':1,'DSX':1}#spwnumber

sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcelocs[source]
nupperpath=sourcepath+'nuppers/'
stdpath=sourcepath+'errorimgs/std/'
slabpath=sourcepath+'spectralslabs/km_s/'
mom0path=sourcepath+'mom0/'
rotdiagpath=sourcepath+'pixelwiserotationaldiagrams/'
figpath=sourcepath+'figures/'

overleafpath="/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/"

picklepath=sourcepath+'ch3ohlinesdict.obj'

if source == 'DS10' or source == 'DS11':
    contpath=sourcepath+'reprojectedcontinuum.fits'
else:
    origsourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+origsourcelocs[source]
    contpath=origsourcepath+'reprojectedcontinuum.fits'

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
    
if source not in origsourcelocs.keys():
    if os.path.exists(sourcepath+'reprojectedcontinuum.fits'):
        pass
    else:
        dummycube=sc.read(datacubes[0])
        dummyslab=dummycube.spectral_slab(216.65*u.GHz,216.7*u.GHz)
        dummymom0=dummyslab.moment0()
        dummyslab.write(sourcepath+'dummyslab_H2S.fits',overwrite=True)
        dummymom0.write(sourcepath+'dummymom0_H2S.fits',overwrite=True)
        
        cntminfile='/orange/adamginsburg/sgrb2/2017.1.00114.S/imaging_results/Sgr_B2_DS_B6_uid___A001_X1290_X46_continuum_merged_12M_robust2_selfcal4_finaliter_feathered_with_bolocam.fits'
        notreproj_cntmimage=fits.open(cntminfile)
        notreproj_restfreq=notreproj_cntmimage[0].header['RESTFRQ']*u.Hz
        cntmbeam=radio_beam.Beam.from_fits_header(notreproj_cntmimage[0].header)
        
        newcont,footprint=reproject.reproject_interp((notreproj_cntmimage[0].data.squeeze(),WCS(notreproj_cntmimage[0].header).celestial),dummymom0.header)
        newcontprimaryHDU=fits.PrimaryHDU(newcont)
        newcontprimaryHDU.header=dummymom0.header
        newcontprimaryHDU.header['BUNIT']='Jy/beam'
        newcontprimaryHDU.header['BTYPE']='1mm continuum flux'
        newcontprimaryHDU.header['RESTFRQ']=notreproj_restfreq.to('Hz').value
        newcontprimaryHDU.header['BMAJ']=cntmbeam.major.to('deg').value
        newcontprimaryHDU.header['BMIN']=cntmbeam.minor.to('deg').value
        #newcontprimaryHDU.header['BPA']=cntmbeam.
        
        newcontHDUList=fits.HDUList([newcontprimaryHDU])
        
        newcontpath=sourcepath+'reprojectedcontinuum.fits'
        newcontHDUList.writeto(newcontpath,overwrite=True)
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

doublet=['15_6+-16_5+vt1','15_6--16_5-vt1']

excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1','15_6-15_7E1vt1','9_6-9_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','8_6-8_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','15_6-15_7E1vt1','16_6-16_7E1vt1'],'DSii':['7_6-7_7E1vt1','9_6-9_7E1vt1','14_6-14_7E1vt1','10_6-10_7E1vt1','13_6-13_7E1vt1','11_6-11_7E1vt1','23_5-22_6E1vt0'],'DSiii':'','DSiv':'','DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1'],'DSVII':'','DSVIII':'','DSIX':'','DS10':'','DS11':'','DSX':''}
restfreq_representativeline={'SgrB2S':217.88650400*u.GHz,'DSi':220.07856100*u.GHz,'DSii':220.07856100*u.GHz,'DSiii':231.28111000*u.GHz,'DSiv':217.88650400*u.GHz,'DSv':220.07856100*u.GHz,'DSVI':220.07856100*u.GHz,'DSVII':220.07856100*u.GHz,'DSVIII':220.07856100*u.GHz,'DSIX':220.07856100*u.GHz,'DS10':231.28111000*u.GHz,'DS11':220.07856100*u.GHz,'DSX':220.07856100*u.GHz}#All taken from Splatalogue;  oldS 218.44006300
representative_filename_base=sourcepath+representativelines[source]+'repline_'
rep_mom1=representative_filename_base+'mom1.fits'
rep_fwhm=representative_filename_base+'fwhm.fits'
rep_slab=representative_filename_base+'slab.fits'
rep_maskedslab=representative_filename_base+'maskedslab.fits'

regiondims=regionparams[3]

if os.path.isfile(rep_mom1):
    print(f'{source} representative line objects already exist.')
    print(f'Grabbing spectralslab from {rep_slab}')
    spectralslab_rep_3sigma=fits.getdata(rep_maskedslab)*u.K
    print(f'Grabbing mom1 from {rep_mom1}')
    velocityfield_representative=fits.getdata(rep_mom1)*u.km/u.s
    print(f'Grabbing fwhm from {rep_fwhm}\n')
    fwhm_representative=fits.getdata(rep_fwhm)*u.km/u.s
else:
    print(f'{representativelines[source]} representative line objects will be computed for {source}.\n')
    pass

targetworldcrds={'SgrB2S':[[0,0,0],[266.8351718,-28.3961210, 0]], 'DSi':[[0,0,0],[266.8316149,-28.3972040,0]], 'DSii':[[0,0,0],[266.8335363,-28.3963158,0]],'DSiii':[[0,0,0],[266.8332758,-28.3969269,0]],'DSiv':[[0,0,0],[266.8323834, -28.3954424,0]],'DSv':[[0,0,0],[266.8321331, -28.3976585, 0]],'DSVI':[[0,0,0],[266.8380037, -28.4050741,0]],'DSVII':[[0,0,0],[266.8426074, -28.4094401,0]],'DSVIII':[[0,0,0],[266.8418408, -28.4118242, 0]],'DSIX':[[0,0,0],[266.8477371, -28.4311386,0]],'DS10':[[0,0,0],[266.8373798, -28.4009340,0]],'DS11':[[0,0,0],[266.8374572, -28.3996894, 0]],'DSX':[[0,0,0],[266.8452950, -28.4282608,0]]}

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
    stdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)#43:52 selects the region size set by the region variables in the cube>core section of the code
    stddata=stdimage[0].data*u.K

    reprojcontfits=fits.open(contpath)
    reprojcont=reprojcontfits[0].data*u.Jy
    reprojcontrestfreq=225*u.GHz#manual addition 11/9/2022, wiggle room w/i GHz
    cntmbeam=radio_beam.Beam.from_fits_header(reprojcontfits[0].header)
    reprojcont_K=reprojcont.to('K',cntmbeam.jtok_equiv(reprojcontrestfreq))
    
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
    
    cube_w=cube.wcs
    stdwcs=WCS(stdimage[0].header)
    
    targetworldcrd=targetworldcrds[source]
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
    
    '''Create representative line to use as masking template'''
    if not os.path.isfile(rep_mom1):
        print(f'Creating {representativelines[source]} representative line data objects for {source}')
        print(f'Opening file {datacubes[representativecubes[source]]}')
        reffreq_repline=restfreq_representativeline[source]/(1+z)
        repcube=sc.read(datacubes[representativecubes[source]])
        
        assert reffreq_repline > min(repcube.spectral_axis) and reffreq_repline < max(repcube.spectral_axis), f'Incorrect representative cube chosen.\nCube limits: {min(repcube.spectral_axis)},{max(repcube.spectral_axis)}\nRepresentative line reference frequency: {reffreq_repline}'
        
        repstdmain=fits.open(stdhome+f'spw{representativecubes[source]}minimize.image.pbcor_noise.fits')
        repstdmain_data=repstdmain[0].data*u.K
        repstdmain_wcs=WCS(repstdmain[0].header)
        repstdmain_cellsize=(np.abs(repstdmain[0].header['CDELT1']*u.deg)).to('arcsec')
        repstdmain_pixcrd=repstdmain_wcs.wcs_world2pix(targetworldcrd,1,ra_dec_order=True)
        repstdmain_xcrd,repstdmain_ycrd=int(round(repstdmain_pixcrd[1][0])),int(round(repstdmain_pixcrd[1][1]))
        
        assert repstdmain_xcrd and repstdmain_ycrd > 0, 'Negative representative std pixel coords'
        
        repstdcutoutsize=round(((float(regiondims)*u.deg)/stdcellsize).to('').value)
        repstdcutout=Cutout2D(repstdmain_data,(repstdmain_xcrd,repstdmain_ycrd), repstdcutoutsize)
        upperfreq=reffreq_repline+velocitytofreq(representativelws[source],reffreq_repline)
        lowerfreq=reffreq_repline-velocitytofreq(representativelws[source],reffreq_repline)
        print(f'Creating spectral slab between {lowerfreq} and {upperfreq}')
        spectralslab_representative=repcube.spectral_slab(lowerfreq,upperfreq)
        spectralslab_representative=spectralslab_representative.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=reffreq_repline)
        rep_3sigmamask=spectralslab_representative >= (3*repstdcutout.data)
        spectralslab_rep_3sigma=spectralslab_representative.with_mask(rep_3sigmamask)
        
        print('\nComputing moment1')
        velocityfield_representative=spectralslab_rep_3sigma.moment1()
        print('\nComputing fwhm')
        fwhm_representative=spectralslab_rep_3sigma.linewidth_fwhm()
        #pdb.set_trace()
        print('Writing objects to file')
        spectralslab_representative.write(rep_slab)
        spectralslab_rep_3sigma.write(rep_maskedslab)
        velocityfield_representative.write(rep_mom1)
        fwhm_representative.write(rep_fwhm)
        print('Continuing to line loops\n')
    else:
        pass

    #pdb.set_trace()
    
    freq_min=freqs[0]*(1+z)
    freq_max=freqs[(len(freqs)-1)]*(1+z)
    
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
    qns=maintable['QNs']
    euks=maintable['EU_K']*u.K
    eujs=euks*k
    degeneracies=sparetable['Upper State Degeneracy']
    log10aijs=maintable['log10_Aij']
    aijs=10**log10aijs*u.Hz
    
    singlecmpntwidth=(0.00485/8)*u.GHz
    linewidth=fwhm_representative[pixycrd,pixxcrd]#representativelws[source]#10*u.km/u.s#8*u.MHz
    linewidth_freq=velocitytofreq(linewidth,restfreq_representativeline[source])
    linewidth_sigma=linewidth_freq/(2*np.sqrt(2*np.log(2)))
    oldwideslabwidth=(15.15*u.MHz)
    originallinewidth=(11231152.36688232*u.Hz/2)#0.005*u.GHz####0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
    nu_offset=oldwideslabwidth-originallinewidth
    linewidth_vel=vradio(singlecmpntwidth,spwrestfreq)#(singlecmpntwidth*c.to(u.km/u.s)/spwrestfreq).to('km s-1')#vradio(linewidth,spw1restfreq)
    #slicedqns=[]
    
    pixeldict={}
    transitiondict={}
    linelooplte(lines,linewidth,len(lines),qns)
    spwdict.update([(images[imgnum],transitiondict)])
    tempkeys=list(spwdict[images[imgnum]].keys())
    
    kstdimgpath=stdpath+f'{images[imgnum]}fluxstd.fits'
    if os.path.isfile(kstdimgpath):
        print(f'{images[imgnum]} brightness std image already exists')
        spwstdarray=fits.getdata(kstdimgpath)*u.K#stdcutout.data
        #kkmsstdarray=fits.getdata(kkmsstdimgpath)*u.K*u.km/u.s
        #print(f'Retrieved integrated intensity std data from {kkmsstdimgpath}\n')
    else:
        print(f'Start {images[imgnum]} std calculations')
        spwstdarray=stdcutout.data
        print('Set Primary HDU')
        khdu=fits.PrimaryHDU(spwstdarray.value)
        '''This transmoment0 file has intensity (K km/s) units'''
        if len(tempkeys) == 0:
            print(f'No transitions detected in this spw ({images[imgnum]})')
            transmomslab=cube.spectral_slab((lines[0]-linewidth_freq),(lines[0]+linewidth_freq))
            transmoment0=transmomslab.moment0()
            transmom0header=transmoment0.header
            print(f'Set transmoment0 to moment0 from {(lines[0]+linewidth_freq).to("GHz")} to {(lines[0]-linewidth_freq).to("GHz")}')
        else:
            transmoment0=fits.open(spwdict[images[imgnum]][tempkeys[0]]['filename'])
            transmom0header=transmoment0[0].header
            print(f'Set header from {spwdict[images[imgnum]][tempkeys[0]]["filename"]}')
        khdu.header=transmom0header
        khdu.header['BUNIT']='K'
        print('Wrapping Primary HDU in HDUList')
        khdul=fits.HDUList([khdu])
        print(f'Writing to {kstdimgpath}')
        khdul.writeto(kstdimgpath,overwrite=True)
        print(f'{images[imgnum]} std calculations complete.\n')
    #transitiondict.update({'restfreq':spwrestfreq})
    #,('pixel_0',(pixycrd,pixxcrd))])
    kstddict.update([(images[imgnum],spwstdarray)])
    #kkmsstddict.update([(images[imgnum],kkmsstdarray)])
    print(f'Finished loop for {images[imgnum]}\n')
    #masterqns.append(slicedqns)
    #pdb.set_trace()

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
intensities,t_brights=brightnessTandintensities(spwdict)

print(intensityerror)

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
master_transkeys=[]
spws_with_detections=[]
for key in spwdictkeys:
    transdict=spwdict[key]#Dictionary of transitions in a given spw
    #print(f'transdict: {transdict}')
    transitionkeys=list(spwdict[key])#Trandition QNs
    master_transkeys.append(transitionkeys)
    if len(transitionkeys) > 0:
        spws_with_detections.append(transdict)
    else:
        pass
    for transkey in range(len(transitionkeys)):
        print(f'Transition: {transitionkeys[transkey]}/Nupper array z-coord: {pixelzcoord_nupper}')
        nupperimage_filepath=nupperpath+'CH3OH~'+transitionkeys[transkey]+'data.fits'
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
                        tempnupper=lte_molecule.nupper_of_kkms(intensities[transitionkeys[transkey]][y,x],transdict[transitionkeys[transkey]]['freq'],transdict[transitionkeys[transkey]]['aij'])#kkmsstddict[key][y,x])
                        tempnuerr=tempnupper*(transdict[transitionkeys[transkey]]['mom0err'][y,x]/intensities[transitionkeys[transkey]][y,x])
                        #pdb.set_trace()
                        n_us.append((tempnupper.value)/transdict[transitionkeys[transkey]]['degen'])
                        n_uerr.append((tempnuerr.value)/transdict[transitionkeys[transkey]]['degen'])
                        #pdb.set_trace()
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

print('Setting up and executing model fit')
texmap=np.empty((testyshape,testxshape))
ntotmap=np.empty((testyshape,testxshape))

ntoterrmap=np.empty((testyshape,testxshape))
texerrormap=np.empty((testyshape,testxshape))

texsigclipmap=np.empty((testyshape,testxshape))
ntotsigclipmap=np.zeros((testyshape,testxshape))
texsnrmap=np.empty((testyshape,testxshape))
ntotsnrmap=np.zeros((testyshape,testxshape))
numtransmap=np.empty((testyshape,testxshape))
degensforfit=[]
snr=3

fitdict={}

for y in range(testyshape):
    print(f'Start Row {y} fitting')
    for x in range(testxshape):
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        nuperrors=[]
        qnsfitted=[]
        excludedqns=[]
        excludednuppers=[]
        excludedeuks=[]
        for zed in range(testzshape):
            if nugsmap[y,x,zed] <= 0 or np.isnan(nugsmap[y,x,zed]):
                continue
            elif nugsmap[y,x,zed]/nugserrormap[y,x,zed] == 1:
                #print(f'Excluded line detected: {masterqns[z]}')
                #print('Appending to exclusion lists')
                #excludedqns.append(masterqns[zed])
                #excludednuppers.append(nugsmap[y,x,zed])
                #excludedeuks.append(mastereuks[zed])
                pass
            else:
                nupperstofit.append(nugsmap[y,x,zed])
                eukstofit.append(mastereuks[zed])
                nuperrors.append(nugserrormap[y,x,zed])
                qnsfitted.append(masterqns[zed])
                degensforfit.append(ordereddegens[zed])
        numtransmap[y,x]=len(nupperstofit)        
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texsnrmap[y,x]=np.nan
            texsigclipmap[y,x]=obsTex
            texerrormap[y,x]=np.nan
            ntoterrmap[y,x]=np.nan
            ntotsigclipmap[y,x]=np.nan
            ntotsnrmap[y,x]=np.nan
        else:
            #log10nuerr=[]
            errstofit=[]
            log10variances=[]
            
            for num in range(len(nupperstofit)):
                templog10=(1/nupperstofit[num])*nuperrors[num]
                temperrfit=1/templog10
                #log10nuerr.append(templog10)
                errstofit.append(temperrfit)
                log10variances.append(templog10**2)
            
            #pdb.set_trace()
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit),weights=errstofit)
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            obsTrot=-np.log10(np.e)/(fit_lin.slope)
            qrot_partfunc=qrot(obsTrot)
            
            obsNtot=qrot_partfunc*10**(fit_lin.intercept)
            
            A=np.stack((eukstofit,np.ones_like(eukstofit)),axis=1)
            C=np.diagflat(log10variances)
            atc_1a=np.dot(np.dot(A.T, np.linalg.inv(C)), A)
            if np.linalg.det(atc_1a) == 0:
                print(f'Singular C matrix detected in pixel {y,x}')
                m_unc=np.nan
            else:
                covmat = np.linalg.inv(atc_1a)
                m_unc = covmat[0,0]**0.5
                b_unc = covmat[1,1]**0.5
            
            dobsTrot=np.abs(np.abs(m_unc/fit_lin.slope)*obsTrot*u.K)
            dobsNtot=np.abs(qrot_partfunc*10**(fit_lin.intercept)*(np.log(10)*b_unc))*u.cm**-2
            sigTrot=(obsTrot*u.K/dobsTrot).to('')
            sigNtot=(obsNtot*u.cm**-2/dobsNtot).to('')
            
            texmap[y,x]=obsTrot
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTrot.to('K').value
            ntoterrmap[y,x]=dobsNtot.value
            texsnrmap[y,x]=sigTrot
            ntotsnrmap[y,x]=sigNtot
            
            if sigTrot >= snr:
                texsigclipmap[y,x]=obsTrot
            else:
                texsigclipmap[y,x]=np.nan
                
            if sigNtot >= snr:
                '''
                if obsNtot >= 1e29 or dobsNtot.value <= 1:
                    ntotsigclipmap[y,x]=np.nan
                else:
                '''
                ntotsigclipmap[y,x]=obsNtot
            #elif np.isnan(dobsNtot) or dobsNtot == 0:
            #     ntotsigclipmap[y,x]=np.nan
            else:
                ntotsigclipmap[y,x]=np.nan

detectnum=0
transmaskarr=np.ma.masked_where(numtransmap<=detectnum,texsigclipmap)
transmasktexmap=transmaskarr.filled(fill_value=np.nan)#np.array(np.ma.masked_where(numtransmap<detectnum,texsigclipmap))
            
transmoment0=fits.open(spws_with_detections[0][master_transkeys[0][0]]['filename'])
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
print(f'Saving raw ntot map at {sourcepath+"ntotmap_allspw_withnans_weighted_useintercept.fits"}\n')
hdulntot.writeto(sourcepath+'ntotmap_allspw_withnans_weighted_useintercept.fits',overwrite=True)

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

primaryhduntoterr=fits.PrimaryHDU(ntoterrmap)
primaryhduntoterr.header=transmom0header
primaryhduntoterr.header['BTYPE']='Total column density error'
primaryhduntoterr.header['BUNIT']='cm-2'
hdulntoterr=fits.HDUList([primaryhduntoterr])
ntoterrpath=sourcepath+f"ntoterrmap_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntoterr map at {ntoterrpath}\n')
hdulntoterr.writeto(ntoterrpath,overwrite=True)

primaryhduntotsig=fits.PrimaryHDU(ntotsigclipmap)
primaryhduntotsig.header=transmom0header
primaryhduntotsig.header['BTYPE']='Total column density'
primaryhduntotsig.header['BUNIT']='cm-2'
hdulntotsig=fits.HDUList([primaryhduntotsig])
ntotsigpath=sourcepath+f"ntotmap_allspw_withnans_weighted_useintercept_{snr}sigma.fits"
print(f'Saving sigmaclip ntot map at {ntotsigpath}\n')
hdulntotsig.writeto(ntotsigpath,overwrite=True)

primaryhduntotsnr=fits.PrimaryHDU(ntotsnrmap)
primaryhduntotsnr.header=transmom0header
primaryhduntotsnr.header['BTYPE']='Total column density SNR'
primaryhduntotsnr.header['BUNIT']='cm-2/cm-2'
hdulntotsnr=fits.HDUList([primaryhduntotsnr])
ntotsnrpath=sourcepath+f"ntotmap_snr_allspw_withnans_weighted_useintercept.fits"
print(f'Saving ntot snr map at {ntotsnrpath}\n')
hdulntotsnr.writeto(ntotsnrpath,overwrite=True)

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
np.savetxt(sourcepath+'mastereuksqnsfreqsdegens.txt',eukqns,fmt='%s',header=f'Methanol transitions, excitation temperatures, and degeneracies used in this folder. Temperatures in units of K, frequencies are redshifted ({z}/{(z*c).to("km s-1")}) and in Hz.\nExcluded lines: {excludedlines[source]}')

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

saveimghome=overleafpath+f'{source}/{sourcelocs[source]}/'
saveimgpath=saveimghome+f'texmap_{detectnum}transmask_{snr}sigma_allspw_withnans_weighted.png'

if not os.path.exists(saveimghome):
    os.makedirs(saveimghome)
else:
    print(f"Figure directory {saveimghome} already exists.")
    pass

plottexwcs=WCS(plottexhdu)

sliced=['x','y']
ax=plt.subplot(projection=plottexwcs,slices=sliced)
plt.rcParams['figure.dpi'] = 150

ra=ax.coords[0]
dec=ax.coords[1]

vmaxdict={'SgrB2S':525,'DSi':320,'DSii':224,'DSiii':300,'DSiv':312,'DSv':280,'DSVI':377,'DSVII':248,'DSVIII':225,'DSIX':215}
if source in vmaxdict.keys():
    sourcevmax=vmaxdict[source]
    plottedtex=ax.imshow(plottexhdu.data,vmax=sourcevmax,cmap=colormap)#,vmax=1000,vmin=10)
else:
    plottedtex=ax.imshow(plottexhdu.data,cmap=colormap)


scaledict={'SgrB2S':5000*u.AU,'DSi':5000*u.AU,'DSii':2000*u.AU,'DSiii':2000*u.AU,'DSiv':2000*u.AU,'DSv':2000*u.AU,'DSVI':5000*u.AU,'DSVII':5000*u.AU,'DSVIII':5000*u.AU,'DSIX':5000*u.AU,'DS10':5000*u.AU,'DSX':5000*u.AU}
scale=scaledict[source]
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
elif source == 'DSii':
    print(f'Scalebar source: DSii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.1115 -28:23:47.596', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.1087 -28:23:48.634', unit=(u.hour,u.deg), frame='icrs')
    pass
elif source == 'DSiii':
    print(f'Scalebar source: DSiii')
    make_scalebar(ax, coordinates.SkyCoord('17:47:20.0577 -28:23:49.912', unit=(u.hour, u.deg), 
                                       frame='icrs'),
              length=lenn, label=f'{int(scale.value)} AU', fontsize=12, 
              text_offset=0.02*u.arcsec)
    crd=coordinates.SkyCoord('17:47:20.0549 -28:23:49.950', unit=(u.hour,u.deg), frame='icrs')
    
else:
    print(f'Need to make {source} scalebar!!!')
    pass

ra.set_axislabel('RA (J2000)',fontsize=14,minpad=0.9)
ra.set_ticklabel(exclude_overlapping=True)
dec.set_axislabel('Dec (J2000)',fontsize=14,minpad=-0.7)
ax.tick_params(fontsize=14)
cbar1=plt.colorbar(plottedtex,pad=0,label='K')#plt.colorbar()
plt.savefig(saveimgpath,overwrite=True)
plt.show()

print('Cube-->core-->Texmap complete.')
