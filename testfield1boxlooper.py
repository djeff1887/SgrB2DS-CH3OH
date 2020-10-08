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

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

sanitytable1 = utils.minimize_table(Splatalogue.query_lines(200*u.GHz, 300*u.GHz, chemical_name=' HCN ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=['JPL'],
                                    show_upper_degeneracy=True))

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

print('Setting input LTE parameters')
testT=500*u.K
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
    targetspecK_stddev=targetspec_K.std()#np.nanstd(targetspec_K)
    transitionbeamlist=[]
    transitionfluxlist=[]
    for i in range(iterations):
        print(f'Start {quantum_numbers[i]} moment0 procedure')
        temptransdict={}
        line=line_list[i]#*u.Hz
        nu_upper=line+line_width
        nu_lower=line-line_width
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        slabstart=time.time()
        slab=cube.spectral_slab(nu_upper,nu_lower)
        slabend=time.time()-slabstart
        print(f'{quantum_numbers[i]} spectral slab done in {time.strftime("%H:%M:%S", time.gmtime(slabend))}')
        
        slabbeams=(slab.beams.value)*u.sr/u.beam
        #print(f'slabbeams: {slabbeams}')
        slab_K=slab[:,pixycrd,pixxcrd]#JybeamtoK(slabbeams,)
        #print(f'slab_K: {slab_K}')
        mulu2=(mulu(aijs[i],line)).to('cm5 g s-2')
        linewidth_vel=vradio(singlecmpntwidth,line)
        tbthick=Tbthick(testntot,line,linewidth_vel,mulu2,degeneracies[i],qrot_partfunc,eujs[i],testT).to('K')
        print('LTE params calculated')
        print(f'tbthick: {tbthick} targetspecK_stddev: {targetspecK_stddev}')
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        moment0filename=home+'CH3OH~'+transition+'.fits'
        slabfilename=slabpath+'CH3OH~'+transition+'_slab.fits'
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(moment0filename):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(moment0filename)*u.K*u.km/u.s
            #isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(moment0filename)
            temptransdict.update([('freq',line),('flux',isfilemom0),('stddev',np.nanstd(isfilemom0)),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename)])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            print('\nDictionaries populated for this transition.')
            if os.path.isfile(slabfilename):
                print('Proceeding...\n')
                pass
            else:
                slab.write(slabfilename)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            pass
        elif tbthick >= targetspecK_stddev:#*u.K:
            print('Commence moment0')
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=lines[i])#spwrestfreq)
            momstart=time.time()
            slabmom0=slab.moment0()
            momend=time.time()-momstart
            print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
            print('\nSaving...')
            #name='test'+str(i)
            slabmom0.write((moment0filename),overwrite=True)
            moment0beam=slabmom0.beam.value*u.sr
            targetpixflux=slabmom0[pixycrd,pixxcrd]
            temptransdict.update([('freq',line),('flux',slabmom0),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i]),('filename',moment0filename)])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i])
            masterstddevs.append(targetspecK_stddev)
            print(f'{quantum_numbers[i]} calculations complete.\n')
            if os.path.isfile(slabfilename):
                print(f'Spectral slab {slabfilename} already exists.\nProceeding...\n')
                pass
            else:
                slab.write(slabfilename)
                print(f'Slab written to {slabfilename}. Proceeding...\n')
            pass
        else:
            print('LTE Model max brightnessTemp below 1sigma threshold')
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
    
def N_u(nu,Aij,velocityintegrated_intensity_K,velint_intK_err):#(ntot,qrot,gu,eu_J,T_ex):
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
     
qrot_partfunc=Q_rot_asym(testT).to('')

sanitytable2 = utils.minimize_table(Splatalogue.query_lines(200*u.GHz, 300*u.GHz, chemical_name=' HCN ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=['JPL'],
                                    show_upper_degeneracy=True))

datacubes=glob.glob('/blue/adamginsburg/d.jeff/imaging_results/field1core1box2/*.fits')
images=[]#'spw0','spw2','spw1','spw3']
for files in datacubes:
    images.append(files[58:62])#[57:61])
assert 'spw1' in images, f'image name list does not match spw# format'
sourcepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/'
nupperpath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/nuppers/'
stdpath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/errorimgs/std/'
slabpath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/spectralslabs/km_s/'
picklepath='/blue/adamginsburg/d.jeff/SgrB2DSreorg/field1/CH3OH/SgrB2S/testbox2dict.obj'

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
    targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    
    pixxcrd,pixycrd=int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1]))
    print(f'x: {pixxcrd}/y: {pixycrd}')
    #data=cube_unmasked[:,368,628]#[:,383,649]#Jy*km/s
    #test=cube_unmasked[:,383:413,649:679]
    
    #print(np.shape(data))
    
    #sum=np.sum(test[:,0:,0:])
    #print(np.shape(sum))
    
    #spec=np.stack((freqs,data),axis=1)
    #print(np.shape(spec))
    
    #plt.plot(spec[:,0],spec[:,1])
    #plt.show()
    
    z=0.000234806#0.0002333587
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
    lines=maintable['Freq']*10**9*u.Hz/(1+z)#Redshifted to Sgr B2
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
    linewidth=0.005*u.GHz#11231152.36688232*u.Hz#0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
    linewidth_vel=vradio(singlecmpntwidth,spwrestfreq)#(singlecmpntwidth*c.to(u.km/u.s)/spwrestfreq).to('km s-1')#vradio(linewidth,spw1restfreq)
    #slicedqns=[]
    
    pixeldict={}
    transitiondict={}
    linelooplte(lines,linewidth,len(lines),qns)
    spwdict.update([(images[imgnum],transitiondict)])
    tempkeys=list(spwdict[images[imgnum]].keys())
    
    
    stdfitsimgpath=stdpath+f'{images[imgnum]}intensitystd.fits'
    if os.path.isfile(stdfitsimgpath):
        print(f'{images[imgnum]} brightness std image already exists')
        spwstdarray=fits.getdata(stdfitsimgpath)*u.K
        kkmsstdarray=spwstdarray*linewidth_vel
        print(f'Retrieved brightness std data from {stdfitsimgpath}\n')
    else:
        print(f'Start {images[imgnum]} std calculations')
        spwstdarray,kkmsstdarray=pixelwisestd(cube)
        print('Set Primary HDU')
        hdu=fits.PrimaryHDU(kkmsstdarray.value)
        transmoment0=fits.open(spwdict[images[imgnum]][tempkeys[0]]['filename'])
        transmom0header=transmoment0[0].header
        print(f'Set header from {spwdict[images[imgnum]][tempkeys[0]]["filename"]}')
        hdu.header=transmom0header
        #hdu.header['BUNIT']='K km s-1'
        print('Wrapping Primary HDU in HDUList')
        hdul=fits.HDUList([hdu])
        print(f'Writing to {stdfitsimgpath}')
        hdul.writeto(stdfitsimgpath)
        print(f'{images[imgnum]} std calculations complete.\n')
    #transitiondict.update({'restfreq':spwrestfreq})
    #,('pixel_0',(pixycrd,pixxcrd))])
    kstddict.update([(images[imgnum],spwstdarray)])
    kkmsstddict.update([(images[imgnum],kkmsstdarray)])
    print(f'Finished loop for {images[imgnum]}')
    #masterqns.append(slicedqns)



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
testyshape=60
testxshape=60
testzshape=len(masterslicedqns)
nugsmap=np.empty(shape=(testyshape,testxshape,testzshape))
nugserrormap=np.empty(shape=(testyshape,testxshape,testzshape))
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
      
print('pixels looped.')  

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

fitdict={}
pdb.set_trace()
for y in range(testyshape):
    print(f'Start Row {y} fitting')
    for x in range(testxshape):
        linemod=models.Linear1D(slope=1.0,intercept=14)
        fit=fitting.LinearLSQFitter()
        nupperstofit=[]
        eukstofit=[]
        for z in range(testzshape):
            if nugsmap[y,x,z] <= 0:# or np.isnan(nugsmap[y,x,z]):
                continue
            else:
                nupperstofit.append(nugsmap[y,x,z])
                eukstofit.append(mastereuks[z])
        if len(nupperstofit)==0:
            obsTex=np.nan
            obsNtot=np.nan
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
        else:
            fit_lin=fit(linemod,eukstofit,np.log10(nupperstofit))
            linemod_euks=np.linspace(min(eukstofit),max(mastereuks),100)
            #print('Model fit complete')
            #print('Compute obsTex and obsNtot')
            obsTex=-np.log10(np.e)/(fit_lin.slope)
            obsNtot=qrot_partfunc*10**(np.log10(nugsmap[y,x,0])+fit_lin.slope*eukstofit[0])
            dobsTex=(eukstofit[0]*u.K*np.log(10)*np.log(np.e))/(np.log(nupperstofit[0]/spwdict['spw2']['10_2--9_3-vt0']['degen'])-np.log(obsNtot/qrot_partfunc))**2
            sigTex=(obsTex/dobsTex).to('')
            
            texmap[y,x]=obsTex
            ntotmap[y,x]=obsNtot
            texerrormap[y,x]=dobsTex.to('K').value
            texsnrmap[y,x]=sigTex
            if sigTex >= 3:
                texsigclipmap[y,x]=obsTex
            else:
                texsigclipmap[y,x]=np.nan
            
transmoment0=fits.open(transdict[transitionkeys[transkey]]['filename'])
transmom0header=transmoment0[0].header
            
primaryhdutex=fits.PrimaryHDU(texmap)
primaryhdutex.header=transmom0header
primaryhdutex.header['BTYPE']='Excitation temperature'
primaryhdutex.header['BUNIT']='K'
hdultex=fits.HDUList([primaryhdutex])
hdultex.writeto(sourcepath+'texmap_allspw.fits',overwrite=True)

primaryhduntot=fits.PrimaryHDU(ntotmap)
primaryhduntot.header=transmom0header
primaryhduntot.header['BTYPE']='Total column density'
primaryhduntot.header['BUNIT']='cm-2'
hdulntot=fits.HDUList([primaryhduntot])
hdulntot.writeto(sourcepath+'ntotmap_allspw.fits',overwrite=True)

primaryhdutexerr=fits.PrimaryHDU(texerrormap)
primaryhdutexerr.header=transmom0header
primaryhdutexerr.header['BTYPE']='Excitation temperature'
primaryhdutexerr.header['BUNIT']='K'
hdultexerror=fits.HDUList([primaryhdutexerr])
hdultexerror.writeto(sourcepath+'texmap_error_allspw.fits',overwrite=True)

primaryhdutexclip=fits.PrimaryHDU(texsigclipmap)
primaryhdutexclip.header=transmom0header
primaryhdutexclip.header['BTYPE']='Excitation temperature'
primaryhdutexclip.header['BUNIT']='K'
hdultexclip=fits.HDUList([primaryhdutexclip])
hdultexclip.writeto(sourcepath+'texmap_3sigma_allspw.fits',overwrite=True)

primaryhdutexsnr=fits.PrimaryHDU(texsnrmap)
primaryhdutexsnr.header=transmom0header
primaryhdutexsnr.header['BTYPE']='Excitation temperature SNR'
primaryhdutexsnr.header['BUNIT']=''
hdultexsnr=fits.HDUList([primaryhdutexsnr])
hdultexsnr.writeto(sourcepath+'texmap_snr_allspw.fits',overwrite=True)

print('Finished.')

'''
log10nuerr=[]
for num in range(len(n_us)):
    templog10=(1/n_us[num])*n_uerr[num]
    log10nuerr.append(templog10)
plt.clf()
print('Begin plotting')
plt.errorbar(mastereuks,np.log10(n_us),yerr=log10nuerr,fmt='o')
plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {obsTex*u.K}\nobsNtot: {obsNtot/u.cm**2}'))
plt.title(f'field1 {spwdict.keys()} CH$_3$OH Rotational Diagram')
plt.xlabel(r'E$_u$ (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
plt.legend()
plt.savefig((fieldpath+'rotdiag.png'),dpi=100,overwrite=True)
plt.show()
print('cubes loopered.')
'''

