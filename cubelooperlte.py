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

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

sanitytable=Splatalogue.query_lines(100*u.GHz, 700*u.GHz, chemical_name='HCN',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=['JPL'],
                                    show_upper_degeneracy=True, get_query_payload=True)
                                    
print(sanitytable)

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
    print('\nlinelooperLTE...')
    print('Grab cube and reference pixel')
    targetspec_K=cube[:,pixycrd,pixxcrd]
    #targetpixjybeam=targetpixjybeam.mask_channels(np.isnan(targetpixjybeam)==False)
    cubebeams=(cube.beams.value)*u.sr/u.beam
    print('Convert from Jy/beam to K')
    #targetspec_K=targetpixjybeam.to(u.K)#JybeamtoK(cubebeams,targetpixjybeam)
    #pdb.set_trace()
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
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(moment0filename):
            print(f'{moment0filename} already exists.')
            isfilemom0=fits.getdata(moment0filename)*u.K*u.km/u.s
            isfilepixflux=isfilemom0[pixycrd,pixxcrd]
            isfilebeam=beamer(moment0filename)
            
            temptransdict.update([('freq',line),('flux',isfilepixflux),('stddev',targetspecK_stddev),('beam',isfilebeam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i])])
            transitiondict.update({transition:temptransdict})
            masterslicedqns.append(transition)
            mastereuks.append(euks[i].value)
            masterstddevs.append(targetspecK_stddev)
            #transitionfluxlist.append(isfilepixflux)
            #masterfluxes.append(isfilepixflux)
            #transitionbeamlist.append(isfilebeam)
            #masterbeams.append(isfilebeam)
            #masterslicedqns.append(transition)
            #masterlines.append(line)
            #mastereuks.append(euks[i])
            #mastereujs.append(eujs[i])
            #masterdegens.append(degeneracies[i])
            #masterlog10aijs.append(log10aijs[i])
            #masteraijs.append(aijs[i])
            print('\nMaster lists populated for this transition. Proceeding...\n')
            pass
        elif tbthick >= targetspecK_stddev:#*u.K:
            print('Commence moment0')
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
            momstart=time.time()
            slabmom0=slab.moment0()
            momend=time.time()-momstart
            print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
            print('\nSaving...')
            #name='test'+str(i)
            slabmom0.write((home+'CH3OH~'+transition+'.fits'),overwrite=True)
            moment0beam=slabmom0.beam.value*u.sr
            targetpixflux=slabmom0[pixycrd,pixxcrd]#/u.beam#Unsure about the division here, will have to discuss later on
            temptransdict.update([('freq',line),('flux',targetpixflux),('stddev',targetspecK_stddev),('beam',moment0beam),('euk',euks[i]),('eujs',eujs[i]),('degen',degeneracies[i]),('aij',aijs[i])])
            transitiondict.update({transition:temptransdict})
            mastereuks.append(euks[i])
            masterstddevs.append(targetspecK_stddev)
            #transitionbeamlist.append(moment0beam)
            #masterbeams.append(moment0beam)
            #transitionfluxlist.append(targetpixflux)
            #masterfluxes.append(targetpixflux)
            #masterslicedqns.append(transition)
            #masterlines.append(line)
            #mastereuks.append(euks[i])
            #mastereujs.append(eujs[i])
            #masterdegens.append(degeneracies[i])
            #masterlog10aijs.append(log10aijs[i])
            #masteraijs.append(aijs[i])
            print(f'{quantum_numbers[i]} calculations complete.\n')
        else:
            print('LTE Model max brightnessTemp below 1sigma threshold')
            print(f'{quantum_numbers[i]} skipped\n')
            pass
    spectraKdict.update({images[imgnumber]:targetspec_K})
    #masterdict[images[imgnum]]={'beam':transitionbeamlist,'fluxes':transitionfluxlist}
    print('lines looped.\n')

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
        print(temptransdictkeys)
        for i in range(len(temptransdictkeys)):
            if 'restfreq' in temptransdictkeys[i]:
                continue
            else:
                velflux_T=temptransdict[temptransdictkeys[i]]['flux']
                intensitydict.update({temptransdictkeys[i]:velflux_T})
                temp=velflux_T/linewidth_vel
                t_bright.update({temptransdictkeys[i]:temp})
                d_velfluxT=(temptransdict[temptransdictkeys[i]]['stddev']/temp)*velflux_T
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

datacubes=glob.glob('/blue/adamginsburg/d.jeff/imaging_results/field1core1box/*.fits')
images=[]#'spw0','spw2','spw1','spw3']
for files in datacubes:
    images.append(files[57:61])
assert 'spw0' in images, f'image name list does not match spw# format'
fieldpath=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/Mom0/field1core1testbox/'

spwdict={}
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

for imgnumber in range(len(datacubes)):
    print(f'Accessing data cube {datacubes[imgnumber]}')
    assert images[imgnumber] in datacubes[imgnumber], f'{images[imgnumber]} not in filename {datacubes[imgnumber]}'
    home=fieldpath+f'{images[imgnumber]}/'#Make sure to include slash after path
    readstart=time.time()
    cube=sc.read(datacubes[imgnumber])
    readelapsed=time.time()-readstart
    print(f'Cube read in {time.strftime("%H:%M:%S", time.gmtime(readelapsed))}')
    
    #cube=cube.rechunk(save_to_tmp_dir=True)
    header=fits.getheader(datacubes[imgnumber])
    
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
    table = utils.minimize_table(Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=[linelist],
                                    show_upper_degeneracy=True))
    '''Needed for upper state degeneracies'''                                
    sparetable=Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ',
                                    energy_max=1840, energy_type='eu_k',
                                    line_lists=[linelist],
                                    show_upper_degeneracy=True)
                                    
    
    print('Gathering Splatalogue table parameters')    
    lines=table['Freq']*10**9*u.Hz/(1+z)
    #masterlines.append(lines)
    #vel_lines=vradio(lines,spw1restfreq)
    qns=table['QNs']
    euks=table['EU_K']*u.K
    eujs=[]
    for eupper_K in euks:
        eujs.append(KtoJ(eupper_K))
    degeneracies=sparetable['Upper State Degeneracy']
    log10aijs=table['log10_Aij']
    aijs=10**log10aijs*u.Hz
    
    '''
    for i in range(len(test)):
        plt.axvline(x=test[i],color='red')
    plt.show()
    '''
    singlecmpntwidth=(0.00485/8)*u.GHz
    linewidth=11231152.36688232*u.Hz#0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
    linewidth_vel=vradio(singlecmpntwidth,spwrestfreq)#(singlecmpntwidth*c.to(u.km/u.s)/spwrestfreq).to('km s-1')#vradio(linewidth,spw1restfreq)
    #slicedqns=[]
    
    transitiondict={}
    linelooplte(lines,linewidth,len(lines),qns)
    #transitiondict.update({'restfreq':spwrestfreq})
    spwdict.update({images[imgnumber]:transitiondict})
    #masterqns.append(slicedqns)
    
######

#print(qns)
#files=glob.glob(fieldpath+'*/*.fits')
#print(files)

#print(beamlist)
#fluxes=fluxvalues(int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1])),files)*u.Jy*u.km/u.s#/u.sr
#print(fluxes)
'''
unscrambledqns,unscrambledfreqs,unscrambledeuks,unscrambleddegs,unscrambledaijs=unscrambler(files,masterqns,masterlines)
print(f'files: {files}')
print(f'unscrambledqns: {unscrambledqns}')
print(f'unscrambledfreqs: {unscrambledfreqs}')
'''

'''Places Splatalogue table params, fluxes, and beams into dictionar
datadict={}
for i in range(len(masterfluxes)):
    datadict[i]={'qns':unscrambledqns[i],'freq':unscrambledfreqs[i],'beam':beamlist[i],'flux':fluxes[i],'E_u(K)':unscrambledeuks[i],'degen':unscrambleddegs[i],'aij':unscrambledaijs[i]}
'''

print('Computing K km/s intensities and K brightness temperatures')
intensityerror=[]
intensities,t_brights=brightnessTandintensities(spwdict)#JybeamtoKkms(spwdict)

print(t_brights)
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
n_us=[]
n_uerr=[]
spwdictkeys=spwdict.keys()
for key in spwdictkeys:
    transdict=spwdict[key]
    transitionkeys=list(spwdict[key])
    for transkey in range(len(transitionkeys)):
        tempnupper,tempnuerr=N_u(transdict[transitionkeys[transkey]]['freq'],transdict[transitionkeys[transkey]]['aij'],intensities[transitionkeys[transkey]],intensityerror[transkey])
        n_us.append((tempnupper.to('cm-2')*u.cm**2).to('')/transdict[transitionkeys[transkey]]['degen'])
        n_uerr.append((tempnuerr.to('cm-2')*u.cm**2).to('')/transdict[transitionkeys[transkey]]['degen'])

print('Setting up and executing model fit')
linemod=models.Linear1D(slope=1.0,intercept=14)
fit=fitting.LinearLSQFitter()
fit_lin=fit(linemod,mastereuks,np.log10(n_us))
linemod_euks=np.linspace(min(mastereuks),max(mastereuks),100)
print('Model fit complete')

print('Compute obsTex and obsNtot')
obsTex=-np.log10(np.e)/(fit_lin.slope)
obsNtot=qrot_partfunc*10**(np.log10(n_us[0])+fit_lin.slope*mastereuks[0])

log10nuerr=[]
for num in range(len(n_us)):
    templog10=(1/n_us[num])*n_uerr[num]
    log10nuerr.append(templog10)
    
dobsTex=(mastereuks[0]*u.K*np.log(10)*np.log(np.e))/(np.log(n_us[0]/spwdict['spw2']['10_2--9_3-vt0']['degen'])-np.log(obsNtot/qrot_partfunc))**2
print(f'dobsTex: {dobsTex}')
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


