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
    
qrot_partfunc=Q_rot_asym(testT).to('')

images=['spw0','spw2','spw1','spw3']
imgnum=0



print('Accessing data cube')
home=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/Mom0/field1/{images[imgnum]}stat/'#Make sure to include slash after path
fname=f'/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_{images[imgnum]}_cube_minimize_hasbeams.image_line.fits'
readstart=time.time()
cube=sc.read(fname)
readelapsed=time.time()-readstart
print('Cube read in')
print(time.strftime("%H:%M:%S", time.gmtime(readelapsed)))

#cube=cube.rechunk(save_to_tmp_dir=True)
header=fits.getheader(fname)

print('Acquiring cube rest frequency and computing target pixel coordinates')
spwrestfreq=header['RESTFRQ']*u.Hz 
freqs=cube.spectral_axis#Hz
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
freq_max=freqs[0]*(1+z)#215*u.GHz
#print(freq_max)
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
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
                                
'''Loop through a given list of lines (in Hz), computing and saving moment0 maps of the entered data cube'''
def linelooplte(line_list,line_width,iterations,quantum_numbers):
    print('\nlinelooperLTE...')
    print('Grab cube and reference pixel')
    targetpixjybeam=cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]
    #targetpixjybeam=targetpixjybeam.mask_channels(np.isnan(targetpixjybeam)==False)
    cubebeams=(cube.beams.value)*u.sr/u.beam
    print('Convert from Jy/beam to K')
    targetspec_K=JybeamtoK(cubebeams,targetpixjybeam)
    print('Compute cube brightness temperature stddev')
    targetspecK_stddev=np.nanstd(targetspec_K)
    #pdb.set_trace()
    for i in range(iterations):
        line=line_list[i]#*u.Hz
        print(f'Start {quantum_numbers[i]} moment0 procedure')

        nu_upper=line+line_width
        nu_lower=line-line_width
        print(f'Make spectral slab between {nu_lower} and {nu_upper}')
        slabstart=time.time()
        slab=cube.spectral_slab(nu_upper,nu_lower)
        slabend=time.time()-slabstart
        print(f'{quantum_numbers[i]} spectral slab done in {time.strftime("%H:%M:%S", time.gmtime(slabend))}')
        
        slabbeams=(slab.beams.value)*u.sr/u.beam
        #print(f'slabbeams: {slabbeams}')
        slab_K=JybeamtoK(slabbeams,slab[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))])
        #print(f'slab_K: {slab_K}')
        mulu2=(mulu(aijs[i],line)).to('cm5 g s-2')
        linewidth_vel=vradio(singlecmpntwidth,line)
        tbthick=Tbthick(testntot,line,linewidth_vel,mulu2,degeneracies[i],qrot_partfunc,eujs[i],testT).to('K')
        print('LTE params calculated')
        print(f'tbthick: {tbthick} targetspecK_stddev: {targetspecK_stddev}')
        print('Slicing quantum numbers')
        transition=qn_replace(quantum_numbers[i])
        slicedqns.append(transition)
        moment0filename=home+'CH3OH~'+transition+'.fits'
        #print('Done')
        #print('Moment 0')
        if os.path.isfile(moment0filename):
            print(f'{moment0filename} already exists.\nSkipping transition {quantum_numbers[i]}\n')
            #stddevs=np.append(targetspecK_stddev)
            pass
        elif tbthick >= targetspecK_stddev*u.K:
            print('Commence moment0')
            slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spwrestfreq)
            momstart=time.time()
            slabmom0=slab.moment0()
            momend=time.time()-momstart
            print(f'{quantum_numbers[i]} elapsed time: {time.strftime("%H:%M:%S", time.gmtime(momend))}')
            print('\nSaving...')
            #name='test'+str(i)
            slabmom0.write((home+'CH3OH~'+transition+'.fits'),overwrite=True)
            stddevs.append(targetspecK_stddev)
            print('Done\n')
        else:
            print('LTE Model max brightnessTemp below 1sigma threshold')
            print(f'{quantum_numbers[i]} skipped\n')
            pass
    print('lines looped.')

'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    return string

print('Gathering Splatalogue table parameters')    
lines=table['Freq']*10**9*u.Hz/(1+z)
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
slicedqns=[]

stddevs=[]

linelooplte(lines,linewidth,len(lines),qns)

######

#print(qns)
files=glob.glob(home+'*.fits')
#print(files)

def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix),(xpix)])#Corrects for different pixel counting procedures
    return vals

'''Gathers beam data from cube headers'''    
def beamer(fitsfiles):
    beams=[]
    print(f'in beamer function: {fitsfiles}')
    for filename in fitsfiles:
        hdu=fits.getheader(filename)
        temp=radio_beam.Beam.from_fits_header(hdu).value
        beams.append(temp)
    return beams

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
            comp=(sliced_qns[j]==tempfiles[i][83:])
            if comp==True:
                print(f'{sliced_qns[j]} == {tempfiles[i][55:]}\comp==True')
                unscrambled_qns.append(sliced_qns[j])
                unscrambled_freqs.append(linelist[j])
                unscrambled_eus.append(euks[j]/u.K)
                unscrambled_degs.append(degeneracies[j])
                unscrambled_aijs.append(aijs[j])
                #unscrambled_stddevs.append()
                break
            else: 
                print(f'{sliced_qns[j]} != {tempfiles[i][55:]}')
    return unscrambled_qns,unscrambled_freqs,unscrambled_eus,unscrambled_degs,unscrambled_aijs
beamlist=beamer(files)*u.sr
#print(beamlist)
fluxes=fluxvalues(int(round(targetpixcrd[1][0])),int(round(targetpixcrd[1][1])),files)*u.Jy*u.km/u.s#/u.sr
#print(fluxes)
unscrambledqns,unscrambledfreqs,unscrambledeuks,unscrambleddegs,unscrambledaijs=unscrambler(files,slicedqns,lines)
print(f'files: {files}')
print(f'unscrambledqns: {unscrambledqns}')
print(f'unscrambledfreqs: {unscrambledfreqs}')

'''Places Splatalogue table params, fluxes, and beams into dictionary'''
datadict={}
for i in range(len(fluxes)):
    datadict[i]={'qns':unscrambledqns[i],'freq':unscrambledfreqs[i],'beam':beamlist[i],'flux':fluxes[i],'E_u(K)':unscrambledeuks[i],'degen':unscrambleddegs[i],'aij':unscrambledaijs[i]}

'''Rough estimate: 0.75 arcsec/15pixels >>> 0.05 arsec/pixel >>> 2.424e-7rad/pixel
Taken from DS9 tradiation region'''

'''Compute Kkm/s intensity from datadict'''
def kkms(beams,data_dict):
    intensitylist=[]
    t_bright=[]
    for i in range(len(fluxes)):
        temp=(data_dict[i]['flux']/linewidth_vel).to('Jy')
        #print(temp)
        equiv=u.brightness_temperature(data_dict[i]['freq'])
        #print(equiv)
        jy_sr=temp/beams[i]
        #print(jy_sr)
        conversion=jy_sr.to(u.K,equivalencies=equiv)
        t_bright.append(conversion)
        #print(conversion)
        velflux_T=conversion*linewidth_vel
        #print(velflux_T)
        #print('\n')
        intensitylist.append(velflux_T)
    return intensitylist,t_bright

intensities,t_brights=kkms(beamlist,datadict)

print(t_brights)
    
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

def N_u(nu,Aij,velocityintegrated_intensity_K):#(ntot,qrot,gu,eu_J,T_ex):
    return ((8*np.pi*k*nu**2)/(h*c**3*Aij))*velocityintegrated_intensity_K#ntot/(qrot*np.exp(eu_J/(k*T_ex)))
    
def S_j(j_upper,k_upper):#Works for symmetric tops
    return (j_upper**2-k_upper**2)/(j_upper*(2*j_upper+1))
    
def KtoJ(T):
    return (3/2)*k*T
    
def Ntot_rj_thin_nobg(nu,s,g,q,eu_J,T_ex,vint_intensity):
    #nu=nu
    #T_ex=T_ex
    #T_r=T_r
    return ((3*k)/(8*np.pi**3*nu*mu_a**2*s*R_i))*(q/g)*np.exp((eu_J/(k*T_ex)))*((vint_intensity))#((nu+templatewidth)-(nu-templatewidth)))
    
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
'''Compute N_uppers'''
n_us=[]
for i in range(len(intensities)):
    temp=N_u(datadict[i]['freq'],datadict[i]['aij'],intensities[i])
    n_us.append((temp.to('cm-2')*u.cm**2)/unscrambleddegs[i])
    
linemod=models.Linear1D(slope=1.0,intercept=14)
fit=fitting.LinearLSQFitter()
fit_lin=fit(linemod,unscrambledeuks,np.log10(n_us))
linemod_euks=np.linspace(min(unscrambledeuks),max(unscrambledeuks),100)

obsTex=-np.log10(np.e)/(fit_lin.slope)
obsNtot=qrot_partfunc*10**(np.log10(n_us[0])+fit_lin.slope*unscrambledeuks[0])
#jupper,kupper=jupperfinder(unscrambledqns)
'''
qrots=[]
s_j=[]
for i in range(len(texs)):
    qrots.append(Q_rot_asym(texs[i]))
    s_j.append(S_j(jupper[i],kupper[i]))'''

#testntot=Ntot_rj_thin_nobg(datadict[0]['freq'],s_j[0],datadict[0]['degen'],qrots[0],
#KtoJ(datadict[0]['E_u']),texs[0],intensities[0])
#testnuoverg=N_u(testntot,qrots[0],datadict[0]['degen'],KtoJ(datadict[0]['E_u']),texs[0])
print(n_us[i])
#ntot=Ntot_rj_thin_nobg(lines[0],linewidth,s,9,KtoJ(tex[0])
home2=home.replace('spw0/','')
plt.clf()
plt.scatter(unscrambledeuks,np.log10(n_us))
plt.plot(linemod_euks,fit_lin(linemod_euks),label=(f'obsTex: {obsTex*u.K}\nobsNtot: {obsNtot/u.cm**2}'))
plt.title(f'statcont --continuum {images[imgnum]} CH$_3$OH Rotational Diagram')
plt.xlabel(r'E$_u$ (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
plt.legend()
plt.savefig((home2+'trialrotdiag2.png'),dpi=100,overwrite=True)
plt.show()


