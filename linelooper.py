import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam

'''Collect constants and CH3OH-specific quantum parameters'''
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
        
home='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS-CH3OH/Mom0/field1/spw2'#Make sure to include slash after path
fname='/blue/adamginsburg/d.jeff/imaging_results/SgrB2DS_field1_spw2_cube_medsub.image.fits'
cube=sc.read(fname)
header=fits.getheader(fname)


spw1restfreq=header['RESTFRQ']*u.Hz 
freqs=cube.spectral_axis#Hz
velcube=cube.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spw1restfreq)
#print(velcube.spectral_axis)
cube_unmasked=velcube.unmasked_data

cube_w=cube.wcs
targetworldcrd=[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]
targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)

#data=cube_unmasked[:,368,628]#[:,383,649]#Jy*km/s
#test=cube_unmasked[:,383:413,649:679]

#print(np.shape(data))

#sum=np.sum(test[:,0:,0:])
#print(np.shape(sum))

#spec=np.stack((freqs,data),axis=1)
#print(np.shape(spec))

#plt.plot(spec[:,0],spec[:,1])
#plt.show()

z=0.0002333587
freq_max=freqs[0]*(1+z)#215*u.GHz
#print(freq_max)
freq_min=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
#print(freq_min)
linelist='JPL'

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
def lineloop(line_list,line_width,iterations,quantum_numbers):
    for i in range(iterations):
        line=line_list[i]#*u.Hz
        print(f'Computing {quantum_numbers[i]} moment0')
        nu_upper=line+line_width
        nu_lower=line-line_width
        #print('Done')
        #print('Make spectral slab')
        slab=cube.spectral_slab(nu_upper,nu_lower)
        slab=slab.with_spectral_unit((u.km/u.s),velocity_convention='radio',rest_value=spw1restfreq)
        #print('Done')
        #print('Moment 0')
        slabmom0=slab.moment0()
        print('Saving...')
        transition=qn_replace(quantum_numbers[i])
        slicedqns.append(transition)
        #name='test'+str(i)
        slabmom0.write((home+'CH3OH~'+transition+'.fits'),overwrite=True)
        print('Done')

'''Replaces unwanted characters from the QN table for use in file names'''
def qn_replace(string):
    string=string.replace('=','')
    string=string.replace('(','_')
    string=string.replace(')','')
    return string
    
'''Converts given line list in frequency to radio velocity'''
def vradio(frequency,rest_freq):
    velocity_list=c.to(u.km/u.s)*((rest_freq-frequency)/rest_freq)
    return velocity_list.to('km s-1')
    
lines=table['Freq']*10**9*u.Hz/(1+z)
#vel_lines=vradio(lines,spw1restfreq)
qns=table['QNs']
eus=table['EU_K']*u.K
degeneracies=sparetable['Upper State Degeneracy']
log10aijs=table['log10_Aij']

'''
for i in range(len(test)):
    plt.axvline(x=test[i],color='red')
plt.show()
'''
linewidth=0.5*0.0097*u.GHz#from small line @ 219.9808GHz# 0.0155>>20.08km/s 
linewidth_vel=(linewidth*c.to(u.km/u.s)/spw1restfreq).to('km s-1')#vradio(linewidth,spw1restfreq)
slicedqns=[]

print('\nlinelooper...')
lineloop(lines,linewidth,len(lines),qns)
print('lines looped.\n')

######

#print(qns)
files=glob.glob(home+'*.fits')
#print(files)

def fluxvalues(xpix,ypix,filenames):
    vals=[]
    for i in filenames:
        data=fits.getdata(i)
        vals.append(data[(ypix-1),(xpix-1)])#Corrects for different pixel counting procedures
    return vals

'''Gathers beam data from cube headers'''    
def beamer(fitsfiles):
    beams=[]
    for i in fitsfiles:
      print(i)
      hdu=fits.getheader(i)
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
        #print(f'filename: {filenames[i]}')
        tempfiles[i]=tempfiles[i].replace('.fits','')
        for j in range(len(sliced_qns)):
            #print(f'sliced_qns: {sliced_qns[j]}')
            #print(f'comparison qns: {tempfiles[i][55:]}')
            comp=(sliced_qns[j]==tempfiles[i][55:])
            if comp==True:
                print('comp==True')
                unscrambled_qns.append(sliced_qns[j])
                unscrambled_freqs.append(linelist[j])
                unscrambled_eus.append(eus[j]/u.K)
                unscrambled_degs.append(degeneracies[j])
                unscrambled_aijs.append((10**log10aijs[j])*u.Hz)
    return unscrambled_qns,unscrambled_freqs,unscrambled_eus,unscrambled_degs,unscrambled_aijs

beamlist=beamer(files)*u.sr
#print(beamlist)
fluxes=fluxvalues(round(targetpixcrd[1][0],round(targetpixcrd[1][1]),files)*u.Jy*u.km/u.s#/u.sr
#print(fluxes)
unscrambledqns,unscrambledfreqs,unscrambledeus,unscrambleddegs,unscrambledaijs=unscrambler(files,slicedqns,lines)
print(f'files: {files}')
print(f'unscrambledqns: {unscrambledqns}')
print(f'unscrambledfreqs: {unscrambledfreqs}')

'''Places Splatalogue table params, fluxes, and beams into dictionary'''
datadict={}
for i in range(len(fluxes)):
    datadict[i]={'qns':unscrambledqns[i],'freq':unscrambledfreqs[i],'beam':beamlist[i],'flux':fluxes[i],'E_u(K)':unscrambledeus[i],'degen':unscrambleddegs[i],'aij':unscrambledaijs[i]}

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
    
def Q_rot_asym(T):
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)
    
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
home2=home.replace('VelocityMoms3/','spw1ltemodelspec/')
plt.clf()
plt.scatter(unscrambledeus,np.log10(n_us))
plt.title(r'spw1 CH$_3$OH Rotational Diagram')
plt.xlabel(r'E$_u$ (K)')
plt.ylabel(r'log$_{10}$(N$_u$/g$_u$)')
plt.savefig((home2+'trialrotdiag2.png'),dpi=100,overwrite=True)
plt.show()


