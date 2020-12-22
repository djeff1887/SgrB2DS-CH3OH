import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.ticker as mtick
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
import radio_beam
import regions
import math
import matplotlib as mpl

plt.close('all')
files=glob.glob('/blue/adamginsburg/d.jeff/imaging_results/*.fits')

mpl.interactive(True)

z=0.000186431
#z=0.00017594380066803095#DSii/iii
#z=0.000186431#DSi
#z=0.0002306756533745274#<<average of 2 components of 5_2-4_1 transition using old redshift(0.000236254)#0.000234806#0.000236254#0.0002333587#SgrB2S
c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1

testT=500*u.K
#n_totes=[1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20,1e21,1e22,1e23]*u.cm**-2
n_total=1e17*u.cm**-2

contaminants=[' CH3OCHO ',' HOONO ',' C3H6O2 ',' g-CH3CH2OH ',' HNCO ']
colors=cm.rainbow(np.linspace(0,1,len(contaminants)))

linelist='JPL'
linelistlist=['JPL','CDMS','SLAIM']

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in files:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

mdict={}
contamdata={}
'''
imgnames=['spw1','spw3','spw2','spw0']

assert imgnames[0] in files[0], 'Files out of order'
'''

def specmaker(plot,x,y,xmin,xmax,center,trans,ymax,ymin,moddata,thickmoddata):
    plot.set_xlim(xmin.value,xmax.value)
    plot.axvline(x=center,color='green',linestyle='--',linewidth=2.0,label='CH3OH')
    plot.set_ylim(ymin,ymax)
    plot.plot(x,y,drawstyle='steps')
    plot.plot(x,moddata,color='brown', label=(r'$\tau<<1$'))
    plot.plot(x,thickmoddata,color='cyan',label=(r'$\tau>1'))
    plot.set_title(trans)
    '''
    for mols in range(len(contaminants)):
        contamlabel=0
        linelistcheck=0
        for lis in linelistlist:
            if linelistcheck > 0:
                #print(contaminants[mols]+' already plotted.')
                break
            else:
                contamtable=Splatalogue.query_lines((mins[col+rowoffset]*(1+z)), (maxs[col+rowoffset]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[mols], line_lists=[lis],show_upper_degeneracy=True)
                if len(contamtable)==0:
                    print('No '+contaminants[mols]+' lines in '+lis+' frequency range '+str(mins[col+rowoffset])+'-'+str(maxs[col+rowoffset])+'.')
                    continue
                else:
                    linelistcheck+=1
                    print('('+lis+') '+contaminants[mols]+' contaminants identified for CH3OH '+mqns[col+rowoffset]+' at '+str(mins[col+rowoffset]+linewidth)+' GHz.')
                    table = utils.minimize_table(contamtable)
                    line=(table['Freq']*10**9)/(1+z)#Redshifted
                    qns=table['QNs']
                    for g in range(len(table)):
                        if g==0 and contamlabel==0:
                            contamline=ax[row,col].axvline(x=line[g],color=colors[mols],label=contaminants[mols])
                            print(contaminants[mols])
                            contamlabel+=1
                        else:
                            ax[row,col].axvline(x=line[g],color=colors[mols])
    '''
    '''
    print(f'ymin: {y.min()}')
    print(f'yvaluemin: {y.value.min()}')
    print(f'tempymin/ymin: {tempymin/y.value.min()}')
    '''

def gauss(x,A,mu,sig):
    return A*np.exp((-1/2)*((x-mu)/sig)**2)

def Q_rot_asym(T):#Eq 58, (Magnum & Shirley 2015); sigma=1, defined in Table 1 of M&S 2015
    return np.sqrt(m*np.pi*((k*T)/(h*b_0))**3)

def mulu(aij,nu):#Rearranged from Eq 11 (Magnum & Shirley 2015), returns product in units of cm5 g s-2
    return (3*h*c**3*aij)/(64*np.pi**4*nu**3)
    
def rjequivtemp(nu,T_ex):
    return ((h*nu)/k)/(np.exp((h*nu)/(k*T_ex))-1)
    
q=Q_rot_asym(testT).to('')
    
def vradio(frequency,rest_freq):
    velocity=c.to(u.km/u.s)*(1-((rest_freq-frequency)/rest_freq))
    return velocity.to('cm s-1')
    
def KtoJ(T):
    return (3/2)*k*T
    
def qngrabber(nums):
    temp=nums.split('(')
    temp2=temp[1].split(',')
    jupper=int(temp[0])
    if linelist == 'JPL':
        temp3=temp2[0].split(')')
        kupper=temp3[0]
        if 'a' in kupper:#What are these things?
            kupper=0
        else:
            kupper=int(temp3[0])
    else:
        kupper=int(temp2[0])
    
    return jupper, kupper

def Tb3(ntot,nu,line_width,mulu_2,s,g,q,eu_J,T_ex):#Rearranged from Eq 82, M&S 2015
    print(f'ntot: {ntot} nu: {nu} line_width: {line_width} mulu_2: {mulu_2} g: {g} q: {q} eu_J: {eu_J} Tex: {T_ex}')
    return ((8*np.pi**3*nu*mulu_2*R_i*g*f)/(3*k*q*np.exp(eu_J/(k*T_ex))*line_width))*ntot
    
def Tbthick(ntot,nu,line_width,mulu_2,g,q,eu_J,T_ex):
    return (1-np.exp(((-8*np.pi**3*mulu_2*R_i*g)/(3*h*q*line_width))*((np.exp((h*nu)/(k*T_ex))-1)/np.exp((eu_J)/(k*T_ex)))*ntot))*(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))
    
def opticaldepth(Tr,nu,T_ex):
    return -np.log(1-(Tr/(f*(rjequivtemp(nu,T_ex)-rjequivtemp(nu,Tbg)))))
    
def N_u(ntot,qrot,gu,eu_J,T_kin):#Rearranged from Eq 31, M&S 2015
    return ntot/((qrot/gu)*np.exp(eu_J/(k*T_kin)))

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

def contamlines(plot,contamlinelist):
    return
    
pixelcoords=[]
for i in range(len(files)):
    print('Getting ready - '+datacubes[i])
    cube=sc.read(datacubes[i],use_dask=True)
    header=fits.getheader(datacubes[i])
    
    cube_w=cube.wcs
    targetworldcrd=[[0,0,0],[266.8323912,-28.3954383,0]]#DSiv
    #[[0,0,0],[266.8332640,-28.3969259,0]]#DSii/iii
    #[266.8316149,-28.3972040,0]]#DSi
    #[[0,0,0],[2.66835339e+02, -2.83961660e+01, 0]]#SgrB2S
    targetpixcrd=cube_w.all_world2pix(targetworldcrd,1,ra_dec_order=True)
    pixelcoords.append(targetpixcrd[1])
    
    assert targetpixcrd[1,0] > 0, 'Negative pixel coords'
    
    cubebeams=(cube.beams.value)*u.sr/u.beam
    targetpixspec=cube[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]
    targetpixspec_K=JybeamtoK(cubebeams,targetpixspec)
    targetpixK_std=np.nanstd(targetpixspec_K)
    
    freqs=cube.spectral_axis
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    linewidth=0.00485*u.GHz#Half of original 0.0097GHz
    lw2=linewidth/8
            
    '''Generate methanol table for contaminant search''' 
    mtable1=Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelist], show_upper_degeneracy=True)
    methanol_table= utils.minimize_table(mtable1)
        
    mdict[i]={datacubes[i]:methanol_table}
    mlines=(methanol_table['Freq']*10**9)/(1+z)
    mqns=methanol_table['QNs']
    meuks=methanol_table['EU_K']*u.K
    mlog10aijs=np.array(methanol_table['log10_Aij'])
    maijs=10**mlog10aijs*u.s**-1
    meujs=[]
    for euk in meuks:
        meujs.append(KtoJ(euk))
    mdegs=mtable1['Upper State Degeneracy']
    
    mins=[]
    maxs=[]
    
    yoffset=0.5#K
    yoffset2=5#K
    
    opticaldepths={}
    opticaldepthlist=[]
    
    print('Setting figure and ax variables')
    numcols=5
    numrows=math.ceil(len(mlines)/numcols)
    fig,ax=plt.subplots(numrows,numcols,sharey=True)
    print('Number of rows: ', numrows)
    
    print('Gathering mlines and and plot widths')
    for line in mlines:
        centroid=line*u.Hz
        minfreq=centroid-linewidth
        maxfreq=centroid+linewidth
        mins.append(minfreq)
        maxs.append(maxfreq)
        
    print('Begin figure plot loops')
    rowoffset=0
    preymax=-100
    preymin=100
    for row in range(numrows):
        print('Start Row '+str(row)+'.')
        for col in range(numcols):
            if col+rowoffset >= len(mlines):
                #handles, labels = ax[row,col].get_legend_handles_labels()
                break
            f1,f2 = maxs[col+rowoffset],mins[col+rowoffset]
            if f1 > f2:
                f1,f2 = f2,f1
            sub=cube.spectral_slab(f1,f2)
            spw=sub[:,int(round(targetpixcrd[1][1])),int(round(targetpixcrd[1][0]))]
            beamlist=spw.beams
            beamlist=(beamlist.value)*u.sr/u.beam
            spwtbs=JybeamtoK(beamlist,spw)
            
            spwtbs_stddev=np.std(spwtbs)
            
            lw2vel=vradio(lw2,mlines[col+rowoffset]*u.Hz)
            J,K=qngrabber(mqns[col+rowoffset])
            s_j=(J**2-K**2)/(J*(2*J+1))#Eq 58, M&S 2015
            n_upper=N_u(n_total,q,mdegs[col+rowoffset],meujs[col+rowoffset],testT).to('cm-2')
            mulu2=(mulu(maijs[col+rowoffset],mlines[col+rowoffset]*u.Hz)).to('cm5 g s-2')#u.statC*u.cm.to('cm(3/2) g(1/2) s-1 cm')
            print(f'n_upper: {n_upper}')
            tbright=Tb3(n_total,mlines[col+rowoffset]*u.Hz,lw2vel,mulu2,s_j,mdegs[col+rowoffset],q,meujs[col+rowoffset],testT).to('K')#Tb2(mlines[col+rowoffset]*u.Hz,lw2vel,s_j,n_upper).to('K')
            tbthick=Tbthick(n_total,mlines[col+rowoffset]*u.Hz,lw2vel,mulu2,mdegs[col+rowoffset],q,meujs[col+rowoffset],testT).to('K')
            print(f'targetpixK_std: {targetpixK_std}')
            if tbthick.value >= targetpixK_std:
                tau=opticaldepth(tbthick,mlines[col+rowoffset]*u.Hz,testT)
                print(f'tau: {tau}')
                
                print(f'Tbthick: {tbthick}')
                modeltbs=[]
                thickmodeltbs=[]
                
                for hz in spw.spectral_axis:
                    modeltbs.append((gauss(hz,tbright,mlines[col+rowoffset]*u.Hz,lw2)/u.K))    
                    thickmodeltbs.append((gauss(hz,tbthick,mlines[col+rowoffset]*u.Hz,lw2)/u.K))
                tempymax=max(spwtbs)
                tempymin=min(spwtbs)

                if row*col > numrows*numcols:
                    break
                
                if tempymax > preymax:
                    reymax=tempymax#+yoffset2
                    #print('new max: ',reymax)
                else:
                    reymax=preymax

                if tempymin < preymin:
                    reymin=tempymin#-yoffset2
                    #print('new min: ',reymin)
                else:
                    reymin=preymin
                '''
                print(f'row: {row} col:{col}')
                print(f'tempymax: {tempymax} spw max: {spw.max().to("mJy/beam")}')
                print(f'tempymin: {tempymin} spw min: {spw.min().to("mJy/beam")}')
                print(f'reymax: {reymax} reymin: {reymin}')
                '''

                specmaker(ax[row,col],spw.spectral_axis,spwtbs,mins[col+rowoffset],maxs[col+rowoffset], mlines[col+rowoffset], mqns[col+rowoffset],reymax,reymin,modeltbs,thickmodeltbs)
                ax[row,col].xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
                preymax=reymax
                preymin=reymin
                '''
                for mols in range(len(contaminants)):
                    contamlabel=0
                    linelistcheck=0
                    for lis in linelistlist:
                        if linelistcheck > 0:
                            #print(contaminants[mols]+' already plotted.')
                            break
                        else:
                            contamtable=Splatalogue.query_lines((mins[col+rowoffset]*(1+z)), (maxs[col+rowoffset]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[mols], line_lists=[lis],show_upper_degeneracy=True)
                            if len(contamtable)==0:
                                print('No '+contaminants[mols]+' lines in '+lis+' frequency range '+str(mins[col+rowoffset])+'-'+str(maxs[col+rowoffset])+'.')
                                continue
                            else:
                                linelistcheck+=1
                                print('('+lis+') '+contaminants[mols]+' contaminants identified for CH3OH '+mqns[col+rowoffset]+' at '+str(mins[col+rowoffset]+linewidth)+' GHz.')
                                table = utils.minimize_table(contamtable)
                                line=(table['Freq']*10**9)/(1+z)#Redshifted
                                qns=table['QNs']
                                for g in range(len(table)):
                                    if g==0 and contamlabel==0:
                                        ax[row,col].axvline(x=line[g],color=colors[mols],label=contaminants[mols])
                                        contamlabel+=1
                                    else:
                                        ax[row,col].axvline(x=line[g],color=colors[mols])
                '''
            else:
                print('Line below 1sigma threshold')
                pass        
        rowoffset+=5
        
        '''
        if a == 0:
            ax.axvline(x=centroid.value,color='green',label='CH3OH')
        else:
            ax.axvline(x=centroid.value,color='green')
        ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
    print('Begin plotting contaminant lines')
    for j in range(len(contaminants)):
        print('Checking'+contaminants[j])
        dum=0
        for d in range(len(mins)):
            contamtable=Splatalogue.query_lines((mins[d]*(1+z)), (maxs[d]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[j], line_lists=[linelist],show_upper_degeneracy=True)
            if len(contamtable)==0:
                print('No '+contaminants[j]+' lines in frequency range '+str(mins[d])+'-'+str(maxs[d])+'.')
            else:
                print(contaminants[j]+' contaminants identified for CH3OH '+mqns[d]+' at '+str(mins[d]+linewidth)+' GHz.')
                table = utils.minimize_table(contamtable)
                line=(table['Freq']*10**9)/(1+z)#Redshifted
                qns=table['QNs']
                for g in range(len(table)):
                    if g==0 and dum==0:
                        ax.axvline(x=line[g],color=colors[j],label=contaminants[j])
                        print('hiii')
                        dum+=1
                    else:
                        ax.axvline(x=line[g],color=colors[j])
        '''
    fig.suptitle(f'{datacubes[i]}, Tkin: {testT}, N_total: {n_total}')
    plt.legend(loc=0,bbox_to_anchor=(1.7,2.12))
    fig.subplots_adjust(wspace=0.2,hspace=0.55)
print('Plotting complete. plt.show()')
plt.show()
'''    
    elif i >= 2:
        print('Setting figure and ax variables')
        numcols=5
        numrows=math.ceil(len(mlines)/numcols)
        fig,ax=plt.subplots(numrows,numcols,sharey=True)
        print('Number of rows: ', numrows)

        
'''
'''        
        plt.plot(freqs,spw.value,drawstyle='steps')
        plt.ylabel('Jy/beam')
        plt.xlabel('Frequency (Hz)')
        plt.title((imgnames[i]+' '+'Contaminant-labeled Spectra'))
        ax=plt.subplot(111)
'''
'''
        print('Gathering mlines and plot widths')
        for line in mlines:
            centroid=line*u.Hz
            minfreq=centroid-(linewidth*1.5)
            maxfreq=centroid+(linewidth*1.5)
            mins.append(minfreq)
            maxs.append(maxfreq)
            
        print('Begin figure plot loops')
        rowoffset=0
        preymax=-100
        preymin=100
        for row in range(numrows):
            print('Start Row '+str(row)+'.')
            for col in range(numcols):
                if col+rowoffset >= len(mlines):
                    #handles, labels = ax[row,col].get_legend_handles_labels()
                    break
                f1,f2 = maxs[col+rowoffset],mins[col+rowoffset]
                if f1 > f2:
                    f1,f2 = f2,f1
                sub=cube.spectral_slab(f1,f2)
                spw=sub[:,649,383]
                beamlist=spw.beams
                beamlist=(beamlist.value)*u.sr/u.beam
                spwtbs=JybeamtoK(beamlist,spw)
                J,K=qngrabber(mqns[col+rowoffset])
                s_j=(J**2-K**2)/(J*(2*J+1))#Eq 58, M&S 2015
                lw2vel=vradio(lw2,mlines[col+rowoffset]*u.Hz)
                n_upper=N_u(n_total,q,mdegs[col+rowoffset],meujs[col+rowoffset],testT).to('cm-2')
                print(f'n_upper: {n_upper}\n')
                print(f'aij: {maijs[col+rowoffset]} lines: {mlines[col+rowoffset]}')
                mulu2=(mulu(maijs[col+rowoffset],mlines[col+rowoffset]*u.Hz)).to('cm5 g s-2')#u.statC*u.cm.to('cm(3/2) g(1/2) s-1 cm')
                print(f'mulu2: {mulu2}')
                tbright=Tb3(n_total,mlines[col+rowoffset]*u.Hz,lw2vel,mulu2,s_j,mdegs[col+rowoffset],q,meujs[col+rowoffset],testT).to('K')#Tb2(mlines[col+rowoffset]*u.Hz,lw2vel,s_j,n_upper).to('K')
                tbthick=Tbthick(n_total,mlines[col+rowoffset]*u.Hz,lw2vel,mulu2,mdegs[col+rowoffset],q,meujs[col+rowoffset],testT).to('K')
                tau=opticaldepth(tbthick,mlines[col+rowoffset]*u.Hz,testT)
                print(f'tau: {tau}')
                
                modeltbs=[]
                thickmodeltbs=[]
                for hz in spw.spectral_axis:
                    modeltbs.append(gauss(hz,tbright,mlines[col+rowoffset]*u.Hz,lw2)/u.K)   
                    thickmodeltbs.append(gauss(hz,tbthick,mlines[col+rowoffset]*u.Hz,lw2)/u.K) 
                
                print(f'Tbmax: {max(modeltbs)*u.K}')
                tempymax=max(spwtbs)
                tempymin=min(spwtbs)
                
                if tempymax > preymax:
                    reymax=tempymax+yoffset
                    #print('new max: ',reymax)
                else:
                    reymax=preymax

                if tempymin < preymin:
                    reymin=tempymin-yoffset
                    #print('new min: ',reymin)
                else:
                    reymin=preymin

                print(f'row: {row} col:{col}')
'''
'''
                print(f'tempymax: {tempymax} spw max: {spw.max().to("mJy/beam")}')
                print(f'tempymin: {tempymin} spw min: {spw.min().to("mJy/beam")}')
                print(f'reymax: {reymax} reymin: {reymin}')
'''
'''

                specmaker(ax[row,col],spw.spectral_axis,spwtbs,mins[col+rowoffset],maxs[col+rowoffset], mlines[col+rowoffset], mqns[col+rowoffset],reymax,reymin,modeltbs,thickmodeltbs)
                preymax=reymax
                preymin=reymin
                
'''
'''
                if row == 0:
                    specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col],maxs[col],mlines[col],mqns[col])
                    continue
                if row == 1:
                    specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+5],maxs[col+5], mlines[col+5],mqns[col+5])
                    continue
                if row == 2:
                    if col >= int(len(mlines)/numrows):
                        break
                    else:
                        specmaker(ax[row,col],freqs,spw.to('mJy/beam'),mins[col+10],maxs[col+10], mlines[col+10], mqns[col+10])
                        continue
'''
'''
            rowoffset+=5
        fig.subplots_adjust(wspace=0.2,hspace=0.55) 
        fig.suptitle(f'{imgnames[i]}, Tkin: {testT}, N_total: {n_total}')
        labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        lines=[]
        fig.legend(labels, loc='upper right', bbox_to_anchor=(1.7,2.12))
        #plt.legend(loc=0,bbox_to_anchor=(1.7,2.12))
        print('Plotting complete. plt.show()')
        plt.show()
'''
'''
            if b == 0:
                ax.axvline(x=centroid.value,color='green',label='CH3OH')
            else:
                ax.axvline(x=centroid.value,color='green')
            if (freqs[0]-freqs[1])<0:
                ax.plot(freqs[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],spw.value[cube.closest_spectral_channel(minfreq):cube.closest_spectral_channel(maxfreq)],drawstyle='steps',color='orange')
            else:
                ax.plot(freqs[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],spw.value[cube.closest_spectral_channel(maxfreq):cube.closest_spectral_channel(minfreq)],drawstyle='steps',color='orange')
'''
'''
        print('Begin plotting contaminant lines')
        for k in range(len(contaminants)):
            print('Checking'+contaminants[k]+'...')
            dummy=0
            for c in range(len(mins)):
                contamtable=Splatalogue.query_lines((mins[c]*(1+z)), (maxs[c]*(1+z)),energy_max=1840, energy_type='eu_k', chemical_name=contaminants[k], line_lists=[linelist],show_upper_degeneracy=True)
                if len(contamtable)==0:
                    print('No '+contaminants[k]+' lines in frequency range '+str(mins[c])+'-'+str(maxs[c])+'.')
                    continue
                else:
                    print(contaminants[k]+' contaminants identified for CH3OH '+mqns[c]+' in frequency range '+str(mins[c])+'-'+str(maxs[c])+'.')
                    table = utils.minimize_table(contamtable)
                    line=(table['Freq']*10**9)/(1+z)#Redshifted
                    qns=table['QNs']
                    dummy+=1
                    for f in range(len(table)):
                        if f == 0 and dummy == 1:
                            ax.axvline(x=line[f],color=colors[k],label=contaminants[k])
                        else:
                            ax.axvline(x=line[f],color=colors[k])
        plt.legend()
        plt.show()
'''