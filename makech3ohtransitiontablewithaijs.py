import pickle
from astropy.table import QTable
import numpy as np
from utilities import *
import pdb
from astropy.units.quantity import Quantity
import glob
from astropy.io import fits
import sys

sources=sourcedict.keys()
obj='ch3ohlinesdict.obj'

#excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','15_6-15_7E1vt1'],'DSii':'','DSiii':'','DSiv':'','DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1'],'DSVII':'','DSVIII':'','DSIX':'','DSX':''}

pixdict={'SgrB2S':(73,54),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35)}#SgrB2S-73,54#69,58

exceptions=['DSv','DSiii']

for s in sources:
    fnum=fields[s]
    sourcehome=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{s}'
    sourcepath=sourcehome+sourcedict[s]
    objectpath=sourcepath+obj
    predataset=open(objectpath,'rb')
    dataset=pickle.load(predataset)
    nopeakpix=True
    #pdb.set_trace()
    spws=dataset.keys()
    excludedqns=excludedlines[s]
    masterfreqs=[]
    mastereuks=[]
    mastereujs=[]
    masterdegens=[]
    masteraijs=[]
    masterflux=[]
    mastererrflux=[]
    masterqns=[]
    masterfwhm=[]
    masterntot=[]
    for spw in spws:
        params=dataset[spw]
        lineqns=params.keys()
        for line in lineqns:
            if nopeakpix:
                if s in exceptions:
                    peakfluxpix=pixdict[s]
                    peakpixy=peakfluxpix[0]
                    peakpixx=peakfluxpix[1]
                    peakpix=dataset['spw1']['4_2-3_1E1vt0']['flux'][peakpixy,peakpixx]
                    nopeakpix=False
                else:
                    peakpix=np.nanmax(dataset['spw1']['4_2-3_1E1vt0']['flux'].value)
                    peakfluxpix=np.where(dataset['spw1']['4_2-3_1E1vt0']['flux'].value==peakpix)
                    peakpixy=int(peakfluxpix[0])
                    peakpixx=int(peakfluxpix[1])
                    nopeakpix=False
            else:
                pass
            #peakfluxpix=np.where(params[line]['flux']==pf)
            pf=params[line]['flux'][peakpixy,peakpixx]
            epf=params[line]['mom0err'][peakpixy,peakpixx]
            if pf <= 3*epf:
                continue
            elif np.isnan(pf):
                continue
            elif line in excludedqns:
                continue
            else:
                masterqns.append(line)
                euk=mastereuks.append(params[line]['euk'])
                euj=mastereujs.append(params[line]['eujs'])
                degen=masterdegens.append(params[line]['degen'])
                aij=masteraijs.append(params[line]['aij'])
                restfreq=masterfreqs.append(params[line]['freq'])#shift_freq is the observed frequency
                err_peakflux=mastererrflux.append(epf)
                if type(pf) != Quantity:
                    pf=Quantity(pf)
                peakflux=masterflux.append(pf)
                mom2s=glob.glob(sourcepath+'mom2/*fwhm.fits')
                ntot=np.squeeze(fits.getdata(sourcepath+'bootstrap_ntot_intstd_boostrap1000_nonegativeslope.fits'))#Selected the non-bootstrap image to not have to find a new peak pixel, the areas over the peaks are never masked anyway
                masterntot.append(ntot[peakfluxpix[0],peakfluxpix[1]]*u.cm**-2)
                for mom2 in mom2s:
                    if line in mom2:
                        mom2img=fits.getdata(mom2)
                        fwhm=mom2img[peakfluxpix[0],peakfluxpix[1]]*u.km/u.s
                        #pdb.set_trace()
                        masterfwhm.append(fwhm)
                    else:
                        pass
        
        #pdb.set_trace()
    columns=['QNs','Freq','EU(K)','EU(J)','g','Aij','Peak Velocity Integrated Flux','Peak Flux Error','Line Width','Ntot']
    sourcelist=[masterqns,masterfreqs,mastereuks,mastereujs,masterdegens,masteraijs,masterflux,mastererrflux,masterfwhm,masterntot]
    sourcetable=QTable(sourcelist,names=columns)
    '''
    if s != 'SgrB2S':
        sys.exit()
    else:
        pass
    '''
    paramtablepath=datadir+f'OpticalDepthTables/{s}_ntot_4-3peak.fits'#{s}_contpeak_nothiiregion.fits'
    if not os.path.exists(datadir+'OpticalDepthTables/'):
        print(f'Creating save directory at {datadir}OpticalDepthTables/')
        os.mkdir(datadir+'OpticalDepthTables/')
    #pdb.set_trace()
    print(f'Saving CH3OH paramtable at {paramtablepath}')
    sourcetable.write(paramtablepath,overwrite=True)
    #print(f'{s}\n')
    #print(sourcetable)
    #pdb.set_trace()

print('Done')
