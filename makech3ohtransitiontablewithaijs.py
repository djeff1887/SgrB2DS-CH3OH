import pickle
from astropy.table import QTable
import numpy as np
from utilities import *
import pdb
from astropy.units.quantity import Quantity
import glob
from astropy.io import fits

sources=sourcedict.keys()
obj='ch3ohlinesdict.obj'

excludedlines={'SgrB2S':['7_6-7_7E1vt1','14_6-14_7E1vt1','11_6-11_7E1vt1'],'DSi':['11_6-11_7E1vt1','25_3-24_4E1vt0','14_6-14_7E1vt1','7_6-7_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','15_6-15_7E1vt1'],'DSii':'','DSiii':'','DSiv':'','DSv':'','DSVI':["6_1--7_2-vt1",'14_6-14_7E1vt1','10_6-10_7E1vt1','9_6-9_7E1vt1','11_6-11_7E1vt1','13_6-13_7E1vt1','12_6-12_7E1vt1','13_3--14_4-vt2','13_3+-14_4+vt2','7_6-7_7E1vt1','16_6-16_7E1vt1','8_6-8_7E1vt1'],'DSVII':'','DSVIII':'','DSIX':'','DSX':''}

for s in sources:
    sourcepath=sourcedict[s]
    objectpath=sourcepath+obj
    predataset=open(objectpath,'rb')
    dataset=pickle.load(predataset)
    
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
    for spw in spws:
        params=dataset[spw]
        lineqns=params.keys()
        for line in lineqns:
            pf=np.nanmax(params[line]['flux'])
            peakfluxpix=np.where(params[line]['flux']==pf)
            epf=params[line]['mom0err'][peakfluxpix[0],peakfluxpix[1]]
            if pf <= 3*epf:
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
                for mom2 in mom2s:
                    if line in mom2:
                        mom2img=fits.getdata(mom2)
                        fwhm=mom2img[peakfluxpix[0],peakfluxpix[1]]*u.km/u.s
                        #pdb.set_trace()
                        masterfwhm.append(fwhm)
                    else:
                        pass
        
        #pdb.set_trace()
    columns=['QNs','Freq','EU(K)','EU(J)','g','Aij','Peak Velocity Integrated Flux','Peak Flux Error','Line Width']
    sourcelist=[masterqns,masterfreqs,mastereuks,mastereujs,masterdegens,masteraijs,masterflux,mastererrflux,masterfwhm]
    sourcetable=QTable(sourcelist,names=columns)
    
    paramtablepath=f'OpticalDepthTables/{s}.fits'
    print(f'Saving CH3OH paramtable at {paramtablepath}')
    sourcetable.write(paramtablepath,overwrite=True)
    #print(f'{s}\n')
    #print(sourcetable)
    #pdb.set_trace()

print('Done')
