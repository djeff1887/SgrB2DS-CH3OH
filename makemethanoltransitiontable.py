from utilities import *
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import astropy.units as u
import astropy.constants as cnst
import os

c=cnst.c

source = 'SgrB2S'
sourcepath=sourcedict[source]

sgrb2sz=0.000234806

mom0home=sourcepath+'mom0/*masked.fits'

txtfile=sourcepath+'mastereuksqnsfreqsdegens.txt'
array=np.genfromtxt(txtfile,dtype=str)
transitions=array[:,1]
#copy_transitions=np.copy(transitions)
freqs=(array[:,2].astype(np.float)*u.Hz).to('GHz')
freqs=[round((x.value*(1+sgrb2sz)),5) for x in freqs]
euppers=[int(float(y)) for y in array[:,0]]

'''String slicing'''
newtransitions=[]
for trans in transitions:
	trans=''.join(('$',trans))
	trans=trans.replace('(','_{')
	if 'vt=0' in trans:
		trans=trans.replace('vt=0','')
	trans=trans.replace('vt',' $v_t$')
	if ')--' in trans:
		trans=trans.replace(')--','-}-')
	if ')+-' in trans:
		trans=trans.replace(')+-','+}-')
	if ')+' in trans:
		trans=trans.replace(')+','+}$')
	if ')- ' in trans:
		trans=trans.replace(')-','-}')
	
	if 'E1' in trans:
		trans=trans.replace('E1',' E1')
	if 'E2' in trans:
		trans=trans.replace('E2',' E2')
	
	if ') ' in trans:
		trans=trans.replace(') ', '}$ ')
	if '} ' in trans:
		trans=trans.replace('} ', '}$ ')
	badints=[1,2,3,4,8,7,9]
	for x in badints:
		bad=f')-{x}'
		if bad in trans:#or ')-2' in trans or ')-8' in trans or ')-7' in trans:
			trans=trans.replace(')-','}-')
	if ')-' in trans:
		trans=trans.replace(')-','-}$')
	if '$15_{6-}-16_{5-}$ $v_t$=1' in trans:
		trans=trans.replace('-}','\\pm}')
	
	newtransitions.append(trans)


stack=hstack([newtransitions,freqs,euppers])
table=Table(stack,names=['Transition','Frequency','$E_U$'],units=['','(GHz)','(K)'])
for row in table:
	if '$15_{6+}-16_{5+}$ $v_t$=1' in row:
		table.remove_row(int(np.where(table['Transition']=='$15_{6+}-16_{5+}$ $v_t$=1')[0]))
#table.write('methanoltransitiontable.fits',overwrite=True)
table.write('methanoltransitiontable.tex',overwrite=True)
#mom0s=glob.glob(mom0home)

print(transitions)
#print(table)
os.system('cat methanoltransitiontable.tex')

