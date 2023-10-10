from utilities import *
from astropy.table import Table, vstack, hstack
import glob
import numpy as np
import astropy.units as u
import astropy.constants as cnst
import os

c=cnst.c

'''This script creates the transition name table that's used in makemultipanelmom0.py. It only uses one source (SgrB2S) because the frequencies in this table need to be unredshifted back to their rest values.'''

source = 'SgrB2S'
fnum=fields[source]
sourcepath=sourcepath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcedict[source]#sourcedict[source]

sgrb2sz=dopplershifts[source]

txtfile=sourcepath+'mastereuksqnsfreqsdegens.txt'
array=np.genfromtxt(txtfile,dtype=str)
transitions=array[:,1]
#copy_transitions=np.copy(transitions)
freqs=(array[:,2].astype(np.float64)*u.Hz).to('GHz')
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
	'''
    if '$15_{6-}-16_{5-}$ $v_t$=1' in trans:
		trans=trans.replace('-}','\\pm}')
	'''
	newtransitions.append(trans)


stack=hstack([newtransitions,transitions,freqs,euppers])
table=Table(stack,names=['Transition','OldTransition','Frequency','$E_U$'],units=['','','(GHz)','(K)'])
'''#Used to remove the second entry of the line doublet and replace it with \\pm, but I think that's wrong/not what we want to do here
for row in table:
	if '$15_{6+}-16_{5+}$ $v_t$=1' in row:
		table.remove_row(int(np.where(table['Transition']=='$15_{6+}-16_{5+}$ $v_t$=1')[0]))
'''
#table.write('methanoltransitiontable.fits',overwrite=True)
datatblpath=datadir+'multiplot_methanoltransitiontable.fits'
table.write(datatblpath,overwrite=True)
print(f'Transition table saved at {datatblpath}')

print(transitions)
textable=Table(table)
textable.remove_column('OldTransition')
textblpath=datadir+'methanoltransitiontable.tex'
textable.write(textblpath,overwrite=True)
#os.system(f'cat {textblpath}')
print(f'Tex\'d version of table saved at {textblpath}')

