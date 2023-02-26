from astropy.table import Table,QTable,Column,vstack,hstack
import numpy as np
import astropy.units as u
from math import log10, floor
import pdb
import os

def transpose_table(tab_before, id_col_name='Core'):#source: https://gist.github.com/PBarmby/61c6addda169568f4a0d
    '''Returns a copy of tab_before (an astropy.Table) with rows and columns interchanged
        id_col_name: name for optional ID column corresponding to
        the column names of tab_before'''
    # contents of the first column of the old table provide column names for the new table
    # TBD: check for duplicates in new_colnames & resolve
    new_colnames=tuple(tab_before[tab_before.colnames[0]])
    # remaining columns of old table are row IDs for new table 
    new_rownames=tab_before.colnames[1:]
    # make a new, empty table
    tab_after=QTable(names=new_colnames)
    # add the columns of the old table as rows of the new table
    for r in new_rownames:
        tab_after.add_row(tab_before[r])
    if id_col_name != '':
        # add the column headers of the old table as the id column of new table
        tab_after.add_column(Column(new_rownames, name=id_col_name),index=0)
    return tab_after

pretablepath='contsanitycheck_bootmasked_t180_compositehotcoresummarytable.fits'
'''
inputtablepath='bootstrap_t150_compositetabletranspose.fits'

if not os.path.exists(inputtablepath):
    print(f'Input table at {inputtablepath} not detected')
    pretable=Table.read(pretablepath)
    transpose_pretable=transpose_table(pretable)
    print(transpose_pretable)
    pdb.set_trace()
    print(f'Saving transposed input table at {inputtablepath}')
    transpose_pretable.write(inputtablepath)
'''
sumtable=Table.read(pretablepath)
densitytable=Table.read('contsanitycheck_densityslopes_bootmasked.fits')
#cores=np.array(['DS1','DS2','DS3','DS4','DS5','DS6','DS7','DS8','DS9','SgrB2S'])
#cores_t=cores.reshape([10,1])
#core_table=QTable(cores_t)
#newdensitytable=hstack([core_table,densitytable])
#transdensitytable=transpose_table(newdensitytable)
newstack=hstack([sumtable,densitytable])
newsumtable=Table(newstack)
#newsumtable['Core'][18]=densitytable.keys()[0]
#newsumtable['Core'][19]=densitytable.keys()[1]
#newsumtable.remove_column(' ')
print(newsumtable)
pdb.set_trace()

newsumtable.write('contsanitycheck_t180_compositedensitytable.fits')




