from astropy.table import Table,QTable,Column,vstack,hstack
import numpy as np
import astropy.units as u
from math import log10, floor
import pdb

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

def round_to_1(x):
    return round(x, -int(floor(log10(abs(x))))) 

sumtable=QTable.read('bootstrap_t150_compositetransposedensitytable.fits')
sumtable_transpose=sumtable#transpose_table(sumtable)
transtable=sumtable_transpose#Table.read('compositetable_transpose_preservesnrsforsgrb2s_wip.fits')
tableclip=transtable

#transtable.write('t150_compositetabletranspose.fits')
#pdb.set_trace()

trotstrlist=[]
massstrlist=[]
lumstrlist=[]
abunstrlist=[]
radiuslist=[]
a1list=[]
a2list=[]
dindexlist=[]
breaklist=[]

for key in tableclip.keys():
    if key == 'Core':
        trotstrlist.append('$T_{peak}$ [K]')
        radiuslist.append('$R_{150}$ (AU)')
        massstrlist.append('$M_{core}$ [M$_{\odot}$]')
        lumstrlist.append('$L$ [L$_{\odot}$]')
        abunstrlist.append('X(\methanol)$_{peak}$ ($10^{-8}$)')
        a1list.append('$\\alpha_1$')
        a2list.append('$\\alpha_2$')
        dindexlist.append('$p$')
        breaklist.append('$r_{break}$')
    if key != 'Core':
        masserr=str(int(tableclip[key][5]))
        troterr=str(int(tableclip[key][1]))
        lumerr=tableclip[key][7]
        raderr=str(int(tableclip[key][9]))
        a1err=tableclip[key][11]
        a2err=tableclip[key][13]
        breakerr=int(tableclip[key][15])
        derr=tableclip[key][19]
        templumerr=np.inf
        if lumerr < 1:
            templumerr=round_to_1(lumerr)
            if templumerr==1:
                lumerr=str(int(templumerr))
            else:
                lumerr=templumerr
        else:
            lumerr=str(int(lumerr))
            pass
        abunerr=tableclip[key][17]/1e-8
        tempabunerr=np.inf
        if abunerr < 1:
            tempabunerr=round_to_1(abunerr)
            abunerr=str(tempabunerr)
        else:
            abunerr=str(int(abunerr))
            pass
        mass=str(round(tableclip[key][4]))
        trot=str(round(tableclip[key][0]))
        lum=tableclip[key][6]
        if lum < 1:
            lum=str(round(lum,(len(str(lumerr))-2)))
        else:
            if templumerr >= 1:
                lum=str(round(lum))
            else:
                lum=str(round(lum,(len(str(lumerr))-2)))
                #lum=str(int(lum))
                pass
        abun=tableclip[key][16]/1e-8
        if tempabunerr < 1:
            abun=str(round(abun,(len(str(abunerr))-2)))
        else:
            abun=str(round(abun))
        rad=str(int(tableclip[key][8]))
        alpha1=str(tableclip[key][10])
        alpha2=str(tableclip[key][12])
        dindex=str(round(tableclip[key][18],(len(str(derr))-2)))
        rbreak=str(int(tableclip[key][14]))

        massstr=f'${mass}\pm{masserr}$'
        trotstr=f'${trot}\pm{troterr}$'
        lumstr=f'${lum}\pm{lumerr}$'
        abunstr=f'${abun}\pm{abunerr}$'
        radstr=f'${rad}\pm{raderr}$'
        a1str=f'${alpha1}\pm{a1err}$'
        a2str=f'${alpha2}\pm{a2err}$'
        dstr=f'${dindex}\pm{derr}$'
        breakstr=f'{rbreak}\pm{breakerr}'
        
        massstrlist.append(massstr)
        trotstrlist.append(trotstr)
        lumstrlist.append(lumstr)
        abunstrlist.append(abunstr)
        radiuslist.append(radstr)
        a1list.append(a1str)
        a2list.append(a2str)
        dindexlist.append(dstr)
        breaklist.append(breakstr)
finaltable=Table(rows=[trotstrlist,radiuslist,massstrlist,lumstrlist,abunstrlist,a1list,a2list,breaklist,dindexlist],names=tableclip.keys())

print(finaltable)

finaltable.write('bootstrap_t150_density_summarytable_latex.csv',overwrite=True)
