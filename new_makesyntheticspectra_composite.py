import numpy as np
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
from astroquery.linelists.cdms import CDMS
import astropy.units as u
import matplotlib.pyplot as plt
from astroquery.splatalogue import utils, Splatalogue
import scipy.constants as cnst
from astropy.io import fits
import glob
from astropy.modeling import models#Fittable1DModel, Parameter, fitting
from utilities import *#Q_rot_asym,mulu,vradio,t_rad,nupper_estimated,opticaldepth,qngrabber
import matplotlib as mpl
import pdb
import sys

Splatalogue.QUERY_URL= 'https://splatalogue.online/c_export.php'

mpl.interactive(True)

plt.close('all')

def lineprofile(sigma,nu_0,nu):
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(nu-nu_0)**2/(2*sigma**2))

def cdms_get_molecule_name(my_molecule_name, **kwargs):
    basename = dict(CDMS.query_lines(min_frequency=1*u.GHz, max_frequency=500*u.GHz, molecule=my_molecule_name, parse_name_locally=True, get_query_payload=True, **kwargs))['Molecules']
    return " ".join(basename.split(" ")[1:])

'''Collect constants for N_tot and N_upper calculations'''

source='DS10'

c=cnst.c*u.m/u.s
k=cnst.k*u.J/u.K
h=cnst.h*u.J*u.s
sigma_sb=cnst.sigma*u.W/((u.m)**(2)*(u.K)**(4))
b_0=24679.98*u.MHz
a_0=127484*u.MHz
c_0=23769.70*u.MHz
m=b_0**2/(a_0*c_0)
Tbg=2.7355*u.K

trotdict={'SgrB2S':300*u.K,'DSi':300*u.K,'DSii':150*u.K,'DSiii':150*u.K,'DSiv':150*u.K,'DSv':100*u.K,'DSVI':300*u.K,'DSVII':200*u.K,'DSVIII':215*u.K,'DSIX':150*u.K,'DS10':150*u.K}

testT=trotdict[source]
#qrot_partfunc=partfunc(testT)#Q_rot_asym(testT).to('')

R_i=1
kappa=((2*b_0)-a_0-c_0)/(a_0-c_0)
f=1


dopplershifts={'SgrB2S':0.000228,'DSi':0.0001865,'DSii':0.000163,'DSiii':0.00017500261911843952,'DSiv':0.00018225233186845314,'DSv':0.0001838576164010067,'DSVI':0.0001661613132158407,'DSVII':0.00016320118280935546,'DSVIII':0.0001662062062062062,'DSIX':0.00015453732389175085,'DS10':0.00015794099431186572}#:0.000190713}/old doppler S: 0.0002306756533745274/0.00015954965399894244/0.00016236367659115043

s_othermol_dshift_v={' CH3CHO ':67.45330305*u.km/u.s,' C2H5OH ':67.45330305*u.km/u.s,' CH3OCHO ':67.45330305*u.km/u.s,' C(18)O ':69.551850256*u.km/u.s,' 13CH3OH ':67.5*u.km/u.s,' SO ':70.5*u.km/u.s}#' CH3OH ':68352.680424
ds2_othermol_dshift_v={' CH3OCHO ':49*u.km/u.s,' CH3CHO ':49*u.km/u.s,' C2H5OH ':49.3*u.km/u.s}#47831.782945392486 m / s
ds5_othermol_dshift_v={}
othermol_dopplershift={' CH3CHO ':0.000225,' C2H5OH ':0.000225,' CH3OCHO ':0.000225,' C(18)O ':0.000232}
ds9_othermol_dshift_v={}

sourceothers={'SgrB2S':s_othermol_dshift_v,'DSi':{},'DSii':ds2_othermol_dshift_v,'DSiii':{},'DSiv':{},'DSv':ds5_othermol_dshift_v,'DSVI':{},'DSVII':{},'DSVIII':{},'DSIX':ds9_othermol_dshift_v,'DS10':{}}
othermol_dshift_v=sourceothers[source]

z=dopplershifts[source]
z_vel=z*c

sourcelocs={'SgrB2S': r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/SgrB2S/OctReimage_K','DSi':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSi/field10originals_K','DSii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSii/field10originals_K','DSiii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiii/field10originals_K','DSiv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSiv/field10originals_K','DSv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSv/field10originals_K','DSVI':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVI/field2originals_K','DSVII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVII/field3originals_K','DSVIII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSVIII/field3originals_K','DSIX':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DSIX/field7originals_K','DS10':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Spectra/files/DS10/OctReimage_K'}

texlocs={'SgrB2S':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/SgrB2S/nov2022continuumsanitycheck_limitvt1lines_centeronlinepeak_repline20-20/','DSi':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS1/nov2022continuumsanitycheck/','DSii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS2/nov2022continuumsanitycheck/','DSiii':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DSiii/Kfield10originals_noexclusions/','DSiv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DSiv/Kfield10originals_noexclusions/','DSv':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS5/Kfield10originals_noexclusions_include4-3_150K_trial2/','DSVI':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS6/nov2022continuumsanitycheck/','DSVII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS7/Kfield3originals_trial1_noexclusions/','DSVIII':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS8/nov2022contniuumsanitycheck/','DSIX':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS9/Kfield7originals_150K_trial1_noexclusions/','DS10':r'C:/Users/desmond/Dropbox/Research/SgrB2DS/Sources/DS10/march2023discovery_5kmslw/'}

arabicswitch={'DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7','DSVIII':'DS8','DSIX':'DS9','DS10':'DS10'}

if source in arabicswitch.keys():
    texlocs[source]=texlocs[source].replace(source,arabicswitch[source])
    texmappath=texlocs[source]+'texmap_3sigma_allspw_withnans_weighted.fits'
else:
    texmappath=texlocs[source]+'texmap_3sigma_allspw_withnans_weighted.fits'

texmapdata=fits.getdata(texmappath)*u.K

pixdict={'SgrB2S':(70,59),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,32),'DS10':(35,35)}#SgrB2S:61,64, DSIX:(34,35)

targetpix=pixdict[source]

testT=texmapdata[targetpix[0],targetpix[1]]#350*u.K

#pdb.set_trace()
sourcepath=sourcelocs[source]

print(f'Collecting spectra from {sourcepath}')
incubes=glob.glob(sourcepath)
inspecs=glob.glob(sourcepath+'/*.txt')
instds=glob.glob(texlocs[source]+'errorimgs/std/*.fits')
images=['spw0','spw1','spw2','spw3']

spectra=[]

stds=[]

for spew in images:
    for f1 in inspecs:
        if spew in f1:
            if str(targetpix[0]) in f1 and str(targetpix[1]) in f1:
                spectra.append(f1)
                continue
            else:
                continue
    for f2 in instds:
        if spew in f2:
            stds.append(f2)
            continue
        else:
            continue
    
assert 'spw0' in spectra[0] and 'spw0' in stds[0], 'List out of order'

print('Spectra and stds are sequential order')

#imgnum=0
testline=0

linewidth=2.5*u.km/u.s#2.5 km/s is ideal for DSVI
print(f'Absolute model line width: {linewidth}\n')

specieslist=[' CH3OH ',' CH3OCHO ',' HOONO ',' HNCO ',' DCN ',' H2CO ',' C2H5OH ',' CH3CHO ',' CH3COOH ',' CH3NH2 ', ' CH3OCH3 ', ' HC3N ',' NH2CHO ', ' NH2CN ',' NH2D ',' SO2 ',' SO ',' t-HCOOH ',' a-H2CCHOH ',' s-H2CCHOH ',' H2CCNH ',' CH3CH2CHO ',' HCCCHO ',' SiS ',' CH2DOH ',' C(18)O ',' HDO ',' CH2CHCN ',' CH3CH2CN ',' c-H2COCH2 ', ' c-HCCCH ',' CCS ',' CH2NH ',"Ethylene Glycol",' cis-CH2OHCHO ','Acetone', ' CH3CN ',' CH2CHCN ',' CH3O13CHO ', ' SiO ', ' OCS ', 'N2D+', ' CH3C(15)N ',' CH3CCH ',' CH3SH ',' 13CS ', ' H2S ', ' SO ',' CH3(18)OH ', ' 13CH3OH ']

linelistdict={' CH3OH ':'JPL',' CH3OCHO ':'JPL',' CH3CHO ':'JPL',' C2H5OH ':'CDMS',' CH3OCH3 ':'JPL',' DCN ':'JPL',' OCS ':'CDMS',' 13CH3OH ':'CDMS',' H2CO ':'CDMS',' HC3N ':'CDMS',' C(18)O ':'CDMS',' 13CS ':'CDMS',' SO2 ':'CDMS',' NH2CHO ':'JPL',' HNCO ':'CDMS',' SO ':'CDMS', ' SiO ':'CDMS',' H2S ':'CDMS',' c-HCCCH ':'CDMS', 'HC3N v7=1':'CDMS',' H213CO ':'CDMS',' 13CH3CN ':'CDMS',' CH3COOH ':'CDMS',' t-HCOOH ':'CDMS',' CH3O13CHO ':'TopModel',' HNO3 ':'JPL','CH3O13CHO, vt = 0, 1':'CDMS',' NH2CN ':'JPL',' CH2CHCN ':'CDMS','CH3OCHO v=1':'JPL',' 18OCS ':'CDMS',' CH3NCO, vb = 0 ':'CDMS'}

jplnamelist={}
cdmsnamelist={' 13CH3OH ':'C-13-H3OH, vt=0,1',' C(18)O ':'CO-18',' 13CS ':'C-13-S, v=0,1',' NH2CHO ':'HC(O)NH2, v=0',' c-HCCCH ':'c-C3H2','HC3N v7=1':'HC3N, v7=1',' H213CO ':'H2C-13-O',' 13CH3CN ':'C-13-H3CN, v=0',' 18OCS ':'O-18-CS',}#' CH3NCO, vb=0 ':'CH3NCO, vb=0'}#' H2S ':'H2S'}#'CO-18'}#' HC3N ':'HC3N, v=0',' C2H5OH ':'C2H5OH,v=0',' CH3OCH3 ':'CH3OCH3, v=0',' OCS ':'060503 OCS, v=0',

sgrb2scolumns={' CH3OH ':1.7e18*u.cm**-2,' CH3OCHO ':3e16*u.cm**-2, ' CH3CHO ':1.5e16*u.cm**-2,' C2H5OH ':9e16*u.cm**-2,' CH3OCH3 ':9e15*u.cm**-2,' DCN ':3.5e16*u.cm**-2, ' OCS ':6e17*u.cm**-2,' 13CH3OH ':7e17*u.cm**-2,' H2CO ':7e17*u.cm**-2,' HC3N ':1e16*u.cm**-2, ' C(18)O ':1.3e19*u.cm**-2,' 13CS ':3e16*u.cm**-2,' SO2 ':2e17*u.cm**-2,' NH2CHO ':9e15*u.cm**-2,' HNCO ':3e17*u.cm**-2,' SO ':5e17*u.cm**-2,' SiO ':1e15*u.cm**-2,' H2S ':2e18*u.cm**-2,' c-HCCCH ':5e15*u.cm**-2, 'HC3N v7=1':5e15*u.cm**-2,' H213CO ':7e16*u.cm**-2,' 13CH3CN ':1e16*u.cm**-2,' CH2CHCN ':4e15*u.cm**-2,' 18OCS ':3e17*u.cm**-2,' CH3NCO, vb = 0 ':1e16*u.cm**-2,}#'CH3OCHO v=1':1e17*u.cm**-2}#' CH3O13CHO ':1e14*u.cm**-2,' H2CCO ':1e16*u.cm**-2,}#' H2CS ':1e18*u.cm**-2,' CH3(18)OH ':2.5e16*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' c-HCCCH ':2.5e15*u.cm**-2, 'Acetone':6e13*u.cm**-2,' CH3C(15)N ':3e13*u.cm**-2,' SiN ':2e15*u.cm**-2, ' CH3NH2 ':9e15*u.cm**-2,}#' HOONO ':5e15*u.cm**-2,' CH3COOH ':2e15*u.cm**-2,
#CDMS - ' CH3OH ':1.2e17*u.cm**-2,

dsicolumns={' CH3OH ':3e17*u.cm**-2,' CH3OCHO ':5e14*u.cm**-2,' CH3CHO ':8e14*u.cm**-2,' C2H5OH ':7e15*u.cm**-2,' CH3OCH3 ':6e14*u.cm**-2,' DCN ':7e14*u.cm**-2,' OCS ':9e16*u.cm**-2,' 13CH3OH ':1.5e16*u.cm**-2,' H2CO ':4e16*u.cm**-2,' HC3N ':9e14*u.cm**-2,' C(18)O ':1.5e19*u.cm**-2,' 13CS ':3e15*u.cm**-2,' SO2 ':2.5e15*u.cm**-2,' NH2CHO ':2e15*u.cm**-2,' HNCO ':3e16*u.cm**-2,' SO ':3e16*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e17*u.cm**-2,' c-HCCCH ':7e14*u.cm**-2, 'HC3N v7=1':5e15*u.cm**-2,' H213CO ':6e15*u.cm**-2,' 13CH3CN ':3e14*u.cm**-2,' CH2CHCN ':1e14*u.cm**-2,' 18OCS ':6e15*u.cm**-2,' CH3NCO, vb = 0 ':3e14*u.cm**-2,}#{' CH3OH ':1e17*u.cm**-2,' CH3OCHO ':3e14*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':0.7e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':2e14*u.cm**-2,' CH3COOH ':0.7e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':6e13*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':0.8e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':8e14*u.cm**-2, ' SO ':7e11*u.cm**-2, ' t-HCOOH ':5e14*u.cm**-2,' SiS ':5e13*u.cm**-2,' C(18)O ':5e17*u.cm**-2,' CH2DOH ':5e15*u.cm**-2, ' CH2NH ':1e13*u.cm**-2,"Ethylene Glycol":1e13*u.cm**-2,'Acetone':0.25e13*u.cm**-2, ' SiO ':7e13*u.cm**-2, ' OCS ':5.25e15*u.cm**-2, ' 13CS ':5e14*u.cm**-2, ' 13CH3OH ':2e15*u.cm**-2}

ds2columns={' CH3OH ':1.2e18*u.cm**-2,' CH3OCHO ':3e15*u.cm**-2,' CH3CHO ':3e15*u.cm**-2,' C2H5OH ':3e16*u.cm**-2,' CH3OCH3 ':3e15*u.cm**-2,' DCN ':5e15*u.cm**-2,' OCS ':3e17*u.cm**-2,' 13CH3OH ':1.5e17*u.cm**-2,' H2CO ':4e17*u.cm**-2,' HC3N ':9e15*u.cm**-2,' C(18)O ':6e19*u.cm**-2,' 13CS ':1e16*u.cm**-2,' SO2 ':2.5e16*u.cm**-2,' NH2CHO ':2e15*u.cm**-2,' HNCO ':5e16*u.cm**-2,' SO ':1.2e17*u.cm**-2, ' SiO ':5e15*u.cm**-2,' H2S ':1.2e18*u.cm**-2,' c-HCCCH ':3e15*u.cm**-2, 'HC3N v7=1':6e16*u.cm**-2,' H213CO ':2e16*u.cm**-2,' 13CH3CN ':2e15*u.cm**-2,' CH2CHCN ':5e14*u.cm**-2,' 18OCS ':6e15*u.cm**-2,' CH3NCO, vb = 0 ':3e14*u.cm**-2,}#'CH3O13CHO, vt = 0, 1':5e15*u.cm**-2' CH3(18)OH ':2.5e15*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' (13)CN ':1.5e15*u.cm**-2  NH2CN ':1e13*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,  ' SiS ':4e15*u.cm**-2,}' HOONO ':1e16*u.cm**-2,'Carbon Dioxide':5e16*u.cm**-2,' NaCl ':1e16*u.cm**-2,' CH3COOH ':7e14*u.cm**-2,' HDCO ':1e18*u.cm**-2,

ds5columns={' CH3OH ':7e15*u.cm**-2,' CH3OCHO ':5e14*u.cm**-2,' CH3CHO ':1e14*u.cm**-2,' C2H5OH ':1e15*u.cm**-2,' CH3OCH3 ':9e13*u.cm**-2,' DCN ':4e14*u.cm**-2,' OCS ':3e15*u.cm**-2,' 13CH3OH ':7e15*u.cm**-2,' H2CO ':9e15*u.cm**-2,' HC3N ':2e14*u.cm**-2,' C(18)O ':1.5e18*u.cm**-2,' 13CS ':8e14*u.cm**-2,' SO2 ':9e14*u.cm**-2,' NH2CHO ':1e14*u.cm**-2,' HNCO ':5e14*u.cm**-2,' SO ':3e15*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e16*u.cm**-2,' c-HCCCH ':9e13*u.cm**-2, ' HC3N v7=1':5e14*u.cm**-2,' H213CO ':7e14*u.cm**-2,' 13CH3CN ':6e13*u.cm**-2,}#{' CH3OH ':1e15*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1.75e15*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e13*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':0.3e14*u.cm**-2}

ds6columns={' CH3OH ':3e17*u.cm**-2,' CH3OCHO ':5e14*u.cm**-2,' CH3CHO ':8e14*u.cm**-2,' C2H5OH ':7e15*u.cm**-2,' CH3OCH3 ':6e14*u.cm**-2,' DCN ':7e14*u.cm**-2,' OCS ':9e16*u.cm**-2,' 13CH3OH ':1.5e16*u.cm**-2,' H2CO ':4e16*u.cm**-2,' HC3N ':9e14*u.cm**-2,' C(18)O ':1.5e19*u.cm**-2,' 13CS ':3e15*u.cm**-2,' SO2 ':2.5e15*u.cm**-2,' NH2CHO ':2e15*u.cm**-2,' HNCO ':3e16*u.cm**-2,' SO ':3e16*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e17*u.cm**-2,' c-HCCCH ':7e14*u.cm**-2, ' HC3N v7=1':5e15*u.cm**-2,' H213CO ':6e15*u.cm**-2,' 13CH3CN ':3e14*u.cm**-2,' CH2CHCN ':1e14*u.cm**-2,' 18OCS ':6e15*u.cm**-2,' CH3NCO, vb = 0 ':3e14*u.cm**-2,}#{' CH3OH ':1e17*u.cm**-2,' CH3OCHO ':5e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1.3e16*u.cm**-2,' DCN ':3e15*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':5e16*u.cm**-2, ' CH3CHO ':3e15*u.cm**-2,' CH3COOH ':5e15*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':1e13*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' SiO ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' 13CH3OH ':1e17*u.cm**-2}

ds7columns={' CH3OH ':3e17*u.cm**-2,' CH3OCHO ':5e14*u.cm**-2,' CH3CHO ':8e14*u.cm**-2,' C2H5OH ':7e15*u.cm**-2,' CH3OCH3 ':6e14*u.cm**-2,' DCN ':7e14*u.cm**-2,' OCS ':9e16*u.cm**-2,' 13CH3OH ':1.5e16*u.cm**-2,' H2CO ':4e16*u.cm**-2,' HC3N ':9e14*u.cm**-2,' C(18)O ':1.5e19*u.cm**-2,' 13CS ':3e15*u.cm**-2,' SO2 ':2.5e15*u.cm**-2,' NH2CHO ':2e15*u.cm**-2,' HNCO ':3e16*u.cm**-2,' SO ':3e16*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e17*u.cm**-2,' c-HCCCH ':7e14*u.cm**-2, ' HC3N v7=1':5e15*u.cm**-2,' H213CO ':6e15*u.cm**-2,' 13CH3CN ':3e14*u.cm**-2,' CH2CHCN ':1e14*u.cm**-2,' 18OCS ':6e15*u.cm**-2,' CH3NCO, vb = 0 ':3e14*u.cm**-2,}#{' CH3OH ':1e16*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':1e14*u.cm**-2}

ds9columns={' CH3OH ':2.5e17*u.cm**-2,' CH3OCHO ':9e14*u.cm**-2,' CH3CHO ':1e14*u.cm**-2,' C2H5OH ':7e15*u.cm**-2,' CH3OCH3 ':9e14*u.cm**-2,' DCN ':3e15*u.cm**-2,' OCS ':1e17*u.cm**-2,' 13CH3OH ':1.5e16*u.cm**-2,' H2CO ':9e15*u.cm**-2,' HC3N ':9e14*u.cm**-2,' C(18)O ':9e18*u.cm**-2,' 13CS ':4e15*u.cm**-2,' SO2 ':2.5e15*u.cm**-2,' NH2CHO ':3e14*u.cm**-2,' HNCO ':5e15*u.cm**-2,' SO ':2e16*u.cm**-2, ' SiO ':6e14*u.cm**-2,' H2S ':3.5e17*u.cm**-2,' c-HCCCH ':8e14*u.cm**-2, 'HC3N v7=1':5e15*u.cm**-2,' H213CO ':5e15*u.cm**-2,' 13CH3CN ':5e14*u.cm**-2,}#{' CH3OH ':1e15*u.cm**-2,' CH3OCHO ':1e15*u.cm**-2,' HOONO ':1e15*u.cm**-2,' HNCO ':1e15*u.cm**-2,' DCN ':1e14*u.cm**-2,' H2CO ':1e16*u.cm**-2,' C2H5OH ':3e15*u.cm**-2, ' CH3CHO ':1e14*u.cm**-2,' CH3COOH ':1e14*u.cm**-2, ' CH3NH2 ':1e15*u.cm**-2, ' CH3OCH3 ':1e14*u.cm**-2,' HC3N ':1e14*u.cm**-2, ' NH2CHO ':1e12*u.cm**-2,' NH2CN ':1e14*u.cm**-2, ' NH2D ':1e15*u.cm**-2, ' SO2 ':1e15*u.cm**-2, ' SO ':6e11*u.cm**-2, ' t-HCOOH ':1e14*u.cm**-2,' SiS ':1e14*u.cm**-2, ' OCS ':3e15*u.cm**-2, ' SiO ':1e14*u.cm**-2}

ds10columns={' CH3OH ':6e17*u.cm**-2,' CH3OCHO ':3e15*u.cm**-2,' CH3CHO ':9e14*u.cm**-2,' C2H5OH ':9e15*u.cm**-2,' CH3OCH3 ':2e15*u.cm**-2,' DCN ':3e15*u.cm**-2,' OCS ':1.5e17*u.cm**-2,' 13CH3OH ':1e17*u.cm**-2,' H2CO ':6e16*u.cm**-2,' HC3N ':5e15*u.cm**-2,' C(18)O ':3e19*u.cm**-2,' 13CS ':1e16*u.cm**-2,' SO2 ':6e15*u.cm**-2,' NH2CHO ':5e14*u.cm**-2,' HNCO ':5e15*u.cm**-2,' SO ':5e16*u.cm**-2, ' SiO ':4e15*u.cm**-2,' H2S ':4e17*u.cm**-2,' c-HCCCH ':9e14*u.cm**-2, 'HC3N v7=1':5e15*u.cm**-2,' H213CO ':1e16*u.cm**-2,' 13CH3CN ':1.5e15*u.cm**-2,}

weeds=[' CH3OCHO ', ' CH3CHO ']
cdmsproblemchildren=['OCS','13CS','C(18)O','HNCO','SO','HC3N','CH3NCO, vb=0']
problemchildren2=['CH3NCO, vb=0']

sourcecolumns={'SgrB2S':sgrb2scolumns,'DSi':{}, 'DSii':ds2columns,'DSiii':{},'DSiv':{},'DSv':ds5columns,'DSVI':{},'DSVII':{},'DSVIII':{},'DSIX':ds9columns,'DS10':ds10columns}

columndict=sourcecolumns[source]

plt.rcParams['figure.dpi'] = 150
plt.figure(1, figsize=(30,10))
molcolors=['red','cyan','orange','brown','deepskyblue','darkviolet','yellow','pink','darkviolet','darkkhaki','silver','blue','lime','magenta','grey','plum','fuchsia','darkcyan','magenta','deeppink','gold','palegreen','goldenrod','indigo']
spwmoldict={}
dummylist=[]
p1firstmolline={}#list(np.ones(len(columndict.keys())))
p2firstmolline={}
p3firstmolline={}
p4firstmolline={}
firstmolline=p1firstmolline
plotspecpad=0.005*u.GHz
n=1

for m in columndict.keys():
    p1firstmolline.update({m:1})
    p2firstmolline.update({m:1})
    p3firstmolline.update({m:1})
    p4firstmolline.update({m:1})

for spectrum, img, stdimage in zip(spectra,images,stds):
    
    print('Getting ready - '+img)
    plt.rcParams['figure.dpi'] = 150
    plt.figure(n, figsize=(30,10))
    n+=1
    if img == 'spw1':
        firstmolline=p2firstmolline
    if img == 'spw2':
        firstmolline=p3firstmolline
    if img == 'spw3':
        firstmolline=p4firstmolline
        '''
        plt.xlim(xmin=(p1minfreq-plotspecpad).value,xmax=(p1maxfreq+plotspecpad).value)
        plt.xlabel(r'$\nu$ (Hz)',fontsize=16)
        plt.ylabel('T$_b$ (K)',fontsize=16)
        plt.ylim(ymax=100)
        plt.tick_params(labelsize=13)
        plt.legend()
        plt.tight_layout()
        plt.show()
        '''
    
    spec=np.genfromtxt(spectrum)
    error=fits.getdata(stdimage)[targetpix[0],targetpix[1]]*u.K

    freqs=(spec[:,0]*u.MHz).to('GHz')#cube.spectral_axis
    data=spec[:,1]*u.K
    freqflip=False
    if freqs[0] > freqs[1]:
        freqs=freqs[::-1]
        data=data[::-1]
        freqflip=True
        print('Corrected decreasing frequency axis')
    else:
        pass
    
    freq_min=freqs[0]*(1+z)#215*u.GHz
    freq_max=freqs[(len(freqs)-1)]*(1+z)#235*u.GHz
    
    assert freq_max > freq_min, 'Decreasing frequency axis'
    
    print('Plotting model spectra')
    plt.plot(freqs.value,data.value,drawstyle='steps-mid',color='black')
    
    '''Generate methanol table for use during contaminant search'''
    Jfreqs, Jaij, Jdeg, JEU, qrot = get_molecular_parameters('CH3OH',
                                                         catalog='JPL',
                                                         fmin=freq_min,
                                                         fmax=freq_max)
    qrot_partfunc=qrot(testT)
    methanol_table=Splatalogue.query_lines(freq_min, freq_max, chemical_name=' CH3OH ', energy_max=1840, energy_type='eu_k', line_lists=[linelistdict[' CH3OH ']], show_upper_degeneracy=True)
    minmethtable=utils.minimize_table(methanol_table)
    mlines=((minmethtable['Freq']*10**9)/(1+z)*u.Hz).to('GHz')
    mqns=minmethtable['QNs']
    meuks=minmethtable['EU_K']*u.K
    meujs=[]
    for euk in meuks:
        meujs.append(KtoJ(euk))
    mdegs=methanol_table['Upper State Degeneracy']
    mlog10aijs=minmethtable['log10_Aij']
    maijs=10**mlog10aijs*u.s**-1
    
    '''Create background model for the methanol lines and other species'''
    baseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    #mbaseline=models.Linear1D(slope=(0*(u.K/u.Hz)),intercept=0*u.K)
    baseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    #mbaseline.bounding_box=(freqs[0],freqs[(len(freqs)-1)])
    #modelspec=baseline
    methmodelspec=baseline
    plot=np.linspace(freqs[0],freqs[(len(freqs)-1)],np.shape(spec)[0]).to('GHz')
    modeldict={}
    
    for molecule,hue,first in zip(list(columndict.keys())[1:], molcolors,list(firstmolline.keys())[1:]):
        '''Generate species table for contaminant search'''
        modelspec=baseline
        species_table= Splatalogue.query_lines(freq_min, freq_max, chemical_name=molecule, energy_max=1840, energy_type='eu_k', line_lists=[linelistdict[molecule]], show_upper_degeneracy=True)
        if len(species_table['Chemical Name']) == 0:
            print(f'No transitions for {molecule} in {img}. Continue')
            continue
        else:
            pass
        #a=cdms_get_molecule_name('CH3NCO, vb=0')
        ll=linelistdict[molecule]
        param_molname=molecule.replace(' ','')
        if ll == 'CDMS':
            #if molecule in cdmsnamelist.keys():
            if molecule in cdmsnamelist.keys():
                molname=cdmsnamelist[molecule]
            elif molecule == ' CH3NCO, vb = 0 ':
                molname=cdms_get_molecule_name('CH3NCO, vb=0')
            else:
                molname=cdms_get_molecule_name(param_molname)#cdmsnamelist[molecule]
            #else:
                #molname=param_molname
                
            if param_molname in cdmsproblemchildren:
                cCfreqs, cCaij, cCdeg, cCEU, c_qrot = get_molecular_parameters(molname,catalog='CDMS',
                                                                               fmin=freq_min, fmax=(freq_max+50*u.GHz),)
            else:
                cCfreqs, cCaij, cCdeg, cCEU, c_qrot = get_molecular_parameters(molname,
                                                                               catalog='CDMS',fmin=freq_min,fmax=freq_max,)
        if ll == 'JPL':
            if molecule in jplnamelist.keys():
                molname=cdmsnamelist[molecule]
            else:
                molname=param_molname
            cJfreqs, cJaij, cJdeg, cJEU, c_qrot = get_molecular_parameters(molname,catalog='JPL',
                                                         fmin=freq_min,
                                                         fmax=freq_max,)
        c_qrot_partfunc=c_qrot(testT)
        
        minchemtable=utils.minimize_table(species_table)
        if molecule in othermol_dshift_v.keys():
            otherz=othermol_dshift_v[molecule]/c
            clines=((minchemtable['Freq']*10**9)/(1+otherz)*u.Hz).to('GHz')
        else:
            clines=((minchemtable['Freq']*10**9)/(1+z)*u.Hz).to('GHz')
        
        cqns=minchemtable['QNs']
        ceuks=minchemtable['EU_K']*u.K
        ceujs=[]
        for euk in ceuks:
            ceujs.append(KtoJ(euk))
        cdegs=species_table['Upper State Degeneracy']
        clog10aijs=minchemtable['log10_Aij']
        caijs=10**clog10aijs*u.s**-1
        cntot=columndict[molecule]
        print(f'Begin model loops for {molecule}')
        
        linedetections=[]
        for line,deg,euj,aij,qn in zip(clines,cdegs,ceujs,caijs,cqns):
            #print(f'Transition: {qn} @ {line.to("GHz")}')
            if molecule in othermol_dshift_v.keys():
                restline=line*(1+otherz)
            else:
                restline=line*(1+z)
            est_nupper=nupper_estimated(cntot,deg,c_qrot_partfunc,euj,testT).to('cm-2')
            modlinewidth=velocitytofreq(linewidth,line)
            phi_nu=lineprofile(modlinewidth,restline,restline)
            intertau=lte_molecule.line_tau(testT, cntot, qrot_partfunc, deg, restline, euj, aij) #opticaldepth(aij,restline,testT,est_nupper,modlinewidth).to('')
            est_tau=(intertau*phi_nu).to('')
            #print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
            trad=t_rad(f,est_tau,restline,testT).to('K')
            if trad >= 3*error:
                #print(f'Estimated brightness: {"{:.3f}".format(trad)}')
                #modlinewidth=velocitytofreq(linewidth,line)
                #print(f'Model linewidth (Hz): {modlinewidth}')
                modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
                #modelgaus+=modelline
                modelspec+=modelline
                linedetections.append(True)
            else:
                #print(f'{qn} line brightness ({trad}) below 3sigma threshold ({3*error})')
                linedetections.append(False)
                continue
        if molecule == ' CH3CHO ':
            dummylist.append((freqs,modelspec(freqs)))#spwmoldict.update({img:(freqs,modelspec(freqs))})
        if firstmolline[first]:
            plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color=hue,label=molecule)
            firstmolline[first]=0
        else:
            plt.plot(freqs,modelspec(freqs),drawstyle='steps-mid',color=hue)
    if ' CH3OH ' in columndict:
        print('Begin CH3OH modeling\n')
        mdetections=[]
        for line,deg,euj,aij,qn in zip(mlines,mdegs,meujs,maijs,mqns):
            #print(f'Transition: {qn} @ {line.to("GHz")}')
            restline=line*(1+z)
            modlinewidth=velocitytofreq(linewidth,line)
            phi_nu=lineprofile(modlinewidth,restline,restline)
            
            methntot=columndict[' CH3OH ']
            est_nupper=nupper_estimated(methntot,deg,qrot_partfunc,euj,testT).to('cm-2')
            intertau=lte_molecule.line_tau(testT, methntot, qrot_partfunc, deg, restline, euj, aij)#opticaldepth(aij,restline,testT,est_nupper,originallinewidth).to('')
            est_tau=(intertau*phi_nu).to('')
            #print(f'Estimated tau: {"{:.3f}".format(est_tau)}')
            trad=t_rad(f,est_tau,restline,testT).to('K')
            if trad >= 3*error:
                #print(f'Estimated brightness: {"{:.3f}".format(trad)}')
                #print(f'Model linewidth (Hz): {modlinewidth}')
                modelline=models.Gaussian1D(mean=line, stddev=modlinewidth, amplitude=trad)
                #modelgaus+=modelline
                methmodelspec+=modelline
                mdetections.append(True)
            else:
                #print(f'{qn} line brightness ({"{:.3f}".format(trad)}) below 3sigma threshold ({3*error})')
                mdetections.append(False)
                continue

        #compositespec=modelspec+methmodelspec
        if firstmolline[' CH3OH ']:
            plt.plot(freqs,methmodelspec(freqs),drawstyle='steps-mid',linestyle='--',color='green',label=' CH3OH ')
            firstmolline[' CH3OH ']=0
            print('yay')
        else:
            plt.plot(freqs,methmodelspec(freqs),drawstyle='steps-mid',linestyle='--',color='green')
            print('yayy')
    '''
    print('Overplotting axvlines and transition annotations')
    for line,qn,detected in zip(clines,cqns,linedetections):
        if detected:
            plt.axvline(x=line.value,linestyle='--',color='yellow',ymin=0.25)
            plt.annotate(qn, (line.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    
    for mline,mqn,detected in zip(mlines,mqns,mdetections):
        if detected:
            plt.axvline(x=mline.value,linestyle='--',color='pink',ymin=0.25)
            plt.annotate(mqn, (mline.value, 0), (line.value-0.002,40),rotation=90)
        else:
            continue
    '''
    '''
    if img == 'spw3':
        p2maxfreq=max(freqs)
        plt.xlim(xmin=(p2minfreq-plotspecpad).value,xmax=(p2maxfreq+plotspecpad).value)
    '''
    plt.xlabel(r'$\nu$ (Hz)',fontsize=16)
    plt.ylabel('T$_b$ (K)',fontsize=16)
    plt.xlim(xmin=(min(freqs)-plotspecpad).value,xmax=(max(freqs)+plotspecpad).value)
    plt.ylim(ymax=100)
    plt.tick_params(labelsize=13)
    plt.tight_layout()
    plt.legend()
    #plt.savefig(fr'C:/Users/desmond/Desktop/CH3OHTemps/CompositeSpectra/bootstrap_{source}_{img}_qrotfix_compositespectra.png')
    plt.show()
