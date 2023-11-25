from utilities import *
from astropy.table import Table
import sys
import glob
from spectral_cube import SpectralCube as sc
import pdb
import math
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle
from pyspeckit.spectrum.models import lte_molecule
from astropy.modeling import models
import matplotlib as mpl

mpl.interactive(True)

minicubedirs={1:'OctReimage_K/',10:'field10originals_K/',2:'field2originals_K/',3:'field3originals_K/',7:'field7originals_K/',8:'field8originals_K/'}

plotlabeltablepath=datadir+'multiplot_methanoltransitiontable.fits'

plotlabeltable=Table.read(plotlabeltablepath)

excludelineformatstrings={}

truesubplotwidth=10*u.km/u.s

reprojcontfits=fits.open(contpath)
reprojcont=reprojcontfits[0].data*u.Jy
reprojcontrestfreq=225*u.GHz#manual addition 11/9/2022, wiggle room w/i GHz
cntmbeam=radio_beam.Beam.from_fits_header(reprojcontfits[0].header)
reprojcont_K=reprojcont.to('K',cntmbeam.jtok_equiv(reprojcontrestfreq))

for source in list(sourcedict.keys()):
    print(f'Source:{source}')
    fnum=fields[source]
    pixcoords=pixdict[source]
    minicubedir=minicubedirs[fnum]
    z=dopplershifts[source]
    
    cubelocs=f'/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/{source}/{minicubedir}*.fits'
    reorgpath=f'/blue/adamginsburg/d.jeff/SgrB2DSreorg/field{fnum}/CH3OH/{source}/'+sourcedict[source]
    mastertxttablepath=reorgpath+'mastereuksqnsfreqsdegens.fits'
    fwhmpath=glob.glob(reorgpath+'*fwhm*')[0]
    nch3ohpath=reorgpath+'bootstrap_ntot_intstd_boostrap1000_nonegativeslope.fits'
    picklepath=reorgpath+'ch3ohlinesdict.obj'
    trotmappath=reorgpath+'bootstrap_texmap_3sigma_allspw_withnans_weighted.fits'
    savefigbase=f'/blue/adamginsburg/d.jeff/repos/CH3OHTemps/figures/{source}'
    savefighome=savefigbase+sourcedict[source]
    savefigpath=savefighome+'multilineplot.png'
    #pdb.set_trace()
    if os.path.exists(mastertxttablepath):
        print(f'Retrieving mastertxttable from {mastertxttablepath}')
        mastertxttable=Table.read(mastertxttablepath)
        print('Done\n')
    else:
        print('Converting txt and void construct to astropy table')
        mastertxtpath=reorgpath+'mastereuksqnsfreqsdegens.txt'
        mastertxtdata=np.genfromtxt(mastertxtpath,dtype=None)
        mastertxttable=Table(mastertxtdata,names=['EuK','OldTransition','RedshiftedFrequency','Degeneracy'],units=[u.K,'',u.Hz,''])
        print(f'Saving mastertxttable at {mastertxttablepath}')
        mastertxttable.write(mastertxttablepath,overwrite=True)
        print('Done\n')
    
    safetycopy_mastertxttable=Table(mastertxttable)
    orderedlabels=[]
    i=0
    for line in safetycopy_mastertxttable:
        for label in plotlabeltable:
            if line['OldTransition']==label['OldTransition']:
                if qn_replace(line['OldTransition']) in excludedlines[source]:
                    #print(mastertxttable)
                    #pdb.set_trace()
                    exclindex=np.where(mastertxttable['OldTransition']==line['OldTransition'])[0][0]
                    mastertxttable.remove_row(exclindex)
                    #print(mastertxttable)
                    #pdb.set_trace()
                else:
                    #print(line)
                    #pdb.set_trace()
                    orderedlabels.append(label['Transition'])
            elif line['OldTransition'] not in plotlabeltable['OldTransition'] and line['OldTransition'] in mastertxttable['OldTransition']:
                exclindex=np.where(mastertxttable['OldTransition']==line['OldTransition'])[0][0]
                mastertxttable.remove_row(exclindex)
    #if i != len(safetycopy_mastertxttable):
    #    pdb.set_trace()
    mastertxttable.add_column(orderedlabels,name='Transition')
    sourcefwhm=fits.getdata(fwhmpath)*u.km/u.s
    sourcench3oh=np.squeeze(fits.getdata(nch3ohpath))*u.cm**-2
    openfile=open(picklepath,'rb')
    pklch3oh=pickle.load(openfile)
    trotmap=fits.getdata(trotmappath)*u.K
    measTrot=trotmap[pixcoords[0],pixcoords[1]]
    measlinewidth=sourcefwhm[pixcoords[0],pixcoords[1]]
    measQrot=qrot(measTrot)
    #pdb.set_trace()
    incubes=glob.glob(cubelocs)
    images=['spw0','spw1','spw2','spw3']

    cubenames=[]
    datacubes=[]
    
    print('Ordering spectral windows from lowest to highest frequency')
    for spew in images:
        for f1 in incubes:
            if spew in f1:
                cubenames.append(f1)
                datacubes.append(sc.read(f1).with_spectral_unit('GHz'))
                continue

    assert 'spw0' in cubenames[0], 'Cube list out of order'
    
    print('Sorting mastertxttable by EuK')
    mastertxttable.sort('EuK')
    exclfmtqns=[]
    
    print('Assemble plot layout')
    numcols=5
    #numexclusions=len(excludedlines[source])
    numdetections=len(mastertxttable['OldTransition'])#len(mastertxttable['OldTransition'])-numexclusions
    numrows=math.ceil(numdetections/numcols)
    fig,ax=plt.subplots(numrows,numcols,sharey=True,figsize=(12,7))#,use_gridspec=True)
    plt.rcParams['figure.dpi'] = 300
    print(f'Number of rows: {numrows}')
    print(f'Number of columns: {numcols}')
    datatoplot=[]
    modelstoplot=[]
    
    print('Collecting data for plotting')
    for line in mastertxttable:
        #assert len(mastertxttable)==len(plotlabeltable), 'Transition table lengths do not match'
        qn=line['OldTransition']
        labelqn=line['Transition']
        reffreq=(line['RedshiftedFrequency']*u.Hz).to('GHz')
        exclfmtqn=qn_replace(line['OldTransition'])
        tempstr=exclfmtqn
        exclfmtqns.append(tempstr)
        linewasplotted=False
        for cube,spwname in zip(datacubes,images):
            freqs=cube.spectral_axis
            freqmin=freqs.min()
            freqmax=freqs.max()
            if reffreq >= freqmin and reffreq <= freqmax and tempstr not in excludedlines[source]:
                if source == 'SgrB2S':
                    restline=reffreq*(1+z)
                    z_sgrb2s=0.0002223
                    reffreq=restline/(1+z_sgrb2s)
                    #print('Collecting measured spectra from spectral slab')
                    subplotwidth=velocitytofreq(truesubplotwidth,reffreq)
                    lineslab=cube.spectral_slab((reffreq-subplotwidth),(reffreq+subplotwidth))
                    #vel_lineslab=lineslab.with_spectral_unit('km s-1',velocity_convention='radio',rest_value=reffreq)
                    lineflux=lineslab[:,pixcoords[0],pixcoords[1]]
                    linespecax=lineslab.spectral_axis
                    vel_linespecax=lineslab.with_spectral_unit('km s-1',velocity_convention='radio',rest_value=reffreq).spectral_axis#freqtovelocity(linespecax,reffreq)
                    datatoplot.append([labelqn,vel_linespecax,lineflux,reffreq])
                    
                    #print('Collecting values for line model')
                    modlinewidth=velocitytofreq(measlinewidth,reffreq)
                    lineprofilesigma=modlinewidth/(2*np.sqrt(2*np.log(2)))
                    vel_lineprofilesigma=measlinewidth/(2*np.sqrt(2*np.log(2)))
                    phi_nu=lineprofile(sigma=lineprofilesigma,nu_0=restline,nu=restline)
                    methntot=sourcench3oh[pixcoords[0],pixcoords[1]]
                    linedeg=line['Degeneracy']
                    lineeuk=line['EuK']*u.K
                    lineeuj=(lineeuk*k).to('J')
                    lineaij=pklch3oh[spwname][tempstr]['aij']
                    
                    #print('Calculate line model')
                    modnupper=nupper_estimated(methntot,linedeg,measQrot,lineeuj,measTrot).to('cm-2')
                    intertau=lte_molecule.line_tau(measTrot, methntot, measQrot, linedeg, restline, lineeuj, lineaij)
                    est_tau=(intertau*phi_nu).to('')
                    trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=measTrot).to('K')
                    modelline=models.Gaussian1D(mean=(0*u.km/u.s), stddev=vel_lineprofilesigma, amplitude=trad)
                    modelstoplot.append([labelqn,modelline])
                    
                else:
                    #print('Collecting measured spectra from spectral slab')
                    subplotwidth=velocitytofreq(truesubplotwidth,reffreq)
                    lineslab=cube.spectral_slab((reffreq-subplotwidth),(reffreq+subplotwidth))
                    #vel_lineslab=lineslab.with_spectral_unit('km s-1',velocity_convention='radio',rest_value=reffreq)
                    lineflux=lineslab[:,pixcoords[0],pixcoords[1]]
                    linespecax=lineslab.spectral_axis
                    vel_linespecax=lineslab.with_spectral_unit('km s-1',velocity_convention='radio',rest_value=reffreq).spectral_axis#freqtovelocity(linespecax,reffreq)
                    datatoplot.append([labelqn,vel_linespecax,lineflux,reffreq])

                    #print('Collecting values for line model')
                    restline=reffreq*(1+z)
                    modlinewidth=velocitytofreq(measlinewidth,reffreq)
                    lineprofilesigma=modlinewidth/(2*np.sqrt(2*np.log(2)))
                    vel_lineprofilesigma=measlinewidth/(2*np.sqrt(2*np.log(2)))
                    phi_nu=lineprofile(sigma=lineprofilesigma,nu_0=restline,nu=restline)
                    methntot=sourcench3oh[pixcoords[0],pixcoords[1]]
                    linedeg=line['Degeneracy']
                    lineeuk=line['EuK']*u.K
                    lineeuj=(lineeuk*k).to('J')
                    lineaij=pklch3oh[spwname][tempstr]['aij']

                    #print('Calculate line model')
                    modnupper=nupper_estimated(methntot,linedeg,measQrot,lineeuj,measTrot).to('cm-2')
                    intertau=lte_molecule.line_tau(measTrot, methntot, measQrot, linedeg, restline, lineeuj, lineaij)
                    est_tau=(intertau*phi_nu).to('')
                    trad=t_rad(tau_nu=est_tau,ff=f,nu=restline,T_ex=measTrot).to('K')
                    modelline=models.Gaussian1D(mean=(0*u.km/u.s), stddev=vel_lineprofilesigma, amplitude=trad)
                    modelstoplot.append([labelqn,modelline])
                    #sys.exit()
    #sys.exit()
    print('Beginning plotting')
    firstpanel=True
    i=0
    for row in np.arange(numrows):
        for col in np.arange(numcols):
            if i >= numdetections:
                ax[row,col].set_axis_off()
                #break
            else:
                assert datatoplot[i][0]==modelstoplot[i][0], 'Data line and model line do not match'
                if firstpanel:
                    print(f'Plotting {datatoplot[i][0]} at {datatoplot[i][3]}')
                    plotspecax=datatoplot[i][1]
                    im=ax[row,col].plot(plotspecax.value,datatoplot[i][2].value,color='black',drawstyle='steps-mid')
                    ax[row,col].plot(plotspecax.value,modelstoplot[i][1](plotspecax),color='red')
                    #fmax=maxtb#np.nanmax(transition)
                    #fmin=mintb#np.nanmin(transition)
                    #ax[row,col].tick_params(direction='in')
                    ax[row,col].set_xticks([-10,-5,0,5,10])
                    ax[row,col].set_xticklabels([-10,-5,0,5,10])
                    #ax[row,col].set_yticks([])
                    #ax[row,col].set_yticklabels([])
                    ax[row,col].set_title(fr'{datatoplot[i][0]}',fontsize=14)#annotate(f'{datatoplot[i][0]}',xy=(0,0),xytext=((max(datatoplot[i][2])-20*u.K).value,7),fontsize=10)#(datatoplot[i][3]-(subplotwidth/2))
                    #ax[row,col].axvline(x=datatoplot[i][3].value)
                    firstpanel=False
                    i+=1
                    #plt.show()
                    #sys.exit()
                else:
                    print(f'Plotting {datatoplot[i][0]} at {datatoplot[i][3]}')
                    plotspecax=datatoplot[i][1]
                    ax[row,col].plot(plotspecax.value,datatoplot[i][2].value,color='black',drawstyle='steps-mid')
                    ax[row,col].plot(plotspecax.value,modelstoplot[i][1](plotspecax),color='red')
                    #ax[row,col].tick_params(direction='in')
                    ax[row,col].set_xticks([-10,-5,0,5,10])
                    ax[row,col].set_xticklabels([-10,-5,0,5,10])
                    #ax[row,col].set_yticks([])
                    #ax[row,col].set_yticklabels([])
                    ax[row,col].set_title(fr'{datatoplot[i][0]}',fontsize=14)#ax[row,col].annotate(f'{datatoplot[i][0]}',xy=(0,0),xytext=((max(datatoplot[i][2])-20*u.K),(datatoplot[i][3]-(subplotwidth/2))),fontsize=10)#12
                    #ax[row,col].axvline(x=datatoplot[i][3].value)
                    i+=1
                    
    excludelineformatstrings.update({source:exclfmtqns})
    fig.text(0.5, 0.04, r'$v_{offset}$ (km s$^{-1}$)', ha='center',fontsize=14)
    fig.text(0.04, 0.5, r'$T_{b}$ (K)', va='center', rotation='vertical',fontsize=14)
    #fig.supylabel(r'$T_{b}$ (K)')
    #fig.supxlabel(r'$v_{rel}$ (km s$^{-1}$)')
    if source == 'SgrB2S':
        fig.subplots_adjust(hspace=0.55)
    else:
        fig.subplots_adjust(hspace=0.37)
    print(f'Saving figure at {savefigpath}')
    plt.savefig(savefigpath)
    plt.show()
    #sys.exit()