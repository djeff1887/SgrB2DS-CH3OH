from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import os
import pdb
import matplotlib as mpl

mpl.interactive(True)
plt.close('all')

comptable=Table.read('bootstrap_t150_compositetransposedensitytable.fits')#compositetable_transpose_preservesnrsforsgrb2s_wip.fits')

masses=list(comptable[4])[1:]
errormass=np.array(list(comptable[4])[1:])*0.56#list(comptable[5])[1:]
lums=list(comptable[6])[1:]
errorlum=np.array(list(comptable[6])[1:])*0.56#list(comptable[7])[1:]
temps=list(comptable[0])[1:]
errortemps=list(comptable[1])[1:]
abuns=np.array(list(comptable[16]))[1:]
abuns=list(map(float,abuns))
errorabun=np.array(list(comptable[17]))[1:]
errorabun=list(map(float,errorabun))
nh2s=list(comptable[2])[1:]

plt.figure()
plt.errorbar(temps,masses,yerr=errormass,xerr=errortemps,linestyle='')
plt.scatter(temps,masses,c=abuns,cmap='viridis',norm=mpl.colors.LogNorm(),zorder=5)
plt.yscale('log')
plt.xlabel('T$_{peak}$ (K)',fontsize=14)
plt.ylabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.colorbar(pad=0,label='X(CH$_3$OH)$_{peak}$')#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.scatter(temps,abuns,c=nh2s,cmap='bone',norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xlabel('T$_{K}$ (K)',fontsize=14)
plt.ylabel('X(CH$_3$OH)',fontsize=14)
plt.colorbar(pad=0,label='N(H$_2}$) (cm$^{-2}$)')#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.scatter(masses,lums,c=temps,cmap='inferno')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('L$_{core}$ (L$_{\odot}$)',fontsize=14)
plt.colorbar(pad=0,label=('T$_{K}$ (K)'))#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.scatter(temps,lums,c=abuns,cmap='viridis',norm=mpl.colors.LogNorm())
plt.yscale('log')
#plt.xscale('log')
plt.xlabel('T$_{K}$ (K)',fontsize=14)
plt.ylabel('L$_{core}$ (L$_{\odot}$)',fontsize=14)
plt.colorbar(pad=0,label=('X(CH$_3$OH)'))#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.errorbar(masses,abuns,xerr=errormass,yerr=errorabun,linestyle='')
plt.scatter(masses,abuns,c=temps,cmap='inferno',zorder=5)#,norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xscale('log')
plt.xlabel('M$_{core}$ (M$_\odot$)',fontsize=14)
plt.ylabel('X(CH$_3$OH)$_{peak}$',fontsize=14)
plt.colorbar(pad=0,label=('T$_{peak}$ (K)'))#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()

plt.figure()
plt.errorbar(lums,abuns,xerr=errorlum,yerr=errorabun,linestyle='')
plt.scatter(lums,abuns,c=temps,cmap='inferno',zorder=5)#,norm=mpl.colors.LogNorm())
plt.yscale('log')
plt.xscale('log')
plt.xlabel('L$_{core}$ (L$_{\odot}$)',fontsize=14)
plt.ylabel('X(CH$_3$OH)$_{peak}$',fontsize=14)
plt.colorbar(pad=0,label=('T$_{peak}$ (K)'))#'T$_K$ (K)')
#figsavepath=figpath+f'radialavgnh2s_r{r}px_rphys{int(pixtophysicalsize.value)}AU_smoothed.png'
#plt.savefig(figsavepath,bbox_inches='tight',overwrite=True)
plt.show()
