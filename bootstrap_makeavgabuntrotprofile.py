import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from matplotlib.pyplot import cm

fielddict={'SgrB2S':1,'DSi':10,'DSii':10,'DSiii':10,'DSiv':10,'DSv':10,'DSVI':2,'DSVII':3,'DSVIII':3,'DSIX':7}
romansources={'SgrB2S':'Sgr B2(S)','DSi':'DS1','DSii':'DS2','DSiii':'DS3','DSiv':'DS4','DSv':'DS5','DSVI':'DS6','DSVII':'DS7',
              'DSVIII':'DS8','DSIX':'DS9'}
sources=list(fielddict.keys())[:2]

plt.rcParams["figure.dpi"]=150
dsfmt='*'
rainbow=cm.rainbow(np.arange(0,10,len(sources)))
cmaps=['plasma','viridis','cividis','inferno','magma','Greys','Reds','Blues','Greens','PuRd']
sourcecolors=['red','cyan','orange','green','deepskyblue','black','rosybrown','darkviolet','dimgray','olivedrab']
ls=['solid', 'dashed', 'dashdot', 'dotted','solid', 'dashed', 'dashdot', 'dotted','solid','dashed']

for s in sources:
    abuns=np.loadtxt(f'intstd_{s}_abuns.txt')
    errabuns=np.loadtxt(f'intstd_{s}_errabuns.txt')
    tex=np.loadtxt(f'intstd_{s}_tex.txt')
    errtex=np.loadtxt(f'intstd_{s}_errtex.txt')
    #X,Y=np.meshgrid(tex,abuns
    ok=np.isfinite(tex)*np.isfinite(abuns)*np.isfinite(errabuns)
    abuns=np.array(abuns[ok])
    errabuns=np.array(errabuns[ok])
    tex=np.array(tex[ok])
    errtex=np.array(errtex[ok])
    
    texabun=[tex,abuns]
    sortedtexabun=texabun[texabun[:,0].argsort()]
    print(sortedtexabun[0:5])