import numpy as np
from spectral_cube import SpectralCube as sc
import glob
import os

source='DS10'
homedict={'SgrB2S':"/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/SgrB2S/OctReimage_K/",'DSi':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSi/field10originals_K/",'DSii':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSii/field10originals_K/",'DSiii':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiii/field10originals_K/",'DSiv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSiv/field10originals_K/",'DSv':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSv/field10originals_K/",'DSVI':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSVI/field2originals_K/",'DSVII':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSVII/field3originals_K/",'DSVIII':"/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/DSVIII/field3originals_K/",'DSIX':"/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/DSIX/field7originals_K/",'DS10':
"/orange/adamginsburg/sgrb2/2017.1.00114.S/desmond/SgrB2DSminicubes/DS10/OctReimage_K/",'DSX':"/blue/adamginsburg/d.jeff/SgrB2DSminicubes/DSX/field7originals_K/"}
print(f'Source: {source}')

versiondict={'SgrB2S':'OctReimage_K','DSi':'field10originals_K','DSii':'field10originals_K','DSiii':'field10originals_K','DSiv':'field10originals_K','DSv':'field10originals_K','DSVI':'field2originals_K','DSVII':'field3originals_K','DSVIII':'field3originals_K','DS10':'OctReimage_K','DSIX':'field7originals_K','DSX':'field7originals_K'}
version=versiondict[source]
print(f'Data version: {version}')

xclasshome='/blue/adamginsburg/d.jeff/XCLASS2021/files/'

pixdict={'SgrB2S':(70,59),'DSi':(36,42),'DSii':(22,24),'DSiii':(24,24),'DSiv':(32,31),'DSv':(19,19),'DSVI':(62,62),'DSVII':(75,75),'DSVIII':(50,50),'DSIX':(34,35),'DS10':(35,35),'DSX':(60,60)}#SgrB2S uses the sample pixel, not its hotspot
pix=pixdict[source]
print(f'Sample pixel: {pix}')

incubes=glob.glob(homedict[source]+'*.fits')

images=['spw0','spw1','spw2','spw3']

datacubes=[]

for spew in images:
    for f1 in incubes:
        if spew in f1:
            datacubes.append(f1)
            continue
    
assert 'spw0' in datacubes[0], 'Cube list out of order'

outpath=xclasshome+f'{source}/{version}/'

if os.path.exists(outpath):
    print(f'Outpath {outpath} already exists.')
    pass
else:
    print(f'Creating outpath {outpath}')
    os.makedirs(outpath)

for cubepath,spw in zip(datacubes,images):
    print(f'Evaluating {cubepath}')
    cube=sc.read(cubepath)
    cubespec=cube[:,pix[0],pix[1]]
    if cubespec[0] > cubespec[1]:
        cubespec=cubespec[::-1]
    cubespecMHz=cubespec.spectral_axis.to('MHz')
    cubespecstack=np.stack((cubespecMHz.value,cubespec.value),axis=1)
    
    outfilepath=outpath+f'{spw}_{pix[0]}_{pix[1]}.txt'
    print('Saving')
    np.savetxt(outfilepath,cubespecstack)
    
print('Done')
