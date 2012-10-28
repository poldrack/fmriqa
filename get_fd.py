#!/usr/bin/env python
"""
implement the frame displacement measure described by Power et al. (2011)
"""


import numpy as N
import os,sys

# load the datafile

if len(sys.argv)>1:
    subcode=sys.argv[1]
else:
    sys.exit('python get_fd.py <subcode>')

FDthresh=0.5
basedir='/scratch/01329/poldrack/connectome-genome/'
bolddir=os.path.join(basedir,subcode,'epi')

mcfile=os.path.join(bolddir,'%s_epi_mcf.par'%subcode)
motpars=N.loadtxt(mcfile)

# compute absolute displacement
dmotpars=N.zeros(motpars.shape)

dmotpars[1:,:]=N.abs(motpars[1:,:] - motpars[:-1,:])

# convert rotation to displacement on a 50 mm sphere
# mcflirt returns rotation in radians
# from Jonathan Power:
#The conversion is simple - you just want the length of an arc that a rotational
# displacement causes at some radius. Circumference is pi*diameter, and we used a 5
# 0 mm radius. Multiply that circumference by (degrees/360) or (radians/2*pi) to get the 
# length of the arc produced by a rotation.


headradius=50
disp=dmotpars.copy()
disp[:,0:3]=N.pi*headradius*2*(disp[:,0:3]/(2*N.pi))

FD=N.sum(disp,1)

N.savetxt(os.path.join(bolddir,'FD.txt'),FD)
N.savetxt(os.path.join(bolddir,'FDsquared.txt'),FD**2)
badvols=N.where(FD>FDthresh)[0]+1
if len(badvols)>0:
    N.savetxt(os.path.join(bolddir,'scrubvols.txt'),badvols,fmt='%d')
