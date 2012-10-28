#!/usr/bin/env python
"""
create scrubbed dataset ala power et al
"""

import nibabel as nib
import numpy as N
import os,sys



# load the datafile

if len(sys.argv)>1:
    subcode=sys.argv[1]
else:
    subcode='EJ2108'
    #sys.exit('python get_fd.py <subcode>')


basedir='/scratch/01329/poldrack/connectome-genome/'
bolddir=os.path.join(basedir,subcode,'epi')

scrubfile=os.path.join(bolddir,'scrubvols.txt')

