#!/usr/bin/env python
"""
mk_slice_mosaic(imgdata,outfile,title='',ncols=6,colorbar=True)
 - given an image, make a mosaic image and save to a file
"""

import nibabel as nib
import numpy as N
import matplotlib.pyplot as plt
import sys


def mk_slice_mosaic(imgdata,outfile,title,contourdata=[],ncols=6,colorbar=True):
    if imgdata.shape[0]==imgdata.shape[1]==imgdata.shape[2]:
        min_dim=2
    else:
        min_dim=N.where(N.min(imgdata.shape[0:3])==imgdata.shape[0:3])[0]

        slice_dims=N.where(N.min(imgdata.shape[0:3])!=imgdata.shape[0:3])[0]
        print 'min_dim:',min_dim
        print 'slice_dim:',slice_dims
        
    #if not len(imgdata.shape)==3:
    #    sys.stdout.write(__doc__)
    if 1:
        nrows=int(N.ceil(N.float(imgdata.shape[min_dim])/ncols))
        mosaic=N.zeros((nrows*imgdata.shape[slice_dims[0]],ncols*imgdata.shape[slice_dims[1]]))
        if not contourdata==[]:
            contourmosaic=N.zeros((nrows*imgdata.shape[slice_dims[0]],ncols*imgdata.shape[slice_dims[1]]))
        ctr=0
        print mosaic.shape
        
        for row in range(nrows):
            rowstart=row*imgdata.shape[slice_dims[0]]
            rowend=(row+1)*imgdata.shape[slice_dims[0]]
            for col in range(ncols):
                if ctr<imgdata.shape[min_dim]:
                    colstart=col*imgdata.shape[slice_dims[1]]
                    colend=(col+1)*imgdata.shape[slice_dims[1]]
                    print rowstart,rowend,colstart,colend
                    if min_dim==0:
                        imgslice=imgdata[ctr,:,::-1].T
                        if not contourdata==[]:
                            contourslice=contourdata[ctr,:,::-1].T             
                    elif min_dim==1:
                        imgslice=imgdata[:,ctr,::-1].T
                        if not contourdata==[]:
                            contourslice=contourdata[:,ctr,::-1].T       
                    elif min_dim==2:
                        imgslice=imgdata[:,::-1,ctr].T
                        if not contourdata==[]:
                            contourslice=contourdata[:,::-1,ctr].T          
                    print imgslice.shape
                    mosaic[rowstart:rowend,colstart:colend]=imgslice
                    if not contourdata==[]:
                        contourmosaic[rowstart:rowend,colstart:colend]=contourslice
                    ctr+=1
    elif dir=='saggital':
        nrows=int(N.ceil(N.float(imgdata.shape[1])/ncols))
        mosaic=N.zeros((nrows*imgdata.shape[2],ncols*imgdata.shape[1]))
        if not contourdata==[]:
            contourmosaic=N.zeros((nrows*imgdata.shape[0],ncols*imgdata.shape[2]))
        ctr=0
        for row in range(nrows):
            rowstart=row*imgdata.shape[2]
            rowend=(row+1)*imgdata.shape[2]
            for col in range(ncols):
                if ctr<imgdata.shape[1]:
                    colstart=col*imgdata.shape[0]
                    colend=(col+1)*imgdata.shape[0]
                    mosaic[rowstart:rowend,colstart:colend]=N.flipud(imgdata[ctr,:,:].squeeze().T)
                    if not contourdata==[]:
                        contourmosaic[rowstart:rowend,colstart:colend]=contourdata[ctr,::-1,:].T
                    ctr+=1

    plt.figure(figsize=(8,8))
    plt.imshow(mosaic,cmap=plt.cm.gray)
    if not title=='':
        plt.title(title)
    if colorbar:
        plt.colorbar()
    if not contourdata==[]:
        plt.contour(contourmosaic,colors='red')
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()
    
