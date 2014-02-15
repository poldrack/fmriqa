#!/usr/bin/env python
"""
fMRI quality control
- adapted from fsld_raw.R and fBIRN QA tools

USAGE: fmriqa.py bold_mcf.nii.gz <TR>
"""

import ctypes, sys, os

flags = sys.getdlopenflags()
sys.setdlopenflags(flags|ctypes.RTLD_GLOBAL)

import numpy as N
import nibabel as nib
from compute_fd import *
from statsmodels.tsa.tsatools import detrend
import statsmodels.api
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn.cross_validation
from matplotlib.backends.backend_pdf import PdfPages
from mk_slice_mosaic import *
from matplotlib.mlab import psd
from plot_timeseries import plot_timeseries
from MAD import MAD
from mk_report import mk_report

sys.setdlopenflags(flags)

# thresholds for scrubbing and spike detection
FDthresh=0.5
DVARSthresh=0.5
AJKZ_thresh=25
# number of timepoints forward and back to scrub
nback=1
nforward=2


def error_and_exit(msg):
    print msg
    sys.stdout.write(__doc__)
    sys.exit(2)


    
def main():
    verbose=True
    
    if len(sys.argv)>2:
        infile=sys.argv[1]
        TR=float(sys.argv[2])
    else:
        error_and_exit('')
        #infile='/corral-repl/utexas/poldracklab/openfmri/shared2/ds105/sub001/BOLD/task001_run001/bold_mcf.nii.gz'
        #TR=2.5

    qadir=fmriqa(infile,TR,verbose=verbose)
    
def fmriqa(infile,TR,outdir=None,maskfile=None,motfile=None,verbose=False,plot_data=True):
    save_sfnr=True

    if os.path.dirname(infile)=='':
        basedir=os.getcwd()
        infile=os.path.join(basedir,infile)
    elif os.path.dirname(infile)=='.':
        basedir=os.getcwd()
        infile=os.path.join(basedir,infile.replace('./',''))
    else:
        basedir=os.path.dirname(infile)
        
    if outdir==None:
        outdir=basedir
        
    qadir=os.path.join(outdir,'QA')

    if not infile.find('mcf.nii.gz')>0:
        error_and_exit('infile must be of form XXX_mcf.nii.gz')
        
    if not os.path.exists(infile):
        error_and_exit('%s does not exist!'%infile)

    if maskfile==None:
        maskfile=infile.replace('mcf.nii','mcf_brain_mask.nii')
    
    if not os.path.exists(maskfile):
        error_and_exit('%s does not exist!'%maskfile)

    if motfile==None:
        motfile=infile.replace('mcf.nii.gz','mcf.par')
    if not os.path.exists(motfile):
        error_and_exit('%s does not exist!'%motfile)
    
    
    if not os.path.exists(qadir):
        os.mkdir(qadir)
    else:
        print 'QA dir already exists - overwriting!'
       
    
    if verbose:
        print 'infile:',infile
        print 'maskfile:',maskfile
        print 'motfile:',motfile
        print 'outdir:',outdir
        print 'computing image stats'
        
    img=nib.load(infile)
    imgdata=img.get_data()
    
    nslices=imgdata.shape[2]
    ntp=imgdata.shape[3]
    
    maskimg=nib.load(maskfile)
    maskdata=maskimg.get_data()
    maskvox=N.where(maskdata>0)
    nonmaskvox=N.where(maskdata==0)
    if verbose:
        print 'nmaskvox:',len(maskvox[0])
    
    # load motion parameters and compute FD and identify bad vols for
    # potential scrubbing (ala Power et al.)
    
    motpars=N.loadtxt(motfile)
    fd=compute_fd(motpars)
    N.savetxt(os.path.join(qadir,'fd.txt'),fd)
    
    voxmean=N.mean(imgdata,3)
    voxstd=N.std(imgdata,3)
    voxcv=voxstd/N.abs(voxmean)
    voxcv[N.isnan(voxcv)]=0
    voxcv[voxcv>1]=1
    
    # compute timepoint statistics
    
    
    maskmedian=N.zeros(imgdata.shape[3])
    maskmean=N.zeros(imgdata.shape[3])
    maskmad=N.zeros(imgdata.shape[3])
    maskcv=N.zeros(imgdata.shape[3])
    imgsnr=N.zeros(imgdata.shape[3])
    
    for t in range(imgdata.shape[3]):
        tmp=imgdata[:,:,:,t]
        tmp_brain=tmp[maskvox]
        tmp_nonbrain=tmp[nonmaskvox]
        maskmad[t]=MAD(tmp_brain)
        maskmedian[t]=N.median(tmp_brain)
        maskmean[t]=N.mean(tmp_brain)
        maskcv[t]=maskmad[t]/maskmedian[t]
        imgsnr[t]=maskmean[t]/N.std(tmp_nonbrain)
        
 
    # perform Greve et al./fBIRN spike detection
    #1. Remove mean and temporal trend from each voxel.
    #2. Compute temporal Z-score for each voxel.
    #3. Average the absolute Z-score (AAZ) within a each slice and time point separately.
    # This gives a matrix with number of rows equal to the number of slices (nSlices)
    # and number of columns equal to the number of time points (nFrames).
    #4. Compute new Z-scores using a jackknife across the slices (JKZ). For a given time point,
    # remove one of the slices, compute the average and standard deviation of the AAZ across
    # the remaining slices. Use these two numbers to compute a Z for the slice left out
    # (this is the JKZ). The final Spike Measure is the absolute value of the JKZ (AJKZ).
    # Repeat for all slices. This gives a new nSlices-by-nFrames matrix (see Figure 8).
    # This procedure tends to remove components that are common across slices and so rejects motion.
    if verbose:
        print 'computing spike stats'
    
    detrended_zscore=N.zeros(imgdata.shape)
    detrended_data=N.zeros(imgdata.shape)
    
    for i in range(len(maskvox[0])):
        tmp=imgdata[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]
        tmp_detrended=detrend(tmp)
        detrended_data[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]=tmp_detrended
        detrended_zscore[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]=(tmp_detrended - N.mean(tmp_detrended))/N.std(tmp_detrended)
    
    loo=sklearn.cross_validation.LeaveOneOut(nslices)
    AAZ=N.zeros((nslices,ntp))
        
    for s in range(nslices):
        for t in range(ntp):
            AAZ[s,t]=N.mean(N.abs(detrended_zscore[:,:,s,t]))
    
    JKZ=N.zeros((nslices,ntp))
    
    if verbose:
        print 'computing outliers'
        
    for train,test in loo:
        for tp in range(ntp):
            train_mean=N.mean(AAZ[train,tp])
            train_std=N.std(AAZ[train,tp])
            JKZ[test,tp]=(AAZ[test,tp] - train_mean)/train_std
     
    AJKZ=N.abs(JKZ)
    spikes=[]
    if N.max(AJKZ)>AJKZ_thresh:
        print 'Possible spike: Max AJKZ = %f'%N.max(AJKZ)
        spikes=N.where(N.max(AJKZ,0)>AJKZ_thresh)[0]
    if len(spikes)>0:
        N.savetxt(os.path.join(qadir,'spikes.txt'),spikes)
            
    voxmean_detrended=N.mean(detrended_data,3)
    voxstd_detrended=N.std(detrended_data,3)
    voxsfnr=voxmean/voxstd
    meansfnr=N.mean(voxsfnr[maskvox])
    
    # create plots

    #
    #imgdata_flat=imgdata.reshape(N.prod(imgdata.shape))
    #imgdata_nonzero=imgdata_flat[imgdata_flat>0.0]
    
    scaledmean=(maskmean - N.mean(maskmean))/N.std(maskmean)
    mean_running_diff=N.zeros(maskmad.shape)
    mean_running_diff=(maskmean[1:]-maskmean[:-1])/((maskmean[1:]+maskmean[:-1])/2.0)
    DVARS=N.zeros(fd.shape)
    DVARS[1:]=N.sqrt(mean_running_diff**2)*100.0
    N.savetxt(os.path.join(qadir,'dvars.txt'),DVARS)

    badvol_index_orig=N.where((fd>FDthresh)*(DVARS>DVARSthresh))[0]
    #print badvol_index_orig
    badvols=N.zeros(len(DVARS))
    badvols[badvol_index_orig]=1
    badvols_expanded=badvols.copy()
    for i in badvol_index_orig:
        if i>(nback-1):
            start=i-nback
        else:
            start=0
        if i<(len(badvols)-nforward):
            end=i+nforward+1
        else:
            end=len(badvols)
        #print i,start,end
        badvols_expanded[start:end]=1
    badvols_expanded_index=N.where(badvols_expanded>0)[0]
    #print badvols_expanded_index
    if len(badvols_expanded_index)>0:
        N.savetxt(os.path.join(qadir,'scrubvols.txt'),badvols_expanded_index,fmt='%d')

        # make scrubing design matrix - one colum per scrubbed timepoint
        scrubdes=N.zeros((len(DVARS),len(badvols_expanded_index)))
        for i in range(len(badvols_expanded_index)):
            scrubdes[badvols_expanded_index[i],i]=1
        N.savetxt(os.path.join(qadir,'scrubdes.txt'),scrubdes,fmt='%d')

    else:
        scrubdes=[]
        
    # save out complete confound file
    confound_mtx=N.zeros((len(DVARS),14))
    confound_mtx[:,0:6]=motpars
    confound_mtx[1:,6:12]=motpars[:-1,:]-motpars[1:,:] # derivs
    confound_mtx[:,12]=fd
    confound_mtx[:,13]=DVARS
    if not scrubdes==[]:
        confound_mtx=N.hstack((confound_mtx,scrubdes))

    N.savetxt(os.path.join(qadir,'confound.txt'),confound_mtx)
    
    #plot_timeseries(scaledmean,'Mean in-mask signal (Z-scored)',
    #            os.path.join(qadir,'scaledmaskmean.png'),spikes,'Potential spikes')

    datavars={'imgsnr':imgsnr,'meansfnr':meansfnr,'spikes':spikes,'badvols':badvols_expanded_index}

    if plot_data:
        print 'before plot'
        trend=plot_timeseries(maskmean,'Mean signal (unfiltered)',os.path.join(qadir,'maskmean.png'),
                        plottrend=True,ylabel='Mean MR signal')
        print 'after plot'
        datavars['trend']=trend
        plot_timeseries(maskmad,'Median absolute deviation (robust SD)',
                        os.path.join(qadir,'mad.png'),ylabel='MAD')

        plot_timeseries(DVARS,'DVARS (root mean squared signal derivative over brain mask)',
                        os.path.join(qadir,'DVARS.png'),plotline=0.5,ylabel='DVARS')

        plot_timeseries(fd,'Framewise displacement',os.path.join(qadir,'fd.png'),
                        badvols_expanded_index,'Timepoints to scrub (%d total)'%len(badvols),
                        plotline=0.5,ylims=[0,1],ylabel='FD')

        psd=matplotlib.mlab.psd(maskmean,NFFT=128,noverlap=96,Fs=1/TR)

        plt.clf()
        fig=plt.figure(figsize=[10,3])
        fig.subplots_adjust(bottom=0.15)
        plt.plot(psd[1][2:],N.log(psd[0][2:]))
        plt.title('Log power spectrum of mean signal across mask')
        plt.xlabel('frequency (secs)')
        plt.ylabel('log power')
        plt.savefig(os.path.join(qadir,'meanpsd.png'),bbox_inches='tight')
        plt.close()

        plt.clf()
        plt.imshow(AJKZ,vmin=0,vmax=AJKZ_thresh)
        plt.xlabel('timepoints')
        plt.ylabel('slices')
        plt.title('Spike measure (absolute jackknife Z)')
        plt.savefig(os.path.join(qadir,'spike.png'),bbox_inches='tight')
        plt.close()

        if img.shape[0]<img.shape[1] and img.shape[0]<img.shape[2]:
            orientation='saggital'
        else:
            orientation='axial'

        mk_slice_mosaic(voxmean,os.path.join(qadir,'voxmean.png'),'Image mean (with mask)',contourdata=maskdata)
        mk_slice_mosaic(voxcv,os.path.join(qadir,'voxcv.png'),'Image CV')
        mk_slice_mosaic(voxsfnr,os.path.join(qadir,'voxsfnr.png'),'Image SFNR')


        mk_report(infile,qadir,datavars)

 
    # def save_vars(infile,qadir,datavars):
    datafile=os.path.join(qadir,'qadata.csv')
    f=open(datafile,'w')
    f.write('SNR,%f\n'%N.mean(datavars['imgsnr']))
    f.write('SFNR,%f\n'%datavars['meansfnr'])
    #f.write('drift,%f\n'%datavars['trend'].params[1])
    f.write('nspikes,%d\n'%len(datavars['spikes']))
    f.write('nscrub,%d\n'%len(datavars['badvols']))
    f.close()

    if save_sfnr:
        sfnrimg=nib.Nifti1Image(voxsfnr,img.get_affine())
        sfnrimg.to_filename(os.path.join(qadir,'voxsfnr.nii.gz'))

    return qadir

if __name__=='__main__':
    main()
    
