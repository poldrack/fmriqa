"""
make report for QA using reportlab module
"""

from reportlab.pdfgen import canvas
import numpy as N
import time
import os

def mk_report(infile,qadir,datavars):
              #imgsnr,meansfnr,spikes,badvols):
    
    timestamp=time.strftime('%B %d, %Y: %H:%M:%S')
    
    report_header=[]
    report_header.append('QA Report: %s'%timestamp)
    report_header.append('infile: %s'%infile)
    report_header.append('Mean SNR: %f'%N.mean(datavars['imgsnr']))
    report_header.append('Mean SFNR (in brain mask): %f'%datavars['meansfnr'])
    report_header.append('Drift term: %f (p=%f)'%(datavars['trend'].params[1],datavars['trend'].pvalues[1]))
    report_header.append('# potential spikes: %d'%len(datavars['spikes']))
    report_header.append('# timepoints exceeding FD threshold: %d'%len(datavars['badvols']))
    
    c = canvas.Canvas(os.path.join(qadir,"QA_report.pdf"))
    yloc=820
    stepsize=16
    for line in report_header:
        c.drawString(10,yloc,line)
        yloc=yloc-stepsize
    
    timeseries_to_draw=['maskmean.png','meanpsd.png','mad.png','DVARS.png','fd.png']
    
    tsfiles=[os.path.join(qadir,t) for t in timeseries_to_draw]
    
    ts_img_size=[467,140]
    yloc=yloc-ts_img_size[1]
    
    for imfile in tsfiles:
        c.drawImage(imfile, 45,yloc,width=ts_img_size[0],height=ts_img_size[1])
        yloc=yloc-ts_img_size[1]
    
    c.showPage()
    
    yloc=650
    c.drawImage(os.path.join(qadir,'spike.png'),20,yloc,width=500,height=133)
    yloc=330
    images_to_draw=['voxmean.png','voxsfnr.png','voxcv.png']
    imfiles=[os.path.join(qadir,t) for t in images_to_draw]
    c.drawImage(imfiles[0],0,yloc,width=300,height=300)
    c.drawImage(imfiles[1],300,yloc,width=300,height=300)
    yloc=20
    c.drawImage(imfiles[2],0,yloc,width=325,height=325)
    
    c.save()