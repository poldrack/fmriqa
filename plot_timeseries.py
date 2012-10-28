import numpy as N
import matplotlib.pyplot as plt
import statsmodels.api


def plot_timeseries(data,title,outfile,markers=[],markername=[],
                    plottrend=False,ylims=[],plotline=[],xlabel='timepoints',
                    ylabel=[]):
    
    fig=plt.figure(figsize=[10,3])
    fig.subplots_adjust(bottom=0.15)
    plt.plot(data)
    datarange=N.abs(N.max(data)-N.min(data))
    ntp=len(data)
    if ylims==[]:
        axislims=[0,ntp+1,N.min(data) - datarange*0.1,N.max(data) + datarange*0.1]
    else:
        axislims=[0,ntp-1,ylims[0],ylims[1]]
        
    if plottrend:
        X=N.vstack((N.ones(ntp), N.arange(ntp)- N.mean(N.arange(ntp)),N.arange(ntp)**2-N.mean(N.arange(ntp)**2))).T
        model=statsmodels.api.OLS(data,X)
        results=model.fit()
    
    plt.axis(axislims)
    fig.hold('on')
    if len(markers)>0:
        for s in markers:
            lp,=plt.plot([s,s],axislims[2:4],linewidth=2,color='red')
        plt.legend([lp],[markername])
    if plottrend:
        plt.plot(N.dot(X,results.params),color='black')
    if plotline:
        plt.plot([0,ntp],[plotline,plotline])
        
    plt.title(title)
    plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()

    if plottrend:
        return results
    else:
        return []