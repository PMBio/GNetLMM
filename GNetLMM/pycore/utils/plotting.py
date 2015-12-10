import numpy as np
import matplotlib.pylab as plt
import scipy.cluster.hierarchy as sch
import pdb

import utils


def plot_transhits_per_snp(pv,alpha,snp_pos,snp_chrom, gene_start, gene_chrom, gene_stop=None, dist=0, 
                           color=None,ylim=None,label=None,fn=None):
    """
    manhattan plot where each peak says how many gene expression levels are significant for that SNP

    input
    pv         :    pv-association matrix [FxT]
    alpha      :    association threshold
    snp_pos    :    snp positions [F]
    snp_chrom  :    snp chrom [F]
    gene_start :    start gene positions [T]
    gene_chrom :    gene chrom [T]
    gene_stop  :    stop gene positions [T]
    dist       :    distance
    
    color      :    coloring of the line
    ylim       :    sets the upper y-limit of the plot
    label      :    name of the method
    fn         :    filename to save the figure
    """
    if gene_stop is None:
        gene_stop = gene_start

    F,T = pv.shape
    nHits = np.zeros(F)
    
    for t in range(T):
        cis_idx = utils.getSnpIndicesWithinGeneWindow(snp_chrom,snp_pos,gene_chrom[t],gene_start[t],gene_stop[t],window=dist)
        _pv = pv[~cis_idx,t]
        nHits[~cis_idx] += _pv<alpha
  

    fig = plt.figure(figsize=(10,2.5))
    plt.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111)

    # plotting
    snp_chrom = np.array(snp_chrom, dtype=int)
    posCum = utils.getCumPos(snp_chrom,snp_pos)
    plt.plot(posCum,nHits,color=color,label=label)

    # setting ticks, etc...
    chromBounds = utils.getChromBounds(snp_chrom,posCum)
    n_chroms = chromBounds.shape[0]
    for chrom_i in range(0,n_chroms-1,2):
        plt.fill_between(posCum,0,ylim,where=(posCum>chromBounds[chrom_i]) & (posCum<chromBounds[chrom_i+1]),facecolor='LightGray',linewidth=0,alpha=0.5)
        
    if ylim is not None: plt.ylim((0,ylim))
    xticks = np.array([chromBounds[i:i+2].mean() for i in range(chromBounds.shape[0]-1)])
    ax.set_xticks(xticks)
    plt.xticks(fontsize=6)
    plt.xlim(0,posCum.max())
    ax.set_xticklabels(np.arange(1,n_chroms+1))
 
    plt.xlabel('genetic position')
    plt.ylabel('number of trans eQTLs')
    if label is not None:
        plt.legend(frameon=False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    
    if fn is not None:
        plt.savefig(fn)
        plt.close()

        
def plot_corr(corr, fn=None):
    """
    plotting correlation matrix after hierarchically clustering the data

    input:
    corr   :   correlation matrix [NxN]
    fn     :   output filename
    """

    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(corr, method='complete')
    Z = sch.dendrogram(Y, orientation='right',link_color_func = lambda k: '#8A0808')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["left"].set_visible(False)
    ax1.spines["bottom"].set_visible(False)
    

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Z = sch.dendrogram(Y, orientation='right',link_color_func = lambda k: '#8A0808')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    
    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx = Z['leaves']
    im = plt.imshow(1-np.absolute(corr[idx][:,idx]), aspect='auto', origin='lower',vmin=0,vmax=2,cmap= plt.cm.RdGy)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    axmatrix.spines["right"].set_visible(False)
    axmatrix.spines["top"].set_visible(False)
    axmatrix.spines["left"].set_visible(False)
    axmatrix.spines["bottom"].set_visible(False)

    if fn is not None:
        plt.savefig(fn)
        plt.close()


def plotROCcurve(methods,TPR,FPR,xlim,ylim,ncol=2,fn=None):
    """
    plotting ROC curve (True Positive Rate vs. False Positive Rate)
    """
    fig = plt.figure(figsize=(4,4))
    fig.subplots_adjust(top=0.85,bottom=0.2,left=0.2)
    ax = fig.add_subplot(111)


    for method in methods:
        plt.plot(FPR[method],TPR[method],label=method,linewidth=2)

    plt.xlim((xlim))
    plt.ylim((ylim))
    plt.legend(loc='upper center',frameon=False,bbox_to_anchor=(0.5,1.25),ncol=ncol,prop={'size':10})
    plt.grid(True)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if fn is not None:
        plt.savefig(fn)
        plt.close()
