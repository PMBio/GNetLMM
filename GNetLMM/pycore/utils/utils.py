import numpy as np

def computeLinearKernel(X,jitter=1e-3):
    """
    computing linear kernel K=XX^T

    input:
    X      :   genetic markers [NxF]
    jitter :   noise that is added to the diagonal to ensure that the matrix is not ill-conditioned (default: 1e-3)
    """
    N = X.shape[0]
    K = np.dot(X,X.T)
    K/= np.diag(K).mean()
    K+= jitter*np.eye(N)
    return K



def getSnpIndicesWithinGeneWindow(snp_chrom,snp_pos,gene_chrom,gene_start,gene_stop=None,window=0):
    """
    computes for a given gene, a boolean vector indicating which SNPs are close to the gene

    input:
    snp_chrom    :    chromosomal position [F]
    snp_pos      :    position [F]
    gene_chrom   :    gene chromosome (scalar)
    gene_start   :    start of the gene (scalar)
    gene_stop     :    end of the gene (scalar), if not given, gene_start is used
    window       :    window around gene (scalar)
    """
    
    if gene_stop is None: gene_stop = gene_start
    idx = snp_chrom==gene_chrom
    idx*= gene_start - window <= snp_pos
    idx*= snp_pos <= gene_stop + window
    return idx


def getCumPos(chrom,pos):
    """
    getCumulativePosition

    input
    chrom    :    chromosomal position
    pos      :    position
    """
    poscum = np.copy(pos)
    n_chroms = int(chrom.max())
    x = 0
 
    for chrom_i in range(1,n_chroms+1):
        if np.sum(chrom==chrom_i)==0: continue
        I = chrom==chrom_i
        poscum[I]+=x
        x=poscum[I].max()

    return poscum



def getChromBounds(chrom,posCum):
    """
    getChromBounds


    input
    chrom    :   chromosomal position
    posCum   :   cumulatative position
    """
    n_chroms = int(chrom.max())
    chrom_bounds = [0]
    for chrom_i in range(2,n_chroms+1):
        I1 = chrom==chrom_i
        I0 = chrom==chrom_i-1
        _cb = 0.5*(posCum[I0].max()+posCum[I1].min())
        chrom_bounds.append(_cb)

    chrom_bounds = np.array(chrom_bounds)
    return chrom_bounds
