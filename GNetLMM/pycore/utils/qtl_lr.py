import limix.modules.qtl as qtl
import limix.modules.varianceDecomposition as varianceDecomposition
import pdb
import utils
import numpy as np
import scipy.stats as stats
import time

from GNetLMM.pycore.mtSet.gp.gp3kronSum import gp3kronSum
from GNetLMM.pycore.mtSet.mean.mean import mean
import GNetLMM.pycore.mtSet.covariance.covariance as covariance
import GNetLMM.pycore.mtSet.covariance.freeform as freeform
import GNetLMM.pycore.mtSet.optimize.optimize_bfgs as OPT

def test_lmm_lr(G, y, Z, Kbg, Covs=None):
    """
    low-rank lmm

    input:
    G   :   genotypes
    y   :   phenotype
    Z   :   features of low-rank matrix
    Kbg   :   background covariance matrix
    Covs :  fixed effect covariates
    """

    vd = varianceDecomposition.VarianceDecomposition(y)
    if Covs is not None:
        vd.addFixedEffect(Covs)
    vd.addRandomEffect(Kbg)
    Klr = utils.computeLinearKernel(Z)
    vd.addRandomEffect(Klr)
    vd.addRandomEffect(is_noise=True)
    vd.optimize()

    varComps = vd.getVarianceComps()[0]
    Ktotal = varComps[0]*Kbg + varComps[1]*Klr

    lm = qtl.test_lmm(G,y,covs=Covs,K=Ktotal)
    pv = lm.getPv()[0]
    beta = lm.getBetaSNP()[0]
    
    var_snps =  beta**2 * np.var(G,axis=0)
    var_genes = np.zeros(len(beta)) + varComps[1]
    var_covs  = np.zeros(len(beta))
    if Covs is not None: var_covs += np.dot(Covs, vd.getWeights()).var()
   

    return pv, beta, var_snps, var_covs, var_genes
 
def test_lmm_lr_speed(G,y,Z,Kbg,Covs=None,S=None,U=None):
    """
    low-rank lmm

    input:
    G   :   genotypes
    y   :   phenotype
    Z   :   features of low-rank matrix
    Kbg   :   background covariance matrix
    Covs :  fixed effect covariates

    using mtset implementation
    """
    m = mean(y)
    one = np.ones((1,1))

    if Z.shape[1] > G.shape[0]:
        return test_lmm_lr(G, y, Z, Kbg, Covs=Covs)
    
    if Covs is not None:
        m.addFixedEffect(Covs)
        nCovs = Covs.shape[1]

    Cg = freeform(1)
    Cn = freeform(1)

    Z/=np.sqrt(np.mean(np.sum(Z**2, axis=1)))
    #Z/=np.sqrt(Z.shape[1])
    gp = gp3kronSum(m,Cg,Cn,XX=Kbg,Xr=Z,S_XX=S,U_XX=U)

    
    params_rnd = {}
    params_rnd['Cg'] = 1e-4*np.random.randn(1)
    params_rnd['Cn'] = 1e-4*np.random.randn(1)
    params_rnd['Cr'] = 1e-4*np.random.randn(1)
    if Covs is not None:
        params_rnd['mean'] = 1e-6*np.random.randn(nCovs)
    
    conv,info = OPT.opt_hyper(gp,params_rnd)
        
    LML0 = gp.LML()
    params0 = gp.getParams()

    params_rnd = params0.copy()
    if Covs is not None:
        mean0 = params0['mean']
        params_rnd['mean'] = 1e-6*np.random.randn(nCovs+1)
        params_rnd['mean'][:nCovs] = mean0
    else:
        params_rnd['mean'] = 1e-6*np.random.randn(1)
    
    F = G.shape[1]
    LML = np.zeros(F)
    beta = np.zeros(F)

    var_snps = np.zeros(F)
    var_covs = np.zeros(F)
    var_genes = np.zeros(F)    

    for f in xrange(F):
        m.clearFixedEffect()
        if Covs is not None: m.addFixedEffect(Covs)
        m.addFixedEffect(G[:,[f]])
        conv,info = OPT.opt_hyper(gp, params_rnd)
        beta[f] = m.getParams()[-1]
        LML[f] = gp.LML()

        var_genes[f]  = gp.getParams()['Cr']**2
        var_snps[f]   = beta[f]**2 * np.var(G[:,[f]])
        if Covs is not None: 
            var_covs[f] = np.var(np.dot(Covs, m.getParams()[:-1]))
    

    LRT = 2*(LML0-LML)
    pv = stats.chi2.sf(LRT,1)
    return pv, beta, var_snps, var_covs, var_genes
