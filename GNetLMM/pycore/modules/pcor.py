import scipy as SP
import math

from scipy.stats import t
import scipy.linalg as LA
import scipy.stats as STATS

import limix.modules.varianceDecomposition as varianceDecomposition
import limix.modules.qtl as qtl
import limix

import warnings



        
def pcor(X,Y,Z):
    """
    computes the correlation amtrix of X and Y conditioning on Z
    """
    if X.ndim==1: X = X[:,SP.newaxis]
    if Y.ndim==1: Y = Y[:,SP.newaxis]
    
    if Z is None: return STATS.pearsonr(X,Y)

    if Z.ndim==1: Z = Z[:,SP.newaxis]
    nSamples = X.shape[0]
    betaX, _, _, _ = LA.lstsq(Z,X)
    betaY, _, _, _ = LA.lstsq(Z,Y)
    Xres = X - SP.dot(Z,betaX)
    Yres = Y - SP.dot(Z,betaY)
    corr_cond = SP.corrcoef(Xres[:,0],Yres[:,0])[0,1]
    dz = Z.shape[1]  # dimension of conditioning variable
    df = max(nSamples - dz - 2,0)  # degrees of freedom

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        tstat = corr_cond / SP.sqrt(1.0 - corr_cond ** 2)  # calculate t statistic
        
    tstat = math.sqrt(df) * tstat
    pv_cond = 2 * t.cdf(-abs(tstat), df, loc=0, scale=1)  # calculate p value
    return corr_cond,pv_cond



def pcorParallel(X,Z,Y=None):
    """
    computes the correlation matrix between X and Y conditioning on Z
    """
    if Y is None: return pcorParallelSym(X,Z)
    if Z is None: return corrParallel(X,Y)

 
    if Z.ndim==1: Z = Z[SP.newaxis,:]

    X = X.T
    Y = Y.T
    Z = Z.T
    
    beta,_,_,_ = LA.lstsq(Z,Y)
    Yres = Y - SP.dot(Z,beta)

    beta,_,_,_ = LA.lstsq(Z,X)
    Xres = X - SP.dot(Z,beta)

    
    nSamples = Z.shape[0]
    nCovs = Z.shape[1]
    df = max(nSamples  - 2 - nCovs,0)
    
    return corrParallel(Xres.T,Yres.T,df=df)
    
  

def pcorParallelSym(Y,Z):
    """
    computes the correlation matrix of Y conditioning on Z
    """
    if Y.ndim==1: Y = Y[SP.newaxis,:]
    if Z is None: return corrParallel(Y,Z)

    if Z.ndim==1: Z = Z[SP.newaxis,:]
            
    Y = Y.T
    Z = Z.T
    beta,_,_,_ = LA.lstsq(Z,Y)
    Yres = Y - SP.dot(Z,beta)


    nSamples = Z.shape[0]
    nCovs = Z.shape[1]
    df = max(nSamples  - 2 - nCovs,0)

  
    return corrParallelSym(Yres.T,df=df)

 

def corrParallelSym(Y,df=None):
    """
    computes the correlation matric of Y
    """
    nSamples = Y.shape[1]
    corr = SP.corrcoef(Y)
    if df is None:
        df = max(nSamples  - 2,0)  # degrees of freedom

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        tstat = corr / SP.sqrt(1.0 - corr ** 2)  # calculate t statistic

    tstat = math.sqrt(df) * tstat
    pv = 2 * t.cdf(-abs(tstat), df, loc=0, scale=1)  # calculate p value
    return corr,pv

def corrParallel(X,Y=None,df=None):
    """
    computes the mxk correlation matrix between the mxn matrix X and the kxn matrix Z
    """    
    if Y is None:
        return corrParallelSym(X,df=df)

    assert X.shape[1]==Y.shape[1], 'ouch, samples do not match'
    nSamples = X.shape[1]
    
    Xstd = X.T
    Xstd-= Xstd.mean(0)
    Xstd/= Xstd.std(0)

    Ystd = Y.T
    Ystd-= Ystd.mean(0)
    Ystd/= Ystd.std(0)

    corr =  SP.dot(Xstd.T,Ystd)/nSamples
    if df is None:
        df = max(nSamples  - 2,0)  # degrees of freedom

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        tstat = corr / SP.sqrt(1.0 - corr ** 2)  # calculate t statistic

        
  
    tstat = math.sqrt(df) * tstat
    pv = 2 * t.cdf(-abs(tstat), df, loc=0, scale=1)  # calculate p value
    return corr,pv


