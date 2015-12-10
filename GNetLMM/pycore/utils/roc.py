##ROC
#helper class for the evaluation of ROC curves and area under ROC

import scipy as S
import numpy as N
import pdb


def roc(labels, predictions):
    """roc - calculate receiver operator curve
    labels: true labels (>0 : True, else False)
    predictions: the ranking generated from whatever predictor is used"""
    #1. convert to arrays
    labels = S.array(labels).reshape([-1])
    predictions = S.array(predictions).reshape([-1])

    #threshold
    t = labels>0
    
    #sort predictions in descending order
    #get order implied by predictor (descending)
    Ix = S.argsort(predictions)[::-1]
    #reorder truth
    t = t[Ix]

    #compute true positiive and false positive rates
    tpr = S.double(N.cumsum(t))/t.sum()
    fpr = S.double(N.cumsum(~t))/(~t).sum()

    #add end points
    tpr = S.concatenate(([0],tpr,[1]))
    fpr = S.concatenate(([0],fpr,[1]))

    return [fpr,tpr]

def pr(labels, predictions):
    #1. convert to arrays
    labels = S.array(labels).reshape([-1])
    predictions = S.array(predictions).reshape([-1])
    #threshold
    t = labels>0
    Ix = S.argsort(predictions)[::-1]
    #reorder truth
    t = t[Ix]
    pr =  S.double(N.cumsum(t))/(N.cumsum(t)+N.cumsum(~t))
    rr =  S.double(N.cumsum(t))/(N.cumsum(t)+((t).sum()-N.cumsum(t)))
    return [rr,pr]

def auroc(labels=None,predictions=None,tp=None,fp=None):
    """auroc - calculate area under the curve from a given fp/rp plot"""

    if labels is not None:
        [fp, tp] = roc(labels,predictions)
    n = tp.size
    auc = 0.5*((fp[1:n]-fp[0:n-1]) * (tp[1:n]+tp[0:n-1])).sum()
    return auc
    pass

def aupr(labels=None,predictions=None,tp=None,fp=None):
    """auroc - calculate area under the curve from a given precision-recall plot"""
    if labels is not None:
        [rr, prr] = pr(labels,predictions)
    n = prr.size
    auc = 0.5*((rr[1:n]-rr[0:n-1]) * (prr[1:n]+prr[0:n-1])).sum()
    return auc
    pass


def AUC(qv, true_assoc):    
    sens, spec = roc(true_assoc.flatten(), -qv.flatten())
    ydata = S.array(sens)
    xdata = S.array(spec)	    
    
    x = (N.roll(xdata, -1) - xdata)[:-1]
    y = (N.roll(ydata, -1) + ydata)[:-1]/2
    auc = sum(map(lambda x, y: x*y, x, y))
    return auc
    

