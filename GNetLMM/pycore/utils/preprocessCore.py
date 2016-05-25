
import os
import subprocess
import pdb
import sys
import numpy as np
import numpy.linalg as la
import scipy as sp
import warnings
from optparse import OptionParser
import time
import scipy as sp
import warnings
import scipy.sparse.linalg as ssla
import glob
import limix.io.plink as plink



import GNetLMM.pycore.io.bedReader as bedReader
#import GNetLMM.pycore.io.phenoReader as phenoReader
import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.io.writer as writer

    
def computeCovarianceMatrixPlink(plink_path,out_dir,bfile,cfile,sim_type='RRM'):
    """
    computing the covariance matrix via plink
    """
    
    print "Using plink to create covariance matrix"
    cmd = '%s --bfile %s '%(plink_path,bfile)

    if sim_type=='RRM':
        # using variance standardization
        cmd += '--make-rel square '
    else:
        raise Exception('sim_type %s is not known'%sim_type)

    cmd+= '--out %s'%(os.path.join(out_dir,'plink'))
    
    subprocess.call(cmd,shell=True)

    # move file to specified file
    if sim_type=='RRM':
        old_fn = os.path.join(out_dir, 'plink.rel')
        os.rename(old_fn,cfile+'.cov')
        
        old_fn = os.path.join(out_dir, 'plink.rel.id')
        os.rename(old_fn,cfile+'.cov.id')

    if sim_type=='IBS':
        old_fn = os.path.join(out_dir, 'plink.mibs')
        os.rename(old_fn,cfile+'.cov')

        old_fn = os.path.join(out_dir, 'plink.mibs.id')
        os.rename(old_fn,cfile+'.cov.id')

    os.remove(os.path.join(out_dir, 'plink.nosex'))
    os.remove(os.path.join(out_dir, 'plink.log'))                 
    

    
def computeCovarianceMatrixPython(out_dir,bfile,cfile,sim_type='RRM'):
    print "Using python to create covariance matrix. This might be slow. We recommend using plink instead."

    if sim_type is not 'RRM':
        raise Exception('sim_type %s is not known'%sim_type)

    """ loading data """
    data = plink.readBED(bfile)
    iid  = data['iid']
    X = data['snps']
    N = X.shape[1]
    print '%d variants loaded.'%N
    print '%d people loaded.'%X.shape[0]
    
    """ normalizing markers """
    print 'Normalizing SNPs...'
    p_ref = X.mean(axis=0)/2.
    X -= 2*p_ref

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X /= sp.sqrt(2*p_ref*(1-p_ref))
        
    hasNan = sp.any(sp.isnan(X),axis=0)
    print '%d SNPs have a nan entry. Exluding them for computing the covariance matrix.'%hasNan.sum()

    """ computing covariance matrix """
    print 'Computing relationship matrix...'
    K = sp.dot(X[:,~hasNan],X[:,~hasNan].T)
    K/= 1.*N
    print 'Relationship matrix calculation complete'
    print 'Relationship matrix written to %s.cov.'%cfile
    print 'IDs written to %s.cov.id.'%cfile

    """ saving to output """
    np.savetxt(cfile + '.cov', K, delimiter='\t',fmt='%.6f')
    np.savetxt(cfile + '.cov.id', iid, delimiter=' ',fmt='%s')
    




def computeCovarianceMatrix(plink_path,bfile,cfile,sim_type='RRM'):
    """
    compute similarity matrix using plink

    Input:
    plink_path   :   plink path
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    cfile        :   the covariance matrix will be written to cfile.cov and the corresponding identifiers
                         to cfile.cov.id. If not specified, the covariance matrix will be written to cfile.cov and
                         the individuals to cfile.cov.id in the current folder.
    sim_type     :   {IBS/RRM} are supported
    """
 
    try:
        output    = subprocess.check_output('%s --version --noweb'%plink_path,shell=True)
        use_plink = float(output.split(' ')[1][1:4])>=1.9
    except:
        use_plink = False

    assert bfile!=None, 'Path to bed-file is missing.'
    assert os.path.exists(bfile+'.bed'), '%s.bed is missing.'%bfile
    assert os.path.exists(bfile+'.bim'), '%s.bim is missing.'%bfile
    assert os.path.exists(bfile+'.fam'), '%s.fam is missing.'%bfile

    # create dir if it does not exist
    out_dir = os.path.split(cfile)[0]
    if out_dir!='' and (not os.path.exists(out_dir)):
        os.makedirs(out_dir)


    if use_plink:
        computeCovarianceMatrixPlink(plink_path,out_dir,bfile,cfile,sim_type=sim_type)
    else:
        computeCovarianceMatrixPython(out_dir,bfile,cfile,sim_type=sim_type)




def preprocess(options):

    """ computing the covariance matrix """
    if options.compute_cov:
       assert options.bfile!=None, 'Please specify a bfile.'
       assert options.cfile is not None, 'Specify covariance matrix basename'
       print 'Computing covariance matrix'
       t0 = time.time()
       computeCovarianceMatrix(options.plink_path,options.bfile,options.cfile,options.sim_type)
       t1 = time.time()
       print '... finished in %s seconds'%(t1-t0)
  
    
