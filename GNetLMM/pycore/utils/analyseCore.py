
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
import pdb
import glob

import GNetLMM.pycore.io.bedReader as bedReader
import GNetLMM.pycore.io.phenoReaderFile as phenoReaderFile
import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.io.writer as writer

import GNetLMM.pycore.modules.gnetlmm as gnetlmm
import GNetLMM.pycore.modules.assoc_results as assoc_results


def merge_assoc0_scan(assoc0file,nSnps,bfile):
    """
    merging associations files

    input:
    assoc0   :   basename of assoc0 results
    nSnps    :   number of SNPs in each block
    bfile    :   basename of bed file
    """
    greader =  bedReader.BedReader(bfile)
    F = greader.get_nrows()

    fn_beta = []
    fn_pv   = []
    
    f = 0
    while f < F:

        _fn = assoc0file + '.startSnp_%d'%f
        _fn_beta = '%s.beta.matrix'%_fn
        _fn_pv   = '%s.pv.matrix'%_fn
        
        assert os.path.exists(_fn_beta), 'File %s is missing'%(_fn_beta)
        assert os.path.exists(_fn_pv), 'File %s is missing'%(_fn_pv)
        fn_beta.append(_fn_beta)
        fn_pv.append(_fn_pv)
        f += nSnps

    
    merge_files(fn_beta, assoc0file + '.beta.matrix')
    merge_files(fn_pv, assoc0file + '.pv.matrix')


def merge_genecorr(pfile, gfile, nTraits):
    """
    merging gene-gene correlations file 

    input:
    pfile    : basename of the phenotype file
    gfile    : basename of the correlation file
    nTraits  : number of traits in each block
    """
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    T = preader.get_nrows()

    t = 0
    fn_corr = []
    fn_pv = []
    while t < T:

        _fn = gfile + '.startTrait_%d'%t
        _fn_corr = '%s.corr.matrix'%_fn
        _fn_pv   = '%s.pv.matrix'%_fn
        
        assert os.path.exists(_fn_corr), 'File %s is missing'%(_fn_corr)
        assert os.path.exists(_fn_pv), 'File %s is missing'%(_fn_pv)
        fn_corr.append(_fn_corr)
        fn_pv.append(_fn_pv)
        t += nTraits

    
    merge_files(fn_corr, gfile + '.corr.matrix')
    merge_files(fn_pv, gfile + '.pv.matrix')


def merge_files(fns_in, fn_out):
    """
    writing all files fns_ins into fn_out
    """
    with open(fn_out, 'w') as fout:
        for fn_in in fns_in:
            with open(fn_in) as fin:
                for line in fin:
                    fout.write(line)
                    

    

def scan(bfile,pfile,cfile,ffile,vfile,assocfile,startTraitIdx,nTraits):
    """
    running association scan

    input:
    bfile      :   basefilename of plink file
    pfile      :   phenotype file
    cfile      :   covariance file
    ffile      :   fixed effects file
    vfile      :   file containing vstructures
    assocfile  :   file for saving results
    """
    K = None
    if cfile is not None:
        K = np.loadtxt(cfile)

    Covs = None
    if ffile is not None:
        Covs = np.loadtxt(ffile)
        
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader, Covs=Covs,K=K)
    model.load_vstructures(vfile+".csv")
    model.update_associations(startTraitIdx, nTraits)
    model.save_updates(assocfile)

    

def find_vstructures(bfile, pfile,gfile,anchorfile, assoc0file,window,vfile,startTraitIdx,nTraits, corr_thresh, ind_thresh,max_genes):
    """
    running association scan

    input:
    bfile      :   blink file
    pfile      :   phenotype file
  
    gfile      :   correlation file
    anchorfile :   file containing anchors
    assoc0file :   file containing the results from the initial association scan
    window     :   window size
    vfile      :   file containing v-structures

    startTraitIdx   :   index of first trait to be analyzed
    nTraits         :   number of traits to be analyzed

    corr_thresh   : q-value for calling a correlation significant
    ind_thresh    : q-value for calling a correlation not significant
    max_genes     : maximal number of genes in the conditioning set
    """
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    
    model = gnetlmm.GNetLMM(preader, greader, window=window, thresh_corr=corr_thresh, thresh_ind=ind_thresh)

    genecorr_reader = reader.FileReader(gfile + '.pv')
    model.set_genecorr_reader(genecorr_reader)
    
    assoc0Reader = reader.FileReader(assoc0file + '.pv')
    model.set_assoc0_reader(assoc0Reader)
    
    model.load_anchors(anchorfile)
    model.find_vstructures(startTraitIdx, nTraits, max_genes)
    model.save_vstructures(vfile+'.csv')




def initial_scan(bfile, pfile, cfile, ffile, assoc0file, startSnpIdx=0, nSnps=np.inf,
                 memory_efficient = False):
    """
    running initial scan using a standard linear mixed model

    Input:
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    pfile        :   phenotype file
    cfile        :   covariance matrix file
    ffile        :   covariates file

    assoc0file   :   basename of output file 
    """
    K = None
    if cfile is not None:
        K = np.loadtxt(cfile)

    Covs = None
    if ffile is not None:
        Covs = np.loadtxt(ffile)

    preader = phenoReaderFile.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader,K=K,Covs=Covs)
    beta0, pv0 = model.initial_scan(startSnpIdx, nSnps,memory_efficient)

    write = writer.Writer(assoc0file+'.pv')
    write.writeMatrix(pv0, fmt='%.4e')

    write = writer.Writer(assoc0file+'.beta')
    write.writeMatrix(beta0, fmt='%.4f')

   

def marginal_genecorr(pfile, gfile, startTraitIdx=0, nTraits=np.inf):
    """
    running marginal gene-gene correlations

    Input:
    pfile        :   phenotype file
    cfile        :   covariance matrix file
    ffile        :   covariates file

    gfile        :   basename of output file 

    startTraidIdx : first trait to be analyses 
    nTraits       : number of traits to be analysed
    """
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    model = gnetlmm.GNetLMM(preader,None)
    corr, pv = model.marginal_gene_correlations(startTraitIdx=startTraitIdx, nTraits=nTraits)
    write = writer.Writer(gfile+'.pv')
    write.writeMatrix(pv, fmt='%.4e')
    write = writer.Writer(gfile+'.corr')
    write.writeMatrix(corr, fmt='%.4f')
    
  
def gene_has_anchor(bfile, pfile, assoc0, anchor_thresh, anchorfile,window,cis=True):
    """
    tests if a gene has a cis anchor

    input:
    bfile           :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    pfile           :   phenotype file
    assoc0file      :   basefilename of initial association scan
    anchor_thresh   :   thrshold for anchor-associations
    anchor_file     :   filename for saving cis assocaitions
    window          :   maximal distance between cis-snp and gene
    cis             :   if set, look for cis-associations only (default: True)
    """
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    
    model = gnetlmm.GNetLMM(preader, greader, window=window)
    assoc0Reader = reader.FileReader(assoc0 + '.pv')
    model.set_assoc0_reader(assoc0Reader)
    model.gene_has_anchor(anchor_thresh,cis)
    model.save_anchors(anchorfile)


    
def analyse(options):

    """ running initial scan """
    if options.initial_scan:
        assert options.bfile is not None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify a phenotypic file.'
        assert options.assoc0file is not None, 'Please specify an output file for saving the associations results.'
        t0 = time.time()
        print 'Running initial scan'
        initial_scan(options.bfile,options.pfile,options.cfile,options.ffile, options.assoc0file,startSnpIdx=options.startSnpIdx, nSnps=options.nSnps,
                     memory_efficient=options.memory_efficient)
        t1 = time.time()
        print '... finished in %s seconds'%(t1-t0)


    """ merging assoc0 files """
    if options.merge_assoc0_scan:
        assert options.assoc0file is not None, 'Please specify assoc0 basefilename.'
        assert options.nSnps is not None, 'Please specify number of SNPs.'
        assert options.bfile is not None, 'Please specify bed basefilename.'
        t0 = time.time()
        print 'Merge assoc0 files'
        merge_assoc0_scan(options.assoc0file,options.nSnps,options.bfile)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)


    """ computing marginal gene-gene correlations """
    if options.gene_corr:
        assert options.pfile!=None, "Please specify a phenotypic file."
        assert options.gfile!=None, "Please specify an output file for saving the gene-gene correlations"
        t0 = time.time()
        print "Computing marginal gene-gene correlations"
        marginal_genecorr(options.pfile,options.gfile, startTraitIdx=options.startTraitIdx, nTraits=options.nTraits)
        t1=time.time()
        print '... finished in %s seconds'%(t1-t0)

    """ merging correlation files """
    if options.merge_corr:
        assert options.gfile is not None, "Please specify an output file for saving the gene-gene correlations."
        assert options.nTraits is not None, "Please specify the number of traits."
        assert options.pfile is not None, "Please specify a phenotypic file."
        t0 = time.time()
        print "Merging gene-gene correlation"
        merge_genecorr(options.pfile, options.gfile, nTraits=options.nTraits)
        t1 = time.time()
        print '... finished in %s seconds'%(t1-t0)


    """ determining cis anchors """
    if options.compute_anchors:
        assert options.bfile!=None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify the phenotype file'
        assert options.assoc0file is not None, "Please specify assoc0 basefilename."
        assert options.anchor_thresh is not None, "Please specify an anchor threshold."
        assert options.anchorfile is not None, 'Please specify file for saving the anchors.'
        assert (options.cis or options.trans), 'Please specify if you look for cis- or trans-associations.'
        assert options.cis*options.trans==0, 'Please specify if you look for cis- or trans-associations.'
    
        t0 = time.time()
        print "Determining anchors"
        gene_has_anchor(options.bfile,options.pfile,options.assoc0file, options.anchor_thresh,options.anchorfile,options.window,options.cis)
        t1=time.time()
        print '... finished in %s seconds'%(t1-t0)

    """ finding v-structures """
    if options.find_vstructures:
        t0 = time.time()
        print 'Finding v-structures'
        assert options.bfile!=None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify the phenotype file'
        assert options.gfile is not None, 'Please specify the gene-gene correlation file'
        assert options.anchorfile is not None, 'Please specify the cis-anchor file'
        assert options.assoc0file is not None, 'Please specify the assoc0 file'
        assert options.vfile is not None, 'Please specify output file for vstructures'
        find_vstructures(options.bfile, options.pfile, options.gfile, options.anchorfile,options.assoc0file, options.window, options.vfile,options.startTraitIdx,options.nTraits, options.corr_thresh, options.ind_thresh, options.max_genes)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
    """ updating associations """
    if options.update_assoc:
        t0 = time.time()
        print 'Updating associations'
        assert options.bfile!=None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify the phenotype file'
        assert options.vfile is not None, 'Please specify vstructure file'
        assert options.assocfile is not None, 'Please specify an output file'
        scan(options.bfile,options.pfile,options.cfile,options.ffile,options.vfile, options.assocfile,options.startTraitIdx,options.nTraits)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
