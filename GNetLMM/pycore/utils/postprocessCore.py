from optparse import OptionParser
import pdb
import numpy as np
import time
import glob

import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.io.bedReader as bedReader
import GNetLMM.pycore.io.phenoReaderFile as phenoReaderFile
import GNetLMM.pycore.utils.utils as utils
import GNetLMM.pycore.utils.roc as roc
import GNetLMM.pycore.utils.plotting as plotting
import GNetLMM.pycore.modules.assoc_results as assoc_results




def merge_files(fns_in, fn_out):
    """
    writing all files fns_ins into fn_out
    """
    with open(fn_out, 'w') as fout:
        for fn_in in fns_in:
            with open(fn_in) as fin:
                for line in fin:
                    fout.write(line)

                    
def concatenate(infiles, outfile):
    fns_in = glob.glob(infiles)
    merge_files(fns_in, outfile)


def merge_results(assocfile, assoc0file):
    """
    merging association updates with lmm-scan
    """
    results = assoc_results.AssocResultsList()
    results.load_csv(assocfile + '.csv')
    results.save_matrix(assoc0file, assocfile)




def get_groundtruth(Agene, snp_pos, snp_chrom, gene_pos, gene_chrom, window):

    T = len(gene_pos)
    F = len(snp_pos)

    P_cis   = np.zeros((F,T),dtype=bool)
    P_trans = np.zeros((F,T),dtype=bool)
    for t in range(T):
        icis = utils.getSnpIndicesWithinGeneWindow(snp_chrom,snp_pos,gene_chrom[t],gene_pos[t],window=window)         
        P_cis[icis,t]= True
        itrans = np.zeros(F,dtype=bool)
        parents = np.nonzero(Agene[t])[0]
        for ipar in parents:
            if ipar==t: continue
            itrans += utils.getSnpIndicesWithinGeneWindow(snp_chrom,snp_pos,gene_chrom[ipar],gene_pos[ipar],
                                                          window=window)
            
        P_trans[itrans,t] = True

    return P_cis, P_trans
    
def plot_power(bfile,pfile,assoc0file, assocfile, plotfile, window):
    # reading out p-values
    score = {}
    assoc0Reader = reader.FileReader(assoc0file + '.pv')
    score['LMM'] = -np.log10(assoc0Reader.getMatrix())
    assocReader = reader.FileReader(assocfile + '.pv')
    score['GNetLMM'] = -np.log10(assocReader.getMatrix())

    # get network
    Agene = np.loadtxt(pfile + '.Agene')

    # gene info
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    gene_start = preader.getGeneStart()
    gene_end = preader.getGeneEnd()
    gene_chrom = preader.getGeneChrom()
    
    # snp info
    breader = bedReader.BedReader(bfile)
    snp_pos = breader.getSnpPos()
    snp_chrom = breader.getSnpChrom()
    

    
    P_cis, P_trans = get_groundtruth(Agene, snp_pos, snp_chrom, gene_start, gene_chrom, window)

    # delete cis-associations
    for key in score.keys():
        score[key][P_cis] = 0

    # compute receiver operator characteristics
    FPR = {}
    TPR = {}
    for key in score.keys():
        FPR[key], TPR[key] = roc.roc(P_trans, score[key])

    # plotting results
    plotting.plotROCcurve(['LMM','GNetLMM'],TPR,FPR,xlim=(0,0.05),ylim=(0,0.42),fn=plotfile)




def postprocess(options):

    """ concatenate files """
    if options.concatenate:
        t0 = time.time()
        print 'Concatenating files'
        assert options.infiles is not None, 'Please specify the regular expression for the input filename.'
        assert options.outfiles is not None, 'Please specify the output filename.'
        concatenate(options.infiles, options.outfile)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)




    """ merging results """
    if options.merge_assoc:
        t0 = time.time()
        print 'Merging associations'
        assert options.assocfile is not None, 'Please specify an output file'
        assert options.assoc0file is not None, 'Please specify assoc0 results'
        merge_results(options.assocfile, options.assoc0file)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
  


    if options.plot_power:
        assert options.assoc0file is not None, 'Please specify the assoc0 file'
        assert options.assocfile is not None, 'Please specify the assoc0 file'
        assert options.pfile is not None, 'Please specify the pheno file'
        assert options.plotfile is not None, 'Please specify the file for saving the figure'
        assert options.bfile is not None, 'Please specify the geno file'

        t0 = time.time()
        print 'Plotting power'
        plot_power(options.bfile,options.pfile,options.assoc0file, options.assocfile,options.plotfile,options.window)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)