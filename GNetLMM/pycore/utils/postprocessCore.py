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
import GNetLMM.pycore.modules.vstructures as vstructures



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
    fns_in = glob.glob(infiles + ".startTrait_*.csv")
    merge_files(fns_in, outfile + '.csv')


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
    
def plot_power(bfile,pfile,assoc0file, assocfile, plotfile, window, blockfile=None):
    # reading out p-values
    score = {}
    assoc0Reader = reader.FileReader(assoc0file + '.pv')
    score['LMM'] = -np.log10(assoc0Reader.getMatrix())
    assocReader = reader.FileReader(assocfile + '.pv')
    score['GNetLMM'] = -np.log10(assocReader.getMatrix())

    if blockfile is not None:
        assocReader = reader.FileReader(blockfile + '.pv')
        score["BlockLMM"] = -np.log10(assocReader.getMatrix())

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
    plotting.plotROCcurve(score.keys(),TPR,FPR,xlim=(0,0.05),ylim=(0,0.42),fn=plotfile)



def create_nice_output(vfile, bfile, pfile, assoc0file, assocfile, blockfile, outfile):
    """
    creating human readbable output file
    """
    preader = phenoReaderFile.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)

    snp_ids  = greader.getSnpIds()
    gene_ids = preader.getGeneIds()

    vstruct = vstructures.VstructureFile(vfile + '.csv')
    assoc0Reader = reader.FileReader(assoc0file + '.pv')
    assocReader  = assoc_results.AssocResultsList()
    assocReader.load_csv(assocfile + '.csv')

    blockReader = None
    if blockfile is not None: 
        blockReader = assoc_results.AssocResultsList()
        blockReader.load_csv(blockfile + '.csv')

    f = open(outfile + '.csv','w')
    if blockReader is None:
        header = "Anchor Snp\t Anchor Gene\t Focal Gene\t Orthogonal Genes\t pv(LMM) \t pv(GNetLMM)"
    else:
        header = "Anchor Snp\t Anchor Gene\t Focal Gene\t Orthogonal Genes\t pv(LMM) \t pv(GNetLMM) \t pv(BlockLMM)"

    if assocReader.var_snps is not None: header += '\t varSnp(GnetLMM)'
    if assocReader.var_covs is not None: header += '\t varCovs(GnetLMM)'
    if assocReader.var_genes is not None: header += '\t varGenes(GnetLMM)'

    if blockReader is not None:
        if blockReader.var_snps is not None: header += '\t varSnp(BlockLMM)'
        if blockReader.var_covs is not None: header += '\t varCovs(BlockLMM)'
        if blockReader.var_genes is not None: header += '\t varGenes(BlockLMM)'

    header += "\n"
    f.write(header)
    for idx_focal_gene, idx_anchor_snp, idx_orth_gene, idx_anchor_gene in vstruct.iterator(full=True):

        focal_gene = ",".join(gene_ids[idx_focal_gene])
        orth_gene = ",".join(gene_ids[idx_orth_gene])
        anchor_snps = snp_ids[idx_anchor_snp]
        
        for i in range(len(anchor_snps)):

            try:
                anchor_gene = gene_ids[idx_anchor_gene[i]]
            except:
                anchor_gene = ""

        
            pv_lmm = assoc0Reader.getRows([idx_anchor_snp[i]])[:,idx_focal_gene][0,0]

            idx_gnet = np.logical_and(assocReader.focal_gene==idx_focal_gene, assocReader.snp_anchor==idx_anchor_snp[i])
            if blockReader is not None:
                idx_block = np.logical_and(blockReader.focal_gene==idx_focal_gene, blockReader.snp_anchor==idx_anchor_snp[i])

            if blockReader is None:
                line = "%s\t%s\t%s\t%s\t%.4e\t%.4e"%(anchor_snps[i],anchor_gene,focal_gene, orth_gene, pv_lmm, assocReader.pv[idx_gnet])
            else:
                line = "%s\t%s\t%s\t%s\t%.4e\t%.4e\t%.4e"%(anchor_snps[i],anchor_gene,focal_gene, orth_gene, pv_lmm, assocReader.pv[idx_gnet], blockReader.pv[idx_block])

            if assocReader.var_snps is not None: line += '\t%.4e'%assocReader.var_snps[idx_gnet]
            if assocReader.var_covs is not None: line += '\t%.4e'%assocReader.var_covs[idx_gnet]
            if assocReader.var_genes is not None: line += '\t%.4e'%assocReader.var_genes[idx_gnet]

            if blockReader is not None:
                if blockReader.var_snps is not None: line += '\t%.4e'%blockReader.var_snps[idx_block]
                if blockReader.var_covs is not None: line += '\t%.4e'%blockReader.var_covs[idx_block]
                if blockReader.var_genes is not None: line += '\t%.4e'%blockReader.var_genes[idx_block]

            line += "\n"
            f.write(line)


    f.close()


def postprocess(options):

    """ concatenate files """
    if options.concatenate:
        t0 = time.time()
        print 'Concatenating files'
        assert options.infiles is not None, 'Please specify the regular expression for the input filename.'
        assert options.outfile is not None, 'Please specify the output filename.'
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
        

    """ creating nice output file """
    if options.nice_output:
        t0 = time.time()
        print "Creating nice outputfile"
        assert options.bfile is not None, 'Please specify bed file'
        assert options.pfile is not None, 'Please specify phenotype file'
        assert options.vfile is not None, 'Please specify v-structures file'
        assert options.outfile is not None, 'Please specify output file'

        assert options.assocfile is not None, 'Please specify assoc file'
        assert options.assoc0file is not None, 'Please specify assoc0 file'
        create_nice_output(options.vfile, options.bfile, options.pfile, options.assoc0file, options.assocfile,
            options.blockfile, options.outfile)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
  
    """ plotting power on simulated data """
    if options.plot_power:
        assert options.assoc0file is not None, 'Please specify the assoc0 file'
        assert options.assocfile is not None, 'Please specify the assoc0 file'
        assert options.pfile is not None, 'Please specify the pheno file'
        assert options.plotfile is not None, 'Please specify the file for saving the figure'
        assert options.bfile is not None, 'Please specify the geno file'

        t0 = time.time()
        print 'Plotting power'
        plot_power(options.bfile,options.pfile,options.assoc0file, options.assocfile,options.plotfile,options.window, options.blockfile)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
