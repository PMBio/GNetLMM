import scipy as sp
import numpy as np

import os

import GNetLMM.pycore.io.writer as writer
import GNetLMM.pycore.modules.simulator as simulator

def simPheno(options):
    assert options.bfile!=None, 'Please specify a bfile.'
    assert options.pfile!=None, 'Please specify a pfile.'
    
    sp.random.seed(options.seed)
    
    print 'simulating genes'
    if options.networkDesign=='star':
        sim = simulator.HotSpotSimulator(options.bfile)
    elif options.networkDesign=='sparse':
        sim = simulator.SparseSimulator(options.bfile)
    else:
        raise Exception("networkDesign '%s' is not known."%options.networkDesign)
    RV = sim.simulateGenes(T=options.T,varSnp=options.varSnp,varNetwork=options.varNetwork,expN=options.expN,
                           alpha=options.alpha,nConfounder=options.nConfounder,confPerGene=options.confPerGene)

    print 'exporting simulated data'
    outdir = os.path.split(options.pfile)[0]
    if not(os.path.exists(outdir)): os.mkdir(outdir)
    np.savetxt('%s.Aconf'%options.pfile, RV['Aconf'], fmt='%d')
    np.savetxt('%s.Agene'%options.pfile, RV['Agene'], fmt='%d')
    np.savetxt('%s.conf'%options.pfile, RV['H'], fmt='%.6f')

    write = writer.Writer(options.pfile)
    write.writeColumnInfo(data={'fid':sim.genoreader.fam[:,0]})
    write.writeRowInfo(data={'gene_ids':RV['gene_ids'], 'gene_chrom':RV['gene_chrom'], 'gene_start':RV['gene_start'],
                             'causal_rs':RV['rs']})
    write.writeMatrix(RV['Y'], fmt='%.6f')
 
