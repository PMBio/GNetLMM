import unittest
import filecmp
import numpy as np

from GNetLMM.pycore.utils.simPhenoCore import simPheno
from GNetLMM.pycore.utils.preprocessCore import preprocess
from GNetLMM.pycore.utils.analyseCore import analyse
from GNetLMM.pycore.utils.postprocessCore import postprocess



class BlankObject(object):
	pass 

class TestMain(unittest.TestCase):
	
	def test_simpheno(self):
		
		bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		ffile = './../data/1000G_chr22/ones.txt'
		pfile = './out/pheno'
		
		options = BlankObject()
		options.bfile = bfile 
		options.pfile = pfile 
		options.T = 30
		options.varSnp = 0.1 
		options.varNetwork = 0.8
		options.alpha = 1 
		options.nConfounder = 3 
		options.confPerGene = 1 
		options.expN = 1 
		options.seed = 0
		options.networkDesign = 'sparse'
		
		simPheno(options)
		
		assert filecmp.cmp('out/pheno.cols','comp/pheno.cols')
		assert filecmp.cmp('out/pheno.Agene','comp/pheno.Agene')
		assert filecmp.cmp('out/pheno.conf','comp/pheno.conf')
		assert filecmp.cmp('out/pheno.matrix','comp/pheno.matrix')
		assert filecmp.cmp('out/pheno.rows','comp/pheno.rows')
		
	def test_cov(self):
	
		bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		ffile = './../data/1000G_chr22/ones.txt'
		cfile = './out/chrom22'
		
		options = BlankObject()
		options.bfile = bfile 
		options.cfile = cfile 
		options.compute_cov = True 
		options.sim_type = 'RRM'
		options.plink_path = 'plink'
		
		preprocess(options)
		assert filecmp.cmp('out/chrom22.cov','comp/chrom22.cov')
		assert filecmp.cmp('out/chrom22.cov.id','comp/chrom22.cov.id')

		
	def test_init_assoc(self):

		options = BlankObject()
		options.initial_scan = True 
		options.bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		options.pfile = './comp/pheno'
		options.cfile = './comp/chrom22.cov'
		options.ffile = './../data/1000G_chr22/ones.txt'
		options.assoc0file = './out/lmm'
		options.nSnps = 10000
		options.memory_efficient = False
		options.merge_assoc0_scan = False 
		options.gene_corr = False
		options.merge_corr = False 
		options.compute_anchors = False 
		options.find_vstructures = False 
		options.update_assoc = False 
		options.block_assoc = False
		
		i = 0
		while i < 50000:
			options.startSnpIdx = i 
			analyse(options)
			
			assert filecmp.cmp('out/lmm.startSnp_%d.beta.matrix'%i,'comp/lmm.startSnp_%d.beta.matrix'%i)
			assert filecmp.cmp('out/lmm.startSnp_%d.pv.matrix'%i,'comp/lmm.startSnp_%d.pv.matrix'%i)
			i+= 10000
			
		options.initial_scan = False 
		options.merge_assoc0_scan = True 
		analyse(options)
		
		assert filecmp.cmp('out/lmm.beta.matrix','comp/lmm.beta.matrix')
		assert filecmp.cmp('out/lmm.pv.matrix','comp/lmm.pv.matrix')
			

	def test_gene_corr(self):

		options = BlankObject()
		options.bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		options.pfile = './comp/pheno'
		options.cfile = './comp/chrom22.cov'
		options.ffile = './../data/1000G_chr22/ones.txt'
		options.gfile = './out/genes'
		options.initial_scan = False
		options.nSnps = 10000
		options.memory_efficient = False
		options.merge_assoc0_scan = False 
		options.gene_corr = True 
		options.merge_corr = False 
		options.compute_anchors = False 
		options.find_vstructures = False 
		options.update_assoc = False 
		options.block_assoc = False
		options.nTraits = 10 
		
		
		# Compute marginal gene-gene correlations
		i=0 
		while i<=30: 
			options.startTraitIdx = i 
			analyse(options)
			
			assert filecmp.cmp('out/genes.startTrait_%d.corr.matrix'%i,'comp/genes.startTrait_%d.corr.matrix'%i)
			assert filecmp.cmp('out/genes.startTrait_%d.pv.matrix'%i,'comp/genes.startTrait_%d.pv.matrix'%i)
			i+= 10 
			
			
		options.gene_corr = False 
		options.merge_corr = True 
		analyse(options)

		assert filecmp.cmp('out/genes.corr.matrix','comp/genes.corr.matrix')
		assert filecmp.cmp('out/genes.pv.matrix','comp/genes.pv.matrix')
		

	def test_anchors(self):
	
		options = BlankObject()
		options.bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		options.pfile = './comp/pheno'
		options.cfile = './comp/chrom22.cov'
		options.ffile = './../data/1000G_chr22/ones.txt'
		options.gfile = './out/genes'
		options.assoc0file = './comp/lmm'
		options.anchorfile = './out/anchors.txt'
		options.initial_scan = False
		options.nSnps = 10000
		options.memory_efficient = False
		options.merge_assoc0_scan = False 
		options.gene_corr = False
		options.merge_corr = False 
		options.compute_anchors = True
		options.find_vstructures = False 
		options.update_assoc = False 
		options.block_assoc = False
		options.nTraits = 10 
		options.cis = True 
		options.trans = False 
		options.anchor_thresh = 1e-6 
		options.window = 2000
		
		analyse(options)
		assert filecmp.cmp('out/anchors.txt','comp/anchors.txt')
		
	
	def test_find_vstructures(self):
	
		
		options = BlankObject()
		options.bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		options.pfile = './comp/pheno'
		options.cfile = './comp/chrom22.cov'
		options.ffile = './../data/1000G_chr22/ones.txt'
		options.gfile = './comp/genes'
		options.assoc0file = './comp/lmm'
		options.anchorfile = './comp/anchors.txt'
		options.vfile = 'out/vstructures'
		options.initial_scan = False
		options.nSnps = 10000
		options.memory_efficient = False
		options.merge_assoc0_scan = False 
		options.gene_corr = False
		options.merge_corr = False 
		options.compute_anchors = False 
		options.find_vstructures = True 
		options.update_assoc = False 
		options.block_assoc = False
		options.nTraits = 10 
		options.cis = True 
		options.trans = False 
		options.anchor_thresh = 1e-6 
		options.window = 2000
		options.corr_thresh = 0.01 
		options.ind_thresh = 0.1
		options.max_genes = 10 
		
		i = 0 
		while i < 40:
			options.startTraitIdx = i 
			analyse(options)
			
			assert filecmp.cmp('out/vstructures.startTrait_%d.csv'%i,'comp/vstructures.startTrait_%d.csv'%i)
			
			i+= 10 

			
		options = BlankObject()
		options.merge_assoc = False 
		options.concatenate = True 
		options.plot_power = False 
		options.nice_output = False 
		options.infiles = 'out/vstructures'
		options.outfile = 'out/vstructures'

		postprocess(options)
		assert filecmp.cmp('out/vstructures.csv','comp/vstructures.csv')
		
	def test_update_assoc(self):
	
		options = BlankObject()
		options.bfile = './../data/1000G_chr22/chrom22_subsample20_maf0.10'
		options.pfile = './comp/pheno'
		options.cfile = './comp/chrom22.cov'
		options.ffile = './../data/1000G_chr22/ones.txt'
		options.gfile = './comp/genes'
		options.assoc0file = './comp/lmm'
		options.assocfile = './out/gnetlmm'
		options.anchorfile = './comp/anchors.txt'
		options.vfile = 'comp/vstructures'
		options.initial_scan = False
		options.nSnps = 10000
		options.memory_efficient = False
		options.merge_assoc0_scan = False 
		options.gene_corr = False
		options.merge_corr = False 
		options.compute_anchors = False 
		options.find_vstructures = False
		options.update_assoc = True
		options.block_assoc = False
		options.nTraits = 10 
		options.cis = True 
		options.trans = False 
		options.anchor_thresh = 1e-6 
		options.window = 2000
		options.corr_thresh = 0.01 
		options.ind_thresh = 0.1
		options.max_genes = 10 
		
		i = 0 
		while i < 40:
			options.startTraitIdx = i 
			analyse(options)
			
			M1 = np.loadtxt('out/gnetlmm.startTrait_%d.csv'%i)
			M2 = np.loadtxt('comp/gnetlmm.startTrait_%d.csv'%i)
			
			assert M1.shape==M2.shape 
			if M1.ndim>1:
				assert np.allclose(M1[:,:2],M2[:,:2])
		
			i+= 10 


		options = BlankObject()
		options.merge_assoc = False 
		options.concatenate = True 
		options.plot_power = False 
		options.nice_output = False 
		options.infiles = 'comp/gnetlmm'
		options.outfile = 'out/gnetlmm'
		postprocess(options)
		
		M1 = np.loadtxt('out/gnetlmm.csv')
		M2 = np.loadtxt('comp/gnetlmm.csv')
		assert np.allclose(M1[:,:2],M2[:,:2])
			
		options.concatenate = False 
		options.merge_assoc = True 
		options.assoc0file = './comp/lmm'
		options.assocfile = './out/gnetlmm'
		postprocess(options)
		M1 = np.loadtxt('comp/gnetlmm.beta.matrix')
		M2 = np.loadtxt('out/gnetlmm.beta.matrix')
		assert np.allclose(M1,M2)
	

		
if __name__ =='__main__':
	unittest.main()