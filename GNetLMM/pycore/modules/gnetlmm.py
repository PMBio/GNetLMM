
import pdb
import numpy as np
import scipy.stats
import time

import limix.modules.qtl as qtl

import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.utils.utils as utils
import GNetLMM.pycore.utils.qvalue as qvalue
import GNetLMM.pycore.utils.qtl_lr as qtl_lr

import vstructures
import anchors
import assoc_results
import pcor

class GNetLMM:

    def __init__(self, phenoreader=None, genoreader=None,K=None,Covs=None, pv_genes=None, thresh_corr=0.01, thresh_ind=0.1,
                 window=None):
        """
        Constructor

        input:
        Y      : phenotype matrix [NxT]
        K      : covariance matrix [NxN]
        Covs   : fixed effects matrix [NxS]
        """
        self.phenoreader = phenoreader
        self.genoreader  = genoreader
        self.window = window

        if self.genoreader is not None:
            N = self.genoreader.get_ncols()
        if self.phenoreader is not None:
            N = self.phenoreader.get_ncols()
        
        self.K = None
        if K is not None:
            assert K.shape[0]==N
            assert K.shape[1]==N
            self.set_K(K)

        self.Covs = None
        if Covs is not None:
            assert Covs.shape[0]==N
            if Covs.ndim==1: Covs = Covs[:,np.newaxis]
            self.Covs = Covs
     
        if pv_genes is not None:
            assert pv_genes.shape[0]==T
            assert pv_genes.shape[1]==T
            self.pv_genes = pv_genes

        self.thresh_corr = thresh_corr
        self.thresh_ind  = thresh_ind
            

    def set_assoc0_reader(self,assoc0_reader):
        self.assoc0_reader = assoc0_reader

    def set_genecorr_reader(self,genecorr_reader):
        self.genecorr_reader = genecorr_reader

    def set_window(self,window):
        self.window = window
        
    def set_K(self,K=None, X=None, U=None, S=None):
        """
        setting background covariance matrix

        input:
        K   :   covariance matrix [NxN]
        X   :   genetric markers used for constructing the covariance marix [NxF].
        """
        N = self.genoreader.get_ncols()
        
        if K is not None:
            assert K.shape[0]==N
            assert K.shape[1]==N
            self.K = K

        if X is not None:
            assert X.shape[0]==N
            self.K = utils.computeLinearKernel(X)

        if (U is None) or (S is None):
            S,U = scipy.linalg.eigh(K)
        self.S = S
        self.U = U

    def initial_scan(self, startSnpIdx=0, nSnps=np.inf,memory_efficient=False):
        """
        running initial scan using a linear mixed model

        input:
        startSnpIdx   :    index of first snp (default : 0)
        nSnps         :    number of SNPs to use (default: infinite)
        memory_efficient : if turned on (default: false), phenotype are processed sequentially,
                           leading to longer runtime but less memory.
        """
        F = self.genoreader.get_nrows()
        T = self.phenoreader.get_nrows()
        
        N = self.genoreader.get_ncols()

        if ~np.isfinite(nSnps):
            nSnps = F

        nSnps = min(nSnps, F-startSnpIdx)
        
        G = self.genoreader.loadSnpBlock(startSnpIdx, nSnps)

        if memory_efficient:
            pv = np.zeros((nSnps,T))
            beta = np.zeros((nSnps,T))
            for t in range(T):
                y = self.phenoreader.getRows([t]).T
                lm = qtl.test_lmm(snps=G.T, pheno=y, K=self.K, covs=self.Covs)
                pv[:,t] = lm.getPv()[0]
                beta[:,t] = lm.getBetaSNP()[0]
        else:
            Y = self.phenoreader.getMatrix()
            lm = qtl.test_lmm(snps=G.T, pheno=Y.T, K=self.K, covs=self.Covs)
            pv = lm.getPv().T
            beta = lm.getBetaSNP().T
            
        
        self.assoc0_reader = reader.MatrixReader(pv)
        
        return beta, pv


    def marginal_gene_correlations(self, startTraitIdx=0, nTraits=np.inf):
        """
        running marginal gene-gene correlations
        """
        Y = self.phenoreader.getMatrix()

        if startTraitIdx==0 and np.isinf(nTraits):
            corr, pv = pcor.corrParallelSym(Y)
        else:
            corr, pv = pcor.corrParallel(Y[startTraitIdx:startTraitIdx+nTraits],Y)
        
        self.genecorr_reader = reader.MatrixReader(pv)
        
        return corr,pv


    def gene_has_anchor(self, thresh, cis=True):
        """
        computes if a gene has a cis anchor

        input:
        cis_thresh   :    threshold for cis-association
        cis_window   :    max distance between snp and gene
        """
        
        F = self.genoreader.get_nrows()
        T = self.phenoreader.get_nrows()
                
        snp_ids  = self.genoreader.getSnpIds()
        gene_ids = self.phenoreader.getGeneIds()
        
        RV = {'pv':[], 'snp_ids':[], 'gene_ids':[], 'isnp':[], 'igene':[]}

        for f, pv0_f in self.assoc0_reader.getRowIterator():
   
            pv_min = np.min(pv0_f)
            if pv_min > thresh: continue
            idx_anchor = pv0_f==pv_min
            if not(idx_anchor.any()): continue

            if cis:
                idx_anchor[idx_anchor] = self.find_cis_genes(f, idx_anchor)
     
            if idx_anchor.any():
                igenes = np.nonzero(idx_anchor)[0]
                for t in igenes:
                    RV['pv'].append(pv_min)
                    RV['snp_ids'].append(snp_ids[f])
                    RV['gene_ids'].append(gene_ids[t])
                    RV['isnp'].append(f)
                    RV['igene'].append(t)
          
        for key in RV.keys():
            RV[key] = np.array(RV[key])

        self.anchors = anchors.Anchors(F,T,pv=RV['pv'],snp_ids=RV['snp_ids'],gene_ids=RV['gene_ids'],
                                                  igene=RV['igene'],isnp=RV['isnp'])
        
            

    def save_anchors(self,fn):
        """
        save anchors in file
        """
        self.anchors.save(fn)

    def load_anchors(self,fn):
        """
        load anchors from file
        """
        F = self.genoreader.get_nrows()
        T = self.phenoreader.get_nrows()
        self.anchors = anchors.Anchors(F,T)
        
        self.anchors.load(fn)
        

    def find_cis_genes(self,isnp,igene=None):
        """
        finding genes in close proximity to the snp

        isnp         :   snp index
        window       :   defines close proximity
        """
        if igene is None:
            T = self.phenoreader.get_nrows()
            igene = np.ones(T, dtype=bool)
        snp_pos = self.genoreader.getSnpPos()[isnp]
        snp_chrom = self.genoreader.getSnpChrom()[isnp]
        gene_start = self.phenoreader.getGeneStart()[igene]
        gene_end   = self.phenoreader.getGeneEnd()[igene]
        gene_chrom = self.phenoreader.getGeneChrom()[igene]

        isel = gene_chrom==snp_chrom
        isel = np.logical_and(isel, gene_start - self.window <= snp_pos)
        isel = np.logical_and(isel, snp_pos <= gene_end + self.window)
        return isel



    def update_associations(self, startTraitIdx=0, nTraits=np.inf):
        """
        updating association scan
        """
        self.assoc_updates = assoc_results.AssocResultsList()

        focal_gene_prev = None
        
        for focal_gene, snp_anchor, orth_gene in self.vstructures.iterator():
            if focal_gene<startTraitIdx or startTraitIdx+nTraits<=focal_gene:
                continue

            if (focal_gene_prev is None) or (focal_gene!=focal_gene_prev):
                print ".... Updating associations for gene %d"%(focal_gene)
                focal_gene_prev = focal_gene
      
            y_focal  = self.phenoreader.getRows(focal_gene).T
            y_orth   = self.phenoreader.getRows(orth_gene).T
            if y_focal.ndim==1: y_focal = y_focal[:,np.newaxis]
            if y_orth.ndim==1: y_orth = y_orth[:,np.newaxis]
                
            startSnpIdx = np.min(snp_anchor)
            nSnps  = np.max(snp_anchor) - startSnpIdx + 1
            G_anchor = self.genoreader.loadSnpBlock(startSnpIdx, nSnps).T
            G_anchor = G_anchor[:,snp_anchor-startSnpIdx]

            pv, beta = qtl_lr.test_lmm_lr_speed(G_anchor,y_focal, Z=y_orth,Kbg=self.K,Covs=self.Covs, S=self.S, U=self.U)
            #pv, beta = qtl_lr.test_lmm_lr(G_anchor,y_focal, Z=y_orth,Kbg=self.K,Covs=self.Covs)
            self.assoc_updates.add(pv,beta, snp_anchor, focal_gene)


    def save_updates(self,fn):
        self.assoc_updates.save_updates(fn)
        

    def find_incoming_edges(self,t):
        """
        find incoming edges that build a v-structure a-->t<--b with
        a) corr(a,t)
        b) corr(b,t)
        c) ind(a,b)
        d) corr(a,b  t)

        input:
        t   :   index of the gene t
f        """
        # incoming edges are associated with the gene of interest...        
        pv_genes = self.genecorr_reader.getRows([t])[0]
        pv_genes[t] = np.inf # don't count self-correlation in when estimating q-values (always 0)
        qv_genes = np.ones(pv_genes.shape)
        qv_genes[np.isfinite(pv_genes)] = qvalue.estimate(pv_genes[np.isfinite(pv_genes)])
        idx_assoc = qv_genes<self.thresh_corr
        idx_assoc[t] = False
        if not(idx_assoc).any(): return None,None

        # independent of each other
        _idx_assoc = np.nonzero(idx_assoc)[0]
        pv_genes = self.genecorr_reader.getRows(_idx_assoc)[:,idx_assoc]
        idx_vstruct = np.nonzero(idx_assoc)[0]
        pv_genes[~np.isfinite(pv_genes)] = 0
        vstruct     = pv_genes > self.thresh_ind
        idx_ind     = vstruct.any(axis=1)
        idx_vstruct = idx_vstruct[idx_ind]
        vstruct = vstruct[idx_ind][:,idx_ind]
        if not(idx_vstruct).any(): return None,None
         
        # becoming dependent once we condition on the gene under observation
        Yv = self.phenoreader.getRows(idx_vstruct)
        Yt = self.phenoreader.getRows([t])[0]
        _,pv_cond = pcor.pcorParallel(Yv,Yt)

        # take nans into account
        #qv_cond = qvalue.estimate(pv_cond)
        np.fill_diagonal(pv_cond,np.nan)
        qv_cond = np.ones(pv_cond.shape)
        
        # upper triangular matrix
        iu = np.triu_indices(pv_cond.shape[0])
        pv_cond_iu = pv_cond[iu]
        idx_finite = np.isfinite(pv_cond_iu)
        qv_cond_iu = np.ones(pv_cond_iu.shape)
        qv_cond_iu[idx_finite] = qvalue.estimate(pv_cond_iu[idx_finite])
        qv_cond[iu] = qv_cond_iu

        # lower triangular matrix
        il = np.tril_indices(pv_cond.shape[0])
        qv_cond[il] = 0
        qv_cond += qv_cond.T
        
        vstruct*= (qv_cond < self.thresh_corr)
        idx_partcorr = vstruct.any(axis=0)
        if not(idx_partcorr).any(): return None,None
        vstruct = vstruct[idx_partcorr][:,idx_partcorr]
        idx_vstruct = idx_vstruct[idx_partcorr]

        return vstruct, idx_vstruct

    def remove_vstructures_with_related_parents(self,vstruct,idx_orth,idx_snp):
        """
        remove v-structures if the parents are not independent
        from each other

        input:
        vstruct   :   indicator matrix [T x F],

        where T are the orthognal genes and F are the anchor snps
        """
        pv0 = self.assoc0_reader.getRows(idx_snp)
        vstruct *= pv0[:,idx_orth].T  > self.thresh_ind
        vstruct *= ~(self.find_cis_genes(idx_snp[:,np.newaxis], idx_orth)).T
        return vstruct


    def find_anchored_parents(self,vstruct, idx_vstruct):
        """
        subsets vstruct such that only parents with at least one anchored parent remain

        input:
        vstruct     : indicator matrix [TxT]
        idx_vstruct : gene indices,

        where T are all genes building at least one v-structure with the gene t
        """
        idx_orth = np.copy(idx_vstruct)
        idx_anchor = self.anchors.get_gene_idx()
        idx_anchor = np.intersect1d(idx_anchor, idx_vstruct)
        if not(idx_anchor.any()): return None, None,None
            
        idx_y = np.nonzero(idx_vstruct==idx_anchor[:,np.newaxis])[1]
        vstruct = vstruct[:,idx_y]
        idx_x = vstruct.any(axis=1)
        vstruct = vstruct[idx_x]
        idx_orth = idx_orth[idx_x]
        return vstruct, idx_orth, idx_anchor

 

    def remove_snps_related_to_focal_gene(self, t,idx_snp, vstruct):
        """
        remove anchor snps that are physically close to the focal gene

        input:
        t         :   focal gene index
        idx_snp   :   anchor snp indices
        """
        icols = self.find_cis_genes(idx_snp,t)
        if icols.any(): vstruct[:,icols] = False
        return vstruct
                
    def find_vstructures(self,startTraitIdx=0, nTraits=np.inf, max_genes=np.inf):
        """
        returns an iterator over all (snp,orth_gene) pairs where
        snp -> anchor gene -> gene t <- orth gene
        """
        T = self.phenoreader.get_nrows()

        self.vstructures = vstructures.VstructureList()
        for t in range(startTraitIdx, min(startTraitIdx + nTraits,T)):
            print ".... Finding vstructures for gene %d"%t
            for isnps, igenes, ianchors in self.find_vstructures_given_focal_gene(t, max_genes):
                if (isnps is not None) and (igenes is not None):
                    self.vstructures.add(t,isnps,igenes,ianchors)


    def save_vstructures(self,fn):
        self.vstructures.save(fn)


    def load_vstructures(self,fn):
        self.vstructures = vstructures.VstructureFile(fn)


    def find_vstructures_given_focal_gene(self,t, max_genes):
   
        # find incoming edges
        vstruct, idx_vstruct = self.find_incoming_edges(t)
        if vstruct is None:
            yield None, None, None
            return
        
        t1 = time.time()
     
        # find incoming edges that have at least one anchor
        vstruct,idx_orth,idx_anchor = self.find_anchored_parents(vstruct, idx_vstruct)
        if vstruct is None:
            yield None, None, None
            return

        # map anchor genes to anchor snps
        anchor2snp, idx_snp = self.anchors.gene2snp(idx_anchor)
        vstruct = np.dot(vstruct,anchor2snp)
 
        # remove edge if snp is correlated with orthogonal gene
        vstruct = self.remove_vstructures_with_related_parents(vstruct,idx_orth,idx_snp)
 
        # remove anchor snps if close to gene
        vstruct = self.remove_snps_related_to_focal_gene(t,idx_snp, vstruct)
 
        # bundle SNPs that have the same in-coming genes
        np.random.seed(0)
        w_rnd = np.random.randn(vstruct.shape[0])
        _,idx,inv = np.unique(w_rnd.dot(vstruct), return_index=True,return_inverse=True)
        

        for i in np.unique(inv):
            isnps = np.array(idx_snp[inv==i],)
            if vstruct[:,idx[i]].sum()==0: continue
            igenes = np.array(idx_orth[vstruct[:,idx[i]]])

            if len(igenes)>max_genes:
                idx_sorted = np.argsort(self.genecorr_reader.getRows([t])[0,igenes])
                igenes = np.sort(igenes[idx_sorted][:max_genes])

            # map snps back to genes
            ianchors = np.zeros(isnps.shape, dtype=int)
            for j,isnp in enumerate(isnps):
                try:
                    ianchors[j] = idx_anchor[anchor2snp[:,idx_snp==isnp][:,0]]
                except:
                    ianchors[j] = np.inf
         
            yield (isnps,igenes,ianchors)
            
 
