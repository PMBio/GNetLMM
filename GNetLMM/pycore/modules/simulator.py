import scipy as sp
import scipy.stats 
import numpy as np

import GNetLMM.pycore.io.bedReader as bedReader

def genBinormal(dim1,dim2,percSign=0.5,std=1e-1):
    
    rv = (2*(sp.rand(dim1,dim2)>percSign)-1)+std*sp.randn(dim1,dim2)
    return rv


class Simulator:
    def __init__(self,bfile,dist_min=100):
        """
        initializes object

        input:
        bfile      :   basefilename of plink file
        dist_min   :   minimal distance between two causal SNPs
        """
        self.dist_min = dist_min
        self.bfile = bfile

        self.genoreader = bedReader.BedReader(self.bfile)

        self.chrom = self.genoreader.getSnpChrom()
        self.pos   = self.genoreader.getSnpPos()
        self.rs    = self.genoreader.getSnpIds()
      
   
        self.N = self.genoreader.fam.shape[0]
        self.F = self.chrom.shape[0]

    
    def simulateConfoundingFactors(self,nConfounder=1, standardized=True):
        """
        simulate confounding factors

        input:
        nConfounder   :   number of confounding factors
        standardized  :   if True (default), hidden confounders are standardized to zero mean and unit variance
        """
        H = sp.random.randn(self.N,nConfounder)
        if standardized:
            H-= H.mean(axis=0)
            H/= H.std(axis=0)
        return H


    def standardizeGenes(self,Y):
        """
        normalize genes to zero mean and unit variance

        input:
        Y   :   phenotypic matrix [NxT]
        """
        T = Y.shape[0]
        for irow in range(T):
            Y[irow] -= Y[irow].mean()
            Y[irow] /= Y[irow].std()
        return Y

    def simulateSNPindices(self,T):
        """
        select causal SNP indices

        input:
        T   :   number of phenotypes
        """
        iSnp = sp.random.randint(0,self.F,T)
        for irow in range(1,T):
            # take only snps that are not in close proximity to an already selected causal SNP
            dist = sp.absolute(self.pos[iSnp[irow]]-self.pos[iSnp[:irow]])
            dist[self.chrom[:irow]!=self.chrom[iSnp[irow]]] = self.dist_min
            while np.min(dist)<self.dist_min:
                iSnp[irow] = sp.random.randint(0,self.F,1)
        return iSnp


    def simConfounderNetwork(self, T, confPerGene, nConfounder):
        """
        simulate confounding network

        input:
        T             :   number of traits
        confPerGene   :    expected confounders per gene
        nConfounder   :    number of confounders
        """
        pConf = 0
        if nConfounder!=0: pConf = 1.*confPerGene/nConfounder
        Aconf = sp.random.rand(T,nConfounder)
        Aconf[Aconf<pConf] = 1
        Aconf[Aconf<1]     = 0
        Aconf = sp.array(Aconf,dtype=bool)
        return Aconf
    
    def simulateGenes(self,T=100,varSnp=0.1,varNetwork=0.8,alpha=0.8,nConfounder=3,confPerGene=1,**kwargs):
        """
        simulate genes 

        Y = varSnp * Ysnp + varNetwork * Ynetwork + (1-varSnp - varNetwork)*Ynoise

        Ynetwork = alpha*Ygene + (1-alpha)*Yconf

        input:
        T             :    number of traits (default: 100)
        varNetwork    :    variance explained by the network (default: 0.8)
        varSnp        :    variance explained by the snp (default: 0.1)
        alpha         :    ratio between variance explained by the gene-gene network and the gene-confounder network
        nConfounder   :    number of confounders
        confPerGene   :    expected number of confounders per gene
        """
        # simulate adjacency matrix for gene-gene network
        Agene = self.simulateAdjacencyMatrix(T,**kwargs)

        # simulate hidden confounding factors
        H  = self.simulateConfoundingFactors(nConfounder)

        # simulated adjacency matrix for confounder-gene network
        Aconf = self.simConfounderNetwork(T, confPerGene, nConfounder)
        

        varGene = alpha * varNetwork * Agene.any(axis=1)
        varConf = (1-alpha) * varNetwork * Aconf.any(axis=1)
        varNoise= 1 - varSnp - varGene - varConf

        # simulating causal SNPs
        iSnp  = self.simulateSNPindices(T)

        # simulating gene-gene network
        Y      = sp.zeros((T,self.N))
        Ysnp   = sp.zeros((T,self.N))
        Ynoise = sp.zeros((T,self.N))
        Ygene  = sp.zeros((T,self.N))
        for t in range(T):

            G = self.genoreader.loadSnpBlock(start=iSnp[t], nSNPs=1)[0]
            G = (G-G.mean())/G.std()
            
            Ysnp[t]   = sp.sqrt(varSnp)*G
            Ynoise[t] = sp.random.randn(self.N)
            Ynoise[t]*= sp.sqrt(varNoise[t])/Ynoise[t].std()
            if Agene[t].sum()>0:
                nGenes = Agene[t].sum()
                Y_in  = Y[Agene[t]].T
                Y_in -= Y_in.mean()
                Y_in /= Y_in.std()
                wGene = genBinormal(nGenes,1)[:,0] 
                Ygene[t] = sp.dot(Y_in,wGene)
                Ygene[t]*= sp.sqrt(varGene[t])/Ygene[t].std()
            Y[t] = Ysnp[t] + Ygene[t] + Ynoise[t]

        # simulating confounder-gene network
        Yconf  = sp.zeros((T,self.N))
        for t in range(T):
            if Aconf[t].sum()>0:
                nConf = Aconf[t].sum()
                wConf = genBinormal(nConf,1)[:,0]
                Yconf[t] = sp.dot(H[:,Aconf[t]],wConf)
                Yconf[t]*= sp.sqrt(varConf[t])/Yconf[t].std()

        # plugging phenotype together
        Y += Yconf
        Y = self.standardizeGenes(Y)
        
        # set genes
        RV = {}
        RV['rs'] = self.rs[iSnp]
        RV['gene_ids'] = sp.array(['id_%d'%x for x in range(T)])
        RV['gene_chrom'] = self.chrom[iSnp]
        RV['gene_start'] = self.pos[iSnp]
        RV['H'] = H
        RV['Agene'] = Agene
        RV['Aconf'] = Aconf
        RV['Y'] = Y
        return RV


class SparseSimulator(Simulator):
    
    def simulateAdjacencyMatrix(self,T,expN=5,**kwargs):
        """
        simulate sparse adjacency matrix

        input:
        T      :   number of traits
        expN   :   expected number of in-coming and out-going edges (default: 5).
        """
        s = 1.*expN/(T-1) # sparsity level
        A = sp.zeros((T,T),dtype=bool)

        for irow in range(1,T):
            a_row = sp.random.rand(irow)<s
            icol  = sp.nonzero(a_row)[0][:5]
            A[irow,icol] = True

        return A

    
class HotSpotSimulator(Simulator):
    
    def simulateAdjacencyMatrix(self,T,minProbReg=0.1,maxProbReg=0.5,expN=5,**kwargs):
        """
        simulate adjaency matrix with hubs

        input:
        T            :   number of traits
        minProbReg   :   min fraction of genes that are regugulated by an hotspot
        maxProbReg   :   max fraction of genes that are regugulated by an hotspot
        expN         :   expected number of in-coming and out-going edges (default: 5).
        """
        nHotspot = int(sp.ceil(expN / (maxProbReg + minProbReg)))
        A = sp.zeros((T,T),dtype=bool)
        probHotspot = minProbReg + (maxProbReg-minProbReg)*sp.rand(nHotspot)
        for ihotspot in range(nHotspot):
            irow = sp.random.rand(T)
            irow[irow<probHotspot[ihotspot]] = 0
            irow[irow>probHotspot[ihotspot]]  = 1
            irow = sp.array(1-irow,dtype=bool)
            irow[:nHotspot] = False
            A[irow,ihotspot] = True

        self.nHotspot = nHotspot

        return A


 
