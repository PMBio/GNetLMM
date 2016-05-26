import pdb
import numpy as np

class Anchors:

    def __init__(self,F,T,pv=None,snp_ids=None,gene_ids=None,igene=None,isnp=None):
        self.F = F
        self.T = T
        
        self.pv = pv
        self.snp_ids = snp_ids
        self.gene_ids = gene_ids
        self.igene = igene
        self.isnp = isnp


    def save(self,fn):
        header = 'gene_ids\tsnp_ids\tgene_index\tsnp_index\tpv\n'
      
        fout = open(fn,'w')
        fout.write(header)
        for i in range(len(self.gene_ids)):
            line = '%s\t%s\t%d\t%d\t%.4e\n'%(self.gene_ids[i],self.snp_ids[i],self.igene[i],self.isnp[i],self.pv[i])
            fout.write(line)
        fout.close()
    

    def load(self,fn):
        M = np.loadtxt(fn,delimiter='\t',dtype=str,skiprows=1)
        
        if M.shape[0]==0:
            self.gene_ids = []
            self.snp_ids = []
            self.igene = []
            self.isnp = []
            self.pv = []

        else:
            if M.ndim==1: M=M[np.newaxis,:]
            self.gene_ids = M[:,0]
            self.snp_ids = M[:,1]
            self.igene = np.array(M[:,2],dtype=int)
            self.isnp = np.array(M[:,3],dtype=int)
            self.pv   = np.array(M[:,4],dtype=np.float)
        
      
    def get_gene_idx(self):
        return np.unique(self.igene)


    def gene2snp(self,igene):
        """
        creates mapping from genes to snps
        """
        N = len(self.gene_ids)

        isnp  = np.unique(self.isnp)
        T = len(igene)
        F = len(isnp)

        M = np.zeros((T,F),dtype=bool)
        
        for n in range(N):
            irow = self.igene[n]==igene
            if irow.any():
                icol  = self.isnp[n]==isnp
                M[irow,icol]=1
            
        icol = M.any(axis=0)
        M = M[:,icol]
        isnp = isnp[icol]
        
        return M,isnp
