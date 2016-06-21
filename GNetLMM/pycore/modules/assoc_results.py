import numpy as np
import pdb 


class AssocResultsList:

    def __init__(self):
        self.pv = []
        self.beta = []
        self.focal_gene = []
        self.snp_anchor  = []

        self.var_snps = []
        self.var_covs = []
        self.var_genes = []

    def add(self, pv, beta, var_snps, var_covs, var_genes,snp_anchor, focal_gene):

        n = len(pv)
        focal_gene = np.repeat(focal_gene, n)

        self.pv.append(pv)
        self.beta.append(beta)
        self.snp_anchor.append(snp_anchor)
        self.focal_gene.append(focal_gene)
        self.var_snps.append(var_snps)
        self.var_covs.append(var_covs)
        self.var_genes.append(var_genes)
        

    def save_updates(self,fn):
        if len(self.pv)>0:
            pv = np.concatenate(self.pv)[:,np.newaxis]
            beta = np.concatenate(self.beta)[:,np.newaxis]
            focal_gene = np.concatenate(self.focal_gene)[:,np.newaxis]
            snp_anchor = np.concatenate(self.snp_anchor)[:,np.newaxis]
            var_snps = np.concatenate(self.var_snps)[:,np.newaxis]
            var_covs = np.concatenate(self.var_covs)[:,np.newaxis]
            var_genes = np.concatenate(self.var_genes)[:,np.newaxis]
            matrix = np.hstack((focal_gene,snp_anchor,pv,beta, var_snps, var_covs, var_genes))
            np.savetxt(fn+'.csv', matrix, fmt='%d\t%d\t%.4e\t%.4f\t%.4e\t%.4e\t%.4e')
        else:
            f = open(fn + ".csv","w")
            f.close()
  

    def load_csv(self,fn):
        M = np.loadtxt(fn)
        self.var_snps = None 
        self.var_covs = None 
        self.var_genes = None

        if M.shape[0]==0:
            self.focal_gene = np.array([])
            self.snp_anchor =  np.array([])
            self.pv =  np.array([])
            self.beta =  np.array([])
        else:
            if M.ndim==1: M = M[np.newaxis,:]
            self.focal_gene = np.array(M[:,0], dtype=int)
            self.snp_anchor = np.array(M[:,1], dtype=int)
            self.pv = M[:,2]
            self.beta = M[:,3]

            if M.shape>4:
                self.var_snps = M[:,4]
                self.var_covs = M[:,5]
                self.var_genes = M[:,6]


    
    def save_matrix(self,fn0,fn_out):
        # re-arrange such that SNPs are ordered
        idx = np.argsort(self.snp_anchor)
        self.pv = self.pv[idx]
        self.beta = self.beta[idx]
        self.focal_gene = self.focal_gene[idx]
        self.snp_anchor = self.snp_anchor[idx]
        N = len(idx)

        self._save_matrix(self.pv, fn0+'.pv.matrix', fn_out+'.pv.matrix','%.4e')
        self._save_matrix(self.beta, fn0+'.beta.matrix', fn_out+'.beta.matrix','%.4f')


    def _save_matrix(self,x, fn_in, fn_out, fmt):

   
        fin = open(fn_in,'r')
        fout = open(fn_out, 'w')
        N = len(x)
        
        i = 0 # line index of assoc0 file
        j = 0 # line index of sorted list

        

        line = fin.readline()
        while line:
            if j<N and self.snp_anchor[j]==i:
            
                arr = line.split(' ')
                while j<N and self.snp_anchor[j]==i:

            
                    arr[self.focal_gene[j]] = fmt%x[j]
                    j += 1

                
                line = ' '.join(arr)
                if not(line.endswith('\n')): line += '\n'
                    
            fout.write(line)                
            line = fin.readline()
            i += 1

        fin.close()
        fout.close()
        



