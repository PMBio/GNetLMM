import numpy as np
import genoReader 
import GNetLMM.pycore.utils.utils as utils


class GenoReaderMatrix(genoReader.GenoReader):
    """
    geno reader using matrix/hdf file
    """
    
    def __init__(self, M, pos, chrom, ids):
        """
        constructor

        input:
        M       :   genotype matrix [FxN]
        pos     :   positions [F]
        chrom   :   chromosomal positions [F]
        ids     :   snp identifiers [TF]
        """
        self.M = M
        self.pos = pos
        self.chrom = chrom
        self.ids = ids
        

    def getSnpPos(self):
        """
        returns the position of the SNPs (in basepairs)
        """
        return self.pos

    def getSnpChrom(self):
        """
        returns the chromosomal information
        """
        return self.chrom

    def getSnpIds(self):
        """
        returns the snp identifiers
        """
        return self.ids

    def get_nrows(self):
        """
        returns the number of SNPs/ rows
        """
        return self.M.shape[0]

    def get_ncols(self):
        """
        returns the number of samples/ columns
        """
        return self.M.shape[1]


    def loadSnpBlock(self,start = 0, nSNPs = np.inf):
      """
      start           : index of the first SNP to be loaded 
      nSNPs           : load nSNPs (default np.inf, meaning all)
      """
      return self.M[start:start+nSNPs,:]


    def compute_Kbg_excluded_chrom(self,chrom_i):
        """
        computes the background covariance matrix from all SNPs that are not on
        chromosome chrom_i
        """
        idx_bg = self.chrom!=chrom_i
        M = self.M[idx_bg]
        Kbg = utils.computeLinearKernel(M.T)
        return Kbg


    def get_chrom_idx(self,chrom_i):
        """
        returns a boolean vector indicating which SNPs are on chromosome chrom_i
        """
        idx = np.nonzero(self.chrom==chrom_i)[0]
        start = idx[0]
        nSnps = len(idx)
        return start, nSnps
        
        return self.chrom==chrom_i




        
