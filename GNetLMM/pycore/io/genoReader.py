import numpy as np


class GenoReader():
    """
    Interface for genotype reader. Must be implemented from all inherited classes.
    """

    def getSnpPos(self):
        """
        returns the snp position (in basepairs)
        """
        raise Exception('must be implemented from inherited class.')

    def getSnpChrom(self):
        """
        returns the chromosomal information
        """
        raise Exception('must be implemented from inherited class.')

    def getSnpIds(self):
        """
        returns the snp identifiers
        """ 
        raise Exception('must be implemented from inherited class.')

    def get_nrows(self):
        """
        returns the number of SNPs/ rows
        """
        raise Exception('must be implemented from inherited class.')

    def get_ncols(self):
        """
        return the number of samples/ columns
        """
        raise Exception('must be implemented from inherited class.')


    def loadSnpBlock(self,start = 0, nSNPs = np.inf):
      """
      start           : index of the first SNP to be loaded 
      nSNPs           : load nSNPs (default np.inf, meaning all)
      """
      raise Exception('must be implemented from inherited class.')

      
        
