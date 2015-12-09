
import reader
from reader import Reader


class PhenoReader(Reader):
    """
    Interface for phenotype reader. Must be implemented from all inherited classes.
    """

    def getGeneStart(self):
        """
        returns the gene start
        """
        raise Exception('must be implemented from inherited class.')

    def getGeneEnd(self):
        """
        returns the gene ened
        """
        raise Exception('must be implemented from inherited class.')


    def getGeneChrom(self):
        """
        returns the chromosomal information
        """
        raise Exception('must be implemented from inherited class.')

    def getGeneIds(self):
        """
        returns the gene identifiers
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




      





