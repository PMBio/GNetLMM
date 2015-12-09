import reader
import phenoReader

class PhenoReaderMatrix(reader.MatrixReader, phenoReader.PhenoReader):

    def __init__(self, M, gene_start, gene_stop, chrom,ids):
        """
        constructor

        input:
        M       :   phenotype matrix [TxN]
        pos     :   positions [T]
        chrom   :   chromosomal positions [T]
        ids     :   snp identifiers [T]
        """
        self.M = M
        self.gene_start = gene_start
        self.gene_stop  = gene_stop
        self.chrom = chrom
        self.ids = ids

    def getGeneStart(self):
        return self.gene_start

    def getGeneEnd(self):
        return self.gene_stop
        
    def getGeneChrom(self):
        return self.chrom

    def getGeneIds(self):
        return self.ids

    def get_nrows(self):
        return self.M.shape[0]

    def get_ncols(self):
        return self.M.shape[1]
 
