import reader
import phenoReader

class PhenoReaderFile(reader.FileReader, phenoReader.PhenoReader):

    def getGeneStart(self):
        return self.row_info['gene_start']

    def getGeneEnd(self):
        if 'gene_end' in self.row_info:
            gene_end =  self.row_info['gene_end']
        elif 'gene_stop' in self.row_info:
            gene_end =  self.row_info['gene_stop']
        else:
            gene_end =  self.row_info['gene_start']
        return gene_end
        
    def getGeneChrom(self):
        gene_chrom = self.row_info['gene_chrom']
        return gene_chrom

    def getGeneIds(self):
        if 'gene_ids' in self.row_info:
            return self.row_info['gene_ids']

    def get_nrows(self):
        keys = self.row_info.keys()
        return len(self.row_info[keys[0]])

    def get_ncols(self):
        keys = self.col_info.keys()
        return len(self.col_info[keys[0]])
