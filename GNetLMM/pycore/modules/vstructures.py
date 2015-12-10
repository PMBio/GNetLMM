import pdb
import numpy as np
import csv

class VstructureList:
    """
    container for Vstructures
    """

    def __init__(self):
        self.focal_gene = []
        self.snp_anchor = []
        self.orth_gene  = []


    def add(self, focal_gene, snp_anchor, orth_gene):
        """
        adding a v-structure 
        """
        self.focal_gene.append(focal_gene)
        self.snp_anchor.append(snp_anchor)
        self.orth_gene.append(orth_gene)

    def save(self,fn):
        """
        writing v-structures to file
        """
        M = len(self.orth_gene)

        csvfile = fout = open(fn,'w')
        csvwriter = csv.writer(csvfile, delimiter='\t')

        for m in range(M):

            snp_anchor = ','.join([str(x) for x in self.snp_anchor[m]])
            orth_gene  = ','.join([str(x) for x in self.orth_gene[m]])
            csvwriter.writerow([self.focal_gene[m], snp_anchor, orth_gene])

        fout.close()


    def iterator(self):
        M = len(self.orth_gene)

        for m in range(M):
            yield self.focal_gene[m], self.snp_anchor[m], self.orth_gene[m]


class VstructureFile:

    def __init__(self,fn):
        self.fn = fn

    def iterator(self):
        """
        iterating over v-structures
        """
        fin = open(self.fn,'r')

        line = fin.readline()

        while line:
            line = line.rstrip()
            arr = line.split('\t')
            focal_gene = np.array([int(arr[0])])
            snp_anchor = np.array(arr[1].split(','), dtype=int)
            orth_gene  = np.array(arr[2].split(','), dtype=int)
            yield focal_gene, snp_anchor, orth_gene
            line = fin.readline()
            
        fin.close()
