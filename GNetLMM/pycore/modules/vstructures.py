import pdb
import numpy as np
import csv

class VstructureList:
    """
    container for Vstructures
    """

    def __init__(self):
        self.focal_gene = []
        self.anchor_snp = []
        self.orth_gene  = []
        self.anchor_gene = []


    def add(self, focal_gene, anchor_snp, orth_gene, anchor_gene):
        """
        adding a v-structure 
        """
        self.focal_gene.append(focal_gene)
        self.anchor_snp.append(anchor_snp)
        self.orth_gene.append(orth_gene)
        self.anchor_gene.append(anchor_gene)

    def save(self,fn):
        """
        writing v-structures to file
        """
        M = len(self.orth_gene)

        csvfile = fout = open(fn,'w')
        csvwriter = csv.writer(csvfile, delimiter='\t')

        for m in range(M):

            anchor_snp  = ','.join([str(x) for x in self.anchor_snp[m]])
            anchor_gene = ','.join([str(x) for x in self.anchor_gene[m]])
            orth_gene   = ','.join([str(x) for x in self.orth_gene[m]])
            csvwriter.writerow([self.focal_gene[m], anchor_snp, orth_gene, anchor_gene])

        fout.close()


    def iterator(self, full=False):
        M = len(self.orth_gene)

        for m in range(M):
            if full:
                yield self.focal_gene[m], self.snp_anchor[m], self.orth_gene[m], self.anchor_gene[m]
            else:
                yield self.focal_gene[m], self.snp_anchor[m], self.orth_gene[m]


class VstructureFile:

    def __init__(self,fn):
        self.fn = fn

    def iterator(self, full=False):
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
            
            if full:
                anchor_gene = None
                if len(arr)==4:
                    anchor_gene = np.array(arr[3].split(','), dtype=int)
                yield focal_gene, snp_anchor, orth_gene, anchor_gene
            else:
                yield focal_gene, snp_anchor, orth_gene
            line = fin.readline()
            
        fin.close()
