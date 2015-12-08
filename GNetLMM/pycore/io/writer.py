import pdb
import numpy as np
import csv

class Writer:
    """
    basic class for writing out experiments
    """
    def __init__(self, basefile):
        """
        constructor

        input:
        basefile   : name of basefile
        """
        self.basefile = basefile


    def _writeInfo(self,data,which):
        """
        writing out column (which=cols) or row (which=rows) information

        input:
        data   :   dictionary, containing information
        which  :   specifies file ending
        """
        with open(self.basefile + '.' + which, "wb") as outfile:
            csv_writer = csv.writer(outfile,delimiter=' ')
            csv_writer.writerow(data.keys())
            csv_writer.writerows(zip(*data.values()))
            
    def writeColumnInfo(self,data):
        """
        writing out column information
        """
        self._writeInfo(data, 'cols')
       
    def writeRowInfo(self,data):
        """
        writing out row information
        """
        self._writeInfo(data, 'rows')

    def writeMatrix(self,M,**kwargs):
        """
        writing out matrix

        input:
        M   :   data matrix
        """
        np.savetxt(self.basefile + '.matrix', M,**kwargs)
