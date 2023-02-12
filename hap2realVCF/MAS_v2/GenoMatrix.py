import os
import sys
import pandas as pd
import numpy as np

import numpy as np
from scipy.spatial.distance import pdist

def dfun(u, v):
    return np.sqrt(((u-v)**2).sum())


dm = pdist(X, dfun)
class GenotypeDF:
    def __init__(self, hap_file):
        self.hap_file = hap_file
        GenotypeDF.samples,GenotypeDF.alleles,GenotypeDF.gt, GenotypeDF.depth=self.__parse()
     
    def __parse(self):
        
        read.csv(self.hap_file):

	return samples, allels, gt, depth

    def dfun(u,v):
        d=len(set(u).intersection(set(v)))/(set(u).union(set(v)) 
        return d

    def __ld_pair(x,y)

    def ld_multiple( )

    def missing_ind()
    
    def 
    
 
