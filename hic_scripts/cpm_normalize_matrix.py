#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

for i in range(1,len(sys.argv)):
    mat = pd.read_csv(sys.argv[i], sep='\t', header=None)
    mmr = np.sum(mat[2])
    mat[2] = mat[2] / (mmr / 1e6)
    mat.to_csv(sys.argv[i].replace('.matrix', '_cpm.matrix'), sep = "\t", 
               index = False, header = False)
