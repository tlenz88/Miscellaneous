#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

arr = np.genfromtxt(sys.argv[1], delimiter = '\t')

# Extract scalar values
scalars = [item[2] for item in arr]

# Calculate median of scalar values
med_scalar = np.median(scalars)

# Calculate absolute deviations from the median
abs_dev = [np.abs(item[2] - med_scalar) for item in arr]

# Calculate MAD (Median absolute deviation)
mad = np.median(abs_dev)

# Set a threshold for outlier detection (e.g., 2 or 3 times MAD)
threshold_mad = 10

# Identify outliers
ol_mad_rem = [arr[i] for i in range(len(arr)) 
              if abs_dev[i] <= threshold_mad * mad]
ol_mad = [arr[i] for i in range(len(arr)) 
          if abs_dev[i] > threshold_mad * mad]

df_ol_rem = pd.DataFrame(ol_mad_rem)
df_ol = pd.DataFrame(ol_mad)
df_ol_rem[0] = df_ol_rem[0].astype(int)
df_ol_rem[1] = df_ol_rem[1].astype(int)
df_ol[0] = df_ol[0].astype(int)
df_ol[1] = df_ol[1].astype(int)
df_ol_rem.to_csv(sys.argv[1].replace('.matrix', '_outliers_removed.matrix'), 
          sep = '\t', header = False, index = False)
df_ol.to_csv(sys.argv[1].replace('.matrix', '_outliers.matrix'), 
          sep = '\t', header = False, index = False)