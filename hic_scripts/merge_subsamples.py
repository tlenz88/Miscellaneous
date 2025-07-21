import sys
import pandas as pd
from operator import itemgetter

subsample_dict = {}
for input_file in sys.argv[1:]:
	print(input_file)
	with open(input_file, "r") as f:
		for line in f:
			if line[0] not in subsample_dict:
				subsample_dict[line] = 1
			else:
				subsample_dict[line] += 1
output_dict = dict(sorted(subsample_dict.items(), key = itemgetter(1), reverse = True)[:36900589])
with open("/mnt/f/Deitsch_collab/HiC/output_files/A3/A3_subsample.allValidPairs", "w") as outfile:
	for k, v in output_dict.items():
		outfile.write(k)
