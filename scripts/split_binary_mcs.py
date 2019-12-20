#!/usr/bin/env python
import numpy as np
import sys


bin_mcs = sys.argv[1]
out_prefix = sys.argv[2]

mcs = np.genfromtxt(bin_mcs, delimiter=1, dtype=int)
ds = mcs.sum(1)

# wipe files if they already exist
for d in ds:
	open(out_prefix + str(d), 'w').close()

# write binary MCS to the respective file
for row in mcs:
	with open(out_prefix + str(row.sum()), 'a') as f:
		f.write(''.join(row.astype(str)) + '\n')

