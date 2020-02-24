#!/usr/bin/env python

import sys

sys.path.append('/home/jlibiseller/git/mcs_stuff/code/pof_python_implementation')

from PoF9 import PoF_calculator
import numpy as np


mcs_file = sys.argv[1]
d0 = int(sys.argv[2])
if len(sys.argv) > 3:
	threads = int(sys.argv[3])
else:
	threads = 1

mcs = np.genfromtxt(mcs_file, dtype=int, delimiter=1)

calc = PoF_calculator(mcs)

weights, Fs, pofs = calc.get_all_PoFs(max_d=d0, ncores=threads)

for w, f, p in zip(weights, Fs, pofs):
	print(f'{p}\t{f}\t{w}')


