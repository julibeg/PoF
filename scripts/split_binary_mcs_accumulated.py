#!/usr/bin/env python
import numpy as np
from subprocess import Popen, PIPE
import sys


mcsfile = sys.argv[1]
output_prefix = sys.argv[2]

a = np.genfromtxt(mcsfile, delimiter=1, dtype=int)

s_old = a[0].sum()
for i, s in enumerate(a.sum(1)):
	if s > s_old:
		with open(f'{output_prefix}{s_old}', 'w') as f:
			Popen(['head', f'-n{i}', mcsfile], stdout=f)
		print(f'd<={s_old}: done')
		s_old = s
	elif s < s_old:
		print('Decreasing cardinalities encountered. Is the input a valid binary MCS file and sorted properly?')
		break
# else of for loop
else:
	with open(f'{output_prefix}{s_old}', 'w') as f:
		# this essentially just copies the input file
		Popen(['head', f'-n{i+1}', mcsfile], stdout=f)
	print(f'd<={s}: done')
