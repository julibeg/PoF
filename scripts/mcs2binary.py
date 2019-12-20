#!/usr/bin/env python
import sys


mcsfname = sys.argv[1]
rxnsfname = sys.argv[2]
outfname = sys.argv[3]

with open(rxnsfname, 'r') as f:
	rxns = f.read().strip().split()
	rxns = [r.replace('"', '') for r in rxns]

with open(mcsfname, 'r') as mcsfile, open(outfname, 'w') as outfile:
	for line in mcsfile:
		arr = ['0'] * len(rxns)
		sep = ' ' if ' ' in line else ','
		mcs = line.strip().split(sep)
		for rxn in mcs:
			arr[rxns.index(rxn)] = '1'
		outfile.write(''.join(arr) + '\n')


