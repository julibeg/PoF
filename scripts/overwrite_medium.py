#!/usr/bin/env python

import sys
import cobra
import os
import util
import pickle


infname, ext = os.path.splitext(sys.argv[1])

m = util.read_model(sys.argv[1])
medium = pickle.load(open(sys.argv[2], 'rb'))

m.medium = medium

util.write_model(m, f'minimal_medium_{infname}{ext}')
if m.slim_optimize() > m.tolerance:
	print(f'overwritten medium for {infname} - still viable')
else:
	print(f'overwritten medium for {infname} - no longer viable')
	
