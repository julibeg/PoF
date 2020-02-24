#!/usr/bin/env python

import sys
import os 
import cobra


fname, extension = os.path.splitext(sys.argv[1])

if extension == '.json':
	m = cobra.io.load_json_model(sys.argv[1])
else:
	m = cobra.io.read_sbml_model(sys.argv[1])

for key, item in m.medium.items():
	print(fname, key, item)
print()


