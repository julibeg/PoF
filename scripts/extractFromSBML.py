#!/usr/bin/env python
import numpy as np
import cobra
import sys
import os


fname = sys.argv[1]
model_name, extension = os.path.splitext(fname)
if extension == '.json':
	m = cobra.io.load_json_model(fname)
elif extension == '.xml' or extension == '.sbml':
	m = cobra.io.read_sbml_model(fname)
else:
	print(f'ERROR: input file ({fname}) missing matching extension (.json/.xml/.sbml)')
	sys.exit(1)

#write sfile
np.savetxt(f'{model_name}.sfile', 
		cobra.util.array.create_stoichiometric_matrix(m),
		delimiter='\t', fmt='%g')

#write mfile
with open(f'{model_name}.mfile', 'w') as f:
	f.write(' '.join([met.id for met in m.metabolites]))

#write rfile
with open(f'{model_name}.rfile', 'w') as f:
	f.write(' '.join([rxn.id for rxn in m.reactions]))

#write rvfile
with open(f'{model_name}.rvfile', 'w') as f:
	f.write(' '.join(['1' if rxn.reversibility \
			else '0' for rxn in m.reactions]))


