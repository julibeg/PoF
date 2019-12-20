#!/usr/bin/env python

import cobra
from cobra.flux_analysis.variability import find_blocked_reactions
import sys


model_fname = sys.argv[1]

model_orig = cobra.io.read_sbml_model(model_fname)
model = model_orig.copy()

model_outfname = sys.argv[2] if len(sys.argv) > 2 else f'{model.id}_all_C_sources.xml'

medium_orig = model.medium

exchanges = []
for r in model.exchanges:
	for met in r.metabolites.keys():
		if 'C' in met.formula:
			exchanges.append(r)
			break

c_sources = []
for ex in exchanges:
	medium = medium_orig.copy()
	medium['EX_glc__D_e'] = 0 
	medium[ex.id] = 10
	with model as m:
		m.medium = medium
		if m.slim_optimize() > m.tolerance:
			c_sources.append(ex)

medium = medium_orig.copy()
for c_source in c_sources:
	medium[c_source.id] = 10

model_orig.medium = medium

cobra.io.write_sbml_model(model_orig, model_outfname)
print(f'{model_orig.id}: {len(c_sources)} C-sources')

