#!/usr/bin/env python

import sys
import cobra
from cobra.flux_analysis.variability import find_blocked_reactions
import os


infname = sys.argv[1]
model_name, ext = os.path.splitext(infname)

m = cobra.io.read_sbml_model(infname)

fba_orig = m.slim_optimize()
rxns_orig = len(m.reactions)
mets_orig = len(m.metabolites)
tol_orig = m.tolerance

m.remove_reactions(find_blocked_reactions(m))
m, _ = cobra.manipulation.delete.prune_unused_metabolites(m)

fba_cons = m.slim_optimize()
rxns_cons = len(m.reactions)
mets_cons = len(m.metabolites)

if abs(fba_orig - fba_cons) > (fba_orig * tol_orig):
	print(f'ERROR: {infname}: difference in FBA objective is too large')
	sys.exit(1)

cobra.io.write_sbml_model(m, f'consistent_{model_name}{ext}')
print(f'{model_name}:\n{rxns_orig} rxns, {mets_orig} mets, obj: {fba_orig} --> {rxns_cons} rxns, {mets_cons} mets, obj: {fba_cons}\n')
