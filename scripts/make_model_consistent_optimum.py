#!/usr/bin/env python

import sys
import cobra
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra.flux_analysis import flux_variability_analysis as FVA
import os


infname = sys.argv[1]
fva_proc = int(sys.argv[2]) if len(sys.argv) > 2 else 1
model_name, ext = os.path.splitext(infname)

m = cobra.io.read_sbml_model(infname)

fba_orig = m.slim_optimize()
rxns_orig = len(m.reactions)
mets_orig = len(m.metabolites)
tol_orig = m.tolerance

fva = FVA(m, processes=fva_proc, fraction_of_optimum=0.99)
blocked_rxns = list(fva.index[fva.abs().sum(1) < tol_orig])

m.remove_reactions(blocked_rxns)
m, _ = cobra.manipulation.delete.prune_unused_metabolites(m)

fba_cons = m.slim_optimize()
rxns_cons = len(m.reactions)
mets_cons = len(m.metabolites)

if abs(fba_orig - fba_cons) > tol_orig:
	print(f'ERROR: {infname}: difference in FBA objective is too large')
	print(f'FBA orig: {fba_orig}\tFBA cons: {fba_cons}')
	sys.exit(1)

cobra.io.write_sbml_model(m, f'consistent_{model_name}{ext}')
print(f'{model_name}:\n{rxns_orig} rxns, {mets_orig} mets, obj: {fba_orig} --> {rxns_cons} rxns, {mets_cons} mets, obj: {fba_cons}\n')
