import os
import cobra
import sys


def read_model(infname):
	model_name, ext = os.path.splitext(infname)
	if ext == '.json':
		m = cobra.io.load_json_model(infname)
	elif ext == '.xml' or ext == '.sbml':
		m = cobra.io.read_sbml_model(infname)
	else:
		print(f'ERROR: input file {infname} lacking viable extension (.xml/.sbml/.json)')
		sys.exit(1)
	return m

def write_model(model, fname):
	ext = os.path.splitext(fname)[1]
	if ext == '.json':
		cobra.io.save_json_model(model, fname)
	else:
		cobra.io.write_sbml_model(model, fname)
