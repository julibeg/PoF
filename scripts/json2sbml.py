#!/usr/bin/env python

import sys
import cobra
import os


fname, extension = os.path.splitext(sys.argv[1])

m = cobra.io.load_json_model(sys.argv[1])

cobra.io.write_sbml_model(m, f'{fname}.sbml')
