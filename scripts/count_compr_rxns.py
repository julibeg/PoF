#!/usr/bin/env python

import numpy as np
import sys


infname = sys.argv[1]
outfname = sys.argv[2]

rx = np.genfromtxt(infname, dtype=str)
rx_cnt = np.char.count(rx, '%') + 1

np.savetxt(outfname, rx_cnt.reshape(1, -1), fmt='%g', delimiter=' ')
