#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/EXP00/TIDY/ARCHIVE/'

# --------------------------------------------------------------------------- #

field01 = np.squeeze(np.load(''.join([homedir, nemodir,'POST/',
                                      'vosal_ANN_1997-1997.npy'])))

field02 = np.squeeze(np.load(''.join([homedir, nemodir,'POST/',
                                      'vosal_ANN_1998-1998.npy'])))

field03 = np.squeeze(np.load(''.join([homedir, nemodir,'POST/',
                                      'vosal_ANN_1999-1999.npy'])))

avefield = (field01 + field02 + field03) / 3.0

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'POST/vosal_ANN_1997-1999.npy']), avefield)

# --------------------------------------------------------------------------- #
