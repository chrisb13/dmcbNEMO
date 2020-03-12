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

# Specify the field to average
fieldname1 = 'vomece3v'
fieldname2 = 'e3v'

# Specify the years to do the averaging over and the season to average.
init_year = 1997
num_years = 3

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant mask.
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to make the require seasonal
# average.
field = nemo.calc_annual_ratio(fieldname1, fieldname2, 
                               vmask, 'V', 1.0, homedir, nemodir, 'd05/V/',
                               prefix='ORCH0083-LIM3_', suffix='_V_d05.nc',
                               start_year=init_year, no_years=num_years)

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'POST/', fieldname1[:-3], '_', 'ANN_',
                 str(init_year), '-', str(init_year + num_years - 1)]), field)

# --------------------------------------------------------------------------- #
