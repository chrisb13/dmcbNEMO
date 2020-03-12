#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55FIL-ORCH0083-LIM3/EXP00/TIDY/ARCHIVE/'

# Specify the field to average
fieldname1 = 'sossussu'
fieldname2 = 'sozotaux'

# Specify the years to do the averaging over and the season to average.
init_year = 2003
num_years = 5

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to make the require seasonal
# average.
field = nemo.calc_annual_product(fieldname1, fieldname2, 
                               umask, 'U', homedir, nemodir, 'd01/U/',
                               prefix='ORCH0083-LIM3_', suffix='_U_d01.nc',
                               start_year=init_year, no_years=num_years)

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'POST/', fieldname1, '_', fieldname2, '_', 'ANN_',
                 str(init_year), '-', str(init_year + num_years - 1)]), field)

# --------------------------------------------------------------------------- #
