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
fieldname = 'e3u'

# Specify the years to do the averaging over and the season to average.
init_year = 2003
num_years = 5

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to make the require seasonal
# average.
field = nemo.calc_annual_ave(fieldname, umask, 'U', homedir, nemodir, 'm01/U/',
                             prefix='ORCH0083-LIM3_', suffix='_U_m01.nc',
                             start_year=init_year, no_years=num_years)

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to extract the required files.
# field = nemo.extract_field(fieldname, tmask, 'T', homedir, nemodir, 'T/',
#                           prefix='ORCH0083-LIM3_', suffix='_T_d05.nc',
#                           start_year=init_year, no_years=num_years)

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'POST/', fieldname, '_', 'ANN_',
                 str(init_year), '-', str(init_year + num_years - 1)]), field)

# --------------------------------------------------------------------------- #
