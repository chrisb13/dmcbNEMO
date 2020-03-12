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
fieldname = 'siconc'

# Specify the years to do the averaging over and the season to average.
month_to_ave = range(12)
# month_to_ave = [ 0, 7 ]
init_year = 2003
num_years = 5

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Loop over the number of years to make the required monthly average.
for k in month_to_ave:
    field = nemo.calc_monthly_ave(fieldname, tmask, 'I', homedir, nemodir,
                                  'd01/I/',
                                  prefix='ORCH0083-LIM3_', suffix='_I_d01.nc',
                                  start_year=init_year, no_years=num_years,
                                  month=k+1)

    # Save the averaged field to a numpy matrix file for later plotting.
    np.save(''.join([homedir, nemodir, 'POST/', fieldname, '_', str(k+1).rjust(2, '0'), '_',
                     str(init_year), '-', str(init_year + num_years - 1)]),
            field)

# --------------------------------------------------------------------------- #
