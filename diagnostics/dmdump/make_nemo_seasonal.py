#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/TIDY/ARCHIVE/'

# Specify the field to average.
fieldname = 'siconc'

# Specify the years to do the averaging over and the season to average.
season_to_ave = ('JFM', 'AMJ', 'JAS', 'OND')
init_year = 1960
num_years = 1

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to make the require seasonal
# average.
for k in season_to_ave:
    field = nemo.calc_seasonal_ave(fieldname, tmask, 'T', homedir, nemodir,
                                   'd01/I/',
                                   prefix='ORCH0083-LIM3_', suffix='_I_d01.nc',
                                   start_year=init_year, no_years=num_years,
                                   season=k)

    # Save the averaged field to a numpy matrix file for later plotting.
    if k is not 'DJF':
        np.save(''.join([homedir, nemodir, 'POST/', fieldname, '_', k, '_', str(init_year),
                         '-', str(init_year + num_years - 1)]), field)
    else:
        np.save(''.join([homedir, nemodir, 'POST/', fieldname, '_', k, '_', str(init_year-1),
                         '-', str(init_year + num_years - 1)]), field)

# --------------------------------------------------------------------------- #
