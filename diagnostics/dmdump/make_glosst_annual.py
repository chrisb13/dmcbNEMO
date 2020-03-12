#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import glosst

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
glosstdir = 'MO-GLO-SST/'

# Specify the field to average
fieldname = 'analysed_sst'

# Specify the years to do the averaging over and the season to average.
init_year = 2014
num_years = 3

# Specify what power to raise the field to.
power = 1

# --------------------------------------------------------------------------- #

# Loop over the number of years and months to make the require seasonal
# average.
field = glosst.calc_annual_ave(fieldname, homedir, glosstdir,
                               'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_', '.nc',
                               start_year=init_year, no_years=num_years,
                               power_to_raise=power)

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
if power == 1:
    np.save(''.join([fieldname, '_', 'ANN_', str(init_year), '-',
                     str(init_year + num_years - 1)]), field.data)
elif power == 2:
    np.save(''.join([fieldname, '_sq_', 'ANN_', str(init_year), '-',
                     str(init_year + num_years - 1)]), field.data)

# --------------------------------------------------------------------------- #
