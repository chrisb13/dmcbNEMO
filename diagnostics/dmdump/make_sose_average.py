#!/opt/local/bin/python
#
# --------------------------------------------------------------------------- #

# My handrolled modules.
import sose

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/SOSE/'
sosedir = 'RAWDATA/'

# Specify the names of the different files that I want to load from.
tfilename = 'TIDY/ARCHIVE/1949/T/ORCH0083-LIM3_19491227_T_d05.nc'

# Specify the field to average
fieldname = 'ETAN'

# Specify what power to raise the field to.
power = 2

# Specify the indices to begin and end on.
start = 73     # First index at end of year 1/start of year 2 is 73 (0-based).
final = 437    # Final index should be ~< 437, as there are 438 5 day means.

# --------------------------------------------------------------------------- #

field = sose.calc_annual_average(fieldname, homedir, sosedir,
                                 iteration_number=100,
                                 start_index=start, final_index=final,
                                 power_to_raise=power, fieldtype='T')

# --------------------------------------------------------------------------- #
# Save the averaged field to a numpy matrix file for later plotting.
if power == 1:
    np.save(''.join([fieldname, '_', 'ANN_', str(start+1).rjust(3, '0'), '-',
                     str(final+1).rjust(3, '0')]), field)
elif power == 2:
    np.save(''.join([fieldname, '_sq_', 'ANN_', str(start+1).rjust(3, '0'),
                     '-', str(final+1).rjust(3, '0')]), field)

# --------------------------------------------------------------------------- #
