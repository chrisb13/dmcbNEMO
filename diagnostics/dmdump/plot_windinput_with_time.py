#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
import matplotlib
matplotlib.rcParams['backend'] = 'Agg'
import matplotlib.pyplot as plt

plt.ioff()

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/COAREJRA-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/m01/U/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load volumes.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamu = np.squeeze(nemo.load_field('glamu', homedir, nemodir, gridfile, 'U'))
gphiu = np.squeeze(nemo.load_field('gphiu', homedir, nemodir, gridfile, 'U'))

# Load the grid spacings.
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2u = nemo.load_field('e2u', homedir, nemodir, gridfile, 'U')

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculate the ice volume.

if save_output:
    # Preallocate the output variable.
    wind = np.ndarray(shape=[len(ufiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        print(ufiles[k])
        # Load & mask the current ice fraction field.
        zonalwind = nemo.load_field('sozotaux', '', '', ufiles[k], 'U')
        zonalwind = nemo.mask_field(zonalwind, umask)

        # Calculate the area of the grid box covered in ice.
        total = (e1u[0:nx-1,:,:] * e2u[0:nx-1,:,:] * zonalwind[0:nx-1,:,:]).sum()

        # Store the total area for output.
        wind[k, 0] = ufiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_U', 1)[0]
        wind[k, 1] = total

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/zonalwindinput_m01']), wind)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    wind = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/zonalwindinput.npy']))

# --------------------------------------------------------------------------- #
