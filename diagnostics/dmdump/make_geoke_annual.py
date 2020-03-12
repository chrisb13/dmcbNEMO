#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
#homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2TAU-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/1972/d01/T/'
udir = 'TIDY/ARCHIVE/19972/d01/U/'
vdir = 'TIDY/ARCHIVE/19972/d01/V/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the years to do the averaging over.
init_year = 1972 # Sets the output file names.
num_years = 1

# Set physical constants.
g = 9.80665 # gravitational acceleration.

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the relevant masks.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Load the grid spacings.
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2v = nemo.load_field('e2v', homedir, nemodir, gridfile, 'V')

# Load the Coriolis parameter.
f = nemo.load_field('ff_f', homedir, nemodir, gridfile, 'Z')

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to the geostrophic wind work for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))
ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #

# Initialise a counter to zero and arrays to hold utaux & vtauy.
counter = 0.0
ssu = 0.0
ssv = 0.0
ssk = 0.0

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(tfiles), 1):
    print(tfiles[k])

    # Load the SSH field and remask it.
    ssh = nemo.load_field('sossheig', '', '', tfiles[k], 'T')
    ssh = nemo.mask_field(ssh, tmask)

    # Calculate the geostrophic velocities.
    ssu = -2.0 * g * nemo.grad_dtdy(ssh.data, e2v, vmask, nx, ny, 1) / (f[0:nx, :, :] + f[1:nx+1, :, :])
    ssu = ssu * vmask
    ssv = 2.0 * g * nemo.grad_dtdx(ssh.data, e1u, umask, nx, ny, 1)/ (f[:, 0:ny, :] + f[:, 1:ny+1, :])
    ssv = ssv * umask

    # Calculate the geostrophic KE for the current geostrophic velocity fields
    ssk = ssk + nemo.calc_ke(ssv, ssu, nx, ny, 1, tmask)

    # Increment the counter by 1.
    counter = counter + 1

# --------------------------------------------------------------------------- #

# Now the loop over all the files is complete, divide the accumulations by counter.
ssk = ssk / counter

# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/geoke_ANN_', str(init_year), '-', str(init_year + num_years - 1)]), ssk.data)

# --------------------------------------------------------------------------- #
