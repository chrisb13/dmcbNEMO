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
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55FIL-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/200[34567]/m01/T/*'
udir = 'TIDY/ARCHIVE/200[34567]/m01/U/*'
vdir = 'TIDY/ARCHIVE/200[34567]/m01/V/*'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the years to do the averaging over.
init_year = 2003 # Sets the output file names.
num_years = 5

# Set physical constants.
g = 9.80665 # gravitational acceleration.

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant masks.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))[:, :, 0:1]
fmask = np.squeeze(nemo.load_field('fmask', homedir, nemodir, gridfile, 'Z'))[:, :, 0:1]

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
utaux = 0.0
vtauy = 0.0

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(tfiles), 1):
    print(tfiles[k])

    # Load the SSH field and remask it.
    ssh = nemo.load_field('sossheig', '', '', tfiles[k], 'T')
    ssh = nemo.mask_field(ssh, tmask)

    # Take the X & Y gradients of the SSH field.
    dhdx = umask * nemo.grad_dtdx(ssh.data, e1u, umask, nx, ny, 1)
    dhdy = vmask * nemo.grad_dtdy(ssh.data, e2v, vmask, nx, ny, 1)

    # Average the gradients on to the other velocity point to allow for
    # geostrophic velocity calculation.
    dhdx = vmask * nemo.ave_u_onto_v(dhdx, nx, ny, 1, vmask)
    dhdy = umask * nemo.ave_v_onto_u(dhdy, nx, ny, 1, umask)

    # Convert the gradients into geostrophic velocity estimates.
    dhdx = 2.0 * g * dhdx * vmask / (f[0:nx, :, :] + f[1:nx+1, :, :])
    dhdy = -2.0 * g * dhdy * umask / (f[:, 0:ny, :] + f[:, 1:ny+1, :])

    # Load the tau fields and remask them.
    taux = nemo.load_field('sozotaux', '', '', ufiles[k], 'U')
    taux = nemo.mask_field(taux, umask)
    tauy = nemo.load_field('sometauy', '', '', vfiles[k], 'V')
    tauy = nemo.mask_field(tauy, vmask)

    # Multiply the tau? by dhd? & accumulate them onto their averages.
    utaux = utaux + taux * dhdy
    vtauy = vtauy + tauy * dhdx

    # Increment the counter by 1.
    counter = counter + 1

# --------------------------------------------------------------------------- #

# Now the loop over all the files is complete, divide the accumulations by counter.
utaux = utaux / counter
vtauy = vtauy / counter

# Save the averaged field to a numpy matrix file for later plotting.
np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/tauxdhdy_m01_ANN_',
                 str(init_year), '-', str(init_year + num_years - 1)]), utaux.data)
np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/tauydhdx_m01_ANN_',
                 str(init_year), '-', str(init_year + num_years - 1)]), vtauy.data)


# --------------------------------------------------------------------------- #
