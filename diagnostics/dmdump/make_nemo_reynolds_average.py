#!/opt/local/bin/python

#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np
from six.moves import cPickle as pickle

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Set physical constants.
cp = 3991.86795711963 # Specific heat capacity (J/kg).
rho0 = 1026.000000000 # Boussinesq reference density (kg/m^3)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55FIL-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
years = '2003-2007'

# Setup the filenames from where data will be loaded.
e3tfilename = ''.join(['e3t_ANN_', years, '.npy'])
e3vfilename = ''.join(['e3v_ANN_', years, '.npy'])

tfilename = ''.join(['voteme3t_ANN_', years, '.npy'])

vfilename = ''.join(['vomece3v_ANN_', years, '.npy'])
vtfilename = ''.join(['vomtem3d_ANN_', years, '.npy'])

gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the relevant masks.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))

# Load the coordinates.
gphiv = np.squeeze(nemo.load_field('gphiv', homedir, nemodir, gridfile, 'V'))
gphiv = nemo.mask_field(gphiv, vmask[:,:,0])

# Load the zonal grid spacing.
e1v = np.squeeze(nemo.load_field('e1v', homedir, nemodir, gridfile, 'V'))
e1v = nemo.mask_field(e1v, vmask[:,:,0])

# --------------------------------------------------------------------------- #

# Load the e3t field.
e3t = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/', e3tfilename])))
e3t = nemo.mask_field(e3t, tmask)

# Load the e3v field.
e3v = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/', e3vfilename])))
e3v = nemo.mask_field(e3v, vmask)

# --------------------------------------------------------------------------- #

# Load the VT field & factor out the vertical grid spacing.
vt = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/', vtfilename])))
vt = nemo.mask_field(vt, vmask) / e3v / rho0 / cp / e1v[:, :, None] # Diagnostic has factor of rho0*cp*e1v already.

# --------------------------------------------------------------------------- #

# Load the theta field & factor out the vertical grid spacing.
t = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/', tfilename])))
t = nemo.mask_field(t, tmask) / e3t

# Average theta onto the V points.
t = nemo.ave_t_onto_v(t, nx, ny, nz, vmask)

# --------------------------------------------------------------------------- #

# Load the V field & factor out the vertical grid spacing.
v = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/', vfilename])))
v = nemo.mask_field(v, vmask) / e3v

# --------------------------------------------------------------------------- #

# Calculate the total heat transport and Reynolds average it into different components.
heat_transport = nemo.reynolds_average_vt(vt, v, t, e1v, e3t, e3v)

# --------------------------------------------------------------------------- #

# Save the heat_transport dictionary to a pickle file.
with open(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/heat_transport_', years, '.pkl']), 'wb') as f:
    pickle.dump(heat_transport, f)

# --------------------------------------------------------------------------- #
