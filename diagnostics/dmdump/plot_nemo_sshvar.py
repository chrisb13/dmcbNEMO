#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the dense colormap.
plt.register_cmap(name='dense', cmap=cmocean.cm.dense)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP02/'

# Specify the names of the different files that I want to load from.
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Load the velocity fields and remask them.
ssh2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosssh2_ANN_1952-1954.npy'])))[:, :, None]
ssh2 = nemo.mask_field(ssh2, tmask)
ssh = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sossheig_ANN_1952-1954.npy'])))[:, :, None]
ssh = nemo.mask_field(ssh, tmask)

# --------------------------------------------------------------------------- #

# Calculate the SSH variance.
sshvar = ssh2 - ssh*ssh

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=-150, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25, color='w')
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25,
                  color='w')
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_ssz = map.contourf(glamt.T, gphit.T, np.log10(1.E4*np.squeeze(sshvar)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('../figs/EXP02/nemo_ss_sshvar.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
