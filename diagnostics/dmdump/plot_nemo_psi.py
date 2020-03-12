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
plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

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

# Load the grid spacings.
e2u = nemo.load_field('e2u', homedir, nemodir, gridfile, 'U')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))

# Load the velocity fields and remask them.
uvele3u = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                                      'vozoce3u_ANN_1952-1954.npy'])))
uvele3u = nemo.mask_field(uvele3u, umask)

# --------------------------------------------------------------------------- #

# Calculate the barotropic streamfunction.
psi = uvele3u.sum(axis=2)
psi = np.cumsum(psi.data[:, :, None]*umask[:, :, 0:1]*e2u[:, :, 0:1], axis=1)
psi = nemo.mask_field(psi, umask[:, :, 0:1])

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
cs_ssz = map.contourf(glamt.T, gphit.T, np.squeeze(psi[0:-1, :, :]).T/1.E6,
                      np.arange(-50., 261., 5.),
                      latlon='true', cmap='balance', vmin=-50., vmax=260.,
                      extend='both')
foo = map.contour(glamt.T, gphit.T, np.squeeze(psi[0:-1, :, :]).T/1.E6,
                  np.arange(-260., 261., 25.), colors='k',
                  latlon='true')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-50., 261., 50.), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('../figs/EXP02/nemo_ss_psi.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
