#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the thermal colormap.
plt.register_cmap(name='ice', cmap=cmocean.cm.ice)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
filename = 'TIDY/ARCHIVE/1948/I/ORCH0083-LIM3_19481227_I_d05.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the coordinates.
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Load the field to plot and remask it.
ice = nemo.load_field('siconc', homedir, nemodir, filename, 'T')
ice = nemo.mask_field(ice, tmask)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-10,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)

# Draw the mask to fill in the continents.
cs_tmask = map.pcolormesh(glamt.T, gphit.T, 1.0-np.squeeze(tmask[:, :, 0].T),
                          latlon='true', cmap='Greys')

# Draw the contour plot.
cs_mxl = map.contourf(glamt.T, gphit.T, np.squeeze(ice).T,
                      np.arange(0.1, 1.025, 0.025),
                      latlon='true', cmap='ice')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0.1, 1.225, 0.225), pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('ice1948.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
