#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cmocean

# Register the ice colormap.
plt.register_cmap(name='ice', cmap=cmocean.cm.ice)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/TIDY/ARCHIVE/'

# Specify the names of the different files that I want to load from.
gridfile = 'MESH/mesh_mask.nc'

# Specify the filename.
filename = 'siconc_DJF_1951-1954.npy'

# --------------------------------------------------------------------------- #

# Load the coordinates.
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

ice = np.load(''.join([homedir, nemodir, 'POST/', filename]))
ice = nemo.mask_field(ice, tmask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_ice = 1.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spstere', boundinglat=-35, lon_0=-150, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.5, color='w')
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True],
                  linewidth=0.5, color='w')
map.drawmapboundary(fill_color='black')
map.bluemarble()

# Draw the mask to fill in the continents - replaced by using drawmapboundary.
# cs_tmask = map.contourf(nemo.pad_field(glamt, 0.0083).T,
#                        nemo.pad_field(gphit, 0.).T,
#                        1.0-np.squeeze(nemo.pad_ocean_field(tmask)).T,
#                        latlon='true', cmap='Greys')

# Draw the contour plot.
cs_ice = map.contourf(nemo.pad_field(glamt, 0.0083).T,
                      nemo.pad_field(gphit, 0.).T,
                      np.squeeze(nemo.pad_ocean_field(ice)).T,
                      np.arange(0., max_ice+0.01, max_ice/51),
                      latlon='true', cmap='ice', vmin=0., vmax=max_ice,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ice, 'right',
                    ticks=np.arange(0., max_ice+0.01, max_ice/4), pad=0.75)

# Save a hires version of the figure
plt.savefig(''.join([filename[:-3], 'png']), bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
