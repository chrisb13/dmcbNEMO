#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the deep colormap.
plt.register_cmap(name='deep', cmap=cmocean.cm.deep)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55ABS-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
filename = 'TIDY/ARCHIVE/1952/d01/T/ORCH0083-LIM3_19521231_T_d01.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the coordinates.
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Load the field to plot and remask it.
mxl = nemo.load_field('somxlr103', homedir, nemodir, filename, 'T')
mxl = nemo.mask_field(mxl, tmask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 100.

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
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(glamt.T, gphit.T, np.squeeze(mxl).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='deep', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('mxl1954.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
