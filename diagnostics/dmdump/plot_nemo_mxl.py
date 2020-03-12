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
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

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

# Load the mixed layer fields remask them.
somxlr103_01 = np.squeeze(np.load(''.join([homedir, nemodir,
                                           'TIDY/ARCHIVE/POST/',
                                           'somxlr103_01_1952-1952.npy'])))
somxlr103_01 = nemo.mask_field(somxlr103_01[:, :, None], tmask)
somxlr103_08 = np.squeeze(np.load(''.join([homedir, nemodir,
                                           'TIDY/ARCHIVE/POST/',
                                           'somxlr103_08_1952-1952.npy'])))
somxlr103_08 = nemo.mask_field(somxlr103_08[:, :, None], tmask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 100.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(glamt.T, gphit.T, np.squeeze(somxlr103_01).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='deep', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_ss_mxl_01.png']), 
            bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 500.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(glamt.T, gphit.T, np.squeeze(somxlr103_08).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='deep', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_ss_mxl_08.png']), 
            bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #