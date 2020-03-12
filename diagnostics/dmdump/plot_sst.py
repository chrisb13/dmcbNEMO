#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
# import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)

# Register the thermal colormap.
# plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
filename = 'TIDY/ARCHIVE/1951/d01/T/ORCH0083-LIM3_19511231_T_d01.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the coordinates.
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Load the field to plot and remask it.
sst = nemo.load_field('sosstsst', homedir, nemodir, filename, 'T')
sst = nemo.mask_field(sst, tmask)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the mask to fill in the continents.
# cs_tmask = map.pcolormesh(glamt.T, gphit.T, 1.0-np.squeeze(tmask[:, :, 0].T),
#                           latlon='true', cmap='Greys')

# Draw the contour plot.
cs_sst = map.contourf(glamt.T, gphit.T, np.squeeze(sst).T,
                      np.arange(-2., 22.1, 0.1),
                      latlon='true', cmap='viridis', vmin=-2., vmax=22.,
                      extend='both')

# cs_sst.cmap.set_under([0.015556013331540799, 0.13824424546464084,
#                       0.20181088645583051, 1.0])
# cs_sst.cmap.set_over([0.90904184166740365, 0.98215740632167059,
#                      0.35550780642995311, 1.0])

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-2., 22., 4.), pad=0.5)

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/sst.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
