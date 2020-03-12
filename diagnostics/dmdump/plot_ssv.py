#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the balance colormap.
plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
filename = 'TIDY/ARCHIVE/1948/V/ORCH0083-LIM3_19481227_V_d05.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamv = np.squeeze(nemo.load_field('glamv', homedir, nemodir, gridfile, 'V'))
gphiv = np.squeeze(nemo.load_field('gphiv', homedir, nemodir, gridfile, 'V'))

# Load the relevant mask.
vmask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Load the field to plot and remask it.
ssv = nemo.load_field('sossvssv', homedir, nemodir, filename, 'V')
ssv = nemo.mask_field(ssv, vmask)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure()

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-10,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)

# Draw the mask to fill in the continents.
cs_vmask = map.pcolormesh(glamv.T, gphiv.T, 1.0-np.squeeze(vmask[:, :, 0].T),
                          latlon='true', cmap='Greys')

# Draw the contour plot.
cs_ssu = map.contourf(glamv.T, gphiv.T, np.squeeze(ssv).T,  np.arange(-2., 2.0625, 0.0625),
                      latlon='true', cmap='balance')

# Add a colour bar.
cbar = map.colorbar(cs_ssu, 'right', ticks=np.arange(-2., 2.5, 0.5), pad=0.5)

# Show the plot.
plt.show()

# --------------------------------------------------------------------------- #
