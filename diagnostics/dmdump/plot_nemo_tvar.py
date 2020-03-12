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
plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tfilename = 'TIDY/ARCHIVE/1952/ORCH0083-LIM3_1952_T_d01.nc'
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

# Load the temperature fields and remask them.
sst2 = nemo.load_field('sosstsqu', homedir, nemodir, tfilename, 'T')
sst2 = nemo.mask_field(sst2, tmask)
sst = nemo.load_field('sosstsst', homedir, nemodir, tfilename, 'T')
sst = nemo.mask_field(sst, tmask)

# --------------------------------------------------------------------------- #

# Calculate the temperature variance.
tvar = sst2 - sst*sst

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
cs_sst = map.contourf(glamt.T, gphit.T, np.log10(np.squeeze(tvar)).T,
                      np.arange(-1.5, 1.51, 0.025),
                      latlon='true', cmap='thermal', vmin=-1.5, vmax=1.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-1.5, 1.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_sss_tvar.png']),
            bbox_inches='tight', dpi=600)

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
cs_sst = map.contourf(glamt.T, gphit.T, np.squeeze(sst).T,
                      np.arange(-2., 30.1, 1.),
                      latlon='true', cmap='thermal', vmin=-2., vmax=30.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-2., 30.1, 4.), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_sss_sst.png']),
            bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
