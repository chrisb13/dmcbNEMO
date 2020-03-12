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
plt.register_cmap(name='haline', cmap=cmocean.cm.haline)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
sfilename = 'TIDY/ARCHIVE/1952/ORCH0083-LIM3_1952_T_d01.nc'
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

# Load the salinity fields and remask them.
sss2 = nemo.load_field('sosalsqu', homedir, nemodir, sfilename, 'T')
sss2 = nemo.mask_field(sss2, tmask)
sss = nemo.load_field('sosaline', homedir, nemodir, sfilename, 'T')
sss = nemo.mask_field(sss, tmask)

# --------------------------------------------------------------------------- #

# Calculate the temperature variance.
svar = sss2 - sss*sss

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
cs_sst = map.contourf(glamt.T, gphit.T, np.log10(np.squeeze(svar)).T,
                      np.arange(-3.5, 1.01, 0.025),
                      latlon='true', cmap='haline', vmin=-3.5, vmax=1.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-3.5, 1.01, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_ss_svar.png']),
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
cs_sst = map.contourf(glamt.T, gphit.T, np.squeeze(sss).T,
                      np.arange(34., 36.01, 0.0125),
                      latlon='true', cmap='haline', vmin=34., vmax=36.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(34., 36.01, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/nemo_ss_sss.png']),
            bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
