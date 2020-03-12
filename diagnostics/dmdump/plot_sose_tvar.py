#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import sose

# Import required modules.
import cmocean
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the speed colormap.
plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/SOSE/'
sosedir = 'RAWDATA/'
postdir = '/Users/munday/Documents/Data/SOSE/POST/'

# --------------------------------------------------------------------------- #

# Load the grid.
sosegrid = sose.load_grid(homedir)

# Make the relevant masks.
tmask = sosegrid['hFacC'][:, :, 0:1]
tmask[tmask > 0] = 1

# Load the velocity fields and remask them.
sst2 = spio.loadmat(''.join([postdir,
                             'CT_TEOS10_sq_ANN_074-438.mat']))['average'][:, :, 0:1]
sst2 = np.ma.masked_array(sst2, 1.0-tmask)

sst = spio.loadmat(''.join([postdir,
                            'CT_TEOS10_ANN_074-438.mat']))['average'][:, :, 0:1]
sst = np.ma.masked_array(sst, 1.0-tmask)

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
cs_sst = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T,
                      np.log10(np.squeeze(tvar)).T,
                      np.arange(-1.5, 1.51, 0.025),
                      latlon='true', cmap='thermal', vmin=-1.5, vmax=1.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-1.5, 1.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_ss_tvar.png', bbox_inches='tight', dpi=600)

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
cs_sst = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T, np.squeeze(sst).T,
                      np.arange(-2., 30.1, 1.),
                      latlon='true', cmap='thermal', vmin=-2., vmax=30.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-2., 30.1, 4.), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_ss_sst.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
