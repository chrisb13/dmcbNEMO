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
plt.register_cmap(name='haline', cmap=cmocean.cm.haline)

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
sss2 = spio.loadmat(''.join([postdir,
                             'SA_TEOS10_sq_ANN_074-438.mat']))['average'][:, :, 0:1]
sss2 = np.ma.masked_array(sss2, 1.0-tmask)

sss = spio.loadmat(''.join([postdir,
                            'SA_TEOS10_ANN_074-438.mat']))['average'][:, :, 0:1]
sss = np.ma.masked_array(sss, 1.0-tmask)

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
cs_sst = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T,
                      np.log10(np.squeeze(svar)).T,
                      np.arange(-3.5, 1.01, 0.025),
                      latlon='true', cmap='haline', vmin=-3.5, vmax=1.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-3.5, 1.01, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_ss_svar.png', bbox_inches='tight', dpi=600)

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
cs_sst = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T, np.squeeze(sss).T,
                      np.arange(34., 36.01, 0.0125),
                      latlon='true', cmap='haline', vmin=34., vmax=36.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(34., 36.01, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_ss_sss.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
