#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import sose

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/SOSE/'
sosedir = ''
postdir = '/Users/munday/Documents/Data/SOSE/POST/'

month_to_plot = 1

# --------------------------------------------------------------------------- #

# Load the grid.
sosegrid = sose.load_grid(homedir)

# Make the relevant masks.
tmask = sosegrid['hFacC'][:, :, 0:1]
tmask[tmask > 0] = 1
tmask = np.squeeze(tmask)

# Load the field to plot.
mxl = sose.load_field('MLD_mnthlyBar', homedir, sosedir,
                      iteration_number=100,
                      records=[month_to_plot+11, month_to_plot+23,
                               month_to_plot+35, month_to_plot+47,
                               month_to_plot+59],
                      levels=(),
                      fieldtype='T').mean(axis=2)
mxl = sose.mask_field(-1.0*mxl, tmask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 100.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spstere', boundinglat=-35, lon_0=-150, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.5, color='w')
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True],
                  linewidth=0.5, color='w')
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T, np.squeeze(mxl).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='seismic', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_mxl_01_2006-2010.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
