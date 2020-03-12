#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import mxl

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the dense colormap.
plt.register_cmap(name='deep', cmap=cmocean.cm.deep)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
mxldir = 'DEBOYER-MXL/'

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(mxl.load_axis('lat', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
lon = np.squeeze(mxl.load_axis('lon', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
lon = np.concatenate((lon, lon[-1:]+2), axis=0)
lon, lat = np.meshgrid(lon[:, None], lat[:, None])
lon = lon.T
lat = lat.T

# Load the mixed layer depth fields and remask them.
mask = np.squeeze(mxl.load_field('mask', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
mld = np.squeeze(mxl.load_field('mld', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))

# Expand the mask out to have monthly values.
mask = np.repeat(mask[:, :, None], 12, axis=2)

# Make negative values masked.
mask[mld < 0] = 0

# Mask the mixed layer depth.
mld = mxl.mask_field(mld, mask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 100.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(lon.T, lat.T, np.squeeze(mld[:,:,1]).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='deep', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('deboyer_mxl_01.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 500.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(lon.T, lat.T, np.squeeze(mld[:,:,8]).T,
                      np.arange(0., max_mxl+1, max_mxl/20),
                      latlon='true', cmap='deep', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('deboyer_mxl_01.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
