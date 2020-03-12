#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
from mpl_toolkits.basemap import Basemap

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
datadir = '/Users/munday/Documents/Data/BYU-IB/'

# Load the iceberg lats and lons.
lat = spio.loadmat(''.join([datadir, 'lat.mat']))
lon = spio.loadmat(''.join([datadir, 'lon.mat']))

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(12.0, 12.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spstere', boundinglat=-35, lon_0=-150, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.5, color='w')
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True],
                  linewidth=0.5, color='w')
map.drawmapboundary(fill_color='black')
map.bluemarble()

map.plot(lon['lon'], lat['lat'], 'r.', latlon='true', markersize=1)

# Save a hires version of the figure
# plt.savefig(''.join([filename[:-3], 'png']), bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
