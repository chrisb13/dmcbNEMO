#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import glosst

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the dense colormap.
plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
glosstdir = 'MO-GLO-SST/'

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(glosst.load_field('lat', homedir, glosstdir, 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20141231.nc'))
lon = np.squeeze(glosst.load_field('lon', homedir, glosstdir, 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20141231.nc'))
lon, lat = np.meshgrid(lon[:, None], lat[:, None])
lon = lon.T
lat = lat.T

# Load the temperature fields and remask them.
sst2 = np.squeeze(np.load(''.join([homedir, glosstdir, 'POST/',
                          'analysed_sst_sq_ANN_2014-2016.npy'])))
sst = np.squeeze(np.load(''.join([homedir, glosstdir, 'POST/',
                         'analysed_sst_ANN_2014-2016.npy'])))
mask = np.zeros(sst.shape)
mask[:] = sst
mask[mask > 0] = 1
sst2 = glosst.mask_field(sst2, mask)
sst = glosst.mask_field(sst, mask)

# --------------------------------------------------------------------------- #

# Calculate the temperature variance.
tvar = sst2 - sst*sst
sst = sst - 273.15

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
cs_sst = map.contourf(lon.T, lat.T, np.log10(np.squeeze(tvar)).T,
                      np.arange(-1.5, 1.51, 0.025),
                      latlon='true', cmap='thermal', vmin=-1.5, vmax=1.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-1.5, 1.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('glosst_ss_tvar.png', bbox_inches='tight', dpi=600)

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
cs_sst = map.contourf(lon.T, lat.T, np.squeeze(sst).T,
                      np.arange(-2., 30.1, 1.),
                      latlon='true', cmap='thermal', vmin=-2., vmax=30.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(-2., 30.1, 4.), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('glosst_ss_sst.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
