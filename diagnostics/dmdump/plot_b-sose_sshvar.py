#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import sose

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

# Register the speed colormap.
plt.register_cmap(name='dense', cmap=cmocean.cm.dense)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/B-SOSE/'
sosedir = 'POST/'

# --------------------------------------------------------------------------- #

# Open the grid file for read-only access.
ncfile = Dataset(''.join([homedir, 'grid.nc']), mode='r')

# Load the required fields from the specified file.
XC = ncfile.variables['XC'][:]
YC = ncfile.variables['YC'][:]
hFacC = ncfile.variables['hFacC'][0:1, :, :]

# Close the file, for the sake of form.
ncfile.close()

# Make a temperature mask.
tmask = np.zeros(hFacC.shape)
tmask[:] = hFacC
tmask[tmask > 0] = 1

# --------------------------------------------------------------------------- #

# Load the SSH variance and remask it.

# Open the file for read-only access.
ncfile = Dataset(''.join([homedir, sosedir,
                          'bsose_i105_2009to2012_4year_SSH.nc']), mode='r')

# Load the SSH variance field.
sshvar = ncfile.variables['SSHvar'][:]

# Close the file, for the sake of form.
ncfile.close()

sshvar = np.ma.masked_array(sshvar, 1.0-tmask)

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
cs_ssz = map.contourf(XC.T, YC.T,
                      np.log10(1.E4*np.squeeze(sshvar)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('b-sose_ss_sshvar.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
