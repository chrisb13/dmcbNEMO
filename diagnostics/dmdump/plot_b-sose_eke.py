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

# Load the velocity fields and remask them.
ssu2 = np.load(''.join([postdir, 'UVEL_sq_ANN_074-438.npy']))[:, :, 0:1]
ssu2 = np.ma.masked_array(ssu2, umask)

ssu = np.load(''.join([postdir, 'UVEL_ANN_074-438.npy']))[:, :, 0:1]
ssu = np.ma.masked_array(ssu, umask)

ssv2 = np.load(''.join([postdir, 'VVEL_sq_ANN_074-438.npy']))[:, :, 0:1]
ssv2 = np.ma.masked_array(ssv2, vmask)

ssv = np.load(''.join([postdir, 'VVEL_ANN_074-438.npy']))[:, :, 0:1]
ssv = np.ma.masked_array(ssv, vmask)

# --------------------------------------------------------------------------- #

# Calculate the Total Kinetic Energy.
tke = sose.calc_eke(np.zeros(ssu.data.shape), ssu2.data,
                    np.zeros(ssv.data.shape), ssv2.data, 2160, 320, 1, tmask)

# Calculate the Eddy Kinetic Energy.
eke = sose.calc_eke(ssu.data, ssu2.data, ssv.data, ssv2.data,
                    2160, 320, 1, tmask)

# Calculate the Mean Kinetic Energy.
mke = tke - eke

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
cs_ssz = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T,
                      np.log10(1.E4*np.squeeze(tke)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
# plt.savefig('sose_ss_tke.png', bbox_inches='tight', dpi=600)

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
cs_ssz = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T,
                      np.log10(1.E4*np.squeeze(eke)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
# plt.savefig('sose_ss_eke.png', bbox_inches='tight', dpi=600)

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
cs_ssz = map.contourf(sosegrid['XC'].T, sosegrid['YC'].T,
                      np.log10(1.E4*np.squeeze(mke)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
# plt.savefig('sose_ss_mke.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
