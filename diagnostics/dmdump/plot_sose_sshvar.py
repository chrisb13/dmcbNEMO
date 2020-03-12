#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import sose

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the speed colormap.
plt.register_cmap(name='dense', cmap=cmocean.cm.dense)

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
ssh2 = np.load(''.join([postdir, 'ETAN_sq_ANN_074-438.npy']))[:, :, None]
ssh2 = np.ma.masked_array(ssh2, 1.0-tmask)

ssh = np.load(''.join([postdir, 'ETAN_ANN_074-438.npy']))[:, :, None]
ssh = np.ma.masked_array(ssh, 1.0-tmask)

# --------------------------------------------------------------------------- #

# Calculate the SSH variance.
sshvar = ssh2 - ssh*ssh

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
                      np.log10(1.E4*np.squeeze(sshvar)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('sose_ss_sshvar.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
