#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import aviso

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the dense colormap.
plt.register_cmap(name='dense', cmap=cmocean.cm.dense)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
avisodir = 'AVISO/'

# Specify the names of the different files that I want to load from.
gridfile = ('dataset-duacs-rep-global-merged-allsat-phy-l4-v3_' +
            '20140101-20141231_adt.nc')

# Specify the number of grid boxes.
nx = 1440
ny = 221

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(aviso.load_field('longitude', homedir, avisodir, gridfile))
gphit = np.squeeze(aviso.load_field('latitude', homedir, avisodir, gridfile))
glamt, gphit = np.meshgrid(glamt[:, None], gphit[:, None])
glamt = glamt.T
gphit = gphit.T

# Load the velocity fields and remask them.
ssh2 = np.load(''.join([homedir, avisodir, 'POST/',
                        'adt_sq_ANN_2014-2016.npy']))
ssh = np.load(''.join([homedir, avisodir, 'POST/', 'adt_ANN_2014-2016.npy']))

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
cs_ssz = map.contourf(glamt.T, gphit.T, np.log10(1.E4*np.squeeze(sshvar)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='dense', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.75)

# Show the plot.
# plt.show()

# Save a hires version of the figure
# plt.savefig('aviso_ss_sshvar.png', bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #
