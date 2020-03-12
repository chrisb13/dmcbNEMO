#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the balance colormap.
plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

# --------------------------------------------------------------------------- #

# Turn interactive plotting off
plt.ioff()

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
ufilename = 'TIDY/ARCHIVE/1952/d01/U/ORCH0083-LIM3_19521231_U_d01.nc'
vfilename = 'TIDY/ARCHIVE/1952/d01/V/ORCH0083-LIM3_19521231_V_d01.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamf = np.squeeze(nemo.load_field('glamf', homedir, nemodir, gridfile, 'Z'))
gphif = np.squeeze(nemo.load_field('gphif', homedir, nemodir, gridfile, 'Z'))

# Load the grid spacings.
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2v = nemo.load_field('e2u', homedir, nemodir, gridfile, 'V')
e1f = nemo.load_field('e1f', homedir, nemodir, gridfile, 'Z')
e2f = nemo.load_field('e2f', homedir, nemodir, gridfile, 'Z')

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]
fmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'Z'))[:, :, 0:1]

# Load the velocity fields and remask them.
ssu = nemo.load_field('sossussu', homedir, nemodir, ufilename, 'U')
ssu = nemo.mask_field(ssu, umask)
ssv = nemo.load_field('sossvssv', homedir, nemodir, vfilename, 'V')
ssv = nemo.mask_field(ssv, vmask)

# --------------------------------------------------------------------------- #

ssz = nemo.calc_xi(ssu.data, ssv.data, nx, ny, 1, e1u, e2v, e1f, e2f, fmask)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('ssz.png', bbox_inches='tight', dpi=1200)

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-10,
              lon_0=180, round='true')
# map.shadedrelief()
# map.fillcontinents(color='black', lake_color='white')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_ssz = map.contourf(glamf.T, gphif.T, np.squeeze(ssz).T,
                      np.arange(-3.05E-5, 3.15E-5, 1.E-6),
                      latlon='true', cmap='balance', vmin=-3.E-5, vmax=3.E-5,
                      extend='both')

cs_ssz.cmap.set_under([0.093176301801157851, 0.11117332947760272,
                       0.26151238855305475, 1.0])
cs_ssz.cmap.set_over([0.23605636466461405, 0.035297479946040287,
                      0.069437442394125581, 1.0])

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-3.E-5, 4.5E-5, 1.5E-5),
                    pad=0.5, format='%.0e')

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('ssz.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
