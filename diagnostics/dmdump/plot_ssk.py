#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# from netCDF4 import Dataset
# import time

# Register the speed colormap.
plt.register_cmap(name='speed', cmap=cmocean.cm.speed)

# --------------------------------------------------------------------------- #

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
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Load the velocity fields and remask them.
ssu = nemo.load_field('sossussu', homedir, nemodir, ufilename, 'U')
ssu = nemo.mask_field(ssu, umask)
ssv = nemo.load_field('sossvssv', homedir, nemodir, vfilename, 'V')
ssv = nemo.mask_field(ssv, vmask)

# --------------------------------------------------------------------------- #

ssk = nemo.calc_ke(ssu.data, ssv.data, nx, ny, 1, tmask)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

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
cs_ssz = map.contourf(glamt.T, gphit.T, np.log10(1.E4*np.squeeze(ssk)).T,
                      np.arange(-0.5, 3.51, 0.025),
                      latlon='true', cmap='speed', vmin=-0.5, vmax=3.5,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5), pad=0.5)

# Show the plot.
# plt.show()

# Save a hires version of the figure
plt.savefig('ssk.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
# Output the ssk field to a netcdf file.

# Create the dataset.
# dataset = Dataset('surfaceke.nc', 'w', format='NETCDF4_CLASSIC')
# ny = dataset.createDimension('ny', 2000)
# nx = dataset.createDimension('nx', 4320)

# Create the coordinate variables.
# nav_lat = dataset.createVariable('nav_lat', np.float32, ('ny', 'nx'))
# nav_lon = dataset.createVariable('nav_lon', np.float32, ('ny', 'nx'))

# Create the actual variable
# surfacekineticenergy = dataset.createVariable('ssk', np.float32, ('ny', 'nx'))

# Variable the dimensions and the variables.
# surfacekineticenergy[:] = np.squeeze(ssk).T
# nav_lat[:] = gphit.T
# nav_lon[:] = glamt.T

# Global Attributes to record a little metadata.
# dataset.description = 'Netcdf output of plot_ssk.py python script'
# dataset.history = 'Created ' + time.ctime(time.time())
# dataset.source = ufilename + ' ' + vfilename

# Variable Attributes
# surfacekineticenergy.long_name = 'Surface kinetic energy'
# surfacekineticenergy.units = 'm2/s2'
# nav_lat.standard_name = 'latitude'
# nav_lat.long_name = 'Latitude'
# nav_lat.units = 'degrees_north'
# nav_lon.standard_name = 'longitude'
# nav_lon.long_name = 'Longitude'
# nav_lon.units = 'degrees_east'
# nav_lon.nav_model = 'grid_U'

# Write and close the dataset.
# dataset.close()

# --------------------------------------------------------------------------- #
