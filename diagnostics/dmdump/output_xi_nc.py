#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
import netCDF4 as nc4
import time

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/1985/d05/T/'
udir = 'TIDY/ARCHIVE/1985/d05/U/'
vdir = 'TIDY/ARCHIVE/1985/d05/V/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load KE values.
save_output = 0

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile))
depth = np.squeeze(nemo.load_field('nav_lev', homedir, nemodir, gridfile))

# Load the grid spacings.
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2v = nemo.load_field('e2u', homedir, nemodir, gridfile, 'V')
e1f = nemo.load_field('e1f', homedir, nemodir, gridfile, 'Z')
e2f = nemo.load_field('e2f', homedir, nemodir, gridfile, 'Z')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))
fmask = np.squeeze(nemo.load_field('fmask', homedir, nemodir, gridfile, 'Z'))

# --------------------------------------------------------------------------- #
# Find the U and V files.

tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))
ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(len(ufiles)):
    print ufiles[k]
    # Load & mask the current U velocity field.
    U = nemo.load_field('vozoce3u', '', '', ufiles[k], 'U')
    U = nemo.mask_field(np.squeeze(U), umask)
    e3u = nemo.load_field('e3u', '', '', ufiles[k], 'U')
    e3u = nemo.mask_field(np.squeeze(e3u), umask)
    U = U / e3u
    del e3u

    # Load and mask the current V velocity field.
    V = nemo.load_field('vomece3v', '', '', vfiles[k], 'V')
    V = nemo.mask_field(np.squeeze(V), vmask)
    e3v = nemo.load_field('e3v', '', '', vfiles[k], 'V')
    e3v = nemo.mask_field(np.squeeze(e3v), vmask)
    V = V / e3v
    del e3v

    # Calculate the KE for the current velocity fields
    XIe3t = nemo.calc_xi(U.data, V.data, nx, ny, nz, e1u, e2v, e1f, e2f, fmask)
    del U, V

    e3t = nemo.load_field('e3t', '', '', tfiles[k], 'T')
    e3t = nemo.mask_field(np.squeeze(e3t), tmask)


    # Average the XI field on to T points.
    XIe3t = 0.25*(XIe3t[0:4320,0:2000,:] + XIe3t[1:4321,0:2000,:] + XIe3t[0:4320,1:2001,:] + XIe3t[1:4321,1:2001,:]) * e3t

    # Add extra columns to match NEMO arrays.
    XIe3t = np.concatenate((XIe3t[-1:,:,:], XIe3t[:,:,:], XIe3t[0:1,:,:]))

    # Output the XI field to a netcdf file.

    # Create the dataset.
    filename=''.join([homedir, nemodir, udir[0:-2], 'Z/ORCH0083-LIM3_', ufiles[k][-17:-9], '_Z_d05.nc'])

    dataset = nc4.Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    t = dataset.createDimension('time_counter', None)
    z = dataset.createDimension('deptht', 75)
    y = dataset.createDimension('y', 2000)
    x = dataset.createDimension('x', 4320+2)

    # Create the coordinate variables.
    nav_lat = dataset.createVariable('nav_lat', np.float32, ('y', 'x'))
    nav_lon = dataset.createVariable('nav_lon', np.float32, ('y', 'x'))
    nav_lev = dataset.createVariable('nav_lev', np.float32, 'deptht')

    # Create the actual variable
    relativevorticity = dataset.createVariable('XIe3t', np.float32, ('time_counter', 'deptht', 'y', 'x'), zlib=True, complevel=9)

    # Variable the dimensions and the variables.
    relativevorticity[:] = XIe3t[:,:,:,None].T
    nav_lat[:] = gphit.T
    nav_lon[:] = glamt.T

    # Global Attributes to record a little metadata.
    dataset.description = 'Netcdf output of plot_ssk.py python script'
    dataset.history = 'Created ' + time.ctime(time.time())
    dataset.source = ufiles[k] + ' ' + vfiles[k]

    # Variable Attributes
    relativevorticity.long_name = 'Vertical component of relative vorticity * e3t'
    relativevorticity.units = '/s'
    nav_lat.standard_name = 'latitude'
    nav_lat.long_name = 'Latitude'
    nav_lat.units = 'degrees_north'
    nav_lon.standard_name = 'longitude'
    nav_lon.long_name = 'Longitude'
    nav_lon.units = 'degrees_east'
    nav_lon.nav_model = 'grid_U'

    # Write and close the dataset.
    dataset.close()

# --------------------------------------------------------------------------- #
