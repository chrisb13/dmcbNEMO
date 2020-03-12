#!/opt/local/bin/python
#
# --------------------------------------------------------------------------- #
# This python script will step through the interpolation of David Hutchinson's
# CM2.1 configuration onto an MITgcm grid.
# --------------------------------------------------------------------------- #

# Import handrolled ECCO-handling and interpolation modules.
import ecco4
import mitgcm
import nemo

# Import required modules.
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# Specify the location of the ECCO v4 data set.

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
eccodir = 'ECCO4/interp_climatology/'

# Specify the file and variable name to interpolate.
filename = 'SALT.0001.nc'
varname = 'SALT'

# Specifu the number of grid boxes in the donor mesh to make life easier.
ecco = {'nx': 720, 'ny': 360, 'nz': 50, 'nt': 12}

# --------------------------------------------------------------------------- #
# Load the field on the donor mesh and the grid box coordinates.

lat = mitgcm.pad_field(ecco4.load_field('lat', homedir, eccodir, filename),
                       0.5)
lon = mitgcm.pad_field(ecco4.load_field('lon', homedir, eccodir, filename),
                       0.0)
dep = ecco4.load_field('dep', homedir, eccodir, filename)
field = mitgcm.pad_ocean_field(np.squeeze(
                                          ecco4.load_field(varname, homedir,
                                                           eccodir, filename)))

# --------------------------------------------------------------------------- #
# Load the longitudes and latitudes for the target mesh (NEMO).

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
gridfile = 'INPUTS/EN4.1.1_75L_monthly_19952014_reg1d_C_0_TEOS-10.nc'

# Specify the number of grid boxes.
nemogrid = {'nx': 360, 'ny': 180, 'nz': 75}

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('nav_lon', homedir, nemodir, gridfile))-180.
gphit = np.squeeze(nemo.load_field('nav_lat', homedir, nemodir, gridfile))
glamt, gphit = np.meshgrid(glamt[:, None], gphit[:, None])
glamt = glamt.T
gphit = gphit.T
gdept_0 = np.squeeze(nemo.load_field('deptht', homedir, nemodir, gridfile))

# --------------------------------------------------------------------------- #
# ECCO is masked with NaNs, so use a nearest neighbour interpolation to flood
# fill them.

# Construct a mask that is TRUE in LAND.
mask = np.squeeze(field[:, :, :, 0:1].copy())
mask[~np.isnan(mask)] = 0.0  # Ocean is 0.0.
mask[np.isnan(mask)] = 1.0   # Land is 1.0.
mask = mask.astype(bool)

print 'At floodfill of field'

# Floodfill the 3D ocean field as a series of 2D operations.
for l in range(ecco['nt']):
    print 'l =', l
    for k in range(ecco['nz']):
        field[mask[:, :, k], k, l] = mitgcm.floodfill_field(lon[~mask[:, :, k]],
                                                            lat[~mask[:, :, k]],
                                                            field[~mask[:, :, k], k, l],
                                                            lon[mask[:, :, k]],
                                                            lat[mask[:, :, k]])

# --------------------------------------------------------------------------- #
# Do a series of 2D interpolations from the ECCO v4 1/2o grid to the NEMO
# 1/12o grid.

# Predefine the output array.
nemofield = np.zeros((nemogrid['nx'], nemogrid['ny'], ecco['nz'], ecco['nt']))

print 'At interpolation of nemofield'

# Interpolate the 3D ocean field as a series of 2D operations.
for l in range(ecco['nt']):
    print 'l =', l
    for k in range(ecco['nz']):
        nemofield[:, :, k, l] = mitgcm.interpolate_field(lon, lat, field[:, :, k, l],
                                                         ecco['nx']+2, ecco['ny'],
                                                         glamt, gphit,
                                                         nemogrid['nx'],
                                                         nemogrid['ny'])

# --------------------------------------------------------------------------- #
# Load the depths of the grid boxes to allow for interpolation in the vertical.

finalfield = np.zeros((nemogrid['nx'], nemogrid['ny'], nemogrid['nz'], ecco['nt']))

print 'At interpolation of finalfield'

for l in range(ecco['nt']):
    print 'l =', l
    for j in range(nemogrid['ny']):
        for i in range(nemogrid['nx']):
            f = sp.interpolate.interp1d(np.concatenate(([0], dep), axis=0),
                                        np.concatenate((nemofield[i, j, 0:1, l],
                                                        nemofield[i, j, :, l]),
                                                       axis=0))
            finalfield[i, j, :, l] = f(gdept_0)

# --------------------------------------------------------------------------- #

np.save(''.join(['ECCOv4_', varname, '_EN4.1']), finalfield)
sp.io.savemat(''.join(['ECCOv4_', varname, '_EN4.1']), {varname: finalfield})

# --------------------------------------------------------------------------- #
