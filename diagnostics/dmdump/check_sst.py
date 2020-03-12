#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/POST/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask)

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, 'sosstsst*'])))

# --------------------------------------------------------------------------- #

# Preallocate the output variable.
meansst = np.ndarray(shape=[len(tfiles), 1])

for k in range(len(tfiles)):
        print tfiles[k]
        # Load & mask the current SST field.
        sst = np.load(tfiles[k])
        sst = nemo.mask_field(sst, tmask)

        # Calculate the average value and save in the meansst array.

        # Calculate the area average surface temperature.
        meansst[k, 0] = (sst*area).sum() / area.sum()

# --------------------------------------------------------------------------- #
