#!/opt/local/bin/python
#
# --------------------------------------------------------------------------- #

# My handrolled modules.
import sys
sys.path.append('/nerc/n01/n01/munday/python/nemo')
import nemo

# Import required modules.
import numpy as np

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/VOLCHECK/EXP05/'

# Specify the names of the different files that I want to load from.
tfilename = 'OUTPUTS/ORCH0083-LIM3_1h_19480101_19480101_grid_T_1948010101-1948010101.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')
e3t_0 = np.squeeze(nemo.load_field('e3t_0', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))

# Mask the e3's.
e3t_0 = nemo.mask_field(e3t_0, tmask)

# --------------------------------------------------------------------------- #

e3t = np.squeeze(nemo.load_field('e3t', homedir, nemodir, tfilename, 'T'))
e3t = nemo.mask_field(e3t, tmask)

# --------------------------------------------------------------------------- #

vt_0 = e1t * e2t * e3t_0
vt_0 = vt_0.sum()
print 'vt_0 =', vt_0

vt = e1t * e2t * e3t
vt = vt.sum()
print 'vt =', vt

# --------------------------------------------------------------------------- #

e3t = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'e3t_ANN_1952-1952.npy'])))
e3t = nemo.mask_field(e3t, tmask)

sivolu = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'sivolu_ANN_1952-1952.npy'])))
sivolu = nemo.mask_field(sivolu, np.squeeze(tmask[:, :, 0:1]))

vt1952 = e1t * e2t * e3t
vt1952 = vt1952.sum()
sivolu1952 = np.squeeze(e1t * e2t) * sivolu
sivolu1952 = sivolu1952.sum()
print 'vt1952 =', vt1952
print 'sivolu1952 =', sivolu1952
print 'vt1952+sivolu1952 =', vt1952+sivolu1952

# --------------------------------------------------------------------------- #

e3t = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'e3t_ANN_1953-1953.npy'])))
e3t = nemo.mask_field(e3t, tmask)

sivolu = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'sivolu_ANN_1953-1953.npy'])))
sivolu = nemo.mask_field(sivolu, np.squeeze(tmask[:, :, 0:1]))

vt1953 = e1t * e2t * e3t
vt1953 = vt1953.sum()
sivolu1953 = np.squeeze(e1t * e2t) * sivolu
sivolu1953 = sivolu1953.sum()
print 'vt1953 =', vt1953
print 'sivolu1953 =', sivolu1953
print 'vt1953+sivolu1953 =', vt1953+sivolu1953

# --------------------------------------------------------------------------- #

e3t = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'e3t_ANN_1954-1954.npy'])))
e3t = nemo.mask_field(e3t, tmask)

sivolu = np.squeeze(np.load(''.join([homedir, 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/', 'TIDY/ARCHIVE/POST/',
                          'sivolu_ANN_1954-1954.npy'])))
sivolu = nemo.mask_field(sivolu, np.squeeze(tmask[:, :, 0:1]))

vt1954 = e1t * e2t * e3t
vt1954 = vt1954.sum()
sivolu1954 = np.squeeze(e1t * e2t) * sivolu
sivolu1954 = sivolu1954.sum()
print 'vt1954 =', vt1954
print 'sivolu1954 =', sivolu1954
print 'vt1954+sivolu1954 =', vt1954+sivolu1954

# --------------------------------------------------------------------------- #
