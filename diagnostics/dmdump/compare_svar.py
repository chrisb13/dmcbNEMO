#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo
import sose

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio

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
sss2 = spio.loadmat(''.join([postdir,
                             'SA_TEOS10_sq_ANN_074-438.mat']))['average'][:, :, 0:1]
sss2 = np.ma.masked_array(sss2, 1.0-tmask)

sss = spio.loadmat(''.join([postdir,
                            'SA_TEOS10_ANN_074-438.mat']))['average'][:, :, 0:1]
sss = np.ma.masked_array(sss, 1.0-tmask)

# Calculate the temperature variance.
svar = sss2 - sss*sss

# Integrate the temperature variance.
sose_svar = (svar * sosegrid['DXC'][:, :, None] * tmask).sum(axis=0) / \
            (sosegrid['DXC'][:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del tmask, sss2, sss, svar

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the zonal grid spaing and latitudes
e1t = np.squeeze(nemo.load_field('e1t', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosalsqu_ANN_1952-1954.npy'])))
sss2 = nemo.mask_field(sss2[:, :, None], tmask)
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_ANN_1952-1954.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Calculate the temperature variance.
svar = sss2 - sss*sss

# Integrate the temperature variance.
nemo_svar_exp00 = (svar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP01/'

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosalsqu_ANN_1952-1954.npy'])))
sss2 = nemo.mask_field(sss2[:, :, None], tmask)
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_ANN_1952-1954.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Calculate the temperature variance.
svar = sss2 - sss*sss

# Integrate the temperature variance.
nemo_svar_exp01 = (svar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP02/'

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosalsqu_ANN_1952-1954.npy'])))
sss2 = nemo.mask_field(sss2[:, :, None], tmask)
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_ANN_1952-1954.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Calculate the temperature variance.
svar = sss2 - sss*sss

# Integrate the temperature variance.
nemo_svar_exp02 = (svar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del e1t, tmask, sss2, sss, nx, ny, nz, svar

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(sosegrid['YC'][0, :], sose_svar, '-b',
             gphit.mean(axis=0), nemo_svar_exp00, '-g',
             gphit.mean(axis=0), nemo_svar_exp01, '-r',
             gphit.mean(axis=0), nemo_svar_exp02, '-m',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/0.5)


plt.axis([-70.0, -30.0, 0, 0.5])
plt.xticks(np.linspace(-70.0, -30.0, 5))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean Salinity Variance ($\mathrm{g}^2\mathrm{/kg}^2$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2], p[3]],
           labels=['SOSE', 'NEMO - EXP00', 'NEMO - EXP01', 'NEMO - EXP02'])

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')

# --------------------------------------------------------------------------- #
# Save the figure to a pdf.

plt.savefig('svar_comparison.pdf', bbox_inches='tight')

# --------------------------------------------------------------------------- #
