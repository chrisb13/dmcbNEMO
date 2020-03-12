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
sss = spio.loadmat(''.join([postdir,
                            'SA_TEOS10_ANN_074-438.mat']))['average'][:, :, 0:1]
sss = np.ma.masked_array(sss, 1.0-tmask)

# Integrate the temperature variance.
sose_sss = (sss * sosegrid['DXC'][:, :, None] * tmask).sum(axis=0) / \
            (sosegrid['DXC'][:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del tmask, sss

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP02/'

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
# Load the field to plot and remask it.
sss = nemo.load_field('vosaline', homedir, nemodir, 'INIT/vosaline.1949.nc', 'T')[:, :, 0:1, 0]
sss = nemo.mask_field(sss, tmask)

# Integrate the temperature variance.
nemo_sss_init = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1948-1948.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1948 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1949-1949.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1949 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1950-1950.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1950 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1951-1951.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1951 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1952-1952.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1952 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1953-1953.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1953 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sss = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosaline_08_1954-1954.npy'])))
sss = nemo.mask_field(sss[:, :, None], tmask)

# Integrate the temperature variance.
nemo_sss_1954 = (sss * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del e1t, tmask, sss, nx, ny, nz

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(sosegrid['YC'][0, :], sose_sss, '--k',
             gphit.mean(axis=0), nemo_sss_init, '-k',
             gphit.mean(axis=0), nemo_sss_1948, '-b',
             gphit.mean(axis=0), nemo_sss_1949, '-g',
             gphit.mean(axis=0), nemo_sss_1950, '-r',
             gphit.mean(axis=0), nemo_sss_1951, '-m',
             gphit.mean(axis=0), nemo_sss_1952, '-y',
             gphit.mean(axis=0), nemo_sss_1953, '--b',
             gphit.mean(axis=0), nemo_sss_1954, '--g',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/2.0)


plt.axis([-70.0, -30.0, 34.0, 36.0])
plt.xticks(np.linspace(-70.0, -30.0, 5))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean Salinity ($\mathrm{g}\mathrm{/kg}$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2], p[3],
                    p[4], p[5], p[6], p[7], p[8]], labels=['SOSE',
                                                           'NEMO - INIT',
                                                           'NEMO - 1948',
                                                           'NEMO - 1949',
                                                           'NEMO - 1950',
                                                           'NEMO - 1951',
                                                           'NEMO - 1952',
                                                           'NEMO - 1953',
                                                           'NEMO - 1954'])

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

plt.savefig('check_sss_exp02_comparison.pdf', bbox_inches='tight')

# --------------------------------------------------------------------------- #
