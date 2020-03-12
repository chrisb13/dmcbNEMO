#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import glosst
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
sst2 = spio.loadmat(''.join([postdir,
                             'CT_TEOS10_sq_ANN_074-438.mat']))['average'][:, :, 0:1]
sst2 = np.ma.masked_array(sst2, 1.0-tmask)

sst = spio.loadmat(''.join([postdir,
                            'CT_TEOS10_ANN_074-438.mat']))['average'][:, :, 0:1]
sst = np.ma.masked_array(sst, 1.0-tmask)

# Calculate the temperature variance.
tvar = sst2 - sst*sst

# Integrate the temperature variance.
sose_tvar = (tvar * sosegrid['DXC'][:, :, None] * tmask).sum(axis=0) / \
            (sosegrid['DXC'][:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del tmask, sst2, sst, tvar

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

# Load the zonal grid spacing and latitudes
e1t = np.squeeze(nemo.load_field('e1t', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sst2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosstsqu_ANN_1952-1954.npy'])))
sst2 = nemo.mask_field(sst2[:, :, None], tmask)
sst = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosstsst_ANN_1952-1954.npy'])))
sst = nemo.mask_field(sst[:, :, None], tmask)

# Calculate the temperature variance.
tvar = sst2 - sst*sst

# Integrate the temperature variance.
nemo_tvar_exp00 = (tvar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP01/'

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sst2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosstsqu_ANN_1952-1954.npy'])))
sst2 = nemo.mask_field(sst2[:, :, None], tmask)
sst = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosstsst_ANN_1952-1954.npy'])))
sst = nemo.mask_field(sst[:, :, None], tmask)

# Calculate the temperature variance.
tvar = sst2 - sst*sst

# Integrate the temperature variance.
nemo_tvar_exp01 = (tvar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP02/'

# --------------------------------------------------------------------------- #

# Load the temperature fields and remask them.
sst2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'sosstsqu_ANN_1952-1954.npy'])))
sst2 = nemo.mask_field(sst2[:, :, None], tmask)
sst = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'sosstsst_ANN_1952-1954.npy'])))
sst = nemo.mask_field(sst[:, :, None], tmask)

# Calculate the temperature variance.
tvar = sst2 - sst*sst

# Integrate the temperature variance.
nemo_tvar_exp02 = (tvar * e1t[:, :, None] * tmask).sum(axis=0) / \
                  (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del e1t, tmask, sst2, sst, nx, ny, nz, tvar

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
glosstdir = 'MO-GLO-SST/'

# Specify the names of the different files that I want to load from.
gridfile = 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20141231.nc'

# Specify the number of grid boxes.
nx = 7200
ny = 1101

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(glosst.load_field('lat', homedir, glosstdir, gridfile))

# Load the temperature fields and remask them.
sst2 = np.squeeze(np.load(''.join([homedir, glosstdir, 'POST/',
                          'analysed_sst_sq_ANN_2014-2016.npy'])))
sst = np.squeeze(np.load(''.join([homedir, glosstdir, 'POST/',
                         'analysed_sst_ANN_2014-2016.npy'])))
mask = np.zeros(sst.shape)
mask[:] = sst
mask[mask > 0] = 1
sst2 = glosst.mask_field(sst2, mask)
sst = glosst.mask_field(sst, mask)

tvar = sst2 - sst*sst
tvar.mask[tvar < 0] = 1

# Calculate the temperature variance.
glosst_tvar = tvar.mean(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del sst2, sst, nx, ny, tvar

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(sosegrid['YC'][0, :], sose_tvar, '-b',
             lat, glosst_tvar, '-k',
             gphit.mean(axis=0), nemo_tvar_exp00, '-g',
             gphit.mean(axis=0), nemo_tvar_exp01, '-r',
             gphit.mean(axis=0), nemo_tvar_exp02, '-m',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/6.0)


plt.axis([-70.0, -30.0, 0, 6.0])
plt.xticks(np.linspace(-70.0, -30.0, 5))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean Temperature Variance ($\mathrm{K}^2$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2], p[3], p[4]],
           labels=['SOSE', 'MO-GLO-SST', 'NEMO - EXP00', 'NEMO - EXP01', 'NEMO - EXP02'])

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

plt.savefig('tvar_comparison.pdf', bbox_inches='tight')

# --------------------------------------------------------------------------- #
