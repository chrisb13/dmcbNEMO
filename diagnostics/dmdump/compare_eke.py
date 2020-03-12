#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import aviso
import nemo
import sose

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt

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

umask = sosegrid['hFacW'][:, :, 0:1]
umask[umask > 0] = 1
umask = sose.pad_ocean_field(umask, 'U')

vmask = sosegrid['hFacS'][:, :, 0:1]
vmask[vmask > 0] = 1
vmask = sose.pad_ocean_field(vmask, 'V')

# Load the velocity fields and remask them.
ssu2 = np.load(''.join([postdir, 'UVEL_sq_ANN_074-438.npy']))[:, :, 0:1]
ssu2 = np.ma.masked_array(ssu2, umask)

ssu = np.load(''.join([postdir, 'UVEL_ANN_074-438.npy']))[:, :, 0:1]
ssu = np.ma.masked_array(ssu, umask)

ssv2 = np.load(''.join([postdir, 'VVEL_sq_ANN_074-438.npy']))[:, :, 0:1]
ssv2 = np.ma.masked_array(ssv2, vmask)

ssv = np.load(''.join([postdir, 'VVEL_ANN_074-438.npy']))[:, :, 0:1]
ssv = np.ma.masked_array(ssv, vmask)

# --------------------------------------------------------------------------- #

# Calculate the Eddy Kinetic Energy.
eke = sose.calc_eke(ssu.data, ssu2.data, ssv.data, ssv2.data,
                    2160, 320, 1, tmask)

# Integrate the Eddy Kinetic Energy.
sose_eke = (sosegrid['DXC'][:, :, None] * tmask * eke).sum(axis=0) / \
           (sosegrid['DXC'][:, :, None] * tmask).sum(axis=0)

# Clean up.
del tmask, umask, vmask, ssu2, ssu, ssv2, ssv, eke

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
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Load the velocity fields and remask them.
ssu2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'uocesq_e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu2 = nemo.mask_field(ssu2, umask)
ssu = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vozocrtx_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu = nemo.mask_field(ssu, umask)

ssv2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'vocesq_e3v_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv2 = nemo.mask_field(ssv2, vmask)
ssv = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vomecrty_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv = nemo.mask_field(ssv, vmask)

e3t = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3t_ANN_1952-1954.npy'])))[:, :, 0:1]
e3u = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
e3v = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3v_ANN_1952-1954.npy'])))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Calculate the Eddy Kinetic Energy.
eke = nemo.calc_eke(ssu.data, ssu2.data, ssv.data, ssv2.data,
                    e3t, e3u, e3v, nx, ny, 1, tmask)

# Integrate the Eddy Kinetic Energy.
nemo_eke_exp00 = (eke * e1t[:, :, None] * e3t * tmask).sum(axis=0) / \
                 (e1t[:, :, None] * e3t * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP01/'

# --------------------------------------------------------------------------- #

# Load the velocity fields and remask them.
ssu2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'uocesq_e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu2 = nemo.mask_field(ssu2, umask)
ssu = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vozocrtx_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu = nemo.mask_field(ssu, umask)

ssv2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'vocesq_e3v_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv2 = nemo.mask_field(ssv2, vmask)
ssv = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vomecrty_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv = nemo.mask_field(ssv, vmask)

e3t = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3t_ANN_1952-1954.npy'])))[:, :, 0:1]
e3u = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
e3v = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3v_ANN_1952-1954.npy'])))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Calculate the Eddy Kinetic Energy.
eke = nemo.calc_eke(ssu.data, ssu2.data, ssv.data, ssv2.data,
                    e3t, e3u, e3v, nx, ny, 1, tmask)

# Integrate the Eddy Kinetic Energy.
nemo_eke_exp01 = (eke * e1t[:, :, None] * e3t * tmask).sum(axis=0) / \
                 (e1t[:, :, None] * e3t * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
nemodir = 'trunk/NEMOGCM/CONFIG/DRM-ORCH0083-LIM3/EXP02/'

# --------------------------------------------------------------------------- #

# Load the velocity fields and remask them.
ssu2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'uocesq_e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu2 = nemo.mask_field(ssu2, umask)
ssu = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vozocrtx_ANN_1952-1954.npy'])))[:, :, 0:1]
ssu = nemo.mask_field(ssu, umask)

ssv2 = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                          'vocesq_e3v_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv2 = nemo.mask_field(ssv2, vmask)
ssv = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'vomecrty_ANN_1952-1954.npy'])))[:, :, 0:1]
ssv = nemo.mask_field(ssv, vmask)

e3t = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3t_ANN_1952-1954.npy'])))[:, :, 0:1]
e3u = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3u_ANN_1952-1954.npy'])))[:, :, 0:1]
e3v = np.squeeze(np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/',
                         'e3v_ANN_1952-1954.npy'])))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Calculate the Eddy Kinetic Energy.
eke = nemo.calc_eke(ssu.data, ssu2.data, ssv.data, ssv2.data,
                    e3t, e3u, e3v, nx, ny, 1, tmask)

# Integrate the Eddy Kinetic Energy.
nemo_eke_exp02 = (eke * e1t[:, :, None] * e3t * tmask).sum(axis=0) / \
                 (e1t[:, :, None] * e3t * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

# Clean up.
del e1t, tmask, umask, vmask, ssu2, ssu, ssv2, ssv, e3t, e3u, e3v, nx, ny, nz
del eke

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
avisodir = 'AVISO/'

# Specify the names of the different files that I want to load from.
gridfile = ('dataset-duacs-rep-global-merged-allsat-phy-l4-v3_' +
            '20140101-20141231_uvgos.nc')

# Specify the number of grid boxes.
nx = 1440
ny = 221

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(aviso.load_field('latitude', homedir, avisodir, gridfile))

# Load the velocity fields and remask them.
ssu2 = np.load(''.join([homedir, avisodir, 'POST/',
                        'ugos_sq_ANN_2014-2016.npy']))
ssu = np.load(''.join([homedir, avisodir, 'POST/', 'ugos_ANN_2014-2016.npy']))

ssv2 = np.load(''.join([homedir, avisodir, 'POST/',
                        'vgos_sq_ANN_2014-2016.npy']))
ssv = np.load(''.join([homedir, avisodir, 'POST/', 'vgos_ANN_2014-2016.npy']))

# --------------------------------------------------------------------------- #

# Calculate the Eddy Kinetic Energy.
eke = aviso.calc_eke(ssu, ssu2, ssv, ssv2)

mask = np.zeros(eke.shape)
mask[:] = eke
mask[mask > 0] = 1

# Integrate the Eddy Kinetic Energy.
aviso_eke = (np.pi * 0.25 * 6371.0E3 * np.cos(np.pi * lat[None, :] / 180.0) *
             eke / 180.0).sum(axis=0) / \
            (np.pi * 0.25 * 6371.0E3 * np.cos(np.pi * lat[None, :] / 180.0) *
             mask / 180.0).sum(axis=0)

# Clean up.
del ssu2, ssu, ssv2, ssv, nx, ny, mask, eke

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(sosegrid['YC'][0, :], 1.E4*sose_eke, '-b',
             lat, 1.E4*aviso_eke, '-k',
             gphit.mean(axis=0), 1.E4*nemo_eke_exp00, '-g',
             gphit.mean(axis=0), 1.E4*nemo_eke_exp01, '-r',
             gphit.mean(axis=0), 1.E4*nemo_eke_exp02, '-m',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/300.)


plt.axis([-70.0, -30.0, 0, 300.])
plt.xticks(np.linspace(-70.0, -30.0, 5))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean EKE ($\mathrm{cm}^2/s^{-1}$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2], p[3], p[4]],
           labels=['SOSE', 'AVISO', 'NEMO - EXP00', 'NEMO - EXP01', 'NEMO - EXP02'])

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

plt.savefig('eke_comparison.pdf', bbox_inches='tight')

# --------------------------------------------------------------------------- #
