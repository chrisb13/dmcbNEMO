#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import mxl
import nemo

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
mxldir = 'DEBOYER-MXL/'

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(mxl.load_axis('lat', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
lon = np.squeeze(mxl.load_axis('lon', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
lon = np.concatenate((lon, lon[-1:]+2), axis=0)
lon, lat = np.meshgrid(lon[:, None], lat[:, None])
lon = lon.T
lat = lat.T

# Load the mixed layer depth fields and remask them.
mask = np.squeeze(mxl.load_field('mask', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))
mld = np.squeeze(mxl.load_field('mld', homedir, mxldir, 'mld_DR003_c1m_reg2.0.nc'))

# Expand the mask out to have monthly values.
mask = np.repeat(mask[:, :, None], 12, axis=2)

# Make negative values masked.
mask[mld < 0] = 0

# Mask the mixed layer depth.
mld = mxl.mask_field(mld, mask)

# --------------------------------------------------------------------------- #

# Zonally average the mixed layer depth.
mld = mld.mean(axis=0)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

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

somxlr103_01 = np.squeeze(np.load(''.join([homedir, nemodir,
                                           'TIDY/ARCHIVE/POST/',
                                           'somxlr103_01_1952-1952.npy'])))
somxlr103_01 = nemo.mask_field(somxlr103_01[:, :, None], tmask)
somxlr103_08 = np.squeeze(np.load(''.join([homedir, nemodir,
                                           'TIDY/ARCHIVE/POST/',
                                           'somxlr103_08_1952-1952.npy'])))
somxlr103_08 = nemo.mask_field(somxlr103_08[:, :, None], tmask)

# --------------------------------------------------------------------------- #

# Zonallly average the Eddy Kinetic Energy.
somxlr103_01 = (somxlr103_01 * e1t[:, :, None] * tmask).sum(axis=0) / \
                 (e1t[:, :, None] * tmask).sum(axis=0)
somxlr103_08 = (somxlr103_08 * e1t[:, :, None] * tmask).sum(axis=0) / \
                 (e1t[:, :, None] * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(lat[0, :].T, mld[:, 0], '-b',
             lat[0, :].T, mld[:, 7], '-g',
             gphit.mean(axis=0), somxlr103_01, '-r',
             gphit.mean(axis=0), somxlr103_08, '-m',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/250.)


plt.axis([-70.0, -30.0, 0, 250.])
plt.xticks(np.linspace(-70.0, -30.0, 5))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean MLD ($\mathrm{m}$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2], p[3]],
           labels=['de Boyer - Jan', 'de Boyer - Aug',
           'NEMO - Jan', 'NEMO - Aug'])

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
