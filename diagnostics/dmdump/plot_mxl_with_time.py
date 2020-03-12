#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
import matplotlib
matplotlib.rcParams['backend'] = 'Agg'
import matplotlib.pyplot as plt

plt.ioff()

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/COAREJRA-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/????/d01/T/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load MXL values.
save_output = 1

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

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask)

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average MXL depth.

if save_output:
    # Preallocate the output variable.
    mxl = np.ndarray(shape=[len(tfiles), 4])

    # Loop over the T and load the mixed layer depths.
    for k in range(len(tfiles)):
        print(tfiles[k])
        # Load & mask the current MXL field.
        somxl010 = nemo.load_field('somxlr101', '', '', tfiles[k], 'T')
        somxl010 = nemo.mask_field(somxl010, tmask)
        somxl030 = nemo.load_field('somxlr103', '', '', tfiles[k], 'T')
        somxl030 = nemo.mask_field(somxl030, tmask)
        somxlkara = nemo.load_field('somxlkara', '', '', tfiles[k], 'T')
        somxlkara = nemo.mask_field(somxlkara, tmask)

        # Calculate the area average surface KE.
        mxl[k, 0] = tfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_T', 1)[0]
        mxl[k, 1] = (somxl010*area).sum() / area.sum()
        mxl[k, 2] = (somxl030*area).sum() / area.sum()
        mxl[k, 3] = (somxlkara*area).sum() / area.sum()

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/mxl']), mxl)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    mxl = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/mxl.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(1, len(mxl), len(mxl))/365.,
             mxl[:, 1],
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/200.0)

plt.axis([0, 40.0, 0, 200.0])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.arange(0., 201.0, 50.))
ax.set_yticks(np.arange(0, 201.0, 10.), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. MXL depth ($\mathrm{m}$)',
           fontname='arial', fontsize=30, fontweight='bold')

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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/mxl010_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(1, len(mxl), len(mxl))/365.,
             mxl[:, 2],
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/200.0)

plt.axis([0, 40.0, 0, 200.0])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.arange(0., 201.0, 50.))
ax.set_yticks(np.arange(0, 201.0, 10.), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. MXL depth ($\mathrm{m}$)',
           fontname='arial', fontsize=30, fontweight='bold')

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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/mxl030_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(1, len(mxl), len(mxl))/365.,
             mxl[:, 3],
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/200.0)

plt.axis([0, 40.0, 0, 200.0])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.arange(0., 201.0, 50.))
ax.set_yticks(np.arange(0, 201.0, 10.), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. MXL depth ($\mathrm{m}$)',
           fontname='arial', fontsize=30, fontweight='bold')

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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/mxlkara_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
