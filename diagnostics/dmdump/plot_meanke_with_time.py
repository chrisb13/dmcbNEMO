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
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2TAU-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/????/d05/U/'
vdir = 'TIDY/ARCHIVE/????/d05/V/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load KE values.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
    vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    meanke = np.ndarray(shape=[len(ufiles), nz+1])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        print ufiles[k]
        # Load & mask the current U velocity field.
        u = np.squeeze(nemo.load_field('vozoce3u', '', '', ufiles[k], 'U'))
        u = nemo.mask_field(u, umask)
        e3u = np.squeeze(nemo.load_field('e3u', '', '', ufiles[k], 'U'))
        e3u = nemo.mask_field(e3u, umask)
        u = u / e3u
        e3u = None

        # Load and mask the current V velocity field.
        v = np.squeeze(nemo.load_field('vomece3v', '', '', vfiles[k], 'V'))
        v = nemo.mask_field(v, vmask)
        e3v = np.squeeze(nemo.load_field('e3v', '', '', vfiles[k], 'V'))
        e3v = nemo.mask_field(e3v, vmask)
        v = v / e3v
        e3v = None

        # Calculate the KE for the current velocity fields
        ke = nemo.calc_ke(u.data, v.data, nx, ny, nz, tmask)

        # Calculate the area average surface temperature.
        meanke[k, 0] = ufiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_U', 1)[0]
        meanke[k, 1:] = np.sum(ke*area, axis=(0, 1)) / np.sum(area, axis=(0, 1))

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meankine']), meanke)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    meanke = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meankine.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(2.5, 5.0*len(meanke)-2.5, len(meanke))/365.0,
             meanke[:, 1], '-b',
             np.linspace(2.5, 5.0*len(meanke)-2.5, len(meanke))/365.0,
             meanke[:, 17], '-g',
             np.linspace(2.5, 5.0*len(meanke)-2.5, len(meanke))/365.0,
             meanke[:, 30], '-r',
             np.linspace(2.5, 5.0*len(meanke)-2.5, len(meanke))/365.0,
             meanke[:, 50], '-m',
             linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/0.025)

plt.axis([0, 40.0, 0, 0.025])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.linspace(0., 0.0250, 5))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. surface KE ($\mathrm{m^2s^{-1}}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/meankine_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
