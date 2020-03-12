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
idir = 'TIDY/ARCHIVE/????/d01/I/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load volumes.
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

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
ifiles = sorted(glob.glob(''.join([homedir, nemodir, idir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the ice volume.

if save_output:
    # Preallocate the output variable.
    thickness = np.ndarray(shape=[len(ifiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ifiles)):
        print(ifiles[k])
        # Load & mask the current U velocity field.
        icethick = nemo.load_field('sithic', '', '', ifiles[k], 'T')
        icethick = nemo.mask_field(icethick, tmask)
        iceconc = nemo.load_field('siconc', '', '', ifiles[k], 'T')
        iceconc = nemo.mask_field(iceconc, tmask)

        # Calculate the KE for the current velocity fields
        thick = (e1t * e2t * icethick * iceconc).sum() / (e1t * e2t * tmask).sum()

        # Calculate the area average surface KE.
        thickness[k, 0] = ifiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_I', 1)[0]
        thickness[k, 1] = thick

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/icethick']), thickness)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    thickness = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/icethick.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, len(thickness), len(thickness)+1)/365.0,
             np.concatenate(([0.], thickness[:, 1]), axis=0),
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/0.25)


plt.axis([0, 40.0, 0, 0.25])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.arange(0., 0.3, 0.05))
ax.set_yticks(np.arange(0, 0.26, 0.01), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. Ice Thickness ($\mathrm{m}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/icethick_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
