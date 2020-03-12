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
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/????/m01/T/'
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

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    meansalty = np.ndarray(shape=[len(tfiles), nz+1])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(tfiles)):
        print tfiles[k]
        # Load & mask the current theta field.
        salty = np.squeeze(nemo.load_field('vosale3t', '', '', tfiles[k], 'T'))
        salty = nemo.mask_field(salty, tmask)

        # Load & mask the current e3t field.
        e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T'))
        e3t = nemo.mask_field(e3t, tmask)

        # Calculate the area average surface temperature.
        meansalty[k, 0] = tfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_T', 1)[0]
        meansalty[k, 1:] = np.sum(salty*area, axis=(0, 1)) / np.sum(e3t*area, axis=(0, 1))

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meansalty_m01']), meansalty)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    meansalty = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meansalty_m01.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(15.0, 30.0*len(meansalty)-15.0, len(meansalty))/365.0,
             meansalty[:, 1], '-b',
             np.linspace(15.0, 30.0*len(meansalty)-15.0, len(meansalty))/365.0,
             meansalty[:, 17], '-g',
             np.linspace(15.0, 30.0*len(meansalty)-15.0, len(meansalty))/365.0,
             meansalty[:, 30], '-r',
             np.linspace(15.0, 30.0*len(meansalty)-15.0, len(meansalty))/365.0,
             meansalty[:, 50], '-m',
             linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/0.25)

plt.axis([0, 40.0, 34.75, 35.0])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.linspace(34.75, 35.0, 5))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel(r'Ave. S ($\mathrm{g/kg}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/meansalty_m01_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
