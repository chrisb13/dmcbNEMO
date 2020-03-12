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

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/19??/d05/T/'
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
e3t_0 = np.squeeze(nemo.load_field('e3t_0', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #

if save_output:
    # Preallocate the output variable.
    vvl = np.ndarray(shape=[len(tfiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(tfiles)):
        print tfiles[k]

        # Load & mask the current e3t field.
        e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T'))
        e3t = nemo.mask_field(e3t, tmask)

        # Calculate the volume of the ocean.
        volume = (e1t * e2t * e3t).sum()

        # Calculate the area average surface KE.
        vvl[k, 0] = tfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_T', 1)[0]
        vvl[k, 1] = volume

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/vvl']), vvl)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    vvl = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/vvl.npy']))

# --------------------------------------------------------------------------- #
# Load the annual average e3t fields.

# Find the number of files in the directory that we want to calculate volume for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir,
                                   'TIDY/ARCHIVE/VVL/', 'e3t_ANN_*'])))

# Preallocate the output variable.
ann_vvl = np.ndarray(shape=[len(tfiles), 2])

# Loop over the e3t files and load the volume.
for k in range(len(tfiles)):
    print tfiles[k]
    # Load & mask the current U velocity field.
    e3t = np.squeeze(np.load(''.join(['', '', tfiles[k]])))
    e3t = nemo.mask_field(e3t, tmask)

    # Calculate the volume of the ocean.
    volume = (e1t * e2t * e3t).sum()

    # Calculate the total annual volume.
    ann_vvl[k, 0] = k
    ann_vvl[k, 1] = volume

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(5, 5*len(vvl), len(vvl))/365.0,
             vvl[:, 1], '-b',
             (ann_vvl[:, 0] + 0.5), ann_vvl[:, 1], '.k',
             linewidth=1, markersize=20)
ax = plt.gca()
ax.set_aspect(40.0/15000000000000.0)


plt.axis([0, 40.0, 5.0567E17, 5.05685E17])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Days', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Domain volume from e3t ($\mathrm{m}^{3}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/vvl_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #