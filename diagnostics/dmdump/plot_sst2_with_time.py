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
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/19??/d01/T/'
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
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask)

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    meansst2 = np.ndarray(shape=[len(tfiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(tfiles)):
        print tfiles[k]
        # Load & mask the current SST field.
        sst = nemo.load_field('sosstsst', '', '', tfiles[k], 'T')
        sst = nemo.mask_field(sst, tmask)

        # Calculate the area average surface temperature squared.
        meansst2[k, 0] = tfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_T', 1)[0]
        meansst2[k, 1] = (sst*sst*area).sum() / area.sum()

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meansst2']), meansst2)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    meansst2 = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/meansst2.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, 5*len(meansst2), len(meansst2)+1),
             np.concatenate(([0.], meansst2[:, 1]), axis=0),
             '-b', linewidth=3)
ax = plt.gca()
ax.set_aspect(5*len(meansst2)/250.0)


plt.axis([0, 5*len(meansst2), 0, 250.0])
plt.xticks(np.arange(0, 5*len(meansst2)+1, 365))
plt.yticks(np.linspace(0., 250.0, 5))
plt.xlabel('Days', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Ave. SST^2 ($\mathrm{^{\circ}C}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/meansst2_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
