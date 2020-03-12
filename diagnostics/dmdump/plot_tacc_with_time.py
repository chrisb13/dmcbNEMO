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

# Turn interactive plotting off
plt.ioff()

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/COAREJRA-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/????/d05/U/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load KE values.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the grid spacings.
e2u = nemo.load_field('e2u', homedir, nemodir, gridfile, 'T'
                      )[2599:2600, :, None]

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U')
                   )[2599:2600, :, :, None]

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    tacc = np.ndarray(shape=[len(ufiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        print(ufiles[k])
        bounds = [2599, 2600, 0, 1999, 0, 74]
        # Load & mask the current U velocity field.
        uoce = nemo.load_field('vozoce3u', '', '', ufiles[k], 'U')[2599:2600,:,:]
        uoce = nemo.mask_field(uoce, umask)

        # Calculate the Drake Passage transport.
        tacc[k, 0] = ufiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_U', 1)[0]
        tacc[k, 1] = (uoce*e2u*umask).sum()/1.E6

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/tacc']), tacc)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    tacc = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/tacc.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, 5*len(tacc), len(tacc)+1)/365.,
             np.concatenate(([0.], tacc[:, 1]), axis=0),
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
plt.ylabel('$T_{ACC}$ ($\mathrm{Sv}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/tacc_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
