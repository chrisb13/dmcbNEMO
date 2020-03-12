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
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the names of the different files that I want to load from.
vdir = 'TIDY/ARCHIVE/1987/m01/V/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Choose whether to save/load KE values.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the relevant mask.
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))

# Load the grid spacings.
e2v = nemo.load_field('e2v', homedir, nemodir, gridfile, 'V')
e2v = nemo.mask_field(e2v, vmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    amoc = np.ndarray(shape=[len(vfiles), 2])

    # Loop over the V files and load Ve3v field.
    for k in range(len(vfiles)):
        print vfiles[k]

        # Load and mask the current V velocity field.
        ve3v = np.squeeze(nemo.load_field('vomece3v', '', '', vfiles[k], 'V'))
        ve3v = nemo.mask_field(ve3v, vmask)

        # Integrate the V velocity zonally and then cumsum from the bottom up.
        ve3v_integrated = (ve3v[2650:3679, :, :] * e2v[2650:3679:, :, :]).sum(axis=0).cumsum(axis=1)

        # Calculate the area average surface KE.
        amoc[k, 0] = vfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_V', 1)[0]
        amoc[k, 1] = ve3v_integrated[1768:1769,45:46]

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/amoc']), amoc)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    amoc = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/amoc.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(15.0, 30.0*len(amoc)-15.0, len(amoc)),
             amoc[:, 1]/1.E6, '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(360.0/30.0)


plt.axis([0, 360.0, 0, 30.0])

plt.yticks(np.arange(0, 31.0, 5.))
ax.set_yticks(np.arange(0, 31.0, 1.), minor=True)

plt.xticks(np.linspace(0., 360.0, 13.0))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Days', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('AMOC @ ~26$^{\circ}$S',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/amoc_26oS_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
