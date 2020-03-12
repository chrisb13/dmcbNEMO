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
nemodir = 'trunk/NEMOGCM/CONFIG/COAREJRA-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/????/d01/U/'
vdir = 'TIDY/ARCHIVE/????/d01/V/'
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
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask)

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
    totalssp = np.ndarray(shape=[len(ufiles), 2])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        print(ufiles[k])
        # Load & mask the current U velocity field.
        ssut = nemo.load_field('sossutau', '', '', ufiles[k], 'U')
        ssut = nemo.mask_field(ssut, umask)

        # Load and mask the current V velocity field.
        ssvt = nemo.load_field('sossvtau', '', '', vfiles[k], 'V')
        ssvt = nemo.mask_field(ssvt, vmask)

        # Calculate the KE for the current velocity fields
        ssp = nemo.calc_uv_on_t(ssut.data, ssvt.data, nx, ny, 1, tmask)

        # Calculate the area average surface KE.
        totalssp[k, 0] = ufiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_U', 1)[0]
        totalssp[k, 1] = (ssp*area).sum()

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/totalssp']), totalssp)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    totalssp = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/totalssp.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, len(totalssp), len(totalssp)+1)/365.0,
             np.concatenate(([0.], totalssp[:, 1]), axis=0)/1.E12,
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(40.0/2.0)


plt.axis([0, 40.0, 0, 2.0])

plt.xticks(np.arange(0, 41.0, 5.))
ax.set_xticks(np.arange(0, 41.0, 1.), minor=True)

plt.yticks(np.linspace(0., 2.0, 5))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Power Input ($\mathrm{m^2s^{-1}}$)',
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/totalssp_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
