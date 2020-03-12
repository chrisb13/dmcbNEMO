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
tdir = 'TIDY/ARCHIVE/????/d01/T/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Set physical constants.
g = 9.80665 # gravitational acceleration.

# Choose whether to save/load KE values.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2v = nemo.load_field('e2v', homedir, nemodir, gridfile, 'V')

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask)

# --------------------------------------------------------------------------- #
# Load the Coriolis parameter.
f = nemo.load_field('ff_f', homedir, nemodir, gridfile, 'Z')

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    geoke = np.ndarray(shape=[len(tfiles), 2])

    # Loop over the T files and load the SSH for geostrophic velocity calc.
    for k in range(len(tfiles)):
        print(tfiles[k])
        # Load & mask the current SSH fields.
        ssh = nemo.load_field('sossheig', '', '', tfiles[k], 'T')
        ssh = nemo.mask_field(ssh, tmask)

        # Calculate the geostrophic velocities.
        ssu = -2.0 * g * nemo.grad_dtdy(ssh.data, e2v, vmask, nx, ny, 1) / (f[0:nx, :, :] + f[1:nx+1, :, :])
        ssu = ssu * vmask
        ssv = 2.0 * g * nemo.grad_dtdx(ssh.data, e1u, umask, nx, ny, 1)/ (f[:, 0:ny, :] + f[:, 1:ny+1, :])
        ssv = ssv * umask

        # Calculate the geostrophic KE for the current geostrophic velocity fields
        ssk = nemo.calc_ke(ssv, ssu, nx, ny, 1, tmask)

        # Calculate the area average surface geostrophic KE.
        geoke[k, 0] = tfiles[k].split('/ORCH0083-LIM3_', 1)[1].split('_T', 1)[0]
        geoke[k, 1] = (ssk*area).sum() / area.sum()

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/geoke']), geoke)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    geoke = np.load(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/geoke.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, len(geoke), len(geoke)+1)/365.0,
             np.concatenate(([0.], geoke[:, 1]), axis=0),
             '-b', linewidth=1)
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

plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/geoke_timeseries.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
