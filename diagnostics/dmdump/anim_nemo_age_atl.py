#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)

plt.switch_backend('agg')

# --------------------------------------------------------------------------- #

# Turn interactive plotting off
plt.ioff()

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/195[123456789]/d05/T/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Specify an offset so that the animation can start partway through a run.
no_years = 0.0
offset = 73.0 * no_years

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
gdept = np.squeeze(nemo.load_field('gdept_1d', homedir, nemodir, gridfile))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T')[3250:3251, :, :])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(tfiles)+1, 1):
    print tfiles[k]

    # Load the velocity fields and remask them.
    age = np.squeeze(nemo.load_field('voagee3t', '', '', tfiles[k], 'T')[3250:3251, :, :])
    age = nemo.mask_field(age, tmask)
    e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T')[3250:3251, :, :])
    e3t = nemo.mask_field(e3t, tmask)

    # Factor out the layer thickness.
    age = age / e3t
    e3t = None

    # Create a new figure window.
    plt.figure(figsize=(4.0, 4.0))

    # Draw the contour plot.
    cs_age = plt.contourf(np.repeat(gphit[3250:3251, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      age,
                      np.arange(0., (k+offset+2.0) / 73.0, 1.0 / 73.0),
                      cmap='viridis', vmin=0.0, vmax=(k+offset+1.0)/73.0)

    # Add a colour bar.
    cbar = plt.colorbar(cs_age)

    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(1)
        tick.label1.set_fontname('arial')
        tick.label1.set_fontweight('bold')

    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/',
                         nemodir[-29:-1], '/anim/age_atl/age_atl_', str(k).zfill(5), '.png']),
                bbox_inches='tight', dpi=300)

    # Close the current plot window.
    plt.close()

# --------------------------------------------------------------------------- #
