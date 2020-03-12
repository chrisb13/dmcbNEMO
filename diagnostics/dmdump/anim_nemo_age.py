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
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55ABS-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/????/d05/T/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Specify an offset so that the animation can start partway through a run.
no_years = 40.0
offset = 73.0 * no_years

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T')[:, :, 30:31])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(tfiles)+1, 1):
    print tfiles[k]

    # Load the velocity fields and remask them.
    age = np.squeeze(nemo.load_field('voagee3t', '', '', tfiles[k], 'T')[:, :, 30:31])
    age = nemo.mask_field(age, tmask)
    e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T')[:, :, 30:31])
    e3t = nemo.mask_field(e3t, tmask)

    # Factor out the layer thickness.
    age = age / e3t
    e3t = None

    # Create a new figure window.
    plt.figure(figsize=(4.0, 4.0))

    # Specify some things to make the plot look nice.
    map = Basemap(projection='spaeqd',
                  boundinglat=-35, lon_0=180, round='true')

    # Draw grid lines and label the longitudes.
    map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
    map.drawmeridians(np.arange(-180, 180, 30), fontsize=4,
                      labels=12*[True], linewidth=0.25)
    map.drawmapboundary(fill_color='black')

    # Draw the contour plot.
    cs_age = map.contourf(glamt.T, gphit.T, np.squeeze(age).T,
                          np.arange(0., (k+offset+2.0) / 73.0, 1.0 / 73.0),
                          latlon='true', cmap='viridis', vmin=0.0, vmax=(k+offset+1.0)/73.0)

    # Add a colour bar.
    # cbar = map.colorbar(cs_age, 'right', pad=0.5)

    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(1)
        tick.label1.set_fontname('arial')
        tick.label1.set_fontweight('bold')

    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/',
                         nemodir[-29:-1], '/anim/age_nocb/age_', str(k).zfill(5), '.png']),
                bbox_inches='tight', dpi=300)

    # Close the current plot window.
    plt.close()

# --------------------------------------------------------------------------- #
