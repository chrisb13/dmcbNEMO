#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# Register the dense colormap.
plt.register_cmap(name='deep', cmap=cmocean.cm.deep)

plt.switch_backend('agg')

# --------------------------------------------------------------------------- #

# Turn interactive plotting off
plt.ioff()

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55TAU-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
tdir = 'TIDY/ARCHIVE/????/d01/T/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(tfiles)+1, 5):
    print(tfiles[k])

    # Load the velocity fields and remask them.
    mxl = nemo.load_field('somxlkara', '', '', tfiles[k], 'T')
    mxl = nemo.mask_field(mxl, tmask)

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
    cs_ice = map.contourf(glamt.T, gphit.T, np.squeeze(mxl).T,
                          np.arange(0., 510.0, 10.0),
                          latlon='true', cmap='deep', extend='max')

    # Add a colour bar.
    cbar = map.colorbar(cs_ice, 'right', ticks=np.arange(0., 550.0, 50.0),
                        pad=0.5)
    cbar.ax.tick_params(labelsize=4)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(1)
        tick.label1.set_fontname('arial')
        tick.label1.set_fontweight('bold')

    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/',
                         nemodir[-29:-1], '/anim/mxl/mxl_', str(k).zfill(5), '.png']),
                bbox_inches='tight', dpi=300)

    # Close the current plot window.
    plt.close()

# --------------------------------------------------------------------------- #
