#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# import colormaps as cmaps
import cmocean

# Register the viridis colormap.
# plt.register_cmap(name='viridis', cmap=cmaps.viridis)

# Register the speed colormap.
plt.register_cmap(name='speed', cmap=cmocean.cm.speed)

# Switch the backend for some reason I forget.
plt.switch_backend('agg')

# --------------------------------------------------------------------------- #

# Turn interactive plotting off
plt.ioff()

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55TAU-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
udir = 'TIDY/ARCHIVE/19??/d01/U/'
vdir = 'TIDY/ARCHIVE/19??/d01/V/'
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
tmask = np.squeeze(nemo.load_field('tmask',
                                   homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(ufiles)+1, 5):
    print(ufiles[k])

    # Load the velocity fields and remask them.
    ssu2 = nemo.load_field('sossusqu', '', '', ufiles[k], 'U')
    ssu2 = nemo.mask_field(ssu2, umask)
    ssv2 = nemo.load_field('sossvsqu', '', '', vfiles[k], 'V')
    ssv2 = nemo.mask_field(ssv2, vmask)

    # Calculate the surface vorticity.
    ssk = 0.5 * nemo.calc_uv_on_t(ssu2.data, ssv2.data, nx, ny, 1, tmask)

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
    cs_ssz = map.contourf(glamt.T, gphit.T, np.log10(1.E4*np.squeeze(ssk)).T,
                          np.arange(-0.5, 3.51, 0.1),
                          latlon='true', cmap='speed', vmin=-0.5, vmax=3.5,
                          extend='both')

    # Add a colour bar.
    cbar = map.colorbar(cs_ssz, 'right', ticks=np.arange(-0.5, 3.6, 0.5),
                        pad=0.5)
    cbar.ax.tick_params(labelsize=4) 
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(1)
        tick.label1.set_fontname('arial')
        tick.label1.set_fontweight('bold')

    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/',
                         nemodir[-29:-1], '/anim/ke/ke_', str(k).zfill(5), '.png']),
                bbox_inches='tight', dpi=300)

    # Close the current plot window.
    plt.close()

# --------------------------------------------------------------------------- #
