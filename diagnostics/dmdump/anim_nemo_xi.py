#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import cmocean

# Register the balance colormap.
plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

# Change the backend for a reason I forget.
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
glamf = np.squeeze(nemo.load_field('glamf', homedir, nemodir, gridfile, 'Z'))
gphif = np.squeeze(nemo.load_field('gphif', homedir, nemodir, gridfile, 'Z'))

# Load the grid spacings.
e1u = nemo.load_field('e1u', homedir, nemodir, gridfile, 'U')
e2v = nemo.load_field('e2u', homedir, nemodir, gridfile, 'V')
e1f = nemo.load_field('e1f', homedir, nemodir, gridfile, 'Z')
e2f = nemo.load_field('e2f', homedir, nemodir, gridfile, 'Z')

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask',
                                   homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'V'))[:, :, 0:1]
fmask = np.squeeze(nemo.load_field('vmask',
                                   homedir, nemodir, gridfile, 'Z'))[:, :, 0:1]

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*'])))
vfiles = sorted(glob.glob(''.join([homedir, nemodir, vdir, '*'])))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and load the surface velocity.
for k in range(0, len(ufiles)+1, 5):
    print(ufiles[k])

    # Load the velocity fields and remask them.
    ssu = nemo.load_field('sossussu', '', '', ufiles[k], 'U')
    ssu = nemo.mask_field(ssu, umask)
    ssv = nemo.load_field('sossvssv', '', '', vfiles[k], 'V')
    ssv = nemo.mask_field(ssv, vmask)

    # Calculatet the surface vorticity.
    ssz = nemo.calc_xi(
        ssu.data, ssv.data, nx, ny, 1, e1u, e2v, e1f, e2f, fmask)

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
    cs_ssz = map.contourf(glamf.T, gphif.T, np.squeeze(ssz).T,
                          np.arange(-3.05E-5, 3.05E-5, 1.0E-6),
                          latlon='true', cmap='balance',
                          vmin=-3.E-5, vmax=3.E-5, extend='both')

    # Add a colour bar.
    cbar = map.colorbar(cs_ssz, 'right',
                        ticks=np.arange(-3.E-5, 4.E-5, 1.E-5),
                        pad=0.5, format='%.0e')
    cbar.ax.tick_params(labelsize=4) 
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(1)
        tick.label1.set_fontname('arial')
        tick.label1.set_fontweight('bold')

    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/',
                         nemodir[-29:-1], '/anim/xi/xi_', str(k).zfill(5), '.png']),
                bbox_inches='tight', dpi=300)

    # Close the current plot window.
    plt.close()

# --------------------------------------------------------------------------- #
