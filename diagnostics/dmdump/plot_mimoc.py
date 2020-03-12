#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# import cmocean

# Register the thermal colormap.
# plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/'
nemodir = 'Data/MIMOC_ML_v2.2_CT_SA/'

# Specify the names of the different files that I want to load from.
filename = 'MIMOC_ML_v2.2_CT_SA_MLP_month08.nc'

# --------------------------------------------------------------------------- #
gphit = np.expand_dims(np.squeeze(
        nemo.load_field('LATITUDE', homedir, nemodir, filename)), 1)
glamt = np.expand_dims(np.squeeze(
        nemo.load_field('LONGITUDE', homedir, nemodir, filename)), 1)


# Load the coordinates.
gphit = np.expand_dims(np.squeeze(
        nemo.load_field('LATITUDE', homedir, nemodir, filename)), 1)
glamt = np.expand_dims(np.squeeze(
        nemo.load_field('LONGITUDE', homedir, nemodir, filename)), 1)

glamt, gphit = np.meshgrid(glamt, gphit)

# Load the field to plot and remask it.
mxl = nemo.load_field('DEPTH_MIXED_LAYER', homedir, nemodir, filename)
mxlmask = np.isfinite(mxl)
mxl[np.isnan(mxl)] = 0.0
mxl = nemo.mask_field(mxl, mxlmask)

# --------------------------------------------------------------------------- #

# Set the max mxl that I want to plot.
max_mxl = 500.

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))

# Specify some things to make the plot look nice.
map = Basemap(projection='spstere', boundinglat=-35, lon_0=-150, round='true')
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.5, color='w')
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True],
                  linewidth=0.5, color='w')
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_mxl = map.contourf(glamt, gphit, mxl.T,
                      np.arange(0., max_mxl+1, max_mxl/51),
                      latlon='true', cmap='seismic', vmin=0., vmax=max_mxl,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_mxl, 'right', ticks=np.arange(0., max_mxl+1, max_mxl/4),
                    pad=0.75)

# Save a hires version of the figure
# plt.savefig(''.join([filename[:-2], 'png']), bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
