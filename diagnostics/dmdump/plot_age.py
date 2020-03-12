#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
# import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.set_cmap(cmaps.viridis)

# Register the thermal colormap.
# plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
# homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
filename = 'TIDY/ARCHIVE/1952/d05/T/ORCH0083-LIM3_19521227_T_d05.nc'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# --------------------------------------------------------------------------- #

# Load the coordinates.
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gdept = np.squeeze(nemo.load_field('gdept_1d', homedir, nemodir, gridfile))

# Load the relevant mask & e3 field.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))

# Load the field to plot and remask it.
nage = np.squeeze(nemo.load_field('vonage3t', homedir, nemodir, filename, 'T'))
nage = nemo.mask_field(nage, tmask)
sage = np.squeeze(nemo.load_field('voagee3t', homedir, nemodir, filename, 'T'))
sage = nemo.mask_field(sage, tmask)
e3t = np.squeeze(nemo.load_field('e3t', homedir, nemodir, filename, 'T'))

# Divide by the cell depth to get units right.
nage = nage / e3t
sage = sage / e3t

# Set the maximum age of the water.
max_age = 2.0

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_sst = map.contourf(glamt.T, gphit.T, np.squeeze(sage[:, :, 30:31]).T,
                      np.arange(0., max_age+max_age/40, max_age/40),
                      latlon='true', cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(0., max_age+max_age/4., max_age/4.), pad=0.5)

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/ssage.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[-1:, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(sage[-1:, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/sage_ind.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[1750:1751, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(sage[1750:1751, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/sage_pac.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[3250:3251, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(sage[3250:3251, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/sage_atl.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)

# Specify some things to make the plot look nice.
map = Basemap(projection='spaeqd', boundinglat=-35,
              lon_0=180, round='true')

# Draw grid lines and label the longitudes.
map.drawparallels(np.arange(-80, 0, 20), linewidth=0.25)
map.drawmeridians(np.arange(-180, 180, 30), labels=12*[True], linewidth=0.25)
map.drawmapboundary(fill_color='black')

# Draw the contour plot.
cs_sst = map.contourf(glamt.T, gphit.T, np.squeeze(nage[:, :, 30:31]).T,
                      np.arange(0., max_age+max_age/40, max_age/40),
                      latlon='true', cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = map.colorbar(cs_sst, 'right', ticks=np.arange(0., max_age+max_age/4., max_age/4.), pad=0.5)

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/nsage.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[-1:, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(nage[-1:, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/nage_ind.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[1750:1751, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(nage[1750:1751, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/nage_pac.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #

# Create a new figure window.
plt.figure(figsize=(8.0, 8.0))
# Use: plt.savefig('sst.png', bbox_inches='tight', dpi=1200)


# Draw the contour plot.
cs_sst = plt.contourf(np.repeat(gphit[3250:3251, :].T, 75, axis=1),
                      -np.repeat(gdept.T[:, None], 2000, axis=1).T,
                      np.squeeze(nage[3250:3251, :, :]),
                      np.arange(0., max_age+max_age/40, max_age/40),
                      cmap='viridis', vmin=0., vmax=2.,
                      extend='both')

# Add a colour bar.
cbar = plt.colorbar(cs_sst, ticks=np.arange(0., max_age+max_age/4., max_age/4.))

# Save a hires version of the figure
plt.savefig('../figs/CORE2NYF-ORCH0083-LIM3/EXP00/nage_atl.png', bbox_inches='tight', dpi=1200)

# --------------------------------------------------------------------------- #
