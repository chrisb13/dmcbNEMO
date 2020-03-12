#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt

# Register the speed colormap.
plt.register_cmap(name='balance', cmap=cmocean.cm.balance)

# Allow plotting without a display.
plt.switch_backend('agg')

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
#! homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP01/'

# Specify the years to do the averaging over and the season to average.
init_year = 1983 # Sets the output file names.
num_years = 5

# Specify the names of the different files that I want to load from.
ufilename = ''.join(['TIDY/ARCHIVE/ORCH0083-LIM3_', str(init_year), '-', str(init_year + num_years - 1), '_U_d01.nc'])
vfilename = ''.join(['TIDY/ARCHIVE/ORCH0083-LIM3_', str(init_year), '-', str(init_year + num_years - 1), '_V_d01.nc'])
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# --------------------------------------------------------------------------- #

# Load the coordinates.
glamt = np.squeeze(nemo.load_field('glamt', homedir, nemodir, gridfile, 'T'))
gphit = np.squeeze(nemo.load_field('gphit', homedir, nemodir, gridfile, 'T'))
glamu = np.squeeze(nemo.load_field('glamu', homedir, nemodir, gridfile, 'U'))
gphiu = np.squeeze(nemo.load_field('gphiu', homedir, nemodir, gridfile, 'U'))
glamv = np.squeeze(nemo.load_field('glamv', homedir, nemodir, gridfile, 'V'))
gphiv = np.squeeze(nemo.load_field('gphiv', homedir, nemodir, gridfile, 'V'))

# Set the longitudes to go from -180o to 180o starting on the left.
glamt[glamt<=glamt[-1,1000]] = glamt[glamt<=glamt[-1,1000]] + 360.
glamt = glamt-glamt.min()
glamt = glamt - 180.

glamu[glamu<=glamu[-1,1000]] = glamu[glamu<=glamu[-1,1000]] + 360.
glamu = glamu-glamu.min()
glamu = glamu - 180.

glamv[glamv<=glamv[-1,1000]] = glamv[glamv<=glamv[-1,1000]] + 360.
glamv = glamv-glamv.min()
glamv = glamv - 180.

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T'))[:, :, 0:1]
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U'))[:, :, 0:1]
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V'))[:, :, 0:1]

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T')
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T')
e1u = nemo.load_field('e1t', homedir, nemodir, gridfile, 'U')
e2u = nemo.load_field('e2t', homedir, nemodir, gridfile, 'U')
e1v = nemo.load_field('e1t', homedir, nemodir, gridfile, 'V')
e2v = nemo.load_field('e2t', homedir, nemodir, gridfile, 'V')

# --------------------------------------------------------------------------- #

# Load the initial 1d level thicknesses.
gdept = np.squeeze(nemo.load_field('gdept_1d', homedir, nemodir, gridfile))

# Load the array with the index containing the deepest level.
mbathy = np.squeeze(nemo.load_field('mbathy', homedir, nemodir, gridfile, 'T'))

# Preallocate the depth array.
depth = np.ndarray((nx, ny))

# Create the bathymetry map.
for i in np.arange(nx):
    for j in np.arange(ny):
        depth[i, j] = gdept[mbathy[i, j]].sum()

# --------------------------------------------------------------------------- #

# Load the power input fields and remask them.
utaux = nemo.load_field('sossutau', homedir, nemodir, ufilename, 'U')
utaux = nemo.mask_field(utaux, umask)
vtauy = nemo.load_field('sossvtau', homedir, nemodir, vfilename, 'V')
vtauy = nemo.mask_field(vtauy, vmask)

# Load the velocity fields and remask them.
ssu = nemo.load_field('sossussu', homedir, nemodir, ufilename, 'U')
ssu = nemo.mask_field(ssu, umask)
ssv = nemo.load_field('sossvssv', homedir, nemodir, vfilename, 'V')
ssv = nemo.mask_field(ssv, vmask)

# Load the wind fields and remask them.
taux = nemo.load_field('sozotaux', homedir, nemodir, ufilename, 'U')
taux = nemo.mask_field(taux, umask)
tauy = nemo.load_field('sometauy', homedir, nemodir, vfilename, 'V')
tauy = nemo.mask_field(tauy, vmask)

# --------------------------------------------------------------------------- #

total_power, mean_power, eddy_power = nemo.calc_wind_power(utaux.data,
                                                           vtauy.data,
                                                           taux.data,
                                                           tauy.data,
                                                           ssu.data, ssv.data,
                                                           nx, ny, 1, tmask)

# --------------------------------------------------------------------------- #

# Plot the total power input.
ax, cbar = nemo.plot_tfield(np.squeeze(total_power), np.squeeze(gphit), np.squeeze(glamt), -0.05, 0.05, 20, 'balance', 'both', gridcolor='black')

# Save a hires version of the figure
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_total_', str(init_year), '-', str(init_year + num_years - 1), '.png']), bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

# Plot the mean power input.
ax, cbar = nemo.plot_tfield(np.squeeze(mean_power), np.squeeze(gphit), np.squeeze(glamt), -0.05, 0.05, 20, 'balance', 'both', gridcolor='black')

# Save a hires version of the figure
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_mean_', str(init_year), '-', str(init_year + num_years - 1), '.png']), bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

# Plot the eddy power input.
ax, cbar = nemo.plot_tfield(np.squeeze(eddy_power), np.squeeze(gphit), np.squeeze(glamt), -0.01, 0.01, 20, 'balance', 'both', gridcolor='black')

# Save a hires version of the figure
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_eddy_', str(init_year), '-', str(init_year + num_years - 1), '.png']), bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

# Plot the mean zonal wind stress.
ax, cbar = nemo.plot_tfield(np.squeeze(taux[0:nx,:]), np.squeeze(gphiu[0:nx,:]), np.squeeze(glamt[0:nx,:]), -0.2, 0.2, 20, 'balance', 'both', gridcolor='black')

# Save a hires version of the figure
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_taux_', str(init_year), '-', str(init_year + num_years - 1), '.png']), bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

# Plot the mean meridional wind stress.
ax, cbar = nemo.plot_tfield(np.squeeze(tauy[:,0:ny]), np.squeeze(gphit), np.squeeze(glamt), -0.2, 0.2, 20, 'balance', 'both', gridcolor='black')

# Save a hires version of the figure
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_tauy_', str(init_year), '-', str(init_year + num_years - 1), '.png']), bbox_inches='tight', dpi=600)

# --------------------------------------------------------------------------- #

print ('')
print ('Integral over entire domain')
print (' ')

# Integrate the power input over the surface of the ocean.
print ((total_power * e1t * e2t).sum()/1.E12)
print ((mean_power * e1t * e2t).sum()/1.E12)
print ((eddy_power * e1t * e2t).sum()/1.E12)

# Setup the new mask.
mask = tmask.copy()
mask[depth <= 1000.0] = 0.0

print (' ')
print ('Integral over depths > 1000m')
print (' ')

# Integrate the power input over the surface of the ocean.
print ((total_power * e1t * e2t * mask).sum()/1.E12)
print ((mean_power * e1t * e2t * mask).sum()/1.E12)
print ((eddy_power * e1t * e2t * mask).sum()/1.E12)

# --------------------------------------------------------------------------- #

# Integrate the power inputs zonally.
total_power = (total_power * e1t * tmask).sum(axis=0) / \
                 (e1t * tmask).sum(axis=0)
mean_power = (mean_power * e1t * tmask).sum(axis=0) / \
                 (e1t * tmask).sum(axis=0)
eddy_power = (eddy_power * e1t * tmask).sum(axis=0) / \
                 (e1t * tmask).sum(axis=0)

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12.0, 12.0))

plt.rc('axes', linewidth=3)

p = plt.plot(
             gphit.mean(axis=0), total_power, '-b',
             gphit.mean(axis=0), mean_power, '-g',
             gphit.mean(axis=0), eddy_power, '-r',
             linewidth=3)
ax = plt.gca()
ax.set_aspect(40.0/0.035)


plt.axis([-70.0, -30.0, -0.01, 0.02])
plt.xticks(np.arange(-70.0, -30.0, 5))
plt.yticks(np.arange(-0.01, 0.025, 0.005))
plt.xlabel('Latitude', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('Zonal Mean Power Input ($\mathrm{W}m^{-2}$)',
           fontname='arial', fontsize=30, fontweight='bold')

plt.legend(handles=[p[0], p[1], p[2]], labels=['Total', 'Mean', 'Eddy'])

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
plt.savefig(''.join([homedir, 'figs/', nemodir[-29:], 'ssp_wind_power_comparison_', str(init_year), '-', str(init_year + num_years - 1), '.pdf']), bbox_inches='tight')

# --------------------------------------------------------------------------- #
