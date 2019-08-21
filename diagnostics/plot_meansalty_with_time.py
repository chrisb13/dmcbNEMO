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
import xarray as xr

plt.ioff()

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n02/shared/chbull/'
nemodir = 'NEMO_JRAspinup4DAVECHK/'

# Specify the names of the different files that I want to load from.
# tdir = 'u-bl504/onm.nc.file/'
tdir = 'u-bc337-bl504_links/'
gridfile = 'mesh_mask_eORCA025-GO7.nc'

# Specify the number of grid boxes.
nx = 1440
ny = 1207
nz = 75

# Choose whether to save/load KE values.
save_output = 0

# --------------------------------------------------------------------------- #

# Load the grid spacings.
e1t = nemo.load_field('e1t', homedir, nemodir, gridfile, 'T', nx, ny)
e2t = nemo.load_field('e2t', homedir, nemodir, gridfile, 'T', nx, ny)

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T', nx, ny))

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    # tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, 'nemo_bl504o_1m_20??????-20??????_grid-T.nc'])))
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-T.nc'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    meansalty = np.ndarray(shape=[len(tfiles), nz])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(tfiles)):
        print(tfiles[k])
        # Load & mask the current theta field.
        salty = np.squeeze(nemo.load_field('vosaline', '', '', tfiles[k], 'T', nx, ny))
        salty = nemo.mask_field(salty, tmask)

        # Load & mask the current e3t field.
        e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T', nx, ny))
        e3t = nemo.mask_field(e3t, tmask)

        # Calculate the area average salt
        meansalty[k, :] = np.sum(salty*e3t*area, axis=(0, 1)) / np.sum(e3t*area, axis=(0, 1))

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'post/meansalty_m01']), meansalty)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    meansalty = np.load(''.join([homedir, nemodir, 'post/meansalty_m01.npy']))

# --------------------------------------------------------------------------- #
ifile=xr.open_dataset(sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-T.nc'])))[0])
depths=ifile['deptht'].values.tolist()
ifile.close()

plt.rc('axes', linewidth=3)

plt.close('all')
fig=plt.figure(figsize=(12, 12))
ax=fig.add_subplot(1, 1,1)

idx=np.linspace(15.0, 30.0*len(meansalty)-15.0, len(meansalty))/365.0

ax.plot(idx,meansalty[:, 0], '-b',label= str(np.round(depths[0],1))+'m',linewidth=1)
ax.plot(idx,meansalty[:, 16], '-g',label=str(np.round(depths[16],1))+'m',linewidth=1)
ax.plot(idx,meansalty[:, 29], '-r',label=str(np.round(depths[29],1))+'m',linewidth=1)
ax.plot(idx,meansalty[:, 49], '-m',label=str(np.round(depths[49],1))+'m',linewidth=1)
ax.plot(idx,meansalty[:, 60], '-y',label=str(np.round(depths[60],1))+'m',linewidth=1)
ax.plot(idx,meansalty[:, 74], '-k',label=str(np.round(depths[74],1))+'m',linewidth=1)

# ax = plt.gca()
ax.set_aspect(40.0/1.0)

ax.set_xlim([0, 65.0])
ax.set_ylim([34.0, 35.00])
ax.legend()

plt.xticks(np.arange(0, 66.0, 5.))
ax.set_xticks(np.arange(0, 66.0, 1.), minor=True)

plt.yticks(np.linspace(34.0, 35.0, 5))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel(r'Ave. S ($\mathrm{g/kg}$)',
           fontname='arial', fontsize=30, fontweight='bold')

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(18)
    tick.label1.set_fontname('arial')
    tick.label1.set_fontweight('bold')

ax.grid(True)
# --------------------------------------------------------------------------- #
# Save the figure to a pdf.

fig.savefig(''.join([homedir,nemodir,'figs/meansalty_m01_timeseries.pdf']), bbox_inches='tight')
fig.savefig(''.join([homedir,nemodir,'figs/meansalty_m01_timeseries.png']),dpi=300,bbox_inches='tight')
print("plot in: "+''.join([homedir,nemodir,'figs/meansalty_m01_timeseries.pdf']))
print("plot in: "+''.join([homedir,nemodir,'figs/meansalty_m01_timeseries.png']))

# --------------------------------------------------------------------------- #
