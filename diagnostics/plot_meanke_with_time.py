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

import os
import xarray as xr


plt.ioff()

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n02/shared/chbull/'
nemodir = 'NEMO_JRAspinup4DAVECHK/'

# Specify the names of the different files that I want to load from.
# tdir = 'u-bl504/onm.nc.file/'
# udir = 'u-bl504/onm.nc.file/'
# vdir = 'u-bl504/onm.nc.file/'
tdir = 'u-bc337-bl504_links/'
udir = 'u-bc337-bl504_links/'
vdir = 'u-bc337-bl504_links/'

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
umask = np.squeeze(nemo.load_field('umask', homedir, nemodir, gridfile, 'U', nx, ny))
vmask = np.squeeze(nemo.load_field('vmask', homedir, nemodir, gridfile, 'V', nx, ny))

# Calculate the surface area of the T grid boxes & mask.
area = np.ma.masked_array(e1t * e2t, mask=1.0-tmask[:, :, 0:1])

# --------------------------------------------------------------------------- #

# Find the number of files in the directory that we want to calculate KE for.
if save_output:
    # tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, 'nemo_bl504o_1m_20??????-20??????_grid-T.nc'])))
    # ufiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, 'nemo_bl504o_1m_20??????-20??????_grid-U.nc'])))
    # vfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, 'nemo_bl504o_1m_20??????-20??????_grid-V.nc'])))
    tfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-T.nc'])))
    ufiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-U.nc'])))
    vfiles = sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-V.nc'])))

# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.

if save_output:
    # Preallocate the output variable.
    meanke = np.ndarray(shape=[len(ufiles), nz])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        print(ufiles[k])
        # Load & mask the current UU velocity field.
        uu = np.squeeze(nemo.load_field('vozocrtx2', '', '', ufiles[k], 'U', nx, ny))
        uu = nemo.mask_field(uu, umask)

        # Load and mask the current VV velocity field.
        vv = np.squeeze(nemo.load_field('vomecrty2', '', '', vfiles[k], 'V', nx, ny))
        vv = nemo.mask_field(vv, vmask)

        # Calculate the KE for the current velocity fields
        ke = nemo.calc_ke(uu.data, vv.data, nx, ny, nz, tmask)

        # Load & mask the current e3t field.
        e3t = np.squeeze(nemo.load_field('e3t', '', '', tfiles[k], 'T', nx, ny))
        e3t = nemo.mask_field(e3t, tmask)

        # Calculate the area average surface temperature.
        meanke[k, :] = np.sum(ke*e3t*area, axis=(0, 1)) / np.sum(e3t*area, axis=(0, 1))

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'post/meankine']), meanke)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    meanke = np.load(''.join([homedir, nemodir, 'post/meankine.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

#was good for 5-day output (DM's beast I think has 5-day output?)
idx=np.linspace(2.5, 5.0*len(meanke)-2.5, len(meanke))/365.0

#for monthly output..
idx=np.linspace(15.0, 30.0*len(meanke)-15.0, len(meanke))/365.0

ifile=xr.open_dataset(sorted(glob.glob(''.join([homedir, nemodir, tdir, '*_grid-U.nc'])))[0])
depths=ifile['depthu'].values.tolist()
ifile.close()

plt.close('all')
fig=plt.figure(figsize=(12, 12))
ax=fig.add_subplot(1, 1,1)

ax.plot(idx,meanke[:, 1], '-b',label= str(np.round(depths[1],1))+'m',linewidth=1)
ax.plot(idx,meanke[:, 17], '-g',label=str(np.round(depths[17],1))+'m',linewidth=1)
ax.plot(idx,meanke[:, 30], '-r',label=str(np.round(depths[30],1))+'m',linewidth=1)
ax.plot(idx,meanke[:, 50], '-m',label=str(np.round(depths[50],1))+'m',linewidth=1)
# ax = plt.gca()
ax.set_aspect(40.0/0.025)

ax.set_xlim([0, 65.0])
ax.set_ylim([0, 0.025])
ax.legend()

plt.xticks(np.arange(0, 96.0, 5.))
ax.set_xticks(np.arange(0, 96.0, 1.), minor=True)

plt.yticks(np.linspace(0., 0.0250, 5))

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

ax.set_xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
ax.set_ylabel('Ave. surface KE ($\mathrm{m^2s^{-1}}$)',fontname='arial', fontsize=30, fontweight='bold')
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
fig.savefig(''.join([homedir,nemodir,'figs/meankine_timeseries.pdf']), bbox_inches='tight')
fig.savefig(''.join([homedir,nemodir,'figs/meankine_timeseries.png']),dpi=300,bbox_inches='tight')
print("plot in: "+''.join([homedir,nemodir,'figs/meankine_timeseries.pdf']))
print("plot in: "+''.join([homedir,nemodir,'figs/meankine_timeseries.png']))

# --------------------------------------------------------------------------- #
