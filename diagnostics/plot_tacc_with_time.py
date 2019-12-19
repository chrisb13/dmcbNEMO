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

# --------------------------------------------------------------------------- #

# Turn interactive plotting off
plt.ioff()

#########
#  JRA  #
#########

# Specify where the input data lives.
homedir = '/nerc/n02/shared/chbull/'
nemodir = 'NEMO_JRAspinup4DAVECHK/'

# Specify the names of the different files that I want to load from.
udir = 'u-bl504/onm.nc.file/'
udir = 'u-bc337-bl504_links/'
gridfile = 'mesh_mask_eORCA025-GO7.nc'

##########
#  CORE  #
##########

# Specify where the input data lives for CORE.
homedir = '/nerc/n02/n02/chbull/RawData/'
nemodir = 'mi-an854/'

# Specify the names of the different files that I want to load from.
udir = 'onm.nc.file/'
gridfile = 'mesh_mask_eORCA025-GO7.nc'

# Specify the number of grid boxes.
nx = 1440
ny = 1207
nz = 75

# make output dir
nemo.mkdir(homedir+nemodir+'figs/')
nemo.mkdir(homedir+nemodir+'post/')

# Choose whether to save/load tacc values.
save_output = 1

# --------------------------------------------------------------------------- #

# Load the grid spacings.
e2u = nemo.load_field('e2u', homedir, nemodir, gridfile, 'U')[885:886, 315:421, None]

# Load the relevant mask.
umask = np.squeeze(nemo.load_field('umask',homedir, nemodir, gridfile, 'U'))[885:886, 315:421, None]

# --------------------------------------------------------------------------- #

if save_output:
    # Find the number of files in the directory that we want to calculate KE for.
    # ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, 'nemo_bl504o_1m_20??????-20??????_grid-U.nc'])))
    # ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*_grid-U.nc'])))

    #for CORE
    ufiles = sorted(glob.glob(''.join([homedir, nemodir, udir, '*_????????_grid_U.nc'])))


# --------------------------------------------------------------------------- #
# If we're not loading the data, then loop over all the available files and
# calculatet the average KE.


if save_output:
    # Preallocate the output variable.
    tacc = np.ndarray(shape=[len(ufiles), 1])

    # Loop over the U/V files and load the surface velocity.
    for k in range(len(ufiles)):
        # Load & mask the current U velocity field.
        uoce = nemo.load_field('vozocrtx', '', '', ufiles[k], 'U')[885:886, 315:421, :]
        e3u = nemo.load_field('e3u', '', '', ufiles[k], 'U')[885:886, 315:421, :]
        # uoce = nemo.mask_field(uoce*e3u, umask)
        uoce = uoce*e3u

        # found that the multiplication with the umask was leading to crazy values - ~300 Sv; so turned it off

        # Calculate the Drake Passage transport.
        tacc[k, 0] = (uoce*e2u).sum()/1.E6
        # print((uoce*e2u).sum()/1.E6)
        print(ufiles[k],str(np.round(tacc[k, 0],1)))

# --------------------------------------------------------------------------- #
# Spit the numbers out to file, if requested.
if save_output:
    np.save(''.join([homedir, nemodir, 'post/tacc']), tacc)

# --------------------------------------------------------------------------- #
# If we're not saving the number, we're loading them.
if not save_output:
    tacc = np.load(''.join([homedir, nemodir, 'post/tacc.npy']))

# --------------------------------------------------------------------------- #

plt.figure(figsize=(12, 12))

plt.rc('axes', linewidth=3)

p = plt.plot(np.linspace(0, 30.0*len(tacc), len(tacc)+1)/365.,
             np.concatenate(([0.], tacc[:,0]), axis=0),
             '-b', linewidth=1)
ax = plt.gca()
ax.set_aspect(50.0/300.0)

plt.axis([0, 50.0, 0, 350.0])

plt.xticks(np.arange(0, 91.0, 5.))
ax.set_xticks(np.arange(0, 91.0, 1.), minor=True)

plt.yticks(np.arange(0., 351.0, 50.))
ax.set_yticks(np.arange(0, 351.0, 10.), minor=True)

ax.tick_params(which='major', length=10, width=2, direction='in')
ax.tick_params(which='minor', length=5, width=2, direction='in')

plt.xlabel('Years', fontname='arial', fontsize=30, fontweight='bold')
plt.ylabel('$T_{ACC}$ ($\mathrm{Sv}$)',
           fontname='arial', fontsize=30, fontweight='bold')

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

# plt.savefig(''.join(['../figs/tacc_timeseries.pdf']), bbox_inches='tight')

plt.savefig(''.join([homedir,nemodir,'figs/tacc_timeseries.pdf']), bbox_inches='tight')
plt.savefig(''.join([homedir,nemodir,'figs/tacc_timeseries.png']),dpi=300,bbox_inches='tight')
print("plot in: "+''.join([homedir,nemodir,'figs/tacc_timeseries.pdf']))
print("plot in: "+''.join([homedir,nemodir,'figs/tacc_timeseries.png']))

# --------------------------------------------------------------------------- #
