#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import glosst
import aviso

# Import required modules.
import cmocean
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Register the dense colormap.
plt.register_cmap(name='thermal', cmap=cmocean.cm.thermal)

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
glosstdir = 'MO-GLO-SST/BACKUP/'

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(glosst.load_field('lat', homedir, glosstdir, 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20140131.nc'))
lon = np.squeeze(glosst.load_field('lon', homedir, glosstdir, 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20140131.nc'))
lon, lat = np.meshgrid(lon[:, None], lat[:, None])
lon = lon.T
lat = lat.T

# --------------------------------------------------------------------------- #

sst = np.squeeze(glosst.load_field('analysed_sst', homedir, glosstdir, 'METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_20140101-20140131.nc'))

mask[:] = sst
mask[mask > 0] = 0

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/Users/munday/Documents/Data/'
glosstdir = 'AVISO/'

# --------------------------------------------------------------------------- #

# Load the coordinates.
lat = np.squeeze(glosst.load_field('lat', homedir, glosstdir, 'dataset-duacs-rep-global-merged-allsat-phy-l4-v3_20140101-20141231_adt.nc'))
lon = np.squeeze(glosst.load_field('lon', homedir, glosstdir, 'dataset-duacs-rep-global-merged-allsat-phy-l4-v3_20140101-20141231_adt.nc'))
lon, lat = np.meshgrid(lon[:, None], lat[:, None])
lon = lon.T
lat = lat.T

# --------------------------------------------------------------------------- #

adt = np.squeeze(aviso.load_field('adt', homedir, glosstdir, 'dataset-duacs-rep-global-merged-allsat-phy-l4-v3_20140101-20141231_adt.nc'))

mask[:] = sst
mask[mask > 0] = 0

# --------------------------------------------------------------------------- #
