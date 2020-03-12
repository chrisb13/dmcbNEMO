#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from scipy import ndimage
from skimage import measure

# Turn interactive plotting off
plt.switch_backend('agg')
plt.ioff()

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
nemodir = 'trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/EXP00/'

# Specify the names of the different files that I want to load from.
idir = 'TIDY/ARCHIVE/201[34567]/d01/I/'
gridfile = 'TIDY/ARCHIVE/MESH/mesh_mask.nc'

# Specify the number of grid boxes.
nx = 4320
ny = 2000
nz = 75

# Bound on the Marginal Ice Zone (MIZ).
miz_bound = 0.8

# --------------------------------------------------------------------------- #

# Load the relevant mask.
tmask = np.squeeze(nemo.load_field('tmask', homedir, nemodir, gridfile, 'T')[:, :, 0])

# Load the grid spacings.
e1t = np.squeeze(nemo.load_field('e1t', homedir, nemodir, gridfile, 'T'))
e2t = np.squeeze(nemo.load_field('e2t', homedir, nemodir, gridfile, 'T'))

# --------------------------------------------------------------------------- #

# Find the files that I want to classify the ice cover for.
ifiles = sorted(glob.glob(''.join([homedir, nemodir, idir, '*'])))

# --------------------------------------------------------------------------- #

# Setup the colormap - this makes sure that the colorbar shows 8 discrete colours
#                      instead of interpolating between them. It seems like 
#                      pointless effort for something that just works on my laptop.
base = plt.cm.get_cmap('Set2',8)
basemap = base(range(8))
basemap[0, :] = 0.0
basemap = ListedColormap(basemap)

# --------------------------------------------------------------------------- #

# Preallocate an array to store the area of each category..
areas = np.zeros((len(ifiles), 9))

# Loop over the U/V files and calculate the turning angle..
for k in range(len(ifiles)):
    print(ifiles[k])
 
    # Load the U velocity field and remask it.
    siconc = np.squeeze(nemo.load_field('siconc', '', '', ifiles[k], 'T'))
    siconc = nemo.mask_field(siconc, tmask)
       
    #fig = plt.figure(figsize=(12.0,12.0))
    #plt.subplot(211)
    #ax1 = plt.pcolormesh(np.squeeze(siconc).T)
    #plt.clim(0,1)
    #plt.axis([0, 500, 950, 1200])
    #cbar = fig.colorbar(ax1, ticks=[0.0, 0.15, 0.8, 1.0, 1.2])

    # Preallocate the ice cateogry array.
    ice_category = np.zeros((nx,ny), dtype=int)

    # Preliminary categorisation.
    ice_category[tmask < 0.001] = 0                           # land
    ice_category[siconc == 0.0] = 1                           # ocean
    ice_category[(siconc < 0.15) & (siconc > 0)] = 2          # loose ice
    ice_category[(siconc >= 0.15) & (siconc < miz_bound)] = 3 # MIZ
    ice_category[siconc >= miz_bound] = 4                     # pack ice
    
    # Preallocate the masks for each category.
    category_mask = np.zeros((nx,ny,5), dtype=int)
    
    # Derive the masks from ice_category.
    category_mask[(ice_category == 0) | (ice_category == 2) | (ice_category == 3), 0] = 1 # land           / Lland
    category_mask[(ice_category == 2), 1] = 1                                             # JUST loose ice mask
    category_mask[(ice_category == 1) | (ice_category == 2) | (ice_category == 3), 2] = 1 # loose ice mask / Ll
    category_mask[ice_category == 3, 3] = 1                                               # MIZ mask       / Lm
    category_mask[ice_category == 4, 4] = 1                                               # pack ice mask  / Lp
    
    #fig = plt.figure(figsize=(12.0,12.0))
    #ax2 = plt.pcolormesh(np.squeeze(ice_category).T, cmap='Set2')
    #plt.clim(0,8)
    #plt.axis([0, 500, 950, 1200])
    #cbar = fig.colorbar(ax2, ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    #cbar.ax.set_yticklabels(['land', 'ocean', 'loose ice', 'MIZ', 'pack ice', 'open pack ice', 'coastal polynyas', 'inner open water'])

    # Preallocate labelling array.
    labelled_array = np.zeros((nx,ny,5), dtype=int)
    num_features = np.zeros((5,1), dtype=int)
    
    # Use scipy.ndimage.label to find the features.
    for l in np.arange(0,5):
        labelled_array[:,:,l], num_features[l] = ndimage.label(category_mask[:,:,l], [[1,1,1],[1,1,1],[1,1,1]])

    # Find the perimeter length for the second stage of classification - test ocean + loose ice + MIZ.
    label_perimeter = measure.regionprops(labelled_array[:,:,2])
    label_perimeter = [prop.perimeter for prop in label_perimeter]
    
    perimeter = np.zeros((nx,ny,4))
    for l in np.arange(num_features[2]):
        perimeter[labelled_array[:,:,2]==l+1,0] = label_perimeter[l]

    # Find the perimeter length for the second stage of classification - test loose ice ONLY.
    label_perimeter = measure.regionprops(labelled_array[:,:,1])
    label_perimeter = [prop.perimeter for prop in label_perimeter]
    
    for l in np.arange(num_features[1]):
        perimeter[labelled_array[:,:,1]==l+1,1] = label_perimeter[l]

    # Find the perimeter length for the second stage of classification - test MIZ ONLY.
    label_perimeter = measure.regionprops(labelled_array[:,:,3])
    label_perimeter = [prop.perimeter for prop in label_perimeter]
    
    for l in np.arange(num_features[3]):
        perimeter[labelled_array[:,:,3]==l+1,3] = label_perimeter[l]
             
    # Refine ice_category based on the labelling.
    ice_category[(ice_category == 3) & (labelled_array[:,:,3] > 1) &
                 (perimeter[:,:,0] <= 10000) &
                 (labelled_array[:,:,2] > 1) & (labelled_array[:,:,0] > 1)] = 5 # open pack ice                 
                 
    ice_category[((ice_category == 1) | (ice_category == 2) | (ice_category == 3)) & 
                 (((perimeter[:,:,0] <= 10000) | (perimeter[:,:,3] <= 2500)) & (perimeter[:,:,1] <= 10000)) &
                 (labelled_array[:,:,0] == 1) & ((labelled_array[:,:,1] > 2) | (labelled_array[:,:,2] > 1))] = 6 # coastal polynya

    ice_category[(ice_category == 2) & (labelled_array[:,:,2] > 0) & 
                 (labelled_array[:,:,3] > 1)] = 7 # inner open water

    # Mask ice_category for plotting.
    ice_category = nemo.mask_field(ice_category, tmask)

    #plt.subplot(212)
    fig = plt.figure(figsize=(12.0,6.0))
    ax3 = plt.pcolormesh(np.squeeze(ice_category).T, cmap=basemap)
    ax3.set_clim(0,8)
    plt.axis([0, 4320, 0, 2000])
    cbar = fig.colorbar(ax3, ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    cbar.ax.set_yticklabels(['land', 'ocean', 'loose ice', 'MIZ', 'pack ice', 'open pack ice', 'coastal polynyas', 'inner open water'])
        
    # Save a hires version of the figure
    plt.savefig(''.join(['../figs/', nemodir[-29:-1], '/anim/ice_classification/ice_classification.', str(k).zfill(4), '.png']), bbox_inches='tight', dpi=600)

    # Calculate the area of each class.
    for l in np.arange(0,7+1):
        areas[k,l] = (e1t * e2t * (ice_category == l) * tmask).sum()
        
    areas[k,8] = (e1t * e2t * (ice_category > 1) * tmask).sum()
    
# Spit the numbers out to file.
np.save(''.join([homedir, nemodir, 'TIDY/ARCHIVE/POST/ice_classification_areas']), areas)
