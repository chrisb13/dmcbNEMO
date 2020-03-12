#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# My handrolled modules.
import nemo

# Import required modules.
import copy
import glob
import netCDF4 as nc
import numpy as np
import os
import shutil

# --------------------------------------------------------------------------- #

# Specify where the input data lives.
homedir = '/nerc/n01/n01/munday/ORCHESTRA/'
# homedir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/'
jra55dir = 'trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA55-do/v1-3/'

# Specify where to output the modified forcing.
outdir = 'trunk/NEMOGCM/CONFIG/JRA55TAU-ORCH0083-LIM3/JRA/JRA55-do/v1-3-tau/'

# --------------------------------------------------------------------------- #

# Specify the amount to scale the wind fields by.
scale_factor = 1.1

# --------------------------------------------------------------------------- #

# Find the files that I want to rescale the U & V winds in.
ufiles = sorted(glob.glob(''.join([homedir, jra55dir, 'u_10.*'])))
vfiles = sorted(glob.glob(''.join([homedir, jra55dir, 'v_10.*'])))

# --------------------------------------------------------------------------- #

# Specify the number of grid boxes in JRA55.
nx = 640
ny = 320

# Load the lats & lons from the first ufile.
lat = np.repeat(nemo.load_field('latitude', '', '', ufiles[0])[:,None], 640, axis=1).T
lon = np.repeat(nemo.load_field('longitude', '', '', ufiles[0])[:,None], 320, axis=1)

# --------------------------------------------------------------------------- #

# Construct the mask for the wind scaling based upon latitude.
mask = copy.deepcopy(lat)
mask[mask<-70] = 0.
mask[mask>-30] = 0.
mask[mask<0] = 1.

for j in np.arange(320):
    if lat[0,j] > -70. and lat[0,j] < -60.:
        mask[:,j] = 0.5*(1.0 - np.cos(np.pi*(lat[:,j]+70.)/10.))
    elif lat[0,j] > -40. and lat[0,j] < -30.:
        mask[:,j] = 0.5*(1.0 + np.cos(np.pi*(lat[:,j]+40.)/10.))

# --------------------------------------------------------------------------- #

# Loop over the U/V files and calculate the turning angle..
for k in range(len(ufiles)):
    
    # Print the name of the current input file to screen.
    print('infile name =', ufiles[k])
 
    # Load the U wind field.
    uwind = np.squeeze(nemo.load_field('uas_10m', '', '', ufiles[k]))

    # Scale the U wind.
    uwind = uwind + (scale_factor-1.0) * mask[:, :, None] * uwind

    # Construct the name of the output file.
    outfilename = ''.join([homedir, outdir, ufiles[k].split('/')[-1]])

    # Print the current output file name.
    print('outfile name =', outfilename)

    # If the output file doesn't exist, copy the original to it.
    if not os.path.isfile(outfilename):
        shutil.copy2(ufiles[k], outfilename)

    # Open the duplicated netcdf file.
    outfile = nc.Dataset(outfilename, 'r+', format="NETCDF3_CLASSIC")

    # Replace the variable with the right bit of the fitered data.
    outfile['uas_10m'][:] = uwind.T

    # Close the nc file, which writes it out to disk.
    outfile.close()
    
    # Print the name of the current input file to screen.
    print('infile name =', vfiles[k])

    # Load the V wind field.
    vwind = np.squeeze(nemo.load_field('vas_10m', '', '', vfiles[k]))

    # Scale the U wind.
    vwind = vwind + (scale_factor-1.0) * mask[:, :, None] * vwind

    # Construct the name of the output file.
    outfilename = ''.join([homedir, outdir, vfiles[k].split('/')[-1]])

    # Print the current output file name.
    print('outfile name =', outfilename)

    # If the output file doesn't exist, copy the original to it.
    if not os.path.isfile(outfilename):
        shutil.copy2(vfiles[k], outfilename)

    # Open the duplicated netcdf file.
    outfile = nc.Dataset(outfilename, 'r+', format="NETCDF3_CLASSIC")

    # Replace the variable with the right bit of the fitered data.
    outfile['vas_10m'][:] = vwind.T

    # Close the nc file, which writes it out to disk.
    outfile.close()

# --------------------------------------------------------------------------- #
