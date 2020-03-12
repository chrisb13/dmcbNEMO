#!/opt/local/bin/python

# --------------------------------------------------------------------------- #

# Import required modules.
import glob
import netCDF4 as nc
import numpy as np
import os
import shutil

# --------------------------------------------------------------------------- #

# Specify where the original JRA data lives.
jradir = '/nerc/n01/n01/munday/ORCHESTRA/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA55-do/v1-3/'
# jradir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA-do/v1-3/'

# Specify where the filtered JRA output data should go.
outdir = '/nerc/n01/n01/munday/ORCHESTRA/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA55-do/v1-3-fil/'
# outdir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA-do/v1-3-fil/'

# --------------------------------------------------------------------------- #

# Create a list of filename prefixes that we're going to use to search
# for matching files.
# prefix = ['cb_rainsnowfix/totalprecip', 'q_10', 'rlds', 'rsds', 'runoff_all', 'slp', 'snow', 't_10', 'u_10', 'v_10']
prefix = ['rsds']

# Create a list of variable names that we want to load from the files.
# varnam = ['prrn', 'huss_10m', 'rlds', 'rsds', 'friver', 'psl', 'prsn', 'tas_10m', 'uas_10m', 'vas_10m']
varnam = ['rsds']

# Set the date stamp to put in the filename.
datestamp = '12Jun2019'

# Set the number of divisions in X & Y to avoid memory overruns.
npx = 8
npy = 4

# Calculate the number of grid boxes in each tile.
nsx = 640 / npx
nsy = 320 / npy

# --------------------------------------------------------------------------- #

# Set the prefix index.
iindex = 0

# Loop over the entries in prefix.
for i in prefix:

    # Find the files that match the current prefix.
    files = sorted(glob.glob(''.join([jradir, i, '*.nc'])))

    # Set the indices for loops over X & Y tiles.
    minter = 0
    ninter = 0

    # Loop over the X & Y indice in each file.
    for m in np.arange(npy):
        for n in np.arange(npx):

            # Initialise the array that's going to hold the data.
            variable = 0.0

            # Initialise the array that's going to hold the size of time.
            entries = np.ndarray((len(files), 1))

            # Initialise the files index.
            jindex = 0

            # Loop over the files that were found.
            for j in files:

                # Tell me where I am:
                print 'm =', m, 'n =', n

                # Print the current file name.
                print 'inputfilename = ', j.split('/')[-1]

                # Open the dataset.
                ncfile = nc.Dataset(j)

                # Load the required variable from the current file.
                if jindex == 0:
                    variable = ncfile.variables[varnam[iindex]][:, minter*nsx:(minter+1)*nsx, ninter*nsy:(ninter+1)*nsy]
                else:
                    variable = np.concatenate((variable, ncfile.variables[varnam[iindex]][:, minter*nsx:(minter+1)*nsx, ninter*nsy:(ninter+1)*nsy]), axis=0)

                # Record how big the record dimension is for each file.
                entries[jindex] = ncfile.variables['time'].size

                # Close the file.
                ncfile.close()

                # Increment the files index.
                jindex = jindex + 1

                # End loop over files.

            # Find the minimum value of the unfiltered variable & print to screen.
            varmin = variable.min()
            print 'og minimum value =', varmin

            # Take the FFT.
            print 'Taking fft...'
            variable_fft = np.fft.rfft(variable, axis=0)

            # Work out which index of the fft corresponds to frequency of 1/(2 days).
            index_to_eliminate = next(index for index, value in enumerate(np.fft.rfftfreq(variable.shape[0], d=3.0*3600)) if value > 1.0/172800) - 1

            # Eliminate all frequencies equal or higher than 1/(2 days).
            variable_fft[index_to_eliminate:, :, :] = 0.0

            # Invert the FFT to get the time series back.
            print 'Inverting fft...'
            variable = np.fft.irfft(variable_fft, n=variable.shape[0], axis=0)

            # Find the minimum value of the filtered variable & print to screen.
            filmin = variable.min()
            print 'filtered minimum value =', filmin

            # If the original minimum value is >= 0 & the filtered minimum
            #  < 0, set all filtered values <= 0 to the original minimum.
            if varmin >= 0 and filmin < 0:
                variable[variable < varmin] = varmin

            # Print the new minimum value of the filtered variable to screen
                print 'output minimum value =', variable.min()

            # Reset the files index.
            jindex = 0

            # Second loop over files for output.
            for j in files:

                # Construct the name of the output file.
                outfilename = ''.join([outdir, i.split('.')[0], '.', j.split('/')[-1].split('.')[1], '.', datestamp, '.nc'])

                # Print the current output file name.
                print 'outfile name =', outfilename.split('/')[-1]

                # If the output file doesn't exist, copy the original to it.
                if not os.path.isfile(outfilename):
                    shutil.copy2(j, outfilename)

                # Open the duplicated netcdf file.
                outfile = nc.Dataset(outfilename, 'r+', format="NETCDF3_CLASSIC")

                # Replace the variable with the right bit of the fitered data.
                outfile[varnam[iindex]][:, minter*nsx:(minter+1)*nsx, ninter*nsy:(ninter+1)*nsy] = variable[int(entries[0:jindex].sum()):int(entries[0:jindex+1].sum()), :, :]

                # Close the nc file, which writes it out to disk.
                outfile.close()

                # Increment the files index.
                jindex = jindex + 1

            # Increase the index for the X tiles.
            ninter = ninter + 1

            # Initialise the fft array back to nothing.
            variable_fft = None

        # Increase the index for the Y tiles.
        minter = minter + 1

        # Initialise the index for the X tiles.
        ninter = 0


    # Increment the prefix index by 1.
    iindex = iindex + 1

# --------------------------------------------------------------------------- #
