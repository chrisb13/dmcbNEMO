#!/opt/local/bin/python
"""
Dave Munday's script for filtering JRA fields using FFTs

#    Shared/written by Dave Munday (see email: date: Jul 19, 2019, 1:59 PM "JRA forcing set")
#    re-factored / argpassed by Christopher Bull (19/07/2019)

Usage:
    jra55_fourier_filter.py VARIABLE_NUM [ --lowcut COL --highcut COH --bandcut COLH --mkpbs PYPATH ]

Arguments:
    VARIABLE_NUM   variable number to extract (options 0-9)

Options:
    -h,--help                   : show this help message
    --lowcut COL                : low-pass cut-off (days)
    --highcut COH               : high-pass cut-off (days)
    --bandcut COLH              : bandpass cut offs low and high, where format is:  low,high
    --mkpbs PYPATH              : mkpbs script, PYPATH is where present script is located (for sending off to the scheduler)

Examples:
    python jra55_fourier_filter.py 8 --lowcut 400 --mkpbs /nerc/n02/n02/chbull/repos/nemo_wed_analyss/g07/configs/jra55_fourier_filter.py
    python jra55_fourier_filter.py 2 --bandcut 100,200 
    python jra55_fourier_filter.py 8 --highcut 5 


Notes:
    -Submit when logged into the DAC
    -bandcut cuts off everything *inbetween* low,high
    -cut offs should be natural numbers, no floats
    -if mkpbs is specified, VARIABLE_NUM is ignored

Here's a way to check that the runs worked (or got to the last tile):
    cd /nerc/n02/n02/chbull/RawData/JRA55-do_v1-3_400low
    for x in `ls *mkpbs*.sh.o*`; do echo $x; grep 'm = ' $x|tail -1; echo ''; done
"""
# --------------------------------------------------------------------------- #

# Import required modules.
from docopt import docopt
arguments = docopt(__doc__)
# print(arguments)
import glob
import netCDF4 as nc
import numpy as np
import os
import shutil
import sys
import contextlib as ctx
import time


pbsheader=\
"""
#!/bin/bash --login
#PBS -l select=serial=true:ncpus=1
#PBS -l walltime=96:00:00
#PBS -A n02-FISSA

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

export PYTHONPATH=~/anaconda2/pkgs;export PATH=~/anaconda2/bin:$PATH;source activate root
# export PYTHONPATH=/work/n02/n02/chbull/anaconda2/pkgs;export PATH=/work/n02/n02/chbull/anaconda2/bin:$PATH;source activate root
# export PYTHONPATH=~/anaconda2/pkgs;export PATH=~/anaconda2/bin:/opt/torque/4.2.6.1/gcc/4.4.7/bin:/opt/maui/3.3.1/gcc/4.4.7/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/ibutils/bin:/nerc/n02/n02/chbull/xclip/:/nerc/n02/n02/chbull/bin;source activate root

"""


def mkdir(p):
    """make directory of path that is passed"""
    try:
       os.makedirs(p)
       print("output folder: "+p+ " does not exist, we will make one.")
    except OSError as exc: # Python >2.5
       import errno
       if exc.errno == errno.EEXIST and os.path.isdir(p):
          pass
       else: raise

# --------------------------------------------------------------------------- #

# Specify where the original JRA data lives.
jradir = '/nerc/n01/n01/munday/ORCHESTRA/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA55-do/v1-3/'
# jradir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA-do/v1-3/'

jradir='/nerc/n02/n02/chbull/RawData/MonsoonForcings/JRA55-do/v1-3/' # CB
jradir='/nerc/n02/n02/chbull/RawData/MonsoonForcings/JRA55-do/v1-3/dave_filter_links/' # CB

# Specify where the filtered JRA output data should go.
outdir = '/nerc/n01/n01/munday/ORCHESTRA/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA55-do/v1-3-fil/'
# outdir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA-do/v1-3-fil/'


# --------------------------------------------------------------------------- #

#Error traps/mungin"
if int(arguments['VARIABLE_NUM']) not in np.arange(10):
    sys.exit('ERROR: Please provide a variable_num b/w 0-9')

if arguments['--lowcut'] is not None:
    cutter=float(arguments['--lowcut'])
    print('doing a --lowcut with '+str(cutter))
    if arguments['--bandcut'] is not None or arguments['--highcut'] is not None:
        sys.exit('ERROR: not allowed')

elif  arguments['--bandcut'] is not None:
    cutter= [float(s) for s in arguments['--bandcut'].split(',')]
    print('doing a --bandcut with '+str(cutter))

    if arguments['--lowcut'] is not None or arguments['--highcut'] is not None:
        sys.exit('ERROR: not allowed')

elif  arguments['--highcut'] is not None:
    cutter=float(arguments['--highcut'])
    print('doing a --highcut with '+str(cutter))

    if arguments['--lowcut'] is not None or arguments['--bandcut'] is not None:
        sys.exit('ERROR: not allowed')
else:
    print('error')
    __import__('pdb').set_trace()

# outdir='/nerc/n02/n02/chbull/RawData/JRA55-do_v1-3_400low/'  #CB

if arguments['--lowcut'] is not None:
    outdir='/nerc/n02/n02/chbull/RawData/JRA55-do_v1-3_'+str(int(cutter))+'_lowcut/'
elif  arguments['--bandcut'] is not None:
    outdir='/nerc/n02/n02/chbull/RawData/JRA55-do_v1-3_'+str(int(cutter[0]))+'-'+str(int(cutter[1]))+'_bandcut/'
elif  arguments['--highcut'] is not None:
    outdir='/nerc/n02/n02/chbull/RawData/JRA55-do_v1-3_'+str(int(cutter))+'_highcut/'
else:
    print('error')
    __import__('pdb').set_trace()

mkdir(outdir)

if arguments['--mkpbs'] is not None:
    if arguments['--lowcut'] is not None:
        run='--lowcut '+arguments['--lowcut']
    elif  arguments['--bandcut'] is not None:
        run='--bandcut '+arguments['--bandcut']
    elif  arguments['--highcut'] is not None:
        run='--highcut '+arguments['--highcut']

    print("Copy and paste the following to get the party started...")

    print("cd "+outdir)
    for file_num in [str(r) for r in np.arange(10)]:
        with ctx.closing(open(outdir+'mkpbs'+file_num.zfill(2)+'.sh','w')) as handle:
             handle.write(pbsheader+"\n")
             handle.write('python '+arguments['--mkpbs'] +' '+str(file_num)+' '+run +"\n")

        print("qsub "+outdir+'mkpbs'+file_num.zfill(2)+'.sh')
    sys.exit("Created all the PBS files, now exiting.")

# --------------------------------------------------------------------------- #

# Create a list of filename prefixes that we're going to use to search
# for matching files.
# prefix = ['cb_rainsnowfix/totalprecip', 'q_10', 'rlds', 'rsds', 'runoff_all', 'slp', 'snow', 't_10', 'u_10', 'v_10']
# prefix = ['cb_rainsnowfix/totalprecip']
prefix = [['cb_rainsnowfix/totalprecip', 'q_10', 'rlds', 'rsds', 'runoff_all', 'slp', 'snow', 't_10', 'u_10', 'v_10'][int(arguments['VARIABLE_NUM'])] ]

# Create a list of variable names that we want to load from the files.
# varnam = ['prrn', 'huss_10m', 'rlds', 'rsds', 'friver', 'psl', 'prsn', 'tas_10m', 'uas_10m', 'vas_10m']
varnam = ['prrn']
varnam = [['prrn', 'huss_10m', 'rlds', 'rsds', 'friver', 'psl', 'prsn', 'tas_10m', 'uas_10m', 'vas_10m'][int(arguments['VARIABLE_NUM'])]]

# Set the date stamp to put in the filename.
datestamp = '19Jul2019'

print('doing a filter with  '+str(prefix) + ' and var '+str(varnam)+ ' with datestamp '+ datestamp)

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

    assert(files!=[]),"glob didn't find anything when trying: "+''.join([jradir, i, '*.nc'])

    # Set the indices for loops over X & Y tiles.
    minter = 0
    ninter = 0

    # Loop over the X & Y indice in each file.
    for m in np.arange(npy):
        for n in np.arange(npx):

            start=time.time()

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

            cutter_sec=cutter*60*60*24
            if  arguments['--highcut'] is not None or arguments['--lowcut'] is not None:
                # cutter_sec=172800
                # Work out which index of the fft corresponds to frequency of 1/(2 days).
                index_to_eliminate = next(index for index, value in enumerate(np.fft.rfftfreq(variable.shape[0], d=3.0*3600)) if value > 1.0/cutter_sec) - 1
            elif  arguments['--bandcut'] is not None:
                index_to_eliminate1 = next(index for index, value in enumerate(np.fft.rfftfreq(variable.shape[0], d=3.0*3600)) if value > 1.0/cutter_sec[0]) - 1
                index_to_eliminate2 = next(index for index, value in enumerate(np.fft.rfftfreq(variable.shape[0], d=3.0*3600)) if value > 1.0/cutter_sec[1]) - 1

            if  arguments['--highcut'] is not None:
                # Eliminate all frequencies equal or higher than 1/(2 days).
                variable_fft[index_to_eliminate:, :, :] = 0.0
            elif arguments['--lowcut'] is not None:
                #retain mean (1:)
                variable_fft[1:index_to_eliminate, :, :] = 0.0
            elif arguments['--bandcut'] is not None:
                #want to keep: variable_fft[index_to_eliminate1:index_to_eliminate2, :, :], so..
                # variable_fft[index_to_eliminate2:, :, :] = 0.0
                # variable_fft[:index_to_eliminate1, :, :] = 0.0
                variable_fft[index_to_eliminate1:index_to_eliminate2, :, :] = 0.0

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
                    mkdir(os.path.dirname(outfilename)+'/')
                    shutil.copy2(j, outfilename)

                # Open the duplicated netcdf file.
                outfile = nc.Dataset(outfilename, 'r+', format="NETCDF3_CLASSIC")

                # Replace the variable with the right bit of the filtered data.
                outfile[varnam[iindex]][:, minter*nsx:(minter+1)*nsx, ninter*nsy:(ninter+1)*nsy] = variable[int(entries[0:jindex].sum()):int(entries[0:jindex+1].sum()), :, :]

                # Close the nc file, which writes it out to disk.
                outfile.close()

                # Increment the files index.
                jindex = jindex + 1

            # Increase the index for the X tiles.
            ninter = ninter + 1

            # Initialise the fft array back to nothing.
            variable_fft = None

            finish=time.time()
            print("time for tile: "+str(finish-start))

        # Increase the index for the Y tiles.
        minter = minter + 1

        # Initialise the index for the X tiles.
        ninter = 0

    # Increment the prefix index by 1.
    iindex = iindex + 1

    print("Filtering finished.")
# --------------------------------------------------------------------------- #
