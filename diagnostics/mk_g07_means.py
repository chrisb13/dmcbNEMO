#!/opt/local/bin/python
"""
Script to creat time-mean files for G07 runs

Usage:
    mk_g07_means.py FTYPE START END [ --indir NCFILE --output PYPATH ]

Arguments:
    FTYPE   file type to extract e.g. (grid-T, grid-U, etc)
    START   year to start the extraction
    END     last year to end the extraction

Options:
    -h,--help                   : show this help message
    --indir NCFILE               : folder to get G07 files from
    --output PYPATH              : output folder to put files

Examples:
    python mk_g07_means.py grid-T 2023 2024
    python mk_g07_means.py grid-T 2023 2023
    python mk_g07_means.py grid-T 2023 2024 --indir /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/ --ouput /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/g07_means/

Notes:
    -Requires nco to be loaded
    -Dates can be just years, or include month/day
    -Assumes one file per month

    - (Xarray method was abandoned b/c) has a problem with some versions of xarray where the slicing doesn't work, e.g.
        TypeError: cannot compare cftime._cftime.DatetimeNoLeap(2022, 9, 16, 0, 0, 0, 0, 0, 259) and '2026'

#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Tue, 13 Aug 2019 10:41:42
#   Machine created on: SB2Vbox
"""
# --------------------------------------------------------------------------- #

# Import required modules.
from docopt import docopt
arguments = docopt(__doc__)
# print(arguments)
import glob
import numpy as np
import os
import shutil
import sys
import contextlib as ctx
import time
import glob
import pandas as pd

# import xarray as xr

from netCDF4 import Dataset
from netCDF4 import num2date
import collections

import sys,os
# sys.path.insert(1,os.path.expanduser('~/hdrive/repos/cms_analysis/'))
from cb2logger import *

if __name__ == "__main__": 
    LogStart('',fout=False)
    #put useful code here!

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')


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

def get_timemean(infol,tstart,tend,ftype,outputpath):
    """@todo: Docstring for get_timemean
    
    :filelist: @todo
    :returns: @todo
    """

    fout=tstart+'_'+tend+'_'+ftype+'.nc'
    if not os.path.exists(outputpath+fout):
        lg.info("Time-mean has not been created, making one with NCO now...")

        lg.info("Creating index of time stamps available")
        ifiles=sorted(glob.glob(infol+'*'+ftype+'.nc' ))
        assert(ifiles!=[]),"glob didn't find anything!"
        lg.info("Found " + str(len(ifiles)) + ' files')

        dates=collections.OrderedDict()
        for f in ifiles:
            # print f
            ifile=Dataset(f, 'r')

            # this gave screwy files..
            # dates[f]=num2date(ifile['time_counter'][:],units=ifile['time_counter'].units).tolist()[0]

            dates[f]=num2date(ifile['time_instant'][:],units=ifile['time_counter'].units).tolist()[0]
            ifile.close()

        # Abanadoned b/c it was unreliable on some environments
        # fls=xr.open_mfdataset(ifiles,concat_dim='time_counter')
        # df=pd.DataFrame({'fls':ifiles},index=fls['time_counter'].values)

        df=pd.DataFrame({'fls':ifiles},index=dates.values())
        lg.info("Dates available, from " + str(df.index[0])+', to '+ str(df.index[-1]) )
        cut=df[tstart:tend]

        # for debugging
        # realfiles=[os.path.basename(f) for f in df[tstart:tend]['fls']]

        # lg.info(outputpath+fout)
    
        lg.info("")
        lg.info("Now executing: ")
        lg.info('ncra ' + ' '.join(cut['fls'].tolist()) + ' ' +outputpath+fout)
        lg.info("")

        assert(len(cut)!=0),"Error with your time slice!"

        # print fout
        subprocess.call('ncra ' + ' '.join(cut['fls'].tolist()) + ' ' +outputpath+fout,shell=True)
        if not os.path.exists(outputpath+fout):
            lg.error("File NOT created in: " + outputpath+fout)
            sys.exit()
        else:
            lg.info("File created in: " + outputpath+fout)

        # fls.close()
    else:
        lg.info("File already existed in: " + outputpath+fout)
    return outputpath+fout


# --------------------------------------------------------------------------- #

# Specify where the original JRA data lives.
indir= '/nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/'

if arguments['--indir'] is not None:
    indir=arguments['--indir']
    lg.info("Indir netcdf files specified: " + indir)

# Specify where the filtered JRA output data should go.
outdir = indir + 'g07_means/'
# outdir = '/Users/munday/Documents/Projects/ORCHESTRA/NEMO/trunk/NEMOGCM/CONFIG/JRA55IAF-ORCH0083-LIM3/JRA/JRA-do/v1-3-fil/'

if arguments['--output'] is not None:
    outdir=arguments['--output']
mkdir(outdir)
lg.info("Files being created in: " + outdir)

get_timemean(indir,arguments['START'],arguments['END'],arguments['FTYPE'],outdir)
