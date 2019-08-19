#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Mon, 19 Aug 2019 16:48:27
#   Machine created on: SB2Vbox
#
"""
Quick script to make links for experiments u-bc337 and bl504 to create a continuous time-series of files.
"""
import sys,os
from cb2logger import *
import glob
def mkdir(p):
    """make directory of path that is passed"""
    try:
       os.makedirs(p)
       lg.info("output folder: "+p+ " does not exist, we will make one.")
    except OSError as exc: # Python >2.5
       import errno
       if exc.errno == errno.EEXIST and os.path.isdir(p):
          pass
       else: raise

if __name__ == "__main__": 
    LogStart('',fout=False)
    bc337='/nerc/n02/n02/chbull/RawData/u-bc337/onm.nc.file/'
    bl504='/nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/'
    odir='/nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bc337-bl504_links/'

    mkdir(odir)

    bcif=sorted(glob.glob(bc337+'*.nc'))
    assert(bcif!=[]),"glob didn't find anything!"

    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20180101-20180201_grid-T.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20180101-20180201_grid-U.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20180101-20180201_grid-V.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20180101-20180201_grid-W.nc')

    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20171201-20180101_grid-T.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20171201-20180101_grid-U.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20171201-20180101_grid-V.nc')
    bcif.remove(os.path.dirname(bc337)+'/'+'nemo_bc337o_1m_20171201-20180101_grid-W.nc')

    blif=sorted(glob.glob(bl504+'nemo_bl504o_1m_20??????-20??????_grid-?.nc'))
    assert(blif!=[]),"glob didn't find anything!"

    for f in (bcif+blif):
        lg.info("Creating link: " +' nemospinny_'+os.path.basename(f)[15:]+ ' pointing towards: ' + f)
        os.symlink(f,odir+'nemospinny_'+os.path.basename(f)[15:])

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
