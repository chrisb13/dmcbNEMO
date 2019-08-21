#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Tue, 13 Aug 2019 17:36:10
#   Machine created on: SB2Vbox
#

"""
Quick script to wrap PyChart
"""
import sys,os
import matplotlib
matplotlib.rcParams['backend'] = 'Agg'
from cb2logger import *
import glob
import re
import shareme as sm

def getdate(string):
    """@todo: Docstring for getdate
    :returns: @todo
    """
    exp = re.search(r'[0-9]*_[0-9]*', string) 
    exp=exp.group()
    return exp


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
    infol='/nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/g07_means/'
    infol='/nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bc337/g07_means/'
    mesh='/nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc'

    outfol=infol+'plots/'
    mkdir(outfol)

    ufiles=sorted(glob.glob(infol + '*grid-U.nc'))
    vfiles=sorted(glob.glob(infol + '*grid-V.nc'))
    assert(ufiles!=[]),"glob didn't find anything!"
    assert(vfiles!=[]),"glob didn't find anything!"

    mkdir(infol+'bsf/')
    for uf,vf in zip(ufiles,vfiles):
        #really handonline regex finder: https://regex101.com/#python
        udate=getdate(os.path.basename(uf))
        vdate=getdate(os.path.basename(vf))
        assert(udate==vdate),"file dates don't match, udate: "+ udate + ". vdate is: "+vdate

        bsff=infol+'bsf/'+udate+'_'+'cdfpsi'+'.nc'
        if not os.path.exists(bsff):
            bsff=sm.get_bsf(uf,vf,mesh,udate,infol+'bsf/')

        # go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+'bsf/'+' --mapf '+os.path.basename(bsff)+ ' --mapv sobarstf --cblvl -100000000 200000000 20000000 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p ant -o '+outfol+'sobarstf_'+os.path.basename(bsff) + '_grid-T' + ' --ft' + ' '+os.path.basename(bsff) + '_sobarstf '
        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+'bsf/'+' --mapf '+os.path.basename(bsff)+ ' --mapv sobarstf --cblvl -100000000 200000000 20000000 --cbn magma_r --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'sobarstf_'+os.path.basename(bsff) + '_grid-T' + ' --ft' + ' '+os.path.basename(bsff) + '_sobarstf '
        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

    # sys.exit('DONE all the BSF, stopping...')
    ifiles=sorted(glob.glob(infol + '*grid-T.nc'))
    assert(ifiles!=[]),"glob didn't find anything!"

    for f in ifiles:
        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+' --mapf '+os.path.basename(f)+ ' --mapv sossheig --cblvl -1.1 0.1 0.1 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'sossheig_'+os.path.basename(f) + '_grid-T' + ' --ft' + ' '+os.path.basename(f) + '_sossheig '
        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+' --mapf '+os.path.basename(f)+ ' --mapv votemper --cblvl -2 5 0.4 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'votemper_'+os.path.basename(f) + '_grid-T' + ' --ft' + ' '+os.path.basename(f) + '_votemper '
        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+' --mapf '+os.path.basename(f)+ ' --mapv soicecov --cblvl 0 1 0.1 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'soicecov_'+os.path.basename(f) + '_grid-T'+ ' --ft' + ' '+os.path.basename(f) + '_soicecov '

        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+' --mapf '+os.path.basename(f)+ ' --mapv somxl010 --cblvl 0 200 10 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'somxl010_'+os.path.basename(f) + '_grid-T'+ ' --ft' + ' '+os.path.basename(f) + '_somxl010 '
        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

        go='python /nerc/n02/n02/chbull/repos/PyChart/plot_maps.py --dir '+infol+' --mapf '+os.path.basename(f)+ ' --mapv vosaline --cblvl 33 35 0.1 --mesh /nerc/n02/n02/chbull/nemo/bld_configs/input_WED12_mathiot_eORCA025-GO7/mesh_mask_eORCA025-GO7.nc  -p south_ocean -o '+outfol+'vosaline_'+os.path.basename(f) + '_grid-T'+ ' --ft' + ' '+os.path.basename(f) + '_vosaline '
        lg.info("Running: " + go)
        subprocess.call(go,shell=True)

    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')



