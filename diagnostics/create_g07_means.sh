###### alternatively for multiple procs
##!/bin/bash --login
##PBS -l select=serial=true:ncpus=8
##PBS -l walltime=03:30:00
##PBS -A n02-FISSA

## Make sure any symbolic links are resolved to absolute path
#export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

## Change to the directory that the job was submitted from
#cd $PBS_O_WORKDIR

#export PYTHONPATH=/nerc/n02/n02/chbull/anaconda3/pkgs;export PATH=/nerc/n02/n02/chbull/anaconda3/bin:$PATH;source activate root

#module unload cray-netcdf-hdf5parallel/4.4.1.1
#module unload cray-hdf5-parallel/1.10.0.1
#module load nco

#python mk_g07_means.py grid-V 2018 2018 &
#python mk_g07_means.py grid-V 2019 2019 &
#python mk_g07_means.py grid-V 2020 2020 &
#python mk_g07_means.py grid-V 2021 2021 &
#python mk_g07_means.py grid-V 2022 2022 &
#python mk_g07_means.py grid-V 2023 2023 &
#python mk_g07_means.py grid-V 2024 2024 &
#python mk_g07_means.py grid-V 2025 2025 &
#python mk_g07_means.py grid-V 2026 2026 &
#cd /nerc/n02/n02/chbull/repos/nemo_wed_analysis/g07/diagnostics
#python mk_g07_means.py grid-V 2027 2027 &
#python mk_g07_means.py grid-V 2028 2028 &
#python mk_g07_means.py grid-V 2029 2029 &
#python mk_g07_means.py grid-V 2030 2030 &

#python mk_g07_means.py grid-U 2027 2027 &
#python mk_g07_means.py grid-U 2028 2028 &
#python mk_g07_means.py grid-U 2029 2029 &
#python mk_g07_means.py grid-U 2030 2030 &
wait

YEARi=2027
YEARf=2030


source /nerc/n02/n02/chbull/.vim/common_bashfunctions
for YEAR in $(seq $YEARi $YEARf)
do

    echo "Creating mean of YEAR : ${YEAR}"
    cd /nerc/n02/n02/chbull/repos/nemo_wed_analysis/g07/diagnostics;qsubdac "python mk_g07_means.py grid-T ${YEAR} ${YEAR} --indir /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/ --output /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/g07_means/" 
    cd /nerc/n02/n02/chbull/repos/nemo_wed_analysis/g07/diagnostics;qsubdac "python mk_g07_means.py grid-U ${YEAR} ${YEAR} --indir /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/ --output /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/g07_means/" 
    cd /nerc/n02/n02/chbull/repos/nemo_wed_analysis/g07/diagnostics;qsubdac "python mk_g07_means.py grid-V ${YEAR} ${YEAR} --indir /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/ --output /nerc/n02/shared/chbull/NEMO_JRAspinup4DAVECHK/u-bl504/onm.nc.file/g07_means/" 
    echo ""
done


