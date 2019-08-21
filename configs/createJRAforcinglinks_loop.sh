#   Author: Christopher Bull. 
#   Affiliation:  British Antarctic Survey
#                 Cambridge, UK
#   Contact: chbull@bas.ac.uk
#   Date created: Tue, 06 Aug 2019 16:05:12
#   Machine created on: SB2Vbox
#

###########################################################
# This is a quick script to create some sym links for JRA #
###########################################################

#UPDATED: hacked the original to now handle looping on the JRA forcing after 2017 is finished.
#NB: YEAR is model year to be used by ROSE/NEMO
#  : NYEAR is the JRA year (so loops back anytime after/including 2018)
#  : YEAR0 needs to be 1958 for this to work

YEAR0=1958
YEAR_MAX=2317
FORCING=/projects/jmmp/chrbu/RawData/MonsoonForcings/JRA55-do/v1-3

#used to be different on ARCHER
WORKDIR=/projects/jmmp/chrbu/RawData/MonsoonForcings/JRA55-do/v1-3

echo "Linking surfacefields for forcing: "
CNT=1
JRACYC=1
for YEAR in $(seq ${YEAR0} ${YEAR_MAX})
    do

    if   [ ${CNT} -eq 61  ]; then
        echo ''
        echo ''
        JRACYC=2
    elif [ ${CNT} -eq 121  ]; then
        echo ''
        echo ''
        JRACYC=3
    elif [ ${CNT} -eq 181  ]; then
        echo ''
        echo ''
        JRACYC=4
    elif [ ${CNT} -eq 241  ]; then
        echo ''
        echo ''
        JRACYC=5
    elif [ ${CNT} -eq 301  ]; then
        echo ''
        echo ''
        JRACYC=6
    elif [ ${CNT} -eq 361  ]; then
        echo ''
        echo ''
        JRACYC=7
    fi 

    if   [ ${JRACYC} -eq 1  ]; then
        NYEAR=${YEAR}
    elif   [ ${JRACYC} -eq 2  ]; then
        NYEAR=`expr ${YEAR} - 60`
    elif   [ ${JRACYC} -eq 3  ]; then
        NYEAR=`expr ${YEAR} - 120`
    elif   [ ${JRACYC} -eq 4  ]; then
        NYEAR=`expr ${YEAR} - 180`
    elif   [ ${JRACYC} -eq 5  ]; then
        NYEAR=`expr ${YEAR} - 240`
    elif   [ ${JRACYC} -eq 6  ]; then
        NYEAR=`expr ${YEAR} - 300`
    elif   [ ${JRACYC} -eq 7  ]; then
        NYEAR=`expr ${YEAR} - 360`
    fi 

    if   [ ${NYEAR} = "2017" ]; then
	JRA_DATE=14Jan2018
    #elif [ ${YEAR} = "2018" ]; then
	#JRA_DATE=07Feb2018
    else
	JRA_DATE=18Oct2017
    fi 

    #just as a check..
    #echo $CNT,${YEAR},${JRA_DATE},'mooo',${NYEAR},${JRACYC}

    ln -s ${FORCING}/u_10.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/u10_JRA_y${YEAR}.nc
    ln -s ${FORCING}/v_10.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/v10_JRA_y${YEAR}.nc   
    ln -s ${FORCING}/rsds.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/rsds_JRA_y${YEAR}.nc     
    ln -s ${FORCING}/rlds.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/rlds_JRA_y${YEAR}.nc     
    ln -s ${FORCING}/t_10.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/t10_JRA_y${YEAR}.nc      
    ln -s ${FORCING}/q_10.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/q10_JRA_y${YEAR}.nc      

    ##totalprecip.${YEAR}.${JRA_DATE}.nc includes both rain+snow (done with CDO offline)
    ln -s ${FORCING}/cb_rainsnowfix/totalprecip.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/rain_JRA_y${YEAR}.nc     
    ln -s ${FORCING}/snow.${NYEAR}.${JRA_DATE}.nc ${WORKDIR}/inputs_loop/snow_JRA_y${YEAR}.nc     

    CNT=`expr ${CNT} + 1`
done
echo "Done"
