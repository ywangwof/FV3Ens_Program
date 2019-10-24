#!/bin/ksh

#set -x

export TOP_DIR=${FV3DIR:-/home1/01540/ctong/sar-fv3/fv3gfs}
export WRK_DIR=${WORKDIR:-/scratch/05294/tg845932/SARFV3}
export CDATE="${eventdate}"            # format yyyymmddhh yyyymmddhh ...
export machine=odin                      # WCOSS_C,WCOSS,THEIA

export CDAS=gfs                          # gfs or gdas
export ENDHOUR=${ENDHOUR:-6}             # integration end hours
export INTHOUR=${INTHOUR:-3}             # Interval hour of external dataset

export NODES=1

#
# the following exports can all be set or just will default to what is in global_chgres_driver.sh
#
export gtype=regional            # grid type = uniform, stretch, nest, or regional
export CASE=${CASE:-C768}        # resolution of tile: 48, 96, 192, 384, 768, 1152, 3072
export LEVS=64
export LSOIL=4
export ictype=opsgfs               # opsgfs for q3fy17 gfs with new land datasets; oldgfs for q2fy16 gfs.
export nst_anl=.false.             # false or true to include NST analysis
export OMP_NUM_THREADS_CH=48       #default for openMP threads
export APRUNC=""

#

export HOMEgfs=$TOP_DIR
export fix_gsm_dir="${TOP_DIR}/fix_am"
export template="${TOP_DIR}/templates/chgres_driver.job"
export ymd=`echo $CDATE | cut -c 1-8`
export gfs_dir=${WRK_DIR}/GDAS

ulimit -s unlimited
ulimit -a

CRES=`echo $CASE | cut -c 2-`
#export OUTDIR=${WRK_DIR}/chgres_${CASE}_$ymd
export OUTDIR=${WRK_DIR}/INPUT
export FIXfv3=${WRK_DIR}/grid_orog

if [[ ! -r $OUTDIR ]]; then
  mkdir -p $OUTDIR
fi

if [ $gtype = regional ] ; then
  export REGIONAL=1
  export HALO=4
  export bchour="000"
  #
  # set the links to use the 4 halo grid and orog files
  # these are necessary for creating the boundary data
  #
  ln -sf $FIXfv3/${CASE}_grid.tile7.halo4.nc     $FIXfv3/${CASE}_grid.tile7.nc
  ln -sf $FIXfv3/${CASE}_oro_data.tile7.halo4.nc $FIXfv3/${CASE}_oro_data.tile7.nc
else
  #
  # for gtype = uniform, stretch or nest
  #
  export REGIONAL=0
fi

#
#execute the chgres driver for initial condition
#
#$TOP_DIR/ush/global_chgres_driver.sh
export TMPDIR=${WRK_DIR}/tmp

if [[ ! -d $TMPDIR ]]; then
  mkdir -p $TMPDIR
fi

cd ${WRK_DIR}

#
#execute the chgres driver for boundary conditions
#

hour=3
end_hour=${ENDHOUR}
while (test "$hour" -le "$end_hour")
do
  hour_name=$(printf %03d $hour)

    #
    #for now on theia run the BC creation sequentially
    #
    export REGIONAL=2
    export HALO=4
    export bchour=$hour_name
    #${TOP_DIR}/ush/global_chgres_driver.sh

    cat <<EOFT > sed_file
        s|WRKDIR|${WRK_DIR}|g
        s|TMPDIR|${TMPDIR}|g
        s|AAAAAA|${REGIONAL}|g
        s|BBBBBB|${CDAS}|g
        s|CCCCCC|${CRES}|g
        s|DDDDDD|${gfs_dir}|g
        s|EEEEEE|${OUTDIR}|g
        s|FFFFFF|${CDATE}|g
        s|GGGGGG|${FIXfv3}|g
        s|HHHHHH|${gtype}|g
        s|IIIIII|${fix_gsm_dir}|g
        s|JJJJJJ|${TOP_DIR}|g
        s|KKKKKK|${bchour}|g
        s|LLLLLL|${OMP_NUM_THREADS_CH}|g
        s|MMMMMM|${APRUNC}|g
        s|NNNNNN|${HALO}|g
EOFT

    runscript="chgres_C${CRES}_${CDATE}_${bchour}.job"
    cat $template | sed -f sed_file > $runscript
    rm -f sed_file

    sbatch $runscript

    hour=`expr $hour + ${INTHOUR}`
done

exit
