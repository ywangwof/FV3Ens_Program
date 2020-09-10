#/bin/bash

HH=$1     # forecast starting hour
NN=${2-40}     # Number of ensember members

if [[ ! $1 =~ ^[0-9]{10}$ ]]; then
  echo "$0 YYYYMMDDHH [NN]"
  exit
fi

DATE=${1:0:8}
HH=${1:8:2}

#WRKDIR=/scratch/ywang/EPIC/GDAS
WRKDIR=/work/00315/tg455890/stampede2/EPIC/GDAS

cd $WRKDIR

if [[ "${HH}" < "18" ]]; then
  RDATE=$(date -d "$DATE 1 day ago" +%Y%m%d)
else
  RDATE=$DATE
fi

#GDADDIR="/oldscratch/nyussouf/GDAS/20190520/$RDATE"
GDADDIR="/scratch/06809/jpark217/GFS"

for m in $(seq 1 1 $NN); do
  mem=$(printf "%03d" $m)

  if [[ ! -d ${DATE}${HH}_mem${mem} ]]; then
    mkdir ${DATE}${HH}_mem${mem}
  fi

  srcdir="$GDADDIR/enkf.$1"

  cd ${DATE}${HH}_mem${mem}

  #
  # Analysis
  #
  #tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.anl.tar ./gdas.t${HH}z.ratmanl.mem${mem}.nemsio ./gdas.t${HH}z.sfcanl.mem${mem}.nemsio
  #ln -s gdas.t${HH}z.ratmanl.mem${mem}.nemsio gfs.t${HH}z.atmanl.nemsio
  #ln -s gdas.t${HH}z.sfcanl.mem${mem}.nemsio  gfs.t${HH}z.sfcanl.nemsio
  ln -s $srcdir/gdas.t${HH}z.ratmanl.mem${mem}.nemsio .
  ln -s $srcdir/gdas.t${HH}z.sfcanl.mem${mem}.nemsio  .

  #
  # forecast hour 03, 09
  #
  for FH in "03" "09"; do
    #tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.fcs${FH}.tar ./gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio ./gdas.t${HH}z.sfcf0${FH}.mem${mem}.nemsio
    #ln -s gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio gfs.t${HH}z.atmf0${FH}.nemsio
    ln -s $srcdir/gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio .
  done

  #
  # forecast hour 06
  #
  FH="06"
  #tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.fcs.tar ./gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio ./gdas.t${HH}z.sfcf0${FH}.mem${mem}.nemsio
  #ln -s gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio gfs.t${HH}z.atmf0${FH}.nemsio
  ln -s $srcdir/gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio .

  cd ../
done

exit 0
