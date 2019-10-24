#/bin/bash

HH=$1     # forecast starting hour
NN=$2     # Number of ensember members

DATE=20190520

if [[ "${HH}" < "18" ]]; then
  DATE=$(date -d "$DATE 1 day" +%Y%m%d)
fi

GDADDIR="/oldscratch/nyussouf/GDAS"

for m in $(seq 1 1 $2); do
  mem=$(printf "%03d" $m)

  if [[ ! -d ${DATE}${HH}_mem${mem} ]]; then
    mkdir ${DATE}${HH}_mem${mem}
  fi

  cd ${DATE}${HH}_mem${mem}

  #
  # Analysis
  #
  tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.anl.tar ./gdas.t${HH}z.ratmanl.mem${mem}.nemsio ./gdas.t${HH}z.sfcanl.mem${mem}.nemsio
  ln -s gdas.t${HH}z.ratmanl.mem${mem}.nemsio gfs.t${HH}z.atmanl.nemsio
  ln -s gdas.t${HH}z.sfcanl.mem${mem}.nemsio  gfs.t${HH}z.sfcanl.nemsio

  #
  # forecast hour 03, 09
  #
  for FH in "03" "09"; do
    tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.fcs${FH}.tar ./gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio ./gdas.t${HH}z.sfcf0${FH}.mem${mem}.nemsio
    ln -s gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio gfs.t${HH}z.atmf0${FH}.nemsio
  done

  #
  # forecast hour 06
  #
  FH="06"
  tar xvf ${GDADDIR}/gpfs_hps_nco_ops_com_gfs_prod_enkf.${DATE}_${HH}.fcs.tar ./gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio ./gdas.t${HH}z.sfcf0${FH}.mem${mem}.nemsio
  ln -s gdas.t${HH}z.atmf0${FH}s.mem${mem}.nemsio gfs.t${HH}z.atmf0${FH}.nemsio

  cd ../
done

exit 0
