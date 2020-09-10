#!/bin/bash

if [[ ! $1 =~ ^[0-9]{10}$ ]]; then
  echo "$0 YYYYMMDDHH [test_runs|test_mp|test_spp|test_mspp] [NS] [NN]"
  exit
fi

casedate=${1-2019052018}
wrkcase=${2-"test_runs"}
numsta=${3-1}
numens=${4-40}
username="yunheng.wang"
maxjobs=40

wrkroot="/scratch/ywang/EPIC"

caseHH=${casedate:8:2}
caseDT=$(date -d "${casedate:0:8} ${caseHH}:00 1 hours ago" +%Y%m%d)

rootdir="/oldscratch/ywang/EPIC/program.git"

wrkdir="${wrkroot}/${wrkcase}/${casedate}"
gdasdir="/scratch/ywang/EPIC/GDAS"

#
# Link grid and orog files
#
if [[ ! -e $wrkdir/grid_orog ]]; then
  griddir="/scratch/ywang/comFV3SAR/test_runs/caps_cntl_${caseDT}"

  mkdir -p $wrkdir/grid_orog

  cd $wrkdir/grid_orog

  fn=$(ls $griddir/grid/C*_grid.tile7.halo3.nc)
  fm=${fn##$griddir/grid/}
  CCASE=${fm%%_grid.tile7.halo3.nc}

  ln -sf $griddir/grid/${CCASE}_grid.tile7.halo*.nc     .
  ln -sf $griddir/orog/${CCASE}_oro_data.tile7.halo*.nc .
  ln -sf $griddir/grid/${CCASE}_mosaic.nc               .

  rename ${CCASE} C768 ${CCASE}_*
#else
#  fn=$(ls $wrkdir/grid_orog/C*_grid.tile7.halo3.nc)
#  fm=${fn##$wrkdir/grid_orog/}
#  CCASE=${fm%%_grid.tile7.halo3.nc}
fi

for imn in $(seq $numsta 1 $numens); do
  ensmemid=$(printf %3.3i $imn)
  memdir="$wrkdir/mem_$ensmemid"
  #
  # Link GDAS
  #
  if [[ ! -e $memdir/GDAS ]]; then
    mkdir -p $memdir/GDAS
    cd $memdir/GDAS
    ln -s $gdasdir/${casedate}_mem$ensmemid/gdas.t${caseHH}z.ratmanl.mem$ensmemid.nemsio gfs.t${caseHH}z.atmanl.nemsio
    ln -s $gdasdir/${casedate}_mem$ensmemid/gdas.t${caseHH}z.sfcanl.mem$ensmemid.nemsio gfs.t${caseHH}z.sfcanl.nemsio
  fi

  cd $memdir

  cat > sed_chgres_ic << EOF
        s|BASE_GSM=.*|BASE_GSM=${rootdir}|g
        s|CDATE=.*|CDATE=${casedate}|g
        s|INIDIR=.*|INIDIR=$memdir/GDAS|g
        s|FIXfv3=.*|FIXfv3=$wrkdir/grid_orog|g
        s|DATA=.*|DATA=$memdir/tmp/wrk.chgres|g
EOF
  sed -f sed_chgres_ic ${rootdir}/scripts/driver_chgres_regional.c768 > run_chgres_ic.ksh

  sbatch run_chgres_ic.ksh

  rm -f sed_chgres_ic

  #sleep $((imn*10))
  numjobs=$(squeue -u $username | wc -l)
  while [[ $numjobs -gt $maxjobs ]]; do
    sleep 10
    numjobs=$(squeue -u $username | wc -l)
  done
done

