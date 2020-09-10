#!/bin/bash

if [[ ! $1 =~ ^[0-9]{10}$ ]]; then
  echo "$0 YYYYMMDDHH [test_runs|test_mp|test_spp|test_mspp] [NS] [NN]"
  exit
fi

casedate=${1-2019052018}
wrkcase=${2-"test_runs"}
numsta=${3-1}
numens=${4-40}
username="tg455890"
maxjobs=25


#wrkroot="/scratch/ywang/EPIC"           # Odin
wrkroot="/scratch/00315/tg455890/EPIC"   # stampede

caseHH=${casedate:8:2}

#rootdir="/oldscratch/ywang/EPIC/Program"              # Odin
rootdir="/work/00315/tg455890/stampede2/EPIC/Program"  # Stampede

wrkdir="${wrkroot}/${wrkcase}/${casedate}"
gdasdir="/work/00315/tg455890/stampede2/EPIC/GDAS"

for imn in $(seq $numsta 1 $numens); do
  ensmemid=$(printf %3.3i $imn)
  memdir="$wrkdir/mem_$ensmemid"
  #
  # Link GDAS
  #
  if [[ ! -e $memdir/GDAS ]]; then
    mkdir -p $memdir/GDAS
  fi
  cd $memdir/GDAS
  ln -sf $gdasdir/${casedate}_mem$ensmemid/gdas.t${caseHH}z.atmf003s.mem$ensmemid.nemsio gfs.t${caseHH}z.atmf003.nemsio
  ln -sf $gdasdir/${casedate}_mem$ensmemid/gdas.t${caseHH}z.atmf006s.mem$ensmemid.nemsio gfs.t${caseHH}z.atmf006.nemsio
  ln -sf $gdasdir/${casedate}_mem$ensmemid/gdas.t${caseHH}z.atmf009s.mem$ensmemid.nemsio gfs.t${caseHH}z.atmf009.nemsio

  cd $memdir

  cat > sed_chgres_lbc << EOF
        s|TOP_DIR=.*|TOP_DIR=${rootdir}|g
        s|CDATE=.*|CDATE=${casedate}|g
        s|WRK_DIR=.*|WRK_DIR=${memdir}|g
        s|FIXfv3=.*|FIXfv3=$wrkdir/grid_orog|g
EOF
  sed -f sed_chgres_lbc ${rootdir}/scripts/driver_chgres_lbc.ksh > run_chgres_lbc.sh

  chmod 755 run_chgres_lbc.sh
  ./run_chgres_lbc.sh

  rm -f sed_chgres_lbc

  #sleep $((imn*10))
  numjobs=$(squeue -u $username | wc -l)
  while [[ $numjobs -gt $maxjobs ]]; do
    sleep 10
    numjobs=$(squeue -u $username | wc -l)
  done
done

