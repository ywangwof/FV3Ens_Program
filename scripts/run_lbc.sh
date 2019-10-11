#!/bin/bash
casedate=${1-2019052018}
numens=${2-4}

caseHH=${casedate:8:2}

rootdir="/oldscratch/ywang/EPIC/Program"

wrkdir="/oldscratch/ywang/EPIC/test_runs/${casedate}"
gdasdir="/oldscratch/ywang/EPIC/GDAS"

griddir="$wrkdir/grid_orog"

for imn in $(seq 1 $numens); do
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
  sed -f sed_chgres_lbc ${rootdir}/scripts/run_chgres.sh > run_chgres_lbc.sh

  chmod 755 run_chgres_lbc.sh
  ./run_chgres_lbc.sh

  rm -f sed_chgres_lbc
done

