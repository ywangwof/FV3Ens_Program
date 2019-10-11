#!/bin/bash
casedate=${1-2019052018}
numens=${2-4}

caseHH=${casedate:8:2}

rootdir="/oldscratch/ywang/EPIC/Program"

wrkdir="/oldscratch/ywang/EPIC/test_runs/${casedate}"
gdasdir="/oldscratch/ywang/EPIC/GDAS"
griddir="/scratch/ywang/comFV3SAR/test_runs/caps_cntl"

#
# Link grid and orog files
#
if [[ ! -e $wrkdir/grid_orog ]]; then
  mkdir -p $wrkdir/grid_orog

  cd $wrkdir/grid_orog

  fn=$(ls $griddir/grid/C*_grid.tile7.halo3.nc)
  fm=${fn##$griddir/grid/}
  CCASE=${fm%%_grid.tile7.halo3.nc}

  ln -sf $griddir/grid/${CCASE}_grid.tile7.halo*.nc     .
  ln -sf $griddir/orog/${CCASE}_oro_data.tile7.halo*.nc .
  ln -sf $griddir/grid/${CCASE}_mosaic.nc               .

  rename ${CCASE} C768 ${CCASE}_*
else
  fn=$(ls $wrkdir/grid_orog/C*_grid.tile7.halo3.nc)
  fm=${fn##$wrkdir/grid_orog/}
  CCASE=${fm%%_grid.tile7.halo3.nc}
fi

for imn in $(seq 1 $numens); do
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
  sed -f sed_chgres_ic ${rootdir}/scripts/run_chgres_regional.c768 > run_chgres_ic.ksh

  sbatch run_chgres_ic.ksh

  rm -f sed_chgres_ic
done

