#!/bin/bash

if [[ ! $1 =~ ^[0-9]{10}$ ]]; then
  echo "$0 YYYYMMDDHH [NN] [test_run|test_mp|test_spp|test_mspp]"
  exit
fi

casedate=${1-2019052018}
numens=${2-40}
wrkcase=${3-"test_runs"}

wrkroot="/scratch/ywang/EPIC"

caseHH=${casedate:8:2}
eventymd=$(date -d "${casedate:0:8} ${caseHH}:00 1 hours ago" +%Y%m%d)

rootdir="/oldscratch/ywang/EPIC/Program"
wrkdir="${wrkroot}/${wrkcase}/${casedate}"

griddir="$wrkdir/grid_orog"
ENDHOUR="6"
quilting=".true."

templates=${rootdir}/templates
FIX_AM=${rootdir}/fix_am
CO2DIR=$FIX_AM/fix_co2_proj

FIXDIR=${wrkdir}/grid_orog

ymdh=${casedate}
ymd=`echo $ymdh |cut -c 1-8`
yyy=`echo $ymdh |cut -c 1-4`
mmm=`echo $ymdh |cut -c 5-6`
ddd=`echo $ymdh |cut -c 7-8`
hhh=`echo $ymdh |cut -c 9-10`

CASE="C768"
if [[ "$wrkcase" =~ "test_runs" || "$wrkcase" =~ "test_spp" ]]; then
  suitelen=1
elif [[ "$wrkcase" =~ "test_mp" || "$wrkcase" =~ "test_mspp" ]]; then
  suitelen=6           # number of suites to be used
else
  echo "Unknown wrkcase = \"$wrkcase\"."
  exit 1
fi

SUITES=("cntl" "lsm1" "sfcl1" "pbl1" "pbl2" "mp1")

hout_min=10                # minutes output
hout_hour=$(python3 <<<"print(${hout_min}/60.)")
hout_second=$((${hout_min}*60))

#-----------------------------------------------------------------------
#
# find equivalent grid RES
#
#-----------------------------------------------------------------------

fn=$(realpath $wrkdir/grid_orog/${CASE}_oro_data.tile7.halo0.nc)
if [[ ! -e $fn ]]; then
  echo "file $wrkdir/grid_orog/${CASE}_oro_data.tile7.halo0.nc not exit."
  exit 0
fi
fm=$(basename $fn)
JCASE=${fm%%_oro_data.tile7.halo0.nc}

#-----------------------------------------------------------------------
#
# Guss process number and block size
#
#-----------------------------------------------------------------------

lonx=$(ncdump -h $fn | grep "lon =" | tr -dc '0-9')
latx=$(ncdump -h $fn | grep "lat =" | tr -dc '0-9')

npx=$((lonx+1))
npy=$((latx+1))

layout_x="14"
layout_y="10"
nquilt="28"
io_layout="1,1"

if [[ "${quilting}" =~ ".true." ]]; then
  npes=$((layout_x * layout_y+nquilt))
else
  npes=$((layout_x * layout_y))
fi

nodes=$((npes/24))
m=$((npes%24))
if [[ $m -gt 0 ]]; then
  nodes=$((nodes+1))
fi

nx=$((lonx/layout_x))
ny=$((latx/layout_y))
npoints=$((nx*ny))

#
# find a suitable block size
#
for block in $(seq 60 -1 1); do
  a=$((npoints%block))
  if [[ $a -eq 0 ]]; then
    break
  fi
done

echo "npx      = $npx, npy = $npy, blocksize = $block"
echo "layout_x = $layout_x, layout_y = $layout_x, nquilt = ${nquilt} "
m=$((npes/nodes))
echo "npes     = ${npes}, nodes = $nodes, ncore_per_node = $m"

#-----------------------------------------------------------------------
#
# find centeral lat/lon
#
#-----------------------------------------------------------------------
fn="$wrkdir/grid_orog/${CASE}_grid.tile7.nc"

function ncattget {
  #ncks -m -M $1 | grep -E "^Global attribute [0-9]+: (CEN_LAT|CEN_LON|TRUELAT1|TRUELAT2|STAND_LON)" | cut -f 4,11 -d ' '
  /scratch/software/Odin/python/anaconda2/bin/ncks -x -M $1 | grep -E "(plat|plon)"
}

domains=$(ncattget $fn)

IFS=$'\n' domelements=($domains)
for var in ${domelements[@]}; do
  IFS='= ' keyval=(${var##Global attribute *:})
  wrfkey=${keyval[0]%%,}
  val=${keyval[-1]}

  case $wrfkey in
    plat)
      key="tlat"
      ;;
    plon)
      key="tlon"
      ;;
    *)
      key=${wrfkey}
      ;;
  esac
  declare "$key=$val"
  echo "$key -> $val"
done

unset IFS

echo "===================================="

#-----------------------------------------------------------------------
#
# run fv3sar forecast for each member
#
#-----------------------------------------------------------------------

l=0
for imn in $(seq 1 1 $numens); do

  ensmemid=$(printf "%03d" $imn)
  memdir="$wrkdir/mem_$ensmemid"
  suite=${SUITES[$l]}
  l=$(( (l+1)%suitelen ))

  if [[ "$wrkcase" =~ "spp" ]]; then
    EXEPRO="${rootdir}/exec/fv3_32bit_${suite}_spp.exe"
  else
    EXEPRO="${rootdir}/exec/fv3_32bit_${suite}.exe"
  fi

  cd $memdir

  if [[ -e tmp ]]; then
    rm -rf tmp
  fi

  mkdir -p RESTART
  ln -sf $EXEPRO fv3_gfs.x

  cp $templates/data_table .
  if [[ "$suite" == "mp1" ]]; then
    cp $templates/field_table.${suite} field_table
  else
    cp $templates/field_table .
  fi

  cp $templates/input.nml.$suite input.nml
  if [[ "$wrkcase" =~ "spp" ]]; then
    echo $wrkcase
    #cp $templates/input.nml.spp.v0  input.nml

    sed -i "/do_sppt/s/.F./.T./" input.nml
    sed -i "/do_shum/s/.F./.T./" input.nml
    sed -i "/do_skeb/s/.F./.T./" input.nml

    cat > spp.v0 <<SPPT
     new_lscale=.true.,
     sppt=1.0,
     sppt_tau = 2.16E4,
     sppt_lscale=150.E3,
     USE_ZMTNBLCK=.FALSE.,
     SPPT_LOGIT=.TRUE.,
     SPPT_SFCLIMIT=.FALSE.,
     ISEED_SPPT=$((casedate*1000 + imn*10 + 3)),
     SHUM=0.006,
     SHUM_TAU=21600,
     SHUM_LSCALE=150.E3,
     ISEED_SHUM=$((casedate*1000 + imn*10 + 2)),
     SKEBNORM=1,
     SKEB=0.5,
     SKEB_VDOF=10,
     SKEB_TAU=21600.,
     SKEB_LSCALE=150.E3,
     ISEED_SKEB=$((casedate*1000 + imn*10 + 1)),
     skebint=3600,
     spptint=3600,
     shumint=3600,
SPPT

    sed -i "/&nam_stochy/r spp.v0" input.nml

    #iseed_skeb=$((casedate*1000 + imn*10 + 1))
    #iseed_shum=$((casedate*1000 + imn*10 + 2))
    #iseed_sppt=$((casedate*1000 + imn*10 + 3))
    #sed -i "/ISEED_SKEB/s/[0-9]\+/$iseed_skeb/" input.nml
    #sed -i "/ISEED_SHUM/s/[0-9]\+/$iseed_shum/" input.nml
    #sed -i "/ISEED_SPPT/s/[0-9]\+/$iseed_sppt/" input.nml

    cp $templates/diag_table.spp diag_table
  else
    #cp $templates/input.nml.$suite input.nml
    cp $templates/diag_table .
  fi
  cp $templates/model_configure_${eventymd} model_configure
  cp $templates/nems.configure .
  cp $templates/CCN_ACTIVATE.BIN .
  cp $templates/global_h2oprdlos.f77 .
  cp $templates/global_o3prdlos.f77 .
  cp $templates/suite_CAPS_${suite}.xml .

  cp -p $FIX_AM/global_solarconstant_noaa_an.txt  solarconstant_noaa_an.txt
  cp -p $FIX_AM/global_o3prdlos.f77               global_o3prdlos.f77
  cp -p $FIX_AM/global_sfc_emissivity_idx.txt     sfc_emissivity_idx.txt
  cp -p $FIX_AM/global_co2historicaldata_glob.txt co2historicaldata_glob.txt
  cp -p $FIX_AM/co2monthlycyc.txt                 co2monthlycyc.txt
  cp -p $FIX_AM/global_climaeropac_global.txt     aerosol.dat

  for file in $(ls $CO2DIR/global_co2historicaldata*); do
    cp $file $(echo $(basename $file) |sed -e "s/global_//g")
  done

  #-----------------------------------------------------------------------
  # copy tile data and orography
  #-----------------------------------------------------------------------
  ln -sf $FIXDIR/${CASE}_grid.tile7.halo3.nc INPUT/
  ln -sf $FIXDIR/${CASE}_grid.tile7.halo4.nc INPUT/
  ln -sf $FIXDIR/${CASE}_mosaic.nc INPUT/

  ln -sf $FIXDIR/${CASE}_oro_data.tile7.halo0.nc INPUT/
  ln -sf $FIXDIR/${CASE}_oro_data.tile7.halo4.nc INPUT/

  cd INPUT
  ln -sf ${CASE}_grid.tile7.halo3.nc     ${JCASE}_grid.tile7.nc
  ln -sf ${CASE}_grid.tile7.halo4.nc     grid.tile7.halo4.nc
  ln -sf ${CASE}_oro_data.tile7.halo4.nc ${CASE}_oro_data.tile7.nc
  ln -sf ${CASE}_oro_data.tile7.halo4.nc oro_data.tile7.halo4.nc
  ln -sf ${CASE}_oro_data.tile7.halo4.nc oro_data.tile7.nc
  ln -sf ${CASE}_oro_data.tile7.halo4.nc oro_data.nc
  ln -sf ${CASE}_mosaic.nc               grid_spec.nc
  #
  # copy GFS datasets
  #
  ln -sf gfs_data.tile7.nc gfs_data.nc
  ln -sf sfc_data.tile7.nc sfc_data.nc

  cd $memdir

  #-----------------------------------------------------------------------
  # Modify local run-time files
  #-----------------------------------------------------------------------

  sed -i -e "/YYYYMMDD/s/YYYYMMDD/$ymd/" diag_table
  sed -i -e "/YYYY/s/YYYY/$yyy/"         diag_table
  sed -i -e "/MM/s/MM/$mmm/"             diag_table
  sed -i -e "/DD/s/DD/$ddd/"             diag_table
  sed -i -e "/HH/s/HH/$hhh/"             diag_table

  sed -i -e "/YYYY/s/YYYY/$yyy/" model_configure
  sed -i -e "/MM/s/MM/$mmm/"     model_configure
  sed -i -e "/DD/s/DD/$ddd/"     model_configure
  sed -i -e "/HH/s/HH/$hhh/"     model_configure

  sed -i -e "/NTASK/s/NTASK/$nquilt/;/quilting/s/.QUILT./$quilting/"  model_configure
  sed -i -e "/NPES/s/NPES/${npes}/;/nhours_fcst/s/ENDHOUR/${ENDHOUR}/"  model_configure
  sed -i -e "/HOUTMIN/s/HOUTMIN/$hout_min/"  model_configure

  sed -i -e "s/HOUTHRF/$hout_hour/;s/HOUTSEC/$hout_second/"  input.nml

  sed -i -e "/LAYOUT/s/LAYOUTX/${layout_x}/;s/LAYOUTY/${layout_y}/" input.nml
  sed -i -e "/LAYOUT/s/LAYOUTIO/${io_layout}/"                      input.nml
  sed -i -e "/npx/s/NPX/${npx}/;/npy/s/NPY/${npy}/;/blocksize/s/30/${block}/" input.nml

  sed -i -e "/FIX_AM/s#FIX_AM#${FIX_AM}#"      input.nml

  sed -i -e "/CENLAT/s/CENLAT/${tlat}/;s/CENLON/${tlon}/" input.nml

  #-----------------------------------------------------------------------
  # Create job file and submit it as needed
  #-----------------------------------------------------------------------

  cat > sed_sarfv3 << EOF
        s|TTTDDD|${wrkdir}|g
        s|MMMMMM|${ensmemid}|g
        s|NPEEEE|${npes}|g
        s|NODESN|${nodes}|g
EOF

  jobscript="run_fv3_${casedate}.job"

  sed -f sed_sarfv3 ${rootdir}/templates/run_fv3_on_Odin.job > $jobscript

  echo -n "Run fv3sar for memeber $ensmemid .... "
  sbatch $jobscript
  sleep 120

  rm -f sed_sarfv3
done

