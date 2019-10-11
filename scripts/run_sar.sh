#!/bin/bash
casedate=${1-2019052018}
numens=${2-4}

caseHH=${casedate:8:2}

rootdir="/oldscratch/ywang/EPIC/Program"
wrkdir="/oldscratch/ywang/EPIC/test_runs/${casedate}"
EXEPRO="${rootdir}/exec/fv3_32bit_cntl.exe"

griddir="$wrkdir/grid_orog"
ENDHOUR="6"

fn=$(realpath $wrkdir/grid_orog/C768_oro_data.tile7.halo0.nc)
if [[ ! -e $fn ]]; then
  echo "file $wrkdir/grid_orog/C768_oro_data.tile7.halo0.nc not exit."
  exit 0
fi
fm=$(basename $fn)
JCASE=${fm%%_oro_data.tile7.halo0.nc}

#
# Guss process number and block size
#
lonx=$(ncdump -h $fn | grep "lon =" | tr -dc '0-9')
latx=$(ncdump -h $fn | grep "lat =" | tr -dc '0-9')

npx=$((lonx+1))
npy=$((latx+1))

layout_x="14"
layout_y="10"
npes=$((${layout_x} * ${layout_y}+28))

n=$((npes/24))
m=$((npes%24))
if [[ $m -gt 0 ]]; then
  n=$((n+1))
fi

nx=$((lonx/layout_x))
ny=$((latx/layout_y))
npoints=$((nx*ny))
for b in $(seq 60 -1 1); do
  a=$((npoints%b))
  if [[ $a -eq 0 ]]; then
    break
  fi
done
#echo $n, $m, $npoints,$b
#echo $npx, $npy
#exit 0

#
# run fv3sar forecast for each member
#
for imn in $(seq 1 $numens); do
  ensmemid=$(printf %3.3i $imn)
  memdir="$wrkdir/mem_$ensmemid"

  cd $memdir

  if [[ -e tmp ]]; then
    rm -rf tmp
  fi

  cat > sed_sarfv3 << EOF
        s|FV3TTT|${rootdir}|g
        s|EXEPPP|$EXEPRO|g
        s|TTTDDD|${wrkdir}|g
        s|WWWDDD|$memdir|g
        s|EEEDDD|${casedate}|g
        s|CCCCCC|C768|g
        s|JJJJJJ|${JCASE}|g
        s|HHHHHH|${ENDHOUR}|g
        s|LYOXXX|${layout_x}|g
        s|LYOYYY|${layout_y}|g
        s|NPXXXX|${npx}|g
        s|NPYYYY|${npy}|g
        s|NPEEEE|${npes}|g
        s|BLOCKI|${b}|g
        s|NODESN|${n}|g
EOF

  jobscript="fv3_${casedate}.job"

  sed -f sed_sarfv3 ${rootdir}/templates/run_on_Odin.job > $jobscript

  echo "sbatch $jobscript"
  sbatch $jobscript

  rm -f sed_sarfv3
done

