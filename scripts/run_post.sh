#!/bin/sh -l

FV3SARDIR=${FV3SARDIR-/oldscratch/ywang/EPIC/Program}  #"/lfs3/projects/hpc-wof1/ywang/regional_fv3/fv3sar.mine"

#-----------------------------------------------------------------------
#
# This script runs the post-processor (UPP) on the NetCDF output files
# of the write component of the FV3SAR model.
#
#-----------------------------------------------------------------------
#
UPPFIX="${FV3SARDIR}/UPP_fix"
#UPPEXE="${FV3SARDIR}/exec"
CRTM_FIX="${FV3SARDIR}/CRTM_v2.2.3_fix"
ENDIAN="Big_Endian"

hostname=$(hostname)
case $hostname in
  odin?)
    template_job="${FV3SARDIR}/templates/run_upp_on_Odin.job"
    ;;
  fe*)
    template_job="${FV3SARDIR}/templates/run_upp_on_Jet.job"
    ;;
  *)
    echo "Unsupported machine: $hostname."
    exit
    ;;
esac

#
#-----------------------------------------------------------------------
#
# Save current shell options (in a global array).  Then set new options
# for this script/function.
#
#-----------------------------------------------------------------------
#
RUNDIR=$1      # $eventdir
casedate=$2

numens=${3-4}
sfhr=${4-0}

nodes1="1"
numprocess="6"
numthread="4"
platppn=$((numprocess/nodes1))
#npes="24"

#-----------------------------------------------------------------------
#
# Create directory (POSTPRD_DIR) in which to store post-processing out-
# put.  Also, create a temporary work directory (FHR_DIR) for the cur-
# rent output hour being processed.  FHR_DIR will be deleted later be-
# low after the processing for the current hour is complete.
#
#-----------------------------------------------------------------------
CDATE=${casedate:0:8}
CHH=${casedate:8:2}

for hr in $(seq $sfhr 1 6); do

  fhr=$(printf "%03d" $hr)

  for mid in $(seq 1 $numens); do
    memid=$(printf "%03d" $mid)
    memdir="$RUNDIR/mem_$memid"

    POSTPRD_DIR="$memdir/postprd"

    dyn_file=${memdir}/dynf${fhr}.nc
    phy_file=${memdir}/phyf${fhr}.nc

    wtime=0
    while [[ ! -f ${dyn_file} ]]; do
      sleep 20
      wtime=$(( wtime += 20 ))
      echo "Waiting ($wtime seconds) for ${dyn_file}"
    done

    while [[ ! -f ${phy_file} ]]; do
      sleep 10
      wtime=$(( wtime += 10 ))
      echo "Waiting ($wtime seconds) for ${phy_file}"
    done

    FHR_DIR="${POSTPRD_DIR}/$fhr"
    if [[ ! -r ${FHR_DIR} ]]; then
      mkdir -p ${FHR_DIR}
    fi

    cd ${FHR_DIR}


#-----------------------------------------------------------------------
#
# Create text file containing arguments to the post-processing executa-
# ble.
#
#-----------------------------------------------------------------------
    vhr=$((hr+CHH))

    POST_TIME=$( date -d "${CDATE} $vhr hours" +%Y%m%d%H%M )
    POST_YYYY=${POST_TIME:0:4}
    POST_MM=${POST_TIME:4:2}
    POST_DD=${POST_TIME:6:2}
    POST_HH=${POST_TIME:8:2}

    cat > itag <<EOF
${dyn_file}
netcdf
grib2
${POST_YYYY}-${POST_MM}-${POST_DD}_${POST_HH}:00:00
FV3R
${phy_file}

&NAMPGB
  KPO=6,PO=1000.,925.,850.,700.,500.,250.,
/
EOF

    rm -f fort.*

#-----------------------------------------------------------------------
#
# Stage files.
#
#-----------------------------------------------------------------------

    ln -sf $UPPFIX/nam_micro_lookup.dat ./eta_micro_lookup.dat
    #ln -s $UPPFIX/postxconfig-NT-fv3sar.txt ./postxconfig-NT.txt
    ln -sf $UPPFIX/postxconfig-NT-fv3sar-hwt2019.txt ./postxconfig-NT.txt
    ln -sf $UPPFIX/params_grib2_tbl_new ./params_grib2_tbl_new

#-----------------------------------------------------------------------
#
# CRTM fix files
#
#-----------------------------------------------------------------------

    spcCoeff_files=(imgr_g15.SpcCoeff.bin imgr_g13.SpcCoeff.bin imgr_g12.SpcCoeff.bin imgr_g11.SpcCoeff.bin \
      amsre_aqua.SpcCoeff.bin tmi_trmm.SpcCoeff.bin \
      ssmi_f13.SpcCoeff.bin ssmi_f14.SpcCoeff.bin ssmi_f15.SpcCoeff.bin ssmis_f16.SpcCoeff.bin \
      ssmis_f17.SpcCoeff.bin ssmis_f18.SpcCoeff.bin ssmis_f19.SpcCoeff.bin ssmis_f20.SpcCoeff.bin \
      seviri_m10.SpcCoeff.bin imgr_mt2.SpcCoeff.bin imgr_mt1r.SpcCoeff.bin \
      imgr_insat3d.SpcCoeff.bin abi_gr.SpcCoeff.bin abi_gr.SpcCoeff.bin )

    for fn in ${spcCoeff_files[@]}; do
      ln -sf ${CRTM_FIX}/SpcCoeff/${ENDIAN}/$fn .
    done

    tauCoeff_files=(imgr_g15.TauCoeff.bin imgr_g13.TauCoeff.bin imgr_g12.TauCoeff.bin imgr_g11.TauCoeff.bin \
        amsre_aqua.TauCoeff.bin tmi_trmm.TauCoeff.bin \
        ssmi_f13.TauCoeff.bin ssmi_f14.TauCoeff.bin ssmi_f15.TauCoeff.bin ssmis_f16.TauCoeff.bin \
        ssmis_f17.TauCoeff.bin ssmis_f18.TauCoeff.bin ssmis_f19.TauCoeff.bin ssmis_f20.TauCoeff.bin \
        seviri_m10.TauCoeff.bin imgr_mt2.TauCoeff.bin imgr_mt1r.TauCoeff.bin \
        imgr_insat3d.TauCoeff.bin abi_gr.TauCoeff.bin abi_gr.TauCoeff.bin)

    for fn in ${tauCoeff_files[@]}; do
      ln -sf ${CRTM_FIX}/TauCoeff/ODPS/${ENDIAN}/$fn .
    done

    cloudAndAerosol_files=(CloudCoeff/${ENDIAN}/CloudCoeff.bin AerosolCoeff/${ENDIAN}/AerosolCoeff.bin \
         EmisCoeff/IR_Land/SEcategory/${ENDIAN}/NPOESS.IRland.EmisCoeff.bin \
         EmisCoeff/IR_Snow/SEcategory/${ENDIAN}/NPOESS.IRsnow.EmisCoeff.bin \
         EmisCoeff/IR_Ice/SEcategory/${ENDIAN}/NPOESS.IRice.EmisCoeff.bin \
         EmisCoeff/IR_Water/${ENDIAN}/Nalli.IRwater.EmisCoeff.bin \
         EmisCoeff/MW_Water/${ENDIAN}/FASTEM6.MWwater.EmisCoeff.bin )

    for fn in ${cloudAndAerosol_files[@]}; do
      ln -sf ${CRTM_FIX}/$fn .
    done

#-----------------------------------------------------------------------
#
# Run the post-processor and move output files from FHR_DIR to POSTPRD_-
# DIR.
#
#-----------------------------------------------------------------------
    jobscript=run_upp_$fhr.job
    cp ${template_job} ${jobscript}
    sed -i -e "s#WWWDDD#${FHR_DIR}#;s#NNNNNN#${nodes1}#;s#MMMMMM#${memid}#;s#PPPPPP#${platppn}#g;s#TTTTTT#${numthread}#g;s#EEEEEE#${FV3SARDIR}#;s#DDDDDD#${CDATE}#;s#CCCHHH#CHH#;s#HHHHHH#${hr}#;" ${jobscript}

    sbatch $jobscript

  done
done

exit

