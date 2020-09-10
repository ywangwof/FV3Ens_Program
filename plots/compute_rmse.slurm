#!/bin/bash
#SBATCH -A smallqueue
#SBATCH -p workq
#SBATCH -J compute_rmse
#SBATCH -N 1 -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH -t 04:30:00
#SBATCH -o /oldscratch/ywang/EPIC/Program/plots/compute_rmse_%j.out
#SBATCH -e /oldscratch/ywang/EPIC/Program/plots/compute_rmse_%j.err

time1=$(date '+%s')
echo "Job Started: $(date). Job Id:  $SLURM_JOBID"
echo " "

cd /oldscratch/ywang/EPIC/Program/plots

export LD_LIBRARY_PATH=/opt/gcc/6.1.0/snos/lib64:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin;
export PYTHONPATH=/home/yunheng.wang/tools/py-ncepbufr/lib/python3.6/site-packages;
source "/scratch/software/Odin/python/anaconda2/etc/profile.d/conda.sh"
conda activate wofs_post

# Set the stack limit as high as we can.
ulimit -s unlimited
#ulimit -a

set echo on

srun -n 1 compute_rmse_agrr.py

set echo off

time2=$(date '+%s')

let diff=time2-time1
let hour=diff/3600
let diff=diff%3600
let min=diff/60
let sec=diff%60

echo -n "Job   Ended: $(date). "
printf 'Job run time:  %02d:%02d:%02d' $hour $min $sec
echo " "

