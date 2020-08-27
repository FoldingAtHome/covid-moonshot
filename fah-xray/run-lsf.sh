#!/bin/bash
#BSUB -P "testing"
#BSUB -J "mpro[1-297]"
#BSUB -n 1
#BSUB -R rusage[mem=3]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W  05:59
#BSUB -m "ls-gpu lg-gpu lt-gpu lp-gpu lg-gpu lu-gpu ld-gpu"
#BSUB -o output/out_%I.stdout
#BSUB -eo output/out_%I.stderr
##BSUB -cwd "/scratch/%U/%J"
#BSUB -L /bin/bash

# quit on first error
set -e

source ~/.bashrc
OPENMM_CPU_THREADS=1

cd $LS_SUBCWD
conda activate perses

# Launch my program.
module load cuda/10.1
env | sort | grep 'CUDA'
python ../scripts/01-prep-xray-for-fah.py --run $LSB_JOBINDEX 
