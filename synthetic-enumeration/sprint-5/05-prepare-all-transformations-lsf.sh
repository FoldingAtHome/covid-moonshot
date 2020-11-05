#!/bin/bash
#
# Prepare all transformations
#

#BSUB -P "testing"
#BSUB -J "mpro[1-1000]" 
#BSUB -n 1
#BSUB -R rusage[mem=3]
#BSUB -R span[hosts=1]
#BSUB -q gpuqueue
#BSUB -sp 1 # low priority. default is 12, max is 25
#BSUB -gpu num=1:j_exclusive=yes:mode=shared
#BSUB -W  02:00
#BSUB -m "ls-gpu lg-gpu lt-gpu lp-gpu lg-gpu lu-gpu ld-gpu"
#BSUB -o output/out_%I.stdout 
#BSUB -eo output/out_%I.stderr
##BSUB -cwd "/scratch/%U/%J"
#BSUB -L /bin/bash

# Make output directory if it doesn't exist
[ ! -d "/path/to/dir" ] && mkdir output

# quit on first error
set -e

source ~/.bashrc
OPENMM_CPU_THREADS=1
NUMEXPR_MAX_THREADS=1

cd $LS_SUBCWD

unset CUDA_OPENMM_COMPILER
conda activate perses

# Launch my program.
module load cuda/10.1
env | sort | grep 'CUDA'
env | sort | grep 'OPENMM'y
export RUN=$(expr $LSB_JOBINDEX-1)
python 05-prepare-single-transformation.py $RUN
