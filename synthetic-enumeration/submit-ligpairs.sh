#!/bin/bash
#BSUB -P "testing"
##BSUB -J "mpro[1344-2688]" 
##BSUB -J "mpro[1-20]" 
#BSUB -J "mpro[2688-3359]" 
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

# quit on first error
set -e

source ~/.bashrc
OPENMM_CPU_THREADS=1
module load cuda/10.2
unset OPENMM_CUDA_COMPILER

cd $LS_SUBCWD

unset CUDA_OPENMM_COMPILER
conda activate perses

# Launch my program.
module load cuda/10.1
env | sort | grep 'CUDA'
export RUN=$(expr $LSB_JOBINDEX)
python run.py $RUN
