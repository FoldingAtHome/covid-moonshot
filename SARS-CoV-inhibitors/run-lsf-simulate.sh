#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 12:00
#BSUB -R "rusage[mem=3]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q gpuqueue
#BSUB -gpu "num=1:mode=shared:mps=no:j_exclusive=yes"
#BSUB -m "lt-gpu ls-gpu lu-gpu lp-gpu ld-gpu"
#BSUB -o %J.sars-simulate.out
#BSUB -J "sars-simulate[1-115]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

export PREFIX="SARS_2003_inhibitors"

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules $PREFIX.csv --index $JOBID --output $PREFIX-docked --simulate
