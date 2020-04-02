#!/bin/bash

# Dock Broad repurposing library


#BSUB -W 0:30
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q gpuqueue
#BSUB -gpu "num=1:mode=shared:mps=no:j_exclusive=yes"
#BSUB -m "lt-gpu ls-gpu lu-gpu lp-gpu ld-gpu"
#BSUB -o %J.repurposing-simulate.out
#BSUB -J "repurposing-simulate[10147]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --molecules broad-repurposing-library-20200324.csv --index $JOBID --output broad-repurposing-docked --simulate --transfer
