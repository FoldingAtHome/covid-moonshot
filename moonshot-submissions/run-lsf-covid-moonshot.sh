#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q gpuqueue
#BSUB -gpu "num=1:mode=shared:mps=no:j_exclusive=yes"
##BSUB -m "lt-gpu ls-gpu lu-gpu lp-gpu ld-gpu"
##BSUB -q cpuqueue
#BSUB -o %J.moonshot.out
#BSUB -J "moonshot[1-1913]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --molecules ../molecules/covid_submissions_03_26_2020.csv --index $JOBID --output docked
