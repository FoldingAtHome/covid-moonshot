#!/bin/bash

# Dock Broad repurposing library


#BSUB -W 6:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q gpuqueue
#BSUB -gpu "num=1:mode=shared:mps=no:j_exclusive=yes"
#BSUB -m "lt-gpu ls-gpu lu-gpu lp-gpu ld-gpu"
#BSUB -o %J.repurposing-simulate.out
#BSUB -J "repurposing-simulate[1-5000]"
##BSUB -J "repurposing-simulate[5001-10147]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

#export PREFIX="14288275.moonshot-aggregate.out"
#export PREFIX="drugbank-5.1.5-2020-01-03"
#export PREFIX="drugset"
export PREFIX="broad-repurposing-library-20200324"

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --molecules $PREFIX.csv --index $JOBID --output $PREFIX --simulate
