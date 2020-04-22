#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 0:30
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.london-dock.out
##BSUB -J "london-dock[1-7000]"
#BSUB -J "london-dock[7001-15278]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

export PREFIX="pyridine_urea"

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules $PREFIX.csv --index $JOBID --output $PREFIX-docked --userfrags
