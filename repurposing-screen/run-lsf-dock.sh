#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.repurposing-dock.out
##BSUB -J "repurposing-dock[1-5000]"
#BSUB -J "repurposing-dock[5001-11335]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

let JOBID=$LSB_JOBINDEX-1
#python ../scripts/02-dock-and-prep.py --molecules broad-repurposing-library-20200324.csv --index $JOBID --output broad-repurposing-docked
python ../scripts/02-dock-and-prep.py --molecules drugbank-5.1.5-2020-01-03.csv --index $JOBID --output drugbank-5.1.5-2020-01-03
