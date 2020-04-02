#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.moonshot-dock.out
#BSUB -J "moonshot-dock[1-172]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --molecules merged-compounds-2020-04-01.csv --index $JOBID --output merged-compounds-2020-04-01
