#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.moonshot-dock.out
##BSUB -J "moonshot-dock[1-2386]"
#BSUB -J "moonshot-dock[1-461]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

let JOBID=$LSB_JOBINDEX-1
#python ../scripts/02-dock-and-prep.py --molecules covid_submissions_03_31_2020.csv --index $JOBID --output covid_submissions_03_31_2020 --userfrags
python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules covid_submissions_with_warhead_info.csv --index $JOBID --output covid_submissions_with_warhead_info-monomer --userfrags --covalent
