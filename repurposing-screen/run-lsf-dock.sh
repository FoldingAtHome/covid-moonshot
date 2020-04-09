#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.repurposing-dock.out
##BSUB -J "repurposing-dock[1-5000]"
##BSUB -J "repurposing-dock[5001-11335]"

#BSUB -o %J.drugbank-dock.out
##BSUB -J "drugbank-dock[1-5000]"
##BSUB -J "drugbank-dock[5001-11335]"
#BSUB -J "drugset-dock[1-7000]"
##BSUB -J "drugset-dock[7001-14670]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

export PREFIX="14288275.moonshot-aggregate.out"
export PREFIX="drugbank-5.1.5-2020-01-03"
export PREFIX="drugset"

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --molecules $PREFIX.csv --index $JOBID --output $PREFIX
