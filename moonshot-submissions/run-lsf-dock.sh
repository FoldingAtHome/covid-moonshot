#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.moonshot-dock.out
##BSUB -J "moonshot-dock[1-2386]"
##BSUB -J "moonshot-dock[1-663]"
##BSUB -J "moonshot-dock[1-2016]"
#BSUB -J "moonshot-dock[1-3719]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

#export PREFIX="covid_submissions_all_info"
#export PREFIX="covalent_warhead_df"
#export PREFIX="nir-london-2020-03-07"
export PREFIX="covid_submissions_all_info"

let JOBID=$LSB_JOBINDEX-1
#python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules $PREFIX.csv --index $JOBID --output $PREFIX-docked --userfrags --covalent
python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules $PREFIX.csv --index $JOBID --output $PREFIX-docked --userfrags
