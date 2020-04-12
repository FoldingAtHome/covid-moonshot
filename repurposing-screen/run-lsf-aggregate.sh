#!/bin/bash

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.repurposing-aggregate.out
#BSUB -J "repurposing-aggregate"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

#export PREFIX="drugbank-5.1.5-2020-01-03"
#export PREFIX="broad-repurposing-library-20200324"
export PREFIX="drugset"

# Extract sorted docking results
python ../scripts/03-aggregate-docking-results.py --docked $PREFIX --output $PREFIX-docked-justscore.csv --clean
python ../scripts/03-aggregate-docking-results.py --docked $PREFIX --output $PREFIX-docked.csv
python ../scripts/03-aggregate-docking-results.py --docked $PREFIX --output $PREFIX-docked.sdf
python ../scripts/03-aggregate-docking-results.py --docked $PREFIX --output $PREFIX-docked.pdb

# Coalesce RUNs for FAH
python ../scripts/04-fah-prep.py --docked $PREFIX --output $PREFIX-fah.csv
