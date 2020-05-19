#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.activity-aggregate.out
#BSUB -J "activity-aggregate"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

conda activate perses

export PREFIX="activity_data"

# Extract sorted docking results
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked-justscore.csv --clean
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.csv
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.sdf
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.pdb

# Generate ROC plots
python ../scripts/06-generate-ROC-plot.py --docked activity_data-docked.sdf 
