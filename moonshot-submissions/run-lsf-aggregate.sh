#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.moonshot-aggregate.out
#BSUB -J "moonshot-aggregate"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

conda activate perses

export PREFIX="covid_submissions_all_info"

# Extract sorted docking results
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked-justscore.csv --clean
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.csv
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.sdf
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output $PREFIX-docked.pdb

# Extract for fragalysis
python ../scripts/03-aggregate-docking-results.py --molecules $PREFIX.csv --docked $PREFIX-docked --output compound-set_ensemble-hybrid-oedocking.sdf --clean --fragalysis "https://discuss.postera.ai/t/ensemble-oedocking-ensemble-hybrid-docking-to-fragment-bound-mpro-structures/1291"

# Prepare RUNs for FAH
python ../scripts/04-fah-prep.py --docked $PREFIX-docked --output $PREFIX-fah.csv

# Compute overlap
python ../scripts/05-score-heavy-atom-overlap.py --docked $PREFIX-docked.sdf --output $PREFIX-docked-overlap --clean --sort
