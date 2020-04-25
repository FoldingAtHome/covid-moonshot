#!/bin/bash

# Dock COVID Moonshot compounds in parallel

#BSUB -W 3:00
#BSUB -R "rusage[mem=2]"
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q cpuqueue
#BSUB -o %J.moonshot-dock.out
#BSUB -J "moonshot-dock[1-4373]"

echo "Job $JOBID/$NJOBS"

echo "LSB_HOSTS: $LSB_HOSTS"

source ~/.bashrc

source activate perses

#export PREFIX="covid_submissions_all_info"
#export PREFIX="covalent_warhead_df"
#export PREFIX="nir-london-2020-03-07"
#export PREFIX="covid_submissions_all_info-2020-04-06"
#export PREFIX="2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April"
#export PREFIX="COVID_MS_final_selection_round_2"
export PREFIX="covid_submissions_all_info-2020-04-23"

let JOBID=$LSB_JOBINDEX-1
python ../scripts/02-dock-and-prep.py --receptors ../receptors/monomer --molecules $PREFIX.csv --index $JOBID --output $PREFIX-docked --userfrags
