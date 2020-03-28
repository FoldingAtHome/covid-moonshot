#!/bin/sh
#PBS -l walltime=24:00:00 # owlsnest::48:00:00, cb2rr::24:00:00
#PBS -N Ms0326_L-RL_1600-1800
#PBS -q large        # owlsnest::large/qcd, cb2rr:gpu
#PBS -l nodes=1:ppn=2 # owlsnest::1:2, cb2rr::1:4
#PBS -o MS0326_L-RL_1600-1800
#PBS

cd $PBS_O_WORKDIR

. ~/.bashrc # removing may cause 'conda init' error
conda activate openff     # activate env with: python==3.7 openmm==7.5.0 openforcefield openmmforcefields
module load cuda/10   # owlsnest::cuda/10, cb2rr::cuda/10.0.130

# ~/path/to/anaconda/env/python simulate.py start end gpu_id

~/packages/anaconda3/envs/openff/bin/python simulate.py 1600 1700 0 &

~/packages/anaconda3/envs/openff/bin/python simulate.py 1700 1800 1 &

wait
