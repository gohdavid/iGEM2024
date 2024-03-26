#!/bin/bash
#PBS -N simulate_protein
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=16gb:ngpus=1:gpu_type=RTX6000

source ~/.bashrc
mamba activate md

/rds/general/user/d21/home/iGEM/scripts/simulate_protein_implicit.py \
-i /rds/general/user/d21/home/iGEM/PDBs/HRT_Single.pdb \
-o $(pwd) \
-t 1000 \
-w 1