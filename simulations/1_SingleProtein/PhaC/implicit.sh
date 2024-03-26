#!/bin/bash
#PBS -N simulate_protein
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=16gb:ngpus=1:gpu_type=RTX6000

mamba activate md

/home/gridsan/ddavid/MD/scripts/simulate_protein_implicit.py\
-i /home/gridsan/ddavid/MD/PDBs/PhaC_Single.pdb\
-o $(pwd)\
-t 1000\
-w 1