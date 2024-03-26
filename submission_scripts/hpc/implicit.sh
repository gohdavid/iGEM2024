#!/bin/bash
#PBS -N simulate_protein
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=16gb:ngpus=1:gpu_type=RTX6000

# /home/gridsan/ddavid/MD/scripts/simulate_protein_implicit.py --i /home/gridsan/ddavid/MD/PDBs/HRT_Single.pdb --o /home/gridsan/ddavid/MD/simulations/1_HRT_I --t 500
/home/gridsan/ddavid/MD/scripts/simulate_protein_implicit.py --i /home/gridsan/ddavid/MD/PDBs/PhaC_Single.pdb --o /home/gridsan/ddavid/MD/simulations/1_PhaC_I --t 500