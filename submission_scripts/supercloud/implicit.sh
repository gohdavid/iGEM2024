#!/bin/bash

#SBATCH --job-name="tunnel"
#SBATCH --time=48:00:00
#SBATCH -p xeon-g6-volta
#SBATCH --gres=gpu:volta:1
#SBATCH --mem-per-cpu=9G

/home/gridsan/ddavid/MD/scripts/simulate_protein_implicit.py --i /home/gridsan/ddavid/MD/PDBs/HRT_Single.pdb --o /home/gridsan/ddavid/MD/simulations/1_HRT_I --t 500
/home/gridsan/ddavid/MD/scripts/simulate_protein_implicit.py --i /home/gridsan/ddavid/MD/PDBs/PhaC_Single.pdb --o /home/gridsan/ddavid/MD/simulations/1_PhaC_I --t 500