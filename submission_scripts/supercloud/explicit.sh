#!/bin/bash

#SBATCH --job-name="tunnel"
#SBATCH --time=48:00:00
#SBATCH -p xeon-g6-volta
#SBATCH --gres=gpu:volta:1
#SBATCH --mem-per-cpu=9G

/home/gridsan/ddavid/MD/scripts/simulate_protein.py --i /home/gridsan/ddavid/MD/PDBs/HRT_Double_Close.pdb --o /home/gridsan/ddavid/MD/simulations/01_Double_Explicit_Close --t 200