#!/bin/bash
#SBATCH --account=def-kagarwal
#SBATCH --mem-per-cpu=24G
#SBATCH --time=3-00:00
#SBATCH --mail-user=simon.bernier@mail.mcgill.ca
#SBATCH --mail-type=ALL

./gap_2dHeis 48 6
