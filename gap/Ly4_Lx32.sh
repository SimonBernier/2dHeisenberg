#!/bin/bash
#SBATCH --account=def-kagarwal
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-00:00
#SBATCH --mail-user=simon.bernier@mail.mcgill.ca
#SBATCH --mail-type=ALL

./gap_2dHeis 32 4
