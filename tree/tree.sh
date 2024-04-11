#!/bin/bash
#SBATCH -A b1059
#SBATCH -p b1059
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="iqtree.oe"
#SBATCH --job-name="iqtree"

# Activate env
source activate	trees

# Run MAFFT on [ORTHO_ID].fa
# Path to [ORTHO_ID].fa] can be full or relative

mafft --auto OG0001886.fa > OG0001886.mafft

# Run IQTREE
# Path to [ORTHO_ID].mafft can be full or relative

iqtree -s OG0001886.mafft -nt 24 -bb 1000
