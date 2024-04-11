#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguestA
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --output="orthofinder.oe"
#SBATCH --job-name="orthofinder"

# $ source activate orthofinder_env
# $ conda install -c bioconda orthofinder
# $ conda install -c bioconda mafft
# $ conda install -c bioconda iqtree

# Activate env
source activate trees

# Run Orthofinder
# Path to [FASTA_DIR] can be full or relative
# [FASTA_DIR] contains all your sample FASTAs + GOI FASTAs

orthofinder -f 4sp_prot/ -og -t 12

# Orthofinder results will be in ~/[FASTA_DIR]/OrthoFinder/Results[DATE]
